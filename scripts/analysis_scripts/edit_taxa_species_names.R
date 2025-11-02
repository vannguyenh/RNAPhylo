install.packages("rentrez")
library(rentrez)
library(ape)
library(phytools)

# accession → organism & taxid
s1 <- entrez_summary(db="nuccore", id="AJ278335.1")
s1$organism; s1$taxid

# taxid → scientific name
xml <- entrez_fetch(db="taxonomy", id="11120", rettype="xml")
sub(".*<ScientificName>([^<]+)</ScientificName>.*", "\\1", xml)

# accession → organism & taxid
s2 <- entrez_summary(db="nuccore", id="L06252.1")
s2$organism; s2$taxid

### test 2
pretty_name_from_acc <- function(acc) {
  s <- entrez_summary(db = "nuccore", id = acc)
  # This is the line you see on the NCBI result card:
  s$title
}

pretty_name_from_acc("AJ278335.1")
# "Avian infectious bronchitis virus, strain D207, 3' UTR"

pretty_name_from_acc("L06252.1")
# "Avian infectious bronchitis virus hypervariable region - related RNA sequence"

pretty_name_from_acc("AY641576.1")

#########################
# ---------------------------
# RF00164 relabelling script
# ---------------------------

# deps
# install.packages(c("ape","rentrez","stringr"), repos="https://cloud.r-project.org")
library(ape)
library(phytools)  # for nodeHeights (to align x-axes nicely)
library(rentrez)
library(stringr)

# (optional) be polite to NCBI:
# Sys.setenv(ENTREZ_EMAIL = "you@anu.edu.au")
# set_entrez_key("YOUR_EUTILS_API_KEY")  # speeds up & raises rate limits

# ---------- helpers ----------
# Your original cleaner: drop odd numeric prefixes / leading underscores,
# then keep up to the first "_" (=> accession + optional coords)
canon <- function(x) {
  x <- sub("^_*[0-9]+\\.?[0-9]*_", "", x)   # drop "64.8_"-style prefix
  x <- sub("^_*", "", x)                    # drop leading underscores
  sub("^([^_]+).*", "\\1", x)               # keep token before first "_"
}

# Extract "accession token" after canon(); then trim coords
# e.g. "Y15936.2/6871-6913" -> "Y15936.2"; "URS000080DEE1" -> "URS000080DEE1"
extract_accession <- function(tip) {
  tok <- canon(tip)
  sub("/.*$", "", tok)
}

# Is this a nuccore-like accession (GenBank/RefSeq)? (we'll only query those)
is_nuccore_like <- function(x) {
  grepl("^([A-Z]{1,2}\\d{5,}|[A-Z]{3}\\d{5,}|[A-Z]{4}\\d{2,}|[A-Z]{2}_\\d+|[A-Z]{4}_\\d+)\\.?\\d*$", x) &
    !grepl("^URS", x, ignore.case = TRUE)
}

# Pull the "taxon chunk" embedded in the original label (fallback for non-nuccore IDs)
# e.g. "_..._synthetic_construct[32630].1" -> "synthetic construct"
label_taxon_guess <- function(tip) {
  t <- str_match(tip, "_([^\\[]+?)\\[")[,2]   # between last "_" and "["
  t <- ifelse(is.na(t), "", t)
  t <- gsub("_+", " ", gsub("\\.+$", "", t))  # underscores -> space; strip trailing dots
  str_trim(t)
}

# Quote labels for Newick (escape inner apostrophes)
quote_newick <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")

# Preferred /source modifiers to append when disambiguating duplicates
PREFERRED_MODS <- c("strain","isolate","serotype","genotype","subtype",
                    "host","country","collection_date","common_name")

# Fetch organism and a preferred modifier from NCBI for a nuccore accession
fetch_label_from_ncbi <- function(acc) {
  s <- tryCatch(entrez_summary(db="nuccore", id=acc), error=function(e) NULL)
  if (is.null(s)) return(list(organism=NA_character_, modifier=NA_character_, title=NA_character_))
  # subname often comes as "a|b|c"; subtype is a parallel vector
  subname <- s$subname; subtype <- s$subtype
  if (!is.null(subname)) subname <- unlist(strsplit(subname, "\\|"))
  modifier <- NA_character_
  if (!is.null(subtype) && !is.null(subname)) {
    df <- data.frame(type=subtype, value=subname, stringsAsFactors = FALSE)
    for (k in PREFERRED_MODS) {
      v <- df$value[df$type == k]
      if (length(v) && nzchar(v[1])) { modifier <- v[1]; break }
    }
  }
  # fallback: try to pull 'strain X' from the Title if no modifier chosen
  if (is.na(modifier) || !nzchar(modifier)) {
    ttl <- s$title %||% NA_character_
    if (!is.na(ttl)) {
      m <- sub("^.*?\\bstrain\\b\\s*([^,;]+).*", "\\1", ttl, perl=TRUE)
      if (!identical(m, ttl)) modifier <- paste("strain", trimws(m))
    }
  }
  list(organism = s$organism %||% NA_character_,
       modifier = modifier,
       title    = s$title %||% NA_character_)
}
`%||%` <- function(a,b) if(!is.null(a)) a else b

# Build a single accession -> pretty label map from a set of tips
# * Queries NCBI only for nuccore-like accessions
# * For others (e.g., URS...), uses the taxon chunk (e.g., "synthetic construct")
build_label_map <- function(tips, sleep_sec = 0.34, always_append_modifier = FALSE) {
  accs <- unique(vapply(tips, extract_accession, "", USE.NAMES = FALSE))
  can_query <- is_nuccore_like(accs)
  
  org <- rep(NA_character_, length(accs))
  mod <- rep(NA_character_, length(accs))
  
  # Prepare a representative original label per accession (for fallback parsing)
  rep_by_acc <- tapply(tips, vapply(tips, extract_accession, "", USE.NAMES = FALSE), `[`, 1)
  
  # 1) Non-nuccore IDs: immediate fallback -> taxon chunk (e.g., "synthetic construct")
  for (i in which(!can_query)) {
    guess <- label_taxon_guess(rep_by_acc[[ accs[i] ]])
    org[i] <- if (nzchar(guess)) guess else accs[i]  # last resort: the accession token itself
    mod[i] <- NA_character_
  }
  
  # 2) Nuccore-like: query NCBI
  for (i in which(can_query)) {
    info <- fetch_label_from_ncbi(accs[i])
    org[i] <- info$organism
    mod[i] <- info$modifier
    if (sleep_sec > 0) Sys.sleep(sleep_sec)  # ~3 req/sec without API key
  }
  
  # 3) Build labels
  labels <- org
  
  if (always_append_modifier) {
    use_mod <- nzchar(ifelse(is.na(mod), "", mod))
    labels[use_mod] <- sprintf("%s, %s", org[use_mod], mod[use_mod])
  } else {
    # append modifier ONLY when organism duplicates exist
    dup_org <- duplicated(org) | duplicated(org, fromLast = TRUE)
    use_mod <- dup_org & nzchar(ifelse(is.na(mod), "", mod))
    labels[use_mod] <- sprintf("%s, %s", org[use_mod], mod[use_mod])
  }
  
  # If still duplicated, append short accession in brackets
  still_dup <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  shortacc  <- sub("\\..*$", "", accs)
  labels[still_dup] <- sprintf("%s [%s]", labels[still_dup], shortacc[still_dup])
  
  # Absolute uniqueness guard (rare)
  if (any(duplicated(labels))) labels <- make.unique(labels, sep = " _")
  
  data.frame(
    accession = accs,
    pretty_label = labels,
    pretty_label_quoted = quote_newick(labels),
    stringsAsFactors = FALSE
  )
}

# Apply the map to a tree
apply_label_map <- function(tr, map_df) {
  acc <- vapply(tr$tip.label, extract_accession, "", USE.NAMES = FALSE)
  idx <- match(acc, map_df$accession)
  newlabs <- ifelse(is.na(idx), tr$tip.label, map_df$pretty_label_quoted[idx])
  tr$tip.label <- newlabs
  tr
}

unquote_newick <- function(x) sub("^'(.*)'$", "\\1", x)
find_tip <- function(tr, target) {
  u <- unquote_newick(tr$tip.label)
  hit <- which(u == target)
  if (!length(hit)) stop(sprintf("Outgroup '%s' not found.", target))
  hit[1]
}
root_by_label <- function(tr, target) {
  tr <- root(tr, outgroup = find_tip(tr, target), resolve.root = TRUE)
  ladderize(tr, right = TRUE)
}
tree_depth <- function(tr) max(nodeHeights(tr))
fmt_support <- function(x, thr = 50) {
  y <- suppressWarnings(as.numeric(x))
  y <- ifelse(is.na(y), NA, ifelse(y <= 1, 100*y, y))
  if (!is.null(thr)) y <- ifelse(is.na(y) | y < thr, NA, y)
  ifelse(is.na(y), "", sprintf("%.0f", y))
}

# compute how much extra x-range we need for full labels (in branch-length units)
label_padding <- function(tr, cex_tip = 0.55, extra_in = 0.2) {
  # width of the longest label (inches)
  labs <- unquote_newick(tr$tip.label)
  w_in <- max(strwidth(labs, units = "inches", cex = cex_tip), na.rm = TRUE) + extra_in
  # convert inches -> user x-units after a temporary plot to get usr/pin
  plot.phylo(tr, type = "phylogram", plot = FALSE)  # sets last_plot.phylo env
  # do a tiny invisible plot to prime par("pin") safely
  op <- par(no.readonly = TRUE); plot.new(); on.exit(par(op), add = TRUE)
  dx_per_in <- diff(par("usr")[1:2]) / par("pin")[1]
  w_in * dx_per_in
}

# plot_with_bootstrap now takes absolute x_end and label_offset_units
plot_with_bootstrap <- function(tr, title_txt, x_end, label_offset_units,
                                cex_tip = 0.52, cex_boot = 0.6, thr = 50) {
  plot.phylo(tr, show.tip.label = TRUE, cex = cex_tip, use.edge.length = TRUE,
             x.lim = c(0, x_end),
             label.offset = label_offset_units,
             align.tip.label = TRUE, no.margin = FALSE)
  title(title_txt)
  if (!is.null(tr$node.label)) {
    # format 0–1 → percent; drop below threshold
    y <- suppressWarnings(as.numeric(tr$node.label))
    y <- ifelse(is.na(y), NA, ifelse(y <= 1, 100*y, y))
    if (!is.null(thr)) y <- ifelse(is.na(y) | y < thr, NA, y)
    labs <- ifelse(is.na(y), "", sprintf("%.0f", y))
    nodelabels(text = labs, frame = "n", cex = cex_boot, adj = c(0.5, -0.2), font = 2)
  }
}

# ---------- main: your files ----------

rfam_path <- '/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/RF00164.rfam.newick'
rna_path  <- '/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/RAxML_bipartitions.RF00164_RNA'
dna_path  <- '/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/RAxML_bipartitions.RF00164_DNA'

rfam_raw <- read.tree(rfam_path)
rna_raw  <- read.tree(rna_path)
dna_raw  <- read.tree(dna_path)

# Build one consistent label map from the union of all tips
all_tips  <- c(rfam_raw$tip.label, rna_raw$tip.label, dna_raw$tip.label)
label_map <- build_label_map(
  all_tips,
  sleep_sec = 0.34,               # set to 0 if you use an API key
  always_append_modifier = FALSE  # TRUE = always add ", strain X" when available
)

# Apply to each tree
rfam_lab <- apply_label_map(rfam_raw, label_map)
rna_lab  <- apply_label_map(rna_raw,  label_map)
dna_lab  <- apply_label_map(dna_raw,  label_map)

# Write relabeled trees
write.tree(rfam_lab, sub("\\.newick$|\\.nwk$|\\.tre$", ".labeled.nwk", rfam_path, ignore.case = TRUE))
write.tree(rna_lab,  paste0(rna_path,  ".labeled.nwk"))
write.tree(dna_lab,  paste0(dna_path,  ".labeled.nwk"))


# Also emit the mapping table for inspection/reuse
map_out <- '/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/label_map.csv'
write.csv(label_map[, c("accession","pretty_label")], map_out, row.names = FALSE)

message("Done.\nRelabeled trees written.\nLabel map: ", map_out)

# --- reroot at URS000080DEE1 ---
outgroup_name <- "URS000080DEE1"
rfam_root <- root_by_label(rfam_lab, outgroup_name)
rna_root  <- root_by_label(rna_lab,  outgroup_name)
dna_root  <- root_by_label(dna_lab,  outgroup_name)

# align x-axes (same base depth), but allow per-panel extra padding for long labels
base_depth <- max(tree_depth(rfam_root), tree_depth(rna_root), tree_depth(dna_root))

# ---------- EXPORT: PNG (wide, no clipping) ----------
png('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/three_trees_RF00164.png',
    width = 5400, height = 1800, res = 300)

par(mfrow = c(1,3), mar = c(2, 1.5, 3.2, 1.5), oma = c(0,0,0,0), xpd = FALSE)
cex_tip <- 0.52

# base depth shared across panels (so x-axes are comparable)
base_depth <- max(tree_depth(rfam_root), tree_depth(rna_root), tree_depth(dna_root))

# compute padding in *branch-length units* using panel width (inches)
panel_w_in <- par("pin")[1]  # width (inches) of ONE panel in current device
max_lab_in <- max(
  strwidth(unquote_newick(rfam_root$tip.label), units = "inches", cex = cex_tip),
  strwidth(unquote_newick(rna_root$tip.label),  units = "inches", cex = cex_tip),
  strwidth(unquote_newick(dna_root$tip.label),  units = "inches", cex = cex_tip),
  na.rm = TRUE
)
dx_per_in  <- base_depth / panel_w_in
pad_units  <- max_lab_in * dx_per_in
label_offset_units <- 0.02 * base_depth
x_end <- base_depth + pad_units + label_offset_units

plot_with_bootstrap(rfam_root, "RFAM seed (RF00164)", x_end, label_offset_units, cex_tip = cex_tip)
plot_with_bootstrap(rna_root,  "RNA tree (RAxML)",    x_end, label_offset_units, cex_tip = cex_tip)
plot_with_bootstrap(dna_root,  "DNA tree (RAxML)",    x_end, label_offset_units, cex_tip = cex_tip)
dev.off()

# ---------- EXPORT: PNG (wide, no clipping) ----------
png('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/three_trees_RF00164.png',
    width = 5400, height = 1800, res = 300)

par(mfrow = c(1,3), mar = c(2, 1.5, 3.2, 1.5), oma = c(0,0,0,0), xpd = FALSE)
cex_tip <- 0.52

# base depth shared across panels (so x-axes are comparable)
base_depth <- max(tree_depth(rfam_root), tree_depth(rna_root), tree_depth(dna_root))

# compute padding in *branch-length units* using panel width (inches)
panel_w_in <- par("pin")[1]  # width (inches) of ONE panel in current device
max_lab_in <- max(
  strwidth(unquote_newick(rfam_root$tip.label), units = "inches", cex = cex_tip),
  strwidth(unquote_newick(rna_root$tip.label),  units = "inches", cex = cex_tip),
  strwidth(unquote_newick(dna_root$tip.label),  units = "inches", cex = cex_tip),
  na.rm = TRUE
)
dx_per_in  <- base_depth / panel_w_in
pad_units  <- max_lab_in * dx_per_in
label_offset_units <- 0.02 * base_depth
x_end <- base_depth + pad_units + label_offset_units

plot_with_bootstrap(rfam_root, "RFAM seed (RF00164)", x_end, label_offset_units, cex_tip = cex_tip)
plot_with_bootstrap(rna_root,  "RNA tree (RAxML)",    x_end, label_offset_units, cex_tip = cex_tip)
plot_with_bootstrap(dna_root,  "DNA tree (RAxML)",    x_end, label_offset_units, cex_tip = cex_tip)
dev.off()


####################################
library(ape)
library(phytools)

# --- helpers: turn pretty labels into species names ---
unquote_newick <- function(x) sub("^'(.*)'$", "\\1", x)
tip_species <- function(tips) {
  # remove ", strain …" and any " [ACCESSION]" suffix
  sp <- unquote_newick(tips)
  sp <- sub(",.*$", "", sp)          # drop anything after first comma (strain/isolate)
  sp <- sub(" \\[.*\\]$", "", sp)    # drop trailing [ACC]
  trimws(sp)
}

# Collapse a tree to one tip per species (choose 1 rep per species),
# reporting which species are NOT monophyletic in the original tree.
collapse_to_species <- function(tr) {
  sp <- tip_species(tr$tip.label)
  split_by_sp <- split(tr$tip.label, sp)
  
  # check monophyly (on the original tree with quoted tip labels)
  non_mono <- names(Filter(function(v) length(v) > 1 && !is.monophyletic(tr, v),
                           split_by_sp))
  
  # pick a representative tip for each species (first one)
  reps <- vapply(split_by_sp, `[`, "", 1, USE.NAMES = TRUE)
  
  # keep only the representatives, then rename them to the species name
  tr_sp <- keep.tip(tr, reps)
  tr_sp$tip.label <- paste0("'", names(reps), "'")  # quote for Newick safety
  tr_sp <- ladderize(tr_sp, right = TRUE)
  
  list(tree = tr_sp, non_mono = non_mono)
}

# --- build species trees + reports ---
rfam_sp <- collapse_to_species(rfam_root)
rna_sp  <- collapse_to_species(rna_root)
dna_sp  <- collapse_to_species(dna_root)

cat("Non-monophyletic species (RFAM):", paste(rfam_sp$non_mono, collapse=", "), "\n")
cat("Non-monophyletic species (RNA): ", paste(rna_sp$non_mono,  collapse=", "), "\n")
cat("Non-monophyletic species (DNA): ", paste(dna_sp$non_mono,  collapse=", "), "\n")

# --- plot/export side-by-side species trees (no clipping) ---
species_png <- '/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/three_species_trees_RF00164.png'
#species_pdf <- '/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/three_species_trees_RF00164.pdf'

# shared x-scale
base_depth_sp <- max(max(nodeHeights(rfam_sp$tree)),
                     max(nodeHeights(rna_sp$tree)),
                     max(nodeHeights(dna_sp$tree)))

plot_species_panel <- function(tr, title_txt, x_end, label_offset_units,
                               cex_tip = 0.65, cex_boot = 0.7) {
  plot.phylo(tr, show.tip.label = TRUE, cex = cex_tip, use.edge.length = TRUE,
             x.lim = c(0, x_end),
             label.offset = label_offset_units,
             align.tip.label = TRUE, no.margin = FALSE)
  title(title_txt)
  # bootstrap labels (formatted as % if 0–1)
  if (!is.null(tr$node.label)) {
    y <- suppressWarnings(as.numeric(tr$node.label))
    y <- ifelse(is.na(y), NA, ifelse(y <= 1, 100*y, y))
    labs <- ifelse(is.na(y), "", sprintf("%.0f", y))
    nodelabels(text = labs, frame = "n", cex = cex_boot, adj = c(0.5, -0.2), font = 2)
  }
}

# open PNG and compute padding once per device
png(species_png, width = 4200, height = 1500, res = 300)
par(mfrow = c(1,3), mar = c(2, 1.5, 3.2, 1.5), oma = c(0,0,0,0))
cex_tip <- 0.65
panel_w_in <- par("pin")[1]
max_lab_in <- max(
  strwidth(unquote_newick(rfam_sp$tree$tip.label), units="inches", cex=cex_tip),
  strwidth(unquote_newick(rna_sp$tree$tip.label),  units="inches", cex=cex_tip),
  strwidth(unquote_newick(dna_sp$tree$tip.label),  units="inches", cex=cex_tip),
  na.rm = TRUE
)
dx_per_in <- base_depth_sp / panel_w_in
pad_units <- max_lab_in * dx_per_in
label_offset_units <- 0.02 * base_depth_sp
x_end <- base_depth_sp + pad_units + label_offset_units

plot_species_panel(rfam_sp$tree, "RFAM seed → species tree", x_end, label_offset_units, cex_tip)
plot_species_panel(rna_sp$tree,  "RNA (RAxML) → species tree", x_end, label_offset_units, cex_tip)
plot_species_panel(dna_sp$tree,  "DNA (RAxML) → species tree", x_end, label_offset_units, cex_tip)
dev.off()


########## RF00740
library(ape)
library(rentrez)
library(stringr)

# Sys.setenv(ENTREZ_EMAIL="you@anu.edu.au")  # recommended
# set_entrez_key("YOUR_API_KEY")             # optional (faster)

# ---- helpers ----
unquote_newick <- function(x) sub("^'(.*)'$", "\\1", x)

extract_taxid <- function(label) {
  m <- str_match(label, "\\[(\\d+)\\]")[,2]
  as.integer(m)
}

is_urs <- function(label) grepl("_URS[0-9A-F]+", label, ignore.case = TRUE) | grepl("^URS[0-9A-F]+", label, ignore.case = TRUE)

# fetch scientific/common name from NCBI taxonomy
fetch_taxon_names <- function(taxids) {
  taxids <- unique(na.omit(taxids))
  out <- data.frame(taxid = taxids, sci = NA_character_, common = NA_character_, stringsAsFactors = FALSE)
  for (i in seq_along(taxids)) {
    s <- tryCatch(entrez_summary(db="taxonomy", id=taxids[i]), error=function(e) NULL)
    if (!is.null(s)) {
      out$sci[i]    <- s$scientificname %||% NA_character_
      out$common[i] <- s$commonname     %||% NA_character_
    }
    Sys.sleep(0.34) # be polite unless you use an API key
  }
  out
}
`%||%` <- function(a,b) if (!is.null(a)) a else b

# build pretty labels for a vector of tip labels
labels_from_urs_or_keep <- function(tips, style=c("sci_common","sci_only")) {
  style <- match.arg(style)
  taxid <- extract_taxid(tips)
  has_urs <- is_urs(tips) & !is.na(taxid)
  
  # fetch names for the taxids we need
  lut <- fetch_taxon_names(taxid[has_urs])
  
  # start with original labels
  new <- tips
  
  # replace URS labels with "Genus species (common name)" or just sci name
  if (nrow(lut)) {
    for (i in which(has_urs)) {
      row <- lut[lut$taxid == taxid[i], ]
      if (nrow(row) == 1 && !is.na(row$sci)) {
        lab <- switch(style,
                      sci_common = if (!is.na(row$common) && nzchar(row$common))
                        sprintf("%s (%s)", row$sci, row$common) else row$sci,
                      sci_only   = row$sci)
        new[i] <- paste0("'", gsub("'", "''", lab, fixed = TRUE), "'")
      }
    }
  }
  
  # disambiguate exact duplicates by appending [taxid]
  uq <- make.unique(new, sep=" _")
  dup <- duplicated(uq) | duplicated(uq, fromLast=TRUE)
  new[dup] <- paste0(unquote_newick(new[dup]), " [", taxid[dup], "]") |> paste0("'") |> gsub("^'(.*)$", "'\\1", .)
  
  new
}

# ---- usage: relabel a tree ----
# tr <- read.tree("your_seed_tree.newick")
# tr$tip.label <- labels_from_urs_or_keep(tr$tip.label, style = "sci_common")
# write.tree(tr, "your_seed_tree_labeled.nwk")

## reroot: URS0000D514AB
# ============================
# Standardise taxon names across RFAM/DNA/RNA trees
# - Parses species names from RFAM labels like ".../…_Genus_species[9606].1"
# - Builds accession -> species-name map
# - Applies to all three trees (adds [ACC] only if needed for uniqueness)
# - Writes: *.labeled.nwk and a CSV mapping
# ============================

# install.packages(c("ape","stringr"), repos="https://cloud.r-project.org")
library(ape)
library(stringr)

# ---- INPUTS: change these ----
rfam_path <- "/path/to/RF00740.rfam.newick"
dna_path  <- "/path/to/RAxML_bipartitionsBranchLabels.RF00740_DNA"
rna_path  <- "/path/to/RAxML_bipartitionsBranchLabels.RF00740_RNA"
out_dir   <- dirname(rfam_path)

# ---- helpers ----

# extract all tip labels (tokens right before ":<number>")
tip_labels_from_newick <- function(txt) {
  m <- gregexpr("(?<=\\(|,)([^():;,]+?)(?=:\\-?\\d)", txt, perl=TRUE)
  labs <- regmatches(txt, m)[[1]]
  trimws(labs)
}

# drop odd leading numeric prefixes and underscores; keep token before first "_"
canon <- function(x) {
  x <- sub("^_*[0-9]+\\.?[0-9]*_", "", x)
  x <- sub("^_*", "", x)
  sub("^([^_]+).*", "\\1", x)
}

# accession token shared across trees (e.g., "URS000075E9DB" or "AY304470.1")
extract_accession <- function(label) {
  tok <- canon(label)
  sub("/.*$", "", tok)       # strip coordinates
}

# try to parse "taxon chunk" from an RFAM-like label
# primary pattern: after last '/' there is "..._<taxon>[", capture <taxon>
rfam_taxon_guess <- function(label) {
  m <- str_match(label, "/[^/_]*_([^[]+)\\[")
  tax <- if (!is.na(m[1,2])) m[1,2] else {
    # fallback: between last "_" and "[" (may capture just epithet; better than nothing)
    lb <- regexpr("\\[", label)
    if (lb > 0) {
      us <- regexpr("_[^_]*\\[", label)
      if (us > 0) substr(label, us + 1, lb - 1) else label
    } else label
  }
  tax <- gsub("_+", " ", tax)
  tax <- sub("\\.+$", "", tax)
  trimws(tax)
}

# quote labels for Newick (escape inner apostrophes)
quote_newick <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")

# Relabel a Newick text using mapping accession -> species
relabel_newick <- function(newick_text, mapping, disambig = TRUE, disambig_tag = c("ACC","taxid")) {
  disambig_tag <- match.arg(disambig_tag)
  tips <- tip_labels_from_newick(newick_text)
  accs <- vapply(tips, extract_accession, "", USE.NAMES = FALSE)
  
  # Build base new labels
  base <- vapply(accs, function(a) if (!is.na(mapping[a])) mapping[a] else NA_character_, "")
  # Fallback: if not in mapping (e.g., tip absent from RFAM), try to parse directly
  missing <- is.na(base) | base == ""
  if (any(missing)) {
    base[missing] <- vapply(tips[missing], rfam_taxon_guess, "")
  }
  
  base <- ifelse(is.na(base) | base == "", tips, base)  # final fallback: original tip
  
  # Disambiguate duplicates
  new_labels <- base
  if (disambig) {
    dups <- duplicated(base) | duplicated(base, fromLast = TRUE)
    if (any(dups)) {
      if (disambig_tag == "ACC") {
        new_labels[dups] <- sprintf("%s [%s]", base[dups], accs[dups])
      } else {
        # try to pull [taxid] from the original tip
        taxid <- str_match(tips, "\\[(\\d+)\\]")[,2]
        new_labels[dups] <- ifelse(!is.na(taxid[dups]),
                                   sprintf("%s [taxid:%s]", base[dups], taxid[dups]),
                                   sprintf("%s [%s]", base[dups], accs[dups]))
      }
    }
  }
  
  # Quote for Newick and splice back into the text
  quoted <- vapply(new_labels, quote_newick, "")
  pieces <- list(); last <- 0L
  it <- gregexpr("(?<=\\(|,)([^():;,]+?)(?=:\\-?\\d)", newick_text, perl=TRUE)[[1]]
  for (i in seq_along(it)) {
    s <- attr(it, "capture.start")[i,1]; e <- s + attr(it, "capture.length")[i,1] - 1
    pieces[[length(pieces)+1]] <- substr(newick_text, last + 1, s - 1)
    pieces[[length(pieces)+1]] <- quoted[i]
    last <- e
  }
  pieces[[length(pieces)+1]] <- substr(newick_text, last + 1, nchar(newick_text))
  paste(pieces, collapse = "")
}

# ---- read files ----
rfam_txt <- paste(readLines(rfam_path, warn=FALSE), collapse="\n")
dna_txt  <- paste(readLines(dna_path,  warn=FALSE), collapse="\n")
rna_txt  <- paste(readLines(rna_path,  warn=FALSE), collapse="\n")

# ---- build accession -> species mapping using RFAM labels ----
rfam_tips <- tip_labels_from_newick(rfam_txt)
accs <- vapply(rfam_tips, extract_accession, "", USE.NAMES = FALSE)
species <- vapply(rfam_tips, rfam_taxon_guess, "", USE.NAMES = FALSE)

acc2name <- setNames(species, accs)

# ---- apply to all three trees ----
rfam_labeled <- relabel_newick(rfam_txt, acc2name, disambig = TRUE, disambig_tag = "ACC")
dna_labeled  <- relabel_newick(dna_txt,  acc2name, disambig = TRUE, disambig_tag = "ACC")
rna_labeled  <- relabel_newick(rna_txt,  acc2name, disambig = TRUE, disambig_tag = "ACC")

# ---- write outputs ----
rfam_out <- file.path(out_dir, sub("\\.newick$|\\.nwk$|\\.tre$", ".labeled.nwk", basename(rfam_path), ignore.case = TRUE))
dna_out  <- file.path(out_dir, paste0(basename(dna_path), ".labeled.nwk"))
rna_out  <- file.path(out_dir, paste0(basename(rna_path), ".labeled.nwk"))
writeLines(rfam_labeled, rfam_out)
writeLines(dna_labeled,  dna_out)
writeLines(rna_labeled,  rna_out)

# mapping CSV (for audit)
map_csv <- file.path(out_dir, "accession_to_species_map.csv")
write.csv(data.frame(accession = names(acc2name), species = unname(acc2name)),
          map_csv, row.names = FALSE)

message("Done.\n  RFAM: ", rfam_out,
        "\n  DNA : ", dna_out,
        "\n  RNA : ", rna_out,
        "\n  Map : ", map_csv)
