# ============================
# Relabel → Reroot (Eptesicus fuscus) → Side-by-side plot (PNG + PDF)
# ============================

# install.packages(c("ape","stringr","phytools"), repos="https://cloud.r-project.org")
library(ape)
library(stringr)
library(phytools)

# ---- INPUTS: change these ----
rfam_path <- "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00740/RF00740.rfam.newick"
dna_path  <- "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00740/RAxML_bipartitions.RF00740_RNA"
rna_path  <- "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00740/RAxML_bipartitions.RF00740_DNA"
out_dir   <- dirname(rfam_path)

## ---- Helpers (base-R, robust) ----

# Grab all tip labels (tokens just before ":<number>")
tip_labels_from_newick <- function(txt){
  m <- gregexpr("(?<=\\(|,)([^():;,]+?)(?=:\\-?\\d)", txt, perl=TRUE)
  trimws(regmatches(txt, m)[[1]])
}

# Clean odd numeric prefixes
canon_head <- function(x){
  x <- sub("^_*[0-9]+\\.?[0-9]*_", "", x)  # drop things like "64.8_"
  sub("^_*", "", x)
}

# Accession token shared across trees
#  - URS: KEEP the full token up to "/" (e.g., "URS000075B235_9606")
#  - non-URS: token before first "_" (e.g., "AY304470.1"), then strip coords
extract_accession <- function(label){
  head <- canon_head(label)
  if (grepl("^URS[0-9A-F]+_", head, ignore.case = TRUE)) {
    sub("/.*$", "", head)
  } else {
    tok <- sub("^([^_]+).*", "\\1", head)
    sub("/.*$", "", tok)
  }
}

# Pull first integer inside square brackets [...] (NCBI taxid)
extract_taxid <- function(label){
  m <- regexpr("\\[(\\d+)\\]", label, perl=TRUE)
  if (m[1] == -1) return(NA_integer_)
  as.integer(substring(label, m[1]+1, m[1]+attr(m, "match.length")-2))
}

# Quote labels for Newick
quote_newick <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
unquote_newick <- function(x) sub("^'(.*)'$", "\\1", x)

# Fetch scientific names from NCBI taxonomy (safe & simple loop)
fetch_scientific_names <- function(taxids, sleep_sec = 0.34){
  taxids <- unique(na.omit(taxids))
  out <- setNames(rep(NA_character_, length(taxids)), taxids)
  if (!length(taxids)) return(out)
  for (id in taxids){
    s <- tryCatch(entrez_summary(db="taxonomy", id=id), error=function(e) NULL)
    if (!is.null(s) && !is.null(s$scientificname)) out[[as.character(id)]] <- s$scientificname
    if (sleep_sec > 0) Sys.sleep(sleep_sec)
  }
  out
}

# Build accession -> preferred label map from RFAM tips using taxid
build_label_map_from_rfam <- function(rfam_tips){
  accs   <- vapply(rfam_tips, extract_accession, "", USE.NAMES = FALSE)
  taxids <- vapply(rfam_tips, extract_taxid,     NA_integer_, USE.NAMES = FALSE)
  
  sci_by_taxid <- fetch_scientific_names(taxids)
  
  # minimal fallback if a tip lacks [taxid]
  guess_species <- function(lbl){
    lb <- regexpr("\\[", lbl, perl=TRUE)[1]
    s  <- if (lb == -1) lbl else substr(lbl, 1, lb-1)
    s  <- sub("^.*?/", "", s)   # after last '/'
    s  <- sub("^[^_]*_", "", s) # after first '_' in that tail
    s  <- gsub("_+", " ", s)
    s  <- sub("\\.+$", "", s)
    trimws(s)
  }
  
  pref <- character(length(rfam_tips))
  for (i in seq_along(rfam_tips)){
    taxid <- taxids[i]
    lbl   <- rfam_tips[i]
    if (!is.na(taxid) && nzchar(sci_by_taxid[[as.character(taxid)]])) {
      pref[i] <- sci_by_taxid[[as.character(taxid)]]
    } else {
      g <- guess_species(lbl)
      pref[i] <- if (nzchar(g)) g else lbl
    }
  }
  
  first_idx <- !duplicated(accs)
  list(
    acc2name  = setNames(pref[first_idx], accs[first_idx]),
    acc2taxid = setNames(taxids[first_idx], accs[first_idx])
  )
}

# Relabel a Newick text using accession -> name map; disambiguate duplicates only if needed
relabel_newick <- function(newick_text, acc2name, acc2taxid=NULL, disambig_tag=c("acc","taxid")){
  disambig_tag <- match.arg(disambig_tag)
  tips <- tip_labels_from_newick(newick_text)
  accs <- vapply(tips, extract_accession, "", USE.NAMES = FALSE)
  
  base <- vapply(accs, function(a) if (!is.na(acc2name[[a]])) acc2name[[a]] else NA_character_, "")
  base[is.na(base) | base==""] <- tips[is.na(base) | base==""]
  
  new_labels <- base
  dups <- duplicated(base) | duplicated(base, fromLast = TRUE)
  if (any(dups)){
    if (disambig_tag == "taxid" && !is.null(acc2taxid)){
      tdx <- vapply(accs, function(a) acc2taxid[[a]], NA_integer_)
      tag <- ifelse(!is.na(tdx), paste0("taxid:", tdx), accs)
      new_labels[dups] <- sprintf("%s [%s]", base[dups], tag[dups])
    } else {
      new_labels[dups] <- sprintf("%s [%s]", base[dups], accs[dups])
    }
  }
  
  quoted <- vapply(new_labels, quote_newick, "")
  pieces <- list(); last <- 0L
  it   <- gregexpr("(?<=\\(|,)([^():;,]+?)(?=:\\-?\\d)", newick_text, perl=TRUE)[[1]]
  caps <- attr(it, "capture.start"); lens <- attr(it, "capture.length")
  for (i in seq_along(it)){
    s <- caps[i,1]; e <- s + lens[i,1] - 1
    pieces[[length(pieces)+1]] <- substr(newick_text, last + 1, s - 1)
    pieces[[length(pieces)+1]] <- quoted[i]
    last <- e
  }
  pieces[[length(pieces)+1]] <- substr(newick_text, last + 1, nchar(newick_text))
  paste(pieces, collapse = "")
}

# find tip by name prefix (works if labels became "Eptesicus fuscus [ACC]")
find_tip_prefix <- function(tr, prefix){
  u <- unquote_newick(tr$tip.label)
  hit <- which(startsWith(u, prefix))
  if (!length(hit)) stop(sprintf("Outgroup '%s' not found.", prefix))
  hit[1]
}

# Common padding-aware plotting
tree_depth <- function(tr) max(nodeHeights(tr))
plot_with_padding <- function(tr, title_txt, base_depth, cex_tip=0.65, cex_boot=0.7, thr=NULL){
  panel_w_in <- par("pin")[1]
  max_lab_in <- max(strwidth(unquote_newick(tr$tip.label), units="inches", cex=cex_tip), na.rm=TRUE)
  dx_per_in  <- base_depth / panel_w_in
  pad_units  <- max_lab_in * dx_per_in
  label_offset_units <- 0.02 * base_depth
  x_end <- base_depth + pad_units + label_offset_units
  
  plot.phylo(tr, show.tip.label=TRUE, cex=cex_tip, use.edge.length=TRUE,
             x.lim=c(0, x_end), label.offset=label_offset_units,
             align.tip.label=TRUE, no.margin=FALSE)
  title(title_txt)
  if (!is.null(tr$node.label) && !all(is.na(tr$node.label))){
    y <- suppressWarnings(as.numeric(tr$node.label))
    y <- ifelse(is.na(y), NA, ifelse(y <= 1, 100*y, y))
    if (!is.null(thr)) y <- ifelse(is.na(y) | y < thr, NA, y)
    labs <- ifelse(is.na(y), "", sprintf("%.0f", y))
    nodelabels(text=labs, frame="n", cex=cex_boot, adj=c(0.5, -0.2), font=2)
  }
}

# ---------- 1) Read trees ----------
rfam_txt <- paste(readLines(rfam_path, warn=FALSE), collapse="\n")
dna_txt  <- paste(readLines(dna_path,  warn=FALSE), collapse="\n")
rna_txt  <- paste(readLines(rna_path,  warn=FALSE), collapse="\n")

# ---------- 2) Build label map from RFAM tips ----------
rfam_tips <- tip_labels_from_newick(rfam_txt)
lut <- build_label_map_from_rfam(rfam_tips)
acc2name  <- lut$acc2name
acc2taxid <- lut$acc2taxid

# sanity: make sure URS tokens for human & orangutan stay distinct
# print(acc2name[c("URS000075B235_9606","URS000075B235_9600")])

# ---------- 3) Relabel all three trees ----------
rfam_labeled_txt <- relabel_newick(rfam_txt, acc2name, acc2taxid, disambig_tag="acc")
dna_labeled_txt  <- relabel_newick(dna_txt,  acc2name, acc2taxid, disambig_tag="acc")
rna_labeled_txt  <- relabel_newick(rna_txt,  acc2name, acc2taxid, disambig_tag="acc")

# write out relabeled trees + mapping
rfam_out <- file.path(out_dir, sub("\\.newick$|\\.nwk$|\\.tre$", ".labeled.nwk",
                                   basename(rfam_path), ignore.case=TRUE))
dna_out  <- file.path(out_dir, paste0(basename(dna_path), ".labeled.nwk"))
rna_out  <- file.path(out_dir, paste0(basename(rna_path), ".labeled.nwk"))
writeLines(rfam_labeled_txt, rfam_out)
writeLines(dna_labeled_txt,  dna_out)
writeLines(rna_labeled_txt,  rna_out)

map_csv <- file.path(out_dir, "accession_to_scientific_name.csv")
write.csv(data.frame(accession=names(acc2name),
                     scientific_name=unname(acc2name),
                     taxid=unname(acc2taxid[names(acc2name)])),
          map_csv, row.names=FALSE)

# ---------- 4) Reroot at Eptesicus fuscus and ladderize ----------
rfam_tr <- read.tree(text=rfam_labeled_txt)
dna_tr  <- read.tree(text=dna_labeled_txt)
rna_tr  <- read.tree(text=rna_labeled_txt)

out_prefix <- "Eptesicus fuscus"
rfam_root <- ladderize(root(rfam_tr, outgroup=find_tip_prefix(rfam_tr, out_prefix), resolve.root=TRUE), right=TRUE)
dna_root  <- ladderize(root(dna_tr,  outgroup=find_tip_prefix(dna_tr,  out_prefix), resolve.root=TRUE), right=TRUE)
rna_root  <- ladderize(root(rna_tr,  outgroup=find_tip_prefix(rna_tr,  out_prefix), resolve.root=TRUE), right=TRUE)

# ---------- 5) Side-by-side export ----------
base_depth <- max(tree_depth(rfam_root), tree_depth(dna_root), tree_depth(rna_root))

png(file.path(out_dir, "RF00740_three_trees_rerooted.png"),
    width=5200, height=1700, res=300)
par(mfrow=c(1,3), mar=c(2,1.5,3.2,1.5), oma=c(0,0,0,0), xpd=FALSE)
plot_with_padding(rfam_root, "RFAM seed (scientific names; outgroup = Eptesicus fuscus)", base_depth, cex_tip=0.65)
plot_with_padding(rna_root,  "RNA (scientific names; rerooted)",                             base_depth, cex_tip=0.65)
plot_with_padding(dna_root,  "DNA (scientific names; rerooted)",                             base_depth, cex_tip=0.65)
dev.off()

# ---------- (Extra) Topology-only export: ignore branch lengths ----------
# Plot helper (no edge lengths). Scales padding to the longest label so nothing truncates.
plot_topology_only <- function(tr, title_txt, cex_tip = 0.65, cex_boot = 0.7, thr = NULL){
  panel_w_in <- par("pin")[1]
  max_lab_in <- max(strwidth(unquote_newick(tr$tip.label), units = "inches", cex = cex_tip), na.rm = TRUE)
  dx_per_in  <- 1 / panel_w_in           # x-range is [0,1] when use.edge.length = FALSE
  pad_units  <- max_lab_in * dx_per_in
  label_offset_units <- 0.02
  x_end <- 1 + pad_units + label_offset_units
  
  plot.phylo(tr, show.tip.label = TRUE, cex = cex_tip,
             use.edge.length = FALSE,           # <- the key bit: ignore branch lengths
             x.lim = c(0, x_end),
             label.offset = label_offset_units,
             align.tip.label = TRUE, no.margin = FALSE)
  title(title_txt)
  
  if (!is.null(tr$node.label) && !all(is.na(tr$node.label))){
    y <- suppressWarnings(as.numeric(tr$node.label))
    y <- ifelse(is.na(y), NA, ifelse(y <= 1, 100*y, y))  # handle 0–1 supports
    if (!is.null(thr)) y <- ifelse(is.na(y) | y < thr, NA, y)
    labs <- ifelse(is.na(y), "", sprintf("%.0f", y))
    nodelabels(text = labs, frame = "n", cex = cex_boot, adj = c(0.5, -0.2), font = 2)
  }
}

# PNG
png(file.path(out_dir, "RF00740_three_trees_rerooted_noscale.png"),
    width = 5200, height = 1700, res = 300)
par(mfrow = c(1,3), mar = c(2, 1.5, 3.2, 1.5), oma = c(0,0,0,0), xpd = FALSE)
plot_topology_only(rfam_root, "RFAM seed (no branch lengths; outgroup = Eptesicus fuscus)")
plot_topology_only(rna_root,  "RNA (no branch lengths)")
plot_topology_only(dna_root,  "DNA (no branch lengths)")
dev.off()
