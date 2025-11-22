# ============================
# RF00740: Relabel -> Reroot (Eptesicus fuscus) -> Side-by-side plots
# ============================

# install.packages(c("ape","phytools","stringr","rentrez"), repos="https://cloud.r-project.org")
library(ape)
library(phytools)
library(stringr)
library(rentrez)

## Optional (be polite to NCBI, and speed up with an API key if you have one):
# Sys.setenv(ENTREZ_EMAIL = "you@anu.edu.au")
# set_entrez_key("YOUR_NCBI_API_KEY")

# ---- INPUTS ----
rfam_path <- "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00740/RF00740.rfam.newick"
dna_path  <- "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00740/RAxML_bipartitions.RF00740_DNA"
rna_path  <- "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00740/RAxML_bipartitions.RF00740_S16"
out_dir   <- dirname(rfam_path)

# ---- Helpers ----

# get all tip labels (token just before :<number>)
tip_labels_from_newick <- function(txt){
  m <- gregexpr("(?<=\\(|,)([^():;,]+?)(?=:\\-?\\d)", txt, perl=TRUE)
  trimws(regmatches(txt, m)[[1]])
}

# strip odd numeric prefixes
canon_head <- function(x){
  x <- sub("^_*[0-9]+\\.?[0-9]*_", "", x)
  sub("^_*", "", x)
}

# accession token shared across trees
# - URS: keep whole "URS..._taxid" up to slash; if no _taxid, the core "URS..."
# - others: token before first "_" (e.g., AY304470.1) then drop coords
extract_accession <- function(label){
  head <- canon_head(label)
  if (grepl("^URS[0-9A-F]+_", head, ignore.case = TRUE)) {
    sub("/.*$", "", head)
  } else if (grepl("^URS[0-9A-F]+$", head, ignore.case = TRUE)) {
    head
  } else {
    tok <- sub("^([^_]+).*", "\\1", head)
    sub("/.*$", "", tok)
  }
}

# taxid inside [12345]
extract_taxid <- function(label){
  m <- regexpr("\\[(\\d+)\\]", label, perl=TRUE)
  if (m[1] == -1) return(NA_integer_)
  as.integer(substring(label, m[1]+1, m[1]+attr(m, "match.length")-2))
}

# fallback species text from the label (keeps genus + species)
guess_species <- function(lbl){
  # take chunk after last '/', drop trailing [taxid], underscores->spaces
  s <- sub(".*?/(.*)", "\\1", lbl)
  s <- sub("\\[\\d+\\].*$", "", s)
  s <- gsub("_+", " ", s)
  trimws(s)
}

# NCBI taxonomy lookup
lookup_taxid_name <- function(tid){
  s <- tryCatch(entrez_summary(db="taxonomy", id=tid), error=function(e) NULL)
  if (is.null(s)) NA_character_ else s$scientificname
}

# quote/unquote for Newick
quote_newick   <- function(x) paste0("'", gsub("'", "''", x, fixed=TRUE), "'")
unquote_newick <- function(x) sub("^'(.*)'$", "\\1", x)

# build accession->name map from RFAM; add URS core aliases
build_label_map_from_rfam <- function(rfam_tips){
  accs   <- vapply(rfam_tips, extract_accession, "", USE.NAMES=FALSE)
  taxids <- vapply(rfam_tips, extract_taxid,     NA_integer_, USE.NAMES=FALSE)
  
  sci <- character(length(rfam_tips))
  for (i in seq_along(rfam_tips)){
    if (!is.na(taxids[i])) sci[i] <- lookup_taxid_name(taxids[i])
    if (!nzchar(sci[i]))   sci[i] <- guess_species(rfam_tips[i])
  }
  
  keep       <- !duplicated(accs)
  acc2name   <- setNames(sci[keep],    accs[keep])
  acc2taxid  <- setNames(taxids[keep], accs[keep])
  
  # add aliases: URSxxxxxxxxxxxx_9606 -> URSxxxxxxxxxxxx
  is_urs_full <- grepl("^URS[0-9A-F]+_\\d+$", names(acc2name))
  if (any(is_urs_full)){
    urs_core <- sub("^(URS[0-9A-F]+)_\\d+$", "\\1", names(acc2name)[is_urs_full])
    acc2name  <- c(acc2name,  setNames(unname(acc2name [is_urs_full]), urs_core))
    acc2taxid <- c(acc2taxid, setNames(unname(acc2taxid[is_urs_full]), urs_core))
  }
  
  list(acc2name=acc2name, acc2taxid=acc2taxid)
}

# relabel a Newick using accession->name map; disambiguate duplicates if needed
relabel_newick <- function(newick_text, acc2name, acc2taxid=NULL, disambig_tag=c("acc","taxid")){
  disambig_tag <- match.arg(disambig_tag)
  tips <- tip_labels_from_newick(newick_text)
  accs <- vapply(tips, extract_accession, "", USE.NAMES=FALSE)
  
  base <- vapply(accs, function(a) if (!is.na(acc2name[[a]])) acc2name[[a]] else NA_character_, "")
  base[is.na(base) | base==""] <- tips[is.na(base) | base==""]
  
  new_labels <- base
  dups <- duplicated(base) | duplicated(base, fromLast=TRUE)
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

# find tip by prefix (works if labels become "Eptesicus fuscus [ACC]")
find_tip_prefix <- function(tr, prefix){
  u <- unquote_newick(tr$tip.label)
  hit <- which(startsWith(u, prefix))
  if (!length(hit)) stop(sprintf("Outgroup '%s' not found.", prefix))
  hit[1]
}

# plotting helpers
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

plot_topology_only <- function(tr, title_txt, cex_tip=0.65, cex_boot=0.7, thr=NULL){
  panel_w_in <- par("pin")[1]
  max_lab_in <- max(strwidth(unquote_newick(tr$tip.label), units="inches", cex=cex_tip), na.rm=TRUE)
  dx_per_in  <- 1 / panel_w_in
  pad_units  <- max_lab_in * dx_per_in
  label_offset_units <- 0.02
  x_end <- 1 + pad_units + label_offset_units
  
  plot.phylo(tr, show.tip.label=TRUE, cex=cex_tip, use.edge.length=FALSE,
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

# ---- 1) Read trees as text ----
rfam_txt <- paste(readLines(rfam_path, warn=FALSE), collapse="\n")
dna_txt  <- paste(readLines(dna_path,  warn=FALSE), collapse="\n")
rna_txt  <- paste(readLines(rna_path,  warn=FALSE), collapse="\n")

# ---- 2) Build label map from RFAM tips (with URS aliasing) ----
rfam_tips <- tip_labels_from_newick(rfam_txt)
lut <- build_label_map_from_rfam(rfam_tips)
acc2name  <- lut$acc2name
acc2taxid <- lut$acc2taxid

# (Optional) quick sanity:
# print(acc2name[c("URS000075E9DB_9796","URS000075E9DB",
#                  "URS000075B235_9606","URS000075B235")])

# ---- 3) Relabel all three trees ----
rfam_labeled_txt <- relabel_newick(rfam_txt, acc2name, acc2taxid, disambig_tag="acc")
dna_labeled_txt  <- relabel_newick(dna_txt,  acc2name, acc2taxid, disambig_tag="acc")
rna_labeled_txt  <- relabel_newick(rna_txt,  acc2name, acc2taxid, disambig_tag="acc")

# write relabeled trees + mapping
rfam_out <- file.path(out_dir, sub("\\.newick$|\\.nwk$|\\.tre$", ".labeled.newick",
                                   basename(rfam_path), ignore.case=TRUE))
dna_out  <- file.path(out_dir, paste0(basename(dna_path), ".labeled.nwk"))
rna_out  <- file.path(out_dir, paste0(basename(rna_path), ".labeled.nwk"))
writeLines(rfam_labeled_txt, rfam_out)
writeLines(dna_labeled_txt,  dna_out)
writeLines(rna_labeled_txt,  rna_out)

map_csv <- file.path(out_dir, "accession_to_scientific_name.csv")
write.csv(
  data.frame(
    accession       = names(acc2name),
    scientific_name = unname(acc2name),
    taxid           = unname(acc2taxid[names(acc2name)])
  ),
  map_csv, row.names = FALSE
)

# ---- 4) Reroot at Eptesicus fuscus & ladderize ----
rfam_tr <- read.tree(text=rfam_labeled_txt)
dna_tr  <- read.tree(text=dna_labeled_txt)
rna_tr  <- read.tree(text=rna_labeled_txt)

out_prefix <- "Eptesicus fuscus"
rfam_root <- ladderize(root(rfam_tr, outgroup=find_tip_prefix(rfam_tr, out_prefix), resolve.root=TRUE), right=TRUE)
dna_root  <- ladderize(root(dna_tr,  outgroup=find_tip_prefix(dna_tr,  out_prefix), resolve.root=TRUE), right=TRUE)
rna_root  <- ladderize(root(rna_tr,  outgroup=find_tip_prefix(rna_tr,  out_prefix), resolve.root=TRUE), right=TRUE)

# ---- 5) Side-by-side export (scaled) ----
base_depth <- max(tree_depth(rfam_root), tree_depth(dna_root), tree_depth(rna_root))

png(file.path(out_dir, "RF00740_three_trees_rerooted.png"),
    width=5200, height=1700, res=300)
par(mfrow=c(1,3), mar=c(2,1.5,3.2,1.5), oma=c(0,0,0,0), xpd=FALSE)
plot_with_padding(rfam_root, "RFAM seed (scientific names; outgroup = Eptesicus fuscus)", base_depth, cex_tip=0.65)
plot_with_padding(rna_root,  "RNA (scientific names; rerooted)",                             base_depth, cex_tip=0.65)
plot_with_padding(dna_root,  "DNA (scientific names; rerooted)",                             base_depth, cex_tip=0.65)
dev.off()

# ---- 6) Side-by-side export (topology only) ----
png(file.path(out_dir, "RF00740_three_trees_rerooted_noscale.png"),
    width=5200, height=1700, res=300)
par(mfrow=c(1,3), mar=c(2,1.5,3.2,1.5), oma=c(0,0,0,0), xpd=FALSE)
plot_topology_only(rfam_root, "RFAM seed (no branch lengths; outgroup = Eptesicus fuscus)")
plot_topology_only(rna_root,  "RNA (no branch lengths)")
plot_topology_only(dna_root,  "DNA (no branch lengths)")
dev.off()
