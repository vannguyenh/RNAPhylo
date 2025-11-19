# ============================
# Relabel taxa in 3 trees (seed, DNA, RNA)
# - No rerooting, no plotting
# - If duplicate taxon names occur, append the sequence range [start-end]
#   (fallback: short accession tag) to make them unique.
# ============================

# install.packages(c("ape","stringr","rentrez"), repos="https://cloud.r-project.org")
library(ape)
library(stringr)
library(rentrez)

# ---------- CONFIG ----------
RNA <- 'RF00167'

# Directory layout (adjust if yours differs)
base_dir    <- file.path("/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison", RNA)
seed_file   <- file.path(base_dir, paste0(RNA, ".seed_tree"))

# For RAxML outputs (rapidBS layout from your previous messages)
raxml_dir   <- file.path(base_dir, "raxml_rapidBS")
dna_file    <- file.path(raxml_dir, paste0("RAxML_bipartitions.", RNA, ".rapidBS.DNA"))
rna_file    <- file.path(raxml_dir, paste0("RAxML_bipartitions.", RNA, ".rapidBS.S6A"))

# Optional: be polite to NCBI (and add your API key if you have one)
# Sys.setenv(ENTREZ_EMAIL = "you@anu.edu.au")
# rentrez::set_entrez_key("YOUR_NCBI_API_KEY")

# ---------- Sanity checks ----------
stopifnot("Seed tree not found" = file.exists(seed_file))
stopifnot("DNA tree not found"  = file.exists(dna_file))
stopifnot("RNA tree not found"  = file.exists(rna_file))

# ---------- Helpers ----------
# Grab all tip labels (tokens just before :<number>)
tip_labels_from_newick <- function(txt){
  m <- gregexpr("(?<=\\(|,)([^():;,]+?)(?=:\\-?\\d)", txt, perl=TRUE)
  trimws(regmatches(txt, m)[[1]])
}

# Clean odd numeric prefixes
canon_head <- function(x){
  x <- sub("^_*[0-9]+\\.?[0-9]*_", "", x)  # drop like "64.8_"
  sub("^_*", "", x)
}

# Accession token shared across trees
# - URS with taxid: keep "URS...._NNNN"
# - bare URS: keep core "URS...."
# - others: first token before "_" (e.g., AY304470.1), then drop coords after "/"
extract_accession <- function(label){
  head <- canon_head(label)
  if (grepl("^URS[0-9A-F]+_\\d+$", head, ignore.case = TRUE)) {
    sub("/.*$", "", head)
  } else if (grepl("^URS[0-9A-F]+$", head, ignore.case = TRUE)) {
    head
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

# Fallback species text from the label (tries to recover "Genus species ...")
guess_species <- function(lbl){
  # chunk after last '/', then drop trailing [taxid], underscores -> spaces
  s <- sub("^.*?/", "", lbl)
  s <- sub("\\[\\d+\\].*$", "", s)
  s <- gsub("_+", " ", s)
  s <- sub("\\.+$", "", s)
  trimws(s)
}

# Lookup scientific name from NCBI taxonomy (by taxid)
lookup_taxid_name <- function(tid){
  s <- tryCatch(rentrez::entrez_summary(db="taxonomy", id=tid), error=function(e) NULL)
  if (is.null(s)) NA_character_ else s$scientificname
}

# Quote labels for Newick
quote_newick   <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
unquote_newick <- function(x) sub("^'(.*)'$", "\\1", x)

# Extract a coordinate range from the original tip
# Supports:
#   ".../197-239_..."            -> 197-239
#   "...:13684-13586_..."        -> 13684-13586
#   ".../2143272-2143175[...]"   -> 2143272-2143175
extract_coords_range <- function(label){
  m1 <- regexpr("/([0-9]+-[0-9]+)", label, perl = TRUE)
  if (m1[1] != -1) return(substring(label, m1[1] + 1, m1[1] + attr(m1, "match.length") - 1))
  m2 <- regexpr(":([0-9]+-[0-9]+)", label, perl = TRUE)
  if (m2[1] != -1) return(substring(label, m2[1] + 1, m2[1] + attr(m2, "match.length") - 1))
  NA_character_
}

# Shorten accession for fallback tags (ABDQ01000007.1 -> ABDQ01000007)
short_acc <- function(acc) sub("\\..*$", "", acc)

# Build accession -> preferred label map from the RFAM seed tips
# We use: scientific name from taxid if available; otherwise parsed fallback.
# Also alias URSxxxx_core to URSxxxx_taxid so DNA/RNA map correctly.
build_label_map_from_seed <- function(seed_tips){
  accs   <- vapply(seed_tips, extract_accession, "", USE.NAMES=FALSE)
  taxids <- vapply(seed_tips, extract_taxid,     NA_integer_, USE.NAMES=FALSE)
  
  sci <- character(length(seed_tips))
  for (i in seq_along(seed_tips)){
    if (!is.na(taxids[i])) sci[i] <- lookup_taxid_name(taxids[i])
    if (!nzchar(sci[i]))   sci[i] <- guess_species(seed_tips[i])  # fallback
  }
  
  keep       <- !duplicated(accs)
  acc2name   <- setNames(sci[keep],    accs[keep])
  acc2taxid  <- setNames(taxids[keep], accs[keep])
  
  # alias URS core
  is_urs_full <- grepl("^URS[0-9A-F]+_\\d+$", names(acc2name))
  if (any(is_urs_full)){
    urs_core <- sub("^(URS[0-9A-F]+)_\\d+$", "\\1", names(acc2name)[is_urs_full])
    acc2name  <- c(acc2name,  setNames(unname(acc2name [is_urs_full]), urs_core))
    acc2taxid <- c(acc2taxid, setNames(unname(acc2taxid[is_urs_full]), urs_core))
  }
  
  list(acc2name=acc2name, acc2taxid=acc2taxid)
}

# Disambiguate duplicates by appending coordinate range (fallback: short accession)
uniquify_with_coords <- function(base_labels, tip_labels_raw){
  out <- base_labels
  dups <- duplicated(base_labels) | duplicated(base_labels, fromLast = TRUE)
  if (!any(dups)) return(out)
  
  coords <- vapply(tip_labels_raw, extract_coords_range, "", USE.NAMES = FALSE)
  acc    <- vapply(tip_labels_raw, extract_accession, "", USE.NAMES = FALSE)
  tag    <- ifelse(nzchar(coords), coords, short_acc(acc))
  
  out[dups] <- sprintf("%s [%s]", base_labels[dups], tag[dups])
  
  # Absolute uniqueness guard
  if (any(duplicated(out))) out <- make.unique(out, sep = " _")
  out
}

# Relabel a Newick using accession->name map; disambiguate duplicates via coords
relabel_newick <- function(newick_text, acc2name){
  tips <- tip_labels_from_newick(newick_text)
  if (!length(tips)) return(newick_text)
  
  accs <- vapply(tips, extract_accession, "", USE.NAMES=FALSE)
  
  # base label from mapping; fallback to parsed species from each tip
  base <- vapply(seq_along(tips), function(i){
    a <- accs[i]
    nm <- acc2name[[a]]
    if (!is.null(nm) && nzchar(nm)) return(nm)
    guess_species(tips[i])
  }, "")
  
  # Make duplicates unique with [start-end] (or short acc)
  unique_labels <- uniquify_with_coords(base, tips)
  
  # Splice back into the Newick
  quoted <- vapply(unique_labels, quote_newick, "")
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

# ---------- MAIN ----------
# Read original Newick text
seed_txt <- paste(readLines(seed_file, warn=FALSE), collapse="\n")
dna_txt  <- paste(readLines(dna_file,  warn=FALSE), collapse="\n")
rna_txt  <- paste(readLines(rna_file,  warn=FALSE), collapse="\n")

# Build label map from the seed tree (drives canonical names)
seed_tips <- tip_labels_from_newick(seed_txt)
lut <- build_label_map_from_seed(seed_tips)
acc2name <- lut$acc2name

# Relabel all three trees
seed_labeled_txt <- relabel_newick(seed_txt, acc2name)
dna_labeled_txt  <- relabel_newick(dna_txt,  acc2name)
rna_labeled_txt  <- relabel_newick(rna_txt,  acc2name)

# Write outputs next to inputs
seed_out <- sub("(\\.newick|\\.nwk|\\.tre|\\.seed_tree)?$", ".labeled.newick", seed_file, ignore.case=TRUE)
dna_out  <- paste0(dna_file, ".labeled.newick")
rna_out  <- paste0(rna_file, ".labeled.newick")

writeLines(seed_labeled_txt, seed_out)
writeLines(dna_labeled_txt,  dna_out)
writeLines(rna_labeled_txt,  rna_out)

message("Done.\n",
        "Seed → ", seed_out, "\n",
        "DNA  → ", dna_out,  "\n",
        "RNA  → ", rna_out,  "\n")
