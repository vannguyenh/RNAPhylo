#install.packages(c("ape","phytools"))  # if needed
library(ape)
library(phytools)

# helper: standardise labels to the accession
# full name: _AAOX01000026.1/47929-47995_Bacillus[313627].21
# output: AAOX01000026.1/47929-47995
canon <- function(x) {
  #x <- gsub("\\[[^]]*\\]", "", x)                 # remove [ ... ] comments
  x <- sub("^_*[0-9]+\\.?[0-9]*_", "", x)         # drop leading "64.8_"-style prefix
  x <- sub("^_*", "", x)                          # drop leading underscores
  sub("^([^_]+).*", "\\1", x)                     # keep up to first '_'
}

#####RF00164
# --- paste your trees as strings (or read from files) ---
rfam_raw <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/RF00164.rfam.newick')
rna_raw  <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/RAxML_bipartitions.RF00164_RNA')
dna_raw <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/RAxML_bipartitions.RF00164_DNA')

# copy trees and clean tip labels
t_rfam <- rfam_raw; t_rfam$tip.label <- canon(t_rfam$tip.label)
t_rna  <- rna_raw;  t_rna$tip.label  <- canon(t_rna$tip.label)
t_dna  <- dna_raw;  t_dna$tip.label  <- canon(t_dna$tip.label)

# (optional) write out clean files you can open in Dendroscope
write.tree(t_rfam, file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/RF00164_rfam_clean.newick")          # simple Newick
write.tree(t_rna,  file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/RF00164_rna_clean.newick")
write.tree(t_dna, file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/RF00164_dna_clean.newick")

### RF

# --- paste your trees as strings (or read from files) ---
rfam_raw <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF01690/RF01690.rfam.newick')
rna_raw  <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF01690/RAxML_bipartitions.RF01690_RNA')
dna_raw <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF01690/RAxML_bipartitions.RF01690_DNA')



# copy trees and clean tip labels
t_rfam <- rfam_raw; t_rfam$tip.label <- canon(t_rfam$tip.label)
t_rna  <- rna_raw;  t_rna$tip.label  <- canon(t_rna$tip.label)
t_dna  <- dna_raw;  t_dna$tip.label  <- canon(t_dna$tip.label)

# keep the common taxa and drop the rest
common <- Reduce(intersect, list(t_rfam$tip.label, t_rna$tip.label, t_dna$tip.label))
t_rfam <- drop.tip(t_rfam, setdiff(t_rfam$tip.label, common))
t_rna  <- drop.tip(t_rna,  setdiff(t_rna$tip.label,  common))
t_dna  <- drop.tip(t_dna,  setdiff(t_dna$tip.label,  common))

# sanity checks
stopifnot(length(t_rfam$tip.label) == length(t_rna$tip.label))
stopifnot(setequal(t_rfam$tip.label, t_rna$tip.label))

stopifnot(length(t_rfam$tip.label) == length(t_dna$tip.label))
stopifnot(setequal(t_rfam$tip.label, t_dna$tip.label))

# (optional) write out clean files you can open in Dendroscope
write.tree(t_rfam, file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF01690/RF01690_rfam_clean.newick")          # simple Newick
write.tree(t_rna,  file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF01690/RF01690_rna_clean.newick")
write.tree(t_dna, file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF01690/RF01690_dna_clean.newick")

set.seed(1)
co <- cophylo(t_rfam, t_rna, rotate = TRUE)  # tries to minimise crossings
plot(co, link.type = "curved", link.lwd = 1, fisze=0.3)


#### 24.09.25
# edit the trees
library(ape)
library(phytools)

# --- read
t_rfam <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF01690/RF01690.rfam.newick')
t_rna  <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF01690/RAxML_bipartitions.RF01690_RNA')
t_dna  <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF01690/RAxML_bipartitions.RF01690_DNA')

# --- keep labels unique: only strip [..] comments + leading underscores
canon <- function(x){
  x <- sub("^_*[0-9]+\\.?[0-9]*_", "", x)         # drop leading "64.8_"-style prefix
  x <- sub("^_*", "", x)                          # drop leading underscores
  sub("^([^_]+).*", "\\1", x)                     # keep up to first '_'
}
t_rfam$tip.label <- canon(t_rfam$tip.label)
t_rna$tip.label  <- canon(t_rna$tip.label)
t_dna$tip.label  <- canon(t_dna$tip.label)

# --- helper: build a one-to-one assoc even if there are duplicate names
make_assoc <- function(labsL, labsR){
  common_names <- intersect(labsL, labsR)
  # index duplicates deterministically: A, A -> A#1, A#2
  make_unique <- function(v){
    k <- ave(seq_along(v), v, FUN = seq_along)
    ifelse(k > 1, paste0(v, "#", k), v)
  }
  uL <- make_unique(labsL); uR <- make_unique(labsR)
  # pair the first min(nL, nR) occurrences for each name
  pairs <- do.call(rbind, lapply(common_names, function(nm){
    iL <- which(labsL == nm); iR <- which(labsR == nm)
    k <- min(length(iL), length(iR))
    if (k) cbind(uL[iL[seq_len(k)]], uR[iR[seq_len(k)]]) else NULL
  }))
  list(uL = uL, uR = uR, assoc = pairs,
       keepL = unique(pairs[,1]), keepR = unique(pairs[,2]))
}

# --- example tanglegram: RFAM vs RNA (keep extras)
pa <- make_assoc(t_rfam$tip.label, t_rna$tip.label)
L  <- t_rfam; L$tip.label <- pa$uL
R  <- t_rna;  R$tip.label <- pa$uR

tipcols_left  <- ifelse(L$tip.label %in% pa$keepL, "black", "grey70")
tipcols_right <- ifelse(R$tip.label %in% pa$keepR, "black", "grey70")

set.seed(1)
co <- cophylo(L, R, assoc = pa$assoc, rotate = TRUE)  # minimises crossings using only linked pairs

# note: it's fsize, not "fisze"
plot(co, link.type = "straight", link.lwd = 0.7, fsize = 0.55,
     tip.color = list(tipcols_left, tipcols_right))


fmt_support <- function(x, thr = NULL){           # returns character vector for plotting
  y <- suppressWarnings(as.numeric(x))
  y <- ifelse(is.na(y), NA, ifelse(y <= 1, 100*y, y))  # 0–1 → %
  if (!is.null(thr)) y <- ifelse(is.na(y) | y < thr, NA, y)
  ifelse(is.na(y), "", sprintf("%.0f", y))
}

# after plot(co, ...):
nodelabels.cophylo(co, which = "left",
                   text = fmt_support(co$trees[[1]]$node.label, thr = 50),
                   frame = "n", adj = c(0.5, -0.2), cex = 0.5)

nodelabels.cophylo(co, which = "right",
                   text = fmt_support(co$trees[[2]]$node.label, thr = 50),
                   frame = "n", adj = c(0.5, -0.2), cex = 0.5)


### reroot

# by outgroup (can be one or several tips)
rfam_raw  <- root(rfam_raw, outgroup = c("AY188187.1/197-239"), resolve.root = TRUE)
dna_raw <- root(dna_raw, outgroup = c("AY188187.1/197-239"), resolve.root = TRUE)
rna_raw <- root(rna_raw, outgroup = c("AY188187.1/197-239"), resolve.root = TRUE)

# or midpoint rooting (no outgroup required)
#t_rna  <- midpoint.root(t_rna)   # phytools
#t_dna  <- midpoint.root(t_dna)
#t_rfam <- midpoint.root(t_rfam)

plot(co, link.type = "curved", gap = 6, link.lwd = 0.6, fsize = 0.55,
     tip.color = list(tipcols_left, tipcols_right))


# ------------- A) View the single trees (side by side) ----------------
fmt_support <- function(x, thr = NULL){
  y <- suppressWarnings(as.numeric(x))
  y <- ifelse(is.na(y), NA, ifelse(y <= 1, 100*y, y))  # 0–1 -> %
  if (!is.null(thr)) y <- ifelse(is.na(y) | y < thr, NA, y)
  ifelse(is.na(y), "", sprintf("%.0f", y))
}

# support helpers
to_pct <- function(x){
  y <- suppressWarnings(as.numeric(x))
  ifelse(is.na(y), NA, ifelse(y <= 1, 100*y, y))  # 0–1 -> %
}
scale_cex <- function(supp, cmin=0.45, cmax=1.2){
  y <- to_pct(supp)
  if (all(is.na(y))) return(rep((cmin+cmax)/2, length(y)))
  rng <- range(y, na.rm = TRUE)
  out <- cmin + (y - rng[1]) / diff(rng) * (cmax - cmin)
  out[!is.finite(out) | out <= 0] <- 1e-6  # never 0 (avoids cex errors)
  out
}

# draw support "bubbles" + numbers offset from the split
add_support <- function(tr, thr = 50, show_numbers = TRUE, bubbles = TRUE){
  lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  nn <- tr$Nnode
  nodes <- (Ntip(tr) + 1):(Ntip(tr) + nn)
  x <- lp$xx[nodes]; y <- lp$yy[nodes]
  sup <- to_pct(tr$node.label)
  
  if (bubbles) {
    points(x, y, pch = 21, bg = "white",
           cex = scale_cex(sup), col = "black", lwd = 0.6)
  }
  if (show_numbers) {
    lab <- ifelse(is.na(sup) | sup < thr, "", sprintf("%.0f", sup))
    dx <- max(node.depth.edgelength(tr)) * 0.02  # small right shift
    text(x + dx, y, labels = lab, cex = 0.45, xpd = NA)
  }
}

# put support numbers slightly to the right of each node (not at the split)
add_support_offset <- function(tr, thr = 50, cex = 0.45, offset_frac = 0.02, font=3, col=black){
  lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  nn <- tr$Nnode
  nodes <- (Ntip(tr) + 1):(Ntip(tr) + nn)
  x <- lp$xx[nodes]; y <- lp$yy[nodes]
  H <- max(node.depth.edgelength(tr))
  dx <- H * offset_frac
  text(x + dx, y, labels = fmt_support(tr$node.label, thr = thr), cex = cex, xpd = NA)
}

op <- par(mfrow = c(1,3), mar = c(1,1,2,1)); on.exit(par(op), add = TRUE)

plot(ladderize(t_rfam), no.margin = TRUE, main = "FastTree (rerooted)", cex = 0.55)
#nodelabels(fmt_support(t_rfam$node.label, thr = 50), frame = "n", cex = 0.45)
add_support_offset(t_rfam, thr = 50, cex = 0.45, font=2)  # <- replaces nodelabels(...)

plot(ladderize(t_dna),  no.margin = TRUE, main = "DNA (rerooted)",  cex = 0.55)
#nodelabels(fmt_support(t_dna$node.label,  thr = 50), frame = "n", cex = 0.45)
add_support_offset(t_dna,  thr = 50, cex = 0.45, font =2)  # <- replaces nodelabels(...)

plot(ladderize(t_rna),  no.margin = TRUE, main = "RNA (S6A - rerooted)",  cex = 0.55)
#nodelabels(fmt_support(t_rna$node.label,  thr = 50), frame = "n", cex = 0.45)
add_support_offset(t_rna,  thr = 50, cex = 0.45, font =2)  # <- replaces nodelabels(...)

#op <- par(mfrow = c(1,3), mar = c(1,1,1,1), xpd = NA); on.exit(par(op), add = TRUE)

#plot(ladderize(t_rfam), no.margin = FALSE, main = "FastTree (rerooted)", cex = 0.55, align.tip.label = TRUE)
#add_support(t_rfam, thr = 50)   # <-- replaces nodelabels(...)

#plot(ladderize(t_dna),  no.margin = FALSE, main = "DNA (rerooted)",    cex = 0.55, align.tip.label = TRUE)
#add_support(t_dna,  thr = 50)   # <-- replaces nodelabels(...)

#plot(ladderize(t_rna),  no.margin = FALSE, main = "RNA (rerooted)",    cex = 0.55, align.tip.label = TRUE)
#add_support(t_rna,  thr = 50)   # <-- replaces nodelabels(...)

#### 07.10.25
# 1) OPEN a PNG device (pick your path/size)
ragg::agg_png(
  "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/RF00164_three_trees.png",
  width = 3000, height = 1600, res = 300, background = "white"
)

# 2) your plotting code (unchanged except nodelabels -> add_support_offset)
op <- par(mfrow = c(1,3), mar = c(1,1,2,1)); on.exit(par(op), add = TRUE)

plot(ladderize(t_rfam), no.margin = TRUE, main = "FastTree (rerooted)", cex = 0.55)
add_support_offset(t_rfam, thr = 50, cex = 0.45, font = 3)

plot(ladderize(t_dna),  no.margin = TRUE, main = "DNA (rerooted)",  cex = 0.55)
add_support_offset(t_dna,  thr = 50, cex = 0.45, font = 3)

plot(ladderize(t_rna),  no.margin = TRUE, main = "RNA (S6A - rerooted)",  cex = 0.55)
add_support_offset(t_rna,  thr = 50, cex = 0.45, font = 3)

# 3) CLOSE the device (this writes the PNG to disk)
dev.off()


### RF00011
# --- paste your trees as strings (or read from files) ---
rfam_raw <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00011/RF00011.rfam.newick')
rna_raw  <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00011/RAxML_bipartitions.RF00011_RNA')
dna_raw <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00011/RAxML_bipartitions.RF00011_DNA')

# copy trees and clean tip labels
t_rfam <- rfam_raw; t_rfam$tip.label <- canon(t_rfam$tip.label)
t_rna  <- rna_raw;  t_rna$tip.label  <- canon(t_rna$tip.label)
t_dna  <- dna_raw;  t_dna$tip.label  <- canon(t_dna$tip.label)

# (optional) write out clean files you can open in Dendroscope
write.tree(t_rfam, file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00011/RF00011_rfam_clean.newick")          # simple Newick
write.tree(t_rna,  file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00011/RF00011_rna_clean.newick")
write.tree(t_dna, file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00011/RF00011_dna_clean.newick")

# keep the common taxa and drop the rest
common <- Reduce(intersect, list(t_rfam$tip.label, t_rna$tip.label, t_dna$tip.label))
t_rfam <- drop.tip(t_rfam, setdiff(t_rfam$tip.label, common))
t_rna  <- drop.tip(t_rna,  setdiff(t_rna$tip.label,  common))
t_dna  <- drop.tip(t_dna,  setdiff(t_dna$tip.label,  common))

# sanity checks
stopifnot(length(t_rfam$tip.label) == length(t_rna$tip.label))
stopifnot(setequal(t_rfam$tip.label, t_rna$tip.label))

stopifnot(length(t_rfam$tip.label) == length(t_dna$tip.label))
stopifnot(setequal(t_rfam$tip.label, t_dna$tip.label))

# by outgroup (can be one or several tips)
t_rna  <- root(t_rna, outgroup = c("AACY023573249.1/971-644"), resolve.root = TRUE)
t_dna  <- root(t_dna, outgroup = c("AACY023573249.1/971-644"), resolve.root = TRUE)
t_rfam <- root(t_rfam, outgroup = c("AACY023573249.1/971-644"), resolve.root = TRUE)

# 1) OPEN a PNG device (pick your path/size)
ragg::agg_png(
  "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/RF00164_three_trees.png",
  width = 3200, height = 2000, res = 300, background = "white"
)

# 2) your plotting code (unchanged except nodelabels -> add_support_offset)
op <- par(mfrow = c(1,3), mar = c(1,1,2,1)); on.exit(par(op), add = TRUE)

plot(ladderize(t_rfam), no.margin = TRUE, main = "FastTree (rerooted)", cex = 0.55)
add_support_offset(t_rfam, thr = 50, cex = 0.45, font = 3)

plot(ladderize(t_dna),  no.margin = TRUE, main = "DNA (rerooted)",  cex = 0.55)
add_support_offset(t_dna,  thr = 50, cex = 0.45, font = 3)

plot(ladderize(t_rna),  no.margin = TRUE, main = "RNA (S6A - rerooted)",  cex = 0.55)
add_support_offset(t_rna,  thr = 50, cex = 0.45, font = 3)

# 3) CLOSE the device (this writes the PNG to disk)
dev.off()


## RF03165
# --- paste your trees as strings (or read from files) ---
rfam_raw <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF03165/RF03165.rfam.newick')
rna_raw  <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF03165/RAxML_bipartitions.RF03165_RNA')
dna_raw <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF03165/RAxML_bipartitions.RF03165_DNA')

# copy trees and clean tip labels
t_rfam <- rfam_raw; t_rfam$tip.label <- canon(t_rfam$tip.label)
t_rna  <- rna_raw;  t_rna$tip.label  <- canon(t_rna$tip.label)
t_dna  <- dna_raw;  t_dna$tip.label  <- canon(t_dna$tip.label)

# (optional) write out clean files you can open in Dendroscope
write.tree(t_rfam, file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF03165/RF03165_rfam_clean.newick")          # simple Newick
write.tree(t_rna,  file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF03165/RF03165_rna_clean.newick")
write.tree(t_dna, file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF03165/RF03165_dna_clean.newick")

# keep the common taxa and drop the rest
common <- Reduce(intersect, list(t_rfam$tip.label, t_rna$tip.label, t_dna$tip.label))
t_rfam <- drop.tip(t_rfam, setdiff(t_rfam$tip.label, common))
t_rna  <- drop.tip(t_rna,  setdiff(t_rna$tip.label,  common))
t_dna  <- drop.tip(t_dna,  setdiff(t_dna$tip.label,  common))

# sanity checks
stopifnot(length(t_rfam$tip.label) == length(t_rna$tip.label))
stopifnot(setequal(t_rfam$tip.label, t_rna$tip.label))

stopifnot(length(t_rfam$tip.label) == length(t_dna$tip.label))
stopifnot(setequal(t_rfam$tip.label, t_dna$tip.label))

# by outgroup (can be one or several tips)
#t_rna  <- root(t_rna, outgroup = c("URS0000D6D221"), resolve.root = TRUE)
#t_dna  <- root(t_dna, outgroup = c("URS0000D6D221"), resolve.root = TRUE)
#t_rfam <- root(t_rfam, outgroup = c("URS0000D6D221"), resolve.root = TRUE)

# 1) OPEN a PNG device (pick your path/size)
ragg::agg_png(
  "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF03165/RF03165_three_trees.png",
  width = 3200, height = 2000, res = 300, background = "white"
)

# 2) your plotting code (unchanged except nodelabels -> add_support_offset)
op <- par(mfrow = c(1,3), mar = c(1,1,2,1)); on.exit(par(op), add = TRUE)

plot(ladderize(t_rfam), no.margin = TRUE, main = "FastTree (rerooted)", cex = 0.55)
add_support_offset(t_rfam, thr = 50, cex = 0.45, font = 3)

plot(ladderize(t_dna),  no.margin = TRUE, main = "DNA (rerooted)",  cex = 0.55)
add_support_offset(t_dna,  thr = 50, cex = 0.45, font = 3)

plot(ladderize(t_rna),  no.margin = TRUE, main = "RNA (S6A - rerooted)",  cex = 0.55)
add_support_offset(t_rna,  thr = 50, cex = 0.45, font = 3)

# 3) CLOSE the device (this writes the PNG to disk)
dev.off()

# (optional) write out clean files you can open in Dendroscope
write.tree(t_rfam, file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF03165/RF03165_rfam_clean.newick")          # simple Newick
write.tree(t_rna,  file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF03165/RF03165_rna_clean.newick")
write.tree(t_dna, file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF03165/RF03165_dna_clean.newick")

############ RF00740
# --- paste your trees as strings (or read from files) ---
rfam_raw <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00740/RF00740.rfam.newick')
rna_raw  <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00740/RAxML_bipartitions.RF00740_RNA')
dna_raw <- read.tree('/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00740/RAxML_bipartitions.RF00740_DNA')

# copy trees and clean tip labels
t_rfam <- rfam_raw; t_rfam$tip.label <- canon(t_rfam$tip.label)
t_rna  <- rna_raw;  t_rna$tip.label  <- canon(t_rna$tip.label)
t_dna  <- dna_raw;  t_dna$tip.label  <- canon(t_dna$tip.label)

# (optional) write out clean files you can open in Dendroscope
write.tree(t_rfam, file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/RF00164_rfam_clean.newick")          # simple Newick
write.tree(t_rna,  file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/RF00164_rna_clean.newick")
write.tree(t_dna, file = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison/RF00164/RF00164_dna_clean.newick")
