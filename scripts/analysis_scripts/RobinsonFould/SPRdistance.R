install.packages(c("ape","TreeDist","TBRDist"))  # run once
library(ape)
library(TreeDist)   # SPRDist() ≈ upper bound (heuristic)
library(TBRDist)    # USPRDist() = exact uSPR for unrooted trees

# ---- helpers ----
read_multi <- function(path) {
  tr <- read.tree(path)
  if (inherits(tr, "phylo")) list(tr) else unclass(tr)  # multiPhylo -> plain list
}

common_tips <- function(tree_list) {
  Reduce(intersect, lapply(tree_list, function(t) t$tip.label))
}

prep_trees <- function(trees, tips) {
  lapply(trees, function(t) {
    t <- keep.tip(t, tips)
    t <- multi2di(t)             # resolve polytomies consistently
    unroot(t)                    # ensure unrooted
  })
}

# Core: compute pairwise uSPR distances between two lists of trees
pairwise_uspr <- function(A, B) {
  # Exact uSPR for all pairs (vector of length length(A)*length(B))
  USPRDist(A, B, allPairs = TRUE, checks = TRUE)
}

pairwise_spr_upper <- function(A, B) {
  # Fast heuristic upper bound on SPR distance
  SPRDist(A, B)  # accepts lists; returns vector / matrix; we coerce below
}

# ---- main function for one family ----
spr_for_family <- function(fileA, fileB, labelA = "setA", labelB = "setB") {
  A0 <- read_multi(fileA)
  B0 <- read_multi(fileB)
  
  # Intersect taxa across all trees in both files
  tipsA <- common_tips(A0); tipsB <- common_tips(B0)
  tips  <- intersect(tipsA, tipsB)
  if (length(tips) < 4L) stop("Fewer than 4 shared taxa after intersection.")
  
  A <- prep_trees(A0, tips)
  B <- prep_trees(B0, tips)
  nA <- length(A); nB <- length(B); n <- length(tips)
  
  # Try exact uSPR; if it errors, compute an upper bound instead
  exact_vec <- tryCatch(pairwise_uspr(A, B), error = function(e) e)
  if (inherits(exact_vec, "error")) {
    message("Exact uSPR failed / too slow; returning an SPR upper bound instead.")
    up <- pairwise_spr_upper(A, B)
    # Coerce to vector if matrix was returned
    up_vec <- as.numeric(up)
    out <- data.frame(
      i = rep(seq_len(nA), each = nB),
      j = rep(seq_len(nB), times = nA),
      n_tips = n,
      uSPR = NA_real_,
      uSPR_norm = NA_real_,
      SPR_upper = up_vec,
      SPR_upper_norm = up_vec / (n - 3),
      setA = labelA,
      setB = labelB
    )
    return(out)
  } else {
    uspr_vec <- as.numeric(exact_vec)
    out <- data.frame(
      i = rep(seq_len(nA), each = nB),
      j = rep(seq_len(nB), times = nA),
      n_tips = n,
      uSPR = uspr_vec,
      uSPR_norm = uspr_vec / (n - 3),
      SPR_upper = NA_real_,
      SPR_upper_norm = NA_real_,
      setA = labelA,
      setB = labelB
    )
    return(out)
  }
}

# ---- example usage ----
# Replace these with your two files for ONE RNA family:
#   e.g. a "DNA-model" file vs a "RNA-structure-model" file
res <- spr_for_family("/Users/u7875558/RNAPhylo/fullAlignment_S6A/outputs/Robinson_Foulds_iqtree3/RF04235/RF04235.raxml", 
                      "/Users/u7875558/RNAPhylo/fullAlignment_S6A/outputs/Robinson_Foulds_iqtree3/RF04235/RF04235.raxmlPi",
                      labelA = "DNA", labelB = "RNA")
print(head(res))

summary(res$uSPR)              # min/median/max number of moves
summary(res$uSPR_norm)         # same, normalized
best <- res[which.min(res$uSPR), ]
best                           # the closest DNA-vs-RNA replicate pair
aggregate(uSPR ~ i, res, min)  # best RNA match for each DNA tree
aggregate(uSPR ~ j, res, min)  # best DNA match for each RNA tree
