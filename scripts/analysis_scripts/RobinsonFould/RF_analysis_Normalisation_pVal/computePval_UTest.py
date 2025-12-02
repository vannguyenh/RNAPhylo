#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mann–Whitney U-test on normalised RF distances for ONE model.

Compares, per RNA:
  - DNA vs DNA     (file: DIR_RF/DNA_extra/<RNA>.rfdist in the DNA_extra folder)
  - DNA vs RNA     (file: DIR_RF/<MODEL>/<RNA>.rfdist in the S* folder)

Outputs:
  1) {out_dir}/{tag}_long.csv   with columns: Model, RNA, n_DNA, n_RNA, U, pvalue
  2) {out_dir}/{tag}_wide.csv   wide RNA × Model with pvalues
"""

import os
from os.path import join, isdir
import numpy as np
import pandas as pd
from Bio import Phylo

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


# -------- PATHS (edit if needed) ---------------------------------------------
DIR_WORKING = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels"
DIR_OUTPUTS = join(DIR_WORKING, "outputs")
DIR_RF      = join(DIR_OUTPUTS, "251121_Robinson_Foulds")          # <MODEL>/<RNA>/<RNA>.rfdist
DIR_DNA     = join(DIR_OUTPUTS, "inferred_trees", "DNA")    # DNA trees for counting taxa
RF_SUFFIX   = ".rfdist"                       # used for both DNA_extra and S* models

ALPHA = 0.05

# ── I/O helpers ────────────────────────────────────────────────────────────────
def read_rfdist_matrix(path: str) -> np.ndarray:
    """Read IQ-TREE .rfdist square matrix to 2D array; empty array if missing/empty."""
    if not os.path.exists(path):
        return np.array([])
    with open(path, "r") as f:
        lines = f.readlines()
    if len(lines) <= 1:
        return np.array([])
    rows = []
    for line in lines[1:]:
        parts = line.strip().split()
        if len(parts) <= 1:
            continue
        rows.append(list(map(float, parts[1:])))  # skip row label
    return np.array(rows, dtype=float) if rows else np.array([])

def find_one_tree_for_taxa(dir_with_trees: str) -> str | None:
    """
    Return a path to any tree file for this RNA to count taxa from.
    """
    if isdir(dir_with_trees):
        for fn in sorted(os.listdir(dir_with_trees)):
            if fn.startswith("RAxML_bestTree"):
                p = join(dir_with_trees, fn)
                if os.path.exists(p):
                    return p
    return None

def count_taxa(tree_path: str) -> int | None:
    """Count terminals in the first tree from a Newick file."""
    if not tree_path or not os.path.exists(tree_path): return None
    try:
        with open(tree_path, "r") as fh:
            tree = next(Phylo.parse(fh, "newick"))
        return len(tree.get_terminals())
    except Exception:
        return None

def normalize_rf_matrix_by_n(mat: np.ndarray, n_taxa: int) -> np.ndarray:
    """
    Standard RF normalization to [0,1]:
      nRF = RF / (2*(n-3))  for unrooted binary trees.
    If max(mat) ≤ 1.0, assume already normalized; return as-is.
    """
    if mat.size == 0 or not n_taxa or n_taxa < 4: return np.array([])
    denom = 2 * (n_taxa - 3)
    return mat / denom

def model_to_category(model):
    return f"DNA vs RNA ({model})"

# ── Core ───────────────────────────────────────────────────────────────────────
# ── Core per-RNA computation ───────────────────────────────────────────────────
def process_one_rna(rna_dir: str, rna: str, model:str) -> dict | None:
    """
    Read, normalize, summarize, and U-test one RNA.
    Returns a dict of results, or None if cannot process.
    """
    # RF matrices: DNA_extra and S* model
    dna_fp   = join(DIR_RF, "DNA_extra", rna, f"{rna}{RF_SUFFIX}") # DNA vs DNA_extra
    rna_fp   = join(DIR_RF, model,       rna, f"{rna}{RF_SUFFIX}") # DNA vs RNA in RNA model
    if not (os.path.exists(dna_fp) and os.path.exists(rna_fp)):
        return None

    # taxon counts from any DNA tree and the model’s DNA trees
    dna_tree = find_one_tree_for_taxa(join(DIR_DNA, rna))
    mod_tree = find_one_tree_for_taxa(join(DIR_OUTPUTS, "inferred_trees", model, rna))
    n_dna = count_taxa(dna_tree)
    n_mod = count_taxa(mod_tree)
    if not n_dna or not n_mod or n_dna < 4 or n_mod < 4: return None
    if n_dna != n_mod:
        # require same n for comparable nRF
        return None
    n = n_dna

    dna_raw = read_rfdist_matrix(dna_fp)
    rna_raw = read_rfdist_matrix(rna_fp)
    if dna_raw.size == 0 or rna_raw.size == 0: return None

    dna_nrf = normalize_rf_matrix_by_n(dna_raw, n)
    rna_nrf = normalize_rf_matrix_by_n(rna_raw, n)
    if dna_nrf.size == 0 or rna_nrf.size == 0: return None

    dna_vals = dna_nrf.flatten().astype(float)   # 100
    rna_vals = rna_nrf.flatten().astype(float)   # 100
    if dna_vals.size == 0 or rna_vals.size == 0: return None

    U, p = mannwhitneyu(dna_vals, rna_vals, alternative="two-sided")

    return {
        "Model": model,
        "RNA": rna,
        "n_taxa": int(n),
        "n_DNA":  int(dna_vals.size),
        "n_RNA":  int(rna_vals.size),
        "DNA_median_nRF": float(np.median(dna_vals)),
        "RNA_median_nRF": float(np.median(rna_vals)),
        "U": float(U),
        "pvalue": float(p),
    }


def main():
    if not isdir(DIR_RF):
        raise SystemExit(f"DIR_RF not found: {DIR_RF}")

    models = [m for m in sorted(os.listdir(DIR_RF)) if isdir(join(DIR_RF, m)) and m.startswith("S")]
    results, med_rows = [], []

    for model in models:
        model_dir = join(DIR_RF, model)
        rnas = [d for d in sorted(os.listdir(model_dir)) if isdir(join(model_dir, d))]
        model_rows = []
        for rna in rnas:
            rec = process_one_rna(model, rna)
            if rec: model_rows.append(rec)
        if not model_rows:
            continue

        # per-model FDR (BH)
        pvals = np.array([r["pvalue"] for r in model_rows], dtype=float)
        rej, qvals, _, _ = multipletests(pvals, alpha=ALPHA, method="fdr_bh")
        for r, q, ok in zip(model_rows, qvals, rej):
            r["p_fdr_model"] = float(q)
            r["sig_fdr_model"] = bool(ok)

        results.extend(model_rows)

        cat = model_to_category(model)
        for r in model_rows:
            med_rows.append({"Model": model, "RNA": r["RNA"], "Category": "DNA vs DNA", "Median nRF": r["DNA_median_nRF"]})
            med_rows.append({"Model": model, "RNA": r["RNA"], "Category": cat,           "Median nRF": r["RNA_median_nRF"]})

    if not results:
        print("No usable data found across models.")
        return None, None, None

    df_long = pd.DataFrame(results).sort_values(["Model","RNA"]).reset_index(drop=True)

    # global FDR (across all model×RNA tests), optional:
    rej_g, q_g, _, _ = multipletests(df_long["pvalue"].values, alpha=ALPHA, method="fdr_bh")
    df_long["p_fdr_global"] = q_g
    df_long["sig_fdr_global"] = rej_g

    df_wide = (
        df_long.pivot_table(index="RNA", columns="Model", values="pvalue", aggfunc="first")
               .sort_index()
    )
    df_medians = pd.DataFrame(med_rows)

    # Write outputs
    out1 = join(DIR_RF, "SEED_Utest_all_models_long.csv")
    out2 = join(DIR_RF, "Utest_all_models_wide.csv")
    out3 = join(DIR_RF, "Median_nRF_all_models.csv")

    df_long.to_csv(out1, index=False)
    df_wide.to_csv(out2, index=False)
    df_medians.to_csv(out3, index=False)

    return df_long, df_wide, df_medians

# ── CLI ────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    df_long, df_wide, df_medians = main()
