#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mann–Whitney U-test on normalised RF distances for ONE model.

Compares, per RNA:
  - DNA vs DNA     (file: <RNA>.raxml.rfdist;     use upper triangle)
  - DNA vs RNA     (file: <RNA>.raxml.raxmlPi.rfdist; use all cells)

Outputs:
  1) {out_dir}/{tag}_long.csv   with columns: Model, RNA, n_DNA, n_RNA, U, pvalue
  2) {out_dir}/{tag}_wide.csv   wide RNA × Model with pvalues
"""

import os
from os.path import join, isdir
import numpy as np
import pandas as pd
from sklearn.preprocessing import normalize
from scipy.stats import mannwhitneyu

# ── PATHS (edit to your machine) ───────────────────────────────────────────────
DIR_WORKING = "/Users/u7875558/RNAPhylo/fullAlignment_S6A"
DIR_OUTPUTS = join(DIR_WORKING, "outputs")
DIR_RF      = join(DIR_OUTPUTS, "Robinson_Foulds_iqtree3")   # contains <RNA>/ files
MODEL       = "S6A"  # label to write in outputs

# Expected file suffixes
DNA_RFDIST_SUFFIX       = ".raxml.rfdist"
DNA_VS_RNA_IPSEU_SUFFIX = ".raxml.raxmlPi.rfdist"

# ── I/O helpers ────────────────────────────────────────────────────────────────
def read_rfdist_matrix(path: str) -> np.ndarray:
    """Read IQ-TREE .rfdist square matrix to 2D array; empty array if missing/empty."""
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

def norm_upper_triangle(mat: np.ndarray) -> np.ndarray:
    """L2-normalize and return upper triangle (k=1) flattened."""
    if mat.size == 0:
        return np.array([])
    if mat.ndim == 1:
        mat = mat.reshape(1, -1)
    n = mat.shape[0]
    norm = normalize(mat, norm="l2")
    iu = np.triu_indices(n, k=1)
    return norm[iu].astype(float)

def norm_full_flat(mat: np.ndarray) -> np.ndarray:
    """L2-normalize and return all entries flattened."""
    if mat.size == 0:
        return np.array([])
    if mat.ndim == 1:
        mat = mat.reshape(1, -1)
    norm = normalize(mat, norm="l2")
    return norm.flatten().astype(float)

# ── Core ───────────────────────────────────────────────────────────────────────
def compute_utest_for_model(model: str, rf_root: str) -> pd.DataFrame:
    """
    Iterate RNA folders under rf_root and compute U-test p-values:
      DNA vs DNA  vs  DNA vs RNA (ignore pseudoknots)
    """
    if not isdir(rf_root):
        return pd.DataFrame(columns=["Model","RNA","n_DNA","n_RNA","U","pvalue"])

    records = []
    for rna in sorted(os.listdir(rf_root)):
        rna_dir = join(rf_root, rna)
        if not isdir(rna_dir):
            continue

        dna_fp   = join(rna_dir, f"{rna}{DNA_RFDIST_SUFFIX}")
        ipseu_fp = join(rna_dir, f"{rna}{DNA_VS_RNA_IPSEU_SUFFIX}")
        if not (os.path.exists(dna_fp) and os.path.exists(ipseu_fp)):
            continue

        dna_vals   = norm_upper_triangle(read_rfdist_matrix(dna_fp))
        ipseu_vals = norm_full_flat(read_rfdist_matrix(ipseu_fp))

        if dna_vals.size == 0 or ipseu_vals.size == 0:
            continue

        res = mannwhitneyu(dna_vals, ipseu_vals, alternative="two-sided")
        U = getattr(res, "statistic", res[0])
        p = getattr(res, "pvalue",    res[1])

        records.append({
            "Model":  model,
            "RNA":    rna,
            "n_DNA":  int(dna_vals.size),
            "n_RNA":  int(ipseu_vals.size),
            "U":      float(U),
            "pvalue": float(p),
        })

    return pd.DataFrame.from_records(records)

def run_single_model(model: str = MODEL,
                     rf_root: str = DIR_RF,
                     out_dir: str = DIR_RF,
                     tag: str = "Utest_pvalues_ignore_pseudoknots"):
    """Compute and write long & wide tables for a single model."""
    df_long = compute_utest_for_model(model, rf_root)
    if df_long.empty:
        print("No usable RNAs found — check paths and files.")
        return df_long, pd.DataFrame()

    df_wide = (
        df_long.pivot(index="RNA", columns="Model", values="pvalue")
               .sort_index()
    )

    os.makedirs(out_dir, exist_ok=True)
    long_fp = join(out_dir, f"{tag}_long.csv")
    wide_fp = join(out_dir, f"{tag}_wide.csv")
    df_long.to_csv(long_fp, index=False)
    df_wide.to_csv(wide_fp)

    print(f"Wrote long table: {long_fp}  (rows={len(df_long)})")
    print(f"Wrote wide table: {wide_fp}  (RNAs={df_wide.shape[0]})")
    return df_long, df_wide

# ── CLI ────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print(f"Running U-test for model: {MODEL}")
    run_single_model()
