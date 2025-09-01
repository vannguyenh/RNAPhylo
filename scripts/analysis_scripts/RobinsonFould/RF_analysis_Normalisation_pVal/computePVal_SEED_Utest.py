"""
Batch U-test on normalized RF distances (single model, many RNA families).

For each RNA directory under DIR_RF:
  • Read DNA↔DNA       matrix: <RNA>.raxml.rfdist         (10×10)
  • Read DNA↔RNA (Pi)  matrix: <RNA>.raxml.raxmlPi.rfdist (10×10)
  • Normalize each cell to nRF = RF / [2*(n-3)]  (auto-skips if already ≤1)
  • DNA↔DNA: use upper triangle (k=1) → 45 values
    DNA↔RNA: use all 100 values
  • Compute two-sided Mann–Whitney U-test
  • Record per-RNA medians for both groups
Outputs:
  1) Utest_long.csv        : Model, RNA, n_DNA, n_RNA, U, pvalue, p_bonf, flags
  2) Utest_wide.csv        : RNA × Model table of (raw) p-values
  3) Median_nRF_long.csv   : long table of per-RNA medians (for plotting)
Notes:
  • n = #taxa is taken from a tree file per RNA. We try, in order:
      <RNA>.raxml  →  <RNA>.raxmlPi  →  any DNA tree in DIR_DNA/<RNA>/RAxML_bestTree.*
"""

import os
from os.path import join, isdir
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.preprocessing import normalize as _ignore   # not used; kept to show change
from Bio import Phylo

# ── PATHS (EDIT THESE) ─────────────────────────────────────────────────────────
DIR_WORKING = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels"
DIR_OUTPUTS = join(DIR_WORKING, "outputs")
DIR_RF      = join(DIR_OUTPUTS, "Robinson_Foulds")   # contains <RNA> subfolders
DIR_DNA     = join(DIR_OUTPUTS, "DNAtrees")                  # fallback place to find a DNA tree

# File suffixes we expect inside each RNA folder
DNA_RFDIST_SUFFIX       = ".raxml.rfdist"        # DNA vs DNA
DNA_VS_RNA_IPSEU_SUFFIX = ".raxml.raxmlPi.rfdist"  # DNA vs RNA (ignore pseudoknots)

# ── OPTIONS ────────────────────────────────────────────────────────────────────
ALPHA = 0.05          # nominal alpha
DO_BONFERRONI = True  # set False to skip adjusted p-values

# ── Helpers ────────────────────────────────────────────────────────────────────
def read_rfdist_matrix(path: str) -> np.ndarray:
    """Read IQ-TREE .rfdist into a 2D array of floats. Returns empty array if nothing."""
    with open(path, "r") as f:
        lines = f.readlines()
    if len(lines) <= 1:
        return np.array([])
    rows = []
    for line in lines[1:]:
        parts = line.strip().split()
        if len(parts) <= 1:
            continue
        rows.append(list(map(float, parts[1:])))  # skip the row label
    return np.array(rows, dtype=float) if rows else np.array([])

def find_one_tree_for_taxa(rna_dir: str, rna: str, dir_dna: str) -> str | None:
    """
    Return a path to any tree file for this RNA to count taxa from.
    Priority: <rna_dir>/<RNA>.raxml → <rna_dir>/<RNA>.raxmlPi → DIR_DNA/<RNA>/RAxML_bestTree.*
    """
    candidates = [
        join(rna_dir, f"{rna}.raxml"),
        join(rna_dir, f"{rna}.raxmlPi"),
    ]
    dna_rna_dir = join(dir_dna, rna)
    if isdir(dna_rna_dir):
        # Pick first DNA tree file found
        for fn in sorted(os.listdir(dna_rna_dir)):
            if fn.startswith("RAxML_bestTree"):
                candidates.append(join(dna_rna_dir, fn))
                break
    for p in candidates:
        if os.path.exists(p):
            return p
    return None

def count_taxa(tree_path: str) -> int | None:
    """Count terminals in the first tree from a Newick file."""
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
    if mat.size == 0:
        return mat
    if np.nanmax(mat) <= 1.0 + 1e-9:
        return mat.copy()
    denom = 2 * (n_taxa - 3)
    if denom <= 0:
        return mat   # degenerate; return raw to avoid division by zero
    return mat / denom

def upper_triangle_values(mat: np.ndarray) -> np.ndarray:
    """Return the upper triangle (k=1) flattened; empty if invalid."""
    if mat.size == 0:
        return np.array([])
    n = mat.shape[0]
    iu = np.triu_indices(n, k=1)
    return mat[iu].astype(float)

def model_to_category(model):
    if model == "DNA_extra":
        return "DNA vs DNA"
    return f"DNA vs RNA ({model})"

# ── Core per-RNA computation ───────────────────────────────────────────────────
def process_one_rna(rna_dir: str, rna: str, model:str) -> dict | None:
    """
    Read, normalize, summarize, and U-test one RNA.
    Returns a dict of results, or None if cannot process.
    """
    dna_fp   = join(rna_dir, f"{rna}{DNA_RFDIST_SUFFIX}")
    ipseu_fp = join(rna_dir, f"{rna}{DNA_VS_RNA_IPSEU_SUFFIX}")
    if not (os.path.exists(dna_fp) and os.path.exists(ipseu_fp)):
        return None

    # Determine n from any available tree
    tree_path = find_one_tree_for_taxa(rna_dir, rna, DIR_DNA)
    if tree_path is None:
        return None
    n_taxa = count_taxa(tree_path)
    if not n_taxa or n_taxa < 4:
        return None

    # Read matrices and normalize by 2*(n-3) (unless already 0..1)
    dna_raw   = read_rfdist_matrix(dna_fp)
    ipseu_raw = read_rfdist_matrix(ipseu_fp)
    if dna_raw.size == 0 or ipseu_raw.size == 0:
        return None

    dna_nrf   = normalize_rf_matrix_by_n(dna_raw,   n_taxa)
    ipseu_nrf = normalize_rf_matrix_by_n(ipseu_raw, n_taxa)

    dna_vals   = upper_triangle_values(dna_nrf)    # 45 values
    rna_vals   = ipseu_nrf.flatten().astype(float) # 100 values
    if dna_vals.size == 0 or rna_vals.size == 0:
        return None

    # U-test (two-sided)
    try:
        res = mannwhitneyu(dna_vals, rna_vals, alternative="two-sided")
        U = float(getattr(res, "statistic", res[0]))
        p = float(getattr(res, "pvalue",    res[1]))
    except Exception:
        return None

    return {
        "Model": model,
        "RNA": rna,
        "n_taxa": int(n_taxa),
        "n_DNA":  int(dna_vals.size),
        "n_RNA":  int(rna_vals.size),
        "DNA_median_nRF": float(np.median(dna_vals)),
        "RNA_median_nRF": float(np.median(rna_vals)),
        "U": U,
        "pvalue": p,
    }

# ── Driver ─────────────────────────────────────────────────────────────────────
def main():
    if not isdir(DIR_RF):
        raise SystemExit(f"DIR_RF not found: {DIR_RF}")

    # discover models as subfolders in DIR_RF
    all_models = [m for m in sorted(os.listdir(DIR_RF)) if isdir(join(DIR_RF, m))]
    results = []
    med_rows=[]

    for model in all_models:
        modelRF_dir = join(DIR_RF, model)
        rna_ids = [d for d in sorted(os.listdir(modelRF_dir)) if isdir(join(modelRF_dir, d))]
        if not rna_ids: 
            continue
        
        model_rows = []
        for rna in rna_ids:
            rec = process_one_rna(join(modelRF_dir, rna), rna, model)
            if rec:
                model_rows.append(rec)

        if not model_rows:
            continue

        # Append long results
        results.extend(model_rows)

        # Per-model Bonferroni
        m = len(model_rows)
        for rec in model_rows:
            rec["p_bonf"] = min(1.0, rec["pvalue"] * m)
            rec["sig_bonf"] = rec["p_bonf"] < ALPHA

        # Build median rows for plotting
        cat_other = model_to_category(model)
        for rec in model_rows:
            med_rows.append({"Model": model, "RNA": rec["RNA"], "Category": "DNA vs DNA",       
                             "Median nRF": rec["DNA_median_nRF"]})
            med_rows.append({"Model": model, "RNA": rec["RNA"], "Category": cat_other,           
                             "Median nRF": rec["RNA_median_nRF"]})

        # quick console heartbeat
        sig_count = sum(1 for r in model_rows if r["sig_bonf"])
        print(f"{model}: {sig_count} / {m} Bonferroni-significant @ alpha={ALPHA}")

    if not results:
        print("No usable data found across models.")
        return

    df_long = pd.DataFrame(results).sort_values(["Model", "RNA"]).reset_index(drop=True)

    # Wide p-value table (RNA × Model)
    df_wide = (
        df_long.pivot_table(index="RNA", columns="Model", values="pvalue", aggfunc="first")
               .sort_index()
    )

    df_medians = pd.DataFrame(med_rows)

    #Write outputs
    out1 = join(DIR_RF, "Utest_all_models_long.csv")
    out2 = join(DIR_RF, "Utest_all_models_wide.csv")
    out3 = join(DIR_RF, "Median_nRF_all_models_long.csv")
    df_long.to_csv(out1, index=False)
    df_wide.to_csv(out2)
    df_medians.to_csv(out3, index=False)
    print(f"Wrote:\n  {out1}\n  {out2}\n  {out3}")


if __name__ == "__main__":
    main()