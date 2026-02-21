#!/usr/bin/env python3
"""
plot_rf_density_by_rna_type.py
==============================
Faceted KDE density plot of normalised RF distances, split by RNA structural
type (rows) and RNA substitution model (columns).

Each panel shows two curves:
  Red  : nRF(T_rna, T_dna)   — RNA model trees vs DNA trees
  Teal : nRF(T_dna, T_dna2)  — DNA baseline (seeds 1–10 vs 11–20)

The separation between curves answers: "does this RNA type benefit from
RNA-specific substitution models?"

RNA types are ordered top-to-bottom by expected structural complexity
(i.e. where RNA models should matter most → least).

Usage
-----
    python plot_rf_density_by_rna_type.py

Dependencies
------------
    pip install numpy matplotlib scipy biopython pandas
"""

import os
import glob
import logging
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde
from Bio import Phylo

# =============================================================================
# CONFIGURATION — edit paths here
# =============================================================================

DIR_WORKING   = os.path.expanduser("~/RNAPhylo/seedAlignment_AllModels")
DIR_OUTPUTS   = os.path.join(DIR_WORKING, "outputs")
DIR_RF        = os.path.join(DIR_OUTPUTS, "260203_Robinson_Foulds")
DIR_DNA_TREES = os.path.join(DIR_OUTPUTS, "inferred_trees", "DNA")

# Rfam metadata files (adjust paths if needed)
RFAM_TBL         = os.path.join(DIR_WORKING, "inputs/Rfam.full.seed.tbl")
RFAM_FAMILY_TXT  = os.path.join(DIR_WORKING, "inputs/family.txt")   # the file you uploaded

OUTPUT_FIG = os.path.join(DIR_OUTPUTS, "rf_density_by_rna_type.pdf")

# --- RNA type groups ----------------------------------------------------------
# Ordered by structural complexity (most structured → least).
# CHANGE 1: added "Other" at the bottom (includes tRNA n=2, CRISPR, unknown)
RNA_TYPE_ORDER = [
    "rRNA",
    "Ribozyme",
    "snRNA",
    "Intron",
    "snoRNA",
    "Cis-reg / Riboswitch",
    "IRES",
    "sRNA / Antisense",
    "miRNA",
    "lncRNA",
    "Other",
]

# Mapping from Rfam 'type' field → display group name
TYPE_MAP = {
    "Gene; rRNA;":                       "rRNA",
    "Gene; tRNA;":                       "Other",         # CHANGE 1: n=2, merged into Other
    "Gene; ribozyme;":                   "Ribozyme",
    "Gene; snRNA; splicing;":            "snRNA",
    "Gene; snRNA;":                      "snRNA",
    "Gene; snRNA; snoRNA; CD-box;":      "snoRNA",
    "Gene; snRNA; snoRNA; HACA-box;":    "snoRNA",
    "Gene; snRNA; snoRNA; scaRNA;":      "snoRNA",
    "Cis-reg; riboswitch;":              "Cis-reg / Riboswitch",
    "Cis-reg; leader;":                  "Cis-reg / Riboswitch",
    "Cis-reg; thermoregulator;":         "Cis-reg / Riboswitch",
    "Cis-reg; frameshift_element;":      "Cis-reg / Riboswitch",
    "Cis-reg;":                          "Cis-reg / Riboswitch",
    "Cis-reg; IRES;":                    "IRES",
    "Gene; miRNA;":                      "miRNA",
    "Gene; sRNA;":                       "sRNA / Antisense",
    "Gene; antisense;":                  "sRNA / Antisense",
    "Gene; antitoxin;":                  "sRNA / Antisense",
    "Intron;":                           "Intron",
    "Gene; lncRNA;":                     "lncRNA",
    "Gene; CRISPR;":                     "Other",
    "Gene;":                             "Other",
}

# Colours
COLOR_RNA = "#E05C5C"   # red  — nRF(T_rna, T_dna)
COLOR_DNA = "#4BBFBF"   # teal — nRF(T_dna, T_dna2)

BW_METHOD  = "scott"
KDE_POINTS = 500

# =============================================================================
# TAXONOMY / TYPE MAP BUILDER
# =============================================================================

logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
log = logging.getLogger(__name__)


def build_rna_type_map() -> dict[str, str]:
    """
    Returns {rfam_acc: rna_group} using family.txt (col 18 = type).
    Falls back to Rfam_full_seed.tbl ID-name heuristics if family.txt missing.
    """
    if os.path.isfile(RFAM_FAMILY_TXT):
        log.info(f"Reading RNA types from {RFAM_FAMILY_TXT} ...")
        df = pd.read_csv(
            RFAM_FAMILY_TXT, sep="\t", header=None,
            encoding="latin1", on_bad_lines="skip",
            usecols=[0, 18]
        )
        df.columns = ["rfam_acc", "type"]
        df["type"] = df["type"].astype(str).str.strip()
        df["group"] = df["type"].map(TYPE_MAP).fillna("Other")
        result = dict(zip(df["rfam_acc"], df["group"]))
        counts = df["group"].value_counts()
        log.info(f"  Type distribution:\n{counts.to_string()}")
        return result

    # Fallback: name heuristics from .tbl file
    import re
    log.warning("family.txt not found — using name heuristics.")
    tbl = pd.read_csv(RFAM_TBL, sep="\t", usecols=["AC", "ID"])
    def heuristic(name):
        n = name.lower()
        if re.search(r'rrna|_rrna|16s|23s|18s|28s|5s_r|ssu_r|lsu_r', n): return "rRNA"
        if re.search(r'ribozyme', n): return "Ribozyme"
        if re.search(r'^u\d|snrna|splicing', n): return "snRNA"
        if re.search(r'snor|snoRNA', n, re.I): return "snoRNA"
        if re.search(r'intron', n): return "Intron"
        if re.search(r'riboswitch|leader|thermometer|frameshift', n): return "Cis-reg / Riboswitch"
        if re.search(r'ires', n): return "IRES"
        if re.search(r'^mir-|let-7|mirna', n): return "miRNA"
        if re.search(r'srna|antisense|csrb|6s\b', n): return "sRNA / Antisense"
        if re.search(r'lncrna', n): return "lncRNA"
        return "Other"
    return {row["AC"]: heuristic(row["ID"]) for _, row in tbl.iterrows()}


# =============================================================================
# RF DISTANCE HELPERS
# =============================================================================

def count_taxa(locus: str) -> int | None:
    pattern = os.path.join(DIR_DNA_TREES, locus, "RAxML_bestTree.*")
    candidates = [p for p in glob.glob(pattern) if os.path.isfile(p)]
    if not candidates:
        return None
    try:
        with open(candidates[0]) as fh:
            tree = next(Phylo.parse(fh, "newick"))
        return len(tree.get_terminals())
    except Exception as e:
        log.debug(f"taxa count failed for {locus}: {e}")
        return None


def read_rfdist(path: str) -> np.ndarray | None:
    try:
        with open(path) as fh:
            lines = fh.readlines()
        rows = []
        for line in lines[1:]:
            parts = line.strip().split()
            if len(parts) > 1:
                rows.append(list(map(float, parts[1:])))
        return np.array(rows, dtype=float) if rows else None
    except Exception as e:
        log.warning(f"cannot read {path}: {e}")
        return None


def normalise(mat: np.ndarray, n_taxa: int) -> np.ndarray | None:
    if n_taxa < 4:
        return None
    denom = 2 * (n_taxa - 3)
    return mat / denom if denom > 0 else None


def collect_by_type(model: str,
                    type_map: dict[str, str]) -> dict[str, list[float]]:
    """
    Returns {rna_group: [nRF values]} for a given model directory,
    with loci split by RNA type. Off-diagonal only (no self-comparisons).
    """
    model_dir = os.path.join(DIR_RF, model)
    loci = sorted(d for d in os.listdir(model_dir)
                  if os.path.isdir(os.path.join(model_dir, d)))

    result: dict[str, list[float]] = {g: [] for g in RNA_TYPE_ORDER + ["Other"]}
    skipped = 0

    for locus in loci:
        rfdist_path = os.path.join(model_dir, locus, f"{locus}.rfdist")
        if not os.path.isfile(rfdist_path):
            skipped += 1
            continue
        mat = read_rfdist(rfdist_path)
        if mat is None:
            skipped += 1
            continue
        n_taxa = count_taxa(locus)
        if n_taxa is None:
            skipped += 1
            continue
        normed = normalise(mat, n_taxa)
        if normed is None:
            skipped += 1
            continue

        # Off-diagonal values only
        mask = ~np.eye(normed.shape[0], normed.shape[1], dtype=bool)
        vals  = normed[mask].tolist()

        group = type_map.get(locus, "Other")
        if group not in result:
            group = "Other"
        result[group].extend(vals)

    log.info(f"  {model}: {len(loci)-skipped}/{len(loci)} loci, "
             f"{skipped} skipped")
    return result


# =============================================================================
# PLOT
# =============================================================================

def plot_kde(ax, vals, color, bw="scott"):
    if len(vals) < 5:
        return False
    arr = np.array(vals, dtype=float)
    kde = gaussian_kde(arr, bw_method=bw)
    x = np.linspace(0, 1, KDE_POINTS)
    y = kde(x)
    ax.plot(x, y, color=color, linewidth=1.8)
    ax.fill_between(x, y, alpha=0.15, color=color)
    return True


def main():
    # ── 1. Build RNA type map ─────────────────────────────────────────────────
    type_map = build_rna_type_map()

    # ── 2. Discover RNA models ────────────────────────────────────────────────
    MODEL_ORDER = [
        "S16", "S16A", "S16B",
        "S7A", "S7B", "S7C", "S7D", "S7E", "S7F",
        "S6A", "S6B", "S6C", "S6D", "S6E",
    ]
    available = {
        d for d in os.listdir(DIR_RF)
        if os.path.isdir(os.path.join(DIR_RF, d)) and d != "DNA_extra"
    }
    rna_models = [m for m in MODEL_ORDER if m in available]
    log.info(f"RNA models: {rna_models}")

    # ── 3. DNA baseline split by type ────────────────────────────────────────
    log.info("Computing DNA baseline (DNA_2)...")
    dna_by_type = collect_by_type("DNA_2", type_map)

    # ── 4. RNA distances split by type ───────────────────────────────────────
    rna_by_type: dict[str, dict[str, list]] = {}
    for model in rna_models:
        log.info(f"Processing {model}...")
        rna_by_type[model] = collect_by_type(model, type_map)

    # ── 5. Determine active rows (types with enough data) ────────────────────
    MIN_VALS = 5
    active_types = [
        g for g in RNA_TYPE_ORDER
        if any(len(rna_by_type[m].get(g, [])) >= MIN_VALS
               for m in rna_models)
    ]
    log.info(f"Active RNA types ({len(active_types)}): {active_types}")

    n_rows = len(active_types)
    n_cols = len(rna_models)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * 3.2, n_rows * 2.8),
        sharex=True, sharey=False
    )
    axes = np.array(axes).reshape(n_rows, n_cols)

    for row_idx, rna_type in enumerate(active_types):
        for col_idx, model in enumerate(rna_models):
            ax = axes[row_idx, col_idx]

            has_dna = plot_kde(ax, dna_by_type.get(rna_type, []),
                               color=COLOR_DNA)
            has_rna = plot_kde(ax, rna_by_type[model].get(rna_type, []),
                               color=COLOR_RNA)

            # ggplot-style panel
            ax.set_facecolor("#F2F2F2")
            ax.grid(True, color="white", linewidth=0.7, linestyle="-")
            ax.spines[:].set_visible(False)
            ax.tick_params(labelsize=11)
            ax.set_xlim(0, 1)
            ax.set_ylim(bottom=0)

            # Column header — top row only
            if row_idx == 0:
                ax.set_title(model, fontsize=13, fontweight="bold", pad=5)

            # Row label — rightmost column only
            if col_idx == n_cols - 1:
                n_loci = len(rna_by_type[model].get(rna_type, [])) // 90  # approx
                ax.annotate(
                    rna_type,
                    xy=(1.02, 0.5), xycoords="axes fraction",
                    fontsize=12, fontweight="bold",
                    va="center", ha="left", rotation=0,
                    annotation_clip=False
                )

            # CHANGE 2: removed per-panel set_ylabel / set_xlabel

    # ── 6. Legend ─────────────────────────────────────────────────────────────
    legend_handles = [
        mpatches.Patch(color=COLOR_RNA, alpha=0.85,
                       label=r"$nRF(T_{\mathrm{rna}},\,T_{\mathrm{dna}})$"),
        mpatches.Patch(color=COLOR_DNA, alpha=0.85,
                       label=r"$nRF(T_{\mathrm{dna}},\,T_{\mathrm{dna2}})$"),
    ]
    fig.legend(
        handles=legend_handles, loc="upper center", ncol=2,
        fontsize=16, frameon=True, framealpha=0.9,
        bbox_to_anchor=(0.5, 1.02)
    )
    fig.suptitle(
        "Normalised RF distance by RNA structural type and substitution model",
        fontsize=20, y=1.05
    )

    plt.tight_layout(rect=[0.03, 0.03, 0.88, 1])

    # CHANGE 2: single shared axis labels
    fig.text(0.44, 0.01, "Normalised RF distance", ha="center", fontsize=16)
    fig.text(0.01, 0.5,  "Density", va="center", rotation=90, fontsize=16)

    # ── 7. Save ───────────────────────────────────────────────────────────────
    fig.savefig(OUTPUT_FIG, dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_FIG.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    log.info(f"Saved: {OUTPUT_FIG}")
    plt.show()


if __name__ == "__main__":
    main()