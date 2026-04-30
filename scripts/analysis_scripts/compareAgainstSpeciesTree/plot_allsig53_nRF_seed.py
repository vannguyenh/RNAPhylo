#!/usr/bin/env python3
"""
Publication-ready figures summarising the nRF tree comparison
for the 53 all-significant RNA families (seed data).

Figures:
  1. Boxplot: nRF for direct tree-to-tree comparisons (all 53 families)
  2. Per-model bar chart: RNA closer vs DNA closer vs Tie (direct vs Rfam)
  3. Scatter: nRF(RNA vs Rfam) vs nRF(DNA vs Rfam), one point per family
  4. Boxplot: nRF vs species tree (42 families with enough unique species)

Usage:
    python plot_allsig53_nRF.py
"""

import os
import sys
from datetime import datetime
from os.path import join
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Shared manuscript-wide plot style
sys.path.insert(0, str(next(
    p for p in Path(__file__).resolve().parents if p.name == "scripts"
)))
from plot_style import (
    apply_style, COLORS, FIG_DOUBLE, diagonal_scatter, save_figure,
)
apply_style()

# ── Config ──────────────────────────────────────────────────────────────────
INPUT_CSV = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/outputs/260306_allsig_species_tree_comparison/nRF_comparison.csv"
OUTPUT_DIR = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/outputs/260306_allsig_species_tree_comparison/figures"

MODEL_ORDER = ['S16', 'S16A', 'S16B',
               'S7A', 'S7B', 'S7C', 'S7D', 'S7E', 'S7F',
               'S6A', 'S6B', 'S6C', 'S6D', 'S6E']

RESOLUTION_THRESHOLD = 0.5
DATASET_TAG = "SEED"
DATE_TAG = datetime.now().strftime("%y%m%d")


def load_data():
    return pd.read_csv(INPUT_CSV)


# ── Figure 1: Boxplot of direct tree-to-tree comparisons per model ────────
def fig_direct_boxplot(df):
    """Per-model boxplot of nRF_rna_vs_dna, nRF_dna_vs_rfam, nRF_rna_vs_rfam."""
    models = [m for m in MODEL_ORDER if m in df['model'].unique()]
    n_fam = df['RNA_family'].nunique()

    fig, axes = plt.subplots(1, 3, figsize=(FIG_DOUBLE, 2.6), sharey=True)
    # Colour each panel by the "non-DNA" tree being compared: RNA-vs-DNA
    # uses the RNA colour, DNA-vs-Rfam uses DNA, RNA-vs-Rfam uses Rfam.
    comparisons = [
        ('nRF_rna_vs_dna',  'RNA vs DNA',  COLORS['RNA']),
        ('nRF_dna_vs_rfam', 'DNA vs Rfam', COLORS['DNA']),
        ('nRF_rna_vs_rfam', 'RNA vs Rfam', COLORS['Rfam']),
    ]

    for ax, (col, title, color) in zip(axes, comparisons):
        data = [df[(df['model'] == m)][col].dropna().values for m in models]
        bp = ax.boxplot(
            data, tick_labels=models, patch_artist=True, widths=0.6,
            showfliers=True,
            flierprops=dict(marker='o', markersize=2, alpha=0.3),
        )
        for patch in bp['boxes']:
            patch.set_facecolor(color)
            patch.set_alpha(0.55)
            patch.set_edgecolor(color)
        for median in bp['medians']:
            median.set_color('black')
            median.set_linewidth(1.2)

        ax.set_title(title)
        ax.set_ylim(-0.02, 1.05)
        ax.tick_params(axis='x', rotation=45)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)

    axes[0].set_ylabel('nRF')
    fig.suptitle(f'Direct tree-to-tree distances per model '
                 f'({n_fam} families, seed)')
    plt.tight_layout()
    _save(fig, 'fig1_direct_boxplot')


# ── Figure 2: Per-model bar chart (RNA vs Rfam closer) ───────────────────
def fig_per_model_rfam_bars(df):
    """Per model: is RNA or DNA tree closer to the Rfam tree?"""
    models = [m for m in MODEL_ORDER if m in df['model'].unique()]
    rna_closer, dna_closer, ties = [], [], []

    for model in models:
        sub = df[df['model'] == model].dropna(subset=['nRF_rna_vs_rfam', 'nRF_dna_vs_rfam'])
        rna_closer.append(int((sub['nRF_rna_vs_rfam'] < sub['nRF_dna_vs_rfam']).sum()))
        dna_closer.append(int((sub['nRF_rna_vs_rfam'] > sub['nRF_dna_vs_rfam']).sum()))
        ties.append(int((sub['nRF_rna_vs_rfam'] == sub['nRF_dna_vs_rfam']).sum()))

    x = np.arange(len(models))
    width = 0.6

    fig, ax = plt.subplots(figsize=(FIG_DOUBLE, 3.0))
    ax.bar(x, rna_closer, width, label='RNA tree closer to Rfam',
           color=COLORS['RNA'], edgecolor='white', linewidth=0.5)
    ax.bar(x, ties, width, bottom=rna_closer, label='Tie',
           color=COLORS['tie'], edgecolor='white', linewidth=0.5)
    ax.bar(x, dna_closer, width,
           bottom=[r + t for r, t in zip(rna_closer, ties)],
           label='DNA tree closer to Rfam',
           color=COLORS['DNA'], edgecolor='white', linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(models, rotation=45, ha='right')
    ax.set_ylabel('Number of RNA families')
    ax.set_title('Which inferred tree is closer to the Rfam tree? '
                 '(per RNA substitution model, 53 families)')
    ax.legend(loc='upper right')
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)

    for i, (rc, dc) in enumerate(zip(rna_closer, dna_closer)):
        if rc > 0:
            ax.text(i, rc / 2, str(rc), ha='center', va='center',
                    fontweight='bold', color='white')
        if dc > 0:
            total = rna_closer[i] + ties[i]
            ax.text(i, total + dc / 2, str(dc), ha='center', va='center',
                    fontweight='bold', color='white')

    plt.tight_layout()
    _save(fig, 'fig2_per_model_vs_rfam')


# ── Layout: 3 rows grouped by model family ───────────────────────────────
# Row 1: S16, S16A, S16B (3 cols)
# Row 2: S7A, S7B, S7C, S7D, S7E, S7F (6 cols)
# Row 3: S6A, S6B, S6C, S6D, S6E (5 cols)
MODEL_ROWS = [
    ['S16', 'S16A', 'S16B'],
    ['S7A', 'S7B', 'S7C', 'S7D', 'S7E', 'S7F'],
    ['S6A', 'S6B', 'S6C', 'S6D', 'S6E'],
]
MAX_COLS = 6  # widest row


PANEL_SIZE = 1.9  # inches per scatter panel (square; nRF is on a 0–1 scale)


def _scatter_grid(df, x_col, y_col, x_label, y_label, suptitle, save_name):
    """Shared scatter grid: 3 rows grouped by model family."""
    fig_w = MAX_COLS * PANEL_SIZE
    fig_h = len(MODEL_ROWS) * PANEL_SIZE + 0.8  # extra room for suptitle
    fig, axes = plt.subplots(len(MODEL_ROWS), MAX_COLS,
                             figsize=(fig_w, fig_h), squeeze=False)

    for row_idx, row_models in enumerate(MODEL_ROWS):
        for col_idx in range(MAX_COLS):
            ax = axes[row_idx][col_idx]
            if col_idx >= len(row_models):
                ax.set_visible(False)
                continue

            model = row_models[col_idx]
            if model not in df['model'].unique():
                ax.set_visible(False)
                continue

            sub = df[df['model'] == model].dropna(subset=[x_col, y_col])

            ax.scatter(sub[x_col], sub[y_col],
                       s=14, alpha=0.65, color=COLORS['RNA'],
                       edgecolor='white', linewidth=0.3, zorder=3)

            lim = 1.08
            if not sub.empty:
                lim = max(sub[x_col].max(), sub[y_col].max(), 0.1) * 1.08
            diagonal_scatter(ax, lim=lim)

            n_below = int((sub[y_col] < sub[x_col]).sum())
            n_above = int((sub[y_col] > sub[x_col]).sum())
            n_eq = int((sub[y_col] == sub[x_col]).sum())

            ax.set_title(f'{model}  (R:{n_below} D:{n_above} T:{n_eq})')

    fig.suptitle(suptitle, y=0.995)
    fig.supxlabel(x_label, y=0.02)
    fig.supylabel(y_label, x=0.005)
    plt.tight_layout(rect=[0.02, 0.04, 1, 0.96])
    _save(fig, save_name)


# ── Figure 3: Scatter nRF(RNA vs Rfam) vs nRF(DNA vs Rfam) per model ─────
def fig_scatter_vs_rfam(df):
    n_fam = df['RNA_family'].nunique()
    _scatter_grid(
        df,
        x_col='nRF_dna_vs_rfam', y_col='nRF_rna_vs_rfam',
        x_label='nRF (DNA vs Rfam)', y_label='nRF (RNA vs Rfam)',
        suptitle=f'Distance to Rfam Tree — per RNA Model ({n_fam} families)',
        save_name='fig3_scatter_vs_rfam',
    )


# ── Figure 4: Scatter nRF(RNA vs NCBI) vs nRF(DNA vs NCBI), coloured by NCBI resolution ──
def fig_scatter_vs_species(df):
    """Scatter grid coloured by NCBI tree resolution.

    Filled circles: families where NCBI taxonomy resolves >= RESOLUTION_THRESHOLD
                    of internal branches (informative reference).
    Open circles:    families where NCBI is largely a polytomy (caveat applies).
    """
    valid = df.dropna(subset=['nRF_dna_vs_sp', 'nRF_rna_vs_sp']).copy()
    if valid.empty:
        return
    n_fam = valid['RNA_family'].nunique()
    fam_res = valid.groupby('RNA_family')['ncbi_resolution'].first()
    n_well = int((fam_res >= RESOLUTION_THRESHOLD).sum())
    n_poor = int((fam_res < RESOLUTION_THRESHOLD).sum())

    fig_w = MAX_COLS * PANEL_SIZE
    fig_h = len(MODEL_ROWS) * PANEL_SIZE + 1.2
    fig, axes = plt.subplots(len(MODEL_ROWS), MAX_COLS,
                             figsize=(fig_w, fig_h), squeeze=False,
                             gridspec_kw={'hspace': 0.45, 'wspace': 0.3})

    for row_idx, row_models in enumerate(MODEL_ROWS):
        for col_idx in range(MAX_COLS):
            ax = axes[row_idx][col_idx]
            if col_idx >= len(row_models):
                ax.set_visible(False)
                continue

            model = row_models[col_idx]
            if model not in valid['model'].unique():
                ax.set_visible(False)
                continue

            sub = valid[valid['model'] == model]
            well = sub[sub['ncbi_resolution'] >= RESOLUTION_THRESHOLD]
            poor = sub[sub['ncbi_resolution'] < RESOLUTION_THRESHOLD]

            # Filled markers = informative reference; open markers = polytomy-rich
            ax.scatter(well['nRF_dna_vs_sp'], well['nRF_rna_vs_sp'],
                       s=18, alpha=0.75,
                       c=COLORS['NCBI'], edgecolor=COLORS['NCBI'],
                       linewidth=0.5, zorder=3)
            ax.scatter(poor['nRF_dna_vs_sp'], poor['nRF_rna_vs_sp'],
                       s=18, alpha=0.7,
                       facecolor='none', edgecolor=COLORS['tie'],
                       linewidth=0.7, zorder=2)

            diagonal_scatter(ax, lim=1.08)

            n_below = int((well['nRF_rna_vs_sp'] < well['nRF_dna_vs_sp']).sum())
            n_above = int((well['nRF_rna_vs_sp'] > well['nRF_dna_vs_sp']).sum())
            n_eq = int((well['nRF_rna_vs_sp'] == well['nRF_dna_vs_sp']).sum())

            ax.set_title(f'{model}  (R:{n_below} D:{n_above} T:{n_eq})')

    pct = int(RESOLUTION_THRESHOLD * 100)
    legend_handles = [
        plt.Line2D([], [], marker='o', linestyle='',
                   markerfacecolor=COLORS['NCBI'],
                   markeredgecolor=COLORS['NCBI'], markersize=5,
                   label=f'NCBI resolution ≥{pct}% (n={n_well})'),
        plt.Line2D([], [], marker='o', linestyle='',
                   markerfacecolor='none',
                   markeredgecolor=COLORS['tie'], markersize=5,
                   label=f'NCBI resolution <{pct}% (n={n_poor})'),
    ]
    # Park the legend in the empty upper-right cells of row 1 (cols 3–5
    # are hidden because row 1 only has S16, S16A, S16B).
    fig.legend(handles=legend_handles, loc='center', ncol=1,
               frameon=True, framealpha=0.9,
               title='NCBI resolution', title_fontsize=9,
               bbox_to_anchor=(0.78, 0.83))

    fig.suptitle((f'Normalised Robinson-Foulds distance to NCBI taxonomy tree '
                  f'— per RNA model ({n_fam} families)'),
                 y=0.99)
    fig.supxlabel('nRF (DNA vs NCBI Taxonomy)', y=0.015)
    fig.supylabel('nRF (RNA vs NCBI Taxonomy)', x=0.005)
    fig.subplots_adjust(top=0.94, bottom=0.08, left=0.05, right=0.99)
    _save(fig, f'{DATE_TAG}_nRF_Taxonomy_{DATASET_TAG}')


def print_stratified_summary(df):
    valid = df.dropna(subset=['nRF_dna_vs_sp', 'nRF_rna_vs_sp']).copy()
    print("\n=== Stratified summary by NCBI resolution (seed) ===")
    for label, mask in [
        ('All', valid['ncbi_resolution'].notna()),
        (f'NCBI resolution ≥ {int(RESOLUTION_THRESHOLD*100)}%',
         valid['ncbi_resolution'] >= RESOLUTION_THRESHOLD),
        (f'NCBI resolution < {int(RESOLUTION_THRESHOLD*100)}%',
         valid['ncbi_resolution'] < RESOLUTION_THRESHOLD),
    ]:
        sub = valid[mask]
        n_fam = sub['RNA_family'].nunique()
        print(f"\n  {label} ({n_fam} families):")
        print(f"    {'Model':<6} {'mean DNA':>9} {'mean RNA':>9} {'R':>4} {'D':>4} {'T':>4}")
        for m in MODEL_ORDER:
            ms = sub[sub['model'] == m]
            if ms.empty:
                continue
            r = int((ms['nRF_rna_vs_sp'] < ms['nRF_dna_vs_sp']).sum())
            d = int((ms['nRF_rna_vs_sp'] > ms['nRF_dna_vs_sp']).sum())
            t = int((ms['nRF_rna_vs_sp'] == ms['nRF_dna_vs_sp']).sum())
            print(f"    {m:<6} {ms['nRF_dna_vs_sp'].mean():>9.4f} "
                  f"{ms['nRF_rna_vs_sp'].mean():>9.4f} {r:>4} {d:>4} {t:>4}")


def _save(fig, name):
    save_figure(fig, join(OUTPUT_DIR, f'{name}.pdf'))
    print(f"  Saved {name}")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = load_data()
    print(f"Loaded {len(df)} rows, {df.RNA_family.nunique()} families\n")

    fig_scatter_vs_species(df)       # nRF vs NCBI (coloured by resolution)
    print_stratified_summary(df)

    print(f"\nAll figures saved to: {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
