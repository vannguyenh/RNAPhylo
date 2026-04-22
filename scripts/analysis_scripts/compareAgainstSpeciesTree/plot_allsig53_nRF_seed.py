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
from os.path import join
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ── Config ──────────────────────────────────────────────────────────────────
INPUT_CSV = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/outputs/260306_allsig_species_tree_comparison/nRF_comparison.csv"
OUTPUT_DIR = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/outputs/260306_allsig_species_tree_comparison/figures"

MODEL_ORDER = ['S16', 'S16A', 'S16B',
               'S7A', 'S7B', 'S7C', 'S7D', 'S7E', 'S7F',
               'S6A', 'S6B', 'S6C', 'S6D', 'S6E']


def load_data():
    return pd.read_csv(INPUT_CSV)


# ── Figure 1: Boxplot of direct tree-to-tree comparisons per model ────────
def fig_direct_boxplot(df):
    """Per-model boxplot of nRF_rna_vs_dna, nRF_dna_vs_rfam, nRF_rna_vs_rfam."""
    models = [m for m in MODEL_ORDER if m in df['model'].unique()]
    n_fam = df['RNA_family'].nunique()

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    comparisons = [
        ('nRF_rna_vs_dna', 'RNA vs DNA', '#CE93D8'),
        ('nRF_dna_vs_rfam', 'DNA vs Rfam', '#90CAF9'),
        ('nRF_rna_vs_rfam', 'RNA vs Rfam', '#A5D6A7'),
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
            patch.set_edgecolor('black')
        for median in bp['medians']:
            median.set_color('black')
            median.set_linewidth(1.5)

        ax.set_title(title, fontsize=11)
        ax.set_ylim(-0.02, 1.05)
        ax.tick_params(axis='x', rotation=45)

    axes[0].set_ylabel('nRF', fontsize=11)
    fig.suptitle(f'Direct Tree-to-Tree Distances per Model\n'
                 f'({n_fam} families, Seed)', fontsize=13)
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

    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.bar(x, rna_closer, width, label='RNA tree closer to Rfam',
           color='#4CAF50', edgecolor='white', linewidth=0.5)
    ax.bar(x, ties, width, bottom=rna_closer, label='Tie',
           color='#BDBDBD', edgecolor='white', linewidth=0.5)
    ax.bar(x, dna_closer, width,
           bottom=[r + t for r, t in zip(rna_closer, ties)],
           label='DNA tree closer to Rfam',
           color='#F44336', edgecolor='white', linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(models, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Number of RNA families', fontsize=11)
    ax.set_title('Which Inferred Tree Is Closer to the Rfam Tree?\n'
                 '(per RNA substitution model, 53 families)', fontsize=12)
    ax.legend(loc='upper right', fontsize=9)

    for i, (rc, dc) in enumerate(zip(rna_closer, dna_closer)):
        if rc > 0:
            ax.text(i, rc / 2, str(rc), ha='center', va='center', fontsize=8,
                    fontweight='bold', color='white')
        if dc > 0:
            total = rna_closer[i] + ties[i]
            ax.text(i, total + dc / 2, str(dc), ha='center', va='center',
                    fontsize=8, fontweight='bold', color='white')

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


def _scatter_grid(df, x_col, y_col, x_label, y_label, suptitle, save_name):
    """Shared scatter grid: 3 rows grouped by model family."""
    fig, axes = plt.subplots(3, MAX_COLS, figsize=(4 * MAX_COLS, 4.5 * 3),
                             squeeze=False)

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
                       s=20, alpha=0.5, color='steelblue', edgecolor='white',
                       linewidth=0.3, zorder=3)

            lim = 1.08
            if not sub.empty:
                lim = max(sub[x_col].max(), sub[y_col].max(), 0.1) * 1.08
            ax.plot([0, lim], [0, lim], 'k--', lw=0.8, alpha=0.4, zorder=1)
            ax.fill_between([0, lim], [0, lim], [lim, lim], alpha=0.03,
                            color='red', zorder=0)
            ax.fill_between([0, lim], [0, 0], [0, lim], alpha=0.03,
                            color='green', zorder=0)

            n_below = int((sub[y_col] < sub[x_col]).sum())
            n_above = int((sub[y_col] > sub[x_col]).sum())
            n_eq = int((sub[y_col] == sub[x_col]).sum())

            ax.set_title(f'{model}  (R:{n_below} D:{n_above} T:{n_eq})',
                         fontsize=10)
            ax.set_xlabel(x_label, fontsize=8)
            ax.set_ylabel(y_label, fontsize=8)
            ax.set_xlim(-0.02, lim)
            ax.set_ylim(-0.02, lim)
            ax.set_aspect('equal')

    fig.suptitle(suptitle, fontsize=13, y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
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


# ── Figure 4: Scatter nRF(RNA vs Sp) vs nRF(DNA vs Sp) per model ─────────
def fig_scatter_vs_species(df):
    valid = df.dropna(subset=['nRF_dna_vs_sp', 'nRF_rna_vs_sp'])
    if valid.empty:
        return
    n_fam = valid['RNA_family'].nunique()
    _scatter_grid(
        valid,
        x_col='nRF_dna_vs_sp', y_col='nRF_rna_vs_sp',
        x_label='nRF (DNA vs NCBI Taxonomy)',
        y_label='nRF (RNA vs NCBI Taxonomy)',
        suptitle=f'Normalised Robinson-Foulds Distance to NCBI Taxonomy Tree '
                 f'— per RNA Model ({n_fam} families with ≥4 unique species)',
        save_name='fig4_scatter_vs_species',
    )


# ── Figure 5: Scatter nIC(RNA vs Sp) vs nIC(DNA vs Sp) per model ─────────
def fig_scatter_ic_vs_species(df):
    """Per-model scatter: nIC of RNA vs species tree vs nIC of DNA vs species tree."""
    if 'nIC_dna_vs_sp' not in df.columns:
        return
    valid = df.dropna(subset=['nIC_dna_vs_sp', 'nIC_rna_vs_sp'])
    if valid.empty:
        return
    n_fam = valid['RNA_family'].nunique()
    _scatter_grid(
        valid,
        x_col='nIC_dna_vs_sp', y_col='nIC_rna_vs_sp',
        x_label='nIC (DNA vs NCBI Taxonomy)',
        y_label='nIC (RNA vs NCBI Taxonomy)',
        suptitle=f'Incompatible Splits vs NCBI Taxonomy Tree '
                 f'— per RNA Model ({n_fam} families)',
        save_name='fig5_scatter_nIC_vs_species',
    )


# ── Figure 6: Scatter nIC(DNA vs Sp) vs nIC(Rfam vs Sp) per model ────────
def fig_scatter_ic_dna_vs_rfam_sp(df):
    """Per-model scatter: is DNA or Rfam tree closer to species tree (nIC)?"""
    if 'nIC_dna_vs_sp' not in df.columns or 'nIC_rfam_vs_sp' not in df.columns:
        return
    valid = df.dropna(subset=['nIC_dna_vs_sp', 'nIC_rfam_vs_sp'])
    if valid.empty:
        return
    n_fam = valid['RNA_family'].nunique()
    _scatter_grid(
        valid,
        x_col='nIC_rfam_vs_sp', y_col='nIC_dna_vs_sp',
        x_label='nIC (Rfam vs Species)', y_label='nIC (DNA vs Species)',
        suptitle=f'Incompatible Splits vs Species Tree: DNA vs Rfam ({n_fam} families)',
        save_name='fig6_scatter_nIC_dna_vs_rfam_sp',
    )


# ── Figure 7: Scatter nIC(RNA vs Sp) vs nIC(Rfam vs Sp) per model ────────
def fig_scatter_ic_rna_vs_rfam_sp(df):
    """Per-model scatter: is RNA or Rfam tree closer to species tree (nIC)?"""
    if 'nIC_rna_vs_sp' not in df.columns or 'nIC_rfam_vs_sp' not in df.columns:
        return
    valid = df.dropna(subset=['nIC_rna_vs_sp', 'nIC_rfam_vs_sp'])
    if valid.empty:
        return
    n_fam = valid['RNA_family'].nunique()
    _scatter_grid(
        valid,
        x_col='nIC_rfam_vs_sp', y_col='nIC_rna_vs_sp',
        x_label='nIC (Rfam vs Species)', y_label='nIC (RNA vs Species)',
        suptitle=f'Incompatible Splits vs Species Tree: RNA vs Rfam ({n_fam} families)',
        save_name='fig7_scatter_nIC_rna_vs_rfam_sp',
    )


def _save(fig, name):
    for ext in ['png', 'pdf']:
        path = join(OUTPUT_DIR, f'{name}.{ext}')
        fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  Saved {name}")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = load_data()
    print(f"Loaded {len(df)} rows, {df.RNA_family.nunique()} families\n")

    # nRF and nIC vs taxonomy tree (per model scatter plots)
    fig_scatter_vs_species(df)       # nRF (from --normalize-dist)
    fig_scatter_ic_vs_species(df)    # nIC

    print(f"\nAll figures saved to: {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
