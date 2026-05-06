"""
Plot the number of significant RNA families per model from AU-test results.

For seed data, reads the pre-computed significance table CSV.
For full data, parses raw CONSEL .pv files to build the table.

Usage:
    python figure_AUtest_significance.py --dataset seed
    python figure_AUtest_significance.py --dataset full
"""

import argparse
import os
import sys
from datetime import datetime
from os.path import join, isdir
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

DATE_TAG = datetime.now().strftime("%y%m%d")

# Shared manuscript-wide plot style
sys.path.insert(0, str(next(
    p for p in Path(__file__).resolve().parents if p.name == "scripts"
)))
from plot_style import apply_style, COLORS, FIG_DOUBLE, save_figure
apply_style()

DATASET_CONFIGS = {
    'seed': {
        'working_dir': '/Users/u7875558/RNAPhylo/seedAlignment_AllModels',
        'au_folder': '260208_AU_Test_RAxML_RNAmodels',
        'sig_table_csv': '260210_AU_significant_table.csv',
    },
    'full': {
        'working_dir': '/Users/u7875558/RNAPhylo/fullAlignment',
        'au_folder': '260208_AU_Test_RAxML_RNAmodels',
        'sig_table_csv': '260210_AU_significant_table.csv',
    },
}

ALPHA = 0.05

# Canonical display order for models (S16 family, then S7 family, then S6 family)
MODEL_ORDER = [
    'S16', 'S16A', 'S16B',
    'S7A', 'S7B', 'S7C', 'S7D', 'S7E', 'S7F',
    'S6A', 'S6B', 'S6C', 'S6D', 'S6E',
]


def parse_pv_file(pv_path):
    """
    Parse a CONSEL .pv file and return (p_DNA, p_RNA) AU test p-values.

    The PV matrix has 2 rows x 8 columns. Column 0 is the AU p-value.
    Row 0 = DNA tree, Row 1 = RNA tree.
    """
    with open(pv_path) as f:
        lines = f.readlines()

    in_pv = False
    skip_next = False
    new_row = False
    row_values = []
    for line in lines:
        stripped = line.strip()
        if stripped == '# PV:':
            in_pv = True
            continue
        if in_pv and stripped.startswith('#!MAT:'):
            skip_next = True
            continue
        if skip_next:
            skip_next = False
            continue
        if in_pv and stripped.startswith('# row:'):
            new_row = True
            continue
        if in_pv and stripped.startswith('#'):
            if stripped == '# SE:':
                break
            continue
        if in_pv and new_row and stripped and not stripped.startswith('#'):
            values = stripped.split()
            row_values.append(float(values[0]))
            new_row = False

    if len(row_values) >= 2:
        return row_values[0], row_values[1]
    return None, None


def build_significance_table_from_pv(dir_au):
    """
    Parse all .pv files under dir_au/{model}/{RNA}/{RNA}_consel.pv
    and build a boolean significance table (RNA_family x models).
    """
    models = sorted(
        m for m in os.listdir(dir_au)
        if isdir(join(dir_au, m)) and m.startswith('S')
    )
    records = {}

    for model in models:
        model_dir = join(dir_au, model)
        for rna in sorted(os.listdir(model_dir)):
            rna_dir = join(model_dir, rna)
            if not isdir(rna_dir):
                continue
            pv_file = join(rna_dir, f"{rna}_consel.pv")
            if not os.path.exists(pv_file):
                continue
            p_dna, p_rna = parse_pv_file(pv_file)
            if p_dna is None:
                continue
            p_min = min(p_dna, p_rna)
            if rna not in records:
                records[rna] = {}
            records[rna][model] = p_min < ALPHA

    df = pd.DataFrame.from_dict(records, orient='index')
    df.index.name = 'RNA_family'
    ordered = [m for m in MODEL_ORDER if m in df.columns]
    remaining = [m for m in df.columns if m not in ordered]
    df = df.reindex(ordered + remaining, axis=1)
    df = df.sort_index()
    df = df.fillna(False)
    return df


def load_significance_table(dataset):
    """Load or build the boolean significance table for a dataset."""
    cfg = DATASET_CONFIGS[dataset]
    dir_au = join(cfg['working_dir'], 'outputs', cfg['au_folder'])
    csv_path = join(dir_au, cfg['sig_table_csv'])

    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path, index_col='RNA_family')
        # Convert string 'True'/'False' to bool
        df = df.apply(lambda col: col.astype(str).str.strip().eq('True'))
        ordered = [m for m in MODEL_ORDER if m in df.columns]
        remaining = [m for m in df.columns if m not in ordered]
        df = df.reindex(ordered + remaining, axis=1)
        print(f"Loaded significance table from {csv_path}")
        return df, dir_au

    print(f"No pre-computed CSV found. Parsing .pv files from {dir_au} ...")
    df = build_significance_table_from_pv(dir_au)
    df.to_csv(csv_path)
    print(f"Saved significance table to {csv_path}")
    return df, dir_au


def plot_significance_bar_chart(sig_table, dataset, output_path):
    """Plot bar chart of significant RNA family counts per model."""
    counts = sig_table.sum().astype(int)
    intersection = int(sig_table.all(axis=1).sum())

    labels = list(counts.index) + ['Intersection']
    values = list(counts.values) + [intersection]
    # Per-model bars in the RNA semantic colour; intersection bar in the
    # highlight colour to set it visually apart.
    colors = [COLORS['RNA']] * len(counts) + [COLORS['highlight']]

    # Width scales with bar count, capped at journal double-column.
    fig_width = min(FIG_DOUBLE, max(FIG_DOUBLE * 0.55, len(labels) * 0.45))
    fig, ax = plt.subplots(figsize=(fig_width, 3.0))
    bars = ax.bar(labels, values, color=colors,
                  edgecolor='white', linewidth=0.5)

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 2,
                str(val), ha='center', va='bottom')

    ax.set_ylabel('Number of significant RNA families')
    ax.set_title(f'AU-test significant RNA families per model — '
                 f'{dataset.upper()}-alignment dataset')
    ax.set_ylim(0, max(values) * 1.12)
    # Tidy axis: drop the box, keep only the bottom and left spines.
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis='x', rotation=45)
    for tick in ax.get_xticklabels():
        tick.set_horizontalalignment('right')
    plt.tight_layout()
    save_figure(fig, output_path)
    print(f"Figure saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Plot AU-test significant RNA families per model.')
    parser.add_argument('--dataset', choices=['seed', 'full'], required=True)
    args = parser.parse_args()

    sig_table, dir_au = load_significance_table(args.dataset)

    print(f"\nSignificant families per model:")
    for col in sig_table.columns:
        sig_families = sorted(sig_table.index[sig_table[col]])
        print(f"  {col}: {len(sig_families)} families")
        print(f"    {sig_families}")

    intersection_mask = sig_table.all(axis=1)
    intersection_families = sorted(sig_table.index[intersection_mask])
    print(f"\n  Intersection (all models): {len(intersection_families)} families")
    print(f"    {intersection_families}")
    print(f"  Total families: {len(sig_table)}")

    output_png = join(dir_au, f'{DATE_TAG}_AUtest_{args.dataset.upper()}.png')
    plot_significance_bar_chart(sig_table, args.dataset, output_png)


if __name__ == '__main__':
    main()
