#!/usr/bin/env python3
"""
AU_pvalues_table.py

Extract AU p-values for the DNA tree (ignore pseudoknots) across all S6*/S16* models,
and assemble a wide table (RNAs × models) of those raw p-values.
"""

import os
from os import listdir
from os.path import join, isdir
import subprocess
import pandas as pd
import numpy as np

# ─── USER PARAMETERS ───────────────────────────────────────────────────────────

DIR_WORKING = '/Users/u7875558/RNAPhylo/seedAlignment_AllModels'
DIR_OUTPUTS = join(DIR_WORKING, 'outputs')
DIR_AU = join(DIR_OUTPUTS, '260208_AU_Test_RAxML_RNAmodels')
DATE = '260210'
SUFFIX_IGNORE = '_consel'  # suffix before .pv


# ─── UTILITY FUNCTIONS ────────────────────────────────────────────────────────

def parse_consel_output(pv_path):
    """
    Run catpv on a .pv file and return [p_DNA, p_RNA] AU p-values.
    """
    p_DNA = None
    p_RNA = None
    proc = subprocess.run(
        ['/Users/u7875558/tools/consel/bin/catpv', pv_path],
        stdout=subprocess.PIPE, text=True
    )
    for line in proc.stdout.splitlines():
        line = line.strip()
        if line.startswith('#') and line[1:].strip()[0].isdigit():
            cols = line[1:].split()
            if cols[1] == '1':  # DNA is item '1'
                p_DNA = float(cols[3])
            elif cols[1] == '2':  # RNA is item '2'
                p_RNA = float(cols[3])
    return [p_DNA, p_RNA]

# ─── MAIN PROCESS ─────────────────────────────────────────────────────────────

def main():
    # Debug: Check directory
    print(f"Looking in: {DIR_AU}")
    print(f"Directory exists: {os.path.exists(DIR_AU)}")

    # Discover S6*/S16*/S7* models
    models = sorted(
        m for m in listdir(DIR_AU)
        if isdir(join(DIR_AU, m)) and (m.startswith('S') or m.startswith('extra'))
    )
    print(f"Models found: {models}")

    records = []

    for model in models:
        consel_model = join(DIR_AU, model)
        for rna in listdir(consel_model):
            # Check both patterns in case of naming inconsistency
            pv_file = join(consel_model, rna, f"{rna}{SUFFIX_IGNORE}.pv")
            if not os.path.exists(pv_file):
                pv_file = join(consel_model, rna, f"{rna}_{SUFFIX_IGNORE}.pv")
            if not os.path.exists(pv_file):
                continue
            p_raw = parse_consel_output(pv_file)
            if len(p_raw) < 2:
                print(f"Warning: Could not parse both p-values from {pv_file}")
                continue
            records.append({
                'Model': model,
                'RNA_family': rna,
                'p_DNA_raw': p_raw[0],
                'p_RNA_raw': p_raw[1]
            })

    # Build DataFrame
    df = pd.DataFrame(records)
    print(f"\nTotal records: {len(df)}")
    print(df.head())

    if df.empty:
        print("ERROR: No records found!")
        return

    # ─── KEY ANALYSIS: Find significant differences ───────────────────────────

    # Minimum p-value (either tree rejected = significant difference)
    df['p_min'] = df[['p_DNA_raw', 'p_RNA_raw']].min(axis=1)

    # Flag significant differences (p < 0.05)
    df['sig_different'] = df['p_min'] < 0.05

    # Which tree was rejected?
    df['rejected_tree'] = np.where(
        df['p_DNA_raw'] < df['p_RNA_raw'], 'DNA', 'RNA'
    )

    # ─── Summary statistics ───────────────────────────────────────────────────

    print("\n" + "=" * 60)
    print("SUMMARY: Significant differences (p < 0.05)")
    print("=" * 60)

    for model in df['Model'].unique():
        model_df = df[df['Model'] == model]
        n_sig = model_df['sig_different'].sum()
        n_total = len(model_df)
        pct = 100 * n_sig / n_total
        print(f"{model}: {n_sig}/{n_total} ({pct:.1f}%) RNA families show significant difference")

    # ─── Save outputs ─────────────────────────────────────────────────────────

    # Full results with all columns
    full_csv = join(DIR_AU, '260208_AU_pValues_both_SEED.csv')
    df.to_csv(full_csv, index=False)

    # Pivot table of minimum p-values
    df_pmin_table = df.pivot(index='RNA_family', columns='Model', values='p_min').sort_index()
    pmin_csv = join(DIR_AU, f'{DATE}_AU_pMin_table.csv')
    df_pmin_table.to_csv(pmin_csv)

    # Pivot table of significance flags
    df_sig_table = df.pivot(index='RNA_family', columns='Model', values='sig_different').sort_index()
    sig_csv = join(DIR_AU, f'{DATE}_AU_significant_table.csv')
    df_sig_table.to_csv(sig_csv)

    # List of significant RNA families per model
    sig_families = df[df['sig_different']][['Model', 'RNA_family', 'p_min', 'rejected_tree']]
    sig_families_csv = join(DIR_AU, f'{DATE}_AU_significant_families.csv')
    sig_families.to_csv(sig_families_csv, index=False)

    print(f"\nSaved outputs:")
    print(f"  Full results: {full_csv}")
    print(f"  p_min table: {pmin_csv}")
    print(f"  Significance table: {sig_csv}")
    print(f"  Significant families: {sig_families_csv}")


if __name__ == '__main__':
    main()