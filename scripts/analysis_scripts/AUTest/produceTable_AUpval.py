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
import numpy as np
import pandas as pd
from Bio import Phylo

# ─── USER PARAMETERS ───────────────────────────────────────────────────────────

#DIR_WORKING     = '/Users/u7875558/Documents/PhD/RNAPhylo/allModels_SEED'
DIR_WORKING = '/Users/u7875558/RNAPhylo/fullAlignment_S6A'
DIR_OUTPUTS     = join(DIR_WORKING, 'outputs')
DIR_DNA         = join(DIR_OUTPUTS, 'DNA')
# CONSEL outputs: outputs/AU_Test_RAxMLNG/<model>/ignore_pseudo/<RNA>/<RNA>_ipseu_consel.pv
DIR_AU = join(DIR_OUTPUTS, 'AU_Test_RAxMLNG')

EXPECTED_SEEDS = {f"{i:02d}" for i in range(1, 11)}
SUFFIX_IGNORE = 'ipseu_consel'  # suffix before .pv for ignore_pseudo

# RNAs to exclude
#ISSUE_RNAS =  ['RF00207', 'RF00976', 'RF01047', 'RF01338', 'RF01380', 'RF03623', 'RF03760', 'RF03969']

# ─── UTILITY FUNCTIONS ────────────────────────────────────────────────────────

def parse_consel_output(pv_path):
    """
    Run catpv on a .pv file and return the AU p-value for the DNA tree.
    """
    result = list()
    proc = subprocess.run(
        ['/Users/u7875558/tools/consel/bin/catpv', pv_path],
        stdout=subprocess.PIPE, text=True
    )
    for line in proc.stdout.splitlines():
        line = line.strip()
        if line.startswith('#') and line[1:].strip()[0].isdigit():
            cols = line[1:].split()
            if cols[1] == '1':  # DNA is item '1'
                result.append(float(cols[3]))
            elif cols[1] == '2': # RNA is item '2'
                result.append(float(cols[3]))
    return result

# ─── MAIN PROCESS ─────────────────────────────────────────────────────────────

# Output
def main():
    # Discover S6*/S16*/S7* models
    models = sorted(
        m for m in listdir(DIR_AU)
        if isdir(join(DIR_AU, m)) and ((m.startswith('S') or (m.startswith('extra'))))
    )

    records = []

    for model in models:
        consel_ignore = join(DIR_AU, model, 'ignore_pseudo')
        for rna in listdir(consel_ignore):
            pv_file = join(consel_ignore, rna, f"{rna}_{SUFFIX_IGNORE}.pv")
            if not os.path.exists(pv_file):
                continue
            p_raw = parse_consel_output(pv_file)
            records.append({'Model': model, 'RNA_family': rna, 
                        'p_DNA_raw': p_raw[0], 'p_RNA_raw': p_raw[1]})

    # Build DataFrame and pivot to wide format
    df = pd.DataFrame(records)
    print(df.head())
    full_csv = join(DIR_AU, 'AU_pValues_FULL_full.csv')
    df.to_csv(full_csv, sep='\t')

    df_table = df.pivot(index='RNA_family', columns='Model', values='p_RNA_raw').sort_index()
    print(df_table.head())
    out_csv = join(DIR_AU, 'AU_pvalues_FULL_RNApVal.csv' )
    df_table.to_csv(out_csv)
    print(f"Saved p-values from AU test of all models to: {full_csv}, {out_csv}")

if __name__ == '__main__':
    main()