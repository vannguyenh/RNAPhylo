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

# ─── USER PARAMETERS ───────────────────────────────────────────────────────────

DIR_WORKING     = '/Users/u7875558/RNAPhylo/seedAlignment_AllModels'
#DIR_WORKING = '/Users/u7875558/RNAPhylo/fullAlignment_S6A'
DIR_OUTPUTS     = join(DIR_WORKING, 'outputs')
# CONSEL outputs: outputs/AU_Test_RAxMLNG/<model>/ignore_pseudo/<RNA>/<RNA>_ipseu_consel.pv
DIR_AU = join(DIR_OUTPUTS, '260208_AU_Test_RAxML_RNAmodels')

EXPECTED_SEEDS = {f"{i:02d}" for i in range(1, 11)}
SUFFIX_IGNORE = 'consel'  # suffix before .pv for ignore_pseudo

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
        if isdir(join(DIR_AU, m)) and (m.startswith('S') or (m.startswith('extra')))
    )

    records = []

    for model in models:
        consel_model= join(DIR_AU, model)
        for rna in listdir(consel_model):
            pv_file = join(consel_model, rna, f"{rna}_{SUFFIX_IGNORE}.pv")
            if not os.path.exists(pv_file):
                continue
            p_raw = parse_consel_output(pv_file)
            records.append({'Model': model, 'RNA_family': rna, 
                        'p_DNA_raw': p_raw[0], 'p_RNA_raw': p_raw[1]})

    # Build DataFrame and pivot to wide format
    df = pd.DataFrame(records)
    print(df.head())
    full_csv = join(DIR_AU, '260210_AU_pValues_both_SEED.csv')
    df.to_csv(full_csv, sep='\t')

    df_table = df.pivot(index='RNA_family', columns='Model', values='p_RNA_raw').sort_index()
    print(df_table.head())
    rna_csv = join(DIR_AU, '260208_AUpvalues_RNApVal_SEED.csv' )
    df_table.to_csv(rna_csv)
    print(f"Saved p-values from AU test of all models to: {full_csv}, {rna_csv}")

    df_table = df.pivot(index='RNA_family', columns='Model', values='p_DNA_raw').sort_index()
    print(df_table.head())
    dna_csv = join(DIR_AU, '260208_AUpvalues_DNApVal_SEED.csv')
    df_table.to_csv(dna_csv)
    print(f"Saved p-values from AU test of all models to: {full_csv}, {dna_csv}")

if __name__ == '__main__':
    main()