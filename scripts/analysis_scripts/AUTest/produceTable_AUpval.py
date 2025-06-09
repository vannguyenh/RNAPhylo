#!/usr/bin/env python3
"""
AU_pvalues_table.py

Extract AU p-values for the DNA tree (ignore pseudoknots) across all S6*/S16* models,
and assemble a wide table (RNAs × models) of those raw p-values.
"""

import os
from os.path import join
import subprocess
import numpy as np
import pandas as pd
from Bio import Phylo

# ─── USER PARAMETERS ───────────────────────────────────────────────────────────

DIR_WORKING     = '/Users/u7875558/Documents/PhD/RNAPhylo/allModels_SEED'
DIR_OUTPUTS     = join(DIR_WORKING, 'outputs')
DIR_DNA         = join(DIR_OUTPUTS, 'DNAtrees')
# CONS EL outputs: outputs/AU_Test/<model>/ignore_pseudo/<RNA>/<RNA>_ipseu_consel.pv
DIR_CONSEL_BASE = join(DIR_OUTPUTS, 'AU_Test')

EXPECTED_SEEDS = {f"{i:02d}" for i in range(1, 11)}
SUFFIX_IGNORE = 'ipseu_consel'  # suffix before .pv for ignore_pseudo

# RNAs to exclude
ISSUE_RNAS =  ['RF00207', 'RF00976', 'RF01047', 'RF01338', 'RF01380', 'RF03623', 'RF03760', 'RF03969']

# ─── UTILITY FUNCTIONS ────────────────────────────────────────────────────────

def check_branch_length(dna_dir, rna):
    """
    Return True if any replicate tree under dna_dir/rna has branch_length > 1.
    """
    folder = join(dna_dir, rna)
    for fn in os.listdir(folder):
        if fn.startswith('RAxML_bestTree') and fn[-2:] in EXPECTED_SEEDS:
            tree = Phylo.read(join(folder, fn), 'newick')
            if any(cl.branch_length and cl.branch_length > 1 for cl in tree.find_clades()):
                return True
    return False


def get_accepted_rnas(dna_dir):
    """
    List RNAs in dna_dir with no overlong branches and not in ISSUE_RNAS.
    """
    accepted = []
    for rna in os.listdir(dna_dir):
        path = join(dna_dir, rna)
        if not os.path.isdir(path) or rna in ISSUE_RNAS:
            continue
        if not check_branch_length(dna_dir, rna):
            accepted.append(rna)
    return sorted(accepted)


def parse_consel_output(pv_path):
    """
    Run catpv on a .pv file and return the AU p-value for the DNA tree.
    """
    proc = subprocess.run(
        ['/Users/u7875558/tools/consel/bin/catpv', pv_path],
        stdout=subprocess.PIPE, text=True
    )
    for line in proc.stdout.splitlines():
        line = line.strip()
        if line.startswith('#') and line[1:].strip()[0].isdigit():
            cols = line[1:].split()
            if cols[1] == '1':  # DNA is item '1'
                return float(cols[3])
    return np.nan

# ─── MAIN PROCESS ─────────────────────────────────────────────────────────────

# Discover S6*/S16*/S7* models
models = sorted(
    m for m in os.listdir(DIR_CONSEL_BASE)
    if os.path.isdir(join(DIR_CONSEL_BASE, m)) and ((m.startswith('S') or (m.startswith('extra'))))
)

accepted_rnas = get_accepted_rnas(DIR_DNA)
records = []

for model in models:
    consel_ignore = join(DIR_CONSEL_BASE, model, 'ignore_pseudo')
    for rna in accepted_rnas:
        pv_file = join(consel_ignore, rna, f"{rna}_{SUFFIX_IGNORE}.pv")
        if not os.path.exists(pv_file):
            continue
        p_raw = parse_consel_output(pv_file)
        records.append({'Model': model, 'RNA': rna, 'p_raw': p_raw})

# Build DataFrame and pivot to wide format
df = pd.DataFrame(records)
df_table = df.pivot(index='RNA', columns='Model', values='p_raw').sort_index()

# Output
def main():
    print(df_table.head())
    out_csv = join(DIR_CONSEL_BASE, 'AU_pvalues_ignore_pseudo_all_models.csv')
    df_table.to_csv(out_csv)
    print(f"Saved p-values from AU test of all models to: {out_csv}")

if __name__ == '__main__':
    main()