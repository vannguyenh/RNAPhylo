"""
Batch U-test on normalized RF distances (single model, many RNA families).

For each RNA directory under DIR_RF:
  • Read DNA↔DNA       matrix: <RNA>.raxml.rfdist         (10×10)
  • Read DNA↔RNA (Pi)  matrix: <RNA>.raxml.raxmlPi.rfdist (10×10)
  • Normalize each cell to nRF = RF / [2*(n-3)]  (auto-skips if already ≤1)
  • DNA↔DNA: use all 100 values (1st 10 DNA trees against 2nd 10 DNA trees) → 100 values
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
from os.path import join, isdir, basename
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.stats import mannwhitneyu
from sklearn.preprocessing import normalize as _ignore
from statsmodels.stats.multitest import multipletests

import subprocess
import re
import logging
from datetime import datetime

# ── PATHS (EDIT THESE) ─────────────────────────────────────────────────────────
LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'

# change to the right pathway of DIR_DNA for SEED data
DIR_WORKING = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels"
DIR_OUTPUTS = join(DIR_WORKING, "outputs")
DIR_TREES = join(DIR_OUTPUTS, 'inferred_trees')
DIR_RF      = join(DIR_OUTPUTS, "251121_Robinson_Foulds")   # contains <RNA> subfolders
os.makedirs(DIR_RF, exist_ok=True)
DIR_RF_LOGS = join(DIR_WORKING, "logs", "RF_distance")
os.makedirs(DIR_RF_LOGS, exist_ok=True)

DIR_DNA=join(DIR_TREES, 'DNA')
DIR_DNAextra=join(DIR_TREES, 'DNA_extra')
RF_ID_RE = re.compile(r'^RF\d{5}$')



# what we will concatenate
#METHOD_SUFFIX = {
#    "raxml":         "raxml",      # DNA set (baseline)
#    "raxmlP_iPseu":  "raxmlPi",    # RNA model set
#    "raxml_extra":   "raxml.extra" # extra DNA set (when MODEL == DNA_extra)
#}

def check_branch_length(diroutput, rna):
    """
    Check out if the tree inferred under DNA model has any branch length larger than 1.
    If yes, there would be an issue with the dataset of this RNA family -- eliminate for the downstream analysis.
    diroutput: the folder containing DNA trees. == DIR_DNA
    """
    expected_files = [f"{i:02d}" for i in range(1, 11)]
    dirRNA = os.path.join(diroutput, rna)
    delete_files = list()

    for file_name in os.listdir(dirRNA):
        for seed in expected_files:
            if file_name.startswith('RAxML_bestTree') and file_name.endswith(seed):
                tree_path = os.path.join(dirRNA, file_name)
                tree = Phylo.read(tree_path, 'newick')
                for clade in tree.find_clades():
                    if clade.branch_length and clade.branch_length > 1:
                        delete_files.append(file_name)
                        logging.info(f"{file_name} has branch length > 1.")

    if len(delete_files) > 0:
        return rna

def extractAnalysedRNAs(log_file):
    # This function produces two sets of accepted RNAs -- RNAs containing pseudoknots and RNAs containing no pseudoknots    
    # extract accepted RNAs which do not have any branch length > 1
    rnas = os.listdir(DIR_DNA)
    
    brlength_larger1=list()
    less4_pat = re.compile(r'Number of sequences of (RF\d{5}) is less than 4\.')
    unpaired_pat = re.compile(r'The secondary structure of (RF\d{5}) has only unpaired bases\.')  # note: "unpaired"
    less4, unpaired = set(), set()
    
    accepted_rnas=list()
    
    with open(log_file, 'r') as lf: #log file used for this full_S6A.log
        for line in lf:
            m1 = less4_pat.search(line)
            if m1: less4.add(m1.group(1))
            m2 = unpaired_pat.search(line)
            if m2: unpaired.add(m2.group(1))
    
    for rna in rnas:
        # consider only the outputs from using DNA
        if check_branch_length(DIR_DNA, rna) == rna:
            brlength_larger1.append(rna)
        elif (rna not in less4 and rna not in unpaired):
            accepted_rnas.append(rna)

    unaccepted_rnas = less4 | unpaired | set(brlength_larger1)

    logging.info(f'{len(accepted_rnas)} can be used for the downstream analysis.')
    logging.info(f'{len(unaccepted_rnas)} cannot be used for the downstream analysis.')
    return accepted_rnas, unaccepted_rnas

def produceCombinedTrees(model, dir_output, rna):
    # 1) DNA set (baseline)
    dna_input = join(DIR_DNA, rna)
    if not isdir(dna_input):
        logging.warning(f"[{model}] Missing DNA folder for {rna}: {dna_input}")
        return False
    
    dir_RFrna = join(dir_output, rna) # RF-RNA folder
    os.makedirs(dir_RFrna, exist_ok=True)

    dna_out = join(dir_RFrna, f"{rna}.dna")
    with open(dna_out, "w") as out:
        files = sorted([f for f in os.listdir(dna_input) if f.startswith('RAxML_bestTree')])
        if not files:
            logging.warning(f"[{model}] No DNA trees for {rna}: {dna_input}")
        for fn in files:
            with open(join(DIR_DNA, rna, fn)) as fh:
                out.write(fh.read())

    if model.endswith('extra'):
        other_root = join(DIR_DNAextra, rna)
        other_suffix = 'dna.extra'
    else:
        other_root = join(DIR_TREES, model, rna)
        other_suffix = 'rna'
    if not isdir(other_root):
        logging.warning(f"[{model}] Missing DNAextra/RNA folder for {rna}: {other_root}")
        return False
    
    other_out = join(dir_RFrna, f"{rna}.{other_suffix}")
    with open(other_out, "w") as out:
        files = sorted([f for f in os.listdir(other_root) if f.startswith("RAxML_bestTree.")])
        for fn in files:
            with open(join(other_root, fn), "r") as fh:
                out.write(fh.read())

def run_command(command):
    try:
        process = subprocess.Popen(command, shell=True)
        process.communicate()
    except Exception as e:
        logging.error(f"Command failed: {command} with error: {e}")

def computeRFdistance_iqtreecmd(dcombine_path, rna):
    dir_combine_rna = os.path.join(dcombine_path, rna)
    for f in os.listdir(dir_combine_rna):
        if f.endswith("dna"):
            DNATree = os.path.join(dir_combine_rna, f)
        elif f.endswith("rna"):
            otherTree = os.path.join(dir_combine_rna, f)
        elif f.endswith("extra"):
            otherTree = os.path.join(dir_combine_rna, f)

    logging.info(f"Compute the RF distances of {rna}.")
    command=f"bash /Users/u7875558/Documents/promotion/projects/projects_code/RNAPhylo/scripts/analysis_scripts/RobinsonFould/computeRFdistance.sh {DNATree} {otherTree} {dir_combine_rna} {rna}"
    run_command(command)

def main():
    MODEL = 'S7F'
    dir_combined = join(DIR_RF, MODEL)
    os.makedirs(dir_combined, exist_ok=True)

    log_filename = join(DIR_RF_LOGS, f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.{MODEL}.log")
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, format=LOG_FORMAT)
    logging.info(f"Running the code with the model {MODEL}.")

    for rna in os.listdir(DIR_DNA):
        if os.path.exists(join(DIR_TREES, MODEL, rna)):
            dna_num_files = len(os.listdir(join(DIR_DNA, rna)))
            rna_num_files = len(os.listdir(join(DIR_TREES, MODEL, rna)))
            if (dna_num_files != 50) & (rna_num_files != 50):
                logging.warning(f"{rna} does not have enough inferred trees in either DNA or {MODEL}")
            else:
                logging.info(f"Working with {rna}.")
                path = join(dir_combined, rna)
                rfdist_file = 0
                if os.path.isdir(path):
                    for f in os.listdir(path):
                        if f.endswith('.rfdist'):
                            rfdist_file += 1

                if rfdist_file != 1: # 3 while using only ramxl and raxmlP_iPseu, 6 while using raxml, raxmlP_iPseu and raxmlP_wPseu
                    #run_command(f"rm -r {path}")
                    produceCombinedTrees(MODEL, dir_combined, rna)
                    computeRFdistance_iqtreecmd(dir_combined, rna)
                else:
                    continue
        else:
            logging.warning(f"{MODEL} does not have {rna}")

if __name__=="__main__":
    main()
