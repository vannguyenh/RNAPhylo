"""
Batch U-test on normalized RF distances (single model, many RNA families).

For each RNA directory under DIR_RF:
  • Read DNA↔DNA       matrix: DNA_extra/<RNA>.rfdist           (10×10)
  • Read DNA↔RNA (Pi)  matrix: S**/<RNA>.rfdist                 (10×10)
  • Normalize each cell to nRF = RF / [2*(n-3)] 

  • Compute two-sided Mann–Whitney U-test
  • Record per-RNA medians for both groups
Outputs:
  1) Utest_long.csv        : Model, RNA, n_DNA, n_RNA, U, pvalue, p_bonf, flags
  2) Utest_wide.csv        : RNA × Model table of (raw) p-values
  3) Median_nRF_long.csv   : long table of per-RNA medians (for plotting)
Notes:
  • n = #taxa is taken from a tree file per RNA that model. We use <>/<RNA>/RAxML_bestTree.*
"""

import os
from os.path import join, isdir
from Bio import Phylo
import pandas as pd 

import subprocess
import logging
from datetime import datetime


LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'

MODEL = 'S6A'
# change to the right pathway of DIR_DNA for FULL or SEED data
DIR_WORKING = "/Users/u7875558/RNAPhylo/fullAlignment_S6A"
DIR_OUTPUTS = join(DIR_WORKING, "outputs")
DIR_TREES = join(DIR_OUTPUTS, 'inferred_trees')
DIR_RF     = join(DIR_OUTPUTS, "260204_Robinson_Foulds")
os.makedirs(DIR_RF, exist_ok=True)
DIR_RF_LOGS = join(DIR_WORKING, "logs", "RF_distance")
os.makedirs(DIR_RF_LOGS, exist_ok=True)

DIR_DNA=join(DIR_TREES, "DNA")
DIR_RNA=join(DIR_TREES, "S6A")

def check_branch_length(diroutput, rna):
    """
    Check out if the tree inferred under DNA model has any branch length larger than 1.
    If yes, there would be an issue with the dataset of this RNA family -- eliminate for the downstream analysis.
    diroutput: the folder containing DNA trees. == DIR_DNA
    """
    issue=False
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
        issue=True
    return issue

def produceCombinedTree(dir_output, rna):
    # 1) DNA set (baseline)
    dna_input = join(DIR_DNA, rna)
    if not isdir(dna_input):
        logging.warning(f"[S6A] Missing DNA folder for {rna}: {dna_input}")
        return False
    
    dir_RFrna = join(dir_output, rna) # RF-RNA folder
    os.makedirs(dir_RFrna, exist_ok=True)

    dna_out = join(dir_RFrna, f"{rna}.dna")
    with open(dna_out, "w") as out:
        files = sorted([f for f in os.listdir(dna_input) if f.startswith('RAxML_bestTree')])
        if not files:
            logging.warning(f"[S6A] No DNA trees for {rna}: {dna_input}")
        for fn in files:
            with open(join(DIR_DNA, rna, fn)) as fh:
                out.write(fh.read())

    # 2) RNA set
    rna_out = join(dir_RFrna, f"{rna}.rna")
    with open(rna_out, "w") as out:
        files = sorted([f for f in os.listdir(join(DIR_RNA, rna)) if f.startswith("RAxML_bestTree.")])
        for fn in files:
            with open(join(DIR_RNA, rna, fn), "r") as fh:
                out.write(fh.read())
        
def run_command(command):
    try:
        process=subprocess.Popen(command, shell=True)
        process.communicate()
    except Exception as e:
        logging.error(f"Command failed: {command} with error: {e}")

def computeRFdistance_iqtreecmd(dcombine_path, rna):
    dir_combine_rna=join(dcombine_path, rna)
    for f in os.listdir(dir_combine_rna):
        if f.endswith('dna'):
            DNATree = join(dir_combine_rna, f)
        elif f.endswith('rna'):
            RNATree = join(dir_combine_rna, f)

    logging.info(f"Compute the RF distance of {rna}.")
    command=f"bash /Users/u7875558/Documents/promotion/projects/projects_code/RNAPhylo/scripts/analysis_scripts/RobinsonFould/computeRFdistance.sh {DNATree} {RNATree} {dir_combine_rna} {rna}"
    run_command(command=command)

def main():
    dir_combined=join(DIR_RF, 'S6A')
    os.makedirs(dir_combined, exist_ok=True)

    log_filename = os.path.join(DIR_RF_LOGS, f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.normalisedRF.iqtree3.S6A.log")
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, format=LOG_FORMAT)
    logging.info(f"Running the code with the model S6A.")

    for rna in os.listdir(DIR_DNA):
        dir_inferred_rna = join(DIR_RNA, rna)

        if os.path.exists(dir_inferred_rna):
            dna_num_files = len(os.listdir(join(DIR_DNA, rna)))
            rna_num_files = len(os.listdir(join(DIR_RNA, rna)))
            if (dna_num_files != 50) & (rna_num_files != 50):
                logging.warning(f"{rna} does not have 50 files in the DNA or RNA inference folder.")
            elif check_branch_length(DIR_DNA, rna):
                logging.warning(f"{rna} has issue with branch length > 1 in the DNA inferred trees.")   
            else:
                logging.info(f"Working with {rna}.")
                path = join(dir_combined, rna)
                rfdist_file = 0
                if os.path.isdir(path):
                    for f in os.listdir(path):
                        if f.endswith('.rfdist'):
                            rfdist_file += 1
                if rfdist_file != 2:
                    produceCombinedTree(dir_combined, rna)
                    computeRFdistance_iqtreecmd(dir_combined, rna)
                else:
                    continue
        else:
            logging.warning(f"{rna} does not have RNA inference folder.")

if __name__=="__main__":
    main()