"""
Combine DNA and RNA (or DNA_extra) trees for each RNA family and compute RF distances.

Usage:
    python combineAllTrees_computeRF.py --dataset seed
    python combineAllTrees_computeRF.py --dataset full
"""

import argparse
import os
from os.path import join, isdir

from Bio import Phylo
import subprocess
import logging
from datetime import datetime

LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

DATASET_CONFIGS = {
    'seed': {
        'working_dir': '/Users/u7875558/RNAPhylo/seedAlignment_AllModels',
        'rf_folder': '260220_RF_distances',
    },
    'full': {
        'working_dir': '/Users/u7875558/RNAPhylo/fullAlignment',
        'rf_folder': '260220_RF_distances',
    },
}


def check_branch_length(diroutput, rna):
    """
    Check if any tree inferred under the DNA model has a branch length > 1.
    Returns True if an issue is found.
    """
    issue = False
    expected_files = [f"{i:02d}" for i in range(1, 11)]
    dirRNA = os.path.join(diroutput, rna)
    delete_files = []

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
        issue = True
    return issue


def produceCombinedTrees(model, dir_output, rna, dir_dna, dir_trees):
    """
    Write all best trees for an RNA family into two files:
      <rna>.dna  — trees from the DNA model
      <rna>.rna  — trees from the RNA model  (or <rna>.extra for DNA_extra)
    Returns True on success, False if any required folder is missing.
    """
    # 1) DNA baseline
    dna_input = join(dir_dna, rna)
    if not isdir(dna_input):
        logging.warning(f"[{model}] Missing DNA folder for {rna}: {dna_input}")
        return False

    dir_RFrna = join(dir_output, rna)
    os.makedirs(dir_RFrna, exist_ok=True)

    dna_out = join(dir_RFrna, f"{rna}.dna")
    with open(dna_out, "w") as out:
        for fn in sorted(f for f in os.listdir(dna_input) if f.startswith('RAxML_bestTree')):
            with open(join(dna_input, fn)) as fh:
                out.write(fh.read())

    # 2) RNA or DNA_extra
    if model.endswith('extra'):
        other_root = join(dir_trees, 'DNA_extra', rna)
        other_suffix = 'extra'
    else:
        other_root = join(dir_trees, model, rna)
        other_suffix = 'rna'

    if not isdir(other_root):
        logging.warning(f"[{model}] Missing folder for {rna}: {other_root}")
        return False

    other_out = join(dir_RFrna, f"{rna}.{other_suffix}")
    with open(other_out, "w") as out:
        for fn in sorted(f for f in os.listdir(other_root) if f.startswith("RAxML_bestTree.")):
            with open(join(other_root, fn)) as fh:
                out.write(fh.read())

    return True


def run_command(command):
    try:
        process = subprocess.Popen(command, shell=True)
        process.communicate()
    except Exception as e:
        logging.error(f"Command failed: {command} with error: {e}")


def computeRFdistance_iqtreecmd(dcombine_path, rna):
    dir_combine_rna = join(dcombine_path, rna)
    dna_tree = other_tree = None

    for f in os.listdir(dir_combine_rna):
        if f.endswith('.dna'):
            dna_tree = join(dir_combine_rna, f)
        elif f.endswith('.rna') or f.endswith('.extra'):
            other_tree = join(dir_combine_rna, f)

    if dna_tree is None or other_tree is None:
        logging.error(f"Missing tree files for {rna} in {dir_combine_rna}")
        return

    logging.info(f"Computing RF distances for {rna}.")
    rf_script = join(os.path.dirname(SCRIPT_DIR), 'computeRFdistance.sh')
    command = f"bash {rf_script} {dna_tree} {other_tree} {dir_combine_rna} {rna}"
    run_command(command)


def main():
    parser = argparse.ArgumentParser(description='Compute RF distances between DNA and RNA/DNA_extra trees.')
    parser.add_argument('--dataset', choices=['seed', 'full'], required=True,
                        help='Which alignment dataset to use: "seed" or "full".')
    args = parser.parse_args()

    cfg = DATASET_CONFIGS[args.dataset]
    DIR_WORKING = cfg['working_dir']
    DIR_OUTPUTS = join(DIR_WORKING, 'outputs')
    DIR_TREES = join(DIR_OUTPUTS, 'inferred_trees')
    DIR_DNA = join(DIR_TREES, 'DNA')
    DIR_RF = join(DIR_OUTPUTS, cfg['rf_folder'])
    os.makedirs(DIR_RF, exist_ok=True)

    DIR_RF_LOGS = join(DIR_WORKING, 'logs', 'RF_distance')
    os.makedirs(DIR_RF_LOGS, exist_ok=True)

    MODEL = input('Model: ')
    dir_RFmodel = join(DIR_RF, MODEL)
    os.makedirs(dir_RFmodel, exist_ok=True)

    log_filename = join(DIR_RF_LOGS, f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.{MODEL}.log")
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, format=LOG_FORMAT)
    logging.info(f"Running with dataset={args.dataset}, model={MODEL}.")

    for rna in os.listdir(DIR_DNA):
        dir_inferred_rna = join(DIR_TREES, MODEL, rna)

        if not os.path.exists(dir_inferred_rna):
            logging.warning(f"{MODEL} does not have a folder for {rna}.")
            continue

        dna_num_files = len(os.listdir(join(DIR_DNA, rna)))
        rna_num_files = len(os.listdir(dir_inferred_rna))

        if dna_num_files != 50 or rna_num_files != 50:
            logging.warning(f"{rna}: expected 50 files each; got DNA={dna_num_files}, {MODEL}={rna_num_files}.")
            continue

        if check_branch_length(DIR_DNA, rna):
            logging.warning(f"{rna}: branch length > 1 detected in DNA trees, skipping.")
            continue

        # Check whether RF distance has already been computed
        path = join(dir_RFmodel, rna)
        rfdist_count = sum(1 for f in os.listdir(path) if f.endswith('.rfdist')) if os.path.isdir(path) else 0

        if rfdist_count == 1:
            logging.info(f"{rna}: RF distance already computed, skipping.")
            continue

        logging.info(f"Working with {rna}.")
        produceCombinedTrees(MODEL, dir_RFmodel, rna, DIR_DNA, DIR_TREES)
        computeRFdistance_iqtreecmd(dir_RFmodel, rna)


if __name__ == "__main__":
    main()
