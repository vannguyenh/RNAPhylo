import argparse
import os
from os.path import join, isfile

from Bio import Phylo
import logging
from datetime import datetime
import subprocess

LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'

DATASET_CONFIGS = {
    'seed': {
        'working_dir': '/Users/u7875558/RNAPhylo/seedAlignment_AllModels',
        'ss_dir': 'ss_files',
        'log_tag': '260208_AU_test_RNAmodels',
    },
    'full': {
        'working_dir': '/Users/u7875558/RNAPhylo/fullAlignment',
        'ss_dir': 'ss_files',
        'log_tag': '260208_AU_test_RNAmodels',
    },
}

SCRIPT_DIR = '/Users/u7875558/Documents/promotion/projects/projects_code/RNAPhylo/scripts/analysis_scripts/AUTest'


def get_fasta_file(rna, dir_inputs):
    """Return the FASTA path for an RNA family.
    - First check subsample (for families >500 taxa)
    - Fall back to fasta (for families ≤500 taxa)
    """
    subsamp = join(dir_inputs, 'subsample', f'{rna}.subsamp.fa')
    if isfile(subsamp):
        return subsamp

    fasta = join(dir_inputs, 'fasta', f'{rna}.nodup.fa')
    if isfile(fasta):
        return fasta

    return None


def check_branch_length(diroutput, rna):
    """
    Check if any tree inferred under the DNA model has a branch length > 1.
    Returns True if an issue is found, False otherwise.
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


def check_inferred_tree(dir_rna):
    if not os.path.isdir(dir_rna):
        logging.warning(f"{dir_rna} does not exist.")
        return False
    if len(os.listdir(dir_rna)) != 50:
        logging.warning(f"{dir_rna} does not have 10 trees.")
        return False
    return True


def extract_highestloglh(dir_path):
    lh = {}
    seeds = [f"{i:02d}" for i in range(1, 11)]

    for file in os.listdir(dir_path):
        for seed in seeds:
            if file.startswith('RAxML_log') and file.endswith(seed):
                with open(os.path.join(dir_path, file), 'r') as f:
                    final_result = f.readlines()[-1].strip()
                    best_value = final_result.split()[-1]
                    lh[seed] = float(best_value)
    return sorted(lh.items(), key=lambda x: x[1], reverse=True)[0]


def extract_highestLH_2Trees_rna_dna2(dir_output, rna, dir_dna):
    dna_path = join(dir_dna, rna)
    rna_path = join(dir_output, rna)
    if check_inferred_tree(dna_path) and check_inferred_tree(rna_path):
        return {
            'DNA': extract_highestloglh(dna_path),
            'RNA_otherDNA': extract_highestloglh(rna_path),
        }
    return None


def combineTreeFiles(dir_combined, rna, tree1, tree2):
    combined_TreeFile = join(dir_combined, f"{rna}.highestLH.trees")
    with open(combined_TreeFile, 'w') as f_combined:
        for tree_file in [tree1, tree2]:
            with open(tree_file, 'r') as f:
                f_combined.write(f.read())


def run_bash(command):
    process = subprocess.Popen(command, shell=True)
    process.communicate()
    if process.returncode != 0:
        logging.error(f"Bash command failed: {command}")
    return process.returncode


def run_command(command):
    process = subprocess.Popen(command, shell=True)
    process.communicate()


def runningCONSEL(rna, fasta_file, ss_file, persite_suffix, dir_AUrna, prefix_consel, model):
    combineTree = join(dir_AUrna, f'{rna}.highestLH.trees')
    output_persite = join(dir_AUrna, f'RAxML_perSiteLLs.{persite_suffix}')
    bash_file = join(SCRIPT_DIR, 'consel_RNA.sh')
    consel_command = f"bash {bash_file} {fasta_file} {combineTree} {ss_file} {persite_suffix} {dir_AUrna} {output_persite} {prefix_consel} {model}"
    run_bash(consel_command)


def has_sitelh(dirpath: str) -> bool:
    """Return True if any file in dirpath starts with RAxML_perSiteLLs."""
    try:
        return any(fn.startswith('RAxML_perSiteLLs') for fn in os.listdir(dirpath))
    except FileNotFoundError:
        return False


def main():
    parser = argparse.ArgumentParser(description='Run AU test via CONSEL for RNA phylogenetics.')
    parser.add_argument('--dataset', choices=['seed', 'full'], required=True,
                        help='Which alignment dataset to use: "seed" or "full".')
    parser.add_argument('--model', type=str, default=None,
                        help='RNA model (e.g., S6A). If omitted, will prompt interactively.')
    parser.add_argument('--skip-existing', action='store_true',
                        help='Skip families that already have .pv files.')
    args = parser.parse_args()

    cfg = DATASET_CONFIGS[args.dataset]
    DIR_WORKING = cfg['working_dir']
    DIR_OUTPUTS = join(DIR_WORKING, 'outputs')
    DIR_TREES = join(DIR_OUTPUTS, 'inferred_trees')
    DIR_DNA = join(DIR_TREES, 'DNA')
    DIR_INPUTS = join(DIR_WORKING, 'inputs')
    DIR_SS = join(DIR_INPUTS, cfg['ss_dir'])
    DIR_AU_LOGS = join(DIR_WORKING, 'logs', cfg['log_tag'])
    os.makedirs(DIR_AU_LOGS, exist_ok=True)

    MODEL = args.model if args.model else input('Model: ')
    DIR_AU = join(DIR_OUTPUTS, '260208_AU_Test_RAxML_RNAmodels', MODEL)
    os.makedirs(DIR_AU, exist_ok=True)

    log_filename = join(DIR_AU_LOGS, f"CONSEL_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.{MODEL}.log")
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, format=LOG_FORMAT)

    # For seed: iterate all inferred DNA trees; for full: iterate existing AU output folders
    rna_list = os.listdir(DIR_DNA if args.dataset == 'seed' else DIR_AU)

    skipped_existing = 0
    for rna in rna_list:
        # Skip families that already have .pv output
        if args.skip_existing:
            pv_path = join(DIR_AU, rna, f'{rna}_consel.pv')
            if isfile(pv_path):
                skipped_existing += 1
                continue

        # Resolve FASTA file (check subsample first for families >500 taxa)
        fasta_file = get_fasta_file(rna, DIR_INPUTS)
        if fasta_file is None:
            logging.warning(f"No FASTA file found for {rna}, skipping.")
            continue

        bestTrees = extract_highestLH_2Trees_rna_dna2(join(DIR_TREES, MODEL), rna, DIR_DNA)
        if bestTrees is None:
            logging.warning(f"{rna} has either under DNA or RNA model no tree inference.")
            continue

        dnaTree = join(DIR_DNA, rna, f"RAxML_bestTree.{rna}.{bestTrees['DNA'][0]}")
        rnaTree = join(DIR_TREES, MODEL, rna, f"RAxML_bestTree.{rna}.{bestTrees['RNA_otherDNA'][0]}")

        persite_path = join(DIR_AU, rna)
        os.makedirs(persite_path, exist_ok=True)

        combineTreeFiles(persite_path, rna, dnaTree, rnaTree)

        persite_suffix = f'{rna}.sitelh'
        prefix_consel = join(persite_path, f'{rna}_consel')
        ss_file = join(DIR_SS, f'{rna}.ss')

        try:
            runningCONSEL(rna, fasta_file, ss_file, persite_suffix, persite_path, prefix_consel, MODEL)
            logging.info(f"CONSEL is running with {rna}.")
        except Exception as e:
            logging.error(f"Error with {rna}: {e}")

        if has_sitelh(persite_path):
            logging.info(f'{rna} ran with CONSEL (found *.sitelh).')
        else:
            # Retry with reduced inputs
            run_command(f"rm -rf {persite_path}")
            os.makedirs(persite_path, exist_ok=True)
            combineTreeFiles(persite_path, rna, dnaTree, rnaTree)

            fasta_red_path = join(DIR_INPUTS, 'fasta', f'{rna}.nodup.fa.reduced')
            fasta_file_red = fasta_red_path if isfile(fasta_red_path) else None
            if fasta_file_red is None:
                logging.warning(f'No reduced FASTA found for {rna}; skipping reduced attempt.')

            ss_red_path = join(DIR_SS, f'{rna}.ss.reduced')
            ss_file_red = ss_red_path if isfile(ss_red_path) else None
            if ss_file_red is None:
                logging.warning(f'No reduced SS found for {rna}; skipping reduced attempt.')

            if fasta_file_red is not None:
                runningCONSEL(rna, fasta_file_red, ss_file_red, persite_suffix, persite_path, prefix_consel, MODEL)
                if has_sitelh(persite_path):
                    logging.info(f"{rna}: reduced attempt produced *.sitelh.")
                else:
                    logging.error(f"{rna}: reduced attempt still did not produce *.sitelh.")


    if args.skip_existing:
        logging.info(f"Skipped {skipped_existing} families with existing .pv files.")
        print(f"Skipped {skipped_existing} families with existing .pv files.")


if __name__ == '__main__':
    main()
