import os
from os.path import join, isfile

from Bio import Phylo
import logging
from datetime import datetime
import time
import subprocess
import re

LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'

DIR_WORKING = '/Users/u7875558/RNAPhylo/seedAlignment_AllModels'
DIR_OUTPUTS = join(DIR_WORKING, 'outputs')
DIR_TREES = join(DIR_OUTPUTS, 'inferred_trees')
DIR_DNA = join(DIR_TREES, 'DNA')
DIR_AU_LOGS = join(DIR_WORKING, 'logs', '260208_AU_test_RNAmodels')
os.makedirs(DIR_AU_LOGS, exist_ok=True)

DIR_INPUTS = join(DIR_WORKING, 'inputs')
DIR_FASTA = join(DIR_INPUTS, 'fasta_files')
DIR_SS = join(DIR_INPUTS, 'ss_files')

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
    
def check_inferred_tree(dir_rna):
    if len(os.listdir(dir_rna)) != 50:
        logging.warning(f"{dir_rna} does not have 10 trees.")
        return False
    else:
        return True
        
def extract_highestloglh(dir_path):
    lh = dict()
    files = os.listdir(dir_path)
    seeds = [f"{i:02d}" for i in range(1, 11)]

    for file in os.listdir(dir_path):
        for seed in seeds:
            if file.startswith('RAxML_log') and file.endswith(seed):
                with open(os.path.join(dir_path, file), 'r') as f:
                    final_result = f.readlines()[-1].strip()
                    best_value = final_result.split()[-1]
                    lh[seed] = float(best_value)
    return sorted(lh.items(), key=lambda x: x[1], reverse=True)[0]
    
def extract_highestLH_2Trees_ipseudo_dnaextra(dir_output, rna):
    dna_path = join(DIR_DNA, rna)
    rna_path = join(dir_output, rna)
    if check_inferred_tree(dna_path) and check_inferred_tree(rna_path):
        bestTreesLH = {'DNA': extract_highestloglh(dna_path),
                       'RNA_otherDNA': extract_highestloglh(rna_path)}
        return bestTreesLH
    else:
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
    # Check process.returncode 
    if process.returncode != 0:
        logging.error(f"Bash command failed: {command}")
    return process.returncode

def runningCONSEL(rna, fasta_file, ss_file, persite_suffix, dir_AUrna, prefix_consel, model):
    combineTree = join(dir_AUrna, f'{rna}.highestLH.trees')
    output_persite = join(dir_AUrna, f'RAxML_perSiteLLs.{persite_suffix}')    
    bash_file = ' /Users/u7875558/Documents/promotion/projects/projects_code/RNAPhylo/scripts/analysis_scripts/AUTest/consel_RNA.sh'
    consel_command = f"bash {bash_file} {fasta_file} {combineTree} {ss_file} {persite_suffix} {dir_AUrna} {output_persite} {prefix_consel} {model}"
    run_bash(consel_command)

def run_command(command):
    process = subprocess.Popen(command, shell=True)
    process.communicate()

def has_sitelh(dirpath: str) -> bool:
    """Return True if any file in dirpath ends with *sitelh (case-insensitive)."""
    try:
        return any(fn.startswith('RAxML_perSiteLLs') for fn in os.listdir(dirpath))
    except FileNotFoundError:
        return False
    
def main():
    MODEL = input('Model: ')
    DIR_AU = join(DIR_OUTPUTS, '260208_AU_Test_RAxML_RNAmodels', MODEL)
    os.makedirs(DIR_AU, exist_ok=True)

    #DIR_COMBINE = join(DIR_AU, 'combine_2trees_highestLH')
    #os.makedirs(DIR_COMBINE, exist_ok=True)

    log_filename = join(DIR_AU_LOGS, f"CONSEL_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.{MODEL}.log")
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
   
    for rna in os.listdir(DIR_DNA):
        bestTrees = extract_highestLH_2Trees_ipseudo_dnaextra(join(DIR_TREES, MODEL), rna)
        if bestTrees is not None:
            dnaTree = join(DIR_DNA, rna, f"RAxML_bestTree.{rna}.{bestTrees['DNA'][0]}")
            rnaTree = join(DIR_TREES, MODEL, rna, 
                           f"RAxML_bestTree.{rna}.{bestTrees['RNA_otherDNA'][0]}")
            
            # create the path of each RNA in the folder of AU test
            persite_path = join(DIR_AU, rna)
            os.makedirs(persite_path, exist_ok=True)
            
            # create a file containing DNA and RNA tree with the best log-likelihood for each RNA family, which would be used as input for CONSEL
            combineTreeFiles(persite_path, rna, dnaTree, rnaTree)

            persite_suffix = f'{rna}.sitelh'
            prefix_consel = join(persite_path, f'{rna}_consel')

            ss_file = join(DIR_SS, f'{rna}.ss')
            fasta_file = join(DIR_FASTA, f'{rna}.nodup.fa')
                
            try:
                runningCONSEL(rna, fasta_file, ss_file, persite_suffix, persite_path, prefix_consel, MODEL)
                logging.info(f"CONSEL is running with {rna}.")
            except Exception as e:
                logging.error(f"Error with {rna}: {e}")

            # If we already have a *.sitelh, move on
            if has_sitelh(persite_path):
                logging.info(f'{rna} ran with CONSEL (found *.sitelh).')
            else:
                # Retry with reduced inputs
                run_command(f"rm -rf {persite_path}")
                os.makedirs(persite_path, exist_ok=True)

                combineTreeFiles(persite_path, rna, dnaTree, rnaTree)
                # Use reduced SS and reduced FASTA if present
                if isfile(join(DIR_FASTA, f'{rna}.nodup.fa.reduced')):
                    fasta_file_red = join(DIR_FASTA, f'{rna}.nodup.fa.reduced')
                else:
                    logging.warning(f'No reduced FASTA found for {rna}; skipping reduced attempt.')
                    fasta_file_red = None

                if isfile(join(DIR_SS, f'{rna}.ss.reduced')):
                    ss_file_red = join(DIR_SS, f'{rna}.ss.reduced')
                else:
                    logging.warning(f'No reduced SS found for {rna}; skipping reduced attempt.')
                    ss_file_red = None

                if fasta_file_red is not None:
                    runningCONSEL(rna, fasta_file_red, ss_file_red, persite_suffix, persite_path, prefix_consel, MODEL)
                    if has_sitelh(persite_path):
                        logging.info(f"{rna}: reduced attempt produced *.sitelh.")
                    else:
                        logging.error(f"{rna}: reduced attempt still did not produce *.sitelh.")
        else:
            logging.warning(f"{rna} has either under DNA or RNA model no tree inference.")
        
if __name__=='__main__':
    main()
