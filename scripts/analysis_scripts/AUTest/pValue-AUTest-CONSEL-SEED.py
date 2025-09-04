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
DIR_AU_LOGS = join(DIR_WORKING, 'logs', 'AU_test_RAxML-NG')
os.makedirs(DIR_AU_LOGS, exist_ok=True)

DIR_INPUTS = join(DIR_WORKING, 'inputs')
DIR_FASTA = join(DIR_INPUTS, 'fasta_files')
DIR_SS = join(DIR_INPUTS, 'ss_files')

LOG_FILE=join(DIR_WORKING, 'logs', 'SEED_AllModels.log')

SUBFOLDERS=['raxml', 'raxmlP_wPseu', 'raxmlP_iPseu']

def check_branch_length(dir_output, rna):
    # return RNAs which either do not have enough 10 trees or have any branch lengths > 1
    expected_files=[f"{i:02d}" for i in range(1,11)]
    dirRNA = join(dir_output, rna)
    delete_files=list()

    if dirRNA.startswith('RF'):
        for file_name in os.listdir(dirRNA):
            for seed in expected_files:
                if file_name.startswith('RAxML_bestTree') and file_name.endswith(seed):
                    tree_path=join(dirRNA, file_name)
                    tree=Phylo.read(tree_path, 'newick')
                    for clade in tree.find_clades():
                        if clade.branch_length and clade.branch_length > 1:
                            delete_files.append(file_name)
                            logging.info(f"{file_name} has branch length > 1.")
    else:
        logging.warning(f"{dirRNA} is not an RNA folder.")
    
    if len(delete_files) > 0:
        return rna
    
def extractAnalysedRNAs(log_file):
    # This function produces two sets of accepted RNAs -- RNAs containing pseudoknots and RNAs containing no pseudoknots    
    # extract accepted RNAs which do not have any branch length > 1
    rnas = os.listdir(DIR_DNA)

    brlength_larger1 =list()
    accepted_rnas=list()

    less4_pat = re.compile(r'Number of sequences of (RF\d{5}) is less than 4\.')
    unpaired_pat = re.compile(r'The secondary structure of (RF\d{5}) has only unpaired bases\.')  # note: "unpaired"
    less4, unpaired = set(), set()

    with open(log_file, 'r') as lf: #log file used for this full_S6A.log
        for line in lf:
            m1 = less4_pat.search(line)
            if m1: less4.add(m1.group(1))
            m2 = unpaired_pat.search(line)
            if m2: unpaired.add(m2.group(1))
    
    #both = less4 & unpaired
    #less_only = sorted(less4 - unpaired)
    #unpaired_only = sorted(unpaired - less4)

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

def extract_highestLH_2Trees_pseudo(dir_output, rna):
    raxTree_path = join(DIR_DNA, rna)
    raxPTree_path = join(dir_output, SUBFOLDERS[1], rna)
    if check_inferred_tree(raxTree_path) and check_inferred_tree(raxPTree_path):
        bestTreesLH = {'DNA': extract_highestloglh(raxTree_path),
                       'RNA considering pseudoknots': extract_highestloglh(raxPTree_path)}
        return bestTreesLH
    else:
        return None
    
def extract_highestLH_2Trees_ipseudo_dnaextra(dir_output, rna):
    raxTree_path = join(DIR_DNA, rna)
    raxPiTree_path = join(dir_output, rna)
    if check_inferred_tree(raxTree_path) and check_inferred_tree(raxPiTree_path):
        bestTreesLH = {'DNA': extract_highestloglh(raxTree_path),
                       'RNA_otherDNA': extract_highestloglh(raxPiTree_path)}
        return bestTreesLH
    else:
        return None

def combineTreeFiles(dir_combined, rna, group, tree1, tree2):
    combined_TreeFile = join(dir_combined, f"{rna}.{group}.highestLH.trees")

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

def runningCONSEL(rna, group, fasta_file, ss_file, persite_suffix, persite_path, prefix_consel, dir_combined):
    combineTree = join(dir_combined, f'{rna}.{group}.highestLH.trees')
    output_persite = join(persite_path, f'RAxML_perSiteLLs.{persite_suffix}')
    
    bash_file = '/Users/u7875558/Documents/RNAPhylo/RNAPhylo/scripts/analysis_scripts/AUTest/consel.sh'
    consel_command = f"bash {bash_file} {fasta_file} {combineTree} {persite_suffix} {ss_file} {persite_path} {output_persite} {prefix_consel}"

    run_bash(consel_command)

def run_command(command):
    process = subprocess.Popen(command, shell=True)
    process.communicate()

def has_sitelh(dirpath: str) -> bool:
    """Return True if any file in dirpath ends with *sitelh (case-insensitive)."""
    try:
        return any(fn.lower().endswith('sitelh') for fn in os.listdir(dirpath))
    except FileNotFoundError:
        return False
    
def main():
    MODEL = input('Model: ')
    DIR_WORKING = join(DIR_OUTPUTS, 'AU_Test_RAxMLNG', MODEL)
    os.makedirs(DIR_WORKING, exist_ok=True)
    DIR_COMBINE = join(DIR_WORKING, 'combine_2trees_highestLH')
    os.makedirs(DIR_COMBINE, exist_ok=True)

    log_filename = join(DIR_WORKING, f"CONSEL_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.{MODEL}.log")
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

    rnas = extractAnalysedRNAs(LOG_FILE)[0]

    # for rna in working_pseu:
    #     bestTrees = extract_highestLH_2Trees_pseudo(join(DIR_OUTPUTS, MODEL), rna)
    #     if bestTrees is not None:
    #         dnaTree = join(DIR_OUTPUTS, MODEL, 'raxml', rna, f"RAxML_bestTree.{rna}.{bestTrees['DNA'][0]}")
    #         rnaTree = join(DIR_OUTPUTS, MODEL, 'raxmlP_wPseu', rna, f"RAxML_bestTree.{rna}.{bestTrees['RNA considering pseudoknots'][0]}")
    #         combineTreeFiles(DIR_COMBINE, rna, 'wPseu', dnaTree, rnaTree)

    #         persite_path = join(DIR_WORKING, 'consider_pseudo', rna)
    #         os.makedirs(persite_path, exist_ok=True)
    #         persite_pseudo_suffix = f'{rna}.wpseu.sitelh'
    #         prefix_pseudo_consel = join(persite_path, f'{rna}_wpseu_consel')
            
    #         if rna in reduced_rnas:
    #             ss_pseudo_file = join(DIR_SS, f'{rna}.brackets.ss.reduced')
    #             fasta_file = join(DIR_FASTA, f'{rna}.nodup.fa.reduced')
    #         else:
    #             ss_pseudo_file = join(DIR_SS, f'{rna}.brackets.ss')
    #             fasta_file = join(DIR_FASTA, f'{rna}.nodup.fa')
                
    #         # running consel
    #         try:
    #             runningCONSEL(rna, MODEL, 'wPseu', fasta_file, ss_pseudo_file, persite_pseudo_suffix, persite_path, prefix_pseudo_consel, DIR_COMBINE)
    #             logging.info(f"CONSEL is running with {rna} having pseudoknots.")
    #         except Exception as e:
    #             logging.error(f"Error with {rna}: {e}")
    #     else:
    #         logging.warning(f"{rna} considering pseudoknots has issue with the tree files.")
    
    for rna in rnas:
        bestTrees = extract_highestLH_2Trees_ipseudo_dnaextra(join(DIR_TREES, MODEL), rna)
        if bestTrees is not None:
            dnaTree = join(DIR_DNA, rna, f"RAxML_bestTree.{rna}.{bestTrees['DNA'][0]}")
            rnaTree = join(DIR_TREES, MODEL, rna, 
                           f"RAxML_bestTree.{rna}.{bestTrees['RNA_otherDNA'][0]}")
            combineTreeFiles(DIR_COMBINE, rna, 'iPseu', dnaTree, rnaTree)
        else:
            logging.warning(f"{rna} has either under DNA or RNA model no tree inference.")
        
        if isfile(join(DIR_COMBINE, f'{rna}.iPseu.highestLH.trees')):
            persite_path = join(DIR_WORKING, 'ignore_pseudo', rna)
            os.makedirs(persite_path, exist_ok=True)
            persite_suffix = f'{rna}.ipseu.sitelh'
            prefix_consel = join(persite_path, f'{rna}_ipseu_consel')

            ss_file = join(DIR_SS, f'{rna}.dots.ss')
            fasta_file = join(DIR_FASTA, f'{rna}.nodup.fa')
                
            try:
                runningCONSEL(rna, MODEL, 'iPseu', fasta_file, ss_file, persite_suffix, 
                              persite_path, prefix_consel, DIR_COMBINE)
                logging.info(f"CONSEL is running with {rna} ignoring pseudoknots.")
            except Exception as e:
                logging.error(f"Error with {rna}: {e}")

            # If we already have a *.sitelh, move on
            if has_sitelh(persite_path):
                logging.info(f'{rna} ran with CONSEL (found *.sitelh).')
            else:
                # Retry with reduced inputs
                run_command(f"rm -rf {persite_path}")
                os.makedirs(persite_path, exist_ok=True)

                # Use reduced SS and reduced FASTA if present
                if isfile(join(DIR_FASTA, f'{rna}.nodup.fa.reduced')):
                    fasta_file_red = join(DIR_FASTA, f'{rna}.nodup.fa.reduced')
                else:
                    logging.warning(f'No reduced FASTA found for {rna}; skipping reduced attempt.')
                    fasta_file_red = None

                if fasta_file_red is not None:
                    runningCONSEL(
                        rna,'iPseu',
                        fasta_file_red, ss_file,
                        persite_suffix, persite_path, prefix_consel, DIR_COMBINE
                    )
                    if has_sitelh(persite_path):
                        logging.info(f"{rna}: reduced attempt produced *.sitelh.")
                    else:
                        logging.error(f"{rna}: reduced attempt still did not produce *.sitelh.")
        else:
            logging.warning(f'{rna} has no combined file.')

    
if __name__=='__main__':
    main()
