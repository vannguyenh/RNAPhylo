import os
from os.path import join

from Bio import Phylo
import logging
from datetime import datetime
import time
import subprocess

DIR_OUTPUTS = '/Users/u7875558/Documents/PhD/RNAPhylo/allModels_SEED/outputs'
MODEL = 'S6A'
DIR_WORKING = join(DIR_OUTPUTS, 'AU_Test', MODEL)
os.makedirs(DIR_WORKING, exist_ok=True)
DIR_COMBINE = join(DIR_WORKING, 'combine_2trees_highestLH')
os.makedirs(DIR_COMBINE, exist_ok=True)

DIR_INPUTS = '/Users/u7875558/Documents/PhD/RNAPhylo/allModels_SEED/nci_inputs'
DIR_FASTA = join(DIR_INPUTS, 'fasta_files')
DIR_SS = join(DIR_INPUTS, 'ss_files')

TREE_OUTPUT = join(DIR_OUTPUTS, 'combinedFiles')
LOG_FILE='/Users/u7875558/Documents/PhD/RNAPhylo/allModels_SEED/logs/2025-03-24_11-55-45.log'

SUBFOLDERS=['raxml', 'raxmlP_wPseu', 'raxmlP_iPseu']

def check_branch_length(diroutput, rna):
    # return RNAs which either do not have enough 10 trees or have any branch lengths > 1
    seeds = [f"{i:02d}" for i in range(1, 11)]
    dir_rna = os.path.join(diroutput, rna)

    delete_files=list()

    for file_name in os.listdir(dir_rna):
        for seed in seeds:
            if file_name.startswith('RAxML_bestTree') and file_name.endswith(seed):
                tree_path = os.path.join(dir_rna, file_name)
                tree = Phylo.read(tree_path, 'newick')
                for clade in tree.find_clades():
                    if clade.branch_length and clade.branch_length > 1:
                        delete_files.append(file_name)
    
    if len(delete_files) > 0:
        return rna
    
def extractAnalysedRNAs(dirmodel, log_file):
    # This function produces two sets of accepted RNAs -- RNAs containing pseudoknots and RNAs containing no pseudoknots    
    # extract accepted RNAs which do not have any branch length > 1
    rnas = os.listdir(join(dirmodel, SUBFOLDERS[0]))
    unaccepted_rnas = list()
    accepted_rnas = list()
    for rna in rnas:
        # consider only the outputs from using DNA
        dir_working = join(dirmodel, SUBFOLDERS[0])
        if check_branch_length(dir_working, rna) == rna:
            #logging.warning(f'{rna} inferred from using DNA has branch length > 1.')
            unaccepted_rnas.append(rna)
        else:
            accepted_rnas.append(rna)

    # extract RNAs containing pseudoknots
    with open(log_file, 'r') as f:
        lines = f.readlines()

    nopseudo_rnas = [line.split()[5] for line in lines if 'not have pseudoknots' in line ]
    accepted_nopseudo=set(accepted_rnas) & set(nopseudo_rnas)
    
    pseudo_rnas=set(rnas)-set(nopseudo_rnas)
    accepted_pseudo=set(accepted_rnas) & set(pseudo_rnas)

    issue_str_rnas = ['RF00207', 'RF00390', 'RF01380', 'RF01338', 'RF01047', 'RF03760', 'RF03969', 'RF00976', 'RF03623']
    working_rnas = set(accepted_rnas) - set(issue_str_rnas)
    working_pseudo = set(accepted_pseudo) - set(issue_str_rnas)

    logging.info(f'{len(accepted_pseudo)} RNAs having pseudoknots.')
    logging.info(f'{len(accepted_nopseudo)} RNAs having no pseudoknots.')
    logging.info(f'RNAs having pseudoknots: {accepted_pseudo}')
    logging.info(f'RNA having no pseudokntos: {accepted_nopseudo}')
    logging.info(f'RNAs having branch length > 1 in DNA trees: {unaccepted_rnas}')
    logging.info(f'Number of working RNA families: {working_rnas}')
    return accepted_rnas, accepted_pseudo, accepted_nopseudo, working_rnas, working_pseudo


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
    print(lh)
    return sorted(lh.items(), key=lambda x: x[1], reverse=True)[0]

def extract_highestLH_2Trees_pseudo(dir_output, rna):
    raxTree_path = join(dir_output, SUBFOLDERS[0], rna)
    raxPTree_path = join(dir_output, SUBFOLDERS[1], rna)
    if check_inferred_tree(raxTree_path) and check_inferred_tree(raxPTree_path):
        bestTreesLH = {'DNA': extract_highestloglh(raxTree_path),
                       'RNA considering pseudoknots': extract_highestloglh(raxPTree_path)}
        return bestTreesLH
    else:
        return None
    
def extract_highestLH_2Trees_ipseudo(dir_output, rna):
    raxTree_path = join(dir_output, SUBFOLDERS[0], rna)
    raxPiTree_path = join(dir_output, SUBFOLDERS[2], rna)
    if check_inferred_tree(raxTree_path) and check_inferred_tree(raxPiTree_path):
        bestTreesLH = {'DNA': extract_highestloglh(raxTree_path),
                       'RNA ignoring pseudoknots': extract_highestloglh(raxPiTree_path)}
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

def runningCONSEL(rna, group, fasta_file, ss_file, persite_suffix, persite_path, prefix_consel):
    combineTree = join(DIR_COMBINE, f'{rna}.{group}.highestLH.trees')
    output_persite = join(persite_path, f'RAxML_perSiteLLs.{persite_suffix}')
    
    bash_file = '/Users/u7875558/Documents/PhD/RNAPhylo/scripts/analysis_scripts/AUTest/consel.sh'
    
    consel_command=f"bash {bash_file} {fasta_file} {combineTree} {persite_suffix} {ss_file} {persite_path} {output_persite} {prefix_consel} {MODEL}"
    run_bash(consel_command)

def main():
    log_filename = join(DIR_WORKING, f"CONSEL_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log")
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

    all_rnas, pseudo_rnas, nopseudo_rnas, working_rnas, working_pseu = extractAnalysedRNAs(join(DIR_OUTPUTS, MODEL), LOG_FILE)

    reduced_rnas = ['RF00001', 'RF00177', 'RF01960', 'RF02345', 'RF02401', 'RF02540', 'RF02842', 'RF02976']
    
    for rna in working_pseu:
        bestTrees = extract_highestLH_2Trees_pseudo(join(DIR_OUTPUTS, MODEL), rna)
        if bestTrees is not None:
            dnaTree = join(DIR_OUTPUTS, MODEL, 'raxml', rna, f"RAxML_bestTree.{rna}.{bestTrees['DNA'][0]}")
            rnaTree = join(DIR_OUTPUTS, MODEL, 'raxmlP_wPseu', rna, f"RAxML_bestTree.{rna}.{bestTrees['RNA considering pseudoknots'][0]}")
            combineTreeFiles(DIR_COMBINE, rna, 'wPseu', dnaTree, rnaTree)

            persite_path = join(DIR_WORKING, 'consider_pseudo', rna)
            os.makedirs(persite_path, exist_ok=True)
            persite_pseudo_suffix = f'{rna}.wpseu.sitelh'
            prefix_pseudo_consel = join(persite_path, f'{rna}_wpseu_consel')
            
            if rna in reduced_rnas:
                ss_pseudo_file = join(DIR_SS, f'{rna}.brackets.ss.reduced')
                fasta_file = join(DIR_FASTA, f'{rna}.nodup.fa.reduced')
            else:
                ss_pseudo_file = join(DIR_SS, f'{rna}.brackets.ss')
                fasta_file = join(DIR_FASTA, f'{rna}.nodup.fa')
                
            # running consel
            try:
                runningCONSEL(rna, 'wPseu', fasta_file, ss_pseudo_file, persite_pseudo_suffix, persite_path, prefix_pseudo_consel)
                logging.info(f"CONSEL is running with {rna} having pseudoknots.")
            except Exception as e:
                logging.error(f"Error with {rna}: {e}")
        else:
            logging.warning(f"{rna} considering pseudoknots has issue with the tree files.")
    
    for rna in working_rnas:
        bestTrees = extract_highestLH_2Trees_ipseudo(join(DIR_OUTPUTS, MODEL), rna)
        if bestTrees is not None:
            dnaTree = join(DIR_OUTPUTS, MODEL, 'raxml', rna, f"RAxML_bestTree.{rna}.{bestTrees['DNA'][0]}")
            rnaTree = join(DIR_OUTPUTS, MODEL, 'raxmlP_iPseu', rna, f"RAxML_bestTree.{rna}.{bestTrees['RNA ignoring pseudoknots'][0]}")
            combineTreeFiles(DIR_COMBINE, rna, 'iPseu', dnaTree, rnaTree)
            
            persite_path = join(DIR_WORKING, 'ignore_pseudo', rna)
            os.makedirs(persite_path, exist_ok=True)
            persite_suffix = f'{rna}.ipseu.sitelh'
            prefix_consel = join(persite_path, f'{rna}_ipseu_consel')

            if rna in reduced_rnas:
                ss_pseudo_file = join(DIR_SS, f'{rna}.dots.ss.reduced')
                fasta_file = join(DIR_FASTA, f'{rna}.nodup.fa.reduced')
            else:
                ss_pseudo_file = join(DIR_SS, f'{rna}.dots.ss')
                fasta_file = join(DIR_FASTA, f'{rna}.nodup.fa')
            
            # running consel 
            try:
                runningCONSEL(rna, 'iPseu', fasta_file, ss_file, persite_suffix, persite_path, prefix_consel)
                logging.info(f"CONSEL is running with {rna} ignoring pseudoknots.")
            except Exception as e:
                logging.error(f"Error with {rna}: {e}")
        else:
            logging.warning(f"{rna} ignoring pseudoknots has issue with the tree files.")
    
if __name__=='__main__':
    main()
