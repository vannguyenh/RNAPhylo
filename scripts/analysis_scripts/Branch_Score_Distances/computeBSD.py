import os
from os.path import join

from Bio import Phylo
import logging
from datetime import datetime
import subprocess

DIR_OUTPUTS = '/Users/u7875558/Documents/PhD/RNAPhylo/allModels_SEED/outputs'
DIR_DNA = '/Users/u7875558/Documents/PhD/RNAPhylo/allModels_SEED/outputs/DNAtrees'

DIR_INPUTS = '/Users/u7875558/Documents/PhD/RNAPhylo/allModels_SEED/inputs'
DIR_FASTA = join(DIR_INPUTS, 'fasta_files')
DIR_SS = join(DIR_INPUTS, 'ss_files')

LOG_FILE='/Users/u7875558/Documents/PhD/RNAPhylo/allModels_SEED/logs/2025-03-24_11-55-45.log'

SUBFOLDERS=['raxmlP_wPseu', 'raxmlP_iPseu']

def check_branch_length(diroutput, rna):
    # return RNAs which either do not have enough 10 trees or have any branch lengths > 1
    seeds = [f"{i:02d}" for i in range(1, 11)]
    dir_rna = os.path.join(diroutput, rna)

    delete_files = list()

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

def extractAnalysedRNAs(log_file):
    # This function produces two sets of accepted RNAs -- RNAs containing pseudoknots and RNAs containing no pseudoknots
    # extract accepted RNAs which do not have any branch length > 1
    rnas = os.listdir(DIR_DNA)
    unaccepted_rnas = list()
    accepted_rnas = list()
    for rna in rnas:
        # consider only the outputs from using DNA
        if check_branch_length(DIR_DNA, rna) == rna:
            #logging.warning(f'{rna} inferred from using DNA has branch length > 1.')
            unaccepted_rnas.append(rna)
        else:
            accepted_rnas.append(rna)

    # extract RNAs containing pseudoknots
    with open(log_file, 'r') as f:
        lines = f.readlines()

    nopseudo_rnas = [line.split()[5] for line in lines if 'not have pseudoknots' in line]
    accepted_nopseudo = set(accepted_rnas) & set(nopseudo_rnas)

    pseudo_rnas = set(rnas)-set(nopseudo_rnas)
    accepted_pseudo = set(accepted_rnas) & set(pseudo_rnas)

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
    raxPTree_path = join(dir_output, SUBFOLDERS[0], rna)
    if check_inferred_tree(raxTree_path) and check_inferred_tree(raxPTree_path):
        bestTreesLH = {'DNA': extract_highestloglh(raxTree_path),
                       'RNA considering pseudoknots': extract_highestloglh(raxPTree_path)}
        return bestTreesLH
    else:
        return None

def extract_highestLH_2Trees_ipseudo(dir_output, rna):
    raxTree_path = join(DIR_DNA, rna)
    raxPiTree_path = join(dir_output, SUBFOLDERS[1], rna)
    if check_inferred_tree(raxTree_path) and check_inferred_tree(raxPiTree_path):
        bestTreesLH = {'DNA': extract_highestloglh(raxTree_path),
                       'RNA ignoring pseudoknots': extract_highestloglh(raxPiTree_path)}
        return bestTreesLH
    else:
        return None

def take_highestLHTrees(dir_output, rna, dna_tree, rna_tree):
    dna_TreeFile = join(dir_output, f"{rna}.highestLH.dnatree")
    rna_TreeFile = join(dir_output, f"{rna}.highestLH.rnatree")

    with open(dna_TreeFile, 'w') as fo_dna:
        with open(dna_tree, 'r') as fi_dna:
            fo_dna.write(fi_dna.read())

    with open(rna_TreeFile, 'w') as fo_rna:
        with open(rna_tree, 'r') as fi_rna:
            fo_rna.write(fi_rna.read())

def run_bash(command):
    process = subprocess.Popen(command, shell=True)
    process.communicate()
    # Check process.returncode
    if process.returncode != 0:
        logging.error(f"Bash command failed: {command}")
    return process.returncode

def runningBSD(bsd_prefix, dnaTree, rnaTree):
    bash_file = '/Users/u7875558/Documents/PhD/RNAPhylo/scripts/analysis_scripts/Branch_Score_Distances/bsd.sh'
    bsd_command = f"bash {bash_file} {dnaTree} {rnaTree} {bsd_prefix}"
    run_bash(bsd_command)

def main():
    MODEL = input('Model: ')
    DIR_WORKING = join(DIR_OUTPUTS, 'BSD', MODEL)
    os.makedirs(DIR_WORKING, exist_ok=True)

    log_filename = join(DIR_WORKING, f"BSD_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.{MODEL}.log")
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

    all_rnas, pseudo_rnas, nopseudo_rnas, working_rnas, working_pseu = extractAnalysedRNAs(LOG_FILE)

    reduced_rnas = ['RF00001', 'RF00177', 'RF01960', 'RF02345', 'RF02401', 'RF02540', 'RF02842']

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

    for rna in working_rnas:
        dir_bsd_output = join(DIR_WORKING, 'ignore_pseudo', rna)
        os.makedirs(dir_bsd_output, exist_ok=True)
        bestTrees = extract_highestLH_2Trees_ipseudo(join(DIR_OUTPUTS, MODEL), rna)
        if bestTrees is not None:
            dnaTree = join(DIR_DNA, rna, f"RAxML_bestTree.{rna}.{bestTrees['DNA'][0]}")
            rnaTree = join(DIR_OUTPUTS, MODEL, 'raxmlP_iPseu', rna, f"RAxML_bestTree.{rna}.{bestTrees['RNA ignoring pseudoknots'][0]}")
            take_highestLHTrees(dir_bsd_output, rna, dnaTree, rnaTree)

            # running the computation of BSD via IQ-TREE
            try:
                runningBSD(f"{dir_bsd_output}/{rna}",
                           join(dir_bsd_output, f"{rna}.highestLH.dnatree"),
                           join(dir_bsd_output, f"{rna}.highestLH.rnatree"))
                logging.info(f"Computation of BSD between highest LH trees is running with {rna} ignoring pseudoknots.")
            except Exception as e:
                logging.error(f"Error with {rna}: {e}")
        else:
            logging.warning(f"{rna} ignoring pseudoknots has issue with the tree files.")

if __name__=='__main__':
    main()
