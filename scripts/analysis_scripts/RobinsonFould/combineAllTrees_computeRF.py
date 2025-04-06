import os
from os.path import join
from Bio import Phylo
import pandas as pd 

import subprocess
import logging
from datetime import datetime

LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'

DIR_WORKING = "/Users/u7875558/Documents/PhD/RNAPhylo/allModels_SEED/"
DIR_OUTPUTS = os.path.join(DIR_WORKING, "outputs")
DIR_RF_LOGS = os.path.join(DIR_WORKING, 'logs', 'RF_distance')
os.makedirs(DIR_RF_LOGS, exist_ok=True)

LOG_FILE = join(DIR_WORKING, 'logs/2025-03-24_11-55-45.log')
SUBFOLDERS = ['raxml', 'raxmlP_wPseu', 'raxmlP_iPseu']


def check_branch_length(diroutput, rna):
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
    
    if len(delete_files) > 0:
        return rna

def extractAnalysedRNAs(diroutput, log_file):
    # This function produces two sets of accepted RNAs -- RNAs containing pseudoknots and RNAs containing no pseudoknots    
    # extract accepted RNAs which do not have any branch length > 1
    rnas = os.listdir(join(diroutput, SUBFOLDERS[0]))
    
    unaccepted_rnas=list()
    accepted_rnas=list()
    
    for rna in rnas:
        # consider only the outputs from using DNA
        dir_working = join(diroutput, SUBFOLDERS[0])
        if check_branch_length(dir_working, rna) == rna:
            unaccepted_rnas.append(rna)
        else:
            accepted_rnas.append(rna)

    # extract RNAs containing pseudoknots
    with open(log_file, 'r') as f:
        lines=f.readlines()

    nopseudo_rnas = [line.split()[5] for line in lines if 'not have pseudoknots' in line ]
    accepted_nopseudo=set(accepted_rnas) & set(nopseudo_rnas)
    
    pseudo_rnas = set(rnas)-set(nopseudo_rnas)
    accepted_pseudo = set(accepted_rnas) & set(pseudo_rnas)
    
    return accepted_rnas, accepted_pseudo, accepted_nopseudo

def produceCombinedTrees(dir_input, dir_output):
    os.makedirs(dir_output, exist_ok=True)
    # endings for output files
    endings={
        "raxml": "raxml",
        "raxmlP_wPseu": "raxmlPw",
        "raxmlP_iPseu" : "raxmlPi"
    }

    # List all files in the input directory
    methods = list(endings.keys())
    rnas = os.listdir(os.path.join(dir_input, methods[0]))
    for method in methods:
        for rna in rnas:
            input_method = os.path.join(dir_input, method, rna)
            input_files = sorted([f for f in os.listdir(input_method) if f.startswith("RAxML_bestTree.")])

            if len(input_files) == 0:
                logging.warning(f"{rna} of {method} does not have tree files.")
            else:
                output_rna=os.path.join(dir_output, rna)
                os.makedirs(output_rna, exist_ok=True)
                output_file=os.path.join(output_rna, f"{rna}.{endings[method]}")
                with open(output_file, "w") as outfile:
                    for filename in input_files:
                        file_path = os.path.join(input_method, filename)
                        with open(file_path, "r") as infile:
                            content = infile.read()
                            outfile.write(content)

def run_command(command):
    try:
        process = subprocess.Popen(command, shell=True)
        process.communicate()
    except Exception as e:
        logging.error(f"Command failed: {command} with error: {e}")

def computeRFdistance_iqtreecmd(dcombine_path, rnas):
    #rnas=os.listdir(dcombine_path)
    for rna in rnas:
        dir_combine_rna = os.path.join(dcombine_path, rna)
        for f in os.listdir(dir_combine_rna):
            if f.endswith("raxml"):
                raxTree = os.path.join(dir_combine_rna, f)
            elif f.endswith("raxmlPw"):
                raxPwPTree = os.path.join(dir_combine_rna, f)
            elif f.endswith("raxmlPi"):
                raxPiPTree = os.path.join(dir_combine_rna, f)
        #prefix=f"{dcombine_path}/{rna}/"
        logging.info(f"Compute the RF distances of {rna}.")
        command=f"bash computeRFdistance.sh {raxTree} {raxPwPTree} {raxPiPTree} {dir_combine_rna} {rna}"
        run_command(command)

def main():
    MODEL = input('Model: ')
    dir_combined = os.path.join(DIR_OUTPUTS, 'Robinson_Foulds', MODEL)

    log_filename = os.path.join(DIR_RF_LOGS, f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log")
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, format=LOG_FORMAT)
    logging.info(f"Running the code with the model {MODEL}.")

    rnas, pseudo_rnas, nopseudo_rnas = extractAnalysedRNAs(join(DIR_OUTPUTS, MODEL), LOG_FILE)

    produceCombinedTrees(os.path.join(DIR_OUTPUTS, MODEL), dir_combined)
    computeRFdistance_iqtreecmd(dir_combined, rnas)

if __name__=="__main__":
    main()


