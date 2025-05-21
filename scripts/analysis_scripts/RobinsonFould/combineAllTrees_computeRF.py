import os
from os.path import join
from Bio import Phylo
import pandas as pd 

import subprocess
import logging
from datetime import datetime

LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'

DIR_WORKING = "/Users/u7875558/Documents/RNAPhylo/allModels_SEED/"
DIR_OUTPUTS = os.path.join(DIR_WORKING, "outputs")
DIR_RF_LOGS = os.path.join(DIR_WORKING, 'logs', 'RF_distance')
DIR_DNA = '/Users/u7875558/Documents/RNAPhylo/allModels_SEED/outputs/DNAtrees'
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
    rnas = os.listdir(DIR_DNA)
    
    unaccepted_rnas=list()
    accepted_rnas=list()
    
    for rna in rnas:
        # consider only the outputs from using DNA
        if check_branch_length(DIR_DNA, rna) == rna:
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

def produceCombinedTrees(dir_input, dir_output, rna):
    os.makedirs(dir_output, exist_ok=True)
    # endings for output files
    endings={
        "raxml": "raxml",
        #"raxmlP_wPseu": "raxmlPw",
        "raxmlP_iPseu" : "raxmlPi"
    }

    # List all files in the input directory
    methods = list(endings.keys())
    #rnas = os.listdir(os.path.join(dir_input, methods[0]))
    for method in methods:
        #for rna in rnas:
        if method == 'raxml':
            input_method = join(DIR_DNA, rna)
        else:
            input_method = os.path.join(dir_input, method, rna)
        input_files = sorted([f for f in os.listdir(input_method) if f.startswith("RAxML_bestTree.")])

        if len(input_files) == 0:
            logging.warning(f"{rna} of {method} does not have tree files.")
        else:
            output_rna = os.path.join(dir_output, rna)
            os.makedirs(output_rna, exist_ok=True)
            output_file = os.path.join(output_rna, f"{rna}.{endings[method]}")
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

def computeRFdistance_iqtreecmd(dcombine_path, rna):
    #rnas=os.listdir(dcombine_path)
    #for rna in rnas:
    dir_combine_rna = os.path.join(dcombine_path, rna)
    for f in os.listdir(dir_combine_rna):
        if f.endswith("raxml"):
            raxTree = os.path.join(dir_combine_rna, f)
        #elif f.endswith("raxmlPw"):
        #    raxPwPTree = os.path.join(dir_combine_rna, f)
        elif f.endswith("raxmlPi"):
            raxPiPTree = os.path.join(dir_combine_rna, f)
    logging.info(f"Compute the RF distances of {rna}.")
    #command=f"bash /Users/u7875558/Documents/PhD/RNAPhylo/scripts/analysis_scripts/RobinsonFould/computeRFdistance.sh {raxTree} {raxPwPTree} {raxPiPTree} {dir_combine_rna} {rna}"
    command=f"bash /Users/u7875558/Documents/RNAPhylo/scripts/analysis_scripts/RobinsonFould/computeRFdistance.sh {raxTree} {raxPiPTree} {dir_combine_rna} {rna}"
    run_command(command)

def main():
    MODEL = input('Model: ')
    dir_combined = os.path.join(DIR_OUTPUTS, 'Robinson_Foulds', MODEL)
    os.makedirs(dir_combined, exist_ok=True)

    log_filename = os.path.join(DIR_RF_LOGS, f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.{MODEL}.log")
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, format=LOG_FORMAT)
    logging.info(f"Running the code with the model {MODEL}.")

    rnas, pseudo_rnas, nopseudo_rnas = extractAnalysedRNAs(join(DIR_OUTPUTS, MODEL), LOG_FILE)
    # issue_str_rnas: RNA families running with DNA model, but not with RNA models.
    issue_rnas = ['RF00207', 'RF00976', 'RF01047', 'RF01338', 'RF01380', 'RF03623', 'RF03760', 'RF03969']
    working_rnas = set(rnas) - set(issue_rnas)

    for rna in working_rnas:
        path = os.path.join(dir_combined, rna)
        rfdist_file = 0
        if os.path.isdir(path):
            for f in os.listdir(path):
                if f.endswith('.rfdist'):
                    rfdist_file += 1

        if rfdist_file != 3: # 3 while using only ramxl and raxmlP_iPseu, 6 while using raxml, raxmlP_iPseu and raxmlP_wPseu
            #run_command(f"rm -r {path}")
            produceCombinedTrees(os.path.join(DIR_OUTPUTS, MODEL), dir_combined, rna)
            computeRFdistance_iqtreecmd(dir_combined, rna)
        else:
            continue

if __name__=="__main__":
    main()
