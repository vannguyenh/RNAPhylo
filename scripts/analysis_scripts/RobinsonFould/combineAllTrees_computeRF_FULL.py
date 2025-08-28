import os
from os.path import join
from Bio import Phylo
import pandas as pd 

import subprocess
import logging
from datetime import datetime

import re, sys, csv, pathlib

LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'

MODEL = 'S6A'
# change to the right pathway of DIR_DNA for FULL or SEED data
DIR_WORKING = "/Users/u7875558/RNAPhylo/fullAlignment_S6A"
DIR_OUTPUTS = join(DIR_WORKING, "outputs")
DIR_RF_LOGS = join(DIR_WORKING, "logs", "RF_distance")
os.makedirs(DIR_RF_LOGS, exist_ok=True)

DIR_DNA=join(DIR_OUTPUTS, "raxml")
DIR_RNA=join(DIR_OUTPUTS, "raxmlP_iPseu")

def check_branch_length(dir_output, rna):
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

def produceCombinedTree( dir_output, rna):
    os.makedirs(dir_output, exist_ok=True)
    # endings for output files
    endings= {
        "raxml":"raxml",
        "raxmlP_iPseu":"raxmlPi"
    }

    # list all types (DNA, RNA) in the output directory
    types = list(endings.keys())
    for type in types:
        if type == 'raxml':
            input_data = join(DIR_DNA, rna)
            suffix = endings[type]
        else:
            input_data = join(DIR_RNA, rna)
            suffix = endings[type]
        
        input_files = sorted([f for f in os.listdir(input_data) if f.startswith('RAxML_bestTree')])

        if len(input_files) == 0:
            logging.warning(f"{rna} of {type} does not have any tree file.")
            continue
        else:
            output_rna = join(dir_output, rna)
            os.makedirs(output_rna, exist_ok=True)
            output_file = join(output_rna, f"{rna}.{suffix}")
            with open(output_file, "w") as outfile:
                for filename in input_files:
                    file_path = join(input_data, filename)
                    with open(file_path, "r") as infile:
                        content = infile.read()
                        outfile.write(content)

def run_command(command):
    try:
        process=subprocess.Popen(command, shell=True)
        process.communicate()
    except Exception as e:
        logging.error(f"Command failed: {command} with error: {e}")

def computeRFdistance_iqtreecmd(dcombine_path, rna):
    dir_combine_rna=join(dcombine_path, rna)
    for f in os.listdir(dir_combine_rna):
        if f.endswith('raxml'):
            raxTree=join(dir_combine_rna, f)
        elif f.endswith('raxmlPi'):
            raxPiPTree=join(dir_combine_rna, f)

    logging.info(f"Compute the RF distance of {rna}.")
    command=f"bash /Users/u7875558/Documents/promotion/projects/projects_code/RNAPhylo/scripts/analysis_scripts/RobinsonFould/computeRFdistance.sh {raxTree} {raxPiPTree} {dir_combine_rna} {rna}"
    run_command(command=command)

def main():
    dir_combined=join(DIR_OUTPUTS, 'Robinson_Foulds_iqtree3')
    os.makedirs(dir_combined, exist_ok=True)

    log_preprocess = join(DIR_WORKING, "logs", "full_S6A.log" )
    log_filename = os.path.join(DIR_RF_LOGS, f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.normalisedRF.iqtree3.{MODEL}.log")
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, format=LOG_FORMAT)
    logging.info(f"Running the code with the model {MODEL}.")

    rnas = extractAnalysedRNAs(log_file=log_preprocess)[0]

    for rna in rnas:
        if os.path.exists(join(DIR_RNA, rna)) and len(os.listdir(join(DIR_RNA, rna))) != 50:
            logging.warning(f"{rna} does not have 50 files in the inference folder.")
            continue
        else:
            path=join(dir_combined, rna)
            rfdist_file=0
            if os.path.isdir(path):
                for f in os.listdir(path):
                    if f.endswith('.rfdist'):
                        rfdist_file += 1
            if rfdist_file != 3:
                produceCombinedTree(dir_combined, rna)
                computeRFdistance_iqtreecmd(dir_combined, rna)
            else:
                logging.info(f"The RF computation of {rna} is done.")

if __name__=="__main__":
    main()