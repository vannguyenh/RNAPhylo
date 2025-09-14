#!/usr/bin/env python3

"""
Script to use all files (fasta files, secondary structure files) which were pre-prepared
to run phylogenetic inference on NCI.

DIR_WORKING = /scratch/dx61/vh5686/tmp/RNAPhylo/fullModels_FULLALIGN
"""

import os
from os.path import join
import pandas as pd
import logging
from datetime import datetime
import subprocess
from typing import Optional, Tuple, List, Dict

def control_fasta_file(fi_fasta):
    qualify = False
    num_seqs = 0
    with open(fi_fasta, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('>'):
                num_seqs += 1

    if num_seqs >= 4:
        qualify = True

    return qualify

def run_command(command: str) -> Optional[subprocess.CompletedProcess]:
    """Runs a shell command and logs the output."""
    process = subprocess.Popen(command, shell=True)
    process.communicate()
    
def main():
    RAXML_EXECUTE = '/scratch/dx61/vh5686/tmp/RNAPhylo/tools/standard-RAxML/raxmlHPC'
    DIR_WORKING = '/scratch/dx61/vh5686/tmp/RNAPhylo/fullModels_FULLALIGN'
    # list out folders for the inputs
    DIR_FASTA = join(DIR_WORKING, 'inputs', 'fasta')
    DIR_SS = join(DIR_WORKING, 'inputs', 'ss_all')
    DIR_SUBSAMP = join(DIR_WORKING, 'inputs', 'subsample')
    # log directory
    DIR_LOG = join(DIR_WORKING, 'logs')
    os.makedirs(DIR_LOG, exist_ok=True)
    # output directory
    MODEL = 'S16'
    DIR_OUTPUT = join(DIR_WORKING, 'outputs', MODEL)
    os.makedirs(DIR_OUTPUT, exist_ok=True)

    log_filename = os.path.join(DIR_LOG, f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.{MODEL}.log")
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info(f"Running code with the model {MODEL}!")

    rnas = [rna.split('.')[0] for rna in os.listdir(DIR_FASTA)]
    rnas_subsamp = [rna.split('.')[0] for rna in os.listdir(DIR_SUBSAMP)]

    for rna in rnas: # sort the rna -- the biggest first, then the smallest later.
        ss_file = join(DIR_SS, f"{rna}.ss")
        if rna in rnas_subsamp:
            fasta_file = join(DIR_SUBSAMP, f"{rna}.subsamp.fa")
        else:
            fasta_file = join(DIR_FASTA, f"{rna}.nodup.fa")

        if control_fasta_file(fasta_file):
            raxml_prefix = join(DIR_OUTPUT, 'raxml', rna)
            raxmlPiPseu_prefix = join(DIR_OUTPUT, 'raxmlP_iPseu', rna)

            # run the raxml command via bash files
            if not os.path.isdir(raxml_prefix) or len(os.listdir(raxml_prefix)) != 50:
                logging.warning(f"{raxml_prefix} has not been run yet or not run enough for 10 trees.")
                raxml_command = (
                    f"qsub -V -N raxml_{rna} -o {log_filename} -e {log_filename} "
                    f"-l ncpus=1 -l mem=4gb -l walltime=48:00:00 -l wd -- " #recommend: including threads or reduce ncpus to 1
                    f"bash /scratch/dx61/vh5686/tmp/RNAPhylo/scripts/bashFiles/raxml.sh {rna} {fasta_file} {raxml_prefix} {RAXML_EXECUTE}"
                )
                run_command(raxml_command)

            if not os.path.isdir(raxmlPiPseu_prefix) or len(os.listdir(raxmlPiPseu_prefix)) != 50:
                logging.warning(f"{raxmlPiPseu_prefix} has not been run yet or not run enough for 10 trees.")
                raxmlPiPseu_command = (
                    f"qsub -V -N raxmlP_{rna} -o {log_filename} -e {log_filename} "
                    f"-l ncpus=1 -l mem=48gb -l walltime=48:00:00 -l wd -- " # default of RAxML 1 core.
                    f"bash /scratch/dx61/vh5686/tmp/RNAPhylo/scripts/bashFiles/raxmlP.sh {rna} {fasta_file} {ss_file} {raxmlPiPseu_prefix} {MODEL} {RAXML_EXECUTE}"  
                )
                run_command(raxmlPiPseu_command)
        else:
            logging.info(f"{rna} does not have more than 3 sequences for the analysis.")

if __name__ == "__main__":
    main()