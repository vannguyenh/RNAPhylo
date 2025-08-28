import os
from os.path import join
from Bio import Phylo
import numpy as np
import pandas as pd
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from datetime import datetime

MODEL = 'S6A'
LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'

# ─── PARAMETERS ────────────────────────────────────────────────────────────────
DIR_WORKING = "/Users/u7875558/RNAPhylo/fullAlignment_S6A"
DIR_OUTPUTS = join(DIR_WORKING, "outputs")
DIR_RF = join(DIR_OUTPUTS, 'Robinson_Foulds')
DIR_RF_LOGS = join(DIR_WORKING, "logs", "RF_distance")
os.makedirs(DIR_RF_LOGS, exist_ok=True)

SUFFIXES     = {
    'DNA vs DNA':  '.raxml.rfdist',
    'DNA vs RNA':  '.raxml.raxmlPi.rfdist',
}

# ─── UTILITY FUNCTIONS ────────────────────────────────────────────────────────

def read_rfdist(path):
    """Read an RF distance matrix file into a NumPy array."""
    with open(path) as f:
        lines = f.readlines()[1:]  # skip header
    mat = [list(map(float, row.strip().split()[1:])) for row in lines]
    return np.array(mat)

def summarize_normalized(mat):
    """L2-normalize a matrix and return its mean and median (flattened)."""
    norm = normalize(mat, norm='l2')
    flat = norm.flatten()
    return float(flat.mean()), float(np.median(flat))

# ─── MAIN PROCESSING ──────────────────────────────────────────────────────────

log_filename = os.path.join(DIR_RF_LOGS, f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.{MODEL}.log")
logging.basicConfig(filename=log_filename, level=logging.DEBUG, format=LOG_FORMAT)
logging.info(f"Running the code with the model {MODEL}.")

all_records = []

for rna in os.listdir(DIR_RF):
    rna_dir = join(DIR_RF, rna)
    
    for category, suffix in SUFFIXES.items():
        file_path = join(rna_dir, f"{rna}{suffix}")

        if not os.path.exists(file_path):
            continue
        mat = read_rfdist(file_path)
        mean_rf, med_rf = summarize_normalized(mat)

        all_records.append({
            'RNA': rna,
            'Category': category,
            'Mean RF': mean_rf,
            'Median RF': med_rf
        })

df = pd.DataFrame(all_records)