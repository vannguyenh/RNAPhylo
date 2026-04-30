#!/usr/bin/env python3
"""
Extract the seed numbers of the best trees (highest log-likelihood) from RAxML log files.

This script reads existing RAxML log files and determines which seed (01-10)
produced the tree with the highest log-likelihood for each RNA family and model.

Output:
    - CSV file with columns: RNA_family, DNA_seed, S16_seed, S16A_seed, ...
    - This can be used to run bootstrap analysis with the correct seed numbers

Usage:
    python extract_best_seeds.py
    python extract_best_seeds.py --rna RF00740
    python extract_best_seeds.py --output best_seeds.csv
"""

import os
from os.path import join, isdir, isfile
import argparse
import csv

# Configuration - adjust these paths as needed
DIR_WORKING = '/Users/u7875558/RNAPhylo/seedAlignment_AllModels'
DIR_OUTPUTS = join(DIR_WORKING, 'outputs')
DIR_TREES = join(DIR_OUTPUTS, 'inferred_trees')
DIR_DNA = join(DIR_TREES, 'DNA')

# All RNA models
RNA_MODELS = ['S16', 'S16A', 'S16B', 'S7A', 'S7B', 'S7C', 'S7D', 'S7E', 'S7F', 'S6A', 'S6B', 'S6C', 'S6D', 'S6E']


def extract_highest_loglh(dir_path):
    """
    Extract the seed number with the highest log-likelihood from RAxML log files.

    Returns: (seed, log_likelihood) e.g., ('04', -1234.56)
             or None if no log files found
    """
    lh = dict()
    seeds = [f"{i:02d}" for i in range(1, 11)]

    if not isdir(dir_path):
        return None

    for file in os.listdir(dir_path):
        for seed in seeds:
            if file.startswith('RAxML_log') and file.endswith(seed):
                filepath = join(dir_path, file)
                try:
                    with open(filepath, 'r') as f:
                        lines = f.readlines()
                        if lines:
                            final_result = lines[-1].strip()
                            best_value = final_result.split()[-1]
                            lh[seed] = float(best_value)
                except Exception as e:
                    print(f"    Warning: Could not read {filepath}: {e}")

    if not lh:
        return None

    # Return the seed with highest (least negative) log-likelihood
    best = sorted(lh.items(), key=lambda x: x[1], reverse=True)[0]
    return best


def extract_all_seeds(rna_family):
    """
    Extract best seed numbers for one RNA family across DNA and all RNA models.

    Returns: dict with keys 'DNA', 'S16', 'S16A', etc. and values as (seed, loglh) tuples
    """
    results = {}

    # DNA
    dna_path = join(DIR_DNA, rna_family)
    dna_result = extract_highest_loglh(dna_path)
    if dna_result:
        results['DNA'] = dna_result

    # RNA models
    for model in RNA_MODELS:
        model_path = join(DIR_TREES, model, rna_family)
        model_result = extract_highest_loglh(model_path)
        if model_result:
            results[model] = model_result

    return results


def main():
    parser = argparse.ArgumentParser(description='Extract best seed numbers from RAxML log files')
    parser.add_argument('--rna', type=str, help='Process only this RNA family (e.g., RF00740)')
    parser.add_argument('--output', type=str, default='best_seeds.csv', help='Output CSV file')
    #parser.add_argument('--trees_dir', type=str, default=DIR_TREES, help='Path to inferred_trees directory')

    args = parser.parse_args()

    #global DIR_TREES, DIR_DNA
    #DIR_TREES = args.trees_dir
    DIR_DNA = join(DIR_TREES, 'DNA')

    # Get list of RNA families to process
    if args.rna:
        rna_families = [args.rna]
    else:
        # Get all RNA families from DNA directory
        if isdir(DIR_DNA):
            rna_families = sorted([d for d in os.listdir(DIR_DNA) if isdir(join(DIR_DNA, d))])
        else:
            print(f"ERROR: DNA directory not found: {DIR_DNA}")
            return

    print(f"Processing {len(rna_families)} RNA families...")
    print()

    # Collect all results
    all_results = []

    for rna in rna_families:
        print(f"Processing {rna}...")
        results = extract_all_seeds(rna)

        if results:
            row = {'RNA_family': rna}

            # Add DNA seed
            if 'DNA' in results:
                seed, loglh = results['DNA']
                row['DNA_seed'] = seed
                row['DNA_loglh'] = f"{loglh:.2f}"
                print(f"  DNA: seed={seed}, logLH={loglh:.2f}")

            # Add RNA model seeds
            for model in RNA_MODELS:
                if model in results:
                    seed, loglh = results[model]
                    row[f'{model}_seed'] = seed
                    row[f'{model}_loglh'] = f"{loglh:.2f}"

            all_results.append(row)
        else:
            print(f"  WARNING: No log files found")

    # Write CSV output
    if all_results:
        # Determine all columns
        columns = ['RNA_family', 'DNA_seed', 'DNA_loglh']
        for model in RNA_MODELS:
            columns.extend([f'{model}_seed', f'{model}_loglh'])

        with open(args.output, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(all_results)

        print()
        print(f"Results saved to: {args.output}")
        print()

        # Also print a simple summary for bash script usage
        print("=" * 60)
        print("Summary (for bootstrap script):")
        print("=" * 60)
        for row in all_results:
            rna = row['RNA_family']
            dna_seed = row.get('DNA_seed', '??')
            print(f"{rna}: DNA_seed={dna_seed}")
            for model in RNA_MODELS:
                model_seed = row.get(f'{model}_seed', '??')
                if model_seed != '??':
                    print(f"  {model}: seed={model_seed}")


if __name__ == '__main__':
    main()