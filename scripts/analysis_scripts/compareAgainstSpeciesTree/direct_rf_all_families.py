#!/usr/bin/env python3
"""
Fast direct tree-to-tree nRF comparison for ALL families (no species tree needed).

Compares RNA vs DNA, DNA vs Rfam, and RNA vs Rfam trees directly.
No API calls, no taxid resolution, no species collapsing.

Usage:
    python direct_rf_all_families.py
    python direct_rf_all_families.py --families all-significant
    python direct_rf_all_families.py --rna RF00740
"""

import argparse
import logging
import os
import sys
import types
from os.path import join, isdir, isfile

import numpy as np
import pandas as pd
import dendropy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from compare_vs_species_tree import (
    load_significance_table,
    parse_newick,
    load_tree_pair,
    compute_rf_iqtree,
    build_rna_type_map,
    DATASET_CONFIGS,
    SEED_MODELS,
    FULL_MODELS,
    MIN_UNIQUE_TAXA,
)
from compare_allsig53_vs_species_tree import (
    load_and_parse_rfam_tree,
    normalize_label,
    get_all_significant_families,
    get_any_significant_families,
    get_all_families_with_trees,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
log = logging.getLogger(__name__)

RFAM_SEED_TREE_DIR = "/Users/u7875558/RNAPhylo/Rfam.seed_tree"


def compute_rf_direct_fast(tree1, tree2):
    """
    Compute RF between two trees. Normalize labels, prune to common set.
    Returns dict with RF, nRF, n_common or None.
    """
    t1 = tree1.clone(depth=2)
    t2 = tree2.clone(depth=2)
    for t in [t1, t2]:
        for leaf in t.leaf_node_iter():
            if leaf.taxon:
                leaf.taxon.label = normalize_label(leaf.taxon.label)

    labels1 = {l.taxon.label for l in t1.leaf_node_iter() if l.taxon}
    labels2 = {l.taxon.label for l in t2.leaf_node_iter() if l.taxon}
    common = labels1 & labels2
    n_common = len(common)

    if n_common < MIN_UNIQUE_TAXA:
        return None

    shared_tns = dendropy.TaxonNamespace()
    nwk1 = t1.as_string(schema="newick")
    nwk2 = t2.as_string(schema="newick")
    t1 = dendropy.Tree.get(data=nwk1, schema="newick", taxon_namespace=shared_tns)
    t2 = dendropy.Tree.get(data=nwk2, schema="newick", taxon_namespace=shared_tns)

    extra1 = labels1 - common
    extra2 = labels2 - common
    if extra1:
        t1.prune_taxa([t for t in shared_tns if t.label in extra1])
    if extra2:
        t2.prune_taxa([t for t in shared_tns if t.label in extra2])

    max_rf = 2 * (n_common - 3)
    if max_rf <= 0:
        return None

    rf = compute_rf_iqtree(t1, t2)
    if rf is None:
        return None

    nrf = rf / max_rf if max_rf > 0 else 0.0
    return {'n_common': n_common, 'RF': rf, 'nRF': nrf, 'max_RF': max_rf}


def process_family(rfam_id, models, dir_au):
    """Process one family: direct RF comparisons only."""
    # Load Rfam tree
    rfam_tree = None
    rfam_path = join(RFAM_SEED_TREE_DIR, f"{rfam_id}.seed_tree")
    if isfile(rfam_path):
        rfam_tree_parsed, _ = load_and_parse_rfam_tree(rfam_id)
        if rfam_tree_parsed is not None:
            rfam_tree = rfam_tree_parsed

    results = []
    dna_vs_rfam_cache = None

    for model in models:
        tree_file = join(dir_au, model, rfam_id, f"{rfam_id}.highestLH.trees")
        if not isfile(tree_file):
            continue

        dna_newick, rna_newick = load_tree_pair(tree_file)
        if dna_newick is None:
            continue

        try:
            dna_tree = parse_newick(dna_newick)
            rna_tree = parse_newick(rna_newick)
        except Exception:
            continue

        n_dna_leaves = len([l for l in dna_tree.leaf_node_iter() if l.taxon])

        # RNA vs DNA
        rna_vs_dna = compute_rf_direct_fast(rna_tree, dna_tree)

        # DNA vs Rfam (cached)
        if dna_vs_rfam_cache is None and rfam_tree is not None:
            dna_vs_rfam_cache = compute_rf_direct_fast(dna_tree, rfam_tree)

        # RNA vs Rfam
        rna_vs_rfam = None
        if rfam_tree is not None:
            rna_vs_rfam = compute_rf_direct_fast(rna_tree, rfam_tree)

        row = {
            'RNA_family': rfam_id,
            'model': model,
            'n_dna_leaves': n_dna_leaves,
        }

        for name, res in [('rna_vs_dna', rna_vs_dna),
                          ('dna_vs_rfam', dna_vs_rfam_cache),
                          ('rna_vs_rfam', rna_vs_rfam)]:
            if res is not None:
                row[f'RF_{name}'] = res['RF']
                row[f'nRF_{name}'] = res['nRF']
                row[f'n_common_{name}'] = res['n_common']
            else:
                row[f'RF_{name}'] = np.nan
                row[f'nRF_{name}'] = np.nan
                row[f'n_common_{name}'] = np.nan

        results.append(row)

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Fast direct tree-to-tree nRF comparison (no species tree).')
    parser.add_argument('--dataset', choices=['seed', 'full'], default='seed')
    parser.add_argument('--families', choices=['all-significant', 'any-significant', 'all'],
                        default='all')
    parser.add_argument('--rna', type=str, nargs='*', default=None)
    args = parser.parse_args()

    cfg = DATASET_CONFIGS[args.dataset]
    dir_au = join(cfg['working_dir'], 'outputs', cfg['au_folder'])
    models = SEED_MODELS if args.dataset == 'seed' else FULL_MODELS

    family_funcs = {
        'all-significant': get_all_significant_families,
        'any-significant': get_any_significant_families,
        'all': get_all_families_with_trees,
    }
    families = family_funcs[args.families](args.dataset)

    if args.rna:
        families = [f for f in families if f in set(args.rna)]

    log.info(f"Processing {len(families)} families × {len(models)} models")

    all_results = []
    for i, rfam_id in enumerate(families):
        if (i + 1) % 100 == 0:
            log.info(f"  [{i+1}/{len(families)}]")

        results = process_family(rfam_id, models, dir_au)
        all_results.extend(results)

    df = pd.DataFrame(all_results)

    # Add RNA type
    type_map = build_rna_type_map()
    if type_map and not df.empty:
        df['rna_type'] = df['RNA_family'].map(type_map).fillna('Other')

    # Save
    output_dir = join(cfg['working_dir'], 'outputs',
                      f'260306_direct_rf_{args.families.replace("-", "_")}')
    os.makedirs(output_dir, exist_ok=True)
    csv_path = join(output_dir, 'direct_nRF_comparison.csv')
    df.to_csv(csv_path, index=False)

    log.info(f"\nDone! {len(df)} rows, {df.RNA_family.nunique()} families")
    log.info(f"Results: {csv_path}")

    # Quick summary
    for col in ['nRF_rna_vs_dna', 'nRF_dna_vs_rfam', 'nRF_rna_vs_rfam']:
        if col in df.columns:
            valid = df[col].dropna()
            log.info(f"  {col}: mean={valid.mean():.4f}, median={valid.median():.4f}, "
                     f"n={len(valid)}, NaN={df[col].isna().sum()}")

    print(f"\nResults saved to: {csv_path}")
    print(f"Families: {df.RNA_family.nunique()}, Rows: {len(df)}")


if __name__ == '__main__':
    main()
