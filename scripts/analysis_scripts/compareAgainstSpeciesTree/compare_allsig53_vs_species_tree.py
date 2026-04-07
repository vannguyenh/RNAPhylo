#!/usr/bin/env python3
"""
Compare Rfam, DNA, and RNA inferred trees against NCBI Taxonomy species trees
for the 53 RNA families significant across ALL 14 AU-test models.

Outputs:
  - Species trees with species name labels (for visualisation)
  - Rfam, DNA, RNA trees with species name labels (for visualisation)
  - CSV with nRF distances between each tree type and species tree
  - Summary statistics

Usage:
    python compare_allsig53_vs_species_tree.py
    python compare_allsig53_vs_species_tree.py --rna RF00740   # single family test
"""

import argparse
import logging
import os
import re
import sys
import types
from os.path import join, isdir, isfile

import numpy as np
import pandas as pd
import dendropy

# ete3 compatibility shim for Python 3.13+ (cgi module removed)
if sys.version_info >= (3, 13):
    sys.modules.setdefault('cgi', types.ModuleType('cgi'))

from ete3 import NCBITaxa

# ---------------------------------------------------------------------------
# Import from sibling scripts (same directory)
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from compare_vs_species_tree import (
    load_significance_table,
    fetch_rfam_tree_taxid_mapping,
    resolve_leaf_taxids,
    get_ncbi_species_tree,
    prune_duplicate_taxa,
    relabel_tree_to_taxids,
    compute_rf_iqtree,
    parse_newick,
    load_tree_pair,
    build_rna_type_map,
    DATASET_CONFIGS,
    SEED_MODELS,
    FULL_MODELS,
    MIN_UNIQUE_TAXA,
)
from getAllTreesforComparison import (
    clean_rfam_tree_for_parsing,
)

# Local Rfam seed trees (no API download needed)
RFAM_SEED_TREE_DIR = "/Users/u7875558/RNAPhylo/Rfam.seed_tree"

logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
log = logging.getLogger(__name__)


# =============================================================================
# FAMILY SELECTION
# =============================================================================

def get_all_significant_families(dataset='seed'):
    """Return sorted list of families significant in ALL models."""
    sig_table = load_significance_table(dataset)
    all_sig = sig_table[sig_table.all(axis=1)]
    families = sorted(all_sig.index.tolist())
    log.info(f"Found {len(families)} families significant in all "
             f"{len(sig_table.columns)} models")
    return families


def get_any_significant_families(dataset='seed'):
    """Return sorted list of families significant in ANY model."""
    sig_table = load_significance_table(dataset)
    any_sig = sig_table[sig_table.any(axis=1)]
    families = sorted(any_sig.index.tolist())
    log.info(f"Found {len(families)} families significant in at least one model")
    return families


def get_all_families_with_trees(dataset='seed'):
    """Return sorted list of all families that have inferred trees."""
    cfg = DATASET_CONFIGS[dataset]
    dir_au = join(cfg['working_dir'], 'outputs', cfg['au_folder'])
    models = SEED_MODELS if dataset == 'seed' else FULL_MODELS
    # Use first model directory to find families
    families = set()
    for model in models:
        model_dir = join(dir_au, model)
        if not isdir(model_dir):
            continue
        for rna in os.listdir(model_dir):
            tree_file = join(model_dir, rna, f"{rna}.highestLH.trees")
            if isfile(tree_file):
                families.add(rna)
    families = sorted(families)
    log.info(f"Found {len(families)} families with inferred trees")
    return families


# =============================================================================
# RFAM TREE HELPERS
# =============================================================================

def clean_seed_tree_for_parsing(raw_text):
    """
    Clean local Rfam seed tree text so dendropy can parse it.

    Local seed tree labels look like:
      _URS000075E9DB_9796/1-75_Equus_caballus_(horse)[9796].1:0.0
      _X06576.1/4-117_Caulobacter_spinosum[80].1:0.10333

    Problem: (horse) uses parentheses which dendropy interprets as subtrees.
    Solution: remove species name, common name, and [taxid].N after coordinates.
    Result:   _URS000075E9DB_9796/1-75:0.0
    """
    # Remove _Species_name_(common)[taxid].N after /start-end coordinates
    # Use [^[\]]* instead of [^:,]* to also match ':' in serotype names
    # (e.g., E. coli O157:H7) that would otherwise confuse the Newick parser
    cleaned = re.sub(r'(/\d+-\d+)[^[\]]*\[\d+\]\.\d*', r'\1', raw_text)
    # Remove trailing % or whitespace
    cleaned = cleaned.rstrip().rstrip('%').rstrip()
    return cleaned


def load_and_parse_rfam_tree(rfam_id):
    """
    Load Rfam seed tree from local file, clean for dendropy, and parse.

    Reads from RFAM_SEED_TREE_DIR/{rfam_id}.seed_tree.
    Returns: (dendropy.Tree, raw_text) or (None, None) on failure.
    """
    tree_path = join(RFAM_SEED_TREE_DIR, f"{rfam_id}.seed_tree")
    if not isfile(tree_path):
        log.warning(f"  Rfam seed tree not found: {tree_path}")
        return None, None

    with open(tree_path) as f:
        raw_text = f.read()

    cleaned = clean_seed_tree_for_parsing(raw_text)

    try:
        tree = dendropy.Tree.get(data=cleaned, schema="newick")
    except Exception as e:
        log.warning(f"  Failed to parse Rfam tree for {rfam_id}: {e}")
        return None, raw_text

    return tree, raw_text


def build_rfam_leaf_to_taxid(rfam_tree, rfam_taxid_map):
    """
    Map cleaned Rfam dendropy leaf labels back to taxon IDs.

    After clean_rfam_tree_for_parsing(), leaf labels look like:
      'URS000075E9DB 9796/1-75 Equus caballus'  (dendropy converts _ to space)
    We extract the accession/coords portion and look it up in rfam_taxid_map.

    Returns: dict {leaf_label: taxon_id}
    """
    label_to_taxid = {}
    for leaf in rfam_tree.leaf_node_iter():
        if leaf.taxon is None:
            continue
        label = leaf.taxon.label
        label_us = label.replace(' ', '_')

        # Try extracting accession/coords
        match = re.search(r'([A-Za-z0-9_.]+/\d+-\d+)', label_us)
        if match:
            acc = match.group(1).lstrip('_')
            if acc in rfam_taxid_map:
                label_to_taxid[label] = rfam_taxid_map[acc]
                continue

        # Fallback: try matching the first part before space/underscore
        for rfam_acc, taxid in rfam_taxid_map.items():
            if label_us.split('/')[0].lstrip('_') in rfam_acc:
                label_to_taxid[label] = taxid
                break

    return label_to_taxid


def relabel_tree_to_species_names(tree, ncbi):
    """
    Clone tree and replace taxon-ID labels with species names.
    If duplicate species names exist, append _TAXID to disambiguate.

    Returns: new dendropy.Tree with species name labels.
    """
    tree_copy = tree.clone(depth=2)

    # Collect all taxon IDs
    taxids = []
    for leaf in tree_copy.leaf_node_iter():
        if leaf.taxon is None:
            continue
        try:
            taxids.append(int(leaf.taxon.label))
        except (ValueError, TypeError):
            pass

    if not taxids:
        return tree_copy

    taxid_to_name = ncbi.get_taxid_translator(taxids)

    # Check for duplicate species names
    name_counts = {}
    for tid in taxids:
        name = taxid_to_name.get(tid, str(tid))
        name_counts[name] = name_counts.get(name, 0) + 1

    for leaf in tree_copy.leaf_node_iter():
        if leaf.taxon is None:
            continue
        try:
            tid = int(leaf.taxon.label)
        except (ValueError, TypeError):
            continue
        name = taxid_to_name.get(tid, str(tid)).replace(' ', '_')
        if name_counts.get(taxid_to_name.get(tid, str(tid)), 1) > 1:
            name = f"{name}_{tid}"
        leaf.taxon.label = name

    return tree_copy


# =============================================================================
# PER-FAMILY PROCESSING
# =============================================================================

def process_family_allsig(rfam_id, models, dir_au, cache_dir, ncbi,
                          trees_dir):
    """
    Process one family across all models.

    Returns: (list_of_result_dicts, skip_reason_or_None)
    Each result dict has: RNA_family, model, nRF_rfam, nRF_dna, nRF_rna, etc.
    """
    log.info(f"Processing {rfam_id}...")

    family_trees_dir = join(trees_dir, rfam_id)
    os.makedirs(family_trees_dir, exist_ok=True)

    # ----- Rfam tree (optional — does not block DNA/RNA comparison) -----
    rfam_result = None
    n_rfam_leaves = 0
    n_rfam_unique = 0
    rfam_tree_original = None  # keep original labels for direct tree comparison

    rfam_tree, rfam_raw = load_and_parse_rfam_tree(rfam_id)
    if rfam_tree is not None:
        rfam_tree_original = rfam_tree.clone(depth=2)  # before any relabeling
        rfam_taxid_map = fetch_rfam_tree_taxid_mapping(rfam_id, cache_dir)
        if rfam_taxid_map:
            rfam_leaf_to_taxid = build_rfam_leaf_to_taxid(rfam_tree, rfam_taxid_map)
            n_rfam_leaves = len([l for l in rfam_tree.leaf_node_iter() if l.taxon])
            n_rfam_mapped = len(rfam_leaf_to_taxid)

            if n_rfam_mapped >= MIN_UNIQUE_TAXA:
                rfam_tree, rfam_dups, rfam_unique_taxids = prune_duplicate_taxa(
                    rfam_tree, rfam_leaf_to_taxid
                )
                n_rfam_unique = len(rfam_unique_taxids)

                if n_rfam_unique >= MIN_UNIQUE_TAXA:
                    relabel_tree_to_taxids(rfam_tree, rfam_leaf_to_taxid)
                    rfam_species_tree, rfam_polytomies, rfam_merged = get_ncbi_species_tree(
                        rfam_unique_taxids, ncbi
                    )
                    if rfam_species_tree is not None:
                        # Apply merged taxid translations (skip if target already exists)
                        if rfam_merged:
                            existing = {l.taxon.label for l in rfam_tree.leaf_node_iter() if l.taxon}
                            for leaf in rfam_tree.leaf_node_iter():
                                if leaf.taxon and leaf.taxon.label in [str(k) for k in rfam_merged]:
                                    new_label = str(rfam_merged[int(leaf.taxon.label)])
                                    if new_label not in existing:
                                        existing.discard(leaf.taxon.label)
                                        leaf.taxon.label = new_label
                                        existing.add(new_label)

                        rfam_sp_copy = rfam_species_tree.clone(depth=2)
                        rfam_copy = rfam_tree.clone(depth=2)
                        rfam_result = compute_rf_and_ic_simple(rfam_sp_copy, rfam_copy)

                        # Save Rfam trees with species name labels
                        rfam_sp_named = relabel_tree_to_species_names(rfam_species_tree, ncbi)
                        save_newick(rfam_sp_named,
                                    join(family_trees_dir, f"{rfam_id}_species_tree_rfam.newick"))
                        rfam_named = relabel_tree_to_species_names(rfam_tree, ncbi)
                        save_newick(rfam_named,
                                    join(family_trees_dir, f"{rfam_id}_rfam_tree.newick"))

    if rfam_result is None:
        log.info(f"  Rfam tree comparison not possible (too few unique species)")

    # ----- DNA / RNA trees (per model) -----
    results = []
    dna_result_cache = None
    dna_vs_rfam_cache = None
    dna_saved = False

    for model in models:
        tree_file = join(dir_au, model, rfam_id, f"{rfam_id}.highestLH.trees")
        if not isfile(tree_file):
            log.debug(f"  {model}: no tree file")
            continue

        dna_newick, rna_newick = load_tree_pair(tree_file)
        if dna_newick is None:
            continue

        try:
            dna_tree = parse_newick(dna_newick)
            rna_tree = parse_newick(rna_newick)
        except Exception as e:
            log.debug(f"  {model}: parse error: {e}")
            continue

        n_dna_leaves = len([l for l in dna_tree.leaf_node_iter() if l.taxon])

        # --- Step 2: Direct tree-to-tree comparisons (no species collapsing) ---
        # RNA vs DNA
        rna_vs_dna = compute_rf_direct(rna_tree, dna_tree)

        # DNA vs Rfam (cached — same DNA tree across models)
        if dna_vs_rfam_cache is None and rfam_tree_original is not None:
            dna_vs_rfam_cache = compute_rf_direct(dna_tree, rfam_tree_original)

        # RNA vs Rfam
        rna_vs_rfam = None
        if rfam_tree_original is not None:
            rna_vs_rfam = compute_rf_direct(rna_tree, rfam_tree_original)

        # --- Step 1: Species tree comparison (requires collapsing to unique species) ---
        label_to_taxid, label_format = resolve_leaf_taxids(
            dna_tree, rfam_id, cache_dir, ncbi
        )

        sp_rfam_result = rfam_result  # from Rfam section above
        sp_dna_result = None
        sp_rna_result = None
        n_unique_taxa = 0

        if len(label_to_taxid) >= MIN_UNIQUE_TAXA:
            # Re-parse trees for species collapsing (originals used above for direct RF)
            dna_sp = parse_newick(dna_newick)
            rna_sp = parse_newick(rna_newick)

            dna_sp, dna_dups, dna_unique = prune_duplicate_taxa(dna_sp, label_to_taxid)
            rna_sp, rna_dups, rna_unique = prune_duplicate_taxa(rna_sp, label_to_taxid)
            common_taxids = dna_unique & rna_unique

            if len(common_taxids) >= MIN_UNIQUE_TAXA:
                relabel_tree_to_taxids(dna_sp, label_to_taxid)
                relabel_tree_to_taxids(rna_sp, label_to_taxid)

                dna_species_tree, dna_polytomies, dna_merged = get_ncbi_species_tree(
                    common_taxids, ncbi
                )

                if dna_species_tree is not None:
                    # Apply merged taxid translations (skip if target already exists)
                    if dna_merged:
                        for tree in [dna_sp, rna_sp]:
                            existing = {l.taxon.label for l in tree.leaf_node_iter() if l.taxon}
                            for leaf in tree.leaf_node_iter():
                                if leaf.taxon and leaf.taxon.label in [str(k) for k in dna_merged]:
                                    new_label = str(dna_merged[int(leaf.taxon.label)])
                                    if new_label not in existing:
                                        existing.discard(leaf.taxon.label)
                                        leaf.taxon.label = new_label
                                        existing.add(new_label)

                    n_unique_taxa = len(common_taxids)

                    if dna_result_cache is None:
                        dna_result_cache = compute_rf_and_ic_simple(
                            dna_species_tree.clone(2), dna_sp.clone(2))

                    sp_dna_result = dna_result_cache
                    sp_rna_result = compute_rf_and_ic_simple(
                        dna_species_tree.clone(2), rna_sp.clone(2))

                    # Save trees with species name labels (once)
                    if not dna_saved:
                        dna_sp_named = relabel_tree_to_species_names(
                            dna_species_tree, ncbi)
                        save_newick(dna_sp_named, join(
                            family_trees_dir, f"{rfam_id}_species_tree_dna.newick"))
                        dna_named = relabel_tree_to_species_names(dna_sp, ncbi)
                        save_newick(dna_named, join(
                            family_trees_dir, f"{rfam_id}_dna_tree.newick"))
                        dna_saved = True

                    rna_named = relabel_tree_to_species_names(rna_sp, ncbi)
                    save_newick(rna_named, join(
                        family_trees_dir, f"{rfam_id}_rna_{model}_tree.newick"))

        # --- Build result row ---
        row = {
            'RNA_family': rfam_id,
            'model': model,
            'n_rfam_leaves': n_rfam_leaves,
            'n_rfam_unique_taxa': n_rfam_unique,
            'n_dna_leaves': n_dna_leaves,
            'n_dna_unique_taxa': n_unique_taxa,
            'label_format': label_format,
        }

        # Step 1: vs species tree
        for prefix, res in [('rfam', sp_rfam_result),
                            ('dna', sp_dna_result),
                            ('rna', sp_rna_result)]:
            if res is not None:
                row[f'RF_{prefix}_vs_sp'] = res['RF']
                row[f'nRF_{prefix}_vs_sp'] = res['nRF']
                row[f'IC_{prefix}_vs_sp'] = res['IC']
                row[f'nIC_{prefix}_vs_sp'] = res['nIC']
            else:
                row[f'RF_{prefix}_vs_sp'] = np.nan
                row[f'nRF_{prefix}_vs_sp'] = np.nan
                row[f'IC_{prefix}_vs_sp'] = np.nan
                row[f'nIC_{prefix}_vs_sp'] = np.nan

        # Step 2: direct tree-to-tree
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

    if not results:
        return [], "no_model_succeeded"

    return results, None


def compute_rf_and_ic_simple(species_tree, inferred_tree):
    """
    Compute RF and IC distances between species tree and inferred tree.

    RF: standard Robinson-Foulds (via IQ-TREE).
    IC: incompatible splits count — only counts bipartitions in the inferred
        tree that CONTRADICT the species tree. Handles polytomies correctly.

    Prunes to common taxa first.
    Returns: dict with RF, nRF, IC, nIC, n_taxa, or None.
    """
    sp_labels = {leaf.taxon.label for leaf in species_tree.leaf_node_iter()
                 if leaf.taxon}
    inf_labels = {leaf.taxon.label for leaf in inferred_tree.leaf_node_iter()
                  if leaf.taxon}

    common = sp_labels & inf_labels
    n_taxa = len(common)

    if n_taxa < MIN_UNIQUE_TAXA:
        return None

    # Re-parse into shared TaxonNamespace
    shared_tns = dendropy.TaxonNamespace()
    sp_newick = species_tree.as_string(schema="newick")
    inf_newick = inferred_tree.as_string(schema="newick")

    species_tree = dendropy.Tree.get(
        data=sp_newick, schema="newick", taxon_namespace=shared_tns)
    inferred_tree = dendropy.Tree.get(
        data=inf_newick, schema="newick", taxon_namespace=shared_tns)

    # Prune to common set
    if sp_labels != inf_labels:
        sp_extra = sp_labels - common
        inf_extra = inf_labels - common
        if sp_extra:
            to_remove = [t for t in shared_tns if t.label in sp_extra]
            species_tree.prune_taxa(to_remove)
        if inf_extra:
            to_remove = [t for t in shared_tns if t.label in inf_extra]
            inferred_tree.prune_taxa(to_remove)
        n_taxa = len(common)

    if n_taxa < MIN_UNIQUE_TAXA:
        return None

    max_rf = 2 * (n_taxa - 3)
    if max_rf <= 0:
        return None

    # RF distance via IQ-TREE
    rf = compute_rf_iqtree(species_tree, inferred_tree)
    if rf is None:
        return None
    nrf = rf / max_rf if max_rf > 0 else 0.0

    # IC: incompatible splits — handles polytomies correctly
    species_tree.encode_bipartitions()
    inferred_tree.encode_bipartitions()

    sp_bipartitions = set()
    for nd in species_tree.internal_nodes():
        if nd.bipartition is not None:
            bp = nd.bipartition.split_bitmask
            complement = nd.bipartition.leafset_bitmask ^ bp
            sp_bipartitions.add(frozenset([bp, complement]))

    ic = 0
    for nd in inferred_tree.internal_nodes():
        if nd.bipartition is None:
            continue
        inf_bp = nd.bipartition.split_bitmask
        inf_complement = nd.bipartition.leafset_bitmask ^ inf_bp

        for sp_bp_pair in sp_bipartitions:
            sp_a, sp_b = sp_bp_pair
            if ((inf_bp & sp_a) and (inf_bp & sp_b) and
                    (inf_complement & sp_a) and (inf_complement & sp_b)):
                ic += 1
                break

    max_ic = n_taxa - 3  # max internal splits in binary inferred tree
    nic = ic / max_ic if max_ic > 0 else 0.0

    return {
        'n_taxa': n_taxa, 'RF': rf, 'nRF': nrf, 'max_RF': max_rf,
        'IC': ic, 'nIC': nic, 'max_IC': max_ic,
    }


def normalize_label(label):
    """Normalize a leaf label: strip whitespace, convert spaces to underscores,
    and remove Rfam numeric prefixes like '744.6_' before accessions."""
    label = label.strip().replace(' ', '_')
    # Remove Rfam numeric prefix (e.g., '744.6_AE016830.1/...' -> 'AE016830.1/...')
    label = re.sub(r'^\d+(\.\d+)?_', '', label)
    return label


def compute_rf_direct(tree1, tree2):
    """
    Compute RF distance between two gene trees directly (no species tree).

    Normalizes leaf labels, prunes to common set, then computes RF.
    Works even when trees have different numbers of leaves.
    Returns: dict with RF, nRF, n_common, or None.
    """
    # Normalize labels on clones
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

    # Re-parse into shared TaxonNamespace
    shared_tns = dendropy.TaxonNamespace()
    nwk1 = t1.as_string(schema="newick")
    nwk2 = t2.as_string(schema="newick")
    t1 = dendropy.Tree.get(data=nwk1, schema="newick", taxon_namespace=shared_tns)
    t2 = dendropy.Tree.get(data=nwk2, schema="newick", taxon_namespace=shared_tns)

    # Prune to common set
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


def save_newick(tree, path):
    """Save a dendropy tree to a Newick file."""
    newick = tree.as_string(schema="newick")
    with open(path, 'w') as f:
        f.write(newick)


# =============================================================================
# ORCHESTRATOR
# =============================================================================

FAMILY_MODES = {
    'all-significant': get_all_significant_families,
    'any-significant': get_any_significant_families,
    'all': get_all_families_with_trees,
}

OUTPUT_SUFFIXES = {
    'all-significant': '260306_allsig53_species_tree_comparison',
    'any-significant': '260306_anysig_species_tree_comparison',
    'all': '260306_all_species_tree_comparison',
}


def process_all(dataset='seed', rna_filter=None, family_mode='all-significant'):
    """Main orchestrator. Returns DataFrame of all results."""
    cfg = DATASET_CONFIGS[dataset]
    dir_au = join(cfg['working_dir'], 'outputs', cfg['au_folder'])
    output_dir = join(cfg['working_dir'], 'outputs',
                      OUTPUT_SUFFIXES[family_mode])
    cache_dir = join(cfg['working_dir'], 'outputs',
                     cfg['output_folder'], 'cache')
    trees_dir = join(output_dir, 'species_trees')
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(cache_dir, exist_ok=True)
    os.makedirs(trees_dir, exist_ok=True)

    models = SEED_MODELS if dataset == 'seed' else FULL_MODELS

    # Select families based on mode
    families = FAMILY_MODES[family_mode](dataset)

    if rna_filter:
        families = [f for f in families if f in rna_filter]

    log.info(f"Processing {len(families)} families × {len(models)} models")

    # Init NCBITaxa
    log.info("Initialising NCBITaxa...")
    ncbi = NCBITaxa()

    all_results = []
    skipped = []

    for i, rfam_id in enumerate(families):
        log.info(f"[{i+1}/{len(families)}] {rfam_id}")

        results, skip_reason = process_family_allsig(
            rfam_id, models, dir_au, cache_dir, ncbi, trees_dir
        )

        if results:
            all_results.extend(results)
        if skip_reason:
            skipped.append({'RNA_family': rfam_id, 'reason': skip_reason})
            log.warning(f"  Skipped {rfam_id}: {skip_reason}")

    log.info(f"\nDone! {len(all_results)} rows from "
             f"{len(families)} families")

    df = pd.DataFrame(all_results)

    # Drop families where all nRF columns are NaN (too few unique species
    # after NCBI taxonomy merges deprecated taxids)
    if not df.empty:
        sp_cols = [c for c in ['nRF_rfam_vs_sp', 'nRF_dna_vs_sp', 'nRF_rna_vs_sp'] if c in df.columns]
        direct_cols = [c for c in ['nRF_rna_vs_dna', 'nRF_dna_vs_rfam', 'nRF_rna_vs_rfam'] if c in df.columns]
        all_nrf_cols = sp_cols + direct_cols
        all_nan_mask = df[all_nrf_cols].isna().all(axis=1)
        dropped_families = df.loc[all_nan_mask, 'RNA_family'].unique()
        if len(dropped_families) > 0:
            df = df[~all_nan_mask]
            for fam in dropped_families:
                skipped.append({'RNA_family': fam, 'reason': 'too_few_species_after_taxid_merge'})
                log.warning(f"  Skipped {fam}: too few unique species after NCBI taxid merge")

    # Add RNA type
    type_map = build_rna_type_map()
    if type_map and not df.empty:
        df['rna_type'] = df['RNA_family'].map(type_map).fillna('Other')

    # Save results
    if not df.empty:
        csv_path = join(output_dir, "nRF_comparison.csv")
        df.to_csv(csv_path, index=False)
        log.info(f"Results saved to {csv_path}")

    # Save skipped
    if skipped:
        skip_df = pd.DataFrame(skipped)
        skip_path = join(output_dir, "skipped_families.csv")

        skip_df.to_csv(skip_path, index=False)
        log.info(f"Skipped families: {skip_path}")

    return df, output_dir


# =============================================================================
# SUMMARY
# =============================================================================

def compute_summary(df):
    """Per-model summary of nRF distances."""
    if df.empty:
        return pd.DataFrame()

    summaries = []
    for model, grp in df.groupby('model'):
        row = {
            'model': model,
            'n_families': len(grp),
        }

        # Step 1: vs species tree
        for col in ['nRF_rfam_vs_sp', 'nRF_dna_vs_sp', 'nRF_rna_vs_sp']:
            if col in grp.columns:
                row[f'mean_{col}'] = grp[col].mean()
                row[f'median_{col}'] = grp[col].median()

        # Step 2: direct tree-to-tree
        for col in ['nRF_rna_vs_dna', 'nRF_dna_vs_rfam', 'nRF_rna_vs_rfam']:
            if col in grp.columns:
                row[f'mean_{col}'] = grp[col].mean()
                row[f'median_{col}'] = grp[col].median()

        summaries.append(row)

    return pd.DataFrame(summaries)


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Compare Rfam/DNA/RNA trees vs species trees.'
    )
    parser.add_argument('--dataset', choices=['seed', 'full'], default='seed')
    parser.add_argument('--families', choices=['all-significant', 'any-significant', 'all'],
                        default='all-significant',
                        help='Which families to process: '
                             'all-significant (AU-test significant in ALL models), '
                             'any-significant (significant in at least one model), '
                             'all (every family with inferred trees)')
    parser.add_argument('--rna', type=str, nargs='*', default=None,
                        help='Process only these families (e.g., RF00740)')
    args = parser.parse_args()

    rna_filter = set(args.rna) if args.rna else None

    df, output_dir = process_all(args.dataset, rna_filter, args.families)

    if df.empty:
        log.warning("No results to report.")
        return

    # Summary
    summary = compute_summary(df)
    summary_path = join(output_dir, "summary.csv")
    summary.to_csv(summary_path, index=False)

    print(f"\n{'=' * 70}")
    print(f"Tree Comparison — {args.families} ({df.RNA_family.nunique()} families)")
    print(f"{'=' * 70}")
    for _, row in summary.iterrows():
        print(f"\n  {row['model']}:")
        print(f"    Families: {row['n_families']}")
        if 'mean_nRF_dna_vs_sp' in row:
            print(f"    vs Species tree — Rfam: {row.get('mean_nRF_rfam_vs_sp', float('nan')):.4f}  "
                  f"DNA: {row.get('mean_nRF_dna_vs_sp', float('nan')):.4f}  "
                  f"RNA: {row.get('mean_nRF_rna_vs_sp', float('nan')):.4f}")
        if 'mean_nRF_rna_vs_dna' in row:
            print(f"    Direct — RNA vs DNA: {row.get('mean_nRF_rna_vs_dna', float('nan')):.4f}  "
                  f"DNA vs Rfam: {row.get('mean_nRF_dna_vs_rfam', float('nan')):.4f}  "
                  f"RNA vs Rfam: {row.get('mean_nRF_rna_vs_rfam', float('nan')):.4f}")

    print(f"\nTrees saved in: {join(output_dir, 'species_trees')}")
    print(f"Results CSV: {join(output_dir, 'nRF_comparison.csv')}")


if __name__ == '__main__':
    main()
