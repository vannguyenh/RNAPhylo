#!/usr/bin/env python3
"""
Compare inferred phylogenetic trees against NCBI Taxonomy species trees.

For each AU-test significant RNA family, computes the Robinson-Foulds (RF)
distance and incompatible-splits (IC) distance between the inferred tree
(DNA-only or RNA model) and the NCBI Taxonomy species tree.

The key question: do RNA substitution models produce trees that are closer
to the true species phylogeny than DNA-only models?

Usage:
    python compare_vs_species_tree.py --dataset seed
    python compare_vs_species_tree.py --dataset full
    python compare_vs_species_tree.py --dataset seed --rna RF00740  # single family test
"""

import argparse
import json
import logging
import os
import re
import subprocess
import sys
import tempfile
import time
import types
from collections import Counter
from os.path import join, isdir, isfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import dendropy

# ete3 compatibility shim for Python 3.13+ (cgi module removed)
if sys.version_info >= (3, 13):
    sys.modules.setdefault('cgi', types.ModuleType('cgi'))

from ete3 import NCBITaxa
import requests

# =============================================================================
# CONFIGURATION
# =============================================================================

DATASET_CONFIGS = {
    'seed': {
        'working_dir': '/Users/u7875558/RNAPhylo/seedAlignment_AllModels',
        'au_folder': '260208_AU_Test_RAxML_RNAmodels',
        'sig_table_csv': '260210_AU_significant_table.csv',
        'output_folder': '260305_species_tree_comparison',
    },
    'full': {
        'working_dir': '/Users/u7875558/RNAPhylo/fullAlignment',
        'au_folder': '260208_AU_Test_RAxML_RNAmodels',
        'sig_table_csv': '260210_AU_significant_table.csv',
        'output_folder': '260305_species_tree_comparison',
    },
}

SEED_MODELS = ['S16', 'S16A', 'S16B', 'S7A', 'S7B', 'S7C', 'S7D', 'S7E',
               'S7F', 'S6A', 'S6B', 'S6C', 'S6D', 'S6E']
FULL_MODELS = ['S16', 'S6A', 'S7A']

MIN_UNIQUE_TAXA = 4   # Minimum unique species for meaningful RF comparison
ALPHA = 0.05
RFAM_API_DELAY = 0.5  # seconds between Rfam API calls
IQTREE_PATH = "/Users/u7875558/tools/build-iqtree3/iqtree3"

# RNA type classification (reused from plot_nRF_by_RNAtypes.py)
RFAM_FAMILY_TXT = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/inputs/family.txt"

RNA_TYPE_ORDER = [
    "miRNA", "sRNA", "snRNA / snoRNA", "lncRNA", "rRNA",
    "Antisense", "Cis-reg", "CRISPR", "Antitoxin", "Intron", "tRNA", "Other",
]

TYPE_MAP = {
    "Gene; miRNA;":                      "miRNA",
    "Gene; sRNA;":                       "sRNA",
    "Gene; snRNA; splicing;":            "snRNA / snoRNA",
    "Gene; snRNA;":                      "snRNA / snoRNA",
    "Gene; snRNA; snoRNA; CD-box;":      "snRNA / snoRNA",
    "Gene; snRNA; snoRNA; HACA-box;":    "snRNA / snoRNA",
    "Gene; snRNA; snoRNA; scaRNA;":      "snRNA / snoRNA",
    "Gene; snRNA; snoRNA;":              "snRNA / snoRNA",
    "Gene; lncRNA;":                     "lncRNA",
    "Gene; rRNA;":                       "rRNA",
    "Gene; antisense;":                  "Antisense",
    "Cis-reg;":                          "Cis-reg",
    "Cis-reg; leader;":                  "Cis-reg",
    "Cis-reg; riboswitch;":             "Cis-reg",
    "Cis-reg; thermoregulator;":        "Cis-reg",
    "Cis-reg; IRES;":                   "Cis-reg",
    "Cis-reg; frameshift_element;":     "Cis-reg",
    "Gene; CRISPR;":                     "CRISPR",
    "Gene; antitoxin;":                  "Antitoxin",
    "Intron;":                           "Intron",
    "Gene; tRNA;":                       "tRNA",
    "Gene; ribozyme;":                   "Other",
    "Gene;":                             "Other",
}

logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
log = logging.getLogger(__name__)


# =============================================================================
# SIGNIFICANCE TABLE LOADING (reused from figure_AUtest_significance.py)
# =============================================================================

def parse_pv_file(pv_path):
    """Parse a CONSEL .pv file and return (p_DNA, p_RNA) AU test p-values."""
    with open(pv_path) as f:
        lines = f.readlines()

    in_pv = False
    row_values = []
    for line in lines:
        stripped = line.strip()
        if stripped == '# PV:':
            in_pv = True
            continue
        if in_pv and stripped.startswith('# row:'):
            continue
        if in_pv and stripped.startswith('#') and not stripped.startswith('# row'):
            if stripped == '# SE:':
                break
            continue
        if in_pv and stripped and not stripped.startswith('#'):
            values = stripped.split()
            row_values.append(float(values[0]))

    if len(row_values) >= 2:
        return row_values[0], row_values[1]
    return None, None


def build_significance_table_from_pv(dir_au):
    """Parse all .pv files under dir_au and build boolean significance table."""
    models = sorted(
        m for m in os.listdir(dir_au)
        if isdir(join(dir_au, m)) and m.startswith('S')
    )
    records = {}
    for model in models:
        model_dir = join(dir_au, model)
        for rna in sorted(os.listdir(model_dir)):
            rna_dir = join(model_dir, rna)
            if not isdir(rna_dir):
                continue
            pv_file = join(rna_dir, f"{rna}_consel.pv")
            if not isfile(pv_file):
                continue
            p_dna, p_rna = parse_pv_file(pv_file)
            if p_dna is None:
                continue
            p_min = min(p_dna, p_rna)
            if rna not in records:
                records[rna] = {}
            records[rna][model] = p_min < ALPHA

    df = pd.DataFrame.from_dict(records, orient='index')
    df.index.name = 'RNA_family'
    df = df.reindex(sorted(df.columns), axis=1)
    df = df.sort_index()
    df = df.fillna(False)
    return df


def load_significance_table(dataset):
    """Load or build the boolean significance table for a dataset."""
    cfg = DATASET_CONFIGS[dataset]
    dir_au = join(cfg['working_dir'], 'outputs', cfg['au_folder'])
    csv_path = join(dir_au, cfg['sig_table_csv'])

    if isfile(csv_path):
        df = pd.read_csv(csv_path, index_col='RNA_family')
        df = df.apply(lambda col: col.astype(str).str.strip().eq('True'))
        log.info(f"Loaded significance table from {csv_path}")
        return df

    log.info(f"No CSV found. Parsing .pv files from {dir_au} ...")
    df = build_significance_table_from_pv(dir_au)
    df.to_csv(csv_path)
    log.info(f"Saved significance table to {csv_path}")
    return df


# =============================================================================
# TAXON ID RESOLUTION
# =============================================================================

def extract_taxids_from_urs_labels(leaf_labels):
    """
    Extract taxon IDs from URS-format labels.

    URS labels: URS0000D67E41_1188252/1-53 → taxon ID 1188252
    Returns: dict {leaf_label: taxon_id} for labels that match.
    """
    pattern = re.compile(r'[A-Z0-9]+_(\d+)/\d+-\d+')
    result = {}
    for label in leaf_labels:
        label_clean = label.replace(' ', '_')
        match = pattern.search(label_clean)
        if match:
            result[label] = int(match.group(1))
    return result


def fetch_rfam_tree_taxid_mapping(rfam_id, cache_dir):
    """
    Download Rfam tree via API and extract {accession_label: taxon_id}.

    Rfam tree labels: _URS000075E9DB_9796/1-75_Equus_caballus_{horse}[9796].1
    Extract [TAXID] from bracket notation.
    Caches result as JSON.
    """
    cache_file = join(cache_dir, f"{rfam_id}_taxid_map.json")
    if isfile(cache_file):
        with open(cache_file) as f:
            return {k: int(v) for k, v in json.load(f).items()}

    url = f"https://rfam.org/family/{rfam_id}/tree"
    log.info(f"    Downloading Rfam tree for {rfam_id}...")
    try:
        response = requests.get(url, timeout=60)
        if response.status_code != 200:
            log.warning(f"    Rfam API returned {response.status_code} for {rfam_id}")
            return {}
    except Exception as e:
        log.warning(f"    Rfam API error for {rfam_id}: {e}")
        return {}

    rfam_tree_text = response.text

    # Pattern: accession/coords ... [TAXID].1
    # Example: _URS000075E9DB_9796/1-75_Equus_caballus_{horse}[9796].1
    # Also handles: CP000002.3/2295936-2295837_Species_name_...[TAXID].1
    mapping = {}

    # Extract accession and taxon ID pairs
    # Match: accession/start-end ... [taxid]
    for match in re.finditer(
        r'([A-Za-z0-9_.]+/\d+-\d+)[^[]*\[(\d+)\]',
        rfam_tree_text
    ):
        acc_raw = match.group(1)
        taxid = int(match.group(2))
        # Clean the accession: remove leading underscore
        acc = acc_raw.lstrip('_')
        mapping[acc] = taxid

    # Save cache
    os.makedirs(cache_dir, exist_ok=True)
    with open(cache_file, 'w') as f:
        json.dump({k: v for k, v in mapping.items()}, f)

    time.sleep(RFAM_API_DELAY)
    return mapping


def resolve_leaf_taxids(tree, rfam_id, cache_dir, ncbi):
    """
    Master function: detect label format and resolve all leaf labels to taxon IDs.

    Returns: dict {leaf_label: taxon_id}
    """
    leaf_labels = [leaf.taxon.label for leaf in tree.leaf_node_iter()
                   if leaf.taxon is not None]

    # Try URS format first
    urs_map = extract_taxids_from_urs_labels(leaf_labels)
    if len(urs_map) >= len(leaf_labels) * 0.5:
        # Mostly URS format
        label_format = "URS"
        result = urs_map
    else:
        # GenBank format — need Rfam API lookup
        label_format = "GenBank"
        rfam_map = fetch_rfam_tree_taxid_mapping(rfam_id, cache_dir)

        result = {}
        for label in leaf_labels:
            label_clean = label.replace(' ', '_')
            # Try direct lookup
            if label_clean in rfam_map:
                result[label] = rfam_map[label_clean]
                continue
            # Try extracting accession part
            match = re.search(r'([A-Za-z0-9_.]+/\d+-\d+)', label_clean)
            if match:
                acc = match.group(1)
                if acc in rfam_map:
                    result[label] = rfam_map[acc]
                    continue
            # Try matching without coordinates
            for rfam_acc, taxid in rfam_map.items():
                if label_clean.split('/')[0] in rfam_acc:
                    result[label] = taxid
                    break

    return result, label_format


# =============================================================================
# TREE MANIPULATION
# =============================================================================

def load_tree_pair(tree_file):
    """Read .highestLH.trees and return (dna_newick, rna_newick)."""
    with open(tree_file) as f:
        lines = [line.strip() for line in f if line.strip()]
    if len(lines) < 2:
        return None, None
    return lines[0], lines[1]


def parse_newick(newick_str):
    """Parse a Newick string into a dendropy Tree."""
    tns = dendropy.TaxonNamespace()
    return dendropy.Tree.get(data=newick_str, schema="newick",
                             taxon_namespace=tns)


def prune_duplicate_taxa(tree, label_to_taxid):
    """
    Remove duplicate taxa (multiple sequences from same species).

    Keeps one leaf per taxon ID (the one with shortest root-to-tip distance).
    Returns: (pruned_tree, n_duplicates_removed, set_of_unique_taxids)
    """
    # Group leaves by taxon ID
    taxid_to_leaves = {}
    unmapped = []
    for leaf in tree.leaf_node_iter():
        if leaf.taxon is None:
            continue
        label = leaf.taxon.label
        if label in label_to_taxid:
            taxid = label_to_taxid[label]
            taxid_to_leaves.setdefault(taxid, []).append(leaf)
        else:
            unmapped.append(leaf)

    # For each taxon ID with duplicates, keep the one with shortest root distance
    labels_to_remove = []
    for taxid, leaves in taxid_to_leaves.items():
        if len(leaves) <= 1:
            continue
        # Compute root-to-tip distance for each leaf
        distances = []
        for leaf in leaves:
            dist = leaf.distance_from_root()
            distances.append((dist, leaf))
        distances.sort(key=lambda x: x[0])
        # Keep the first (shortest distance), remove rest
        for _, leaf in distances[1:]:
            labels_to_remove.append(leaf.taxon)

    # Also remove unmapped leaves
    for leaf in unmapped:
        if leaf.taxon:
            labels_to_remove.append(leaf.taxon)

    n_removed = len(labels_to_remove)

    if labels_to_remove:
        tree.prune_taxa(labels_to_remove)

    unique_taxids = set(taxid_to_leaves.keys())
    return tree, n_removed, unique_taxids


def relabel_tree_to_taxids(tree, label_to_taxid):
    """Replace leaf labels with string taxon IDs for comparison."""
    for leaf in tree.leaf_node_iter():
        if leaf.taxon is None:
            continue
        label = leaf.taxon.label
        if label in label_to_taxid:
            leaf.taxon.label = str(label_to_taxid[label])


def get_ncbi_species_tree(taxon_ids, ncbi):
    """
    Get NCBI taxonomy tree for given taxon IDs using ete3.

    Returns: (dendropy.Tree, n_polytomies, merged_map)
    merged_map: dict mapping old taxids to new taxids for any that were
    merged/translated by NCBITaxa (e.g. deprecated strain IDs).
    """
    if len(taxon_ids) < MIN_UNIQUE_TAXA:
        return None, 0, {}

    try:
        ete_tree = ncbi.get_topology(list(taxon_ids))
    except Exception as e:
        log.warning(f"    NCBITaxa.get_topology failed: {e}")
        return None, 0, {}

    # Detect merged taxids: species tree leaves may differ from input
    sp_leaf_ids = {int(leaf.name) for leaf in ete_tree.get_leaves()}
    merged_map = {}
    for old_id in taxon_ids:
        if old_id not in sp_leaf_ids:
            # Find which species tree leaf this old_id was merged into
            for new_id in sp_leaf_ids:
                if new_id not in taxon_ids:
                    merged_map[old_id] = new_id
                    break

    # Convert ete3 tree to Newick string, then parse with dendropy
    newick_str = ete_tree.write(format=9)  # format 9 = leaf names only

    # ete3 uses taxon IDs as node names — they should already be ints as strings
    tns = dendropy.TaxonNamespace()
    dp_tree = dendropy.Tree.get(data=newick_str, schema="newick",
                                 taxon_namespace=tns)

    # Count polytomies
    n_polytomies = sum(
        1 for node in dp_tree.internal_nodes()
        if len(node.child_nodes()) > 2
    )

    return dp_tree, n_polytomies, merged_map


# =============================================================================
# DISTANCE COMPUTATION
# =============================================================================

def compute_rf_iqtree(species_tree, inferred_tree):
    """
    Compute RF distance using IQ-TREE -rf (consistent with rf_distance_calculator.py).

    Both trees are derooted (trifurcation at root) so IQ-TREE treats them
    as unrooted. Falls back to dendropy for tiny trees where IQ-TREE fails.

    Returns: RF distance (int) or None on failure.
    """
    from dendropy.calculate import treecompare

    # Deroot both trees so IQ-TREE treats them as unrooted
    sp_copy = species_tree.clone(depth=2)
    inf_copy = inferred_tree.clone(depth=2)
    for t in [sp_copy, inf_copy]:
        if t.is_rooted or len(t.seed_node.child_nodes()) == 2:
            t.deroot()
        # Clear root edge length/label — IQ-TREE treats root branch length
        # as evidence of a rooted tree
        t.seed_node.edge.length = None
        if t.seed_node.taxon:
            t.seed_node.taxon = None

    sp_newick = sp_copy.as_string(schema="newick")
    inf_newick = inf_copy.as_string(schema="newick")

    with tempfile.TemporaryDirectory() as tmpdir:
        tree1_path = join(tmpdir, "species.tree")
        tree2_path = join(tmpdir, "inferred.tree")
        prefix = join(tmpdir, "rf_out")

        with open(tree1_path, 'w') as f:
            f.write(sp_newick)
        with open(tree2_path, 'w') as f:
            f.write(inf_newick)

        cmd = f"{IQTREE_PATH} -rf {tree1_path} {tree2_path} -pre {prefix}"
        try:
            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=300
            )
        except subprocess.TimeoutExpired:
            log.warning("    IQ-TREE RF computation timed out")

        rfdist_file = prefix + ".rfdist"
        if isfile(rfdist_file):
            with open(rfdist_file) as f:
                lines = f.readlines()
            for line in lines[1:]:
                parts = line.strip().split()
                if len(parts) >= 2:
                    values = [float(x) for x in parts[1:]]
                    rf_val = max(values)
                    return int(rf_val)

    # Fallback: dendropy symmetric_difference for tiny trees
    log.debug("    IQ-TREE RF failed, falling back to dendropy")
    try:
        # Re-parse into shared TaxonNamespace (required by dendropy)
        shared_tns = dendropy.TaxonNamespace()
        sp_fb = dendropy.Tree.get(data=sp_newick, schema="newick",
                                  taxon_namespace=shared_tns)
        inf_fb = dendropy.Tree.get(data=inf_newick, schema="newick",
                                   taxon_namespace=shared_tns)
        sp_fb.encode_bipartitions()
        inf_fb.encode_bipartitions()
        return treecompare.symmetric_difference(sp_fb, inf_fb)
    except Exception as e:
        log.warning(f"    dendropy RF fallback also failed: {e}")
        return None


def compute_rf_and_ic(species_tree, inferred_tree):
    """
    Compute both standard RF distance (via IQ-TREE) and incompatible-splits (IC).

    RF is computed using IQ-TREE for consistency with the rest of the project.
    IC counts bipartitions in the inferred tree that are incompatible with
    (contradict) bipartitions in the species tree. This handles polytomies
    better than standard RF.

    Returns: dict with RF, nRF, IC, nIC
    """
    # Get leaf labels from both trees
    sp_labels = {leaf.taxon.label for leaf in species_tree.leaf_node_iter()
                 if leaf.taxon}
    inf_labels = {leaf.taxon.label for leaf in inferred_tree.leaf_node_iter()
                  if leaf.taxon}

    common = sp_labels & inf_labels
    n_taxa = len(common)

    if n_taxa < MIN_UNIQUE_TAXA:
        return None

    # Re-parse both trees into a SHARED TaxonNamespace (dendropy requirement for IC)
    shared_tns = dendropy.TaxonNamespace()

    sp_newick = species_tree.as_string(schema="newick")
    inf_newick = inferred_tree.as_string(schema="newick")

    species_tree = dendropy.Tree.get(
        data=sp_newick, schema="newick", taxon_namespace=shared_tns)
    inferred_tree = dendropy.Tree.get(
        data=inf_newick, schema="newick", taxon_namespace=shared_tns)

    # Prune to common set if needed
    if sp_labels != inf_labels:
        sp_extra = sp_labels - common
        inf_extra = inf_labels - common
        if sp_extra:
            sp_taxa_to_remove = [t for t in shared_tns if t.label in sp_extra]
            species_tree.prune_taxa(sp_taxa_to_remove)
        if inf_extra:
            inf_taxa_to_remove = [t for t in shared_tns if t.label in inf_extra]
            inferred_tree.prune_taxa(inf_taxa_to_remove)
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

    # Incompatible splits (IC) via dendropy bipartitions (supplementary metric):
    # Count bipartitions in inferred tree that are INCOMPATIBLE with
    # the species tree (contradict an existing split, not just add resolution).
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

        # Check if this bipartition is incompatible with any species tree split
        is_incompatible = False
        for sp_bp_pair in sp_bipartitions:
            sp_a, sp_b = sp_bp_pair
            # Two bipartitions are incompatible if all 4 intersections are non-empty
            if ((inf_bp & sp_a) and (inf_bp & sp_b) and
                    (inf_complement & sp_a) and (inf_complement & sp_b)):
                is_incompatible = True
                break
        if is_incompatible:
            ic += 1

    nic = ic / (n_taxa - 3) if (n_taxa - 3) > 0 else 0.0

    return {
        'n_taxa': n_taxa,
        'RF': rf,
        'nRF': nrf,
        'IC': ic,
        'nIC': nic,
        'max_RF': max_rf,
    }


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def process_family(rfam_id, model, dir_au, cache_dir, ncbi, dna_cache):
    """
    Process one family-model pair.

    Returns result dict or None on failure.
    dna_cache is a mutable dict used to cache DNA tree results.
    """
    tree_file = join(dir_au, model, rfam_id, f"{rfam_id}.highestLH.trees")
    if not isfile(tree_file):
        return None, "no_tree_file"

    dna_newick, rna_newick = load_tree_pair(tree_file)
    if dna_newick is None:
        return None, "invalid_tree_file"

    # Parse trees
    try:
        dna_tree = parse_newick(dna_newick)
        rna_tree = parse_newick(rna_newick)
    except Exception as e:
        log.debug(f"    Parse error for {rfam_id}/{model}: {e}")
        return None, "parse_error"

    # Resolve taxon IDs (use DNA tree for resolution — same leaves as RNA)
    label_to_taxid, label_format = resolve_leaf_taxids(
        dna_tree, rfam_id, cache_dir, ncbi
    )

    n_leaves = len([l for l in dna_tree.leaf_node_iter() if l.taxon])
    if len(label_to_taxid) < MIN_UNIQUE_TAXA:
        return None, "too_few_taxids_resolved"

    # Prune duplicates from both trees
    dna_tree, dna_dups, dna_unique = prune_duplicate_taxa(dna_tree, label_to_taxid)

    # Re-parse RNA tree for fresh pruning
    rna_tree = parse_newick(rna_newick)
    rna_tree, rna_dups, rna_unique = prune_duplicate_taxa(rna_tree, label_to_taxid)

    # Use intersection of unique taxids from both trees
    common_taxids = dna_unique & rna_unique
    if len(common_taxids) < MIN_UNIQUE_TAXA:
        return None, "too_few_unique_taxa"

    # Relabel trees to taxon IDs
    relabel_tree_to_taxids(dna_tree, label_to_taxid)
    relabel_tree_to_taxids(rna_tree, label_to_taxid)

    # Build species tree (or use cache)
    taxids_key = frozenset(common_taxids)
    species_tree, n_polytomies, merged_map = get_ncbi_species_tree(common_taxids, ncbi)
    if species_tree is None:
        return None, "species_tree_failed"

    # Apply merged taxid translations to inferred trees (skip if target already exists)
    if merged_map:
        for tree in [dna_tree, rna_tree]:
            existing = {l.taxon.label for l in tree.leaf_node_iter() if l.taxon}
            for leaf in tree.leaf_node_iter():
                if leaf.taxon and leaf.taxon.label in [str(k) for k in merged_map]:
                    new_label = str(merged_map[int(leaf.taxon.label)])
                    if new_label not in existing:
                        existing.discard(leaf.taxon.label)
                        leaf.taxon.label = new_label
                        existing.add(new_label)

    # Check if DNA result is cached
    if rfam_id in dna_cache:
        dna_result = dna_cache[rfam_id]
    else:
        # Compute DNA vs species tree (deep copy species tree for safety)
        sp_copy_dna = species_tree.clone(depth=2)
        dna_copy = dna_tree.clone(depth=2)
        dna_result = compute_rf_and_ic(sp_copy_dna, dna_copy)
        if dna_result is not None:
            dna_cache[rfam_id] = dna_result

    if dna_result is None:
        return None, "rf_computation_failed_dna"

    # Compute RNA vs species tree
    sp_copy_rna = species_tree.clone(depth=2)
    rna_copy = rna_tree.clone(depth=2)
    rna_result = compute_rf_and_ic(sp_copy_rna, rna_copy)

    if rna_result is None:
        return None, "rf_computation_failed_rna"

    delta_nic = dna_result['nIC'] - rna_result['nIC']
    delta_nrf = dna_result['nRF'] - rna_result['nRF']

    return {
        'RNA_family': rfam_id,
        'model': model,
        'n_leaves': n_leaves,
        'n_unique_taxa': len(common_taxids),
        'n_dups_removed': dna_dups,
        'label_format': label_format,
        'IC_DNA': dna_result['IC'],
        'IC_RNA': rna_result['IC'],
        'nIC_DNA': dna_result['nIC'],
        'nIC_RNA': rna_result['nIC'],
        'delta_nIC': delta_nic,
        'RF_DNA': dna_result['RF'],
        'RF_RNA': rna_result['RF'],
        'nRF_DNA': dna_result['nRF'],
        'nRF_RNA': rna_result['nRF'],
        'delta_nRF': delta_nrf,
        'n_polytomies': n_polytomies,
    }, None


def process_all(dataset, rna_filter=None):
    """Main orchestrator. Returns DataFrame of all results."""
    cfg = DATASET_CONFIGS[dataset]
    dir_au = join(cfg['working_dir'], 'outputs', cfg['au_folder'])
    output_dir = join(cfg['working_dir'], 'outputs', cfg['output_folder'])
    cache_dir = join(output_dir, 'cache')
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(cache_dir, exist_ok=True)

    models = SEED_MODELS if dataset == 'seed' else FULL_MODELS

    # Load significance table
    sig_table = load_significance_table(dataset)

    # Init NCBITaxa
    log.info("Initialising NCBITaxa (first run downloads ~500MB DB)...")
    ncbi = NCBITaxa()

    results = []
    skipped = []
    dna_cache = {}

    total_pairs = 0
    for model in models:
        if model not in sig_table.columns:
            log.warning(f"Model {model} not in significance table, skipping")
            continue

        sig_families = sorted(sig_table.index[sig_table[model]])

        if rna_filter:
            sig_families = [r for r in sig_families if r in rna_filter]

        log.info(f"\nProcessing model {model}: {len(sig_families)} significant families")

        for i, rfam_id in enumerate(sig_families):
            if (i + 1) % 100 == 0:
                log.info(f"  {model}: {i + 1}/{len(sig_families)} families processed")

            result, skip_reason = process_family(
                rfam_id, model, dir_au, cache_dir, ncbi, dna_cache
            )

            if result is not None:
                results.append(result)
                total_pairs += 1
            else:
                skipped.append({
                    'RNA_family': rfam_id,
                    'model': model,
                    'reason': skip_reason,
                })

    log.info(f"\nDone! Processed {total_pairs} family-model pairs")
    log.info(f"Skipped {len(skipped)} pairs")

    df = pd.DataFrame(results)

    # Save results
    if not df.empty:
        csv_path = join(output_dir, f"species_tree_comparison_{dataset.upper()}.csv")
        df.to_csv(csv_path, index=False)
        log.info(f"Results saved to {csv_path}")

    # Save skipped
    if skipped:
        skip_df = pd.DataFrame(skipped)
        skip_path = join(output_dir, f"skipped_families_{dataset.upper()}.csv")
        skip_df.to_csv(skip_path, index=False)
        log.info(f"Skipped families saved to {skip_path}")

        # Log skip reason counts
        reason_counts = Counter(s['reason'] for s in skipped)
        for reason, count in reason_counts.most_common():
            log.info(f"  Skip reason: {reason} = {count}")

    return df, output_dir


# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

def compute_summary(df):
    """Compute per-model summary statistics based on nRF (primary metric)."""
    if df.empty:
        return pd.DataFrame()

    summaries = []
    for model, group in df.groupby('model'):
        n = len(group)
        n_rna_closer = int((group['delta_nRF'] > 0).sum())
        n_dna_closer = int((group['delta_nRF'] < 0).sum())
        n_ties = int((group['delta_nRF'] == 0).sum())

        summaries.append({
            'model': model,
            'n_families_tested': n,
            'n_rna_closer': n_rna_closer,
            'n_dna_closer': n_dna_closer,
            'n_ties': n_ties,
            'pct_rna_closer': 100 * n_rna_closer / n if n > 0 else 0,
            'median_delta_nRF': group['delta_nRF'].median(),
            'mean_delta_nRF': group['delta_nRF'].mean(),
        })

    return pd.DataFrame(summaries)


# =============================================================================
# RNA TYPE MAP
# =============================================================================

def build_rna_type_map():
    """Build RNA family → RNA type group mapping from Rfam family.txt."""
    if isfile(RFAM_FAMILY_TXT):
        df = pd.read_csv(
            RFAM_FAMILY_TXT, sep="\t", header=None,
            encoding="latin1", on_bad_lines="skip",
            usecols=[0, 18]
        )
        df.columns = ["rfam_acc", "type"]
        df["type"] = df["type"].astype(str).str.strip()
        df["group"] = df["type"].map(TYPE_MAP).fillna("Other")
        return dict(zip(df["rfam_acc"], df["group"]))

    log.warning("family.txt not found — RNA type breakdown will be skipped.")
    return {}


# =============================================================================
# VISUALIZATION
# =============================================================================

def plot_paired_comparison(df, output_dir, dataset):
    """
    Figure 1: Paired dot plot — nRF(DNA) vs nRF(RNA) per family.
    Points below diagonal = RNA tree closer to species tree.
    """
    if df.empty:
        return

    models = sorted(df['model'].unique())
    n_models = len(models)
    ncols = min(4, n_models)
    nrows = (n_models + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4 * nrows),
                             squeeze=False)

    for idx, model in enumerate(models):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        sub = df[df['model'] == model]

        ax.scatter(sub['nRF_DNA'], sub['nRF_RNA'], alpha=0.4, s=12,
                   color='steelblue', edgecolor='none')

        # Diagonal line
        lim = max(sub['nRF_DNA'].max(), sub['nRF_RNA'].max()) * 1.05
        lim = max(lim, 0.1)
        ax.plot([0, lim], [0, lim], 'k--', lw=0.8, alpha=0.5)

        n_below = int((sub['nRF_RNA'] < sub['nRF_DNA']).sum())
        n_total = len(sub)
        ax.set_title(f"{model}  (RNA closer: {n_below}/{n_total})", fontsize=10)
        ax.set_xlabel("nRF (DNA vs species)")
        ax.set_ylabel("nRF (RNA vs species)")
        ax.set_xlim(-0.02, lim)
        ax.set_ylim(-0.02, lim)
        ax.set_aspect('equal')

    # Hide unused subplots
    for idx in range(n_models, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    fig.suptitle(f"Inferred Tree vs Species Tree — {dataset.upper()}", fontsize=13)
    plt.tight_layout()
    path = join(output_dir, f"paired_nRF_{dataset.upper()}")
    plt.savefig(path + ".png", dpi=150)
    plt.savefig(path + ".pdf")
    plt.close()
    log.info(f"Figure saved: {path}.png")


def plot_delta_distribution(df, output_dir, dataset):
    """
    Figure 2: Distribution of delta_nRF per model.
    Positive values = RNA model tree closer to species tree.
    """
    if df.empty:
        return

    models = sorted(df['model'].unique())
    n_models = len(models)
    ncols = min(4, n_models)
    nrows = (n_models + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.5 * nrows),
                             squeeze=False)

    for idx, model in enumerate(models):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        sub = df[df['model'] == model]
        delta = sub['delta_nRF']

        ax.hist(delta, bins=30, color='steelblue', edgecolor='white',
                alpha=0.7, density=True)
        ax.axvline(0, color='red', lw=1.2, ls='--')
        ax.axvline(delta.median(), color='orange', lw=1.2, ls='-',
                   label=f"median={delta.median():.3f}")

        n_pos = int((delta > 0).sum())
        n_neg = int((delta < 0).sum())
        ax.set_title(f"{model}  (RNA+:{n_pos}, DNA+:{n_neg})", fontsize=10)

        ax.set_xlabel("delta nRF (DNA − RNA)\n(+) = RNA closer")
        ax.set_ylabel("Density")
        ax.legend(fontsize=8)

    for idx in range(n_models, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    fig.suptitle(f"Delta nRF Distribution — {dataset.upper()}", fontsize=13)
    plt.tight_layout()
    path = join(output_dir, f"delta_nRF_{dataset.upper()}")
    plt.savefig(path + ".png", dpi=150)
    plt.savefig(path + ".pdf")
    plt.close()
    log.info(f"Figure saved: {path}.png")


def plot_summary_bar_chart(df, output_dir, dataset):
    """
    Figure 3: Per model — count of families where RNA closer / DNA closer / tie.
    """
    if df.empty:
        return

    models = sorted(df['model'].unique())
    rna_closer = []
    dna_closer = []
    ties = []

    for model in models:
        sub = df[df['model'] == model]
        rna_closer.append(int((sub['delta_nRF'] > 0).sum()))
        dna_closer.append(int((sub['delta_nRF'] < 0).sum()))
        ties.append(int((sub['delta_nRF'] == 0).sum()))

    x = np.arange(len(models))
    width = 0.25

    fig, ax = plt.subplots(figsize=(max(10, len(models) * 0.9), 5))
    b1 = ax.bar(x - width, rna_closer, width, label='RNA closer',
                color='#4CAF50', edgecolor='black', linewidth=0.5)
    b2 = ax.bar(x, dna_closer, width, label='DNA closer',
                color='#F44336', edgecolor='black', linewidth=0.5)
    b3 = ax.bar(x + width, ties, width, label='Tie',
                color='#9E9E9E', edgecolor='black', linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(models, rotation=45, ha='right')
    ax.set_ylabel('Number of families')
    ax.set_title(f'Species Tree Proximity: RNA vs DNA Model — {dataset.upper()}')
    ax.legend()

    # Add count labels
    for bars in [b1, b2, b3]:
        for bar in bars:
            h = bar.get_height()
            if h > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, h + 1,
                        str(int(h)), ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    path = join(output_dir, f"summary_bar_{dataset.upper()}")
    plt.savefig(path + ".png", dpi=150)
    plt.savefig(path + ".pdf")
    plt.close()
    log.info(f"Figure saved: {path}.png")


def plot_by_rna_type(df, output_dir, dataset, type_map):
    """
    Figure 4: Summary bar chart broken down by RNA structural type.
    """
    if df.empty or not type_map:
        return

    df = df.copy()
    df['rna_type'] = df['RNA_family'].map(type_map).fillna('Other')

    # Only include types with enough data
    type_counts = df['rna_type'].value_counts()
    active_types = [t for t in RNA_TYPE_ORDER if type_counts.get(t, 0) >= 5]

    if not active_types:
        log.info("Not enough data per RNA type for breakdown figure")
        return

    models = sorted(df['model'].unique())

    fig, axes = plt.subplots(len(active_types), 1,
                             figsize=(max(8, len(models) * 0.8),
                                      2.5 * len(active_types)),
                             squeeze=False)

    for i, rna_type in enumerate(active_types):
        ax = axes[i][0]
        sub = df[df['rna_type'] == rna_type]

        rna_closer = []
        dna_closer = []
        for model in models:
            msub = sub[sub['model'] == model]
            rna_closer.append(int((msub['delta_nRF'] > 0).sum()))
            dna_closer.append(int((msub['delta_nRF'] < 0).sum()))

        x = np.arange(len(models))
        width = 0.35
        ax.bar(x - width / 2, rna_closer, width, label='RNA closer',
               color='#4CAF50', alpha=0.8)
        ax.bar(x + width / 2, dna_closer, width, label='DNA closer',
               color='#F44336', alpha=0.8)

        ax.set_xticks(x)
        ax.set_xticklabels(models, fontsize=8)
        ax.set_ylabel('Count')
        n_total = len(sub[sub['model'] == models[0]]) if models else 0
        ax.set_title(f"{rna_type} (n={len(sub) // max(len(models), 1)} families)",
                     fontsize=10)
        if i == 0:
            ax.legend(fontsize=8)

    fig.suptitle(f"Species Tree Proximity by RNA Type — {dataset.upper()}",
                 fontsize=13)
    plt.tight_layout()
    path = join(output_dir, f"by_rna_type_{dataset.upper()}")
    plt.savefig(path + ".png", dpi=150)
    plt.savefig(path + ".pdf")
    plt.close()
    log.info(f"Figure saved: {path}.png")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Compare inferred trees against NCBI Taxonomy species trees.'
    )
    parser.add_argument('--dataset', choices=['seed', 'full'], required=True)
    parser.add_argument('--rna', type=str, nargs='*', default=None,
                        help='Process only these RNA families (e.g., RF00740)')
    args = parser.parse_args()

    rna_filter = set(args.rna) if args.rna else None

    df, output_dir = process_all(args.dataset, rna_filter)

    if df.empty:
        log.warning("No results to report.")
        return

    # Print summary
    summary = compute_summary(df)
    summary_path = join(output_dir,
                        f"species_tree_summary_{args.dataset.upper()}.csv")
    summary.to_csv(summary_path, index=False)
    log.info(f"\nSummary saved to {summary_path}")

    print(f"\n{'=' * 70}")
    print(f"Species Tree Comparison Summary ({args.dataset.upper()})")
    print(f"{'=' * 70}")
    for _, row in summary.iterrows():
        print(f"\n  {row['model']}:")
        print(f"    Tested: {row['n_families_tested']} families")
        print(f"    RNA closer to species tree: {row['n_rna_closer']} "
              f"({row['pct_rna_closer']:.1f}%)")
        print(f"    DNA closer to species tree: {row['n_dna_closer']}")
        print(f"    Ties: {row['n_ties']}")
        print(f"    Median delta nRF: {row['median_delta_nRF']:.4f}")

    # Generate figures
    print(f"\nGenerating figures...")
    plot_paired_comparison(df, output_dir, args.dataset)
    plot_delta_distribution(df, output_dir, args.dataset)
    plot_summary_bar_chart(df, output_dir, args.dataset)

    type_map = build_rna_type_map()
    if type_map:
        plot_by_rna_type(df, output_dir, args.dataset, type_map)

    print(f"\nAll output in: {output_dir}")


if __name__ == '__main__':
    main()
