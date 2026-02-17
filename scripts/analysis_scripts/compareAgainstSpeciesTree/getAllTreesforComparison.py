#!/usr/bin/env python3
"""
Convert tree labels to species names for Dendroscope tanglegrams.
Processes multiple RNA models (S16, S16A, S16B, S7A, S7B, etc.)

Input structure:
    base_dir/
    ├── S16/RF00740/RF00740.highestLH.trees   (line1=DNA, line2=S16 RNA)
    ├── S16A/RF00740/RF00740.highestLH.trees  (line1=DNA, line2=S16A RNA)
    ├── S6A/RF00740/RF00740.highestLH.trees   (line1=DNA, line2=S6A RNA)
    └── ...

Output:
    output/RF00740/
    ├── RF00740_DNA_tree.newick
    ├── RF00740_S16_tree.newick
    ├── RF00740_S16A_tree.newick
    ├── RF00740_S6A_tree.newick
    ├── ... (all RNA model trees)
    ├── RF00740_species_list.txt      (for TimeTree upload)
    └── RF00740_label_mapping.tsv

Usage:
    python compare_trees_final.py --rfam RF00740 --base_dir /path/to/260208_AU_Test_RAxML_RNAmodels

    # Or process multiple families:
    python compare_trees_final.py --rfam RF00001 RF00002 RF00740 --base_dir /path/to/base
"""

import argparse
import os
import re
import requests

try:
    import dendropy

    HAS_DENDROPY = True
except ImportError:
    HAS_DENDROPY = False
    print("WARNING: dendropy not installed. Install with: pip install dendropy")

# All RNA models
RNA_MODELS = ['S16', 'S16A', 'S16B', 'S7A', 'S7B', 'S7C', 'S7D', 'S7E', 'S7F', 'S6A', 'S6B', 'S6C', 'S6D', 'S6E']


def download_rfam_tree(rfam_id):
    """Download tree from Rfam."""
    url = f"https://rfam.org/family/{rfam_id}/tree"
    print(f"  Downloading Rfam tree for {rfam_id}...")

    try:
        response = requests.get(url, timeout=60)
        if response.status_code == 200:
            return response.text.strip()
    except Exception as e:
        print(f"  Error: {e}")
    return None


def parse_rfam_tree_for_mapping(rfam_tree_text):
    """
    Parse Rfam tree to create mapping: accession → species_name

    Rfam label: _URS000063031F_10116/1-75_Rattus_norvegicus_{Norw..[10116].1
    We want:    URS000063031F_10116/1-75 → Rattus_norvegicus
    """
    acc_to_species = {}

    # Pattern: accession/coords followed by species name (Genus_species)
    pattern = r'_?([A-Z0-9]+_\d+/\d+-\d+)_([A-Z][a-z]+_[a-z]+)'

    for match in re.finditer(pattern, rfam_tree_text):
        accession = match.group(1)
        species = match.group(2)
        acc_to_species[accession] = species

    print(f"  Found {len(acc_to_species)} accession → species mappings")

    return acc_to_species


def clean_rfam_tree_for_parsing(rfam_tree_text):
    """
    Clean Rfam tree so dendropy can parse it.

    Problem: Labels like _URS000075E9DB_9796/1-75_Equus_caballus_{horse}[9796].1
    contain { } which dendropy can't handle.

    Solution: Remove everything from { to the next ] (i.e., remove {horse}[9796].1)
    """
    cleaned = re.sub(r'_{[^:]+', '', rfam_tree_text)
    return cleaned


def convert_tree_labels(tree_newick, acc_to_species, is_rfam_tree=False):
    """
    Convert tree labels from accessions to species names.

    Handles two formats:
    1. Your trees:  URS000063031F_10116/1-75
    2. Rfam trees:  _URS000075E9DB_9796/1-75_Equus_caballus_{horse}[9796].1

    Output label: Rattus_norvegicus
    """
    if not HAS_DENDROPY:
        print("ERROR: dendropy required")
        return None

    try:
        # Clean Rfam tree if needed (remove {common_name}[taxid].1 patterns)
        if is_rfam_tree:
            tree_newick = clean_rfam_tree_for_parsing(tree_newick)

        tree = dendropy.Tree.get(data=tree_newick, schema="newick")

        converted = 0
        failed = []

        for leaf in tree.leaf_node_iter():
            if leaf.taxon is None:
                continue

            label = leaf.taxon.label

            # Dendropy converts underscores to spaces in Newick labels
            # Convert spaces back to underscores for lookup
            label_with_underscores = label.replace(' ', '_')

            # DEBUG: Show first few labels
            if converted == 0 and len(failed) < 3:
                print(f"    DEBUG: label = '{label}'")
                print(f"    DEBUG: label_fixed = '{label_with_underscores}'")
                print(f"    DEBUG: acc_to_species keys sample = {list(acc_to_species.keys())[:3]}")

            # Try direct lookup (with underscores restored)
            if label_with_underscores in acc_to_species:
                leaf.taxon.label = acc_to_species[label_with_underscores]
                converted += 1
                continue

            # Try to extract accession pattern from label (handles both formats)
            # Pattern: ACCESSION_TAXID/start-end
            match = re.search(r'([A-Z0-9]+[_ ]\d+/\d+-\d+)', label_with_underscores)
            if match:
                acc = match.group(1).replace(' ', '_')
                if acc in acc_to_species:
                    leaf.taxon.label = acc_to_species[acc]
                    converted += 1
                    continue

            # For Rfam format: try to extract species name directly from label
            # Format: _ACCESSION_Species_name_{common}[taxid].1
            match = re.search(r'_([A-Z][a-z]+_[a-z]+)', label)
            if match:
                species = match.group(1)
                leaf.taxon.label = species
                converted += 1
                continue

            failed.append(label)

        if failed:
            print(f"    Warning: {len(failed)} labels not converted")

        print(f"    Converted {converted} labels")

        return tree.as_string(schema="newick")

    except Exception as e:
        print(f"  Error: {e}")
        import traceback
        traceback.print_exc()
        return None

    except Exception as e:
        print(f"  Error: {e}")
        import traceback
        traceback.print_exc()
        return None


def process_family(rfam_id, base_dir, output_dir):
    """Process one Rfam family across all RNA models."""

    family_dir = os.path.join(output_dir, rfam_id)
    os.makedirs(family_dir, exist_ok=True)

    print(f"\n{'=' * 60}")
    print(f"Processing {rfam_id}")
    print(f"Output: {family_dir}")
    print(f"{'=' * 60}")

    # Step 1: Download Rfam tree and create mapping
    print("\nStep 1: Download Rfam tree and extract mapping...")
    rfam_tree = download_rfam_tree(rfam_id)
    if not rfam_tree:
        print("ERROR: Failed to download Rfam tree")
        return False

    # Save original Rfam tree
    with open(os.path.join(family_dir, f"{rfam_id}_rfam_tree_original.newick"), 'w') as f:
        f.write(rfam_tree)

    # Parse mapping
    acc_to_species = parse_rfam_tree_for_mapping(rfam_tree)

    # Convert and save Rfam tree with species names
    print("  Converting Rfam tree labels...")
    rfam_converted = convert_tree_labels(rfam_tree, acc_to_species, is_rfam_tree=True)
    if rfam_converted:
        with open(os.path.join(family_dir, f"{rfam_id}_rfam_tree.newick"), 'w') as f:
            f.write(rfam_converted)
        print(f"  Saved: {rfam_id}_rfam_tree.newick")
    else:
        print("  ERROR: Failed to convert Rfam tree labels")

    if not acc_to_species:
        print("ERROR: Could not extract any mappings from Rfam tree")
        return False

    # Save mapping
    with open(os.path.join(family_dir, f"{rfam_id}_label_mapping.tsv"), 'w') as f:
        f.write("accession\tspecies_name\n")
        for acc, species in sorted(acc_to_species.items()):
            f.write(f"{acc}\t{species}\n")

    # Step 2: Export species list for TimeTree
    print("\nStep 2: Export species list for TimeTree...")
    species_names = sorted(set(acc_to_species.values()))

    with open(os.path.join(family_dir, f"{rfam_id}_species_list.txt"), 'w') as f:
        for species in species_names:
            f.write(species.replace('_', ' ') + '\n')

    print(f"  Saved {len(species_names)} species to {rfam_id}_species_list.txt")

    # Step 3: Process each RNA model
    print("\nStep 3: Process trees from each RNA model...")

    dna_tree_saved = False
    models_processed = 0

    for model in RNA_MODELS:
        # Path to the tree file for this model
        tree_file = os.path.join(base_dir, model, rfam_id, f"{rfam_id}.highestLH.trees")

        if not os.path.exists(tree_file):
            print(f"  {model}: File not found, skipping")
            continue

        # Read trees
        with open(tree_file) as f:
            lines = [line.strip() for line in f if line.strip()]

        if len(lines) < 2:
            print(f"  {model}: File has less than 2 lines, skipping")
            continue

        dna_tree_text = lines[0]
        rna_tree_text = lines[1]

        # Save DNA tree (only once, since it's the same for all models)
        if not dna_tree_saved:
            print(f"\n  Converting DNA tree...")
            dna_converted = convert_tree_labels(dna_tree_text, acc_to_species)
            if dna_converted:
                with open(os.path.join(family_dir, f"{rfam_id}_DNA_tree.newick"), 'w') as f:
                    f.write(dna_converted)
                print(f"    Saved: {rfam_id}_DNA_tree.newick")
                dna_tree_saved = True

        # Save RNA tree for this model
        print(f"\n  Converting {model} RNA tree...")
        rna_converted = convert_tree_labels(rna_tree_text, acc_to_species)
        if rna_converted:
            with open(os.path.join(family_dir, f"{rfam_id}_{model}_tree.newick"), 'w') as f:
                f.write(rna_converted)
            print(f"    Saved: {rfam_id}_{model}_tree.newick")
            models_processed += 1

    # Summary
    print(f"\n{'=' * 60}")
    print(f"DONE! Processed {models_processed} RNA models for {rfam_id}")
    print(f"{'=' * 60}")
    print(f"\nOutput files in {family_dir}/:")
    print(f"  - {rfam_id}_DNA_tree.newick")
    for model in RNA_MODELS:
        tree_path = os.path.join(family_dir, f"{rfam_id}_{model}_tree.newick")
        if os.path.exists(tree_path):
            print(f"  - {rfam_id}_{model}_tree.newick")
    print(f"  - {rfam_id}_rfam_tree.newick")
    print(f"  - {rfam_id}_species_list.txt (upload to TimeTree.org)")
    print(f"  - {rfam_id}_label_mapping.tsv")

    print(f"\nNext steps:")
    print(f"  1. Upload {rfam_id}_species_list.txt to https://timetree.org")
    print(f"  2. Download the species tree")
    print(f"  3. Open trees in Dendroscope for tanglegrams!")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Convert tree labels to species names across all RNA models"
    )
    parser.add_argument('--rfam', type=str, nargs='+', required=True,
                        help='Rfam ID(s) (e.g., RF00740 or RF00001 RF00002 RF00740)')
    parser.add_argument('--base_dir', type=str, required=True,
                        help='Base directory containing MODEL/RFAM_ID/files')
    parser.add_argument('--output', type=str, default='output',
                        help='Output directory (default: output)')

    args = parser.parse_args()

    if not HAS_DENDROPY:
        print("ERROR: pip install dendropy")
        return

    # Process each Rfam family
    for rfam_id in args.rfam:
        process_family(rfam_id, args.base_dir, args.output)


if __name__ == "__main__":
    main()