#!/usr/bin/env python3
"""
Robinson-Foulds Distance Calculator
====================================

Compares tree sets between a baseline model (e.g., DNA) and a test model (e.g., DNA_2, S16A).

Workflow:
1. Find common RNA families between baseline and test model
2. Combine 10 trees per family into single files
3. Compute RF distances using IQ-TREE (all trees)
4. Find best tree per model (highest log-likelihood)
5. Compute RF distance between best trees
6. Generate summary tables

Usage:
    python rf_distance_calculator.py --baseline DNA --model DNA_2
    python rf_distance_calculator.py --baseline DNA --model S16A --iqtree /path/to/iqtree2

Output:
    DNA_vs_DNA_2/
    ├── rf_summary.csv              # Summary of all families
    └── RF00001/
        ├── RF00001.DNA             # Combined DNA trees (10)
        ├── RF00001.DNA_2           # Combined DNA_2 trees (10)
        ├── RF00001.rfdist          # RF distances (10x10 matrix)
        ├── RF00001.DNA.best.11.newick      # Best DNA tree
        ├── RF00001.DNA.best.11.info        # Best DNA tree metadata
        ├── RF00001.DNA_2.best.15.newick    # Best DNA_2 tree
        ├── RF00001.DNA_2.best.15.info      # Best DNA_2 tree metadata
        └── RF00001.best.rfdist     # RF distance between best trees
"""

import os
import argparse
import logging
import subprocess
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass
from typing import Optional, List, Set, Tuple
import csv
import re


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class Config:
    """Configuration for RF distance calculation."""

    # Base paths - edit these for your setup
    #DIR_WORKING: str = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels"
    DIR_WORKING = "/Users/u7875558/RNAPhylo//fullAlignment"

    # IQ-TREE executable
    IQTREE_PATH: str = "/Users/u7875558/tools/build-iqtree3/iqtree3"

    # Expected files per family
    EXPECTED_TREES: int = 10
    EXPECTED_FILES: int = 50

    # Derived paths
    DIR_OUTPUTS: str = ""
    DIR_TREES: str = ""
    DIR_RF: str = ""
    DIR_LOGS: str = ""

    def __post_init__(self):
        self.DIR_OUTPUTS = os.path.join(self.DIR_WORKING, "outputs")
        self.DIR_TREES = os.path.join(self.DIR_OUTPUTS, "inferred_trees")
        self.DIR_RF = os.path.join(self.DIR_OUTPUTS, "260220_RF_distances")
        self.DIR_LOGS = os.path.join(self.DIR_WORKING, "logs", "RF_distance")


CONFIG = Config()


# ══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def setup_logging(model: str, config: Config) -> str:
    """Setup logging and return log file path."""
    os.makedirs(config.DIR_LOGS, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file = os.path.join(config.DIR_LOGS, f"{timestamp}_{model}.log")

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return log_file


def get_model_dir(model: str, config: Config) -> str:
    """Get directory path for a model."""
    return os.path.join(config.DIR_TREES, model)


def get_families(model_dir: str) -> Set[str]:
    """Get set of RNA families in a model directory."""
    if not os.path.isdir(model_dir):
        return set()
    return {f for f in os.listdir(model_dir)
            if os.path.isdir(os.path.join(model_dir, f)) and f.startswith("RF")}


def count_trees(family_dir: str) -> int:
    """Count RAxML_bestTree files in a directory."""
    if not os.path.isdir(family_dir):
        return 0
    return len([f for f in os.listdir(family_dir) if f.startswith("RAxML_bestTree")])


def count_files(family_dir: str) -> int:
    """Count all files in a directory."""
    if not os.path.isdir(family_dir):
        return 0
    return len(os.listdir(family_dir))


def parse_raxml_info_likelihood(info_file: str) -> Optional[float]:
    """
    Parse RAxML_info file to extract final log-likelihood.

    Returns log-likelihood value or None if not found.
    """
    if not os.path.exists(info_file):
        return None

    try:
        with open(info_file, "r") as f:
            content = f.read()

        # Look for "Final GAMMA-based Score of best tree" or similar
        # Pattern for RAxML likelihood
        patterns = [
            r"Final GAMMA-based Score of best tree[:\s]+([-\d.]+)",
            r"Final GAMMA  likelihood[:\s]+([-\d.]+)",
            r"final GAMMA-based Likelihood[:\s]+([-\d.]+)",
            r"Likelihood[:\s]+([-\d.]+)",
        ]

        for pattern in patterns:
            match = re.search(pattern, content, re.IGNORECASE)
            if match:
                return float(match.group(1))

        return None
    except Exception as e:
        logging.debug(f"Error parsing {info_file}: {e}")
        return None


def find_best_tree(family_dir: str) -> Optional[Tuple[str, str, float]]:
    """
    Find the best tree (highest log-likelihood) in a family directory.

    Returns tuple of (tree_file_path, seed_number, likelihood) or None.
    """
    if not os.path.isdir(family_dir):
        return None

    best_tree = None
    best_seed = None
    best_likelihood = float('-inf')

    # Find all RAxML_info files
    for f in os.listdir(family_dir):
        if f.startswith("RAxML_info."):
            # Extract seed from filename: RAxML_info.RF00001.11 -> 11
            parts = f.split(".")
            if len(parts) >= 3:
                seed = parts[-1]

                info_path = os.path.join(family_dir, f)
                likelihood = parse_raxml_info_likelihood(info_path)

                if likelihood is not None and likelihood > best_likelihood:
                    best_likelihood = likelihood
                    best_seed = seed

                    # Corresponding tree file
                    rna = parts[1]
                    tree_file = f"RAxML_bestTree.{rna}.{seed}"
                    tree_path = os.path.join(family_dir, tree_file)

                    if os.path.exists(tree_path):
                        best_tree = tree_path

    if best_tree:
        return (best_tree, best_seed, best_likelihood)
    return None


def save_best_tree(tree_path: str, seed: str, likelihood: float,
                   output_dir: str, rna: str, model: str) -> str:
    """
    Save best tree to output directory with metadata.

    Returns path to saved tree file.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Zero-pad seed to 2 digits for consistency
    seed_padded = seed.zfill(2)

    # Save tree: RF00001.DNA.best.01.newick
    output_file = os.path.join(output_dir, f"{rna}.{model}.best.{seed_padded}.newick")

    with open(tree_path, "r") as f:
        tree_content = f.read()

    with open(output_file, "w") as f:
        f.write(tree_content)

    # Save metadata: RF00001.DNA.best.01.info
    meta_file = os.path.join(output_dir, f"{rna}.{model}.best.{seed_padded}.info")
    with open(meta_file, "w") as f:
        f.write(f"model: {model}\n")
        f.write(f"seed: {seed}\n")
        f.write(f"likelihood: {likelihood}\n")
        f.write(f"source: {tree_path}\n")

    return output_file


def compute_best_tree_rf(baseline_tree: str, test_tree: str,
                         output_dir: str, rna: str,
                         iqtree_path: str) -> Optional[str]:
    """
    Compute RF distance between two best trees.

    Returns path to .rfdist file.
    """
    # IQ-TREE needs two separate tree files for -rf comparison
    prefix = os.path.join(output_dir, f"{rna}.best")
    cmd = f"{iqtree_path} -rf {baseline_tree} {test_tree} -pre {prefix}"

    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)

        rfdist_file = prefix + ".rfdist"
        if os.path.exists(rfdist_file):
            return rfdist_file

        # Log error if failed
        if result.returncode != 0:
            logging.debug(f"IQ-TREE best tree RF failed for {rna}: {result.stderr}")

        return None
    except Exception as e:
        logging.error(f"Error computing best tree RF for {rna}: {e}")
        return None


# ══════════════════════════════════════════════════════════════════════════════
# CORE FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def find_common_families(baseline_dir: str, test_dir: str,
                         expected_trees: int = 10) -> Tuple[List[str], dict]:
    """
    Find RNA families present in both models with sufficient trees.

    Returns:
        - List of valid common families
        - Dict with stats about excluded families
    """
    baseline_families = get_families(baseline_dir)
    test_families = get_families(test_dir)

    stats = {
        "baseline_only": [],
        "test_only": [],
        "incomplete_baseline": [],
        "incomplete_test": [],
        "valid": []
    }

    # Families only in one model
    stats["baseline_only"] = sorted(baseline_families - test_families)
    stats["test_only"] = sorted(test_families - baseline_families)

    # Check common families for completeness
    common = baseline_families & test_families

    for rna in sorted(common):
        baseline_count = count_trees(os.path.join(baseline_dir, rna))
        test_count = count_trees(os.path.join(test_dir, rna))

        if baseline_count < expected_trees:
            stats["incomplete_baseline"].append((rna, baseline_count))
        elif test_count < expected_trees:
            stats["incomplete_test"].append((rna, test_count))
        else:
            stats["valid"].append(rna)

    return stats["valid"], stats


def combine_trees(family_dir: str, output_file: str) -> bool:
    """
    Combine all RAxML_bestTree files into a single file.

    Returns True if successful.
    """
    if not os.path.isdir(family_dir):
        return False

    tree_files = sorted([f for f in os.listdir(family_dir)
                         if f.startswith("RAxML_bestTree")])

    if not tree_files:
        return False

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, "w") as out:
        for tf in tree_files:
            tree_path = os.path.join(family_dir, tf)
            with open(tree_path, "r") as f:
                content = f.read().strip()
                if content:
                    out.write(content + "\n")

    return True


def compute_rf_distance(tree_file1: str, tree_file2: str,
                        output_dir: str, rna: str,
                        iqtree_path: str = "iqtree2") -> Optional[str]:
    """
    Compute RF distance between two tree sets using IQ-TREE.

    Returns path to .rfdist file if successful, None otherwise.
    """
    # IQ-TREE command: iqtree2 -rf tree1 tree2 -pre prefix
    prefix = os.path.join(output_dir, rna)

    cmd = f"{iqtree_path} -rf {tree_file1} {tree_file2} -pre {prefix}"

    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )

        if result.returncode != 0:
            logging.error(f"IQ-TREE failed: {result.stderr}")
            return None

        # Find the .rfdist file
        rfdist_file = prefix + ".rfdist"
        if os.path.exists(rfdist_file):
            return rfdist_file

        # Sometimes IQ-TREE names it differently
        for f in os.listdir(output_dir):
            if f.endswith(".rfdist"):
                return os.path.join(output_dir, f)

        return None

    except subprocess.TimeoutExpired:
        logging.error(f"IQ-TREE timed out for {tree_file1}")
        return None
    except Exception as e:
        logging.error(f"Error running IQ-TREE: {e}")
        return None


def parse_rfdist_matrix(rfdist_file: str) -> Optional[List[List[float]]]:
    """Parse IQ-TREE RF distance matrix file."""
    if not os.path.exists(rfdist_file):
        return None

    try:
        with open(rfdist_file, "r") as f:
            lines = f.readlines()

        # Skip header line (number of trees)
        matrix = []
        for line in lines[1:]:
            parts = line.strip().split()
            if len(parts) > 1:
                # First element is tree name, rest are distances
                values = [float(x) for x in parts[1:]]
                matrix.append(values)

        return matrix
    except Exception as e:
        logging.error(f"Error parsing {rfdist_file}: {e}")
        return None


def process_family(rna: str, baseline_dir: str, test_dir: str,
                   output_dir: str, baseline_name: str, test_name: str,
                   config: Config) -> dict:
    """
    Process a single RNA family: combine trees and compute RF distance.

    Returns dict with results.
    """
    result = {
        "rna": rna,
        "status": "unknown",
        "best_status": "unknown",  # Track best tree computation separately
        "baseline_trees": 0,
        "test_trees": 0,
        "rfdist_file": None,
        "best_rfdist_file": None,
        "baseline_best_seed": None,
        "test_best_seed": None,
        "baseline_best_likelihood": None,
        "test_best_likelihood": None,
        "error": None
    }

    family_output_dir = os.path.join(output_dir, rna)
    os.makedirs(family_output_dir, exist_ok=True)

    # Paths - use actual model names
    baseline_family_dir = os.path.join(baseline_dir, rna)
    test_family_dir = os.path.join(test_dir, rna)
    baseline_combined = os.path.join(family_output_dir, f"{rna}.{baseline_name}")
    test_combined = os.path.join(family_output_dir, f"{rna}.{test_name}")

    # Combine baseline trees
    if not combine_trees(baseline_family_dir, baseline_combined):
        result["status"] = "failed"
        result["error"] = "Could not combine baseline trees"
        return result

    result["baseline_trees"] = count_trees(baseline_family_dir)

    # Combine test trees
    if not combine_trees(test_family_dir, test_combined):
        result["status"] = "failed"
        result["error"] = "Could not combine test trees"
        return result

    result["test_trees"] = count_trees(test_family_dir)

    # Check if main RF distance already computed
    main_rfdist = os.path.join(family_output_dir, f"{rna}.rfdist")
    if os.path.exists(main_rfdist):
        result["status"] = "exists"
        result["rfdist_file"] = main_rfdist
    else:
        # Compute RF distance
        rfdist_file = compute_rf_distance(
            baseline_combined,
            test_combined,
            family_output_dir,
            rna,
            config.IQTREE_PATH
        )

        if rfdist_file:
            result["status"] = "success"
            result["rfdist_file"] = rfdist_file
        else:
            result["status"] = "failed"
            result["error"] = "IQ-TREE RF computation failed"

    # ── BEST TREE PROCESSING ──
    # Check if best.rfdist already exists
    best_rfdist_path = os.path.join(family_output_dir, f"{rna}.best.rfdist")
    if os.path.exists(best_rfdist_path):
        result["best_status"] = "exists"
        result["best_rfdist_file"] = best_rfdist_path
        # Still need to get seed info from existing .info files
        for f in os.listdir(family_output_dir):
            if f.endswith(".info") and baseline_name in f:
                info_path = os.path.join(family_output_dir, f)
                with open(info_path, "r") as fh:
                    for line in fh:
                        if line.startswith("seed:"):
                            result["baseline_best_seed"] = line.split(":")[1].strip()
                        if line.startswith("likelihood:"):
                            result["baseline_best_likelihood"] = float(line.split(":")[1].strip())
            elif f.endswith(".info") and test_name in f:
                info_path = os.path.join(family_output_dir, f)
                with open(info_path, "r") as fh:
                    for line in fh:
                        if line.startswith("seed:"):
                            result["test_best_seed"] = line.split(":")[1].strip()
                        if line.startswith("likelihood:"):
                            result["test_best_likelihood"] = float(line.split(":")[1].strip())
        return result

    # Find best tree for baseline model
    baseline_best = find_best_tree(baseline_family_dir)
    if baseline_best:
        baseline_tree_path, baseline_seed, baseline_likelihood = baseline_best
        result["baseline_best_seed"] = baseline_seed
        result["baseline_best_likelihood"] = baseline_likelihood

        # Save best tree
        baseline_best_saved = save_best_tree(
            baseline_tree_path, baseline_seed, baseline_likelihood,
            family_output_dir, rna, baseline_name
        )
    else:
        baseline_best_saved = None

    # Find best tree for test model
    test_best = find_best_tree(test_family_dir)
    if test_best:
        test_tree_path, test_seed, test_likelihood = test_best
        result["test_best_seed"] = test_seed
        result["test_best_likelihood"] = test_likelihood

        # Save best tree
        test_best_saved = save_best_tree(
            test_tree_path, test_seed, test_likelihood,
            family_output_dir, rna, test_name
        )
    else:
        test_best_saved = None

    # Compute RF distance between best trees
    if baseline_best_saved and test_best_saved:
        best_rfdist = compute_best_tree_rf(
            baseline_best_saved, test_best_saved,
            family_output_dir, rna, config.IQTREE_PATH
        )
        result["best_rfdist_file"] = best_rfdist
        if best_rfdist:
            result["best_status"] = "success"
        else:
            result["best_status"] = "failed"
    else:
        result["best_status"] = "failed"

    return result


def generate_summary(results: List[dict], output_file: str):
    """Generate summary CSV from results."""
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "rna", "status", "best_status", "baseline_trees", "test_trees",
            "rfdist_file", "best_rfdist_file",
            "baseline_best_seed", "baseline_best_likelihood",
            "test_best_seed", "test_best_likelihood",
            "error"
        ])
        writer.writeheader()
        writer.writerows(results)

    logging.info(f"Summary written to: {output_file}")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Compute Robinson-Foulds distances between tree sets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Compare DNA (baseline) vs DNA_2
    python rf_distance_calculator.py --baseline DNA --model DNA_2

    # Compare DNA vs S16A with custom IQ-TREE path
    python rf_distance_calculator.py --baseline DNA --model S16A --iqtree /path/to/iqtree2

    # Only show statistics without computing
    python rf_distance_calculator.py --baseline DNA --model DNA_2 --dry-run
        """
    )

    parser.add_argument("--baseline", "-b", required=True,
                        help="Baseline model (e.g., DNA)")
    parser.add_argument("--model", "-m", required=True,
                        help="Test model to compare (e.g., DNA_2, S16A)")
    parser.add_argument("--iqtree", default=None,
                        help="Path to IQ-TREE executable (default: from config)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show statistics without computing RF distances")
    parser.add_argument("--limit", type=int, default=None,
                        help="Limit number of families to process (for testing)")
    parser.add_argument("--force", action="store_true",
                        help="Recompute even if .rfdist files exist")

    args = parser.parse_args()

    # Setup
    config = CONFIG
    if args.iqtree:
        config.IQTREE_PATH = args.iqtree

    log_file = setup_logging(args.model, config)

    logging.info(f"{'=' * 60}")
    logging.info(f"RF Distance Calculator")
    logging.info(f"Baseline: {args.baseline}")
    logging.info(f"Test model: {args.model}")
    logging.info(f"IQ-TREE: {config.IQTREE_PATH}")
    logging.info(f"{'=' * 60}")

    # Get model directories
    baseline_dir = get_model_dir(args.baseline, config)
    test_dir = get_model_dir(args.model, config)

    if not os.path.isdir(baseline_dir):
        logging.error(f"Baseline directory not found: {baseline_dir}")
        return

    if not os.path.isdir(test_dir):
        logging.error(f"Test model directory not found: {test_dir}")
        return

    # Find common families
    valid_families, stats = find_common_families(
        baseline_dir, test_dir, config.EXPECTED_TREES
    )

    # Report statistics
    logging.info(f"\nFamily Statistics:")
    logging.info(f"  Valid (in both, complete): {len(valid_families)}")
    logging.info(f"  Only in {args.baseline}: {len(stats['baseline_only'])}")
    logging.info(f"  Only in {args.model}: {len(stats['test_only'])}")
    logging.info(f"  Incomplete in {args.baseline}: {len(stats['incomplete_baseline'])}")
    logging.info(f"  Incomplete in {args.model}: {len(stats['incomplete_test'])}")

    if stats["baseline_only"]:
        logging.info(
            f"\nFamilies only in {args.baseline}: {stats['baseline_only'][:10]}{'...' if len(stats['baseline_only']) > 10 else ''}")

    if stats["test_only"]:
        logging.info(
            f"\nFamilies only in {args.model}: {stats['test_only'][:10]}{'...' if len(stats['test_only']) > 10 else ''}")

    if args.dry_run:
        logging.info("\nDry run - no RF distances computed.")
        return

    # Create output directory
    output_dir = os.path.join(config.DIR_RF, f"{args.baseline}_vs_{args.model}")
    os.makedirs(output_dir, exist_ok=True)

    # Apply limit if specified
    families_to_process = valid_families
    if args.limit:
        families_to_process = valid_families[:args.limit]
        logging.info(f"\nLimiting to {args.limit} families for testing")

    # Process families
    logging.info(f"\nProcessing {len(families_to_process)} families...")

    results = []
    success_count = 0
    skip_count = 0
    fail_count = 0
    best_success_count = 0
    best_skip_count = 0
    best_fail_count = 0

    for i, rna in enumerate(families_to_process, 1):
        result = process_family(rna, baseline_dir, test_dir, output_dir,
                                args.baseline, args.model, config)
        results.append(result)

        # Track main RF distance status
        if result["status"] == "success":
            success_count += 1
        elif result["status"] == "exists":
            skip_count += 1
        else:
            fail_count += 1

        # Track best tree RF distance status
        if result.get("best_status") == "success":
            best_success_count += 1
        elif result.get("best_status") == "exists":
            best_skip_count += 1
        else:
            best_fail_count += 1

        if i % 50 == 0 or i == len(families_to_process):
            logging.info(f"  Progress: {i}/{len(families_to_process)} "
                         f"(RF: {success_count} new, {skip_count} exist | "
                         f"Best: {best_success_count} new, {best_skip_count} exist)")

    # Generate summary
    summary_file = os.path.join(output_dir, f"rf_summary.csv")
    generate_summary(results, summary_file)

    # Final report
    logging.info(f"\n{'=' * 60}")
    logging.info("COMPLETE")
    logging.info(f"{'=' * 60}")
    logging.info(f"  Processed: {len(families_to_process)}")
    logging.info(f"\n  Main RF distances (10x10):")
    logging.info(f"    New: {success_count}")
    logging.info(f"    Already existed: {skip_count}")
    logging.info(f"    Failed: {fail_count}")
    logging.info(f"\n  Best tree RF distances:")
    logging.info(f"    New: {best_success_count}")
    logging.info(f"    Already existed: {best_skip_count}")
    logging.info(f"    Failed: {best_fail_count}")
    logging.info(f"\nOutput directory: {output_dir}")
    logging.info(f"Summary: {summary_file}")
    logging.info(f"Log: {log_file}")


if __name__ == "__main__":
    main()