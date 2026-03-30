#!/usr/bin/env python3
"""
Unified Rfam Phylogenetic Pipeline
===================================

Handles both Seed Alignment and Full Alignment datasets.
Optimized for NCI/Gadi with job arrays and batching to minimize SU usage.

Features:
- Supports both seed and full alignment datasets
- Job arrays instead of individual qsub submissions
- Batches multiple families per job
- Command-line interface for all settings
- Handles subsampling for large families (full alignment)
- Complete seed processing: Rfam.seed → sections → FASTA + SS

SEED ALIGNMENT WORKFLOW:
    # Step 1: Split Rfam.seed into individual section files
    python rfam_unified_pipeline.py seed --split-seed

    # Step 2: Parse sections to create FASTA and SS files
    python rfam_unified_pipeline.py seed --parse

    # Step 3: Setup and submit PBS jobs
    python rfam_unified_pipeline.py seed --setup --model S16A
    qsub /scratch/.../pbs/S16A_submit.pbs

FULL ALIGNMENT WORKFLOW (files pre-prepared):
    python rfam_unified_pipeline.py full --setup --model DNA_2
    qsub /scratch/.../pbs/DNA_2_submit.pbs
"""

import os
import argparse
import subprocess
import math
from dataclasses import dataclass
from typing import Optional, List
from abc import ABC, abstractmethod


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class GadiConfig:
    """Base Gadi configuration - edit these for your setup."""

    # ── YOUR GADI PROJECT ──
    PROJECT: str = "dx61"

    # ── BASE PATH ──
    DIR_SCRATCH: str = "/scratch/dx61/vh5686/tmp/RNAPhylo"

    # ── RAXML ──
    RAXML_EXECUTABLE: str = "/scratch/dx61/vh5686/tmp/RNAPhylo/tools/standard-RAxML/raxmlHPC"
    RAXML_DNA_SCRIPT: str = "/scratch/dx61/vh5686/tmp/RNAPhylo/scripts/bashFiles/raxml.sh"
    RAXML_RNA_SCRIPT: str = "/scratch/dx61/vh5686/tmp/RNAPhylo/scripts/bashFiles/raxmlP.sh"

    # ── JOB SETTINGS ──
    FAMILIES_PER_JOB: int = 10
    NCPUS_PER_JOB: int = 10  # RAxML default is single-threaded
    MEM_PER_JOB: str = "40GB"
    WALLTIME_PER_JOB: str = "48:00:00"
    QUEUE: str = "normal"

    # ── PROCESSING ──
    MIN_SEQUENCES: int = 4
    EXPECTED_OUTPUT_FILES: int = 50


@dataclass
class SeedAlignmentConfig(GadiConfig):
    """Configuration for Seed Alignment dataset."""

    DATASET: str = "seed"

    # Derived paths (set in __post_init__)
    DIR_WORKING: str = ""
    DIR_INPUTS: str = ""
    DIR_OUTPUTS: str = ""
    DIR_TREES: str = ""
    DIR_SECTION: str = ""
    DIR_FASTA: str = ""
    DIR_SS: str = ""
    DIR_LOGS: str = ""
    DIR_PBS: str = ""
    RFAM_SEED: str = ""

    def __post_init__(self):
        self.DIR_WORKING = os.path.join(self.DIR_SCRATCH, "seedAlignment")
        self.DIR_INPUTS = os.path.join(self.DIR_WORKING, "inputs")
        self.DIR_OUTPUTS = os.path.join(self.DIR_WORKING, "outputs")
        self.DIR_TREES = os.path.join(self.DIR_OUTPUTS, "inferred_trees")
        self.DIR_SECTION = os.path.join(self.DIR_INPUTS, "sections")
        self.DIR_FASTA = os.path.join(self.DIR_INPUTS, "fasta_files")
        self.DIR_SS = os.path.join(self.DIR_INPUTS, "ss_files")
        self.DIR_LOGS = os.path.join(self.DIR_WORKING, "logs")
        self.DIR_PBS = os.path.join(self.DIR_WORKING, "pbs")
        self.RFAM_SEED = os.path.join(self.DIR_INPUTS, "Rfam.seed")


@dataclass
class FullAlignmentConfig(GadiConfig):
    """Configuration for Full Alignment dataset."""

    DATASET: str = "full"

    # Full alignment specific
    DIR_SUBSAMP: str = ""  # Subsampled sequences for large families

    # Derived paths
    DIR_WORKING: str = ""
    DIR_INPUTS: str = ""
    DIR_OUTPUTS: str = ""
    DIR_TREES: str = ""
    DIR_FASTA: str = ""
    DIR_SS: str = ""
    DIR_LOGS: str = ""
    DIR_PBS: str = ""

    def __post_init__(self):
        self.DIR_WORKING = os.path.join(self.DIR_SCRATCH, "fullAlignment")
        self.DIR_INPUTS = os.path.join(self.DIR_WORKING, "inputs")
        self.DIR_OUTPUTS = os.path.join(self.DIR_WORKING, "outputs")
        self.DIR_TREES = os.path.join(self.DIR_OUTPUTS, "inferred_trees")
        self.DIR_FASTA = os.path.join(self.DIR_INPUTS, "fasta")
        self.DIR_SS = os.path.join(self.DIR_INPUTS, "ss_all")
        self.DIR_SUBSAMP = os.path.join(self.DIR_INPUTS, "subsample")
        self.DIR_LOGS = os.path.join(self.DIR_WORKING, "logs")
        self.DIR_PBS = os.path.join(self.DIR_WORKING, "pbs")


def get_config(dataset: str) -> GadiConfig:
    """Factory function to get appropriate config."""
    if dataset == "seed":
        return SeedAlignmentConfig()
    elif dataset == "full":
        return FullAlignmentConfig()
    else:
        raise ValueError(f"Unknown dataset: {dataset}. Use 'seed' or 'full'.")


# ══════════════════════════════════════════════════════════════════════════════
# CONSTANTS (for seed alignment parsing)
# ══════════════════════════════════════════════════════════════════════════════

MATCHING_BRACKETS = {')': '(', ']': '[', '}': '{', '>': '<'}
MATCHING_PSEUDOKNOTS = {'a': 'A', 'b': 'B', 'c': 'C', 'd': 'D'}
OPENING_BRACKETS = set(MATCHING_BRACKETS.values()) | set(MATCHING_PSEUDOKNOTS.values())
UNPAIRED_CHARS = {':': '.', ',': '.', '_': '.', '~': '.', '-': '.'}
PSEUDOKNOT_TO_DOT = {
    'A': '.', 'a': '.', 'B': '.', 'b': '.',
    'C': '.', 'c': '.', 'D': '.', 'd': '.'
}


# ══════════════════════════════════════════════════════════════════════════════
# DATASET HANDLERS
# ══════════════════════════════════════════════════════════════════════════════

class DatasetHandler(ABC):
    """Abstract base class for dataset-specific operations."""

    def __init__(self, config: GadiConfig):
        self.config = config

    @abstractmethod
    def get_family_list(self) -> List[str]:
        """Return list of all family accessions."""
        pass

    @abstractmethod
    def get_fasta_path(self, accession: str) -> Optional[str]:
        """Return path to FASTA file for a family."""
        pass

    @abstractmethod
    def get_ss_path(self, accession: str) -> Optional[str]:
        """Return path to secondary structure file for a family."""
        pass

    def check_sequence_count(self, fasta_path: str) -> bool:
        """Check if FASTA has enough sequences."""
        if not os.path.exists(fasta_path):
            return False

        count = 0
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
                    if count >= self.config.MIN_SEQUENCES:
                        return True
        return False

    def create_directories(self):
        """Create required directories."""
        dirs = [self.config.DIR_OUTPUTS, self.config.DIR_TREES,
                self.config.DIR_LOGS, self.config.DIR_PBS]
        for d in dirs:
            os.makedirs(d, exist_ok=True)


class SeedAlignmentHandler(DatasetHandler):
    """Handler for Seed Alignment dataset.

    Can parse from section files if FASTA/SS don't exist yet.
    """

    def get_family_list(self) -> List[str]:
        """Get families from section files."""
        if not os.path.isdir(self.config.DIR_SECTION):
            return []
        return sorted([
            f.replace('.section', '')
            for f in os.listdir(self.config.DIR_SECTION)
            if f.endswith('.section')
        ])

    def get_fasta_path(self, accession: str) -> Optional[str]:
        """Get FASTA file path."""
        path = os.path.join(self.config.DIR_FASTA, f"{accession}.nodup.fa")
        return path if os.path.exists(path) else None

    def get_ss_path(self, accession: str) -> Optional[str]:
        """Get secondary structure file path."""
        path = os.path.join(self.config.DIR_SS, f"{accession}.dots.ss")
        return path if os.path.exists(path) else None

    def parse_section_if_needed(self, accession: str) -> tuple:
        """Parse section file to create FASTA and SS if they don't exist."""
        fasta_path = self.get_fasta_path(accession)
        ss_path = self.get_ss_path(accession)

        # Already exist
        if fasta_path and ss_path:
            return fasta_path, ss_path

        # Parse from section
        section_file = os.path.join(self.config.DIR_SECTION, f"{accession}.section")
        if not os.path.exists(section_file):
            return None, None

        return self._parse_section(section_file)

    def _parse_section(self, fp_section: str) -> tuple:
        """Parse section file to extract FASTA and SS."""
        with open(fp_section, "r") as f:
            lines = f.readlines()

        # Extract accession
        accession = None
        for line in lines:
            if line.startswith('#=GF AC'):
                parts = line.split()
                if len(parts) > 2:
                    accession = parts[2].strip()
                    break

        if not accession:
            return None, None

        # Extract alignments
        alignments = {}
        for line in lines:
            if line.startswith('#') or line.startswith('/') or not line.strip():
                continue
            parts = line.split()
            if len(parts) == 2:
                name, seq = parts
                if seq not in alignments.values():
                    alignments[name] = seq

        if len(alignments) < self.config.MIN_SEQUENCES:
            return None, None

        # Write FASTA
        fasta_path = os.path.join(self.config.DIR_FASTA, f"{accession}.nodup.fa")
        os.makedirs(os.path.dirname(fasta_path), exist_ok=True)
        with open(fasta_path, "w") as f:
            for name, seq in alignments.items():
                f.write(f">{name}\n{seq}\n")

        # Extract and convert SS
        ss_cons = None
        for line in lines:
            if line.startswith("#=GC SS_cons"):
                ss_cons = line.split(maxsplit=2)[2].strip()
                break

        if not ss_cons:
            return fasta_path, None

        # Convert pseudoknots to dots
        conversion = {**UNPAIRED_CHARS, **PSEUDOKNOT_TO_DOT}
        ss_dots = ''.join(conversion.get(c, c) for c in ss_cons)

        if len(set(ss_dots)) == 1:  # All unpaired
            return fasta_path, None

        ss_path = os.path.join(self.config.DIR_SS, f'{accession}.dots.ss')
        os.makedirs(os.path.dirname(ss_path), exist_ok=True)
        with open(ss_path, 'w') as f:
            f.write(ss_dots + '\n')

        return fasta_path, ss_path


class FullAlignmentHandler(DatasetHandler):
    """Handler for Full Alignment dataset.

    Assumes all files are pre-prepared:
    - FASTA: Either original (.nodup.fa) or subsampled (.subsamp.fa)
    - SS: Secondary structure files (.ss)

    No parsing or conversion needed.
    """

    def get_family_list(self) -> List[str]:
        """Get families from FASTA directory."""
        if not os.path.isdir(self.config.DIR_FASTA):
            return []

        families = set()
        for f in os.listdir(self.config.DIR_FASTA):
            if f.endswith('.fa') or f.endswith('.fasta'):
                # Extract RF accession (handle both .nodup.fa and .fa)
                name = f.replace('.nodup.fa', '').replace('.fa', '').replace('.fasta', '')
                families.add(name)

        # Also check subsample directory
        if os.path.isdir(self.config.DIR_SUBSAMP):
            for f in os.listdir(self.config.DIR_SUBSAMP):
                if f.endswith('.subsamp.fa'):
                    name = f.replace('.subsamp.fa', '')
                    families.add(name)

        return sorted(families)

    def get_fasta_path(self, accession: str) -> Optional[str]:
        """Get FASTA path, preferring subsampled version for large families."""
        # Check subsampled first (these are large families that were downsampled)
        subsamp_path = os.path.join(self.config.DIR_SUBSAMP, f"{accession}.subsamp.fa")
        if os.path.exists(subsamp_path):
            return subsamp_path

        # Fall back to original
        original_path = os.path.join(self.config.DIR_FASTA, f"{accession}.nodup.fa")
        if os.path.exists(original_path):
            return original_path

        # Try without .nodup
        alt_path = os.path.join(self.config.DIR_FASTA, f"{accession}.fa")
        if os.path.exists(alt_path):
            return alt_path

        return None

    def get_ss_path(self, accession: str) -> Optional[str]:
        """Get secondary structure file path."""
        path = os.path.join(self.config.DIR_SS, f"{accession}.ss")
        return path if os.path.exists(path) else None


def get_handler(config: GadiConfig) -> DatasetHandler:
    """Factory function to get appropriate handler."""
    if isinstance(config, SeedAlignmentConfig):
        return SeedAlignmentHandler(config)
    elif isinstance(config, FullAlignmentConfig):
        return FullAlignmentHandler(config)
    else:
        raise ValueError(f"Unknown config type: {type(config)}")


# ══════════════════════════════════════════════════════════════════════════════
# JOB MANAGEMENT
# ══════════════════════════════════════════════════════════════════════════════

def get_pending_families(handler: DatasetHandler, model: str) -> List[str]:
    """Get list of families that need processing."""
    all_families = handler.get_family_list()
    config = handler.config

    pending = []
    model_dir = os.path.join(config.DIR_TREES, model)

    for accession in all_families:
        # Check if already complete
        output_dir = os.path.join(model_dir, accession)
        if os.path.isdir(output_dir):
            num_files = len(os.listdir(output_dir))
            if num_files >= config.EXPECTED_OUTPUT_FILES:
                continue

        # Check if has required input files
        fasta = handler.get_fasta_path(accession)
        ss = handler.get_ss_path(accession)

        if fasta and ss and handler.check_sequence_count(fasta):
            pending.append(accession)

    return pending


def create_job_list(handler: DatasetHandler, model: str) -> Optional[str]:
    """Create file listing families to process, grouped into batches."""
    config = handler.config
    pending = get_pending_families(handler, model)

    if not pending:
        print("All families already processed or no valid families found!")
        return None

    num_batches = math.ceil(len(pending) / config.FAMILIES_PER_JOB)
    job_list_path = os.path.join(config.DIR_PBS, f"{model}_joblist.txt")

    os.makedirs(config.DIR_PBS, exist_ok=True)

    with open(job_list_path, 'w') as f:
        for i in range(num_batches):
            start = i * config.FAMILIES_PER_JOB
            end = min((i + 1) * config.FAMILIES_PER_JOB, len(pending))
            batch = pending[start:end]
            f.write(','.join(batch) + '\n')

    print(f"Created job list: {job_list_path}")
    print(f"  Dataset: {config.DATASET}")
    print(f"  Total families to process: {len(pending)}")
    print(f"  Families per job: {config.FAMILIES_PER_JOB}")
    print(f"  Total jobs (array size): {num_batches}")

    return job_list_path


def create_pbs_script(handler: DatasetHandler, model: str, job_list_path: str) -> str:
    """Generate PBS submission script."""
    config = handler.config

    with open(job_list_path) as f:
        num_batches = sum(1 for _ in f)

    # Determine paths based on dataset type
    if isinstance(config, SeedAlignmentConfig):
        fasta_pattern = "${DIR_FASTA}/${ACCESSION}.nodup.fa"
        ss_pattern = "${DIR_SS}/${ACCESSION}.dots.ss"
        subsamp_logic = ""
    else:
        # Full alignment - check for subsampled version
        fasta_pattern = "$FASTA"
        ss_pattern = "${DIR_SS}/${ACCESSION}.ss"
        subsamp_logic = """
    # Check for subsampled version first (full alignment)
    SUBSAMP="${DIR_SUBSAMP}/${ACCESSION}.subsamp.fa"
    ORIGINAL="${DIR_FASTA}/${ACCESSION}.nodup.fa"
    if [ -f "$SUBSAMP" ]; then
        FASTA="$SUBSAMP"
    else
        FASTA="$ORIGINAL"
    fi
"""

    pbs_script = f"""#!/bin/bash
#PBS -P {config.PROJECT}
#PBS -q {config.QUEUE}
#PBS -l ncpus={config.NCPUS_PER_JOB}
#PBS -l mem={config.MEM_PER_JOB}
#PBS -l walltime={config.WALLTIME_PER_JOB}
#PBS -l storage=scratch/{config.PROJECT}
#PBS -l wd
#PBS -J 1-{num_batches}
#PBS -r y
#PBS -o {config.DIR_LOGS}/{model}_${{PBS_ARRAY_INDEX}}.o
#PBS -e {config.DIR_LOGS}/{model}_${{PBS_ARRAY_INDEX}}.e
#PBS -N rfam_{config.DATASET}_{model}

# ══════════════════════════════════════════════════════════════════════════════
# Rfam Phylogenetic Pipeline - PBS Array Job
# Dataset: {config.DATASET}
# Model: {model}
# ══════════════════════════════════════════════════════════════════════════════

# ── GET THIS JOB'S BATCH ──
JOBLIST="{job_list_path}"
FAMILIES=$(sed -n "${{PBS_ARRAY_INDEX}}p" $JOBLIST)

# ── PATHS ──
MODEL="{model}"
DIR_FASTA="{config.DIR_FASTA}"
DIR_SS="{config.DIR_SS}"
DIR_TREES="{config.DIR_TREES}"
DIR_SUBSAMP="{getattr(config, 'DIR_SUBSAMP', '')}"
RAXML_EXEC="{config.RAXML_EXECUTABLE}"
RAXML_DNA_SCRIPT="{config.RAXML_DNA_SCRIPT}"
RAXML_RNA_SCRIPT="{config.RAXML_RNA_SCRIPT}"
EXPECTED_FILES={config.EXPECTED_OUTPUT_FILES}

# ── LOGGING ──
echo "========================================"
echo "Job: $PBS_JOBID (Array index: $PBS_ARRAY_INDEX)"
echo "Model: $MODEL"
echo "Families: $FAMILIES"
echo "Started: $(date)"
echo "========================================"

# ── PROCESS EACH FAMILY ──
IFS=',' read -ra FAMILY_ARRAY <<< "$FAMILIES"

for ACCESSION in "${{FAMILY_ARRAY[@]}}"; do
    echo ""
    echo "Processing: $ACCESSION"

    OUTPUT="${{DIR_TREES}}/${{MODEL}}/${{ACCESSION}}"

    # Skip if already complete
    if [ -d "$OUTPUT" ]; then
        FILE_COUNT=$(ls -1 "$OUTPUT" 2>/dev/null | wc -l)
        if [ "$FILE_COUNT" -ge "$EXPECTED_FILES" ]; then
            echo "  Already complete ($FILE_COUNT files), skipping"
            continue
        else
            echo "  Incomplete ($FILE_COUNT files), restarting"
            rm -rf "$OUTPUT"
        fi
    fi
    {subsamp_logic}
    # Set file paths
    FASTA="{fasta_pattern}"
    SS="{ss_pattern}"

    # Verify files exist
    if [ ! -f "$FASTA" ]; then
        echo "  ERROR: FASTA not found: $FASTA"
        continue
    fi
    if [ ! -f "$SS" ]; then
        echo "  ERROR: SS not found: $SS"
        continue
    fi

    # Create output directory
    mkdir -p "$OUTPUT"

    # Run RAxML
    echo "  FASTA: $FASTA"
    echo "  SS: $SS"
    echo "  Output: $OUTPUT"

    if [[ "$MODEL" == *"DNA"* ]]; then
        echo "  Running DNA model..."
        bash "$RAXML_DNA_SCRIPT" "$ACCESSION" "$FASTA" "$OUTPUT" "$RAXML_EXEC"
    else
        echo "  Running RNA model: $MODEL"
        bash "$RAXML_RNA_SCRIPT" "$ACCESSION" "$FASTA" "$SS" "$OUTPUT" "$MODEL" "$RAXML_EXEC"
    fi

    # Check result
    if [ -d "$OUTPUT" ]; then
        FILE_COUNT=$(ls -1 "$OUTPUT" 2>/dev/null | wc -l)
        echo "  Completed: $FILE_COUNT files generated"
    else
        echo "  ERROR: Output directory not created"
    fi
done

echo ""
echo "========================================"
echo "Batch $PBS_ARRAY_INDEX complete"
echo "Finished: $(date)"
echo "========================================"
"""

    pbs_path = os.path.join(config.DIR_PBS, f"{model}_submit.pbs")
    with open(pbs_path, 'w') as f:
        f.write(pbs_script)

    os.chmod(pbs_path, 0o755)

    print(f"\nCreated PBS script: {pbs_path}")
    estimate_su_usage(config, num_batches)

    return pbs_path


def estimate_su_usage(config: GadiConfig, num_jobs: int):
    """Estimate Service Unit usage."""
    h, m, s = map(int, config.WALLTIME_PER_JOB.split(':'))
    walltime_hours = h + m / 60 + s / 3600

    queue_factors = {'normal': 1.0, 'express': 3.0, 'copyq': 1.0, 'hugemem': 3.0}
    factor = queue_factors.get(config.QUEUE, 1.0)

    su_per_job = config.NCPUS_PER_JOB * walltime_hours * factor
    total_su = su_per_job * num_jobs

    print(f"\nEstimated SU usage:")
    print(f"  Jobs: {num_jobs}")
    print(f"  CPUs per job: {config.NCPUS_PER_JOB}")
    print(f"  Walltime: {config.WALLTIME_PER_JOB}")
    print(f"  Queue: {config.QUEUE} (factor: {factor}x)")
    print(f"  SU per job: {su_per_job:.1f}")
    print(f"  Maximum total SU: {total_su:.1f}")
    print(f"\n  (Actual usage will be less if jobs finish early)")


def check_status(handler: DatasetHandler, model: str):
    """Check processing status."""
    config = handler.config
    all_families = handler.get_family_list()
    model_dir = os.path.join(config.DIR_TREES, model)

    complete = 0
    incomplete = 0
    pending = 0
    no_input = 0

    for accession in all_families:
        fasta = handler.get_fasta_path(accession)
        ss = handler.get_ss_path(accession)

        if not (fasta and ss and handler.check_sequence_count(fasta)):
            no_input += 1
            continue

        output_dir = os.path.join(model_dir, accession)
        if os.path.isdir(output_dir):
            num_files = len(os.listdir(output_dir))
            if num_files >= config.EXPECTED_OUTPUT_FILES:
                complete += 1
            else:
                incomplete += 1
        else:
            pending += 1

    total = len(all_families)
    processable = total - no_input

    print(f"\n{'=' * 60}")
    print(f"Status: {config.DATASET} alignment, model {model}")
    print(f"{'=' * 60}")
    print(f"  Total families:     {total:5d}")
    print(f"  Missing input:      {no_input:5d}")
    print(f"  Processable:        {processable:5d}")
    print(f"  {'─' * 30}")
    if processable > 0:
        print(f"  Complete:           {complete:5d} ({100 * complete / processable:.1f}%)")
        print(f"  Incomplete:         {incomplete:5d} ({100 * incomplete / processable:.1f}%)")
        print(f"  Pending:            {pending:5d} ({100 * pending / processable:.1f}%)")
    print(f"{'=' * 60}")


def run_local(handler: DatasetHandler, model: str, limit: Optional[int] = None):
    """Run processing locally (for testing)."""
    config = handler.config
    pending = get_pending_families(handler, model)

    if limit:
        pending = pending[:limit]

    if not pending:
        print("No families to process!")
        return

    print(f"Processing {len(pending)} families locally...")

    # Create output directory
    model_dir = os.path.join(config.DIR_TREES, model)
    os.makedirs(model_dir, exist_ok=True)

    for i, accession in enumerate(pending, 1):
        fasta = handler.get_fasta_path(accession)
        ss = handler.get_ss_path(accession)
        output_dir = os.path.join(model_dir, accession)

        print(f"[{i}/{len(pending)}] {accession}")

        if not (fasta and ss):
            print(f"  Skipped: missing files")
            continue

        os.makedirs(output_dir, exist_ok=True)

        if 'DNA' in model:
            cmd = f"bash {config.RAXML_DNA_SCRIPT} {accession} {fasta} {output_dir} {config.RAXML_EXECUTABLE}"
        else:
            cmd = f"bash {config.RAXML_RNA_SCRIPT} {accession} {fasta} {ss} {output_dir} {model} {config.RAXML_EXECUTABLE}"

        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"  ERROR: {result.stderr[:200]}")
        else:
            num_files = len(os.listdir(output_dir)) if os.path.isdir(output_dir) else 0
            print(f"  Done: {num_files} files")


def parse_seed_sections(handler: DatasetHandler):
    """Parse all section files to create FASTA and SS files (seed alignment only)."""
    if not isinstance(handler, SeedAlignmentHandler):
        print("ERROR: --parse is only available for seed alignment")
        return

    config = handler.config

    # Create output directories
    os.makedirs(config.DIR_FASTA, exist_ok=True)
    os.makedirs(config.DIR_SS, exist_ok=True)

    families = handler.get_family_list()

    if not families:
        print(f"No section files found in {config.DIR_SECTION}")
        print(f"Run --split-seed first to create section files from Rfam.seed")
        return

    print(f"Parsing {len(families)} section files...")

    success = 0
    skipped = 0
    failed = 0

    for i, accession in enumerate(families, 1):
        fasta, ss = handler.parse_section_if_needed(accession)

        if fasta and ss:
            success += 1
        elif fasta:
            failed += 1
        else:
            skipped += 1

        if i % 100 == 0 or i == len(families):
            print(f"  Progress: {i}/{len(families)} (success: {success}, skipped: {skipped}, no SS: {failed})")

    print(f"\n{'=' * 50}")
    print(f"Parsing complete:")
    print(f"  Success:  {success}")
    print(f"  Skipped:  {skipped} (< {config.MIN_SEQUENCES} sequences)")
    print(f"  No SS:    {failed}")
    print(f"{'=' * 50}")
    print(f"\nFASTA files: {config.DIR_FASTA}")
    print(f"SS files:    {config.DIR_SS}")


def split_rfam_seed(handler: DatasetHandler, min_sequences: int = 4):
    """
    Split Rfam.seed file into individual section files.

    Workflow: Rfam.seed → section files (one per family)
    """
    if not isinstance(handler, SeedAlignmentHandler):
        print("ERROR: --split-seed is only available for seed alignment")
        return

    config = handler.config
    rfam_seed = config.RFAM_SEED

    if not os.path.exists(rfam_seed):
        print(f"ERROR: Rfam.seed not found at {rfam_seed}")
        print(f"Please download from Rfam and place at: {rfam_seed}")
        return

    # Create section directory
    os.makedirs(config.DIR_SECTION, exist_ok=True)

    print(f"Reading Rfam.seed from: {rfam_seed}")
    print(f"Output directory: {config.DIR_SECTION}")

    with open(rfam_seed, 'r', encoding='utf-8', errors='ignore') as f:
        raw_data = f.read()

    # Split into sections
    sections = raw_data.split('# STOCKHOLM 1.0')
    sections = [('# STOCKHOLM 1.0' + section).strip()
                for section in sections if section.strip()]

    print(f"Found {len(sections)} total sections in Rfam.seed")

    written = 0
    skipped_no_ac = 0
    skipped_no_ss = 0
    skipped_few_seqs = 0

    for i, section in enumerate(sections, 1):
        # Extract accession
        ac = None
        for line in section.split('\n'):
            if line.startswith('#=GF AC'):
                parts = line.split()
                if len(parts) > 2:
                    ac = parts[2].strip()
                    break

        if not ac:
            skipped_no_ac += 1
            continue

        # Check for secondary structure
        has_ss = False
        for line in section.split('\n'):
            if line.startswith('#=GC SS_cons'):
                ss = line.split(maxsplit=2)[2].strip() if len(line.split(maxsplit=2)) > 2 else ""
                # Check if it has actual paired bases
                if any(c in OPENING_BRACKETS for c in ss):
                    has_ss = True
                break

        if not has_ss:
            skipped_no_ss += 1
            continue

        # Count sequences
        seq_count = 0
        for line in section.split('\n'):
            if not line.startswith('#') and not line.startswith('/') and line.strip():
                parts = line.split()
                if len(parts) == 2:
                    seq_count += 1

        if seq_count < min_sequences:
            skipped_few_seqs += 1
            continue

        # Write section file
        section_path = os.path.join(config.DIR_SECTION, f"{ac}.section")
        with open(section_path, 'w') as f_out:
            f_out.write(section)
        written += 1

        if i % 500 == 0:
            print(f"  Progress: {i}/{len(sections)} (written: {written})")

    print(f"\n{'=' * 50}")
    print(f"Split complete:")
    print(f"  Written:          {written}")
    print(f"  No accession:     {skipped_no_ac}")
    print(f"  No structure:     {skipped_no_ss}")
    print(f"  < {min_sequences} sequences:    {skipped_few_seqs}")
    print(f"{'=' * 50}")
    print(f"\nSection files: {config.DIR_SECTION}")


def create_seed_table(handler: DatasetHandler):
    """
    Parse Rfam.seed and create a summary table with family metadata.

    Output: TSV file with columns [AC, ID, NSEQ, NSITES, SS]
    """
    if not isinstance(handler, SeedAlignmentHandler):
        print("ERROR: --create-table is only available for seed alignment")
        return

    config = handler.config
    rfam_seed = config.RFAM_SEED

    if not os.path.exists(rfam_seed):
        print(f"ERROR: Rfam.seed not found at {rfam_seed}")
        return

    print(f"Parsing Rfam.seed: {rfam_seed}")

    families = []

    with open(rfam_seed, "r", encoding="utf-8", errors="ignore") as f:
        ac, id_, sq, si, ss = None, None, None, None, None

        for line in f:
            if line.startswith("#=GF AC"):
                ac = line.split()[2]
            elif line.startswith("#=GF ID"):
                id_ = line.split(maxsplit=2)[2].strip()
            elif line.startswith("#=GF SQ"):
                sq = int(line.split()[2])
            elif line.startswith("#=GC SS_cons"):
                ss = line.split(maxsplit=2)[2].strip()
                if not any(c in OPENING_BRACKETS for c in ss):
                    ss = None
            elif not line.startswith('#') and line.strip():
                parts = line.split()
                if len(parts) == 2:
                    si = len(parts[1])

            if all(v is not None for v in [ac, id_, sq, si, ss]):
                families.append([ac, id_, sq, si, ss])
                ac, id_, sq, si, ss = None, None, None, None, None

    # Create DataFrame and save
    try:
        import pandas as pd
        df = pd.DataFrame(families, columns=["AC", "ID", "NSEQ", "NSITES", "SS"])
        df = df.drop_duplicates().dropna()
        df = df[df['NSEQ'] >= config.MIN_SEQUENCES]

        table_path = os.path.join(config.DIR_INPUTS, "rfam_families.tsv")
        df.to_csv(table_path, sep="\t", index=False)

        print(f"\n{'=' * 50}")
        print(f"Table created: {table_path}")
        print(f"  Total families: {len(df)}")
        print(f"  Columns: {list(df.columns)}")
        print(f"{'=' * 50}")
    except ImportError:
        print("ERROR: pandas is required for --create-table")
        print("Install with: pip install pandas")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Unified Rfam Phylogenetic Pipeline for NCI/Gadi",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # SEED ALIGNMENT WORKFLOW:
  # Step 1: Split Rfam.seed into section files
  python rfam_unified_pipeline.py seed --split-seed

  # Step 2: Parse sections to create FASTA and SS files
  python rfam_unified_pipeline.py seed --parse

  # Step 3: Setup PBS jobs
  python rfam_unified_pipeline.py seed --setup --model S16A

  # Optional: Create summary table
  python rfam_unified_pipeline.py seed --create-table

  # FULL ALIGNMENT (files pre-prepared):
  python rfam_unified_pipeline.py full --setup --model DNA_2

  # OTHER COMMANDS:
  python rfam_unified_pipeline.py seed --status --model S16A
  python rfam_unified_pipeline.py full --local --model S16A --limit 5
        """
    )

    parser.add_argument('dataset', choices=['seed', 'full'],
                        help="Dataset type: 'seed' for seed alignment, 'full' for full alignment")

    parser.add_argument('--model', '-m', default=None,
                        help="Substitution model (e.g., DNA_2, S16A, S16B) - required for --setup, --status, --local")

    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--setup', action='store_true',
                        help="Generate PBS job array scripts")
    action.add_argument('--status', action='store_true',
                        help="Check processing status")
    action.add_argument('--local', action='store_true',
                        help="Run locally (for testing)")
    action.add_argument('--split-seed', action='store_true',
                        help="Split Rfam.seed into section files (seed only)")
    action.add_argument('--parse', action='store_true',
                        help="Parse section files to create FASTA/SS (seed only)")
    action.add_argument('--create-table', action='store_true',
                        help="Create summary table from Rfam.seed (seed only)")
    action.add_argument('--list-families', action='store_true',
                        help="List all families in dataset")

    parser.add_argument('--limit', type=int, default=None,
                        help="Limit number of families (for --local testing)")

    parser.add_argument('--families-per-job', type=int, default=None,
                        help="Override families per job (default: 10)")

    parser.add_argument('--walltime', type=str, default=None,
                        help="Override walltime (e.g., '02:00:00')")

    args = parser.parse_args()

    # Get configuration
    config = get_config(args.dataset)

    # Apply overrides
    if args.families_per_job:
        config.FAMILIES_PER_JOB = args.families_per_job
    if args.walltime:
        config.WALLTIME_PER_JOB = args.walltime

    # Validate --model is provided for actions that need it
    if args.model is None and (args.setup or args.status or args.local):
        parser.error("--model is required for --setup, --status, and --local")

    # Get handler
    handler = get_handler(config)
    handler.create_directories()

    # Execute action
    print(f"\n{'═' * 60}")
    print(f"Rfam Unified Pipeline")
    print(f"Dataset: {args.dataset.upper()} ALIGNMENT")
    if args.model:
        print(f"Model: {args.model}")
    print(f"{'═' * 60}")

    if args.setup:
        job_list = create_job_list(handler, args.model)
        if job_list:
            pbs_script = create_pbs_script(handler, args.model, job_list)
            print(f"\n{'─' * 60}")
            print("To submit to Gadi:")
            print(f"  qsub {pbs_script}")
            print(f"{'─' * 60}")

    elif args.status:
        check_status(handler, args.model)

    elif args.local:
        run_local(handler, args.model, args.limit)

    elif args.split_seed:
        split_rfam_seed(handler, config.MIN_SEQUENCES)

    elif args.parse:
        parse_seed_sections(handler)

    elif args.create_table:
        create_seed_table(handler)

    elif args.list_families:
        families = handler.get_family_list()
        print(f"\nFound {len(families)} families:")
        for f in families[:20]:
            print(f"  {f}")
        if len(families) > 20:
            print(f"  ... and {len(families) - 20} more")


if __name__ == "__main__":
    main()