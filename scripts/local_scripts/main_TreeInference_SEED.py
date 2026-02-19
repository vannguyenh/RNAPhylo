#!/usr/bin/env python3
"""
Rfam Seed Alignment Phylogenetic Pipeline
==========================================

Processes Rfam seed alignments to infer phylogenetic trees under DNA and RNA models.
Parses Stockholm format files, extracts FASTA sequences and secondary structures,
converts pseudoknots, and runs RAxML for tree inference.

Usage:
    python rfam_phylo_pipeline.py

Configuration:
    Edit the CONFIG section below before running.
"""

import os
import subprocess
import logging
from datetime import datetime
from dataclasses import dataclass
from typing import Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION - Edit these settings before running
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class Config:
    """Pipeline configuration settings."""

    # Base directories
    DIR_WORKING: str = "/Users/u7875558/RNAPhylo/seedAlignment_AllModels/"

    # Input files
    RFAM_SEED: str = ""  # Set automatically from DIR_WORKING

    # RAxML settings
    RAXML_EXECUTABLE: str = "/Users/u7875558/tools/standard-RAxML-master/raxmlHPC"
    RAXML_SCRIPT: str = "/Users/u7875558/Documents/promotion/projects/projects_code/RNAPhylo/scripts/nci_scripts/bashFiles/raxml.sh"
    RAXML_RNA_SCRIPT: str = "/Users/u7875558/Documents/promotion/projects/projects_code/RNAPhylo/scripts/nci_scripts/bashFiles/raxmlP.sh"

    # Processing settings
    MIN_SEQUENCES: int = 4  # Minimum sequences required per family
    EXPECTED_OUTPUT_FILES: int = 50  # Expected number of output files per run
    MAX_WORKERS: int = 4  # Number of parallel workers

    DATE = "260218"
    # Derived paths (set in __post_init__)
    DIR_INPUTS: str = ""
    DIR_OUTPUTS: str = ""
    DIR_TREES: str = ""
    DIR_SECTION: str = ""
    DIR_FASTA: str = ""
    DIR_SS: str = ""
    DIR_LOGS: str = ""

    def __post_init__(self):
        """Set derived paths after initialization."""
        self.DIR_INPUTS = os.path.join(self.DIR_WORKING, "inputs")
        self.DIR_OUTPUTS = os.path.join(self.DIR_WORKING, "outputs")
        self.DIR_TREES = os.path.join(self.DIR_OUTPUTS, "inferred_trees")
        self.DIR_SECTION = os.path.join(self.DIR_INPUTS, "sections")
        self.DIR_FASTA = os.path.join(self.DIR_INPUTS, "fasta_files")
        self.DIR_SS = os.path.join(self.DIR_INPUTS, "ss_files")
        self.DIR_LOGS = os.path.join(self.DIR_WORKING, "logs")
        self.RFAM_SEED = os.path.join(self.DIR_INPUTS, "Rfam.seed")

    def create_directories(self):
        """Create all required directories."""
        for path in [self.DIR_INPUTS, self.DIR_OUTPUTS, self.DIR_TREES,
                     self.DIR_SECTION, self.DIR_FASTA, self.DIR_SS,
                     self.DIR_LOGS]:
            os.makedirs(path, exist_ok=True)


# Initialize configuration
CONFIG = Config()

# ══════════════════════════════════════════════════════════════════════════════
# CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════

MATCHING_BRACKETS = {')': '(', ']': '[', '}': '{', '>': '<'}
MATCHING_PSEUDOKNOTS = {'a': 'A', 'b': 'B', 'c': 'C', 'd': 'D'}
OPENING_BRACKETS = set(MATCHING_BRACKETS.values()) | set(MATCHING_PSEUDOKNOTS.values())

# Characters to convert to dots (unpaired bases)
UNPAIRED_CHARS = {':': '.', ',': '.', '_': '.', '~': '.', '-': '.'}
PSEUDOKNOT_TO_DOT = {
    'A': '.', 'a': '.', 'B': '.', 'b': '.',
    'C': '.', 'c': '.', 'D': '.', 'd': '.'
}


# ══════════════════════════════════════════════════════════════════════════════
# DATA CLASSES
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class RfamFamily:
    """Represents an Rfam family with its metadata."""
    accession: str
    identifier: str
    num_sequences: int
    num_sites: int
    secondary_structure: str


@dataclass
class ProcessedFamily:
    """Holds paths to processed files for an Rfam family."""
    accession: str
    fasta_path: Optional[str] = None
    ss_brackets_path: Optional[str] = None
    ss_dots_path: Optional[str] = None
    success: bool = False
    error_message: Optional[str] = None


# ══════════════════════════════════════════════════════════════════════════════
# PARSING FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def parse_rfam_seed_to_table(fp_seedfile: str, fp_table: str,
                             min_sequences: int = 4) -> pd.DataFrame:
    """
    Parse Rfam.seed file and create a summary table.

    Args:
        fp_seedfile: Path to Rfam.seed file
        fp_table: Output path for the table
        min_sequences: Minimum number of sequences required

    Returns:
        DataFrame with family metadata
    """
    families = []

    with open(fp_seedfile, "r", encoding="utf-8", errors="ignore") as f:
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
                # Check if structure has actual paired bases
                if not any(char in OPENING_BRACKETS for char in ss):
                    ss = None
            elif not line.startswith('#') and line.strip():
                parts = line.split()
                if len(parts) == 2:
                    si = len(parts[1])

            # When we have all fields, save the family
            if all(v is not None for v in [ac, id_, sq, si, ss]):
                families.append([ac, id_, sq, si, ss])
                ac, id_, sq, si, ss = None, None, None, None, None

    df = pd.DataFrame(families, columns=["AC", "ID", "NSEQ", "NSITES", "SS"])
    df = df.drop_duplicates().dropna()
    df = df[df['NSEQ'] >= min_sequences]
    df.to_csv(fp_table, sep="\t", index=False)

    logging.info(f"Created table with {len(df)} families at {fp_table}")
    return df


def extract_accession(section: str) -> Optional[str]:
    """Extract RF accession number from a Stockholm section."""
    for line in section.split('\n'):
        if line.startswith('#=GF AC'):
            parts = line.split()
            if len(parts) > 2:
                return parts[2].strip()
    return None


def parse_sections(fp_seedfile: str, rf_list: set, dir_section: str) -> int:
    """
    Split Rfam.seed file into individual family sections.

    Args:
        fp_seedfile: Path to Rfam.seed file
        rf_list: Set of RF accessions to extract
        dir_section: Output directory for section files

    Returns:
        Number of sections successfully written
    """
    with open(fp_seedfile, 'r', encoding='utf-8', errors='ignore') as f:
        raw_data = f.read()

    sections = raw_data.split('# STOCKHOLM 1.0')
    sections = [('# STOCKHOLM 1.0' + section).strip()
                for section in sections if section.strip()]

    count = 0
    for section in sections:
        ac = extract_accession(section)
        if ac in rf_list:
            section_path = os.path.join(dir_section, f"{ac}.section")
            with open(section_path, 'w') as f_out:
                f_out.write(section)
            logging.debug(f"Created section file: {section_path}")
            count += 1

    logging.info(f"Extracted {count} section files to {dir_section}")
    return count


# ══════════════════════════════════════════════════════════════════════════════
# SECONDARY STRUCTURE CONVERSION
# ══════════════════════════════════════════════════════════════════════════════

def convert_pseudoknots_to_dots(rna_structure: str) -> str:
    """
    Convert pseudoknot notation to dots (unpaired bases).

    Converts all pseudoknot characters (A/a, B/b, C/c, D/d) and
    other non-standard characters to dots.
    """
    conversion = {**UNPAIRED_CHARS, **PSEUDOKNOT_TO_DOT}
    return ''.join(conversion.get(char, char) for char in rna_structure)

# ══════════════════════════════════════════════════════════════════════════════
# FILE PROCESSING
# ══════════════════════════════════════════════════════════════════════════════

def parse_fasta_and_ss(fp_section: str, fasta_dir: str,
                       ss_dir: str, min_sequences: int = 4) -> ProcessedFamily:
    """
    Parse a section file to extract FASTA alignment and secondary structure.

    Args:
        fp_section: Path to .section file
        fasta_dir: Output directory for FASTA files
        ss_dir: Output directory for secondary structure files
        min_sequences: Minimum required sequences

    Returns:
        ProcessedFamily with paths to created files
    """
    with open(fp_section, "r") as f:
        lines = f.readlines()

    content = "".join(lines)
    accession = extract_accession(content)

    if not accession:
        return ProcessedFamily(
            accession="unknown",
            success=False,
            error_message="Could not extract accession"
        )

    result = ProcessedFamily(accession=accession)

    # Extract alignments (non-comment, non-empty lines)
    alignments = {}
    for line in lines:
        if line.startswith('#') or line.startswith('/') or not line.strip():
            continue
        parts = line.split()
        if len(parts) == 2:
            name, seq = parts
            # Avoid duplicate sequences
            if seq not in alignments.values():
                alignments[name] = seq

    if len(alignments) < min_sequences:
        result.error_message = f"Only {len(alignments)} sequences (need {min_sequences})"
        logging.info(f'{accession}: {result.error_message}')
        return result

    # Write FASTA file
    fasta_path = os.path.join(fasta_dir, f"{accession}.nodup.fa")
    with open(fasta_path, "w") as f:
        for name, seq in alignments.items():
            f.write(f">{name}\n{seq}\n")
    result.fasta_path = fasta_path

    # Extract secondary structure
    ss_cons = None
    for line in lines:
        if line.startswith("#=GC SS_cons"):
            ss_cons = line.split(maxsplit=2)[2].strip()
            break

    if not ss_cons:
        result.error_message = "No secondary structure found"
        logging.warning(f'{accession}: {result.error_message}')
        return result

    # Convert to dots (pseudoknots as unpaired)
    ss_dots = convert_pseudoknots_to_dots(ss_cons)

    # Check if structure has any paired bases
    if len(set(ss_dots)) == 1:
        result.error_message = "Secondary structure has only unpaired bases"
        logging.warning(f'{accession}: {result.error_message}')
        return result

    # Write dots structure
    ss_dots_path = os.path.join(ss_dir, f'{accession}.dots.ss')
    with open(ss_dots_path, 'w') as f:
        f.write(ss_dots + '\n')
    result.ss_dots_path = ss_dots_path

    result.success = True
    logging.info(f'{accession}: Successfully processed FASTA and SS files')
    return result

# ══════════════════════════════════════════════════════════════════════════════
# RAXML EXECUTION
# ══════════════════════════════════════════════════════════════════════════════

def run_command(command: str, timeout: int = None) -> tuple[int, str, str]:
    """
    Run a shell command and capture output.

    Args:
        command: Shell command to execute
        timeout: Optional timeout in seconds

    Returns:
        Tuple of (return_code, stdout, stderr)
    """
    try:
        result = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True,
            timeout=timeout
        )
        return result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        logging.error(f"Command timed out: {command[:100]}...")
        return -1, "", "Timeout"
    except Exception as e:
        logging.error(f"Command failed: {e}")
        return -1, "", str(e)


def run_raxml_for_family(accession: str, fasta_path: str, ss_path: str,
                         output_prefix: str, model: str, config: Config) -> bool:
    """
    Run RAxML for a single family.

    Args:
        accession: RF accession
        fasta_path: Path to FASTA alignment
        ss_path: Path to secondary structure file
        output_prefix: Output directory prefix
        model: Substitution model name
        config: Pipeline configuration

    Returns:
        True if successful, False otherwise
    """
    # Check if already completed
    if os.path.isdir(output_prefix):
        num_files = len(os.listdir(output_prefix))
        if num_files >= config.EXPECTED_OUTPUT_FILES:
            logging.debug(f'{accession}: Already completed ({num_files} files)')
            return True
        else:
            # Incomplete run - remove and restart
            logging.warning(f'{accession}: Incomplete run ({num_files} files), restarting')
            run_command(f"rm -rf {output_prefix}")

    # Create output directory
    os.makedirs(output_prefix, exist_ok=True)

    # Build and run command
    if 'DNA' in model:
        command = (
            f"bash {config.RAXML_SCRIPT} "
            f"{accession} {fasta_path} {output_prefix} "
            f"{config.RAXML_EXECUTABLE}"
        )
    else:
        command = (
            f"bash {config.RAXML_RNA_SCRIPT} "
            f"{accession} {fasta_path} {ss_path} {output_prefix} "
            f"{model} {config.RAXML_EXECUTABLE}"
        )

    logging.info(f'{accession}: Running RAxML with model {model}')
    returncode, stdout, stderr = run_command(command)

    if returncode != 0:
        logging.error(f'{accession}: RAxML failed with code {returncode}')
        if stderr:
            logging.error(f'{accession}: stderr: {stderr[:500]}')
        return False

    logging.info(f'{accession}: RAxML completed successfully')
    return True


def process_family(args: tuple) -> dict:
    """
    Process a single Rfam family (for parallel execution).

    Args:
        args: Tuple of (section_file, model, config)

    Returns:
        Dictionary with processing results
    """
    section_file, model, config = args

    # Get accession from filename (remove .section extension)
    basename = os.path.basename(section_file)
    accession = basename.replace('.section', '')

    result = {
        'accession': accession,
        'parsed': False,
        'raxml_run': False,
        'error': None
    }

    try:
        # Parse FASTA and secondary structure
        processed = parse_fasta_and_ss(
            section_file,
            config.DIR_FASTA,
            config.DIR_SS,
            config.MIN_SEQUENCES
        )

        if not processed.success:
            result['error'] = processed.error_message
            return result

        result['parsed'] = True

        # Run RAxML if we have all required files
        if processed.fasta_path and processed.ss_dots_path:
            output_prefix = os.path.join(config.DIR_TREES, model, accession)

            success = run_raxml_for_family(
                accession,
                processed.fasta_path,
                processed.ss_dots_path,
                output_prefix,
                model,
                config
            )
            result['raxml_run'] = success
        else:
            result['error'] = "Missing FASTA or SS file"

    except Exception as e:
        result['error'] = str(e)
        logging.exception(f'{accession}: Unexpected error')

    return result


# ══════════════════════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def setup_logging(config: Config) -> str:
    """Configure logging and return log file path."""
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    logpath = os.path.join(config.DIR_LOGS, f"{timestamp}.log")
    os.makedirs(config.DIR_LOGS, exist_ok=True)

    # Configure root logger
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        handlers=[
            logging.FileHandler(logpath),
            logging.StreamHandler()  # Also print to console
        ]
    )

    return logpath


def main():
    """Main pipeline entry point."""
    # Get model from user
    print("\n" + "=" * 60)
    print("Rfam Seed Alignment Phylogenetic Pipeline")
    print("=" * 60 + "\n")

    model = input("Enter substitution model name: ").strip()
    if not model:
        print("Error: Model name required")
        return

    # Setup
    config = CONFIG
    config.create_directories()
    logpath = setup_logging(config)

    logging.info(f"Starting pipeline with model: {model}")
    logging.info(f"Log file: {logpath}")
    logging.info(f"Working directory: {config.DIR_WORKING}")
    logging.info(f"Max workers: {config.MAX_WORKERS}")

    # Get section files
    section_files = [
        os.path.join(config.DIR_SECTION, f)
        for f in os.listdir(config.DIR_SECTION)
        if f.endswith('.section')
    ]

    if not section_files:
        logging.error(f"No section files found in {config.DIR_SECTION}")
        return

    logging.info(f"Found {len(section_files)} section files to process")

    # Create output directory for this model
    model_output_dir = os.path.join(config.DIR_TREES, model)
    os.makedirs(model_output_dir, exist_ok=True)

    # Process families in parallel
    results = {
        'total': len(section_files),
        'parsed': 0,
        'raxml_success': 0,
        'errors': []
    }

    # Prepare arguments for parallel processing
    args_list = [(sf, model, config) for sf in section_files]

    with ProcessPoolExecutor(max_workers=config.MAX_WORKERS) as executor:
        futures = {executor.submit(process_family, args): args[0]
                   for args in args_list}

        for i, future in enumerate(as_completed(futures), 1):
            section_file = futures[future]
            try:
                result = future.result()

                if result['parsed']:
                    results['parsed'] += 1
                if result['raxml_run']:
                    results['raxml_success'] += 1
                if result['error']:
                    results['errors'].append((result['accession'], result['error']))

                # Progress update every 10 families
                if i % 10 == 0 or i == len(section_files):
                    logging.info(f"Progress: {i}/{len(section_files)} "
                                 f"({results['raxml_success']} successful)")

            except Exception as e:
                logging.error(f"Failed to process {section_file}: {e}")
                results['errors'].append((os.path.basename(section_file), str(e)))

    # Summary
    logging.info("\n" + "=" * 60)
    logging.info("PIPELINE COMPLETE")
    logging.info("=" * 60)
    logging.info(f"Total families:     {results['total']}")
    logging.info(f"Successfully parsed: {results['parsed']}")
    logging.info(f"RAxML successful:    {results['raxml_success']}")
    logging.info(f"Errors:              {len(results['errors'])}")

    if results['errors']:
        logging.info("\nErrors encountered:")
        for accession, error in results['errors'][:20]:  # Show first 20
            logging.info(f"  {accession}: {error}")
        if len(results['errors']) > 20:
            logging.info(f"  ... and {len(results['errors']) - 20} more")

    print(f"\nPipeline complete. See log file: {logpath}")


if __name__ == "__main__":
    main()