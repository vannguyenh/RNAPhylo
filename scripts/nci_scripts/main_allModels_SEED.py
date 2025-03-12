import os
import pandas as pd
import re
import subprocess
import logging
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor
from typing import List, Dict, Tuple, Optional

# ============================================================================
# Configuration and Constants
# ============================================================================

# Mapping for converting the structures
MATCHING_BRACKETS: Dict[str, str] = {')': '(', ']': '[', '}': '{', '>': '<'}
MATCHING_PSEUDOKNOTS: Dict[str, str] = {'a': 'A', 'b': 'B', 'c': 'C', 'd': 'D'}

# Logging format
LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'

# ============================================================================
# File and Directory Utilities
# ============================================================================


def create_directories(base_dir: str, sub_dirs: List[str]) -> Dict[str, str]:
    """
    Create subdirectories under a base directory and returns a mapping of subdirectory names to their paths.
    """
    paths = {}
    for sub_dir in sub_dirs:
        path = os.path.join(base_dir, sub_dir)
        os.makedirs(path, exist_ok=True)
        paths[sub_dir] = path
    return paths


def directory_has_enough_files(directory: str, expected_count: int) -> bool:
    """
    Checks if the given directory contains exactly the expected number of entries.
    Uses os.scandir() for efficiency
    """
    try:
        count = sum(1 for _ in os.scandir(directory))
        return count == expected_count
    except Exception as e:
        logging.error(f"Error counting files in directory {directory}: {e}")
        return False

# ============================================================================
# Parsing Functions
# ============================================================================


def parse_data_table(fp_seedfile: str, fp_table: str) -> pd.DataFrame:
    """
    Parses the seed file and generates a table with columns: AC, ID, NSEQ, NSITES, SS_CONS
    :param fp_seedfile: path to the seed file
    :param fp_table: path to the table file
    :return: the table in TSV format
    """
    tdata = []
    # Combine pseudoknots and brackets mappings to determine valid opening symbols
    combined = {**MATCHING_PSEUDOKNOTS, **MATCHING_PSEUDOKNOTS}
    opening_brackets = set(combined.values())

    try:
        with open(fp_seedfile, "r", encoding="utf-8", errors="ignore") as f:
            ac = accession_id = sq = si = ss = None
            for line in f:
                line = line.strip()
                if line.startswith("#=GF AC"):
                    parts = line.split()
                    if len(parts) >= 3:
                        ac = parts[2]
                elif line.startswith("#=GF ID"):
                    parts = line.split(maxsplit=2)
                    if len(parts) >= 3:
                        accession_id = parts[2].strip()
                elif line.startswith("#=GF SQ"):
                    # parts = int(line.split()[2])
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            sq = int(parts[2])
                        except ValueError:
                            sq = None
                elif line.startswith("#=GC SS_cons"):
                    ss = line.split(maxsplit=2)[2].strip()
                    # Only accept if there is at least one base pair present
                    if not set(ss).intersection(opening_brackets):
                        ss = None
                elif not line.startswith('#') and line.strip():
                    parts = line.split()
                    if len(parts) == 2:
                        si = len(parts[1])

                if ac and accession_id and (sq is not None) and (si is not None) and ss:
                    tdata.append([ac, accession_id, sq, si, ss])
                    ac = accession_id = sq = si = ss = None
    except Exception as e:
        logging.error(f"Error parsing data table: {e}")

    df = pd.DataFrame(tdata, columns=["AC", "ID", "NSEQ", "NSITES", "SS"])
    df.drop_duplicates(inplace=True)
    df.dropna()
    df = df[df['NSEQ'] >= 4]
    df.to_csv(fp_table, sep="\t", index=False)
    logging.info(f"Table created at {fp_table}")
    return df


def extract_rf(fp_section: str) -> Optional[str]:
    """
    Extracts the RNA RF from a section
    """
    for line in fp_section.split('\n'):
        if line.startswith('#=GF AC'):
            parts = line.split()
            if len(parts) > 2:
                return parts[2].strip()
    return None


def parse_sections(fp_seedfile: str, rf_list: List[str], po_section: str) -> None:
    """
    Splits the seed file into sections according to the RNA family and write each to a separate file
    :param fp_seedfile: path to the seed file
    :param rf_list: list of all RF values of RNA families
    :param po_section: output path of section files
    """
    try:
        with open(fp_seedfile, 'r', encoding='utf-8', errors='ignore') as f:
            raw_data = f.read()
    except Exception as e:
        logging.error(f"Error reading file {fp_seedfile}: {e}")
        return

    # Split using regex to capture the deliminter line without losing it
    sections = re.split(r'(?=# STOCKHOLM 1\.0)', raw_data)
    sections = [section.strip() for section in sections if sections.strip()]

    for section in sections:
        rf = extract_rf(section)
        if rf in rf_list:
            section_path = os.path.join(po_section, f"{rf}.section")
            try:
                with open(section_path, 'w') as f_out:
                    f_out.write(section)
                logging.info(f"Section file created for {rf} at {section_path}")
            except Exception as e:
                logging.error(f"Error writing section file {section_path}: {e}")


def parse_fasta_and_ss(fp_section: str, po_fasta: str, po_ss: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Parses a section file to extract FASTA alignment andd the consensus secondary structure
    Writes a FASTA file (with duplicates removed) and 2 SS files (ignore pseudokntos and consider pseudoknots)
    :param fp_section: file path to the section file
    :param po_fasta: output path to the fasta file
    :param po_ss: output path to the structure file
    :return: a tuple of file paths: (nodup_fasta_path, ss_considerPseu_path, ss_ignorePseu_path)
    """
    try:
        with open(fp_section, "r") as file:
            lines = file.readlines()
    except Exception as e:
        logging.error(f"Error reading sectional file {fp_section}: {e}")
        return None, None, None

    rna = extract_rf("\n".join(lines))
    alignments: Dict[str, str] = {}

    for line in lines:
        if not line.startswith('#') and not line.startswith('/') and line.strip():
            parts = line.split()
            if len(parts) == 2 and parts[1] not in alignments.values():
                alignments[parts[0]] = parts[1]

    if len(alignments) < 4:
        logging.info(f'Number of sequences for {rna} is less than 4.')
        return None, None, None

    nodup_fasta_path = os.path.join(po_fasta, f"{rna}.noddup.fa")
    try:
        with open(nodup_fasta_path, "w") as fasta_file:
            for name, seq in alignments.items():
                fasta_file.write(f">{name}\n{seq}\n")
    except Exception as e:
        logging.error(f'Error writing FASTA file {nodup_fasta_path}: {e}')
        return None, None, None

    ss_cons = ""
    ss_cons_found = False
    for line in lines:
        if line.startswith("#=GC SS_cons") and not ss_cons_found:
            parts = line.split(maxsplit=2)
            if len(parts) >= 3:
                ss_cons += parts[2].strip()
                ss_cons_found = True

    ss_ignorePseu_path = ss_considerPseu_path = None
    if ss_cons_found:
        ss_ignorePseu_structure = convert_pseudoknots_to_unpaired(ss_cons)
        if len(set(ss_ignorePseu_structure)) != 1:
            ss_ignorePseu_path = os.path.join(po_ss, f"{rna}.ignorePseu.ss")
            try:
                with open(ss_ignorePseu_path, 'w', encoding='utf-8') as ss_file:
                    ss_file.write(ss_ignorePseu_structure + '\n')
            except Exception as e:
                logging.error(f'Error writing the structure (ignore pseudoknots) file {ss_ignorePseu_path}: {e}')
                ss_ignorePseu_path = None

            ss_considerPseu_path = os.path.join(po_ss, f'{rna}.considerPseu.ss')
            ss_considerPseu_structure = convert_pseudoknots_to_brackets(ss_cons, rna)[0]
            if ss_considerPseu_structure:
                try:
                    with open(ss_considerPseu_path, 'w', encoding='utf-8') as ss_file:
                        ss_file.write(ss_considerPseu_structure + '\n')
                except Exception as e:
                    logging.error(f'Error writing the structure (consider pseudoknots) file {ss_considerPseu_path}: {e}')
                    ss_considerPseu_path = None
            else:
                logging.error(f'No conversion occurred for {rna}.')
        else:
            logging.warning(f'{rna}: The secondary structure has only unpaired bases.')
    else:
        logging.warning(f'{rna}: No secondary structure found.')

    if ss_considerPseu_structure == ss_ignorePseu_structure and ss_considerPseu_path and ss_ignorePseu_path:
        logging.info(f'{rna} does not have pseudoknots.')
    logging.info(f'FASTA and secondary structure files for {rna} are created.')


# ============================================================================
# Structure Conversion Functions
# ============================================================================


def convert_pseudoknots_to_unpaired(rna_structure: str) -> str:
    """
    Converts pseudoknots in the structure to unpaired bases using str.translate
    """
    conversion_table = str.maketrans({
            ':': '.', ',': '.', '_': '.', '~': '.', '-': '.',
            'A': '.', 'a': '.', 'B': '.', 'b': '.', 'C': '.', 'c': '.', 'D': '.', 'd': '.'
            })
    return rna_structure.translate(conversion_table)


def convert_pseudoknots_to_brackets(rna_structure: str, rna: str) -> Tuple[str, Dict[str, str]]:
    """
    Converts pseudoknot characters to bracket notation
    Returns the converted structure without any alphabet of pseudokntos and the conversion mapping
    """
    conversion: Dict[str, str] = {':': '.', ',': '.', '_': '.', '~': '.', '-': '.'}
    pos_opening_pseudoknots: Dict[str, List[int]] = {k: [] for k in ['A', 'B', 'C', 'D']}
    pos_closing_pseudoknots: Dict[str, List[int]] = {k: [] for k in ['a', 'b', 'c', 'd']}

    for pos, char in enumerate(rna_structure):
        if char in pos_opening_pseudoknots:
            pos_opening_pseudoknots[char].append(pos)
        elif char in pos_closing_pseudoknots:
            pos_closing_pseudoknots[char].append(pos)

    srt_pos_opening = {k: sorted(v)[0] for k, v in pos_opening_pseudoknots.items() if v}
    srt_pos_closing = {k: sorted(v, reverse=True)[0] for k, v in pos_closing_pseudoknots.items() if v}

    # Convert string to list for in-place modifications:
    converted_list = list(rna_structure)
    for opening in sorted(srt_pos_opening, key=lambda k: srt_pos_opening[k]):
        # Find the maching closing pseudoknot for the opening
        for closing, open_char in MATCHING_PSEUDOKNOTS.items():
            if open_char == opening and closing in srt_pos_closing:
                start_pos = srt_pos_opening[opening]
                end_pos = srt_pos_closing[closing]
                substring = "".join(converted_list[start_pos:end_pos+1])
                used_brackets = set(ch for ch in substring if ch in MATCHING_BRACKETS or ch in MATCHING_BRACKETS.values())
                new_matching = {k: v for k, v in MATCHING_BRACKETS.items() if k not in used_brackets}
                if new_matching:
                    k, v = next(iter(new_matching.items()))
                    conversion[closing] = k
                    conversion[opening] = v
                    for i in range(start_pos, end_pos+1):
                        if converted_list[i] in (closing, opening):
                            converted_list[i] = conversion.get(converted_list[i], converted_list[i])
                else:
                    fully_closed = find_fully_closed_brackets(substring)
                    if fully_closed:
                        k, v = fully_closed[0]
                        conversion[closing] = k
                        conversion[opening] = v
                        for i in range(start_pos, end_pos+1):
                            if converted_list[i] in (closing, opening):
                                converted_list[i] = conversion.get(converted_list[i], converted_list[i])
                            else:
                                logging.error(f'No brackets available for pseudoknots {opening}{closing} in {substring}.')
    converted_structure = "".join(converted_list)
    valid_orig, pos_orig = validate_brackets_positions(rna_structure)
    valid_conv, pos_conv = validate_brackets_positions(converted_structure)

    if pos_orig != pos_conv:
        logging.error(f'{rna}: The converted structure does not have the same bracket positions as the original one.')
        return '', conversion
    return converted_structure, conversion


def validate_brackets_positions(rna_structure: str) -> Tuple[bool, List[Tuple[int, int]]]:
    """
    Validates the brackets positions in the rna structure before and after conversion
    Returns a tuple (is_valid, list of pairing positions)
    """
    stacks = {ch: [] for ch in list(MATCHING_BRACKETS.values()) + list(MATCHING_PSEUDOKNOTS.values())}
    pairing_positions: Dict[str, List[Tuple[int, int]]] = {ch: [] for ch in stacks.keys()}

    for i, char in enumerate(rna_structure):
        if char in MATCHING_BRACKETS.values():
            stacks[char].append(i)
        elif char in MATCHING_BRACKETS.keys():
            if stacks[MATCHING_BRACKETS[char]]:
                opening_index = stacks[MATCHING_BRACKETS[char]].pop()
                pairing_positions[MATCHING_BRACKETS[char]].append((opening_index, i))
            else:
                return False, []
        elif char in MATCHING_PSEUDOKNOTS.values():
            stacks[char].append(i)
        elif char in MATCHING_PSEUDOKNOTS.keys():
            if stacks[MATCHING_PSEUDOKNOTS[char]]:
                opening_index = stacks[MATCHING_PSEUDOKNOTS[char]].pop()
                pairing_positions[MATCHING_PSEUDOKNOTS[char]].append((opening_index, i))
            else:
                return False, []

    all_matched = all(len(stack) == 0 for stack in stacks.values())
    flat_positions = [pos for pairs in pairing_positions.values() for pos in pairs]
    return all_matched, sorted(flat_positions, key=lambda x: x[1])


def find_fully_closed_brackets(substring: str) -> List[Tuple[str, str]]:
    """
    Identifies fully closed brackets pairs in the given structure
    """
    stack = []
    fully_closed = set()
    for char in substring:
        if char in MATCHING_BRACKETS.values():
            stack.append(char)
        elif char in MATCHING_BRACKETS.keys():
            if stack and stack[-1] == MATCHING_BRACKETS[char]:
                stack.pop()
                fully_closed.add((char, MATCHING_BRACKETS[char]))
            else:
                stack.clear()
    return list(fully_closed)

# ============================================================================
# Command Building and Execution
# ============================================================================


def build_dna_command(rf: str, logfile: str, nodup_fasta: str, po_prefix: str, model: str, raxml: str) -> str:
    """
    Constructs the command for the DNA model run.
    """
    return (
            f"qsub -V -N raxml_{rf} -o {logfile} -e {logfile} "
            f"-l ncpus=12 -l mem=48gb -l walltime=48:00:00 -l wd -- "
            f"bash /scratch/dx61/vh5686/tmp/RNAPhylo/scripts/bashFiles/raxml.sh {rf} {nodup_fasta} {po_prefix} {model} {raxml}"
    )


def build_rna_command(rf: str, logfile: str, nodup_fasta: str, po_prefix: str, model: str, raxml: str, rna_structure: str) -> str:
    """
    Constructs the command for the RNA model run.
    """
    return (
        f"qsub -V -N raxmlP_rna_{rf} -o {logfile} -e {logfile} "
        f"-l ncpus=24 -l mem=96gb -l walltime=48:00:00 -l wd -- "
        f"bash /scratch/dx61/vh5686/tmp/RNAPhylo/scripts/bashFiles/raxmlP.sh {rf} {nodup_fasta} {rna_structure} {po_prefix} {model} {raxml}"
    )


def run_command(command: str) -> None:
    """
    Executes an external command using subprocess
    """
    try:
        process = subprocess.Popen(command, shell=True)
        process.communicate()
    except Exception as e:
        logging.error(f"Command failed: {command} with error: {e}")


# ============================================================================
# Main Processing Functions
# ============================================================================


def process_rna_family(rf: str, section_file: str, paths: Dict[str, str], model: str, raxml_exec: str, log_filename: str) -> None:
    """
    Processes one RNA family: extracts FASTA, SS files (ignore and consider pseudoknots), build command line,
    and executes the appropriate external commands.
    """
    nodup_fasta, considerPseu_ss, ignorePseu_ss = parse_fasta_and_ss(section_file, paths['inputs/fasta'],
                                                                     paths['inputs/ss_files'])
    if nodup_fasta and considerPseu_ss and ignorePseu_ss:
        raxml_prefix = os.path.join(paths['outputs'], model, 'raxml', rf)
        raxml_considerPseu_prefix = os.path.join(paths['outputs'], model, 'raxmlP_cPseu', rf)
        raxml_ignorePseu_prefix = os.path.join(paths['outputs'], model, 'raxmlP_iPseu', rf)

        # ensure output directories exist
        os.makedirs(raxml_prefix, exist_ok=True)
        os.makedirs(raxml_considerPseu_prefix, exist_ok=True)
        os.makedirs(raxml_ignorePseu_prefix, exist_ok=True)

        # if not directory_has_enough_files(raxml_prefix, 50):
        #     logging.warning(f"{raxml_prefix} has not been run enough for 10 trees.")
        #     run_command(build_dna_command(rf, log_filename, nodup_fasta, raxml_prefix, model, raxml_exec))

        if considerPseu_ss != ignorePseu_ss:
            if not directory_has_enough_files(raxml_considerPseu_prefix, 50):
                run_command(build_rna_command(rf, log_filename, nodup_fasta, raxml_considerPseu_prefix, model, raxml_exec, considerPseu_ss))
            if not directory_has_enough_files(raxml_ignorePseu_prefix, 50):
                run_command(build_rna_command(rf, log_filename, nodup_fasta, raxml_ignorePseu_prefix, model, raxml_exec, ignorePseu_ss))
        else:
            if not directory_has_enough_files(raxml_ignorePseu_prefix, 50):
                run_command(build_rna_command(rf, log_filename, nodup_fasta, raxml_ignorePseu_prefix, model, raxml_exec, ignorePseu_ss))
    else:
        logging.info(f'RNA family {rf} skipped due to missing input files.')


def main() -> None:
    """
    Main driver function:
        - sets up directories and logging
        - parses the seed file to create a data table
        - splits the seed file into sections
        - processes each RNA family concurrently
    """
    # configure paths and constants
    RAXML_EXECUTE = '/scratch/dx61/vh5686/tmp/RNAPhylo/tools/standard-RAxML/raxmlHPC'
    RFAM_SEED = '/scratch/dx61/vh5686/tmp/RNAPhylo/full_models_SEED/inputs/Rfam.seed'
    DIR_WORKING = '/scratch/dx61/vh5686/tmp/RNAPhylo/full_models_SEED'
    MODEL = 'S16B'
    DIR_SUBS = ['inputs', 'outputs', 'logs',
                'inputs/fasta_files', 'inputs/ss_files', 'inputs/sections']
    paths = create_directories(DIR_WORKING, DIR_SUBS)

    # set up logging with a timestamped log file
    log_filename = os.path.join(paths['logs'], f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log")
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, format=LOG_FORMAT)
    logging.info('Running code with the model S16A!')

    # parse the seed file to generate the data table
    rfam_table = parse_data_table(RFAM_SEED, os.path.join(paths['inputs'], 'Rfam.full.seed.tbl'))
    rf_all = rfam_table["AC"].unique().tolist()

    # extract sections for each RNA entry
    parse_sections(RFAM_SEED, rf_all, paths['inputs/sections'])

    # process each RNA family cocurrently
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = []
        for rf in rf_all:
            section_file = os.path.join(paths['inputs/sections'], f"{rf}.section")
            futures.append(executor.submit(process_rna_family, rf, section_file, paths, MODEL, RAXML_EXECUTE, log_filename))
        for future in futures:
            future.result()  # wait for all tasks to finish


if __name__ == "__main__":
    main()
