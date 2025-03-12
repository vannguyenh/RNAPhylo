#!/usr/bin/env python3
"""
Script to download all full alignment files (.sto) for RNA families,
extract FASTA and secondary structure files, and generate a summary table.
This version downloads all .sto files in one go using wget.
"""

import os
from os.path import join
import pandas as pd
import logging
from datetime import datetime
import subprocess
from typing import Optional, Tuple, List, Dict

# Constants for structure conversion
MATCHING_BRACKETS: Dict[str, str] = {')': '(', ']': '[', '}': '{', '>': '<'}
MATCHING_PSEUDOKNOTS: Dict[str, str] = {'a': 'A', 'b': 'B', 'c': 'C', 'd': 'D'}

# -----------------------------------------------------------------------------
# Directory and Download Utilities
# -----------------------------------------------------------------------------

def create_directory(dir_base: str, section: str, dir_subs: Dict[str, str]) -> None:
    """Creates subdirectories under the specified section within the base directory."""
    for dir_sub in dir_subs.values():
        path = join(dir_base, section, dir_sub)
        os.makedirs(path, exist_ok=True)

def run_command(command: str) -> Optional[subprocess.CompletedProcess]:
    """Runs a shell command and logs the output."""
    logging.info(f"Running: {command}")
    try:
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logging.info(f"Output:\n{result.stdout}")
        if result.returncode != 0:
            logging.error(f"Command failed: {result.stderr}")
        return result
    except Exception as e:
        logging.exception(f"Failed to execute command: {command} - Error: {e}")
        return None

def download_all_sto_files(po_sto: str) -> None:
    """
    Downloads all .sto files from the Rfam full alignments directory using wget.
    All files are saved into the specified po_sto directory.

    This function does not work to download all files automatically without calling out the specific file.
    """
    base_url = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/full_alignments/"
    # wget -r -l1 -np -nd -A "*.sto" -R "index.html*" ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/full_alignments/
    cmd = f'wget -r -l1 -np -nd -P {po_sto} -A "*.sto" -R "index.html*" {base_url}'
    result = run_command(cmd)
    if result is not None and result.returncode == 0:
        logging.info("All .sto files downloaded successfully")
    else:
        logging.error("Error downloading .sto files")

# -----------------------------------------------------------------------------
# File Parsing and Conversion Functions
# -----------------------------------------------------------------------------

def extract_fasta_and_ss(fi_sto: str, po_fasta: str, po_ss_ignore: str, po_ss_consider: str, rf: str) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[int], Optional[int]]:
    """
    Extracts FASTA alignments and secondary structure from the .sto file.
    Writes a nodup FASTA file and secondary structure files (ignore and consider pseudoknots).
    Returns a tuple with file paths, number of sequences, and number of sites.
    """
    try:
        with open(fi_sto, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        logging.error(f"Error reading file {fi_sto}: {e}")
        return None, None, None, None, None

    alignments: Dict[str, str] = {}
    sequences = [l for l in lines if not l.startswith('#') and not l.startswith('/') and l.strip()]

    for line in sequences:
        parts = line.split()
        if len(parts) < 2:
            continue
        seq = parts[1].upper().replace('.', '-')
        if seq not in alignments.values():
            alignments[parts[0]] = seq

    if len(alignments) < 4:
        logging.info(f"Number of sequences of {rf} is less than 4.")
        return None, None, None, None, None

    nodup_fasta_file = join(po_fasta, f'{rf}.nodup.fa')
    try:
        with open(nodup_fasta_file, 'w') as fasta_file:
            for name, seq in alignments.items():
                fasta_file.write(f'>{name}\n{seq}\n')
    except Exception as e:
        logging.error(f"Error writing FASTA file {nodup_fasta_file}: {e}")
        return None, None, None, None, None

    ss_cons = ''
    ss_cons_found = False
    for line in lines:
        if line.startswith('#=GC SS_cons') and not ss_cons_found:
            parts = line.split(maxsplit=2)
            if len(parts) >= 3:
                ss_cons = parts[2].strip()
                ss_cons_found = True

    ss_ignore_file = ss_consider_file = None
    if ss_cons_found:
        ss_unpaired = convert_unpaired_bases(ss_cons)
        if len(set(ss_unpaired)) != 1:
            # Convert structure for ignoring pseudoknots
            ss_ignore = convert_ss_ignore_pseudo(ss_unpaired)
            ss_ignore_file = join(po_ss_ignore, f'{rf}.ss')
            try:
                with open(ss_ignore_file, 'w') as ss_file:
                    ss_file.write(ss_ignore + '\n')
            except Exception as e:
                logging.error(f"Error writing file {ss_ignore_file}: {e}")
                ss_ignore_file = None

            # If pseudoknot characters are present, convert for considering pseudoknots
            if any(pseudo in ss_unpaired for pseudo in MATCHING_PSEUDOKNOTS.values()):
                ss_consider, _ = convert_ss_consider_pseudo(ss_ignore)
                logging.info(f"{rf} has pseudoknots.")
                ss_consider_file = join(po_ss_consider, f'{rf}.pseudo.ss')
                try:
                    with open(ss_consider_file, 'w') as ss_pseudo_file:
                        ss_pseudo_file.write(ss_consider + '\n')
                except Exception as e:
                    logging.error(f"Error writing file {ss_consider_file}: {e}")
                    ss_consider_file = None
            else:
                ss_consider_file = None

            nseq = len(alignments)
            nsite = len(list(alignments.values())[0])
            return nodup_fasta_file, ss_ignore_file, ss_consider_file, nseq, nsite
        else:
            logging.warning(f"The secondary structure of {rf} has only unpaired bases.")
    else:
        logging.warning(f"No secondary structure found for {rf}.")

    return nodup_fasta_file, ss_ignore_file, ss_consider_file, len(alignments), len(list(alignments.values())[0]) if alignments else 0

def convert_unpaired_bases(rna_structure: str) -> str:
    """Converts various unpaired symbols in the RNA structure to a uniform format."""
    conversion = {':': '.', ',': '.', '_': '.', '~': '.', '-': '.'}
    return ''.join(conversion.get(char, char) for char in rna_structure)

def convert_ss_ignore_pseudo(rna_structure: str) -> str:
    """Converts pseudoknot characters to unpaired bases for the 'ignore pseudoknots' scenario."""
    conversion = {'A': '.', 'a': '.', 'B': '.', 'b': '.', 'C': '.', 'c': '.', 'D': '.', 'd': '.'}
    return ''.join(conversion.get(char, char) for char in rna_structure)

def convert_ss_consider_pseudo(rna_structure: str) -> Tuple[str, Dict[str, str]]:
    """
    Converts pseudoknot characters to bracket notation.
    Returns the converted structure and the mapping used for conversion.
    """
    conversion = {':': '.', ',': '.', '_': '.', '~': '.', '-': '.'}
    pos_opening: Dict[str, List[int]] = {k: [] for k in ['A', 'B', 'C', 'D']}
    pos_closing: Dict[str, List[int]] = {k: [] for k in ['a', 'b', 'c', 'd']}

    for pos, char in enumerate(rna_structure):
        if char in pos_opening:
            pos_opening[char].append(pos)
        elif char in pos_closing:
            pos_closing[char].append(pos)

    pos_opening = {k: v for k, v in pos_opening.items() if v}
    pos_closing = {k: v for k, v in pos_closing.items() if v}
    srt_opening = {o: sorted(pos_opening[o])[0] for o in pos_opening}
    srt_closing = {c: sorted(pos_closing[c], reverse=True)[0] for c in pos_closing}

    converted_list = list(rna_structure)
    for opening in sorted(srt_opening, key=lambda k: srt_opening[k]):
        for closing, open_char in MATCHING_PSEUDOKNOTS.items():
            if open_char == opening and closing in srt_closing:
                start_pos = srt_opening[opening]
                end_pos = srt_closing[closing]
                substring = "".join(converted_list[start_pos:end_pos+1])
                used_brackets = {ch for ch in substring if ch in MATCHING_BRACKETS or ch in MATCHING_BRACKETS.values()}
                new_matching = {k: v for k, v in MATCHING_BRACKETS.items() if k not in used_brackets and v not in used_brackets}
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
                        logging.error(f'No brackets available for pseudoknot {opening}{closing} in {substring}')
    converted_structure = "".join(converted_list)
    if validate_brackets_with_positions(rna_structure)[1] != validate_brackets_with_positions(converted_structure)[1]:
        logging.error("The converted structure does not have the same bracket positions as the original one.")
        return '', conversion
    return converted_structure, conversion

def find_fully_closed_brackets(substring: str) -> List[Tuple[str, str]]:
    """Finds fully closed bracket pairs in the given substring."""
    stack = []
    fully_closed = set()
    for char in substring:
        if char in MATCHING_BRACKETS.values():
            stack.append(char)
        elif char in MATCHING_BRACKETS:
            if stack and stack[-1] == MATCHING_BRACKETS[char]:
                stack.pop()
                fully_closed.add((char, MATCHING_BRACKETS[char]))
            else:
                stack.clear()
    return list(fully_closed)

def validate_brackets_with_positions(rna_structure: str) -> Tuple[bool, List[Tuple[int, int]]]:
    """
    Validates the bracket structure of the RNA structure.
    Returns a tuple (is_valid, list of pairing positions).
    """
    stacks = {'(': [], '[': [], '{': [], '<': [], 'A': [], 'B': [], 'C': [], 'D': []}
    pairing_positions = {k: [] for k in stacks.keys()}

    for i, char in enumerate(rna_structure):
        if char in MATCHING_BRACKETS.values():
            stacks[char].append(i)
        elif char in MATCHING_BRACKETS:
            if stacks.get(MATCHING_BRACKETS[char]):
                opening_index = stacks[MATCHING_BRACKETS[char]].pop()
                pairing_positions[MATCHING_BRACKETS[char]].append((opening_index, i))
            else:
                return False, []
        elif char in MATCHING_PSEUDOKNOTS.values():
            stacks[char].append(i)
        elif char in MATCHING_PSEUDOKNOTS:
            if stacks.get(MATCHING_PSEUDOKNOTS[char]):
                opening_index = stacks[MATCHING_PSEUDOKNOTS[char]].pop()
                pairing_positions[MATCHING_PSEUDOKNOTS[char]].append((opening_index, i))
            else:
                return False, []
    all_matched = all(len(stack) == 0 for stack in stacks.values())
    flat_positions = [pos for pairs in pairing_positions.values() for pos in pairs]
    return all_matched, sorted(flat_positions, key=lambda x: x[1])

# -----------------------------------------------------------------------------
# Main Function
# -----------------------------------------------------------------------------

def main() -> None:
    DIR_WORKING = '/Users/u7875558/Documents/PhD/RNAPhylo/fullAlignments_FULL'
    input_subs = {
        'full_align': 'sto',
        'fasta': 'fasta',
        'ss_all': 'ss_all',
        'ss_consider_pseudo': 'ss_pseudo'
    }

    # Setup logging with proper filename and filemode
    log_path = os.path.join(DIR_WORKING, 'logs')
    os.makedirs(log_path, exist_ok=True)
    log_file = join(log_path, f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log")
    logging.basicConfig(filename=log_file, filemode='w', level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s')

    # Create directories (inputs/sto, inputs/fasta, etc.)
    create_directory(DIR_WORKING, 'inputs', input_subs)

    # Download all .sto files using wget in one go.
    sto_dir = join(DIR_WORKING, 'inputs', input_subs['full_align'])

    # Process each RF family by iterating through a known RF list.
    rf_list = [rf.split('.')[0] for rf in os.listdir(sto_dir)]
    tdata = []

    for rf in rf_list:
        sto_file = join(sto_dir, f'{rf}.sto')
        if os.path.isfile(sto_file):
            fasta_file, ss_ignore_file, ss_consider_file, nseq, nsite = extract_fasta_and_ss(
                sto_file,
                join(DIR_WORKING, 'inputs', input_subs['fasta']),
                join(DIR_WORKING, 'inputs', input_subs['ss_all']),
                join(DIR_WORKING, 'inputs', input_subs['ss_consider_pseudo']),
                rf
            )
            tdata.append([rf, nseq, nsite])
        else:
            tdata.append([rf, None, None])

    df = pd.DataFrame(tdata, columns=['RNA', 'NSEQ', 'NSITE'])
    df.drop_duplicates(inplace=True)
    df.dropna(inplace=True)
    output_tbl = join(DIR_WORKING, 'inputs', 'Rfam.fullAlign.tbl')
    df.to_csv(output_tbl, index=False)
    logging.info(f"Summary table of full alignments is created at {output_tbl}")


if __name__ == "__main__":
    main()
