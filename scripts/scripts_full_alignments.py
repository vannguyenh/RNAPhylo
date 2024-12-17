# This script is created to convert full alignment file (.sto format)
# into a fasta file and a file containing the consensus secondary structure
# The two files will be required for the phylogenetic inference later.

# This script works only with the assigned AC (RF) value in Rfam database.

# PATTERN OF .sto FILE:
# #STOCKHOLM 1.0
# #=GF AU Infernal 1.1.4
# #=GS ...: extension for the sequence names as well as the methods to acquire the data
# Sequences (without '#' at the start)
# #=GC SS_cons: the consensus secondary structure

import os
from os.path import join
import subprocess
import logging
from datetime import datetime

MATCHING_BRACKETS = {')': '(', ']': '[', '}': '{', '>': '<'}
MATCHING_PSEUDOKNOTS = {'a': 'A', 'b': 'B', 'c': 'C', 'd': 'D'}
DIR_DATA = '/Users/u7875558/Documents/PhD/RNAPhylo/data/'
DIR_OUTPUTS = '/Users/u7875558/Documents/PhD/RNAPhylo/output_analysis'

"""
def create_directory(dir_base, dir_subs):
    # create new directories to save outputs
    # dir_base: the parent path
    # dir_subs: the subdirectories created under the parent path
    paths = dict()
    for dir_sub in dir_subs:
        path = join(dir_base, dir_sub)
        os.makedirs(path, exist_ok=True)
        paths[dir_sub] = path
    return paths
"""


def create_directory(dir_base, dir_sub):
    path = join(dir_base, dir_sub)
    os.makedirs(path, exist_ok=True)
    return path


def parseFastaAndSS(fi_stoFile, fp_fasta, fp_ss, ac):
    with open(fi_stoFile, 'r') as f:
        lines = f.readlines()

    alignments = dict()
    # sequences lines do not start with any '#' and '/'
    sequences = [l for l in lines if not l.startswith('#') and not l.startswith('/') and l.strip()]

    for line in sequences:
        parts = line.split()
        if len(parts) == 2 and parts[1] not in alignments.values():
            name = parts[0]
            seq = parts[1]
            alignments[name] = seq

    # fasta and ss outputs are created only when the aligments have the size >=4
    if len(alignments) >= 4:
        # produce the fasta file
        nodup_fasta_file = join(fp_fasta, f'{ac}.nodup.fa')
        with open(nodup_fasta_file, 'w') as fasta_file:
            for name, seq in alignments.items():
                fasta_file.write(f'>{name}\n{seq}\n')

        # produce the consensus secondary structure file
        ss_cons = ''
        ss_cons_found = False

        for line in lines:
            if line.startswith('#=GC SS_cons') and not ss_cons_found:
                ss_cons += line.split(maxsplit=2)[2].strip()
                ss_cons_found = True

        if ss_cons_found:
            # ignore pseudoknots and convert all of them into unpaired bases
            ss_convert_pseudo_unpaired = convertPseudoknotstoUnpairedBases(ss_cons)
            if len(set(ss_convert_pseudo_unpaired)) != 1:
                ss_convert_pseudo_unpaired_file = os.path.join(fp_ss, f'{ac}.pseudo.unpaired.ss')
                with open(ss_convert_pseudo_unpaired_file, 'w') as ss_file:
                    ss_file.write(ss_convert_pseudo_unpaired + '\n')

                ss_convert_pseudo_paired_file = os.path.join(fp_ss, f"{ac}.pseudo.paired.ss")
                ss_convert_pseudo_paired = convertPseudoknotstoBrackets(ss_cons)[0]
                if len(ss_convert_pseudo_paired) != 0:
                    with open(ss_convert_pseudo_paired_file, "w") as ss_file:
                        ss_file.write(ss_convert_pseudo_paired + "\n")
                else:
                    logging.error(f'No conversion occurred for {ac}.')
            else:
                logging.warning(f'{ac}: The secondary structure of {ac} has only unpaired bases.')
        else:
            logging.warning(f'{ac}: No secondary structure found.')

        if ss_convert_pseudo_unpaired == ss_convert_pseudo_paired:
            logging.info(f'{ac} does not have pseudoknots.')
        logging.info(f'Fasta and converted secondary structure files of {ac} are created. ')
        return nodup_fasta_file, ss_convert_pseudo_unpaired_file, ss_convert_pseudo_paired_file
    else:
        logging.info(f'Number of sequences of {ac} is less than 4.')
        return None, None, None


def convertPseudoknotstoUnpairedBases(rna_structure):
    conversion = {
        ':': '.', ',': '.', '_': '.', '~': '.', '-': '.',
        'A': '.', 'a': '.', 'B': '.', 'b': '.', 'C': '.', 'c': '.', 'D': '.', 'd': '.'
    }
    return ''.join(conversion.get(char, char) for char in rna_structure)


def convertPseudoknotstoBrackets(rna_structure):
    conversion = {':': '.', ',': '.', '_': '.', '~': '.', '-': '.'}
    pos_opening_pseudoknots = {'A': [], 'B': [], 'C': [], 'D': []}
    pos_closing_pseudoknots = {'a': [], 'b': [], 'c': [], 'd': []}

    for pos, char in enumerate(rna_structure):
        if char in pos_opening_pseudoknots:
            pos_opening_pseudoknots[char].append(pos)
        elif char in pos_closing_pseudoknots:
            pos_closing_pseudoknots[char].append(pos)

    pos_opening_pseudoknots = {k: v for k, v in pos_opening_pseudoknots.items() if len(v) != 0}
    pos_closing_pseudoknots = {k: v for k, v in pos_closing_pseudoknots.items() if len(v) != 0}
    srt_pos_opening_pseudoknots = {o: sorted(pos_opening_pseudoknots[o])[0] for o in pos_opening_pseudoknots}
    srt_pos_closing_pseudoknots = {c: sorted(pos_closing_pseudoknots[c], reverse=True)[0] for c in
                                   pos_closing_pseudoknots}
    srt_pos_opening_pseudoknots = {k: v for k, v in sorted(srt_pos_opening_pseudoknots.items(), key=lambda val: val[1])}

    converted_structure = ''.join(conversion.get(char, char) for char in rna_structure)

    for o in srt_pos_opening_pseudoknots:
        for c in MATCHING_PSEUDOKNOTS.keys():
            if o == MATCHING_PSEUDOKNOTS[c]:
                end_pos = srt_pos_closing_pseudoknots[c]
                start_pos = srt_pos_opening_pseudoknots[o]

                sub = converted_structure[start_pos:end_pos + 1]
                used_brackets = set(
                    [char for char in sub if char in MATCHING_BRACKETS or char in MATCHING_BRACKETS.values()])
                new_matching_bracket = {k: v for k, v in MATCHING_BRACKETS.items() if
                                        k not in used_brackets and v not in used_brackets}

                if len(new_matching_bracket) != 0:
                    k, v = list(new_matching_bracket.items())[0]
                    conversion[c] = k
                    conversion[o] = v
                else:
                    fully_closed_brackets = find_fully_closed_brackets(sub)
                    if fully_closed_brackets:
                        k, v = fully_closed_brackets[0]
                        conversion[c] = k
                        conversion[o] = v
                    else:
                        logging.error(f'No brackets available for pseudoknot {o}{c} in {sub}')

            converted_structure = ''.join(conversion.get(char, char) for char in converted_structure)

    if validateBracketsWithPositions(rna_structure)[1] != validateBracketsWithPositions(converted_structure)[1]:
        logging.error(f'The converted structure does not have the same bracket positions as the original one.')
        return '', conversion
    else:
        return converted_structure, conversion


def validateBracketsWithPositions(rna_structure):
    stacks = {'(': [], '[': [], '{': [], '<': [], 'A': [], 'B': [], 'C': [], 'D': []}
    pairing_positions = {'(': [], '[': [], '{': [], '<': [], 'A': [], 'B': [], 'C': [], 'D': []}

    for i, char in enumerate(rna_structure):
        if char in MATCHING_BRACKETS.values():
            stacks[char].append(i)
        elif char in MATCHING_BRACKETS.keys():
            if stacks[MATCHING_BRACKETS[char]]:
                opening_index = stacks[MATCHING_BRACKETS[char]].pop()
                pairing_positions[MATCHING_BRACKETS[char]].append((opening_index, i))
            else:
                return False, stacks

        elif char in MATCHING_PSEUDOKNOTS.values():
            stacks[char].append(i)
        elif char in MATCHING_PSEUDOKNOTS.keys():
            if stacks[MATCHING_PSEUDOKNOTS[char]]:
                opening_index = stacks[MATCHING_PSEUDOKNOTS[char]].pop()
                pairing_positions[MATCHING_PSEUDOKNOTS[char]].append((opening_index, i))
            else:
                return False, stacks

    all_matched = all(len(stack) == 0 for stack in stacks.values())
    if not all_matched:
        return False, stacks

    flat_pairing_positions = []
    for positions in pairing_positions.values():
        flat_pairing_positions.extend(positions)

    return True, sorted(flat_pairing_positions, key=lambda tup: tup[1])


def find_fully_closed_brackets(substring):
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


def run_command(command):
    process = subprocess.Popen(command, shell=True)
    process.communicate()


def main():
    # IQTREE_EXECUTE = '/Users/u7875558/Documents/PhD/tools/iqtree-2.2.2.6-MacOSX/bin/iqtree2'
    RAXML_EXECUTE = '/Users/u7875558/Documents/PhD/tools/standard-RAxML-master/raxmlHPC'
    #log_filename = os.path.join(paths['logs'], f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log")
    log_path=create_directory(DIR_OUTPUTS, 'logs')
    log_filename=join(log_path, f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log")
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    overlapping_pseudo_rnas = ['RF01833']
    overlapping_nopseudo_rnas = ['RF00740', 'RF00341', 'RF01899', 'RF00958', 'RF00872', 'RF00987', 'RF04090', 'RF02613',
                                 'RF00319', 'RF02451', 'RF00758', 'RF00877', 'RF01002', 'RF04121', 'RF03283', 'RF01007',
                                 'RF01383', 'RF03414', 'RF00723', 'RF02725', 'RF01038', 'RF00870', 'RF04290', 'RF01903',
                                 'RF00119', 'RF04076', 'RF04250', 'RF01014', 'RF00433', 'RF00769', 'RF03029', 'RF00807',
                                 'RF02610', 'RF02724', 'RF00725', 'RF03292', 'RF01902']

    # the .sto files are not all downloaded
    # This section is aimed to download only required files
    dir_sto = create_directory(DIR_DATA, 'full_alignments')

    for rf in overlapping_nopseudo_rnas:
        if not os.path.isfile(join(dir_sto, f'{rf}.sto')):
            url = f'https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/full_alignments/{rf}.sto'
            subprocess.run(['wget', url, '-P', dir_sto])

        sto_file = join(dir_sto, f'{rf}.sto')

        nodup_fasta, paired_pseudo_ss, unpaired_pseudo_ss = parseFastaAndSS(sto_file,
                                                                            create_directory(DIR_DATA, 'fasta'),
                                                                            create_directory(DIR_DATA, 'converted_ss'),
                                                                            rf)

        if nodup_fasta and paired_pseudo_ss and unpaired_pseudo_ss:
            raxml_prefix = join(DIR_OUTPUTS, 'raxml', rf)
            raxml_paired_pseudo_prefix = join(DIR_OUTPUTS, 'raxml_wPseu', rf)
            raxml_unpaired_pseudo_prefix = join(DIR_OUTPUTS, 'raxmlP_iPseu', rf)

            raxml_command = f"bash bashFiles/raxml.sh {rf} {nodup_fasta} {raxml_prefix} {RAXML_EXECUTE}"
            raxml_unpaired_pseudo_command = f"bash bashFiles/raxmlP.sh {rf} {nodup_fasta} {unpaired_pseudo_ss} {raxml_unpaired_pseudo_prefix} {RAXML_EXECUTE}"
            raxml_paired_pseudo_command = f"bash bashFiles/raxmlP.sh {rf} {nodup_fasta} {paired_pseudo_ss} {raxml_paired_pseudo_prefix} {RAXML_EXECUTE}"

            run_command(raxml_command)
            run_command(raxml_unpaired_pseudo_command)
            run_command(raxml_paired_pseudo_command)
        else:
            print(f'{rf} does not qualify for the analysis.')


if __name__ == "__main__":
    main()
