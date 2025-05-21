import os
from os.path import join

import subprocess
import pandas as pd
import logging
from datetime import datetime

MATCHING_BRACKETS = {')': '(', ']': '[', '}': '{', '>': '<'}
MATCHING_PSEUDOKNOTS = {'a': 'A', 'b': 'B', 'c': 'C', 'd': 'D'}

# ── CONFIGURE THESE ONCE ──────────────────────────────────────────────────────
DIR_WORKING  = "/Users/u7875558/Documents/RNAPhylo/allModels_SEED/"
DIR_INPUTS   = os.path.join(DIR_WORKING, "inputs")
DIR_OUTPUTS  = os.path.join(DIR_WORKING, "outputs")

# the only three subfolders we use:
DIR_SECTION  = os.path.join(DIR_INPUTS, "sections")
DIR_FASTA    = os.path.join(DIR_INPUTS, "fasta_files")
DIR_SS       = os.path.join(DIR_INPUTS, "ss_files")

RAXML_EXEC   = "/Users/u7875558/tools/standard-RAxML/raxmlHPC"
RFAM_SEED    = os.path.join(DIR_INPUTS, "Rfam.seed")

RAXML_EXECUTE = '/Users/u7875558/tools/standard-RAxML/raxmlHPC'

# ────────────────────────────────────────────────────────────────────────────────

def create_directories(base_dir, sub_dirs):
    paths = {}
    for sub_dir in sub_dirs:
        path = os.path.join(base_dir, sub_dir)
        os.makedirs(path, exist_ok=True)
        paths[sub_dir] = path
    return paths

def create_directory(base_dir, sub_dir):
    path = join(base_dir, sub_dir)
    os.makedirs(path, exist_ok=True)
    return path

def parseData_Table(fp_seedfile, fp_table):
    tdata = []
    opening_brackets = set(dict(MATCHING_PSEUDOKNOTS.items() | MATCHING_BRACKETS.items()).values())
    with open(fp_seedfile, "r", encoding="utf-8", errors="ignore") as f:
        ac, id, sq, si, ss = None, None, None, None, None
        for line in f:
            if line.startswith("#=GF AC"):
                ac = line.split()[2]
            elif line.startswith("#=GF ID"):
                id=line.split(maxsplit=2)[2].strip()
            elif line.startswith("#=GF SQ"):
                sq=int(line.split()[2])
            elif line.startswith("#=GC SS_cons"):
                ss=line.split(maxsplit=2)[2].strip()
                char_ss=set(ss)
                if len(char_ss.intersection(opening_brackets)) == 0:
                    ss=None
            elif not line.startswith('#') and line.strip():
                parts = line.split()
                if len(parts) == 2:
                    si = len(parts[1])

            if ac and id and sq and si and ss:
                tdata.append([ac, id, sq, si, ss])
                ac, id, sq, si, ss = None, None, None, None, None

    df = pd.DataFrame(tdata, columns=["AC", "ID", "NSEQ", "NSITES", "SS"])
    df.drop_duplicates(inplace=True)
    df.dropna()
    df = df[df['NSEQ'] >= 4]
    df.to_csv(fp_table, sep="\t", index=False)
    logging.info(f"Table created at {fp_table}")
    return df

def parse_sections(fp_seedfile, rf_list, dir_section):
    with open(fp_seedfile, 'r', encoding='utf-8', errors='ignore') as f:
        raw_data = f.read()

    sections = raw_data.split('# STOCKHOLM 1.0')
    sections = [('# STOCKHOLM 1.0' + section).strip() for section in sections if section.strip()]

    for section in sections:
        ac = extract_ac(section)
        if ac in rf_list:
            section_path = os.path.join(dir_section, f"{ac}.section")
            with open(section_path, 'w') as f_out:
                f_out.write(section)
            logging.info(f"Section file created for {ac} at {section_path}")

def extract_ac(section):
    lines = section.split('\n')
    for line in lines:
        if line.startswith('#=GF AC'):
            parts = line.split()
            if len(parts) > 2:
                return parts[2].strip()
    return None

def parseFastaAndSS(fp_section, fasta_dir, ss_dir):
    with open(fp_section, "r") as file:
        lines = file.readlines()
        ac = extract_ac("\n".join(lines))
    alignments = dict()
    sequences = [l for l in lines if not l.startswith('#') and not l.startswith('/') and l.strip()]

    for line in sequences:
        parts = line.split()
        if len(parts) == 2 and parts[1] not in alignments.values():
            alignments[parts[0]] = parts[1]

    if len(alignments) >=4:
        fasta_nodup_path = os.path.join(fasta_dir, f"{ac}.nodup.fa")
        with open(fasta_nodup_path, "w") as fasta_file:
            for name, seq in alignments.items():
                fasta_file.write(f">{name}\n{seq}\n")

        ss_cons = ""
        ss_cons_found = False

        for line in lines:
            if line.startswith("#=GC SS_cons") and not ss_cons_found:
                ss_cons += line.split(maxsplit=2)[2].strip()
                ss_cons_found = True

        if ss_cons_found:
            ss_dots_structure = convertPseudoknotstoUnpairedBases(ss_cons)
            if len(set(ss_dots_structure)) != 1:
                ss_dots_path = os.path.join(ss_dir, f'{ac}.dots.ss')
                with open(ss_dots_path, 'w') as ss_file:
                    ss_file.write(ss_dots_structure + '\n')

                ss_brackets_path = os.path.join(ss_dir, f"{ac}.brackets.ss")
                ss_brackets_structure = convertPseudoknotstoBrackets(ss_cons, ac)[0]
                if len(ss_brackets_structure) != 0:
                    with open(ss_brackets_path, "w") as ss_file:
                        ss_file.write(ss_brackets_structure + "\n")
                else:
                    logging.error(f'No conversion occurred for {ac}.')
            else:
                logging.warning(f'{ac}: The secondary structure of {ac} has only unpaired bases.')
        else:
            logging.warning(f'{ac}: No secondary structure found.')
        if ss_dots_structure == ss_brackets_structure:
            logging.info(f'{ac} does not have pseudoknots.')
        logging.info(f'Fasta and converted secondary structure files of {ac} are created. ')
        return fasta_nodup_path, ss_brackets_path, ss_dots_path
    else:
        logging.info(f'Number of sequences of {ac} is less than 4.')
        return None, None, None

def convertPseudoknotstoUnpairedBases(rna_structure):
    conversion = {
            ':': '.', ',': '.', '_': '.', '~': '.', '-': '.',
            'A': '.', 'a': '.', 'B': '.', 'b': '.', 'C': '.', 'c': '.', 'D': '.', 'd': '.'
            }
    return ''.join(conversion.get(char, char) for char in rna_structure)

def convertPseudoknotstoBrackets(rna_structure, ac):
    conversion = {':': '.', ',': '.', '_': '.', '~': '.', '-': '.'}
    pos_opening_pseudoknots = {'A': [], 'B': [], 'C': [], 'D': []}
    pos_closing_pseudoknots = {'a': [], 'b': [], 'c': [], 'd': []}

    for pos, char in enumerate(rna_structure):
        if char in pos_opening_pseudoknots:
            pos_opening_pseudoknots[char].append(pos)
        elif char in pos_closing_pseudoknots:
            pos_closing_pseudoknots[char].append(pos)

    pos_opening_pseudoknots = {k:v for k,v in pos_opening_pseudoknots.items() if len(v) !=0}
    pos_closing_pseudoknots = {k:v for k,v in pos_closing_pseudoknots.items() if len(v) !=0}
    srt_pos_opening_pseudoknots = {o : sorted(pos_opening_pseudoknots[o])[0] for o in pos_opening_pseudoknots}
    srt_pos_closing_pseudoknots = {c : sorted(pos_closing_pseudoknots[c], reverse=True)[0] for c in pos_closing_pseudoknots}
    srt_pos_opening_pseudoknots = {k:v for k, v in sorted(srt_pos_opening_pseudoknots.items(), key=lambda val:val[1])}

    converted_structure = ''.join(conversion.get(char, char) for char in rna_structure)

    for o in srt_pos_opening_pseudoknots:
        for c in MATCHING_PSEUDOKNOTS.keys():
            if o == MATCHING_PSEUDOKNOTS[c]:
                end_pos = srt_pos_closing_pseudoknots[c]
                start_pos = srt_pos_opening_pseudoknots[o]

                sub = converted_structure[start_pos:end_pos+1]
                used_brackets = set([char for char in sub if char in MATCHING_BRACKETS or char in MATCHING_BRACKETS.values()])
                new_matching_bracket = {k: v for k, v in MATCHING_BRACKETS.items() if k not in used_brackets and v not in used_brackets}

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

    if validatBracketsWithPositions(rna_structure)[1] != validatBracketsWithPositions(converted_structure)[1]:
        logging.error(f'{ac}: The converted structure does not have the same bracket positions as the original one.')
        return '', conversion
    else:
        return converted_structure, conversion

def validatBracketsWithPositions(rna_structure):
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
    MODEL = 'S16B'
    # set up logging
    logpath = os.path.join(DIR_WORKING, "logs",
                           datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + ".missing.log")
    os.makedirs(os.path.dirname(logpath), exist_ok=True)
    logging.basicConfig(filename=logpath, level=logging.INFO,
                        format="%(asctime)s  %(levelname)s  %(message)s")


    #rfam_table = parseData_Table(RFAM_SEED, join(DIR_INPUTS, 'Rfam.full.seed.tbl'))
    #rfam_table = pd.read_table(join(DIR_INPUTS, 'Rfam.full.seed.tbl'))
    #RFs = rfam_table["AC"].unique().tolist()
    RFs = ['RF03132', 'RF03156', 'RF00839', 'RF03479', 'RF01038']
    for rf in RFs:
        section_file = os.path.join(DIR_SECTION, f"{rf}.section")
        nodup_fasta, brackets_ss, dots_ss = parseFastaAndSS(section_file, DIR_FASTA, DIR_SS)

        if nodup_fasta and brackets_ss and dots_ss:
            #iqtree_prefix = os.path.join(paths['outputs'], 'iqtree', rf)
            #raxml_prefix = os.path.join(DIR_OUTPUTS, MODEL, 'raxml', rf)
            #raxml_brackets_prefix = os.path.join(DIR_OUTPUTS, MODEL, 'raxmlP_wPseu', rf)
            raxml_dots_prefix = os.path.join(DIR_OUTPUTS, MODEL, 'raxmlP_iPseu', rf)

            # if not os.path.isdir(raxml_prefix) or len(os.listdir(raxml_prefix)) != 50:
            #     run_command(f"rm -r {raxml_prefix}")
            #     logging.warning(f"{raxml_prefix} has not been run or not run enough for 10 trees.")
            #     raxml_command = (
            #         #f"qsub -V -N raxml_{rf} -o {log_filename} -e {log_filename} "
            #         #f"-l ncpus=12 -l mem=48gb -l walltime=48:00:00 -l wd -- "
            #         f"bash /Users/u7875558/Documents/PhD/RNAPhylo/scripts/nci_scripts/bashFiles/raxml.sh {rf} {nodup_fasta} {raxml_prefix} {RAXML_EXECUTE}"
            #     )
            #     run_command(raxml_command)

            # if not os.path.isdir(raxml_brackets_prefix) or len(os.listdir(raxml_brackets_prefix)) != 50:
            #     run_command(f"rm -r {raxml_brackets_prefix}")
            #     raxml_brackets_command = (
            #         #f"qsub -V -N raxmlP_br_{rf} -o {paths['logs']} -e {paths['logs']} "
            #         #f"-l ncpus=24 -l mem=96gb -l walltime=48:00:00 -l wd -- "
            #         f"bash /Users/u7875558/Documents/PhD/RNAPhylo/scripts/nci_scripts/bashFiles/raxmlP.sh {rf} {nodup_fasta} {brackets_ss} {raxml_brackets_prefix} {MODEL} {RAXML_EXECUTE}"
            #     )
            #     run_command(raxml_brackets_command)
            #
            if not os.path.isdir(raxml_dots_prefix) or len(os.listdir(raxml_dots_prefix)) !=50:
                run_command(f"rm -r {raxml_dots_prefix}")
                raxml_dots_command = (
                    #f"qsub -V -N raxmlP_dot_{rf} -o {paths['logs']} -e {paths['logs']} "
                    #f"-l ncpus=24 -l mem=96gb -l walltime=48:00:00 -l wd -- "
                    f"bash /Users/u7875558/Documents/RNAPhylo/RNAPhylo/scripts/nci_scripts/bashFiles/raxmlP.sh {rf} {nodup_fasta} {dots_ss} {raxml_dots_prefix} {MODEL} {RAXML_EXECUTE}"
               )
                run_command(raxml_dots_command)

if __name__ == "__main__":
    main()
