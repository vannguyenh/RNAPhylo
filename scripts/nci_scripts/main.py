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
from Bio import Phylo
import time
import requests

MATCHING_BRACKETS = {')': '(', ']': '[', '}': '{', '>': '<'}
MATCHING_PSEUDOKNOTS = {'a': 'A', 'b': 'B', 'c': 'C', 'd': 'D'}

# I cannot access the files I used to analyse in my home
# I got the list of RNA families locally and parse the list of RNA families into the code.

def create_directory(dir_base, project, dir_sub):
    path = join(dir_base, project, dir_sub)
    os.makedirs(path, exist_ok=True)
    return path

def create_directories(dir_base, project, dir_subs):
    paths = {}
    for dir_sub in dir_subs:
        path = os.path.join(dir_base, project, dir_sub)
        os.makedirs(path, exist_ok=True)
        paths[dir_sub] = path
    return paths

def parseFastaAndSS(fi_stoFile, fp_fasta, fp_ss, ac):
    with open(fi_stoFile, 'r') as f:
        lines = f.readlines()

    alignments = dict()
    sequences = [l for l in lines if not l.startswith('#') and not l.startswith('/') and l.strip()]

    for line in sequences:
        parts = line.split()
        seq =parts[1].upper().replace('.', '-')
        if len(parts) == 2 and seq not in alignments.values():
            name = parts[0]
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
    """Run shell commands"""
    print(f"Running: {command}")
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"❌ Command failed: {result.stderr}")
    return result.returncode
    #process.communicate()

def check_branch_length(diroutput, rna):
    expected_files = [f"{i:02d}" for i in range(1, 11)]
    dirRNA=os.path.join(diroutput, rna)

    delete_files=list()
    for file_name in os.listdir(dirRNA):
        for seed in expected_files:
            if file_name.startswith('RAxML_bestTree') and file_name.endswith(seed):
                tree_path=os.path.join(dirRNA, file_name)
                tree=Phylo.read(tree_path, 'newick')
                for clade in tree.find_clades():
                    if clade.branch_length and clade.branch_length > 1:
                        delete_files.append(file_name)

    if len(delete_files) > 0:
        return rna


def main():
    RAXML_EXECUTE = '/scratch/dx61/vh5686/tmp/RNAPhylo/tools/standard-RAxML/raxmlHPC'
    DIR_WORKING = '/scratch/dx61/vh5686/tmp/RNAPhylo'
    PROJECT='500_full_alignments'
    DIR_SUBS=['data/sto', 'data/fasta', 'data/converted_ss',
              'logs', 'outputs_analysis', 'logs' ]

    PATHS=create_directories(DIR_WORKING, PROJECT, DIR_SUBS)
    log_filename=join(PATHS['logs'], f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log")
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

    # manually parse 500 RNA families
    analysed_nopseudo_rnas = ['RF01922', 'RF00894', 'RF00712', 'RF02421', 'RF01797', 'RF04204', 'RF01316', 'RF03851', 'RF02539',
                     'RF01025', 'RF01700', 'RF01245', 'RF00226', 'RF03632', 'RF01736', 'RF00273', 'RF00110', 'RF02506',
                     'RF00676', 'RF02896', 'RF02823', 'RF03602', 'RF01202', 'RF02772', 'RF00232', 'RF00839', 'RF00324',
                     'RF04066', 'RF01039', 'RF02449', 'RF03565', 'RF02260', 'RF02611', 'RF00817', 'RF01393', 'RF01042',
                     'RF00624', 'RF02830', 'RF04121', 'RF02827', 'RF00325', 'RF02695', 'RF01548', 'RF01939', 'RF00825',
                     'RF00948', 'RF04061', 'RF02834', 'RF02111', 'RF00786', 'RF00950', 'RF02589', 'RF02709', 'RF00873',
                     'RF00358', 'RF04294', 'RF04080', 'RF02424', 'RF01281', 'RF01278', 'RF00819', 'RF00418', 'RF01168',
                     'RF00128', 'RF02534', 'RF00760', 'RF03833', 'RF02887', 'RF01191', 'RF00125', 'RF03354', 'RF00513',
                    'RF04054', 'RF00301', 'RF01173', 'RF01707', 'RF03212', 'RF01300', 'RF00542', 'RF01211',
                     'RF04171', 'RF00852', 'RF01712', 'RF00196', 'RF00639', 'RF00483', 'RF00565', 'RF01643', 'RF00043',
                     'RF02970', 'RF02808', 'RF01156', 'RF02601', 'RF01001', 'RF03508', 'RF00360', 'RF01468', 'RF00293',
                     'RF00744', 'RF00352', 'RF01781', 'RF02500', 'RF03150', 'RF01943', 'RF02911', 'RF02829', 'RF00629',
                     'RF01216', 'RF03532', 'RF03542', 'RF02349', 'RF00282', 'RF00992', 'RF02697', 'RF01864', 'RF01721',
                     'RF02274', 'RF01690', 'RF01458', 'RF02392', 'RF01918', 'RF01590', 'RF02751', 'RF00620', 'RF01041',
                     'RF01259', 'RF00377', 'RF01776', 'RF01522', 'RF01397', 'RF00166', 'RF00227', 'RF00362', 'RF03313',
                     'RF04254', 'RF00265', 'RF02687', 'RF04194', 'RF00065', 'RF00939', 'RF00112', 'RF01348', 'RF01016',
                     'RF02442', 'RF00425', 'RF02372', 'RF01323', 'RF00224', 'RF00775', 'RF02378', 'RF00397', 'RF02235',
                     'RF02980', 'RF02902', 'RF00394', 'RF00200', 'RF00175', 'RF03017', 'RF01836', 'RF00243', 'RF02801',
                     'RF04189', 'RF03714', 'RF02701', 'RF00969', 'RF04154', 'RF00384', 'RF00814', 'RF01379', 'RF01724',
                     'RF00056', 'RF01288', 'RF01471', 'RF03008', 'RF01413', 'RF00477', 'RF00518', 'RF02447', 'RF03414',
                     'RF01023', 'RF02261', 'RF00549', 'RF00319', 'RF00182', 'RF04297', 'RF00725', 'RF00787', 'RF00554',
                     'RF01902', 'RF00671', 'RF00286', 'RF02779', 'RF02238', 'RF03479', 'RF02715', 'RF02577', 'RF00276',
                     'RF04048', 'RF01988', 'RF00600', 'RF03082', 'RF02385', 'RF01269', 'RF02435', 'RF03662', 'RF00247',
                     'RF01230', 'RF00077', 'RF02718', 'RF01221', 'RF03330', 'RF02359', 'RF00531', 'RF02524', 'RF00393',
                     'RF00692', 'RF00578', 'RF04221', 'RF00091', 'RF01265', 'RF03834', 'RF00857', 'RF03060', 'RF03778',
                     'RF01633', 'RF01243', 'RF00170', 'RF04079', 'RF02017', 'RF03296', 'RF01841', 'RF03792', 'RF00440',
                     'RF00317', 'RF01920', 'RF01483', 'RF00145', 'RF03944', 'RF00217', 'RF04117', 'RF00598', 'RF02454',
                     'RF01656', 'RF03442', 'RF01718', 'RF04097', 'RF04275', 'RF00199', 'RF02782', 'RF00892', 'RF03779',
                    'RF02339', 'RF03997', 'RF01497', 'RF03803', 'RF00212', 'RF00813', 'RF03059', 'RF00474', 'RF00071',
                              'RF00204', 'RF00184', 'RF00561', 'RF00374', 'RF01538', 'RF02355', 'RF01383', 'RF01660',
                              'RF03034', 'RF01280', 'RF02753', 'RF01408', 'RF01756', 'RF02694', 'RF01260', 'RF00400',
                              'RF00220', 'RF02833', 'RF03768', 'RF00341', 'RF02051', 'RF00370', 'RF01815', 'RF00433',
                              'RF01395', 'RF03165', 'RF00471', 'RF02070', 'RF00411', 'RF00099', 'RF00436', 'RF01539', 'RF00451', 'RF04078', 'RF04001', 'RF00788', 'RF01720', 'RF01638', 'RF01769', 'RF03103', 'RF01990', 'RF04268', 'RF01547', 'RF00574', 'RF02068', 'RF01290', 'RF02939', 'RF02846', 'RF03048', 'RF00809', 'RF01859', 'RF04270', 'RF01495', 'RF01198', 'RF02613', 'RF01200', 'RF00300', 'RF02370', 'RF00952', 'RF01167', 'RF02693', 'RF00773', 'RF03042', 'RF01178', 'RF01861', 'RF00263', 'RF00803', 'RF02866', 'RF02946', 'RF04236', 'RF00615', 'RF03600', 'RF03350', 'RF03161', 'RF03326', 'RF04175', 'RF03047', 'RF00280', 'RF00283', 'RF00187', 'RF00777', 'RF02730', 'RF02271', 'RF01830', 'RF01315', 'RF04164', 'RF00081', 'RF02375', 'RF01525', 'RF01004', 'RF04134', 'RF00919', 'RF01655', 'RF00630', 'RF00897', 'RF01925', 'RF00335', 'RF04069', 'RF03021', 'RF01384', 'RF02450', 'RF00784', 'RF02981', 'RF00792', 'RF03981', 'RF03480', 'RF00461', 'RF01647', 'RF01898', 'RF01423', 'RF00972', 'RF01272', 'RF00183', 'RF03242', 'RF00189', 'RF01264', 'RF01487', 'RF00891', 'RF00510', 'RF03058', 'RF00260', 'RF03218', 'RF03093', 'RF02520', 'RF01189', 'RF00399', 'RF03390', 'RF02757', 'RF03199', 'RF03695', 'RF00291', 'RF03255', 'RF00560', 'RF00136', 'RF01250', 'RF01901', 'RF00489', 'RF02412', 'RF03280', 'RF02405', 'RF01653', 'RF02263', 'RF00225', 'RF03053', 'RF02464', 'RF03073', 'RF01268', 'RF00248', 'RF02835', 'RF01232', 'RF01684', 'RF03525', 'RF02950', 'RF00179', 'RF03456', 'RF00996', 'RF01218', 'RF03860', 'RF02374', 'RF00195', 'RF00828', 'RF00714', 'RF01061', 'RF01804', 'RF00410', 'RF00523', 'RF01897', 'RF02269', 'RF00409', 'RF01206', 'RF02971', 'RF02228', 'RF01706', 'RF00697', 'RF01231', 'RF03849', 'RF01419', 'RF01782', 'RF00106', 'RF00641', 'RF00312', 'RF01599', 'RF00997', 'RF00552', 'RF00741', 'RF02009', 'RF02744', 'RF01255', 'RF00014', 'RF02162', 'RF00837', 'RF03101', 'RF01192', 'RF01496', 'RF04072', 'RF01275', 'RF02540', 'RF00417', 'RF02783', 'RF02822', 'RF01999', 'RF01523', 'RF00464', 'RF02377', 'RF00019', 'RF01592', 'RF00749', 'RF00878', 'RF00355', 'RF00084', 'RF03814', 'RF03006', 'RF02923', 'RF01183', 'RF04063', 'RF01819', 'RF00171', 'RF00264', 'RF02225', 'RF01461', 'RF01457', 'RF01307', 'RF03198', 'RF00022', 'RF02411', 'RF00177', 'RF01821', 'RF00194', 'RF02778', 'RF02408', 'RF01742', 'RF01702', 'RF01287', 'RF01912', 'RF02912', 'RF02978', 'RF00240', 'RF03012', 'RF03254', 'RF00928', 'RF01729', 'RF00188', 'RF00753', 'RF01070', 'RF03378', 'RF03402', 'RF01263', 'RF03141', 'RF01900', 'RF00340', 'RF00419', 'RF02639', 'RF00929', 'RF04084', 'RF01464', 'RF00475', 'RF03127', 'RF00737', 'RF01475', 'RF00209', 'RF02360', 'RF01796', 'RF00252', 'RF01317', 'RF00037', 'RF02523', 'RF04181', 'RF00536', 'RF00117', 'RF00062', 'RF01637', 'RF04218', 'RF00991', 'RF00689', 'RF00569', 'RF00392', 'RF03767', 'RF01903', 'RF01820', 'RF03125', 'RF02466', 'RF00608', 'RF02532', 'RF01915', 'RF02831', 'RF01456', 'RF00303', 'RF00681', 'RF00135', 'RF04244', 'RF01251', 'RF02144', 'RF01936', 'RF00599', 'RF00901', 'RF01793', 'RF02429', 'RF01056', 'RF00734', 'RF02596', 'RF00735', 'RF01018', 'RF00215', 'RF00935', 'RF00951', 'RF00158', 'RF01784', 'RF03349', 'RF00915', 'RF00049', 'RF04053', 'RF02344', 'RF00957', 'RF00155', 'RF01234', 'RF01320', 'RF00092', 'RF00137', 'RF04064', 'RF01225', 'RF01771', 'RF00577', 'RF00424', 'RF01394', 'RF01835', 'RF02748', 'RF03193', 'RF00460', 'RF00520', 'RF02496', 'RF01069', 'RF04126', 'RF00650', 'RF03102', 'RF00855', 'RF00408', 'RF04267', 'RF00415', 'RF03883', 'RF04179', 'RF01773', 'RF00610', 'RF00674', 'RF01421', 'RF00296', 'RF00655', 'RF02860', 'RF01285', 'RF00699', 'RF00646', 'RF01010', 'RF01455', 'RF01635', 'RF03939', 'RF00766', 'RF04067', 'RF00118', 'RF02975', 'RF00402', 'RF02227', 'RF03842', 'RF01165', 'RF01557', 'RF02644', 'RF02457', 'RF02517', 'RF00765', 'RF02842', 'RF00533', 'RF00151', 'RF02367', 'RF01390', 'RF02917', 'RF04271', 'RF01417', 'RF02689', 'RF00586', 'RF00921', 'RF00539', 'RF02848', 'RF03947', 'RF01699', 'RF00987', 'RF00484', 'RF02459', 'RF02927', 'RF02095', 'RF01274', 'RF01865', 'RF03843', 'RF01239', 'RF01620', 'RF02093', 'RF01411', 'RF01038', 'RF01588', 'RF00085', 'RF02354', 'RF02713', 'RF01286', 'RF00840', 'RF00582', 'RF01392', 'RF01824', 'RF01743', 'RF01858', 'RF00795', 'RF01170', 'RF00771', 'RF01685', 'RF00421', 'RF00722', 'RF01794', 'RF01625', 'RF00810', 'RF00186', 'RF01426', 'RF03001', 'RF00774', 'RF03608', 'RF01549', 'RF02770', 'RF01550', 'RF01233', 'RF00369', 'RF02234', 'RF02636', 'RF01673', 'RF01544', 'RF03983', 'RF02865', 'RF01696', 'RF00796', 'RF02747', 'RF02790', 'RF00687', 'RF00266', 'RF00144', 'RF01219', 'RF00193', 'RF00133', 'RF01610', 'RF04272', 'RF00082', 'RF00736', 'RF00229', 'RF03091', 'RF00304', 'RF00719', 'RF04151', 'RF02353', 'RF00568', 'RF01491', 'RF03249', 'RF00359', 'RF01158', 'RF03495', 'RF00567', 'RF02226', 'RF04285', 'RF00305', 'RF03297', 'RF03066', 'RF01630', 'RF03052', 'RF01065', 'RF01273', 'RF04249', 'RF01866', 'RF01235', 'RF00778', 'RF00258', 'RF01253', 'RF02422', 'RF00835', 'RF02688', 'RF02876', 'RF01294', 'RF04191', 'RF02547', 'RF03115', 'RF03693', 'RF01517', 'RF01011', 'RF02410', 'RF03395', 'RF00989', 'RF00998', 'RF02893', 'RF00936', 'RF04274', 'RF01792', 'RF02345', 'RF03035', 'RF02966', 'RF00544', 'RF02993', 'RF00448', 'RF02941', 'RF00824', 'RF01552', 'RF01814', 'RF03517', 'RF02974', 'RF00383', 'RF00172', 'RF00441', 'RF02346', 'RF02569', 'RF00096', 'RF00927', 'RF01839', 'RF01238', 'RF03794', 'RF00101', 'RF01732', 'RF00428', 'RF03661', 'RF00073', 'RF01046', 'RF00078', 'RF01790', 'RF00591', 'RF03329', 'RF02433', 'RF00452', 'RF00750', 'RF01318', 'RF00918', 'RF02462', 'RF00447', 'RF03372', 'RF01203', 'RF00427', 'RF02240', 'RF01477', 'RF01791', 'RF02867', 'RF01014', 'RF01693', 'RF01862', 'RF00113', 'RF01811', 'RF00994', 'RF01945', 'RF02384', 'RF02922', 'RF02677', 'RF04095', 'RF03258', 'RF00214', 'RF01402', 'RF02750', 'RF02797', 'RF02755', 'RF00222', 'RF02942', 'RF02954', 'RF02983', 'RF01401', 'RF03436', 'RF02844', 'RF01694', 'RF02581', 'RF00115', 'RF00038', 'RF00089', 'RF00467', 'RF02342', 'RF02362', 'RF00572', 'RF01595', 'RF04219', 'RF00707', 'RF02074', 'RF03445', 'RF01719', 'RF03089', 'RF00758', 'RF03109', 'RF03134', 'RF00553', 'RF01914', 'RF00294', 'RF03166', 'RF00435', 'RF00801', 'RF03039', 'RF01540', 'RF02551', 'RF00845', 'RF00102', 'RF02273', 'RF01415', 'RF00254', 'RF02995', 'RF00088', 'RF02652', 'RF02088', 'RF02544', 'RF00438', 'RF01313', 'RF00495', 'RF03784', 'RF00153', 'RF00181', 'RF01326', 'RF01462', 'RF02952', 'RF03098', 'RF01853', 'RF00121', 'RF00543', 'RF01197', 'RF02061', 'RF02802', 'RF00016', 'RF01982', 'RF01277', 'RF04077', 'RF00757', 'RF01391', 'RF02222', 'RF03379', 'RF01196', 'RF04288', 'RF02724', 'RF02692', 'RF00403', 'RF00558', 'RF00231', 'RF01213', 'RF02961', 'RF00666', 'RF00478', 'RF00961', 'RF00797', 'RF02897', 'RF02803', 'RF00731', 'RF00211', 'RF04124', 'RF00836', 'RF03029', 'RF00018', 'RF00677', 'RF02052', 'RF01896', 'RF00386', 'RF00236', 'RF01723', 'RF04292', 'RF01412', 'RF00332', 'RF02940', 'RF02465', 'RF04163', 'RF00287', 'RF02401', 'RF01825', 'RF02913', 'RF01774', 'RF02714', 'RF01802', 'RF00289', 'RF02232', 'RF02535', 'RF03543', 'RF03741', 'RF03043', 'RF04296', 'RF02885', 'RF04276', 'RF01034', 'RF00139', 'RF01024', 'RF02768', 'RF02510', 'RF00190', 'RF03325', 'RF02828', 'RF02528', 'RF01043', 'RF02675', 'RF01490', 'RF01813', 'RF01494', 'RF00501', 'RF02512', 'RF00690', 'RF01666', 'RF02945', 'RF03692', 'RF01258', 'RF02962', 'RF01407', 'RF03369', 'RF04122', 'RF01919', 'RF04076', 'RF03956', 'RF03070', 'RF00439', 'RF02722', 'RF00482', 'RF00353', 'RF01387', 'RF00203', 'RF00958', 'RF02423', 'RF01215', 'RF01857', 'RF00330', 'RF00691', 'RF02765', 'RF02605', 'RF00649', 'RF00149', 'RF02832', 'RF01787', 'RF04224', 'RF03120', 'RF00551', 'RF01220', 'RF02884', 'RF01227', 'RF00859', 'RF02419', 'RF01581', 'RF02548', 'RF00481', 'RF01668', 'RF04045', 'RF00388', 'RF00670', 'RF02265', 'RF04248', 'RF02900', 'RF03601', 'RF00636', 'RF00614', 'RF02643', 'RF03186', 'RF00965', 'RF00781', 'RF00696', 'RF00806', 'RF01585', 'RF00841', 'RF02376', 'RF00257', 'RF00812', 'RF01867', 'RF02826', 'RF01803', 'RF00116', 'RF00547', 'RF01344', 'RF02895', 'RF01284', 'RF00311', 'RF01308', 'RF01470', 'RF00526', 'RF00584', 'RF02084', 'RF00580', 'RF03718', 'RF00827', 'RF00678', 'RF01398', 'RF00853', 'RF02727', 'RF03155', 'RF00083', 'RF00740', 'RF01262', 'RF01474', 'RF01224', 'RF02804', 'RF00357', 'RF01422', 'RF01321', 'RF00664', 'RF00202', 'RF02092', 'RF01687', 'RF03038', 'RF02869', 'RF01996', 'RF03310', 'RF03533', 'RF03083', 'RF01324', 'RF01335', 'RF01396', 'RF00322', 'RF00072', 'RF03114', 'RF04269', 'RF00769', 'RF04227', 'RF03036', 'RF02343', 'RF02756', 'RF00157', 'RF02444', 'RF00680', 'RF01012', 'RF02673', 'RF01240', 'RF01762', 'RF01400', 'RF00180', 'RF02053', 'RF00437', 'RF02583', 'RF00973', 'RF00856', 'RF00281', 'RF00966', 'RF00721', 'RF01600', 'RF02079', 'RF00762', 'RF04088', 'RF02396', 'RF02796', 'RF00040', 'RF00336', 'RF01938', 'RF00570', 'RF02268', 'RF02746', 'RF00159', 'RF03962', 'RF04263', 'RF00161', 'RF02774', 'RF04234', 'RF01459', 'RF01558', 'RF04082', 'RF04216', 'RF00742', 'RF02798', 'RF00070', 'RF01430', 'RF01554', 'RF04261', 'RF00995', 'RF02371', 'RF00843', 'RF01758', 'RF00613', 'RF00842', 'RF00459', 'RF02031', 'RF00480', 'RF01662', 'RF00150', 'RF01937', 'RF00223', 'RF02738', 'RF02503', 'RF01405', 'RF00887', 'RF03051', 'RF04073', 'RF00705', 'RF01508', 'RF03099', 'RF01044', 'RF00816', 'RF02075', 'RF00623', 'RF00138', 'RF00906', 'RF01309', 'RF03050', 'RF01587', 'RF02030', 'RF02413', 'RF01267', 'RF02728', 'RF04027', 'RF02434', 'RF00338', 'RF00107', 'RF00848', 'RF00888', 'RF04300', 'RF01314', 'RF00728', 'RF03552', 'RF00396', 'RF03149', 'RF00490', 'RF01488', 'RF02773', 'RF01045', 'RF03027', 'RF00329', 'RF02999', 'RF00829', 'RF00693', 'RF04057', 'RF03357', 'RF00126', 'RF02725', 'RF00061', 'RF02243', 'RF00444', 'RF02415', 'RF00218', 'RF00815', 'RF01651', 'RF02781', 'RF01816', 'RF01283', 'RF00954', 'RF00036', 'RF01537', 'RF04295', 'RF02000', 'RF03049', 'RF00672', 'RF02463', 'RF00404', 'RF00349', 'RF02690', 'RF00772', 'RF02752', 'RF02754', 'RF00339', 'RF00718', 'RF00295', 'RF00238', 'RF00198', 'RF00242', 'RF00401', 'RF00527', 'RF01015', 'RF02807', 'RF00423', 'RF02076', 'RF00745', 'RF02682', 'RF01292', 'RF01208', 'RF01291', 'RF01410', 'RF01179', 'RF04150', 'RF00943', 'RF03011', 'RF00550', 'RF00124', 'RF02696', 'RF00785', 'RF03024', 'RF01654', 'RF01927', 'RF00345', 'RF03085', 'RF01337', 'RF00606', 'RF01028', 'RF00528', 'RF02425', 'RF00964', 'RF02928', 'RF01241', 'RF00616', 'RF02448', 'RF00917', 'RF00343', 'RF01226', 'RF02819', 'RF04052', 'RF00279', 'RF00486', 'RF00847', 'RF00375', 'RF00376', 'RF01524', 'RF00962', 'RF01827', 'RF01489', 'RF00492', 'RF00609', 'RF03469', 'RF03995', 'RF00271', 'RF00761', 'RF03094', 'RF02758', 'RF00256', 'RF00955', 'RF00328', 'RF02989', 'RF00798', 'RF00154', 'RF03022', 'RF02431', 'RF00503', 'RF00356', 'RF01157', 'RF03321', 'RF00947', 'RF02976', 'RF02445', 'RF01910', 'RF01923', 'RF00462', 'RF00228', 'RF00865', 'RF04286', 'RF01256', 'RF02509', 'RF03620', 'RF00119', 'RF00453', 'RF01418', 'RF00830', 'RF00879', 'RF03108', 'RF03040', 'RF02402', 'RF03251', 'RF01823', 'RF00316', 'RF00105', 'RF00414', 'RF02883', 'RF03669', 'RF02264', 'RF00302', 'RF00579', 'RF02968', 'RF00512', 'RF04090', 'RF00770', 'RF00063', 'RF01382', 'RF00541', 'RF02691', 'RF01767', 'RF01671', 'RF03033', 'RF00079', 'RF01703', 'RF00621', 'RF02366', 'RF03061', 'RF00463', 'RF01244', 'RF00292', 'RF00656', 'RF00581', 'RF00564', 'RF03687', 'RF00644', 'RF03153', 'RF02469', 'RF04074', 'RF01697', 'RF00272', 'RF03444', 'RF02780', 'RF03731', 'RF00862', 'RF00912', 'RF00337', 'RF01692', 'RF01518', 'RF01271', 'RF02348', 'RF03214', 'RF00468', 'RF01481', 'RF00406', 'RF00849', 'RF02775', 'RF03464', 'RF01270', 'RF00164', 'RF00250', 'RF00861', 'RF03374', 'RF00546', 'RF02386', 'RF00350', 'RF00076', 'RF00420', 'RF01845', 'RF04044', 'RF03668', 'RF04100', 'RF04165', 'RF02855', 'RF00416', 'RF02850', 'RF02931', 'RF00704', 'RF02278', 'RF00885', 'RF00449', 'RF00710', 'RF00642', 'RF01473', 'RF04250', 'RF04225', 'RF03007', 'RF01989', 'RF02580', 'RF02003', 'RF01478', 'RF04000', 'RF02438', 'RF01229', 'RF00479', 'RF00364', 'RF01130', 'RF02552', 'RF00984', 'RF01452', 'RF00368', 'RF02874', 'RF03840', 'RF01480', 'RF04009', 'RF01629', 'RF00493', 'RF03544', 'RF00454', 'RF01535', 'RF00422', 'RF04075', 'RF00525', 'RF01406', 'RF00748', 'RF02100', 'RF02369', 'RF01033', 'RF01832', 'RF00659', 'RF02936', 'RF04226', 'RF00147', 'RF01753', 'RF02716', 'RF01716', 'RF02767', 'RF03023', 'RF00981', 'RF01399', 'RF00872', 'RF04290', 'RF01193', 'RF01186', 'RF00348', 'RF02276', 'RF03025', 'RF00046', 'RF00870', 'RF00267', 'RF03311', 'RF03691', 'RF00764', 'RF04182', 'RF00039', 'RF00385', 'RF03163', 'RF00288', 'RF00611', 'RF04193', 'RF00920', 'RF00445', 'RF01913', 'RF02451', 'RF00069', 'RF02432', 'RF02446', 'RF00686', 'RF03055', 'RF00315', 'RF02531', 'RF02674', 'RF00253', 'RF01325', 'RF03283', 'RF00665', 'RF02891', 'RF01663', 'RF01169', 'RF00275', 'RF02545', 'RF01002', 'RF00877', 'RF00178', 'RF01059', 'RF02504', 'RF02806', 'RF00431', 'RF02083', 'RF00593', 'RF00457', 'RF02556', 'RF02037', 'RF01810', 'RF00985', 'RF01175', 'RF04228', 'RF00432', 'RF00455', 'RF00612', 'RF03111', 'RF01319', 'RF03875', 'RF00723', 'RF00658', 'RF00244', 'RF02877', 'RF00673', 'RF01176', 'RF02065', 'RF03122', 'RF00494', 'RF00207', 'RF00060', 'RF04169', 'RF01546', 'RF02400', 'RF02340', 'RF03838', 'RF02409', 'RF00660', 'RF01609', 'RF00033', 'RF00068', 'RF02530', 'RF00487', 'RF01052', 'RF00529', 'RF00831', 'RF01385', 'RF02067', 'RF03534', 'RF03229', 'RF00709', 'RF01766', 'RF01755', 'RF02077', 'RF00566', 'RF03279', 'RF01433', 'RF04029', 'RF01674', 'RF01013', 'RF03301', 'RF04087', 'RF04220', 'RF02985', 'RF04284', 'RF01162', 'RF03056', 'RF01472', 'RF01214', 'RF00602', 'RF02764', 'RF03175', 'RF00434', 'RF02341', 'RF04282', 'RF02507', 'RF03031', 'RF04156', 'RF02784', 'RF00191', 'RF02953', 'RF03404', 'RF03706', 'RF01924', 'RF00326', 'RF02700', 'RF00603', 'RF00321', 'RF00702', 'RF01626', 'RF00429', 'RF00284', 'RF03110', 'RF01428', 'RF00573', 'RF04280', 'RF00540', 'RF00782', 'RF00562', 'RF02399', 'RF02495', 'RF04085', 'RF01675', 'RF01691', 'RF01228', 'RF00807', 'RF03896', 'RF01531', 'RF03355', 'RF00818', 'RF00911', 'RF00530', 'RF00142', 'RF01027', 'RF00763', 'RF00382', 'RF00698', 'RF01019', 'RF00937', 'RF04266', 'RF00537', 'RF03010', 'RF00854', 'RF01757', 'RF04089', 'RF02992', 'RF01555', 'RF00122', 'RF01007', 'RF00201', 'RF03777', 'RF00730', 'RF00890', 'RF01608', 'RF01006', 'RF00846', 'RF01000', 'RF00733', 'RF02610', 'RF00708', 'RF02081', 'RF00205', 'RF02898', 'RF02453', 'RF03092', 'RF00413', 'RF02233', 'RF02064', 'RF02899', 'RF02541', 'RF02262', 'RF02937', 'RF02723', 'RF03271', 'RF01602', 'RF01808', 'RF01431', 'RF00058', 'RF04205', 'RF00886', 'RF01381', 'RF04229', 'RF00156', 'RF01254', 'RF01899', 'RF01299', 'RF03881', 'RF02595', 'RF04183', 'RF02277', 'RF00443', 'RF00715', 'RF02364', 'RF00900', 'RF03528', 'RF00090', 'RF04291', 'RF01921', 'RF01322', 'RF02982', 'RF03097', 'RF01182', 'RF00627', 'RF00132', 'RF03368', 'RF03531', 'RF01917', 'RF01941', 'RF00470', 'RF01627', 'RF02230', 'RF02786', 'RF02270', 'RF02379', 'RF00703', 'RF01261', 'RF02060', 'RF00097', 'RF02231', 'RF02655', 'RF00405', 'RF02420', 'RF03105', 'RF00021', 'RF00160', 'RF03639', 'RF00511', 'RF00371', 'RF00563', 'RF00910', 'RF04293', 'RF02566', 'RF00407', 'RF00711', 'RF01164', 'RF01959', 'RF02915', 'RF02437', 'RF03002', 'RF03408', 'RF02055', 'RF03470', 'RF01032', 'RF02979', 'RF04214', 'RF01332', 'RF00755', 'RF00206', 'RF02929', 'RF02436', 'RF04217', 'RF03184', 'RF01476', 'RF02266', 'RF00942', 'RF03653', 'RF00756', 'RF01512', 'RF01484', 'RF00426', 'RF02732', 'RF00933', 'RF04056', 'RF01493', 'RF01676', 'RF00389', 'RF01209', 'RF00430', 'RF00601', 'RF00986', 'RF04136', 'RF03538', 'RF02502', 'RF03684', 'RF01926', 'RF03292', 'RF00333', 'RF04158']

    # This section is aimed to download only required files
    dir_sto = PATHS['data/sto']
    dir_output = PATHS['outputs_analysis']
    failed_downloads = []
    for rf in analysed_nopseudo_rnas:
        sto_file=os.path.join(dir_sto, f'{rf}.sto')

        # 🟢 Step 1: Skip download if file exists
        if os.path.isfile(sto_file):
            print(f"✅ {rf}.sto already exists. Skipping download.")
            continue

        # 🟠 Step 2: File does not exist, check URL availability
        url = f"https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/full_alignments/{rf}.sto"
        try:
            response = requests.head(url, timeout=10)
            if response.status_code != 200:
                print(f"⚠️ {rf}.sto does not exist in the Rfam database. Skipping.")
                failed_downloads.append(rf)
                continue
        except requests.RequestException as e:
            print(f"❌ Error checking {url}: {e}")
            failed_downloads.append(rf)
            continue

        # 🔵 Step 3: Download the file
        print(f"⬇️ Downloading {rf}.sto...")
        result = run_command(f"wget {url} -P {dir_sto}")

        # 🛑 Step 4: Check if download was successful
        if result != 0 or not os.path.isfile(sto_file):
            print(f"❌ Failed to download {rf}.sto.")
            failed_downloads.append(rf)
            continue

        print(f"✅ Successfully downloaded {rf}.sto")

        # 🟢 Step 5: Process the file if download succeeded
        nodup_fasta, paired_pseudo_ss, unpaired_pseudo_ss = parseFastaAndSS(
            sto_file,
            PATHS['data/fasta'],
            PATHS['data/converted_ss'],
            rf
        )

        # 🔵 Step 6: Submit jobs if valid files were generated
        if nodup_fasta and paired_pseudo_ss and unpaired_pseudo_ss:
            raxml_prefix = os.path.join(dir_output, 'raxml', rf)
            raxml_unpaired_pseudo_prefix = os.path.join(dir_output, 'raxml_iPseu', rf)

            raxml_command = (
                f"qsub -V -N raxml_{rf} -o {PATHS['logs']} -e {PATHS['logs']} "
                f"-l ncpus=12 -l mem=48gb -l walltime=48:00:00 -l wd -- "
                f"bash bashFiles/raxml.sh {rf} {nodup_fasta} {raxml_prefix} {RAXML_EXECUTE}"
            )
            raxml_unpaired_pseudo_command = (
                f"qsub -V -N raxmlP_{rf} -o {PATHS['logs']} -e {PATHS['logs']} "
                f"-l ncpus=12 -l mem=48gb -l walltime=48:00:00 -l wd -- "
                f"bash bashFiles/raxmlP.sh {rf} {nodup_fasta} {unpaired_pseudo_ss} {raxml_unpaired_pseudo_prefix} {RAXML_EXECUTE}"
            )

            run_command(raxml_command)
            run_command(raxml_unpaired_pseudo_command)

        else:
            print(f"⚠️ {rf} does not qualify for the analysis.")

    # 🔴 Step 7: Print a report of failed downloads
    if failed_downloads:
        print("\n🚨 The following files failed to download or were missing:")
        for rf in failed_downloads:
            print(f" - {rf}")
    else:
        print("✅ All files processed successfully.")
"""
    if os.path.isfile(sto_file):
        nodup_fasta, paired_pseudo_ss, unpaired_pseudo_ss = parseFastaAndSS(sto_file,
                                                                            PATHS['data/fasta'],
                                                                            PATHS['data/converted_ss'],
                                                                            rf)

        if nodup_fasta and paired_pseudo_ss and unpaired_pseudo_ss:
            raxml_prefix = join(dir_output, 'raxml', rf)
            #raxml_paired_pseudo_prefix = join(DIR_OUTPUTS, '500_full_alignments', 'raxml_wPseu', rf)
            raxml_unpaired_pseudo_prefix = join(dir_output, 'raxml_iPseu', rf)

            raxml_command = f"qsub -V -N raxml_{rf} -o {PATHS['logs']} -e {PATHS['logs']} -l ncpus=12 -l mem=48gb -l walltime=48:00:00 -l wd -- bash bashFiles/raxml.sh {rf} {nodup_fasta} {raxml_prefix} {RAXML_EXECUTE}"
            raxml_unpaired_pseudo_command = f"qsub -V -N raxmlP_{rf} -o {PATHS['logs']} -e {PATHS['logs']} -l ncpus=12 -l mem=48gb -l walltime=48:00:00 -l wd -- bash bashFiles/raxmlP.sh {rf} {nodup_fasta} {unpaired_pseudo_ss} {raxml_unpaired_pseudo_prefix} {RAXML_EXECUTE}"
            #raxml_paired_pseudo_command = f"bash bashFiles/raxmlP.sh {rf} {nodup_fasta} {paired_pseudo_ss} {raxml_paired_pseudo_prefix} {RAXML_EXECUTE}"

            run_command(raxml_command)
            run_command(raxml_unpaired_pseudo_command)
            #run_command(raxml_paired_pseudo_command)

        else:
            print(f'{rf} does not qualify for the analysis.')
"""

if __name__ == "__main__":
    main()
