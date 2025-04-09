#!/bin/bash

RAXML='/Users/u7875558/Documents/PhD/tools/standard-RAxML-master/raxmlHPC'
PATH_CONSEL='/Users/u7875558/Documents/PhD/tools/consel/bin/'

fasta_file=$1
combineTree=$2
persite_suffix=$3
ssRNA=$4
persite_path=$5
output_persite=$6
prefix_consel=$7
model=$8

${RAXML} -f G -s ${fasta_file} -m GTRGAMMA -z ${combineTree} -n ${persite_suffix} -S ${ssRNA} -w ${persite_path} -A ${model}
# raxml -f G -s /Users/u7875558/Documents/PhD/Rfam/RfamPhylo/analysis/outputs/results/tmp/inputs/fasta_files/RF00177.nodup.fa 
# -m GTRGAMMA 
# -z ../RF00177.raxml 
# -n RF00177.raxml_inSS.sitelh 
# -S /Users/u7875558/Documents/PhD/Rfam/RfamPhylo/analysis/outputs/results/tmp/inputs/ss_files/RF00177.dots.ss 
# -w /Users/u7875558/Library/CloudStorage/OneDrive-AustralianNationalUniversity/Documents/PhD/Rfam/RfamPhylo/analysis/outputs/tmp_test/test_raxml_persitelh/RF00177/

### CONSEL program
# CONSEL does not work with parsing the directory to save the file -> have to direct the folder for saving - WRONG
# CONSEL works with directory but the name has to be one word or word1_word2 (not word1.word2)

## TODO: prepare a coherent name to make the pipeline reproducible
${PATH_CONSEL}/makermt --puzzle ${output_persite} ${prefix_consel}
${PATH_CONSEL}/consel ${prefix_consel}

# /Users/u7875558/Documents/PhD/tools/consel/bin/makermt 
# --puzzle /Users/u7875558/Library/CloudStorage/OneDrive-AustralianNationalUniversity/Documents/PhD/Rfam/RfamPhylo/analysis/outputs/tmp_test/test_raxml_persitelh/RF02052/persitelh/RAxML_perSiteLLs.RF02052.raxml_inSS.sitelh 
# /Users/u7875558/Library/CloudStorage/OneDrive-AustralianNationalUniversity/Documents/PhD/Rfam/RfamPhylo/analysis/outputs/tmp_test/test_raxml_persitelh/RF02052/consel/persite_test

# /Users/u7875558/Documents/PhD/tools/consel/bin/consel /Users/u7875558/Library/CloudStorage/OneDrive-AustralianNationalUniversity/Documents/PhD/Rfam/RfamPhylo/analysis/outputs/tmp_test/test_raxml_persitelh/RF02052/consel/persite_test

#${PATH_CONSEL}/catpv
