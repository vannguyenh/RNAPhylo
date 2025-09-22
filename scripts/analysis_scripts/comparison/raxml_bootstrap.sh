#!/bin/bash

RAXML='/Users/u7875558/tools/standard-RAxML-master/raxmlHPC'

RNA=$1
SEED_DNA=$2
SEED_RNA=$3

format_seedDNA=$(printf "%02d" $SEED_DNA)
format_seedRNA=$(printf "%02d" $SEED_RNA)

input_path='/Users/u7875558/RNAPhylo/seedAlignment_AllModels/inputs'
inferred_path='/Users/u7875558/RNAPhylo/seedAlignment_AllModels/outputs/inferred_trees'
comparison_path=''
# bootstrap for DNA tree
${RAXML} -m GTRGAMMA -s ${input_path}/fasta_files/${RNA}.nodup.fa -p ${SEED_DNA} -b ${SEED_DNA} -# 100 -n ${comparison_path}/${RNA}/${RNA}.dna.${format_seedDNA}.bstr
${RAXML} -f b -m GTRGAMMA -s ${input_path}/fasta_files/${RNA}.nodup.fa -t ${inferred_path}/DNA/${RNA}/RAxML_bestTree.${RNA}.${format_seedDNA} -z ${comparison_path}/${RNA}/RAxML_bootstrap.${RNA}.dna.${format_seedDNA}.bstr -n ${RNA}_DNA

# bootstrap for RNA tree
${RAXML} -m GTRGAMMA -A S6A -s ${input_path}/fasta_files/${RNA}.nodup.fa -S ${input_path}/ss_files/${RNA}.dots.ss -p ${SEED_RNA} -S -b ${SEED_RNA} -# 100 -n ${comparison_path}/${RNA}/${RNA}.rna.${format_seedRNA}.bstr
${RAXML} -f b -m GTRGAMMA -A S6A -s ${input_path}/fasta_files/${RNA}.nodup.fa -S ${input_path}/ss_files/${RNA}.dots.ss -t ${inferred_path}/S6A/${RNA}/RAxML_bestTree.${RNA}.${format_seedRNA} -z ${comparison_path}/${RNA}/RAxML_bootstrap.${RNA}.rna.${format_seedRNA}.bstr -n ${RNA}_RNA