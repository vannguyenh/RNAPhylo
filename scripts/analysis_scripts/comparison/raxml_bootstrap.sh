#!/bin/bash
set -euo pipefail

RAXML='/Users/u7875558/tools/standard-RAxML-master/raxmlHPC'

RNA="$1"         # e.g. RF00740
SEED_DNA="$2"    # e.g. 10 
SEED_RNA="$3"    # e.g. 6
RNA_MODEL="$4"	 # e.g. S16

format_seedDNA=$(printf "%02d" "$SEED_DNA")
format_seedRNA=$(printf "%02d" "$SEED_RNA")

input_path='/Users/u7875558/RNAPhylo/seedAlignment_AllModels/inputs'
inferred_path='/Users/u7875558/RNAPhylo/seedAlignment_AllModels/outputs/inferred_trees'
comparison_path='/Users/u7875558/RNAPhylo/seedAlignment_AllModels/comparison'   # <-- set this
outdir="${comparison_path}/${RNA}"
mkdir -p "$outdir"


# ---------------- DNA BOOTSTRAP ----------------
# Generate 100 bootstrap replicates
"$RAXML" -m GTRGAMMA \
    -s "${input_path}/fasta_files/${RNA}.nodup.fa" \
    -p "${SEED_DNA}" -b "${SEED_DNA}" -# 100 \
    -w "$outdir" \
    -n "${RNA}.dna.${format_seedDNA}.bstr"

"$RAXML" -f b -m GTRGAMMA \
    -s "${input_path}/fasta_files/${RNA}.nodup.fa" \
    -t "${inferred_path}/DNA/${RNA}/RAxML_bestTree.${RNA}.${format_seedDNA}" \
    -z "${outdir}/RAxML_bootstrap.${RNA}.dna.${format_seedDNA}.bstr" \
    -w "$outdir" \
    -n ${RNA}_DNA

# ---------------- RNA (RNA_MODEL) BOOTSTRAP ----------------
# Generate 100 bootstrap replicates under RNA_MODEL + structure
"$RAXML" -m GTRGAMMA -A "${RNA_MODEL}" \
    -s "${input_path}/fasta_files/${RNA}.nodup.fa" \
    -S "${input_path}/ss_files/${RNA}.dots.ss" \
    -p "${SEED_RNA}" -b "${SEED_RNA}" -# 100 \
    -w "$outdir" \
    -n "${RNA}.rna.${RNA_MODEL}.${format_seedRNA}.bstr"

"$RAXML" -f b -m GTRGAMMA -A "${RNA_MODEL}" \
    -s "${input_path}/fasta_files/${RNA}.nodup.fa" \
    -S "${input_path}/ss_files/${RNA}.dots.ss" \
    -t "${inferred_path}/${RNA_MODEL}/${RNA}/RAxML_bestTree.${RNA}.${format_seedRNA}" \
    -z "${outdir}/RAxML_bootstrap.${RNA}.rna.${RNA_MODEL}.${format_seedRNA}.bstr" \
    -w "$outdir" \
    -n "${RNA}_${RNA_MODEL}"
