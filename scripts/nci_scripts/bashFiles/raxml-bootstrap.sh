#!/bin/bash

RF=$1
FASTA_FILE=$2
PREFIX=$3
RAXML=$4

mkdir -p $PREFIX

for i in {1..10}; do
    formatted_i=$(printf "%02d" $i)
    $RAXML -f a -x $i -p $i -# 100 -s $FASTA_FILE -n ${RF}.${formatted_i} -m GTRCAT -w ${PREFIX}
done

# raxml -f a -x 12345 -p 12345 -# 100 -s fasta/RF00190.nodup.fa -m GTRCAT -n TEST
# -x set the random seed for the boostrap replicates
# -p set the random seed for the parsimony starting tree
