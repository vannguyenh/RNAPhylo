#!/bin/bash

RF=$1
FASTA_FILE=$2
SS_FILE=$3
PREFIX=$4
RAXML=$5

mkdir -p $PREFIX
#mkdir -p $PREFIX_D

for i in {1..10}; do
    formatted_i=$(printf "%02d" $i)
    $RAXML -f a -x $i -p $i -# 100 -s $FASTA_FILE -n ${RF}.${formatted_i} -m GTRGAMMA -p $i -S $SS_FILE -w ${PREFIX}
done
# raxml -f a -x 12345 -p 12345 -#100 -s fasta/RF00190.nodup.fa -S converted_ss/RF00190.pseudo.unpaired.ss -m GTRCAT -n TEST.SS
