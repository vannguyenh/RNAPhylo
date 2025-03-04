#!/bin/bash

RF=$1
FASTA_FILE=$2
PREFIX=$3
MODEL=$4
RAXML=$5

mkdir -p $PREFIX

for i in {1..10}; do
    formatted_i=$(printf "%02d" $i)
    $RAXML -s $FASTA_FILE -n ${RF}.${formatted_i} -m GTRGAMMA -p $i -w ${PREFIX} -A ${MODEL}
done
