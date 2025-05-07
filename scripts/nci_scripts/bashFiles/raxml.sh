#!/bin/bash

RF=$1
FASTA_FILE=$2
PREFIX=$3
RAXML=$4

mkdir -p $PREFIX

for i in {11..20}; do
    formatted_i=$(printf "%02d" $i)
    $RAXML -s $FASTA_FILE -n ${RF}.${formatted_i} -m GTRGAMMA -p $i -w ${PREFIX}
done
