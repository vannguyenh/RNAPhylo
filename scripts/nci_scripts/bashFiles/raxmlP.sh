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
    #$RAXML -s $FASTA_FILE -n ${RF}.${formatted_i} -m GTRGAMMA -p $i -S $SS_BRACKETS_FILE -w ${PREFIX_BR}
    $RAXML -s $FASTA_FILE -n ${RF}.${formatted_i} -m GTRGAMMA -p $i -S $SS_FILE -w ${PREFIX}
done
