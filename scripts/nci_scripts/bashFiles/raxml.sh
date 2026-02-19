#!/bin/bash

RF=$1
FASTA_FILE=$2
PREFIX=$3
RAXML=$4

mkdir -p $PREFIX
# Run all 10 in parallel
for i in {11..20}; do
    formatted_i=$(printf "%02d" $i)
    $RAXML -s $FASTA_FILE -n ${RF}.${formatted_i} -m GTRGAMMA -p $i -w ${PREFIX} & #qsub command line here. -- each command line takes 1 cpus.
done

# Wait for all to finish
wait