#!/bin/bash

IQTREE="/Users/u7875558/tools/iqtree-2.3.6-macOS/bin/iqtree2"

dna_tree=$1
rna_tree=$2
rna=$3

${IQTREE} -t ${dna_tree} -comp ${rna_tree} --prefix ${rna}
