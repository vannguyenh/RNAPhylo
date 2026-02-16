#!/bin/bash

raxtree=$1
#raxPwPtree=$2
raxPiPtree=$2
prefix=$3
rna=$4

IQTREE="/Users/u7875558/tools/build-iqtree3/iqtree3"

${IQTREE} -rf ${raxtree} ${raxPiPtree} --prefix ${prefix}/${rna}
# full alignment
# ${IQTREE} -rf_all ${raxtree} --prefix ${prefix}/${rna}_dna