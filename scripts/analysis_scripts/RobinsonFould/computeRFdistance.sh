#!/bin/bash

raxtree=$1
#raxPwPtree=$2
raxPiPtree=$2
prefix=$3
rna=$4

IQTREE="/Users/u7875558/Documents/promotion/projects/projects_code/iqtree3/build/iqtree3"

${IQTREE} -rf_all ${raxtree}
#${IQTREE} -nt 1 -rf_all ${raxPwPtree}
${IQTREE} -rf_all ${raxPiPtree}

#${IQTREE} -rf ${raxtree} ${raxPwPtree} --prefix ${prefix}/${rna}.raxml.raxmlPw
${IQTREE} -rf ${raxtree} ${raxPiPtree} --prefix ${prefix}/${rna}.raxml.raxmlPi
#${IQTREE} -rf ${raxPwPtree} ${raxPiPtree} --prefix ${prefix}/${rna}.raxmlPw.raxmlPi
