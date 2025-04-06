#!/bin/bash

raxtree=$1
raxPwPtree=$2
raxPiPtree=$3
prefix=$4
rna=$5

IQTREE="/Users/u7875558/Documents/PhD/tools/iqtree-2.3.6-macOS/bin/iqtree2"

${IQTREE} -rf_all ${raxtree}
${IQTREE} -rf_all ${raxPwPtree}
${IQTREE} -rf_all ${raxPiPtree}

${IQTREE} -rf ${raxtree} ${raxPwPtree} --prefix ${prefix}/${rna}.raxml.raxmlPw
#${IQTREE} -rf ${raxPwPtree} ${raxPwPtree} --prefix ${prefix}/${rna}.raxmlPw.raxmlPw
${IQTREE} -rf ${raxtree} ${raxPiPtree} --prefix ${prefix}/${rna}.raxml.raxmlPi
${IQTREE} -rf ${raxPwPtree} ${raxPiPtree} --prefix ${prefix}/${rna}.raxmlPw.raxmlPi