#!/bin/bash

raxtree=$1
raxPtree=$2
#raxPwoPtree=$3
prefix=$3
rna=$4

IQTREE="/Users/u7875558/Documents/PhD/tools/iqtree-2.3.6-macOS/bin/iqtree2"

${IQTREE} -rf_all ${raxtree}
${IQTREE} -rf_all ${raxPtree}
##${IQTREE} -rf_all ${raxPwoPtree}

${IQTREE} -rf ${raxtree} ${raxPtree} --prefix ${prefix}/${rna}.raxml.raxmlP
#${IQTREE} -rf ${raxPtree} ${raxPtree} --prefix ${prefix}/${rna}.raxmlP.raxmlP
#${IQTREE} -rf ${raxtree} ${raxtree} --prefix ${prefix}/${rna}.raxml.raxml
#${IQTREE} -rf ${raxPwPtree} ${raxPwoPtree} --prefix ${prefix}/${rna}.raxmlPw.ramxlPwo
