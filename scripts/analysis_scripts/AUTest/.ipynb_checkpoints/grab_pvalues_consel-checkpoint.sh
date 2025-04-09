#!/bin/bash

PATH_CONSEL='/Users/u7875558/Documents/PhD/tools/consel/bin/'

pv_file=$1

${PATH_CONSEL}/catpv ${pv_file} | awk '{print $5}' | head -4 | tail -1
