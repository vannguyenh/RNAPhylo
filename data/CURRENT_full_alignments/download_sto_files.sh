#!/bin/bash

for i in {1..4310}; do
  formatted_i=$(printf "%05d" $i)
  rf_file=${RF}.{formatted_i}
  wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/full_alignments/${rf_file}.sto
done
