#!/usr/bin/env bash

# bash strict mode:
# -e: exit the script immediately if any command returns a non-zero status
# -u (nounset): treat any unset variable expression as an error and exit
# -o pipefail: the whole pipeline fails if any command fails (not just the last one)
set -euo pipefail

LOG="${1:-}"

if [ -z "$LOG"]; then
    echo "Usage: $0 /Users/u7875558/Documents/promotion/projects/RNAPhylo/fullAlignments_S6A/logs/full_S6A.log" >&2
    exit 1
fi 

# RNAs with "< 4 sequences"
grep -oE 'Number of sequences of RF[0-9]{5} is less than 4' "$LOG" \
  | grep -oE 'RF[0-9]{5}' | sort -u > less4.txt

# RNAs with "only unpaired bases"
grep -oE 'The secondary structure of RF[0-9]{5} has only unpaired bases' "$LOG" \
  | grep -oE 'RF[0-9]{5}' | sort -u > unpaired.txt

# RNAs that are in both lists (inputs to comm must be sorted)
comm -12 less4.txt unpaired.txt > both.txt

# Tiny summary
echo "Wrote:"
echo " - less4.txt        ($(wc -l < less4.txt | tr -d ' ' ) IDs)"
echo " - unpaired.txt     ($(wc -l < unpaired.txt | tr -d ' ' ) IDs)"
echo " - both.txt         ($(wc -l < both.txt | tr -d ' ' ) IDs)"