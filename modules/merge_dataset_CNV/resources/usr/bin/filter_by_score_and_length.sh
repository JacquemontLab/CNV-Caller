#!/bin/bash

# filter_cnvs.sh
#
# Description:
#   Filters CNV calls based on a minimum confidence score and CNV length.
#   Keeps only entries where:
#     - Column 6 (Confidence score) >= score_threshold
#     - Column 8 (CNV length) >= 1000 bp
#
# Usage:
#   ./filter_cnvs.sh input_file score_threshold
#
# Arguments:
#   input_file       Path to the input TSV file containing CNV calls.
#   score_threshold  Minimum confidence score (numeric).
#
# Example:
#   ./filter_cnvs.sh cnvs.tsv 30


set -euo pipefail

input_file="$1"
score_threshold="$2"

if [[ -z "$input_file" || -z "$score_threshold" ]]; then
  echo "Usage: $0 input_file score_threshold"
  exit 1
fi

awk -F'\t' -v score="$score_threshold" '
NR==1 { print; next }  # print header
($6 >= score) && ($8 >= 1000)  # Confidence >= score AND Length >= 1kb
' "$input_file"
