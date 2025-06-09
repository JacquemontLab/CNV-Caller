#!/bin/bash

# Usage: ./filter_cnvs.sh input_file score_threshold

set -euo pipefail

input_file="$1"
score_threshold="$2"

if [[ -z "$input_file" || -z "$score_threshold" ]]; then
  echo "Usage: $0 input_file score_threshold"
  exit 1
fi

awk -F'\t' -v score="$score_threshold" '
NR==1 { print; next }  # print header
($6 >= score) && ($8 > 1)  # Confidence >= score AND Length > 1
' "$input_file"
