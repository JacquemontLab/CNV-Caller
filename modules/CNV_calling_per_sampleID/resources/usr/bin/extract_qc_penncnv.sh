#!/bin/bash

# Florian Bénitière and Zohra Saci
# 11 March 2025
# Extract quality score from multiple PennCNV qc output to produce a TSV file with one sample per line.

# Usage: ./script.sh <input_directory> <output_file>
# Example: ./script.sh pcnv_log QC_per_IID.tsv

set -euo pipefail

# Ensure correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_file>"
    exit 1
fi

# Assign arguments to variables
INPUT_FILE="$1"   # .penncnv.qc file
OUTPUT_FILE="$2" # Output TSV file

# Write the header to the output file
echo -e "SampleID\tLRR_mean\tLRR_median\tLRR_SD\tBAF_mean\tBAF_median\tBAF_SD\tBAF_DRIFT\tWF\tGCWF" > "$OUTPUT_FILE"

# Extract relevant data from compressed log files, format it, and append to output file
echo -e "$(basename "$INPUT_FILE" .penncnv.qc)\t$(zgrep -P "quality summary" "$INPUT_FILE" \
| sed -r 's/.*LRR_mean=//g; s/LRR_median=//g; s/LRR_SD=//g; s/BAF_mean=//g; s/BAF_median=//g; s/BAF_SD=//g; s/BAF_DRIFT=//g; s/GCWF=//g; s/WF=//g' \
| sed -r 's/\s+/\t/g')" >> "$OUTPUT_FILE"


# Print completion message
echo "Processing complete. Results saved in $OUTPUT_FILE"