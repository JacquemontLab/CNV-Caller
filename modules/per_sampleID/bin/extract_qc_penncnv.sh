#!/bin/bash

# Florian Bénitière and Zohra Saci
# 11 March 2025
# Extract quality score from multiple PennCNV log output to produce a TSV file with one sample per line.

# Usage: ./script.sh <input_directory> <output_file>
# Example: ./script.sh pcnv_log QC_per_IID.tsv

# Ensure correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_file>"
    exit 1
fi

# Assign arguments to variables
INPUT_DIR="$1"   # Directory containing .pcnv.log.gz files
OUTPUT_FILE="$2" # Output TSV file

# Write the header to the output file
echo -e "SampleID\tLRR_mean\tLRR_median\tLRR_SD\tBAF_mean\tBAF_median\tBAF_SD\tBAF_DRIFT\tGCWF\tWF" > "$OUTPUT_FILE"

# Extract relevant data from compressed log files, format it, and append to output file
for file in "$INPUT_DIR"/*.pcnv.log.gz; do
    zcat "$file" \
    | grep "quality summary" \
    | sed -r 's/NOTICE: quality summary for IID_//g; s/_vcf.tsvlite://g; s/LRR_mean=//g; s/LRR_median=//g; s/LRR_SD=//g; s/BAF_mean=//g; s/BAF_median=//g; s/BAF_SD=//g; s/BAF_DRIFT=//g; s/GCWF=//g; s/WF=//g' \
    | sed -r 's/\s+/\t/g' >> "$OUTPUT_FILE"
done


# Print completion message
echo "Processing complete. Results saved in $OUTPUT_FILE"