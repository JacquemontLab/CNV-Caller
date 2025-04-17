#!/bin/bash

###############################################################################
# Script Name: count_probes_per_cnv.sh
# Description: This script takes a CNV file and intersects it with a BED file
#              of probe coordinates to count how many probes overlap each CNV.
#              Outputs a tab-delimited file with an added column "Num_Probes".
#
# Usage: ./count_probes_per_cnv.sh <input_file[.gz]> <output_file[.gz]> <probes_coordinates>
#
# Author: Florian Bénitière
# Date: April 2025
###############################################################################

# Check input arguments
if [ "$#" -ne 3 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <input_file[.gz]> <output_file[.gz]> <probes_coordinates.bed>"
    echo ""
    echo "Arguments:"
    echo "  <input_file[.gz]>     : The CNV file to process (either gzipped or uncompressed)."
    echo "                          Must contain at least: SampleID\tChr\tStart\tEnd"
    echo "  <output_file[.gz]>    : TSV output file where results are saved (gzipped or not)."
    echo "  <probes_coordinates>  : BED file containing probe regions to intersect with CNVs."
    exit 1
fi


# Assign input and output from arguments
input_file="$1"
output_file="$2"
probes_coordinates="$3"



# Determine how to read the input (gzipped or not)
if [[ "$input_file" == *.gz ]]; then
    read_cmd="zcat"
    # Optionally fallback to gunzip -c if zcat not available
    if ! command -v zcat &>/dev/null; then
        read_cmd="gunzip -c"
    fi
else
    read_cmd="cat"
fi

# Create temporary files
temp_cnv_bed=$(mktemp)
temp_intersect=$(mktemp)
temp_renamed=$(mktemp)
temp_merged=$(mktemp)
temp_reordered=$(mktemp)
final_output=$(mktemp)

# Extract header from input file
header_update=$($read_cmd "$input_file" | head -n 1)

# STEP 1: Transform input CNV file into BED-like format,  Output: [chr, start, end, full_line]
$read_cmd "$input_file" | tail -n +2 | awk 'BEGIN{OFS="\t"} { print $2, $3, $4, $0 }' | sort -k1,1 -k2,2n > "$temp_cnv_bed"

# STEP 2: Intersect CNV BED with probe BED to find overlapping probes
bedtools intersect -a "$temp_cnv_bed" -b "$probes_coordinates" -wao > "$temp_intersect"

# STEP 3: Rename CNV ID as "SampleID,Chr" for merging later
awk 'BEGIN{OFS="\t"} { $1 = $4 "," $1; print }' "$temp_intersect" | sort -k1,1 -k2,2n > "$temp_renamed"

# STEP 4: Merge by CNV to count number of overlapping probes
bedtools merge -i "$temp_renamed" -c 6 -o count > "$temp_merged"

# STEP 5: Prepare original input file with matching IDs for final merging
$read_cmd "$input_file" | tail -n +2 | awk 'BEGIN{OFS="\t"} { id = $1 "," $2; print id, $3, $4, $0 }' | sort -k1,1 -k2,2n > "$temp_reordered"


# STEP 6: Combine header, CNV content, and number of probes
(
    echo -e "$header_update\tNum_Probes"
    paste "$temp_reordered" <(cut -f4 "$temp_merged") | cut -f4-
) | {
    if [[ "$output_file" == *.gz ]]; then
        gzip > "$final_output"
    else
        cat > "$final_output"
    fi
}

# Move final result to desired output path
mv "$final_output" "$output_file"

# Cleanup temp files
rm -f "$temp_cnv_bed" "$temp_intersect" "$temp_renamed" "$temp_merged" "$temp_reordered"
