#!/bin/bash

###############################################################################
# Script Name: remove_PAR_regions.sh
# Description: This script processes a CNV file and removes entries that overlap 
#              with pseudoautosomal regions (PAR or XTR) for a given genome version.
#              Specifically, it filters out regions with Copy_Number = 2 within PAR/XTR,
#              as such duplications are not considered reliable due to their presence
#              on both the X and Y chromosomes.
#
# Usage: ./remove_PAR_regions.sh <input_file[.gz]> <output_file> <regions_file> <genome_version>
#
# Author: Florian Bénitière
# Date: April 2025
###############################################################################

# Check that at least 2 arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <input_file[.gz]> <output_file> <regions_file> <genome_version>"
    echo ""
    echo "Arguments:"
    echo "  <input_file[.gz]>      : The input CNV file (either gzipped or uncompressed)."
    echo "                          Must contain at least the following columns:"
    echo "                          SampleID\tChr\tStart\tEnd\tCopy_Number"
    echo "  <output_file>          : The TSV output file where the results will be saved (e.g., 'output_results.tsv')."
    echo "  <regions_file>          : Path to the genome regions file in TSV format."
    echo "  <genome_version>       : The genome version to filter by (e.g., GRCh37 or GRCh38)."
    echo "                          Default: GRCh38"
    echo ""
    echo "Example Usage:"
    echo "  $0 input_data output_results.tsv /path/to/regions_file.tsv GRCh38"
    exit 1
fi


# Parse input arguments
input_file="$1"
output_file="$2"
regions_file="$3" # Default PAR regions path
genome_version="${4:-GRCh38}"   # Default to GRCh38 if not provided

# Determine whether the input file is gzipped
if [[ "$input_file" == *.gz ]]; then
    read_cmd="zcat"
    # Optionally fallback to gunzip -c if zcat not available
    if ! command -v zcat &>/dev/null; then
        read_cmd="gunzip -c"
    fi
else
    read_cmd="cat"
fi

# Extract header for later reinsertion
header_update=$($read_cmd "$input_file" | head -n 1)

# Create a temporary BED file with PAR regions matching the genome version
bed_PAR=$(mktemp)
if [ "$regions_file" != "none" ]; then
    awk -v genome="$genome_version" 'BEGIN { OFS="\t" }
    NR > 1 && ($4 == "PAR1" || $4 == "PAR2" || $4 == "XTR") && $5 == genome {
        print $1, $2, $3
    }' "$regions_file" > "$bed_PAR"
fi

# Create temporary files for CNVs in chrX with copy number 2 (possibly in PAR) vs all other CNVs
matched_file_bed=$(mktemp)
non_matched_file_bed=$(mktemp)
# Read input data (excluding header) and split based on criteria
$read_cmd "$input_file" | tail -n +2 | awk -v OFS="\t" -v out1="$matched_file_bed" -v out2="$non_matched_file_bed" '
{
    if ($2 == "chrX" && $5 == 2){
        print $2, $3, $4, $0 > out1;
        }
    else {
        print $2, $3, $4, $0 > out2;
    }
}'


# Create a temporary file to hold final intersect result (PAR regions will be removed)
final_output_bed=$(mktemp)
# Use bedtools to remove CNVs that intersect with PAR regions
bedtools intersect -v -a "$matched_file_bed" -b "$bed_PAR" > "$final_output_bed"

# Combine the header + cleaned CNVs from matched_file + all others
( echo -e "$header_update" && cut -f4- "$final_output_bed" "$non_matched_file_bed") > "$output_file"

# Clean up the temporary files
rm "$bed_PAR" "$final_output_bed" "$matched_file_bed" "$non_matched_file_bed"