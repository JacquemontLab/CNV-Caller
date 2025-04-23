#!/bin/bash

###############################################################################
# Script Name: remove_PAR_regions.sh
# Description: This script processes a CNV file, filters based on the genome version,
#              generates a BED file, removes rows that match PAR regions and with Copy_Number = 2, and
#              outputs the result while cleaning up intermediate files.
#
# Usage: ./remove_PAR_regions.sh <input_file[.gz]> <output_file> <genome_version> <PAR_regions_file>
#
# Author: Florian Bénitière
# Date: April 2025
###############################################################################


# Check that at least 2 arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <input_file[.gz]> <output_file> <genome_version> <PAR_regions_file>"
    echo ""
    echo "Arguments:"
    echo "  <input_file[.gz]>      : The input CNV file (either gzipped or uncompressed)."
    echo "                          Must contain at least the following columns:"
    echo "                          SampleID\tChr\tStart\tEnd\tCopy_Number"
    echo "  <output_file>          : The TSV output file where the results will be saved (e.g., 'output_results.tsv')."
    echo "  <genome_version>       : The genome version to filter by (e.g., GRCh37 or GRCh38)."
    echo "                          Default: GRCh38"
    echo "  <PAR_regions_file>     : Path to the PAR regions file in TSV format (optional)."
    echo "                          Default: /home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/ressources/PAR_regions.tsv"
    echo ""
    echo "Example Usage:"
    echo "  $0 input_data output_results.tsv GRCh38 /path/to/PAR_regions.tsv"
    exit 1
fi


# Parse input arguments
input_file="$1"
output_file="$2"
genome_version="${3:-GRCh38}"   # Default to GRCh38 if not provided
PAR_regions="${4:-/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/per_sampleID/resources/PAR_regions.tsv}" # Default PAR regions path

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
if [ "$PAR_regions" != "none" ]; then
    awk -v genome="$genome_version" 'BEGIN { OFS="\t" }
    NR > 1 && $1 == genome {
        print $3, $4, $5
    }' "$PAR_regions" > "$bed_PAR"
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