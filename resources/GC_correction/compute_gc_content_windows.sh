#!/bin/bash

###############################################################################
# Script Name   : compute_gc_content_windows.sh
# Description   : Computes GC content across fixed-width genomic windows.
#                 Generates chromosome sizes from a compressed genome FASTA file,
#                 creates non-overlapping genomic windows, and computes GC content
#                 for each window using `bedtools nuc`.
#                 Outputs a BED file with chromosome, start, end, and GC fraction.
#
# Usage         : ./compute_gc_content_windows.sh <genome_fasta.bgz> <output_file> <window_width>
# Arguments     :
#   <genome_fasta.bgz>      - Compressed genome FASTA file (bgzipped)
#   <output_file>           - Output BED file with columns: chrom, start, end, GC_fraction
#   <window_width>          - (Optional) Width of genomic windows (default: 999 bp)
#
# Requirements  : bedtools (v2.31+), awk, sed, zcat
#
# Author        : Florian Bénitière
# Date          : 2025-04-16
###############################################################################

set -euo pipefail

# Input arguments
genome_used="$1"
output_file="$2"
width="${3:-999}"  # default to 999 if not specified

# Create temp files with unique names
chr_sizes_tmp=$(mktemp --suffix=_chr_sizes.tsv)
windows_tmp=$(mktemp --suffix=_windows.bed)
nuc_tmp=$(mktemp --suffix=_nuc.bed)

# Clean up temp files on exit
trap "rm -f '$chr_sizes_tmp' '$windows_tmp' '$nuc_tmp'" EXIT


# Step 1: Generate chromosome sizes from genome fasta
zcat "$genome_used" | awk '
    /^>/ {
        split($0, arr, " ");
        chr = arr[1];
        chr = substr(arr[1], 2);
        split($0, arrscnd, ":");
        seq_length = arrscnd[6];
        print chr "\t" seq_length
    }' > "$chr_sizes_tmp"


# Step 2: Create fixed-width genomic windows (999 bp) over the genome
bedtools makewindows -g "$chr_sizes_tmp" -w ${width} -s $((width + 1)) > "$windows_tmp"

# Step 3: Compute GC content in each window
bedtools nuc -fi "$genome_used" -bed "$windows_tmp" > "$nuc_tmp"

# Step 4: Extract relevant columns: chr, start, end, GC fraction
awk -F "\t" -v OFS="\t" 'FNR > 1 {
    gc_num = $7 + $8;                    # G_count + C_count
    atgc_num = $6 + $7 + $8 + $9;        # A + G + C + T
   if (atgc_num == 0) {
        gc_frac = "NA";
    } else {
        gc_frac = gc_num / atgc_num;
    }
    print $1, $2, $3, gc_frac;
}' "$nuc_tmp" > "$output_file"