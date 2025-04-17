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


# Step 1: Generate chromosome sizes from genome fasta
zcat "$genome_used" | awk '
    /^>/ {
        split($0, arr, " ");
        chr = arr[1];
        chr = substr(arr[1], 2);
        split($0, arrscnd, ":");
        seq_length = arrscnd[6];
        print chr "\t" seq_length
    }' > genome_chr_sizes.tsv


# Step 2: Create fixed-width genomic windows (999 bp) over the genome
bedtools makewindows -g genome_chr_sizes.tsv -w ${width} -s $((width + 1)) > genome_${width}bps.bed

# Step 3: Compute GC content in each window
bedtools nuc -fi "$genome_used" -bed genome_${width}bps.bed > genome_nuc_${width}.bed

# Step 4: Extract relevant columns: chr, start, end, GC fraction
awk -F "\t" -v OFS="\t" '{if (FNR>1) {print $1, $2, $3, $5}}' genome_nuc_${width}.bed > "$output_file"

rm genome_nuc_${width}.bed genome_${width}bps.bed genome_chr_sizes.bed