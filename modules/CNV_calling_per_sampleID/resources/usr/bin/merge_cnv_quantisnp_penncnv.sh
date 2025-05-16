#!/bin/bash

###############################################################################
# CNV Integration Pipeline
#
# Author: Florian Bénitière
# Description: This script formats, cleans, merges, and analyzes CNV calls
#              from QuantiSNP and PennCNV, removes PAR regions, counts overlapping
#              probes, and computes overlap statistics.
#
# Requirements:
#   - bash, awk, bedtools, tail, cut, sort
#   - Custom scripts:
#       - format_quantisnp_cnv.sh
#       - format_penncnv_cnv.sh
#       - remove_PAR_regions.sh
#       - count_probes_per_cnv.sh
#       - compute_cnv_overlap_fraction.sh
#
# Usage:
#   ./merge_cnv_quantisnp_penncnv.sh <quantisnp_file> <penncnv_file> <probe_file> <regions_file> <genome_version> <output_file>
#
# Arguments:
#   <quantisnp_file>   : QuantiSNP CNV input file (gzipped QuantiSNP raw CNVs)
#   <penncnv_file>     : PennCNV CNV input file (gzipped PennCNV raw CNVs)
#   <probe_file>       : File containing BAF/LogR probe coordinates
#   <regions_file>     : BED-like file of genome annotation regions (e.g., telomeres)
#   <genome_version>   : Genome version string (e.g., GRCh37, GRCh38)
#   <output_file>      : Output filename for the merged CNV results
#
###############################################################################

set -euo pipefail

# --------------------- Check input arguments --------------------- #
if [ "$#" -ne 6 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <quantisnp_file> <penncnv_file> <probe_file> <regions_file> <genome_version> <output_file>"
    echo ""
    echo "Arguments:"
    echo "  <quantisnp_file>   : Path to QuantiSNP raw CNVs file (.txt.gz)"
    echo "  <penncnv_file>     : Path to PennCNV raw CNVs file (.txt.gz)"
    echo "  <probe_file>       : Path to BAF/LogR probe coordinates file"
    echo "  <regions_file>     : Path to BED-like genome annotation file"
    echo "  <genome_version>   : Genome version string (e.g., GRCh37, GRCh38)"
    echo "  <output_file>      : Path to output file for merged CNV table"
    exit 1
fi


# --------------------- Assign inputs --------------------- #
input_file_quantisnp="$1"
input_file_penncnv="$2"
probe_file="$3"
regions_file="$4"
genome_version="$5"
output="$6"

# Get the directory of the currently running script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

### --------------------- INTERMEDIATE TEMPORARY FILES --------------------- ###

tmp_qs=$(mktemp --suffix=.tsv.gz)
tmp_pc=$(mktemp --suffix=.tsv.gz)
qs_clean=$(mktemp --suffix=.tsv)
pc_clean=$(mktemp --suffix=.tsv)
combined_bed=$(mktemp --suffix=.bed)
merged_bed=$(mktemp --suffix=.bed)
merged_tsv=$(mktemp --suffix=.tsv)
probes_bed=$(mktemp --suffix=.bed)
probes_correct=$(mktemp --suffix=.tsv)
overlap_raw=$(mktemp --suffix=.tsv)
overlap_raw_sorted=$(mktemp --suffix=.bed)
qs_bed=$(mktemp --suffix=.bed)
pc_bed=$(mktemp --suffix=.bed)
overlap_bed=$(mktemp --suffix=.bed)
overlap_sum_bed=$(mktemp --suffix=.bed)
frac_bed=$(mktemp --suffix=.bed)
frac_only=$(mktemp)

cleanup() {
    rm -f "$tmp_qs" "$tmp_pc" "$qs_clean" "$pc_clean" "$combined_bed" "$merged_bed" \
          "$merged_tsv" "$probes_bed" "$probes_correct" "$overlap_raw" \
          "$overlap_raw_sorted" "$qs_bed" "$pc_bed" "$overlap_bed" "$overlap_sum_bed" \
          "$frac_bed" "$frac_only"
}
trap cleanup EXIT

### --------------------- STEP 1: Format QuantiSNP and PennCNV outputs --------------------- ###

# Format raw QuantiSNP and PennCNV outputs into standardized tabular files
"$SCRIPT_DIR"/format_quantisnp_cnv.sh "$input_file_quantisnp" "$tmp_qs"
"$SCRIPT_DIR"/format_penncnv_cnv.sh "$input_file_penncnv" "$tmp_pc"



### --------------------- STEP 2: Remove PAR regions --------------------- ###

# Remove CNVs in pseudoautosomal regions on ChrX with Copy_Number = 2 for each algorithm
"$SCRIPT_DIR"/remove_PAR_regions.sh "$tmp_qs" "$qs_clean" "$regions_file" "${genome_version}"
"$SCRIPT_DIR"/remove_PAR_regions.sh "$tmp_pc" "$pc_clean" "$regions_file" "${genome_version}"



### --------------------- STEP 3: Merge CNVs from both algorithms --------------------- ###

# Combine the formatted PennCNV and QuantiSNP files (excluding headers)
# Transform SampleID and Chr into one field to preserve both
# Then select necessary fields and sort for merging
( tail -n +2 "$qs_clean" && tail -n +2 "$pc_clean" ) | \
    awk 'BEGIN{OFS="\t"} {$1=$1","$2; $2=""; print}' | \
    cut -f1,3,4,5,6,7 | \
    sort -k1,1 -k2,2n > "$combined_bed"

# Merge overlapping CNVs using `bedtools merge`, tracking sample info and metrics, Num_Probes_max, Confidence_max etc.
# Merge overlapping regions in the BED file while aggregating additional columns:
# - Column 4 (Copy_Number): report all distinct values
# - Column 5 (Confidence): report the maximum value
# - Column 6 (Num_Probes): report the maximum value
# - Column 6 again (Num_Probes): count how many records were merged
bedtools merge -i "$combined_bed" -c 4,5,6,6 -o distinct,max,max,count > "$merged_bed"

# Add header and split back SampleID and Chr from field 1
(echo -e "SampleID\tChr\tStart\tEnd\tCopy_Number\tConfidence_max\tNum_Probes_max\tNum_Merged_CNVs" && \
    awk -F'\t' '{split($1, arr, ","); print arr[1] "\t" arr[2] "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' "$merged_bed") > "$merged_tsv"



### --------------------- STEP 4: Correct Probe Count --------------------- ###

# Extract probe coordinates if not NaN in Log R Ratio column, and sort them
awk '$4 != "NaN" {print "chr"$2"\t"$3"\t"$3}' "$probe_file" | tail -n +2 | sort -k1,1 -k2,2n > "$probes_bed"

# Report the total number of input probes and number of probes retained after filtering
original_probe_count=$(wc -l < "$probe_file")
filtered_probe_count=$(wc -l < "$probes_bed")
echo "Original probe count: $original_probe_count"
echo "Filtered BED probe count: $filtered_probe_count"

# Count number of probes overlapping each CNV
"$SCRIPT_DIR"/count_probes_per_cnv.sh "$merged_tsv" "$probes_correct" "$probes_bed"



### --------------------- STEP 5: Overlap With Original CNVs --------------------- ###

# Set variable defining sources and files to overlap with
algo_to_overlap="QuantiSNP:$qs_clean,PennCNV:$pc_clean"

# Compute overlap with original CNV files
"$SCRIPT_DIR"/compute_cnv_overlap_fraction.sh "$probes_correct" "$algo_to_overlap" "$overlap_raw"



### --------------------- STEP 6: Detect Overlap Between Tools --------------------- ###

# Create BED file from overlap result with combined SampleID+Chr as first field
cat "$overlap_raw" | awk 'BEGIN{OFS="\t"} NR>1{$1=$1","$2; $2=""; print}' | cut -f1,3,4 > "$overlap_raw_sorted"
cat "$qs_clean" | awk 'BEGIN{OFS="\t"} NR>1{$1=$1","$2; $2=""; print}' | cut -f1,3,4 | sort -k1,1 -k2,2n > "$qs_bed"
cat "$pc_clean" | awk 'BEGIN{OFS="\t"} NR>1{$1=$1","$2; $2=""; print}' | cut -f1,3,4 | sort -k1,1 -k2,2n > "$pc_bed"

# Intersect raw PennCNV and QuantiSNP CNVs and compute overlap count. Those we obtain Num bp shared by penncnv and quantisnp
bedtools intersect -a <(cut -f1-3 "$qs_bed") -b <(cut -f1-3 "$pc_bed") -wao | bedtools merge -i - -c 7 -o sum > "$overlap_bed"

# Overlap this with the CNV that have been merged, sum overlaps previously computed, and clean up nulls
bedtools intersect -a "$overlap_raw_sorted" -b "$overlap_bed" -wao | sed 's/\./0/g' | bedtools merge -i - -c 7 -o sum > "$overlap_sum_bed"

# Compute overlap fraction for each CNV
awk -F'\t' 'BEGIN{OFS="\t"} {
        size = $3 - $2;
        overlap = $NF;
        frac = overlap / size ;
        print $0, frac
    }' "$overlap_sum_bed" > "$frac_bed"

# Merge fractional overlap info with overlap_raw into final output
paste "$overlap_raw" <(echo -e "Two_Algorithm_Overlap" && cut -f5 "$frac_bed" ) | awk 'BEGIN {
    print "SampleID\tChr\tStart\tEnd\tLength\tCopy_Number\tConfidence_max\tNum_Probes\tNum_Merged_CNVs\tQuantiSNP_Overlap\tPennCNV_Overlap\tTwo_Algorithm_Overlap"
}
NR > 1 {  # Skip the first line (NR is the line number)
    sample = $1
    chr = $2
    start = $3
    end = $4
    len = $4 - $3 + 1
    copynumber = $5
    conf_max = $6
    nb_probe = $9
    nb_cnv_merged = $8
    quantisnp_overlap = $10
    penncnv_overlap = $11
    two_algo_overlap = $12
    
    print sample"\t"chr"\t"start"\t"end"\t"len"\t"copynumber"\t"conf_max"\t"nb_probe"\t"nb_cnv_merged"\t"quantisnp_overlap"\t"penncnv_overlap"\t"two_algo_overlap
}' > "$output"

