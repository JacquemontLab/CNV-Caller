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
#   <regions_file>     : BED-like file of genome annotation regions (e.g., telomeres)
#   <genome_version>   : Genome version string (e.g., GRCh37, GRCh38)
#   <output_file>      : Output filename for the merged CNV results
#
###############################################################################

set -euo pipefail

# --------------------- Check input arguments --------------------- #
if [ "$#" -ne 5 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <quantisnp_file> <penncnv_file> <regions_file> <genome_version> <output_file>"
    echo ""
    echo "Arguments:"
    echo "  <quantisnp_file>   : Path to QuantiSNP raw CNVs file (.txt.gz)"
    echo "  <penncnv_file>     : Path to PennCNV raw CNVs file (.txt.gz)"
    echo "  <regions_file>     : Path to BED-like genome annotation file"
    echo "  <genome_version>   : Genome version string (e.g., GRCh37, GRCh38)"
    echo "  <output_file>      : Path to output file for merged CNV table"
    exit 1
fi


# --------------------- Assign inputs --------------------- #
input_file_quantisnp="$1"
input_file_penncnv="$2"
regions_file="$3"
genome_version="$4"
output_file="$5"

log_step() {
    # Choose a symbol depending on the message
    case "$1" in
        STEP*) icon="🔹" ;;   # blue diamond for pipeline steps
        ERROR*) icon="❌" ;;  # red cross for errors
        WARN*) icon="⚠️ " ;;  # warning sign
        DONE*) icon="✅" ;;   # check mark for done
        *) icon="ℹ️ " ;;      # info icon
    esac
    echo -e "\n[$(date '+%Y-%m-%d %H:%M:%S')] $icon $1"
}

# Warn about Polars
log_step "WARN : Install Polars for Python"

# Try installing polars
if pip install --quiet polars; then
    log_step "DONE : Polars installed successfully"
else
    log_step "ERROR : Failed to install Polars"
    exit 1
fi

# Check if bedtools is available
if command -v bedtools >/dev/null 2>&1; then
    log_step "DONE : bedtools is available"
else
    log_step "ERROR : bedtools not found in PATH"
    exit 1
fi

# Get the directory of the currently running script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

### --------------------- INTERMEDIATE TEMPORARY FILES --------------------- ###

tmp_format_qs=$(mktemp --suffix=.tsv)
tmp_format_pc=$(mktemp --suffix=.tsv)
tmp_filter_qs=$(mktemp --suffix=.tsv)
tmp_filter_pc=$(mktemp --suffix=.tsv)
qs_clean=$(mktemp --suffix=.tsv)
pc_clean=$(mktemp --suffix=.tsv)
combined_bed=$(mktemp --suffix=.bed)
merged_bed=$(mktemp --suffix=.bed)
merged_tsv=$(mktemp --suffix=.tsv)
qs_bed=$(mktemp --suffix=.bed)
pc_bed=$(mktemp --suffix=.bed)
overlap_tsv=$(mktemp --suffix=.bed)
overlap_raw=$(mktemp --suffix=.tsv)
cnv_type=$(mktemp --suffix=.tsv)
overlap_problematic=$(mktemp --suffix=.tsv)

cleanup() {
    rm -f "$tmp_format_qs" "$tmp_format_pc" \
          "$tmp_filter_qs" "$tmp_filter_pc" \
          "$qs_clean" "$pc_clean" \
          "$merged_bed" "$merged_tsv" \
          "$qs_bed" "$pc_bed" \
          "$overlap_tsv" "$overlap_raw" "$cnv_type" "$overlap_problematic"
}
trap cleanup EXIT




### --------------------- STEP 1: Format QuantiSNP and PennCNV outputs --------------------- ###
log_step "STEP 1: Format QuantiSNP and PennCNV outputs"

# Check QuantiSNP file
if [[ "$(head -n 1 "$input_file_quantisnp" | cut -f1)" != "SampleID" ]]; then
  "$SCRIPT_DIR"/format_quantisnp_cnv.sh "$input_file_quantisnp" "$tmp_format_qs"
else
  echo "⚠️ Skipping QuantiSNP formatting: header detected in $input_file_quantisnp"
  cp $input_file_quantisnp "$tmp_format_qs"
fi

# Check PennCNV file
if [[ "$(head -n 1 "$input_file_penncnv" | cut -f1)" != "SampleID" ]]; then
  "$SCRIPT_DIR"/format_penncnv_cnv.sh "$input_file_penncnv" "$tmp_format_pc"
else
  echo "⚠️ Skipping PennCNV formatting: header detected in $input_file_penncnv"
  cp $input_file_penncnv "$tmp_format_pc"
fi 


### --------------------- STEP 2: Pre-filter by length >= 1000 bp and by score (30 for PennCNV and 15 for QuantiSNP) --------------------- ###
log_step "STEP 2: Pre-filter by length >= 1000 bp and by score (30 for PennCNV and 15 for QuantiSNP)"

"$SCRIPT_DIR"/filter_by_score_and_length.sh "$tmp_format_pc" 30 > "$tmp_filter_pc"
"$SCRIPT_DIR"/filter_by_score_and_length.sh "$tmp_format_qs" 15 > "$tmp_filter_qs"




### --------------------- STEP 3: Remove CNVs in pseudoautosomal regions --------------------- ###
log_step "STEP 3: Remove CNVs in pseudoautosomal regions on ChrX with Copy_Number = 2"
"$SCRIPT_DIR"/remove_PAR_regions.sh "$tmp_filter_qs" "$qs_clean" "$regions_file" "${genome_version}"
"$SCRIPT_DIR"/remove_PAR_regions.sh "$tmp_filter_pc" "$pc_clean" "$regions_file" "${genome_version}"


### --------------------- STEP 4: Merge CNVs from both PennCNV and QuantiSNP callers --------------------- ###
log_step "STEP 4: Merge CNVs from both PennCNV and QuantiSNP callers"

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




### --------------------- STEP 5: Compute Overlap With PennCNV, QuantiSNP, and Between Both --------------------- ###
log_step "STEP 5: Compute Overlap With PennCNV, QuantiSNP, and Between Both"

# Create BED file from overlap result with combined SampleID+Chr as first field
cat "$qs_clean" | awk 'BEGIN{OFS="\t"} NR>1{$1=$1","$2; $2=""; print}' | cut -f1,3,4 | sort -k1,1 -k2,2n > "$qs_bed"
cat "$pc_clean" | awk 'BEGIN{OFS="\t"} NR>1{$1=$1","$2; $2=""; print}' | cut -f1,3,4 | sort -k1,1 -k2,2n > "$pc_bed"


bedtools intersect -a <(cut -f1-3 "$qs_bed") -b <(cut -f1-3 "$pc_bed") -wo \
| awk '{start=($2 > $5 ? $2 : $5); end=($3 < $6 ? $3 : $6); if (start < end) print $1"\t"start"\t"end}' \
| awk -F'[\t,]' 'BEGIN {OFS="\t"; print "SampleID", "Chr", "Start", "End"} {print $1, $2, $3, $4}' > "$overlap_tsv"


# Set variable defining sources and files to overlap with
algo_to_overlap="QuantiSNP:$qs_clean,PennCNV:$pc_clean,Two_Algorithm:$overlap_tsv"

# Compute overlap with original CNV files
"$SCRIPT_DIR"/compute_cnv_overlap_fraction.sh "$merged_tsv" "$algo_to_overlap" "$overlap_raw"




### --------------------- STEP 6: Infer CNV Type --------------------- ###
log_step "STEP 6: Infer CNV Type Based on Copy Number"

"$SCRIPT_DIR"/infer_cnv_type.py -i "$overlap_raw" -o "$cnv_type"




### --------------------- STEP 7: Annotate With Problematic Regions Overlap --------------------- ###
log_step "STEP 7: Annotate With Problematic Regions Overlap"


# Create a temporary BED file for telomere
problematicregions_db=$(mktemp --suffix=.bed)
awk -v genome="$genome_version" 'BEGIN {
    OFS="\t"
    print "Chr", "Start", "End"
}
NR > 1 && $4 == "problematic_regions" && $5 == genome {
    print $1, $2, $3
}' "$regions_file" > "$problematicregions_db"


# Format string to pass to overlap computation script
regions_to_overlap="ProblematicRegions:$problematicregions_db"


# Annotate input CNV file with overlap metrics
"$SCRIPT_DIR"/compute_regions_overlap_fraction.sh "$cnv_type" "$regions_to_overlap" "$overlap_problematic"




### --------------------- STEP 8: Add Length Column and Generate Final Output --------------------- ###
log_step "STEP 8: Add Length Column and Generate Final Output"

# Create the Length column and save it in the final output
awk 'BEGIN {
    print "SampleID\tChr\tStart\tEnd\tType\tLength\tCopy_Number\tConfidence_max\tNum_Probes_max\tNum_Merged_CNVs\tQuantiSNP_Overlap\tPennCNV_Overlap\tTwo_Algorithm_Overlap\tProblematicRegions_Overlap"
}
NR > 1 {
    sample = $1
    chr = $2
    start = $3
    end = $4
    type = $5
    len = $4 - $3 + 1
    copynumber = $6
    conf_max = $7
    nb_probe = $8
    nb_cnv_merged = $9
    quantisnp_overlap = $10
    penncnv_overlap = $11
    two_algo_overlap = $12
    problematic_region_overlap = $13
    
    print sample"\t"chr"\t"start"\t"end"\t"type"\t"len"\t"copynumber"\t"conf_max"\t"nb_probe"\t"nb_cnv_merged"\t"quantisnp_overlap"\t"penncnv_overlap"\t"two_algo_overlap"\t"problematic_region_overlap
}' "$overlap_problematic" > "$output_file"

log_step "DONE : Output saved in $output_file"
