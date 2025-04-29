#!/bin/bash

###############################################################################
# Script Name: compute_cnv_overlap_fraction.sh
# Description: This script processes a CNV file (gzipped or uncompressed) and
#              intersects it with multiple region-generated BED files to
#              calculate overlap fractions. The results are stored in a
#              tab-delimited, gzipped TSV file for downstream analysis.
#
# Usage: ./compute_regions_overlap_fraction.sh <input_file[.gz]> <regions_to_overlap> <output_file[.gz]>
#
# Author: Florian Bénitière
# Date: April 2025
###############################################################################

set -euo pipefail

# Check input arguments
if [ "$#" -ne 3 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <input_file[.gz]> <regions_to_overlap[.gz]> <output_file[.gz]>"
    echo ""
    echo "Arguments:"
    echo "  <input_file[.gz]>       : The input CNV file (either gzipped or uncompressed)."
    echo "                            Must contain at least the following columns:"
    echo "                            Chr\tStart\tEnd"
    echo "  <regions_to_overlap> : A comma-separated list of regions with corresponding BED files."
    echo "                            Each entry should be in the format 'first_region_name:bed_file,second_region_name:bed_file'."
    echo "                            Each BED file must contain at least the following columns:"
    echo "                            Chr\tStart\tEnd"
    echo "  <output_file[.gz]>      : The TSV output file where the results will be saved (either gzipped or uncompressed)."
    exit 1
fi



### --------------------- STEP 7: Overlap With Genomes Regions --------------------- ###
cd ~/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/merge_dataset_CNV/test

genome_version=GRCh37
regions_file=/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/CNV_calling_per_sampleID/resources/Genome_Regions_data.tsv

telomere_db=$(mktemp --suffix=.bed)
awk -v genome="$genome_version" 'BEGIN {
    OFS="\t"
    print "Chr", "Start", "End"
}
NR > 1 && $4 == "telomere" && $5 == genome {
    print $1, $2, $3
}' "$regions_file" > "$telomere_db"

centromere_db=$(mktemp --suffix=.bed)
awk -v genome="$genome_version" 'BEGIN {
    OFS="\t"
    print "Chr", "Start", "End"
}
NR > 1 && $4 == "centromere" && $5 == genome {
    print $1, $2, $3
}' "$regions_file" > "$centromere_db"

segmentaldup_db=$(mktemp --suffix=.bed)
awk -v genome="$genome_version" 'BEGIN {
    OFS="\t"
    print "Chr", "Start", "End"
}
NR > 1 && $4 == "segmentaldup" && $5 == genome {
    print $1, $2, $3
}' "$regions_file" > "$segmentaldup_db"

# Set variable defining sources and files to overlap with
regions_to_overlap="telomere:$telomere_db,centromere:$centromere_db,segmentaldup:$segmentaldup_db"

cut -f2- "$output" > reduce.tsv

# Compute overlap with original CNV files
"$SCRIPT_DIR"/compute_regions_overlap_fraction.sh reduce.tsv "$algo_to_overlap" test


# Clean up the temporary files
rm "$bed_PAR" "$telomere_db" "$centromere_db" "$segmentaldup_db"

# Assign input and output from arguments
input_file="$1"
input_file=test.tsv
regions_to_overlap="$2"
output_file="$3"


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

# Pre-process CNV coordinates:
temp_sorted_cnv_bed=$(mktemp)
temp_cnv_bed=$(mktemp)
temp_saved_for_merged_cnv_bed=$(mktemp)

# - Combine columns 1 and 2 into a single column comma separated
# - Sort the data based on columns 1 and 2
$read_cmd "$input_file" | tail -n +2 | sort -k1,1 -k2,2n > "$temp_sorted_cnv_bed"
header_update=$($read_cmd "$input_file" | head -n 1)
# - Extract the chr, start, and end columns
cut -f2,3,4 "$temp_sorted_cnv_bed" > "$temp_cnv_bed"

$read_cmd "$input_file" | tail -n +2 | awk 'BEGIN{OFS="\t"}
{
    print $2, $3, $4, $1
}' | sort -k1,1 -k2,2n > "$temp_cnv_bed"

tail -n +2 "$input_file" | awk 'BEGIN{OFS="\t"} {$1=$1","$2; $2=""; print}' - | sort -k1,1 -k2,2n > "$temp_saved_for_merged_cnv_bed"


# Loop through the regions and their corresponding BED files
for item in ${regions_to_overlap//,/ }; do
    
    # Extract region name (before the colon)
    algo=${item%%:*}
    # Extract BED file name (after the colon)
    tsvfile=${item##*:}

    # Process the BED file for the current region:
    # - Combine columns 1 and 2 into a single column
    # - Extract columns chr, start, and end
    # - Sort the data based on columns 1 and 2
    temp_sorted_algo_bed=$(mktemp)

    # Determine whether the tsvfile file is gzipped
    if [[ "$tsvfile" == *.gz ]]; then
        read_cmd="zcat"
        # Optionally fallback to gunzip -c if zcat not available
        if ! command -v zcat &>/dev/null; then
            read_cmd="gunzip -c"
        fi
    else
        read_cmd="cat"
    fi

    $read_cmd "$tsvfile" | tail -n +2 | sort -k1,1 -k2,2n > "$temp_sorted_algo_bed"
    
    # Use bedtools to find intersections between CNV coordinates and the current region's BED file:
    # - Perform a "wao" (write all overlaps) intersection
    # - Merge the results, summing the overlap column
    temp_non_overlap_bed=$(mktemp)
    temp_overlap_bed=$(mktemp)
    bedtools merge -i "$temp_sorted_algo_bed" > "$temp_non_overlap_bed"
    bedtools intersect -a "$temp_cnv_bed" -b "$temp_non_overlap_bed" -wao > "$temp_overlap_bed"


    temp_sorted_cnv_bed=$(mktemp)
    awk 'BEGIN{OFS="\t"} {$1=$4","$1; $4=""; print}' "$temp_overlap_bed" | sort -k1,1 -k2,2n | bedtools merge -i - -c 8 -o sum > "$temp_sorted_cnv_bed"

    # Calculate the fraction of overlap and append this information:
    # - Calculate the size of the interval (difference between start and end positions)
    # - Calculate the overlap as a fraction of the size
    temp_frac_overlap_bed=$(mktemp)
    awk -F'\t' 'BEGIN{OFS="\t"} {
            size = $3 - $2;
            overlap = $NF;
            frac = overlap / size ;
            print $0, frac
        }' "$temp_sorted_cnv_bed" > "$temp_frac_overlap_bed"

    # Merge the CNV coordinates with the calculated overlap fraction
    temp_tmp_merging_algo_bed=$(mktemp)
    paste "$temp_sorted_cnv_bed" <(cut -f5 "$temp_frac_overlap_bed") > "$temp_tmp_merging_algo_bed"

    # Update the main result file with the merged data
    mv "$temp_tmp_merging_algo_bed" "$temp_sorted_cnv_bed"

    # Append the region name (with tab separation) to the header string
    header_update="$header_update\t${algo}_Overlap"

    # Clean up temporary files created during this iteration
    rm -f "$temp_sorted_algo_bed" "$temp_overlap_bed" "$temp_frac_overlap_bed" "$temp_non_overlap_bed"
done

# Update the final output with the new header and the results
# - Add the region names to the header
# - Merge the CNV data with the newly calculated overlap fractions
final_output=$(mktemp)
( echo -e "$header_update" && \
cat "$temp_sorted_cnv_bed" ) | {
    if [[ "$output_file" == *.gz ]]; then
        gzip > "$final_output"
    else
        cat > "$final_output"
    fi
}

mv "$final_output" "$output_file"