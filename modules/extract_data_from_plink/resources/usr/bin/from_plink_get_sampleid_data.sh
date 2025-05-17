#!/bin/bash

###############################################################################
# Script Name   : from_plink_get_sampleid_data.sh
# Description   : Computes sample call rates and infers sex from PLINK genotype files.
#                 Includes family information from the PLINK .fam file.
#                 Filters probes based on a provided list of BAF/LRR probes.
#                 Outputs a TSV file.
#
# https://www.cog-genomics.org/plink/1.9/formats#sexcheck
# This link refers to documentation on PLINK's .sexcheck output format,
# which is used to validate or infer sex from genetic data.
#
# Usage         : ./from_plink_get_sampleid_data.sh <plink_path_prefix> <baf_lrr_probes_file> <output_file>
# Arguments     :
#   <plink_path_prefix>      - Path and prefix of the PLINK binary files (.bed/.bim/.fam)
#   <baf_lrr_probes_file>    - TSV file listing the probes to retain (with header)
#   <output_file>            - Output TSV file with SampleID, Sex_from_plink, and Call_Rate, FID (Family ID), PID (Paternal ID, Father), MID (Maternal ID, Mother)
#
# Requirements  : PLINK must be available in the system PATH
#                 Bash shell (for process substitution or tmpfile workaround)
#
# Author        : Florian Bénitière
# Date          : 2025-04-16
###############################################################################


set -euo pipefail

# Function to show usage
usage() {
  echo "Usage: $0 <plink_path_prefix> <baf_lrr_probes_file> <output_file>" >&2
  exit 1
}

# Check number of arguments
if [ "$#" -ne 3 ]; then
  usage
fi

# Assign input and output from arguments
plink_path_prefix="$1"
one_baf_lrr_probes_file="$2"
output_file="$3"

# Step 1: Filter PLINK file to include only probes listed in `one_baf_lrr_probes_file`
#plink --bfile "$plink_path_prefix" --extract <(cut -f1 "$one_baf_lrr_probes_file" | tail -n +2) --make-bed --out filtered_plink

# Step 2: Compute missingness and perform sex check on the filtered PLINK file
plink --bfile "$plink_path_prefix" --missing --check-sex --out plink_outputs

# Step 3: Process .imiss file to extract sample call rates.
(echo -e "SampleID\tCall_Rate" && awk 'NR>1 {gsub(/[[:space:]]+/, "\t"); OFS="\t"; print $2"\t"(1-$6)}' plink_outputs.imiss) > plink_outputs.imiss.tsv

# Step 4: Process .sexcheck file to get PLINK-imputed sex
(echo -e "SampleID\tSex_from_checksex" && awk 'NR>1 {
    gsub(/[[:space:]]+/, "\t"); 
    OFS = "\t";
    sex = ($4 == 1 ? "male" : ($4 == 2 ? "female" : "unknown"));
    print $2, sex;
}' plink_outputs.sexcheck ) > plink_outputs.sexcheck.tsv


# Step 5: Extract family information from the .fam file
# This part includes the family information such as FID, SampleID, PID, MID, SEX, and PHENOTYPE
(echo -e "FID\tSampleID\tPID\tMID\tSEX\tPHENOTYPE" && awk '{gsub(/[[:space:]]+/, "\t"); OFS="\t"; print }' "$plink_path_prefix".fam) > filtered_plink.family.tsv


# Step 6: Combine sex info, call rate and family information into a final summary TSV
# (Assumes the order of samples is consistent between files)

# Checking for mismatches between sexcheck and imiss sample IDs
if diff <(cut -f1 plink_outputs.sexcheck.tsv) <(cut -f1 plink_outputs.imiss.tsv) > /dev/null; then
  # Checking for mismatch with family data
  if diff <(cut -f1 plink_outputs.sexcheck.tsv) <(cut -f2 filtered_plink.family.tsv); then
    # Combining the information
    paste plink_outputs.sexcheck.tsv <(cut -f2 plink_outputs.imiss.tsv) > "$output_file"
    # Append family information to the output file
    paste "$output_file" <(cut -f1,3,4 filtered_plink.family.tsv) > temp_output_file
    mv temp_output_file "$output_file"  # Rename temporary file to final output file
    
    else
    rm -f plink_outputs.* filtered_plink.* extracted_probes.txt
    echo "❌ SampleID mismatch between .sexcheck and .fam. Aborting." >&2
  fi
else
  rm -f plink_outputs.* filtered_plink.* extracted_probes.txt
  echo "❌ SampleID mismatch between .sexcheck and .imiss. Aborting." >&2
fi
