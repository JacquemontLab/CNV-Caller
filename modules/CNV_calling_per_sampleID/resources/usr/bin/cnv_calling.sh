#!/bin/bash

# =============================================================================
# CNV Analysis Pipeline Wrapper
# =============================================================================
# Description:
#   This script runs CNV detection on SNP array data using two algorithms:
#   PennCNV and QuantiSNP. It processes both autosomal and sex chromosomes.
#
# Inputs:
#   --sample_id        Sample identifier (required)
#   --BAF_LRR_Probes   Input file containing B Allele Frequency (BAF) and
#                      Log R Ratio (LRR) values (required)
#   --sexfile          Tab-delimited file with sample ID and gender (required)
#
# PennCNV-specific inputs:
#   --pfb              Population Frequency of B Allele file (.pfb)
#   --hmm              Hidden Markov Model file (.hmm)
#   --gcmodel          GC model file for GC correction (.gcmodel)
#
# QuantiSNP-specific inputs:
#   --chr              Chromosome to process (e.g., 1, 2, X)
#   --levels           Levels file for QuantiSNP
#   --config           Configuration file for QuantiSNP
#   --gcdir            Directory containing GC content files for each chromosome
#
# Output:
#   - CNV and LOH files , renamed for clarity
#   - Log files and merged CNV result from both tools
#
# Usage:
#   ./cnv_calling.sh [options]
#
# Example:
#   ./cnv_calling.sh \
#     --sample_id sample123 \
#     --BAF_LRR_Probes probes.tsv \
#     --sexfile sexfile.tsv \
#     --pfb pfb.tsv \
#     --hmm model.hmm \
#     --gcmodel gcmodel.tsv \
#     --chr 1 \
#     --levels levels.dat \
#     --config params.dat \
#     --gcdir /path/to/gcdir
#
# =============================================================================

set -euo pipefail

# Help function
define_usage() {
  echo "Usage: $0 [--sample_id ID] [--BAF_LRR_Probes FILE] [--sexfile FILE] \\
                [--pfb FILE] [--hmm FILE] [--gcmodel FILE] [--chr CHR] \\
                [--levels FILE] [--config FILE] [--gcdir DIR]"
  exit 1
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample_id) sample_id="$2"; shift 2;;
    --BAF_LRR_Probes) BAF_LRR_Probes="$2"; shift 2;;
    --sexfile) sex_file="$2"; shift 2;;
    # ================================
    # PennCNV-specific input parameters
    --pfb) pfb_file="$2"; shift 2;;
    --hmm) hmm_file="$2"; shift 2;;
    --gcmodel) gcmodel_file="$2"; shift 2;;
    # ================================
    # QuantiSNP-specific input parameters
    --chr) chr="$2"; shift 2;;
    --levels) levels="$2"; shift 2;;
    --config) config="$2"; shift 2;;
    --gcdir) gcdir="$2"; shift 2;;
    -h|--help) define_usage;;
    *) echo "Unknown option: $1"; define_usage;;
  esac
done

echo "Process Running: CNV_calling for ${sample_id}"

# -----------------------------------------------------------------------------
# Step 1: Extract Gender
# -----------------------------------------------------------------------------
echo "Extracting gender from sex file..."

sexfile_temp=$(mktemp)
zgrep -P "${sample_id}\t" ${sex_file} | cut -f1,2 > ${sexfile_temp}
gender=$(cut -f2 ${sexfile_temp} | tr '[:upper:]' '[:lower:]')
if [[ -z "$gender" ]]; then
    echo "ERROR: No gender information found for sample ID ${sample_id}."
    exit 1
fi

echo "Gender: $gender"
if [[ "$gender" != "male" && "$gender" != "female" ]]; then
    echo "ERROR: Invalid gender '$gender' for sample ID ${sample_id}. Expected 'male' or 'female'."
    exit 1
fi


# -----------------------------------------------------------------------------
# Step 2: PennCNV (Autosomes)
# -----------------------------------------------------------------------------

# Note:
# - `--minsnp 1` is explicitly set because PennCNV defaults to 3 SNPs per CNV,
#   whereas QuantiSNP applies no such minimum filtering.
echo "Running PennCNV on autosomes..."
/usr/bin/time -v perl /usr/local/bin/detect_cnv.pl --test \
    --conf \
    --minsnp 1 \
    --pfbfile ${pfb_file} \
    --hmmfile ${hmm_file} \
    --logfile ${sample_id}.penncnv.log \
    --output ${sample_id}.penncnv.out \
    --gcmodelfile ${gcmodel_file} \
    ${BAF_LRR_Probes}
    
    
# -----------------------------------------------------------------------------
# Step 3: PennCNV (ChrX)
# -----------------------------------------------------------------------------
echo "Running PennCNV on ChrX..."
/usr/bin/time -v perl /usr/local/bin/detect_cnv.pl --test \
    --conf \
    --minsnp 1 \
    --pfbfile ${pfb_file} \
    --hmmfile ${hmm_file} \
    --logfile ${sample_id}.penncnv.chrx.log \
    --output ${sample_id}.penncnv.chrx.out \
    --gcmodelfile ${gcmodel_file} \
    --sexfile ${sexfile_temp} \
    --chrx \
    ${BAF_LRR_Probes}


# -----------------------------------------------------------------------------
# Step 4: QuantiSNP
# -----------------------------------------------------------------------------
echo "Running QuantiSNP..."
export MCR_CACHE_ROOT=$(mktemp -d -t mcr_cache_XXXXXX)
mkdir -p "$MCR_CACHE_ROOT"

/usr/bin/time -v quantisnp --chr ${chr} \
    --outdir . \
    --sampleid ${sample_id} \
    --gender ${gender} \
    --config ${config} \
    --levels ${levels} \
    --gcdir ${gcdir} \
    --input-files ${BAF_LRR_Probes} \
    --doXcorrect --verbose &> ${sample_id}.quantisnp.log

rm -rf "$MCR_CACHE_ROOT"

# -----------------------------------------------------------------------------
# Step 4: Renaming files
# -----------------------------------------------------------------------------
mv ${sample_id}.cnv ${sample_id}.quantisnp.cnv
mv ${sample_id}.loh ${sample_id}.quantisnp.loh
mv ${sample_id}.qc ${sample_id}.quantisnp.qc

echo "Merging PennCNV results..."
cat ${sample_id}.penncnv.chrx.out ${sample_id}.penncnv.out > ${sample_id}.penncnv.cnv

rm ${sample_id}.penncnv.chrx.out ${sample_id}.penncnv.out

rm ${sexfile_temp}

format_quantisnp_cnv.sh "${sample_id}.quantisnp.cnv" "${sample_id}.quantisnp.cnv.tsv"
format_penncnv_cnv.sh "${sample_id}.penncnv.cnv" "${sample_id}.penncnv.cnv.tsv"
