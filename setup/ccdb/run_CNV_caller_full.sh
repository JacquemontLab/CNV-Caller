#!/bin/bash
#SBATCH --job-name=CNV_caller_full          # Job name
#SBATCH --ntasks=1                        # Single task
#SBATCH --cpus-per-task=16                # Number of CPU cores per task
#SBATCH --mem-per-cpu=3500MB              # Memory per CPU
#SBATCH --time=05:00:00                   # Time limit (hh:mm:ss)
#SBATCH --output=CNV_caller_full_%j.log     # Standard output and error log
#SBATCH --account=rrg-jacquese            # Account name


cd $SLURM_SUBMIT_DIR || exit 1

# ------------------------------------------------------------------------------
# Script: run_CNV_caller_full.sh
#
# Description:
#   SLURM wrapper for running CNV calling (BAF/LRR mode) with Nextflow.
#
# Usage:
#   sbatch run_CNV_caller_full.sh \
#          --git_dir /path/to/git_repo \
#          --dataset_name UKBB_N488k \
#          --genome_version GRCh37 \
#          --plink2samplemetadata_tsv file.tsv \
#          --list_sample_baflrrpath list.tsv \
#          [--batch_size 20] \
#          [--report true|false]
# ------------------------------------------------------------------------------

# Default values
pipeline_mode="pipeline_full"
report="false"
batch_size=10

# -------------------------------
# Parse arguments
# -------------------------------
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --git_dir) git_dir="$2"; shift 2 ;;
        --dataset_name) dataset_name="$2"; shift 2 ;;
        --genome_version) genome_version="$2"; shift 2 ;;
        --plink2samplemetadata_tsv) plink2samplemetadata_tsv="$2"; shift 2 ;;
        --list_sample_baflrrpath) list_sample_baflrrpath="$2"; shift 2 ;;
        --batch_size) batch_size="$2"; shift 2 ;;
        --report) report="$2"; shift 2 ;;
        *) echo "Unknown option $1"; exit 1 ;;
    esac
done


# -------------------------------
# Check required arguments
# -------------------------------
for var in git_dir dataset_name genome_version plink2samplemetadata_tsv list_sample_baflrrpath; do
    if [[ -z "${!var}" ]]; then
        echo "Error: --$var is required"
        exit 1
    fi
done


list_baflrr_path=$(tail -n +2 "$list_sample_baflrrpath" | cut -f2)

# Ensure Nextflow runs offline
export NXF_OFFLINE=true

# Load modules
module load nextflow
module load r

# Resolve main.nf and config based on git_dir
MAIN_NF="${git_dir}/main.nf"
CONFIG_FILE="${git_dir}/setup/ccdb/ccdb.config"

nextflow run "$MAIN_NF" \
    --dataset_name "$dataset_name" \
    --genome_version "$genome_version" \
    --list_sample_baflrrpath "$list_sample_baflrrpath" \
    --list_baflrr_path "$list_baflrr_path" \
    --plink2samplemetadata_tsv "$plink2samplemetadata_tsv" \
    --batch_size $batch_size \
    --report "$report" \
    -c "$CONFIG_FILE" \
    -profile standard \
    -resume


