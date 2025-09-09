#!/bin/bash
#SBATCH --job-name=CNV_caller_light          # Job name
#SBATCH --ntasks=1                        # Single task
#SBATCH --cpus-per-task=16                # Number of CPU cores per task
#SBATCH --mem-per-cpu=3500MB              # Memory per CPU
#SBATCH --time=01:00:00                   # Time limit (hh:mm:ss)
#SBATCH --output=CNV_caller_light_%j.log     # Standard output and error log
#SBATCH --account=rrg-jacquese            # Account name

cd $SLURM_SUBMIT_DIR || exit 1

# -------------------------------
# Usage:
# sbatch run_CNV_caller_light.sh \
#                          --git_dir /path/to/git_repo \
#                          --dataset_name UKBB_N488k \
#                          --genome_version GRCh37 --plink2samplemetadata_tsv file.tsv \
#                          --penncnv_calls_path file.txt --quantisnp_calls_path file.txt \
#                          --penncnv_qc_path file.tsv \
#                          [--report true|false]
# -------------------------------


# Fixed pipeline mode
pipeline_mode="pipeline_partial"

# Default report
report="false"

# Parse arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --git_dir) git_dir="$2"; shift 2 ;;
        --dataset_name) dataset_name="$2"; shift 2 ;;
        --genome_version) genome_version="$2"; shift 2 ;;
        --plink2samplemetadata_tsv) plink2samplemetadata_tsv="$2"; shift 2 ;;
        --penncnv_calls_path) penncnv_calls_path="$2"; shift 2 ;;
        --quantisnp_calls_path) quantisnp_calls_path="$2"; shift 2 ;;
        --penncnv_qc_path) penncnv_qc_path="$2"; shift 2 ;;
        --report) report="$2"; shift 2 ;;
        *) echo "Unknown option $1"; exit 1 ;;
    esac
done

# Check required arguments
for var in git_dir dataset_name genome_version plink2samplemetadata_tsv penncnv_calls_path quantisnp_calls_path penncnv_qc_path; do
    if [[ -z "${!var}" ]]; then
        echo "Error: --$var is required"
        exit 1
    fi
done

# Ensure Nextflow runs offline
export NXF_OFFLINE=true

# Load modules
module load nextflow
module load r

# Resolve main.nf and config based on git_dir
MAIN_NF="${git_dir}/main.nf"
CONFIG_FILE="${git_dir}/setup/ccdb/ccdb.config"

# Run Nextflow
nextflow run "$MAIN_NF" \
    --pipeline_mode "$pipeline_mode" \
    --dataset_name "$dataset_name" \
    --genome_version "$genome_version" \
    --penncnv_calls_path "$penncnv_calls_path" \
    --quantisnp_calls_path "$quantisnp_calls_path" \
    --penncnv_qc_path "$penncnv_qc_path" \
    --plink2samplemetadata_tsv "$plink2samplemetadata_tsv" \
    --report "$report" \
    -c "$CONFIG_FILE" \
    -profile standard \
    -resume

