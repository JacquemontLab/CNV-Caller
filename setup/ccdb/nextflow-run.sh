#!/bin/bash
#SBATCH --job-name=CNV_caller          # Job name
#SBATCH --ntasks=1                        # Single task
#SBATCH --cpus-per-task=16                # Number of CPU cores per task
#SBATCH --mem-per-cpu=3500MB              # Memory per CPU
#SBATCH --time=05:00:00                   # Time limit (hh:mm:ss)
#SBATCH --output=CNV_caller_%j.log     # Standard output and error log
#SBATCH --account=rrg-jacquese            # Account name

dataset_name="UKBB_N488k"
genome_version="GRCh37"

list_sample_baflrrpath=/home/flben/list.tsv

list_baflrr_path=$(tail -n +2 "$list_sample_baflrrpath" | cut -f2)

# Ensure Nextflow runs offline (no internet check)
export NXF_OFFLINE=true

module load nextflow
module load r

nextflow run main.nf \
    --dataset_name "$dataset_name" \
    --list_sample_baflrrpath "$list_sample_baflrrpath" \
    --list_baflrr_path "$list_baflrr_path" \
    --plink2samplemetadata_output "$plink2samplemetadata_output" \
    --genome_version "$genome_version" \
    --batch_size 10 \
    -c setup/ccdb/ccdb.config \
    -profile standard \
    -resume


