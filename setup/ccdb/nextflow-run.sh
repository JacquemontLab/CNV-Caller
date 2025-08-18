#!/bin/bash
#SBATCH --job-name=CNV_caller          # Job name
#SBATCH --ntasks=1                        # Single task
#SBATCH --cpus-per-task=16                # Number of CPU cores per task
#SBATCH --mem-per-cpu=3500MB              # Memory per CPU
#SBATCH --time=05:00:00                   # Time limit (hh:mm:ss)
#SBATCH --output=CNV_caller_%j.log     # Standard output and error log
#SBATCH --account=rrg-jacquese            # Account name

GC_correction_dir=/home/flben/links/projects/rrg-jacquese/flben/WES_CNV_calls/scripts/CNV-Caller/resources/GC_correction/
dataset_name="UKBB_N488k"
genome_version="GRCh37"

penncnv_calls_path=/home/flben/links//projects/rrg-jacquese/LAB_WORKSPACE/RAW_DATA/Genetic/UKBB/PennCNV_raw_calls.txt
quantisnp_calls_path=/home/flben/links//projects/rrg-jacquese/LAB_WORKSPACE/RAW_DATA/Genetic/UKBB/QuantiSNP_raw_calls.txt

# Ensure Nextflow runs offline (no internet check)
export NXF_OFFLINE=true

module load nextflow
module load r

nextflow run main.nf \
    --dataset_name "$dataset_name" \
    --penncnv_calls_path "$penncnv_calls_path" \
    --quantisnp_calls_path "$quantisnp_calls_path" \
    --gc_correction_dir "$GC_correction_dir" \
    --genome_version "$genome_version" \
    --report "true" \
    -c setup/ccdb/ccdb.config \
    -profile standard \
    -resume



    # penncnv_calls_path    = params.penncnv_calls_path
    # quantisnp_calls_path  = params.quantisnp_calls_path
    # list_sample_baflrrpath = params.list_sample_baflrrpath
    # list_baflrr_path = params.list_baflrr_path
    # plink2samplemetadata_tsv = params.plink2samplemetadata_tsv
    # gc_correction_dir = params.gc_correction_dir
    # genome_version = params.genome_version
    # batch_size = params.batch_size
    # dataset_name = params.dataset_name
    # -with-dag flowchart.pdf \
    # -with-timeline timeline.html \
    # -with-report report.html \
    # -with-trace trace.txt \