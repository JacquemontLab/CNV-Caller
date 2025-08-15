#!/bin/bash
#SBATCH --job-name=CNV_caller          # Job name
#SBATCH --ntasks=1                        # Single task
#SBATCH --cpus-per-task=16                # Number of CPU cores per task
#SBATCH --mem-per-cpu=3500MB              # Memory per CPU
#SBATCH --time=05:00:00                   # Time limit (hh:mm:ss)
#SBATCH --output=CNV_caller_%j.log     # Standard output and error log
#SBATCH --account=rrg-jacquese            # Account name


list_sample_baflrrpath=/home/flben/links/projects/rrg-jacquese/flben/WES_CNV_calls/data/baf_lrr_filter_path.tsv
plink2samplemetadata_tsv=/home/flben/links/projects/rrg-jacquese/flben/Plink2SampleMetadata/results/sample_metadata_from_plink.tsv
GC_correction_dir=/home/flben/links/projects/rrg-jacquese/flben/WES_CNV_calls/scripts/CNV-Caller/resources/GC_correction/
batch_size=384
dataset_name="SPARK_WES_N2533"
genome_version="GRCh38"




# Ensure Nextflow runs offline (no internet check)
export NXF_OFFLINE=true

module load nextflow
module load r

cut -f2 $list_sample_baflrrpath > list_baflrr_path.txt

nextflow run main.nf \
    --dataset_name "$dataset_name" \
    --list_sample_baflrrpath "$list_sample_baflrrpath" \
    --list_baflrr_path "list_baflrr_path.txt" \
    --plink2samplemetadata_tsv "$plink2samplemetadata_tsv" \
    --batch_size "$batch_size" \
    --gc_correction_dir "$GC_correction_dir" \
    --genome_version "$genome_version" \
    -c setup/ccdb/ccdb.config \
    -profile standard \
    -resume



    # -with-dag flowchart.pdf \
    # -with-timeline timeline.html \
    # -with-report report.html \
    # -with-trace trace.txt \