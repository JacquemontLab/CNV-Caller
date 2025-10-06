
# Ensure Nextflow runs offline (no internet check)
export NXF_OFFLINE=false
export NXF_CONDA_CACHE=/home/jupyter/.conda_envs
mkdir -p $NXF_CONDA_CACHE


penncnv_calls_path=~/data_raw/merge_raw/PennCNV_raw_calls.txt
quantisnp_calls_path=~/data_raw/merge_raw/QuantiSNP_raw_calls.txt
penncnv_qc_path=~/data_raw/PennCNV_QC.tsv
plink2samplemetadata_tsv=~/plink/Plink2SampleMetadata/results/sample_metadata_from_plink.tsv


nextflow run ~/CNV-Caller/main.nf \
--pipeline_mode pipeline_partial --dataset_name AllOfUs_tierv8_Array_GRCh38 --genome_version GRCh38 \
--penncnv_calls_path $penncnv_calls_path \
--quantisnp_calls_path $quantisnp_calls_path \
--penncnv_qc_path $penncnv_qc_path \
--plink2samplemetadata_tsv $plink2samplemetadata_tsv \
--report true -c ~/CNV-Caller/setup/allofus/allofus.config -profile standard -resume
