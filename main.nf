#!/usr/bin/env nextflow

nextflow.enable.dsl=2
// NXF_OFFLINE=true nextflow run main.nf -resume -c nextflow.config -with-trace -with-report



//Part 1 extract_data_from_plink
include { PLINK_EXTRACTED_DATA } from './modules/extract_data_from_plink'
include { PENNCNV_PARAMS } from './modules/generate_penncnv_params'


params.genome_version = "GRCh37"
params.gcDir = "/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/penncnv_params/resources/" + params.genome_version + "_GCdir/"
params.gc_content_windows = "/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/penncnv_params/resources/gc_content_1k_windows_" + params.genome_version + ".bed"


workflow {
    plink_base_path = "/home/flben/projects/rrg-jacquese/All_user_common_folder/RAW_DATA/Genetic/ALSPAC/PLINK/ALSPAC"
    // directory_BAF_LRR_Probes_by_sample = file("/home/flben/projects/rrg-jacquese/All_user_common_folder/RAW_DATA/Genetic/ALSPAC/BAF_LRR_Probes_by_sample/")

    probes_file = Channel
        .fromPath("/home/flben/projects/rrg-jacquese/All_user_common_folder/RAW_DATA/Genetic/ALSPAC/BAF_LRR_Probes_by_sample/ALSPAC09900260.BAF_LRR_Probes.tsv")


    // from_plink_extracted_data = 
    PLINK_EXTRACTED_DATA(
        plink_base_path,
        probes_file
    )
// .out.from_plink_extracted_data

    // PENNCNV_PARAMS(
    //     directory_BAF_LRR_Probes_by_sample,
    //     from_plink_extracted_data,
    //     file(params.gc_content_windows)
    // )
}

