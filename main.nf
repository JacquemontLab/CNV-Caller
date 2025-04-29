#!/usr/bin/env nextflow

nextflow.enable.dsl=2
// NXF_OFFLINE=true nextflow run main.nf -resume -c nextflow.config -with-trace -with-report


include { PLINK_EXTRACTED_DATA } from './modules/extract_data_from_plink'
include { PENNCNV_PARAMS } from './modules/generate_penncnv_params'
include { CNV_CALLING } from './modules/CNV_calling_per_sampleID'


params.genome_version = "GRCh37"

params.gcDir = "/lustre06/project/6008022/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/generate_penncnv_params/data/" + params.genome_version + "_GCdir/"
params.gc_content_windows = "/lustre06/project/6008022/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/generate_penncnv_params/data/gc_content_1k_windows_" + params.genome_version + ".bed"
params.regions_file = "/lustre06/project/6008022/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/CNV_calling_per_sampleID/data/Genome_Regions_data.tsv"


workflow {
    plink_base_path = "/lustre06/project/6008022/All_user_common_folder/RAW_DATA/Genetic/ALSPAC/PLINK/ALSPAC"

    def directory_BAF_LRR_Probes_by_sample = "/lustre06/project/6008022/All_user_common_folder/RAW_DATA/Genetic/ALSPAC/BAF_LRR_Probes_by_sample/"
    // def selected_sample_ids = ["ALSPAC09897249"]
    def selected_sample_ids = ["ALSPAC09897249", "ALSPAC09902309"]
    
    // Create a channel from the files in the directory
    Channel
        .fromPath("${directory_BAF_LRR_Probes_by_sample}/*")
        .filter { file ->
            def sampleID = file.getName().split("\\.")[0]
            sampleID in selected_sample_ids
        }
        .map { file ->
            def sampleID = file.getName().split("\\.")[0]
            [sampleID, file]
        }
        .set { sample_inputs }

    // Extract the first file from the channel
    sample_inputs
        .map { sampleID, file -> file }
        .first()
        .set { probes_file }


    PLINK_EXTRACTED_DATA(
        plink_base_path,
        probes_file
    )

    PENNCNV_PARAMS(
        directory_BAF_LRR_Probes_by_sample,
        PLINK_EXTRACTED_DATA.out,
        file(params.gc_content_windows)
    )


    sample_inputs
        .combine(PENNCNV_PARAMS.out.pfb_file)
        .combine(PENNCNV_PARAMS.out.gc_model)
        .combine(PLINK_EXTRACTED_DATA.out.from_plink_extracted_data)
        .map { it1, gc_model, plink_data ->
            def (sampleID, file) = it1[0]  // sample_inputs tuple
            def pfb_file = it1[1]          // pfb_file from first combine
            tuple(sampleID, file, pfb_file, gc_model, plink_data, file(params.gcDir)) // The tuple of interest
        }
        .set { cnv_inputs }

    cnv_inputs.view()


    CNV_CALLING(
        cnv_inputs
    )

}

