#!/usr/bin/env nextflow

nextflow.enable.dsl=2
// NXF_OFFLINE=true nextflow run main.nf -resume -c nextflow.config -with-trace -with-report



//Part 1 extract_data_from_plink
include { PLINK_EXTRACTED_DATA } from './modules/extract_data_from_plink'
include { PENNCNV_PARAMS } from './modules/generate_penncnv_params'
include { CNV_CALLING } from './modules/CNV_calling_per_sampleID'


params.genome_version = "GRCh37"
params.gcDir = "/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/penncnv_params/data/" + params.genome_version + "_GCdir/"
params.gc_content_windows = "/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/penncnv_params/data/gc_content_1k_windows_" + params.genome_version + ".bed"
params.regions_file = "/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/CNV_calling_per_sampleID/data/Genome_Regions_data.tsv"


workflow {
    plink_base_path = "/home/flben/projects/rrg-jacquese/All_user_common_folder/RAW_DATA/Genetic/ALSPAC/PLINK/ALSPAC"

    def directory_BAF_LRR_Probes_by_sample = "/home/flben/projects/rrg-jacquese/All_user_common_folder/RAW_DATA/Genetic/ALSPAC/BAF_LRR_Probes_by_sample/"
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

    sample_inputs.view()

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

    CNV_CALLING(
        sample_inputs,
        PENNCNV_PARAMS.out.pfb_file,
        PENNCNV_PARAMS.out.gc_model,
        PLINK_EXTRACTED_DATA.out,
        file(params.gcDir),
        params.genome_version,
        file(params.regions_file)
    )

}

