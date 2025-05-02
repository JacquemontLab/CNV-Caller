#!/usr/bin/env nextflow

nextflow.enable.dsl=2
// NXF_OFFLINE=true nextflow run main.nf -resume -c nextflow.config -with-trace -with-report

//import subworkflows
include { PLINK_EXTRACTED_DATA } from './modules/extract_data_from_plink'
include { PENNCNV_PARAMS } from './modules/generate_penncnv_params'

//import processes
include {call_CNV as CALL_CNVs } from './modules/CNV_calling_per_sampleID'
include { merge_cnv_callers_and_extract_qc as MERGE_CNV_CALLS } from './modules/CNV_calling_per_sampleID'

params.genome_version = "GRCh37"


params.plink_dir = file("test_data/plink")
params.probe_files = file("/lustre06/project/6008022/All_user_common_folder/RAW_DATA/Genetic/ALSPAC/BAF_LRR_Probes_by_sample/*tsv")

workflow {

    // Input definitions
    sample_ch =  Channel.fromPath(params.probe_files)
    plink_dir = file(params.plink_dir)


    // resource definitions
    resource_dir       =  file(projectDir.resolve("resources"))  //resolve is the prefered method for building paths. ProjectDir points to the main dir. 
    gcDir              =  file(resource_dir / params.genome_version / "GCdir" )
    gc_content_windows =  file(resource_dir / params.genome_version / "gc_content_1k_windows.bed") 
    genomic_regions    =  file(resource_dir / "Genome_Regions_data.tsv")


    
    '''
    PREPARE INPUTS for PennCNV
    '''
    // Extract the first file from the channel 
    first_sample = Channel.fromPath(params.probe_files).first() 

    PLINK_EXTRACTED_DATA ( plink_dir,
                           first_sample                 )

    PENNCNV_PARAMS       ( first_sample.parent,
                           PLINK_EXTRACTED_DATA.out,
                           gc_content_windows           )

    test_ch = sample_ch.randomSample( 3 ) //for testing subset


    '''
    CALLING CNVs AND MERGE
    '''
    CALL_CNVs            ( test_ch,
                           PENNCNV_PARAMS.out.pfb_file.first(), //pulling queue channel into value channel using first()
                           PENNCNV_PARAMS.out.gc_model.first(),
                           PLINK_EXTRACTED_DATA.out.first(),
                           gcDir                                )

    MERGE_CNV_CALLS      ( CALL_CNVs.out,
                           genomic_regions,
                           params.genome_version                )

    '''
    COLLECTING FILES FOR OUTPUT
    '''

    MERGE_CNV_CALLS.out.CNVs_tsv.collectFile(keepHeader : true, 
                                             storeDir : projectDir.resolve("results"), 
                                             name : "cohort_CNVs.tsv")

    MERGE_CNV_CALLS.out.PennCNV_QC_tsv.collectFile(keepHeader : true, 
                                             storeDir : projectDir.resolve("results"), 
                                             name : "cohort_QC.tsv")

}

