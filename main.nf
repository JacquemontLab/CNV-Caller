#!/usr/bin/env nextflow

nextflow.enable.dsl=2
// NXF_OFFLINE=true nextflow run main.nf -resume -c conf/ccdb.config -with-trace -with-report

//import subworkflows
include { COLLECT_PLINK_DATA } from './modules/extract_data_from_plink'
include { PREPARE_PENNCNV_INPUTS } from './modules/generate_penncnv_params'

//import processes
include {call_CNV as CALL_CNVs } from './modules/CNV_calling_per_sampleID'
include { merge_cnv_callers_and_extract_qc as MERGE_CNV_CALLS } from './modules/CNV_calling_per_sampleID'

params.genome_version = "GRCh37"
params.cohort_tag = "ALSPAC"
params.outDir = projectDir
params.plink_dir = file("test_data/plink")
params.probe_files = file("/lustre06/project/6008022/All_user_common_folder/RAW_DATA/Genetic/ALSPAC/BAF_LRR_Probes_by_sample/*tsv")

process buildSummary {
    publishDir params.outDir.resolve(params.cohort_tag)

    input:
    val cohort_tag 
    path cnvs // collected cnvs into one file 
    path qc   // collected CNV QC reports into one file

    output:
    path "${cohort_tag}_summary.txt"

    script:
    """
     cat <<EOF > ${cohort_tag}_summary.txt
     CNV_Caller ${cohort_tag} run summary:
       version: ${workflow.manifest.version}
       configs: ${workflow.configFiles}
       workDir: ${workflow.workDir}
       launch_user: ${workflow.userName}
       start_time: ${workflow.start}
       duration: ${workflow.duration}
    """

    stub:
    """ 
    touch summary.txt
    """
}

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

    COLLECT_PLINK_DATA     ( plink_dir,
                             first_sample                 )

    PREPARE_PENNCNV_INPUTS ( first_sample.parent,
                             COLLECT_PLINK_DATA.out,
                             gc_content_windows           )

    test_ch = sample_ch.take( 3 ) //for testing subset


    '''
    CALLING CNVs AND MERGE
    '''
    CALL_CNVs            ( test_ch,
                           PREPARE_PENNCNV_INPUTS.out.pfb_file.first(), //pulling queue channel into value channel using first()
                           PREPARE_PENNCNV_INPUTS.out.gc_model.first(),
                           COLLECT_PLINK_DATA.out.first(),
                           gcDir                                         )

    MERGE_CNV_CALLS      ( CALL_CNVs.out,
                           genomic_regions,
                           params.genome_version                         )

    '''
    COLLECTING FILES FOR OUTPUT
    '''

    collect_cnv = MERGE_CNV_CALLS.out.CNVs_tsv.collectFile(keepHeader : true, 
                                             storeDir : file(params.outDir / params.cohort_tag / "results"),
                                             name : "CNVs.tsv")
                                             

    collect_qc = MERGE_CNV_CALLS.out.PennCNV_QC_tsv.collectFile(keepHeader : true, 
                                                   storeDir : file(params.outDir / params.cohort_tag / "results"), 
                                                   name : "pennCNV_QC.tsv")
    '''
    SUMMARY BUILDER
    '''
                                           
    buildSummary ( params.cohort_tag, collect_cnv, collect_qc )

}

