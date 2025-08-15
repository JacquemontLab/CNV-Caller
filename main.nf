#!/usr/bin/env nextflow

nextflow.enable.dsl=2
nextflow.preview.output = true
nextflow.enable.moduleBinaries = true

//import subworkflows
include { PREPARE_PENNCNV_INPUTS } from './modules/generate_penncnv_params'
include { CALL_CNV_PARALLEL      } from './modules/CNV_calling_per_sampleID'
include { MERGE_CNV_CALLS        } from './modules/merge_dataset_CNV'
include { REPORT_PDF       } from './modules/qc_report_pdf'


process buildSummary {
    tag 'quick'
    
    input:
    val cohort_tag
    val genome_version
    path last_outfile

    output:
    path "launch_report.txt"

    script:
    """
        # Convert workflow start datetime to epoch seconds
        start_sec=\$(date -d "${workflow.start}" +%s)
        # Get current time in epoch seconds
        end_sec=\$(date +%s)

        # Calculate duration in seconds
        duration=\$(( end_sec - start_sec ))

        # Convert duration to minutes and seconds
        minutes=\$(( duration / 60 ))
        seconds=\$(( duration % 60 ))

       cat <<EOF > launch_report.txt
       Calling Pipeline on ${cohort_tag} run summary:
       run name: ${workflow.runName}
       version: ${workflow.manifest.version}
       configs: ${workflow.configFiles}
       workDir: ${workflow.workDir}
       genome_version: ${genome_version}
       launch_user: ${workflow.userName}
       start_time: ${workflow.start}
       duration: \${minutes} minutes and \${seconds} seconds

       Command:
       ${workflow.commandLine}

    
    """

    stub:
    """
    touch launch_report.txt
    """
}


process sexfile_from_plink2samplemetadata {
    tag "generate sexfile from plink2samplemetadata"

    input:
    path plink2samplemetadata_tsv

    output:
    path 'sexfile.tsv'

    script:
    """
    cut -f1,3 ${plink2samplemetadata_tsv} > 'sexfile.tsv'
    """
}


process trio_from_plink2samplemetadata {
    tag "generate trio file from plink2samplemetadata"

    input:
    path plink2samplemetadata_tsv

    output:
    path 'trio.tsv'

    script:
    """
    cut -f1,4,5 ${plink2samplemetadata_tsv} > 'trio.tsv'
    """
}



workflow {
    main:
    // Inputs from params
    penncnv_calls_path    = params.penncnv_calls_path
    quantisnp_calls_path  = params.quantisnp_calls_path
    list_sample_baflrrpath = params.list_sample_baflrrpath
    list_baflrr_path = params.list_baflrr_path
    plink2samplemetadata_tsv = params.plink2samplemetadata_tsv
    gc_correction_dir = params.gc_correction_dir
    genome_version = params.genome_version
    batch_size = params.batch_size
    dataset_name = params.dataset_name
    
    
    has_input_calls = penncnv_calls_path && quantisnp_calls_path

    if ( ! has_input_calls ) {
        '''
        PREPARE INPUTS for PennCNV
        '''
        // PREPARE_PENNCNV_INPUTS ( first_sample.getParent(),
        PREPARE_PENNCNV_INPUTS ( list_sample_baflrrpath,
                                plink2samplemetadata_tsv,
                                gc_correction_dir,
                                genome_version )


        sex_file = sexfile_from_plink2samplemetadata(plink2samplemetadata_tsv)
        trio_file = trio_from_plink2samplemetadata(plink2samplemetadata_tsv)

        '''
        CALLING CNVs AND MERGE
        '''
        CALL_CNV_PARALLEL    ( Channel.fromPath(list_baflrr_path),
                            PREPARE_PENNCNV_INPUTS.out.pfb_file, 
                            PREPARE_PENNCNV_INPUTS.out.gc_model,
                            sex_file,
                            genome_version,
                            batch_size )
        
        penncnv_cnv_raw = CALL_CNV_PARALLEL.out.penncnv_cnv_raw_ch
        quantisnp_cnv_raw = CALL_CNV_PARALLEL.out.quantisnp_cnv_raw_ch

        penncnv_qc = CALL_CNV_PARALLEL.out.penncnv_qc_ch

    } else {
        penncnv_cnv_raw = penncnv_calls_path
        quantisnp_cnv_raw = quantisnp_calls_path
    }

    '''
    PREFILTERING - MERGING
    '''
    MERGE_CNV_CALLS (
                    quantisnp_cnv_raw,
                    penncnv_cnv_raw,
                    file("${projectDir}/resources/Genome_Regions/Genome_Regions_data.tsv"),
                    genome_version
                     )

    merged_cnv = MERGE_CNV_CALLS.out.merged_cnv_ch

    // Run report only if requested
    if (params.report) {
        REPORT_PDF ( 
            dataset_name,
            plink2samplemetadata_tsv,
            trio_file,
            penncnv_qc,
            penncnv_cnv_raw,
            quantisnp_cnv_raw,
            merged_cnv) 
    }

    '''
    SUMMARY BUILDER
    '''
    // Always run summary builder (change if you want this conditional too)
    buildSummary  ( dataset_name,
                    genome_version,
                    params.report ? REPORT_PDF.out.merged_cnv_qc : merged_cnv )


    publish:

    // MERGED CNV DATASET
    merged_cnv = MERGE_CNV_CALLS.out.merged_cnv_ch

    // Before filter results
    penncnv_qc       = !has_input_calls ? CALL_CNV_PARALLEL.out.penncnv_qc_ch : null
    penncnv_cnv      = !has_input_calls ? CALL_CNV_PARALLEL.out.penncnv_cnv_ch : null
    penncnv_cnv_raw  = !has_input_calls ? CALL_CNV_PARALLEL.out.penncnv_cnv_raw_ch : penncnv_cnv_raw
    quantisnp_cnv    = !has_input_calls ? CALL_CNV_PARALLEL.out.quantisnp_cnv_ch : null
    quantisnp_cnv_raw= !has_input_calls ? CALL_CNV_PARALLEL.out.quantisnp_cnv_raw_ch : quantisnp_cnv_raw

    // REPORT outputs only if report was run
    penncnv_unfilter_cnv_qc   = params.report ? REPORT_PDF.out.penncnv_unfilter_cnv_qc : null
    quantisnp_unfilter_cnv_qc = params.report ? REPORT_PDF.out.quantisnp_unfilter_cnv_qc : null
    merged_cnv_qc             = params.report ? REPORT_PDF.out.merged_cnv_qc : null
    sample_qc_report          = params.report ? REPORT_PDF.out.sample_qc_report : null

}


output {
    merged_cnv {
        mode 'copy'
        path "${params.dataset_name}/"
    }
    penncnv_qc {
        mode 'copy'
        path "${params.dataset_name}/calls_unfiltered/penncnv/"
    }
    penncnv_cnv {
        mode 'copy'
        path "${params.dataset_name}/calls_unfiltered/penncnv/"
    }
    penncnv_cnv_raw {
        mode 'copy'
        path "${params.dataset_name}/calls_unfiltered/penncnv/"
    }
    quantisnp_cnv {
        mode 'copy'
        path "${params.dataset_name}/calls_unfiltered/quantisnp/"
    }

    quantisnp_cnv_raw {
        mode 'copy'
        path "${params.dataset_name}/calls_unfiltered/quantisnp/"
    }
    penncnv_unfilter_cnv_qc {
        mode 'copy'
        path "${params.dataset_name}/docs/"
    }

    quantisnp_unfilter_cnv_qc {
        mode 'copy'
        path "${params.dataset_name}/docs/"
    }

    merged_cnv_qc {
        mode 'copy'
        path "${params.dataset_name}/docs/"
    }

    sample_qc_report {
        mode 'copy'
        path "${params.dataset_name}/docs/"
    }

}