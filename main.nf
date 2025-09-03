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
       CNV-Caller ${cohort_tag} run summary:
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


process merge_sample_metadata {
    tag "merge sample level data"

    input:
    path plink2samplemetadata_tsv
    path penncnv_qc

    output:
    path "sampleDB.tsv"

    script:
    """
    duckdb -c "
    COPY (
        SELECT a.*, b.* EXCLUDE SampleID
        FROM read_csv_auto('${plink2samplemetadata_tsv}', sep='\t', header=true) a
        JOIN read_csv_auto('${penncnv_qc}', sep='\t', header=true) b
        USING (SampleID)
    ) TO 'sampleDB.tsv' (HEADER, DELIMITER '\t');
    "
    """
}



workflow {
    
    main:
    // Inputs from params
    penncnv_qc_path    = params.penncnv_qc_path
    penncnv_calls_path    = params.penncnv_calls_path
    quantisnp_calls_path  = params.quantisnp_calls_path
    list_sample_baflrrpath = params.list_sample_baflrrpath
    list_baflrr_path = params.list_baflrr_path
    plink2samplemetadata_tsv = params.plink2samplemetadata_tsv
    genome_version = params.genome_version
    batch_size = params.batch_size
    dataset_name = params.dataset_name
    
    
    if (plink2samplemetadata_tsv) {

        sex_file = sexfile_from_plink2samplemetadata(plink2samplemetadata_tsv)
        trio_file = trio_from_plink2samplemetadata(plink2samplemetadata_tsv)

    }
    
    has_input_calls = penncnv_calls_path && quantisnp_calls_path

    if ( ! has_input_calls ) {
        '''
        PREPARE INPUTS for PennCNV
        '''
        // PREPARE_PENNCNV_INPUTS ( first_sample.getParent(),
        PREPARE_PENNCNV_INPUTS ( list_sample_baflrrpath,
                                plink2samplemetadata_tsv,
                                file("${projectDir}/resources/GC_correction/${genome_version}/gc_content_1k_windows.bed"))

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
        penncnv_cnv_raw   = Channel.fromPath(penncnv_calls_path)
        quantisnp_cnv_raw = Channel.fromPath(quantisnp_calls_path)
        penncnv_qc = Channel.fromPath(penncnv_qc_path)
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


    merge_sample_metadata ( plink2samplemetadata_tsv , penncnv_qc )


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
    sampleDB = merge_sample_metadata.out

    // Before filter results
    penncnv_qc       = !has_input_calls ? CALL_CNV_PARALLEL.out.penncnv_qc_ch : Channel.empty()
    penncnv_cnv      = !has_input_calls ? CALL_CNV_PARALLEL.out.penncnv_cnv_ch : Channel.empty()
    penncnv_cnv_raw  = !has_input_calls ? CALL_CNV_PARALLEL.out.penncnv_cnv_raw_ch : penncnv_cnv_raw
    quantisnp_cnv    = !has_input_calls ? CALL_CNV_PARALLEL.out.quantisnp_cnv_ch : Channel.empty()
    quantisnp_cnv_raw= !has_input_calls ? CALL_CNV_PARALLEL.out.quantisnp_cnv_raw_ch : quantisnp_cnv_raw

    // REPORT outputs only if report was run
    penncnv_unfilter_cnv_qc   = params.report ? REPORT_PDF.out.penncnv_unfilter_cnv_qc : Channel.empty()
    quantisnp_unfilter_cnv_qc = params.report ? REPORT_PDF.out.quantisnp_unfilter_cnv_qc : Channel.empty()
    merged_cnv_qc             = params.report ? REPORT_PDF.out.merged_cnv_qc : Channel.empty()
    sample_qc_report          = params.report ? REPORT_PDF.out.sample_qc_report : Channel.empty()

    report_summary = buildSummary.out
}


output {
    merged_cnv {
        mode 'copy'
        path "${params.dataset_name}/"
    }
    sampleDB {
        mode 'copy'
        path "${params.dataset_name}/"
    }
    report_summary {
        mode 'copy'
        path "${params.dataset_name}/docs/"
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