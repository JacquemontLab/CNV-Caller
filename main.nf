#!/usr/bin/env nextflow

nextflow.enable.dsl=2


//import subworkflows
include { PREPARE_PENNCNV_INPUTS } from './modules/prepare_penncnv_inputs'
include { CALL_CNV               } from './modules/call_cnv'
include { MERGE_CNV_CALLS        } from './modules/merge_cnv_calls'
include { MAKE_REPORT            } from './modules/make_report'


process buildSummary {
    
    input:
    val input_file
    val git_hash
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

    # Convert duration to hours, minutes, seconds
    hours=\$(( duration / 3600 ))
    minutes=\$(( (duration % 3600) / 60 ))
    seconds=\$(( duration % 60 ))


    cat <<EOF > launch_report.txt
    CNV_Caller ${cohort_tag} run summary:
    run name: ${workflow.runName}
    version: ${workflow.manifest.version}
    configs: ${workflow.configFiles}
    workDir: ${workflow.workDir}
    input_file: ${input_file}
    genome_version: ${genome_version}
    launch_user: ${workflow.userName}
    start_time: ${workflow.start}
    duration: \${hours}h \${minutes}m \${seconds}s

    Command:
    ${workflow.commandLine}

    Git hash working version:
    commit ${git_hash}
    """

    stub:
    """
    touch launch_report.txt
    """
}


process mergeSampleMetadata {
    label "duckdb"

    input:
    path plink2samplemetadata
    path penncnv_qc

    output:
    path "sampleDB.tsv"

    script:
    """
    duckdb -c "
    COPY (
        SELECT
             CAST(meta.SampleID AS VARCHAR) AS SampleID,
             meta.* EXCLUDE SampleID,
             pennqc.* EXCLUDE SampleID
         FROM read_csv_auto('$penncnv_qc', sep='\t', header=true) AS pennqc
         LEFT JOIN read_csv_auto('$plink2samplemetadata', sep='\t', header=true) AS meta
         ON CAST(pennqc.SampleID AS VARCHAR) = CAST(meta.SampleID AS VARCHAR)

    )  TO 'sampleDB.tsv' (HEADER, DELIMITER '\t');
    "
    """
}

process formatRawCNV {

    input:
    path penncnv_cnv_raw
    path quantisnp_cnv_raw

    output:
    path "PennCNV_CNV.tsv",     emit: penncnv_cnv
    path penncnv_cnv_raw,       emit: penncnv_cnv_raw
    path "QuantiSNP_CNV.tsv",   emit: quantisnp_cnv
    path quantisnp_cnv_raw,     emit: quantisnp_cnv_raw

    script:
    """
    # PennCNV formatting:

    # Remove quotes from input files
    sed 's/"//g' "$penncnv_cnv_raw" > pc_no_quotes.txt
    format_penncnv_cnv.sh "pc_no_quotes.txt" "PennCNV_CNV.tsv"

    # Quantisnp Formatting:
    sed 's/"//g' "$quantisnp_cnv_raw" > qs_no_quotes.txt
    format_quantisnp_cnv.sh "qs_no_quotes.txt" "QuantiSNP_CNV.tsv"
    """
}


process extractPlink {

    input:
    path plink2samplemetadata

    output:
    path 'sexfile.tsv', emit: sexfile
    path 'trio.tsv', emit: triofile

    script:
    """
    duckdb -c "
    COPY (
        SELECT SampleID, Sex
        FROM read_csv_auto('$plink2samplemetadata', sep='\t', header=true)
    ) TO 'sexfile.tsv' (HEADER, DELIMITER '\\t');

    COPY (
        SELECT SampleID, FatherID, MotherID
        FROM read_csv_auto('$plink2samplemetadata', sep='\t', header=true)
    ) TO 'trio.tsv' (HEADER, DELIMITER '\\t');
    "
    """
}



include { STAGE_FILE as stageFile1 } from './modules'
include { STAGE_FILE as stageFile2 } from './modules'
include { STAGE_FILE as stageFile3 } from './modules'
include { STAGE_FILE as stageFile4 } from './modules'
include { STAGE_FILE as stageFile5 } from './modules'



// Default params
params.pipeline_mode            = "full"
params.plink2samplemetadata = ""
params.report                   = false
params.genome_version           = "GRCh38"
params.dataset_name             = "dataset"
params.git_hash                 = "git -C ${projectDir} rev-parse HEAD".execute().text.trim()
params.skip_hmm_training        = false

// Full specific
params.sample_file              = ""
params.pfb_max_sample_size      = 1000
params.data_type                = "array"
params.batch_size               = 64
params.batch_num                = -1 // for tuning batch sizes: default -1 means take all batches. 
                                     // Any other number restricts the execution to N number of batches  

// Partial specific
params.penncnv_qc_path          = ""
params.penncnv_calls_path       = ""
params.quantisnp_calls_path     = ""


workflow {
    
    main:
    
    plink_ch = params.plink2samplemetadata ?
            extractPlink(params.plink2samplemetadata)
            : channel.empty()

    if (params.pipeline_mode == "full"){

        log.info("Running Full pipeline on ${params.dataset_name}")

        sample_file_ch  = channel.fromPath(params.sample_file)

        batch_ch = sample_file_ch.splitCsv(sep: "\t",  header: ['sampleID', 'path_to_BAF_LRR'])                       
                                 .map {row -> row['path_to_BAF_LRR']}                              // grab filepaths only
                                 .buffer( size : params.batch_size, remainder : true)              // split channel into batches of nextflow lists
                                 .take(params.batch_num)                                           // for tuning batch sizes: default -1 means take all batches
               
        
        
       
        '''
        PREPARE INPUTS FOR PENNCNV

        Will create: 
            - HMM File for given tech
            - PFB File 
            - GC Model 
        '''

        PREPARE_PENNCNV_INPUTS ( sample_file_ch, 
                                 params.plink2samplemetadata,
                                 file("${projectDir}/resources/GC_correction/${params.genome_version}/gc_content_1k_windows.bed"),
                                 params.pfb_max_sample_size,
                                 params.data_type,
                                 params.skip_hmm_training                                                                                        
                                )

        '''
        CALLING CNVs AND MERGE
        '''
        CALL_CNV              ( batch_ch,                                                               // File of paths to baf_lrr files without the sampleID
                                PREPARE_PENNCNV_INPUTS.out.pfb_file,                                    // PFB file
                                PREPARE_PENNCNV_INPUTS.out.hmm_file,                                    // HMM file
                                PREPARE_PENNCNV_INPUTS.out.gc_model,                                    // GC model
                                plink_ch.sexfile,                                                       // Sexfile from metadata input 
                                params.genome_version                                                   // genome version for choosing gc content directory                
                              )
        
        // Collect outputs
        penncnv_cnv_raw     = CALL_CNV.out.penncnv_cnv_raw_ch
                                                    .flatten()
                                                    .collectFile(keepHeader : true,
                                                                 name       :"PennCNV_raw_calls.txt")

        quantisnp_cnv_raw   = CALL_CNV.out.quantisnp_cnv_raw_ch
                                                   .flatten()
                                                   .collectFile(keepHeader : true,
                                                                name       :"QuantiSNP_raw_calls.txt")

        penncnv_cnv         = CALL_CNV.out.penncnv_cnv_ch                                                   
                                                   .flatten()
                                                   .collectFile(keepHeader : true,
                                                                name       :"PennCNV_CNV.tsv")
        quantisnp_cnv       = CALL_CNV.out.quantisnp_cnv_ch                                                   
                                                   .flatten()
                                                   .collectFile(keepHeader : true,
                                                                name       :"QuantiSNP_CNV.tsv")

        penncnv_qc          = CALL_CNV.out.penncnv_qc_ch
                                                   .flatten()
                                                   .collectFile(keepHeader : true,
                                                                name       :"PennCNV_QC.tsv")

        sample_db_ch = mergeSampleMetadata(
                    params.plink2samplemetadata,
                    penncnv_qc
                    )
    
    } else if (params.pipeline_mode == "partial") {

        log.info("Running Partial pipeline on ${params.dataset_name}")

        formatRawCNV(params.penncnv_calls_path, params.quantisnp_calls_path)

        penncnv_cnv    = formatRawCNV.out.penncnv_cnv
        penncnv_cnv_raw   = stageFile1(params.penncnv_calls_path, "PennCNV_raw_calls.txt")

        quantisnp_cnv  = formatRawCNV.out.quantisnp_cnv
        quantisnp_cnv_raw = stageFile2(params.quantisnp_calls_path, "QuantiSNP_raw_calls.txt")

        
        penncnv_qc = params.penncnv_qc_path \
            ? stageFile3(params.penncnv_qc_path, "PennCNV_QC.tsv")
            : channel.empty()

        sample_db_ch = (
            params.penncnv_qc_path && params.plink2samplemetadata
            ) ? mergeSampleMetadata(
                    params.plink2samplemetadata,
                    penncnv_qc
                )
            : params.penncnv_qc_path
                ? stageFile4(params.penncnv_qc_path, "sampleDB.tsv")
            : params.plink2samplemetadata
                ? stageFile5(params.plink2samplemetadata, "sampleDB.tsv")
            : channel.empty()

    }


    '''
    PREFILTERING - MERGING
    '''
    MERGE_CNV_CALLS (   quantisnp_cnv,                                                             
                        penncnv_cnv,
                        file("${projectDir}/resources/Genome_Regions/Genome_Regions_data.tsv"),
                        params.genome_version
                    )

    merged_cnv = MERGE_CNV_CALLS.out.merged_cnv_ch
     

    // Run report only if requested
    if (params.report) {

        MAKE_REPORT (    params.dataset_name,
                         params.plink2samplemetadata,
                         plink_ch.triofile,
                         penncnv_qc,
                         penncnv_cnv_raw,
                         quantisnp_cnv_raw,
                         merged_cnv                     
                    )
    }

    '''
    SUMMARY BUILDER
    '''
    
    buildSummary  ( params.sample_file,
                    params.git_hash,
                    params.dataset_name,
                    params.genome_version,
                    params.report ? MAKE_REPORT.out.merged_cnv_qc : merged_cnv 
                  )


    
    publish:

    // MERGED CNV DATASET
    merged_cnv = MERGE_CNV_CALLS.out.merged_cnv_ch
    sampleDB   = sample_db_ch ?: channel.empty()    //publish even if empty

    // Before filter results
    penncnv_qc        = penncnv_qc ?: channel.empty()
    penncnv_cnv_raw   = penncnv_cnv_raw
    penncnv_cnv       = penncnv_cnv
    quantisnp_cnv_raw = quantisnp_cnv_raw
    quantisnp_cnv     = quantisnp_cnv

    // REPORT outputs only if report was run
    penncnv_unfilter_cnv_qc   = params.report ? MAKE_REPORT.out.penncnv_unfilter_cnv_qc   : channel.empty()
    quantisnp_unfilter_cnv_qc = params.report ? MAKE_REPORT.out.quantisnp_unfilter_cnv_qc : channel.empty()
    merged_cnv_qc             = params.report ? MAKE_REPORT.out.merged_cnv_qc             : channel.empty()
    sample_qc_report          = params.report ? MAKE_REPORT.out.sample_qc_report          : channel.empty()

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
