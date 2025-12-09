#!/usr/bin/env nextflow

nextflow.enable.dsl=2


//import subworkflows
include { PREPARE_PENNCNV_INPUTS } from './modules/prepare_penncnv_inputs'
include { CALL_CNV               } from './modules/call_cnv'
include { MERGE_CNV_CALLS        } from './modules/merge_cnv_calls'
include { MAKE_REPORT            } from './modules/make_report'



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

process extractPlink {
    input:
    path plink2samplemetadata_tsv

    output:
    path 'sexfile.tsv', emit: sexfile
    path 'trio.tsv', emit: triofile

    script:
    """
    cut -f1,3 ${plink2samplemetadata_tsv} > 'sexfile.tsv'
    cut -f1,4,5 ${plink2samplemetadata_tsv} > 'trio.tsv'
    """

}

process mergeSampleMetadata {
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
        SELECT 
            CAST(b.SampleID AS VARCHAR) AS SampleID, 
            b.* EXCLUDE SampleID,
            a.* EXCLUDE SampleID
        FROM read_csv_auto('$penncnv_qc', sep='\t', header=true) b
        LEFT JOIN read_csv_auto('$plink2samplemetadata_tsv', sep='\t', header=true) a
        ON CAST(b.SampleID AS VARCHAR) = CAST(a.SampleID AS VARCHAR)
    )  TO 'sampleDB.tsv' (HEADER, DELIMITER '\t');
    "
    """
}

process formatRawCNV {
    tag "formating Raw CNV output"

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




//Default params
params.pipeline_mode            = "full"
params.list_sample_baflrrpath   = ""
params.penncnv_qc_path          = ""
params.penncnv_calls_path       = ""
params.quantisnp_calls_path     = ""
params.plink2samplemetadata_tsv = ""
params.genome_version           = "GRCh38"
params.batch_size               = 64
params.pfb_max_sample_size      = 1000
params.report                   = false
params.data_type                = "array"
params.dataset_name             = ""



workflow FORMAT_CNV {  
    take :
        penncnv_calls_path
        quantisnp_calls_path
        penncnv_qc_path
        plink2samplemetadata_tsv

    main :

        formatRawCNV(penncnv_calls_path, quantisnp_calls_path)

        if (file(penncnv_qc_path).exists() && plink2samplemetadata_tsv){

            mergeSampleMetadata( plink2samplemetadata_tsv, penncnv_qc_path )   
        }

    emit:
    formatted_penncnv    = formatRawCNV.out.penncnv_cnv
    formatted_quantisnp  = formatRawCNV.out.quantisnp_cnv
    sample_DB            = mergeSampleMetadata.out ?: Channel.empty()
    penncnv_qc           = penncnv_qc_path

}

workflow {
    
    main:

    if (params.pipeline_mode == "full"){
        
        list_sample_baflrrpath   = Channel.fromPath(params.list_sample_baflrrpath)

        batch_ch = list_sample_baflrrpath.splitCsv(sep: "\t",  header: true)                       
                                         .map {row -> row['path_to_BAF_LRR']}                              // grab filepaths only
                                         .buffer( size : params.batch_size, remainder : true)              // split channel into batches  
               
        


        plink_ch = extractPlink(params.plink2samplemetadata_tsv)
       
        '''
        PREPARE INPUTS for PennCNV
        '''

        PREPARE_PENNCNV_INPUTS ( list_sample_baflrrpath, //CHANGE TO FILE CHANNEL OF ALL SAMPLES
                                 params.plink2samplemetadata_tsv,
                                 file("${projectDir}/resources/GC_correction/${params.genome_version}/gc_content_1k_windows.bed"),
                                 params.pfb_max_sample_size,
                                 params.data_type                                                                                        
                                )

        '''
        CALLING CNVs AND MERGE
        '''
        CALL_CNV              ( batch_ch,                                                               // File of paths to baf_lrr files without the sampleID
                                PREPARE_PENNCNV_INPUTS.out.pfb_file.first(),                            // PFB file, passing into value channel using first()
                                PREPARE_PENNCNV_INPUTS.out.hmm_file.first(),                            // HMM file, passing into value channel using first()
                                PREPARE_PENNCNV_INPUTS.out.gc_model.first(),                            // GC model
                                plink_ch.sexfile.first(),                                               // Sexfile from metadata input 
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
        
        //Make SampleDB
        sample_db_ch = mergeSampleMetadata(params.plink2samplemetadata_tsv, penncnv_qc)

    
    } else if (params.pipeline_mode == "partial") {
        
       FORMAT_CNV (    
                       file(params.penncnv_calls_path),
                       file(params.quantisnp_calls_path),
                       file(params.penncnv_qc_path),
                       file(params.plink2samplemetadata_tsv)       
                  )

        //formatting output from partial run to mask variables
        sample_db_ch      = FORMAT_CNV.out.sample_DB
        penncnv_qc        = FORMAT_CNV.out.penncnv_qc
        penncnv_cnv_raw   = file(params.penncnv_calls_path)
        quantisnp_cnv_raw = file(params.quantisnp_calls_path)
        quantisnp_cnv     = FORMAT_CNV.out.formatted_quantisnp
        penncnv_cnv       = FORMAT_CNV.out.formatted_penncnv
    }

    //RUN for both partial and full runs 
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
        trio_file = plink_ch.triofile

        MAKE_REPORT (    params.dataset_name,
                         params.plink2samplemetadata_tsv,
                         trio_file,
                         penncnv_qc,
                         penncnv_cnv_raw,
                         quantisnp_cnv_raw,
                         merged_cnv                     \
                    )
    }

    '''
    SUMMARY BUILDER
    '''
    
    buildSummary  ( params.dataset_name,
                    params.genome_version,
                    params.report ? MAKE_REPORT.out.merged_cnv_qc : merged_cnv )


    
    publish:

    // MERGED CNV DATASET
    merged_cnv = MERGE_CNV_CALLS.out.merged_cnv_ch
    sampleDB   = sample_db_ch ?: Channel.empty()                                                            //publish even if empty

    // Before filter results
    penncnv_qc        = penncnv_qc ?: (CALL_CNV.out.penncnv_qc_ch ?: Channel.empty())              //only emits from the full run 
    penncnv_cnv_raw   = penncnv_cnv_raw
    penncnv_cnv       = penncnv_cnv
    quantisnp_cnv_raw = quantisnp_cnv_raw
    quantisnp_cnv     = quantisnp_cnv

    // REPORT outputs only if report was run
    penncnv_unfilter_cnv_qc   = params.report ? MAKE_REPORT.out.penncnv_unfilter_cnv_qc   : Channel.empty()
    quantisnp_unfilter_cnv_qc = params.report ? MAKE_REPORT.out.quantisnp_unfilter_cnv_qc : Channel.empty()
    merged_cnv_qc             = params.report ? MAKE_REPORT.out.merged_cnv_qc             : Channel.empty()
    sample_qc_report          = params.report ? MAKE_REPORT.out.sample_qc_report          : Channel.empty()

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