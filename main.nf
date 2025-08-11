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
    // Define inputs
    list_path_to_BAF_LRR = params.list_path_to_BAF_LRR
    plink2samplemetadata_tsv = params.plink2samplemetadata_tsv
    gc_correction_dir = params.gc_correction_dir
    genome_version = params.genome_version
    batch_size = params.batch_size
    dataset_name = params.dataset_name
    
    '''
    PREPARE INPUTS for PennCNV
    '''
    // Extract the first file from the file 
    // first_sample = Channel.fromPath(file(list_path_to_BAF_LRR))
    //                             .splitCsv()
    //                             .first() // select first row
    //                             .map { f -> file(f[0].toString())} // convert its content to a file

    baf_lrr_files = Channel.fromPath(file(list_path_to_BAF_LRR))
                .splitText()
                .map { line -> 
                    def cols = line.split('\t')
                    cols[1]  // index 1 = second column
                }


    // PREPARE_PENNCNV_INPUTS ( first_sample.getParent(),
    PREPARE_PENNCNV_INPUTS ( list_path_to_BAF_LRR,
                             plink2samplemetadata_tsv,
                             gc_correction_dir,
                             genome_version )


    sex_file = sexfile_from_plink2samplemetadata(plink2samplemetadata_tsv)
    trio_file = trio_from_plink2samplemetadata(plink2samplemetadata_tsv)

    '''
    CALLING CNVs AND MERGE
    '''
    // CALL_CNV_PARALLEL    ( Channel.fromPath(file(list_path_to_BAF_LRR)),
    CALL_CNV_PARALLEL    ( baf_lrr_files,
                           PREPARE_PENNCNV_INPUTS.out.pfb_file, 
                           PREPARE_PENNCNV_INPUTS.out.gc_model,
                           sex_file,
                           genome_version,
                           batch_size )


    '''
    PREFILTERING - MERGING
    '''
    MERGE_CNV_CALLS ( CALL_CNV_PARALLEL.out.quantisnp_cnv_raw_ch,
                    CALL_CNV_PARALLEL.out.penncnv_cnv_raw_ch,
                    file("${projectDir}/resources/Genome_Regions/Genome_Regions_data.tsv"),
                    genome_version
                     )

    penncnv_qc = CALL_CNV_PARALLEL.out.penncnv_qc_ch
    penncnv_cnv = CALL_CNV_PARALLEL.out.penncnv_cnv_ch
    quantisnp_cnv = CALL_CNV_PARALLEL.out.quantisnp_cnv_ch

    merged_cnv = MERGE_CNV_CALLS.out.merged_cnv_ch

    REPORT_PDF ( 
        dataset_name,
        plink2samplemetadata_tsv,
        trio_file,
        penncnv_qc,
        penncnv_cnv,
        quantisnp_cnv,
        merged_cnv)

    '''
    SUMMARY BUILDER
    '''
                                           
    buildSummary  (dataset_name,
                    genome_version,
                    REPORT_PDF.out.quantisnp_raw_cnv_qc )


    publish:

// MERGED CNV DATASET
    merged_cnv = MERGE_CNV_CALLS.out.merged_cnv_ch

// Before filter results
    penncnv_qc = CALL_CNV_PARALLEL.out.penncnv_qc_ch
    penncnv_cnv = CALL_CNV_PARALLEL.out.penncnv_cnv_ch
    penncnv_cnv_raw = CALL_CNV_PARALLEL.out.penncnv_cnv_raw_ch
    quantisnp_cnv = CALL_CNV_PARALLEL.out.quantisnp_cnv_ch
    quantisnp_cnv_raw = CALL_CNV_PARALLEL.out.quantisnp_cnv_raw_ch 

// REPORT
    penncnv_raw_cnv_qc = REPORT_PDF.out.penncnv_raw_cnv_qc
    quantisnp_raw_cnv_qc = REPORT_PDF.out.quantisnp_raw_cnv_qc
    merged_cnv_qc = REPORT_PDF.out.merged_cnv_qc
    sample_qc_report = REPORT_PDF.out.sample_qc_report

}


output {
    merged_cnv {
        mode 'copy'
    }
    penncnv_qc {
        mode 'copy'
        path "calls_unfiltered/penncnv/"
    }
    penncnv_cnv {
        mode 'copy'
        path "calls_unfiltered/penncnv/"
    }
    penncnv_cnv_raw {
        mode 'copy'
        path "calls_unfiltered/penncnv/"
    }
    quantisnp_cnv {
        mode 'copy'
        path "calls_unfiltered/quantisnp/"
    }
    quantisnp_cnv_raw {
        mode 'copy'
        path "calls_unfiltered/quantisnp/"
    }

    quantisnp_raw_cnv_qc {
        mode 'copy'
        path "docs/"
    }

    penncnv_raw_cnv_qc {
        mode 'copy'
        path "docs/"
    }

    merged_cnv_qc {
        mode 'copy'
        path "docs/"
    }

    sample_qc_report {
        mode 'copy'
        path "docs/"
    }

}



