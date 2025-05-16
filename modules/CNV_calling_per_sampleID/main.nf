#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * This Nextflow pipeline performs Copy Number Variation (CNV) calling on selected samples using PennCNV and QuantiSNP, 
 * merges the results, and extracts quality control (QC) metrics. It filters input samples from a provided list, runs 
 * CNV detection on each sample, and organizes the outputs for downstream analysis.
 */

// NXF_OFFLINE=true nextflow run main.nf -resume

// Define the process to run CNV calling
process callBatchCNVs {

    input:
    //val sample_id
    val BAF_LRR_Probes
    path pfb_file
    path gcmodel_file
    path from_plink_extracted_data
    path gcdir

    output:
    path "*.penncnv.cnv",   emit: penn_cnv
    path "*.penncnv.log",   emit: penn_log
    path "*.quantisnp.cnv", emit: quanti_cnv
    path "batch_list.txt",  emit: batch_list
    
    script:
    """
    #turn the nextflow variable into lines of a file
    echo ${BAF_LRR_Probes} | sed 's/[][]//g' | tr ',' '\\n' | tr -d " " > batch_list.txt

    batch_cnv_call.sh --batch_list batch_list.txt \
                      --plink_data ${from_plink_extracted_data} \
                      --gcmodel ${gcmodel_file} \
                      --gcdir ${gcdir} \
                      --pfb ${pfb_file}
    """
}


// Merge CNV callers results and extract QC metrics
process mergeCNVCallers {
    tag "merge_cnv_callers_and_extract_qc: ${sample_id}"

    input:
    tuple val(sample_id), file(BAF_LRR_Probes), file(quantisnp_file), file(penncnv_file), file(penncnv_logfile)
    path regions_file
    val genome_version

    output:
    path "${sample_id}.PennCNV_QC.tsv", emit: PennCNV_QC_tsv
    path "${sample_id}.CNVs.tsv", emit: CNVs_tsv

    script:
    """
    echo "Process Running: merge_cnv_callers_and_extract_qc for ${sample_id}"

    merge_cnv_quantisnp_penncnv.sh "$quantisnp_file" "$penncnv_file" "$BAF_LRR_Probes" "$regions_file" "$genome_version" ${sample_id}.CNVs.tsv

    extract_qc_penncnv.sh "$penncnv_logfile" ${sample_id}.PennCNV_QC.tsv
    """
}
//Index helper function. This function builds a tuple from a filepath where the first entry is file-prefix and the second is the full filepath.
def index(file_ch){
    return file_ch.flatten() | map {f -> tuple(f.getSimpleName(), f)}
}

workflow  CALL_CNV_PARALLEL {
    take:
    list_sample_file    //file containing a list of filepaths to probe files. Value Channel
    pfb                 //pfb file generated from prepare_penncnv_params
    gc_content_windows  //gc model from resources
    plink               //plink data extracted using extract_plink_data
    gcDir               //resource directory pointing to per-chromosome 1k binned gc content regions

    main:

    //Splitting into groups by splitting csv sample file
    batch_ch = list_sample_file.splitCsv(by : 64)    //batch size of 64
                    .take( 3 )            //debug, take first 3 batches

    //Calling CNVs    
    callBatchCNVs ( batch_ch, 
                    pfb,
                    gc_content_windows,
                    plink,
                    gcDir               )

    
    //turning batch sample list into tuple pairs
    sample_paths_ch = callBatchCNVs.out.batch_list
                        .splitCsv()
                        .flatten()
                        .map { f -> tuple(file(f).getSimpleName(),file(f))}

    
    //building sample-level tuples from batched output via serial joins. Each output 
    //from index() is a tuple(sampleID, filepath) which can be matched with other channels 
    //to create one large tuple.
     sample_cnvs_ch = sample_paths_ch.join( index(callBatchCNVs.out.quanti_cnv ))
                                     .join( index(callBatchCNVs.out.penn_cnv   ))
                                     .join( index(callBatchCNVs.out.penn_log   ))
    emit:
    sample_cnvs_ch //sample-level cnvs from both cnv callers plus QC files
    
}