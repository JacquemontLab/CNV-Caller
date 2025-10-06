#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * This Nextflow pipeline performs Copy Number Variation (CNV) calling on selected samples using PennCNV and QuantiSNP, 
 * merges the results, and extracts quality control (QC) metrics. It filters input samples from a provided list, runs 
 * CNV detection on each sample, and organizes the outputs for downstream analysis.
 */


// Define the process to run CNV calling
process callBatchCNVs {
    label "cnv_calling"

    input:
    val BAF_LRR_Probes
    path pfb_file
    path gcmodel_file
    path sexfile
    val genome_version
    val autosome_only

    output:
    path "*.penncnv.qc",        emit: penncnv_qc_raw
    path "*.PennCNV_QC.tsv",    emit: penncnv_qc
    path "*.penncnv.cnv",       emit: penncnv_cnv_raw
    path "*.penncnv.cnv.tsv",   emit: penncnv_cnv
    path "*.quantisnp.cnv",     emit: quantisnp_cnv_raw
    path "*.quantisnp.cnv.tsv", emit: quantisnp_cnv
    path "batch_list.txt",      emit: batch_list

    script:
    if(autosome_only){
    
    """
    # Turn the nextflow variable into lines of a file
    echo ${BAF_LRR_Probes} | sed 's/[][]//g' | tr ',' '\\n' | tr -d " " > batch_list.txt

    # Default parameters avalable in the docker:
    chr="1:23"
    gcdir=/usr/local/QuantiSNP-2.3/GC_correction/${genome_version}/GCdir/
    hmm_file="/usr/local/PennCNV-1.0.5/lib/wgs.hmm"
    levels="/usr/local/QuantiSNP-2.3/bin/config/levels.dat"
    config="/usr/local/QuantiSNP-2.3/bin/config/params.dat"

    batch_cnv_call.sh --batch_list batch_list.txt \
                    --autosome_only \
                    --sexfile ${sexfile} \
                    --pfb ${pfb_file} \
                    --gcmodel ${gcmodel_file} \
                    --gcdir \$gcdir \
                    --hmm_file \$hmm_file \
                    --levels \$levels \
                    --config \$config \
                    --chr \$chr \
                    --mode taskset \
                    --cpus ${task.cpus}
    """
    } else {
    """
    # Turn the nextflow variable into lines of a file
    echo ${BAF_LRR_Probes} | sed 's/[][]//g' | tr ',' '\\n' | tr -d " " > batch_list.txt

    # Default parameters avalable in the docker:
    chr="1:23"
    gcdir=/usr/local/QuantiSNP-2.3/GC_correction/${genome_version}/GCdir/
    hmm_file="/usr/local/PennCNV-1.0.5/lib/wgs.hmm"
    levels="/usr/local/QuantiSNP-2.3/bin/config/levels.dat"
    config="/usr/local/QuantiSNP-2.3/bin/config/params.dat"

    batch_cnv_call.sh --batch_list batch_list.txt \
                    --sexfile ${sexfile} \
                    --pfb ${pfb_file} \
                    --gcmodel ${gcmodel_file} \
                    --gcdir \$gcdir \
                    --hmm_file \$hmm_file \
                    --levels \$levels \
                    --config \$config \
                    --chr \$chr \
                    --mode taskset \
                    --cpus ${task.cpus}
    """

    }
}


workflow  CALL_CNV_PARALLEL {
    take:
    list_baflrr_path    //file containing a list of filepaths to probe files. Value Channel
    pfb                 //pfb file generated from prepare_penncnv_params
    gc_content_windows  //gc model from resources
    sexfile             //plink data extracted using extract_plink_data
    genome_version               //resource directory pointing to per-chromosome 1k binned gc content regions
    autosome_only       //boolean for skipping x-chromosome calling

    main:

    
    callBatchCNVs ( list_baflrr_path, 
                    pfb,
                    gc_content_windows,
                    sexfile,
                    genome_version,
                    autosome_only               )


                    
    emit:
        penncnv_qc_ch        = callBatchCNVs.out.penncnv_qc
        penncnv_cnv_ch       = callBatchCNVs.out.penncnv_cnv
        penncnv_cnv_raw_ch   = callBatchCNVs.out.penncnv_cnv_raw
        quantisnp_cnv_ch     = callBatchCNVs.out.quantisnp_cnv
        quantisnp_cnv_raw_ch = callBatchCNVs.out.quantisnp_cnv_raw
    
}