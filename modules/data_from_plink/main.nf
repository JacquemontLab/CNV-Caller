#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow {
    // Define the base output path for PLINK dataset
    plink_base_path = "/home/flben/projects/rrg-jacquese/All_user_common_folder/RAW_DATA/Genetic/ALSPAC/PLINK/ALSPAC"
    
    // Extract the base name (prefix) without directory and extension
    plink_prefix = plink_base_path.split('/').last()  // Extract "prefix"

    // Define the input channel for PLINK dataset prefix and probe file
    plink_set = Channel.of(
        tuple(
            plink_prefix,
                [
                    file("${plink_base_path}.bed"),
                    file("${plink_base_path}.bim"),
                    file("${plink_base_path}.fam")
                ]
            )
        )

    Channel
        .fromPath("/home/flben/projects/rrg-jacquese/All_user_common_folder/RAW_DATA/Genetic/ALSPAC/BAF_LRR_Probes_by_sample/ALSPAC09900260.BAF_LRR_Probes.tsv")       // A TSV file with probe IDs in column 1
        .set { probes_file }

    // Run the metadata collection process
    collect_plink_metadata(plink_set, probes_file)
}


// Nextflow process to compute call rate and infer sex from PLINK files
process collect_plink_metadata {
    tag "collect_plink_metadata"  // Optional identifier for tracking/logging

    input:
    tuple val(plink_prefix), path(plink_files)                 // Path to the PLINK binary dataset (prefix of .bed/.bim/.fam files)
    path one_baf_lrr_probes_file // A TSV file with probe IDs in the first column (used for filtering SNPs)

    output:
    path "plink_metadata.tsv"  // Final output file with sample ID, call rate, and imputed sex

    script:
    """
    echo "Process Running: collect_plink_metadata"

    plink_get_sex_callrate.sh "$plink_prefix" "$one_baf_lrr_probes_file" "plink_metadata.tsv"
    """
}

