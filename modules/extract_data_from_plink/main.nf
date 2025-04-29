#!/usr/bin/env nextflow

nextflow.enable.dsl=2
// NXF_OFFLINE=true nextflow run main.nf -resume -c nextflow.config -with-trace -with-report


// Nextflow process to compute call rate, infer sex and family information from PLINK files
process collect_plink_data {
    tag "collect_plink_data"  // Optional identifier for tracking/logging

    input:
    tuple val(plink_prefix), path(plink_files)                 // Path to the PLINK binary dataset (prefix of .bed/.bim/.fam files)
    path one_baf_lrr_probes_file // A TSV file with probe IDs in the first column (used for filtering SNPs)

    output:
    path "from_plink_extracted_data.tsv", emit: from_plink_extracted_data  // Final output file with sample ID, call rate, and imputed sex

    script:
    """
    echo "Process Running: collect_plink_data"

    from_plink_get_sampleid_data.sh "$plink_prefix" "$one_baf_lrr_probes_file" "from_plink_extracted_data.tsv"
    """
}


workflow PLINK_EXTRACTED_DATA {
    take:
        // Define inputs
        plink_base_path
        probes_file

    main:
        // Extract the base name (prefix) without directory and extension
        plink_prefix = plink_base_path.split('/').last()

        // Define the input channel for PLINK dataset prefix and files
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

        // Call the process (or workflow) that uses these inputs
        from_plink_extracted_data = collect_plink_data(plink_set, probes_file)

    emit:
        from_plink_extracted_data
}


workflow {
    // Define the base output path for PLINK dataset
    plink_base_path = "/lustre06/project/6008022/All_user_common_folder/RAW_DATA/Genetic/ALSPAC/PLINK/ALSPAC"
    
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
        .fromPath("/lustre06/project/6008022/All_user_common_folder/RAW_DATA/Genetic/ALSPAC/BAF_LRR_Probes_by_sample/ALSPAC09900260.BAF_LRR_Probes.tsv")       // A TSV file with probe IDs in column 1
        .set { probes_file }

    // Run the metadata collection process
    collect_plink_data(plink_set, probes_file)
}