#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * This Nextflow pipeline concatenates multiple CNV and PennCNV QC files
 * into two unified output tables using a custom merge script.
 */
 
// NXF_OFFLINE=true nextflow run main.nf -resume -c nextflow.config -with-trace -with-report

// Define the first process: concatenate CNV files
process concatenated_CNVs {
    tag "concatenated_CNVs"
    input:
    path files   // Input: A list of files to concatenate

    output:
    path "concatenated_CNVs.tsv"  // Output: The concatenated result

    script:
    """
    echo "Process Running: concatenated_CNVs"
    
    merge_files.sh --directory "." --chunk_size 2 --num_cores 1 --output_file concatenated_CNVs.tsv --file_extension ".CNVs.tsv"
    """
}

process genomic_region_overlap {
    tag "genomic_region_overlap"
    input:
    tuple path(concatenated_CNVs), path(regions_file), val(genome_version)

    output:
    path "CNVs_with_genomic_regions_overlap.tsv"  // Output: The concatenated result

    script:
    """
    echo "Process Running: genomic_region_overlap"
    
    add_regions_overlap.sh "$concatenated_CNVs" "$regions_file" "$genome_version" "CNVs_with_genomic_regions_overlap.tsv"
    """
}


// Define the second process: concatenate PennCNV QC files
process concatenated_penncnv_qc {
    tag "concatenated_penncnv_qc"
    input:
    path files   // Input: A list of files to concatenate

    output:
    path "concatenated_penncnv_qc.tsv"  // Output: The concatenated result

    script:
    """
    echo "Process Running: concatenated_penncnv_qc"

    merge_files.sh --directory "." --chunk_size 2 --num_cores 1 --output_file concatenated_penncnv_qc.tsv --file_extension ".PennCNV_QC.tsv"
    """
}

// Define the workflow section: controls how data moves through processes
workflow TEST{
    // Create a Channel from a list of CNV files
    file_channel = Channel.from([
        file("/lustre06/project/6008022/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/merge_dataset_CNV/test/ALSPAC09897249.CNVs.tsv"),
        file("/lustre06/project/6008022/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/merge_dataset_CNV/test/ALSPAC09902309.CNVs.tsv")
    ])

    file_channel.collectFile(name: "concatenated_CNVs.tsv", keepHeader: true, save: true)

    // Then, emit as a single value
    file_channel = file_channel.collect()

    // Launch the process
    concatenated_CNVs(file_channel)


    // Create a Channel from a list of PennCNV QC files)
    file_channel_qc = Channel.from([
        file("/lustre06/project/6008022/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/merge_dataset_CNV/test/ALSPAC09897249.PennCNV_QC.tsv"),
        file("/lustre06/project/6008022/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/merge_dataset_CNV/test/ALSPAC09902309.PennCNV_QC.tsv")
    ])

    // Run the first process with the CNV files
    file_channel_qc = file_channel_qc.collect()

    // Run the second process with the QC files
    concatenated_penncnv_qc(file_channel_qc)
}
