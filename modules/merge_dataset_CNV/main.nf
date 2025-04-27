#!/usr/bin/env nextflow

nextflow.enable.dsl=2
// NXF_OFFLINE=true nextflow run main.nf -resume -c nextflow.config -with-trace -with-report

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

workflow {
    // Group all files into a single list (tuple)
    file_channel = Channel.from([
        file("/home/flben/projects/rrg-jacqueese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/merge_dataset_CNV/test/ALSPAC09897249.CNVs.tsv"),
        file("/home/flben/projects/rrg-jacqueese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/merge_dataset_CNV/test/ALSPAC09902309.CNVs.tsv")
    ])

    // Then, emit as a single value
    file_channel = file_channel.collect()

    // Launch the process
    concatenated_CNVs(file_channel)


    // Group all files into a single list (tuple)
    file_channel_qc = Channel.from([
        file("/home/flben/projects/rrg-jacqueese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/merge_dataset_CNV/test/ALSPAC09897249.PennCNV_QC.tsv"),
        file("/home/flben/projects/rrg-jacqueese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/merge_dataset_CNV/test/ALSPAC09902309.PennCNV_QC.tsv")
    ])

    // Then, emit as a single value
    file_channel_qc = file_channel_qc.collect()

    // Launch the process
    concatenated_penncnv_qc(file_channel_qc)
}
