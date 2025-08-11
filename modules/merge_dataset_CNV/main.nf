#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
 * Process: mergeCNVcalls
 * Merges CNV calls from PennCNV and QuantiSNP using a shell script.
 * Strips quotes from input files and runs `merge_cnv_quantisnp_penncnv.sh`.
 *
 * Inputs:
 *   - quantisnp_cnv_raw: Raw QuantiSNP CNV file
 *   - penncnv_cnv_raw:   Raw PennCNV CNV file
 *   - regions_file:      Genome regions annotation file
 *   - genome_version:    Genome build (e.g. GRCh37)
 *
 * Output:
 *   - merged_cnv.tsv: Combined CNV call file
 */
process mergeCNVcalls {
    label 'merge_cnv'

    input:
    path quantisnp_cnv_raw
    path penncnv_cnv_raw
    path regions_file
    val  genome_version

    output:
    path "CNV_merged_dataset.tsv",   emit: merged_cnv

    script:
    """
    # Remove quotes from input files
    sed 's/"//g' "$quantisnp_cnv_raw" > qs_no_quotes.txt
    sed 's/"//g' "$penncnv_cnv_raw" > pc_no_quotes.txt

    merge_cnv_quantisnp_penncnv.sh qs_no_quotes.txt pc_no_quotes.txt ${regions_file} ${genome_version} "CNV_merged_dataset.tsv"
    """
}

/*
 * Workflow: MERGE_CNV_CALLS
 * Wraps the CNV merge process.
 *
 * Emits:
 *   - merged_cnv: Merged CNV calls
 */
workflow MERGE_CNV_CALLS {
    take:
    quantisnp_cnv_raw
    penncnv_cnv_raw
    regions_file
    genome_version

    main:
    mergeCNVcalls ( quantisnp_cnv_raw, 
                    penncnv_cnv_raw,
                    regions_file,
                    genome_version)

    merged_cnv_ch = mergeCNVcalls.out.merged_cnv

    emit:
        merged_cnv_ch
}