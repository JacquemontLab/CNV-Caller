#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * STAGE_FILE
 * -----------
 * This process stages (renames or moves) an input file to a specified output name to be published.
 *
 * Purpose:
 * - Ensure that downstream processes receive files with standardized or expected names.
 * - Avoid unnecessary file operations when the input and output filenames are identical.
 */

process STAGE_FILE {
    tag "$outname"

    input:
    path infile
    val  outname

    output:
    path outname

    script:
    """
    if [ "${infile}" != "${outname}" ]; then
        mv "${infile}" "${outname}"
    else
        echo "Input and output are the same, skipping copy."
    fi
    """
}