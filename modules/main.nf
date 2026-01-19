#!/usr/bin/env nextflow

nextflow.enable.dsl=2


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