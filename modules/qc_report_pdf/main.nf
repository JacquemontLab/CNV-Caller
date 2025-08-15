#!/usr/bin/env nextflow

nextflow.enable.dsl=2



process formating_quantisnp{
    tag "formating_quantisnp"

    input:
    path quantisnp_cnv_raw

    output:
    path "QuantiSNP_CNV.tsv"

    script:
    """
    format_quantisnp_cnv.sh "${quantisnp_cnv_raw}" "QuantiSNP_CNV.tsv"
    """
}



process formating_penncnv{
    tag "formating_penncnv"

    input:
    path penncnv_cnv_raw

    output:
    path "PennCNV_CNV.tsv"

    script:
    """
    format_penncnv_cnv.sh "${penncnv_cnv_raw}" "PennCNV_CNV.tsv"
    """
}


process sample_qc_report {
    tag "sample_qc_report"

    input:
    val dataset_name                      // dataset_name, used in the report
    path plink2samplemetadata_tsv   // Metadata: SampleID -> info (LRR_SD, BAF_SD, WF )
    path penncnv_qc                 // QC table: SampleID -> QC metrics
    path sample_qc_rmd  // <- This is your .Rmd file

    output:
    path "sample_qc_report.pdf" // Final rendered report

    script:
    """
    # Create merged header (first column shared, rest from both headers except duplicate ID)
    (echo -e "\$(head -n 1 ${plink2samplemetadata_tsv} | cut -f1)\t\$(head -n 1 ${plink2samplemetadata_tsv} | cut -f2-)\t\$(head -n 1 ${penncnv_qc} | cut -f2-)"; \
     join -t \$'\\t' -1 1 -2 1 \
       <(tail -n +2 ${plink2samplemetadata_tsv} | sort) \
       <(tail -n +2 ${penncnv_qc} | sort)) > qc_input.tsv


    Rscript -e "rmarkdown::render('$sample_qc_rmd', params=list(
        dataset_name='$dataset_name',
        qc_input='qc_input.tsv'),
    output_file='sample_qc_report.pdf'
    )"
    """
}


process denovo_cnv_annotation_type {
    tag "denovo_cnv_annotation and type"

    input:
    path input_cnv
    path trio_file

    output:
    path "denovo_cnv_annotation.tsv"

    script:
    """

    ${workflow.projectDir}/modules/merge_dataset_CNV/resources/bin/infer_cnv_type.py -i "$input_cnv" -o "cnv_type.tsv" || true

    cnv_trio_inheritance.py --cnv "cnv_type.tsv" --pedigree $trio_file --type_col Type --overlap 0.5 --output denovo_cnv_annotation.tsv

    """
}

process unfilter_cnv_qc {
    tag "unfilter_cnv_qc"

    input:
    val dataset_name                      // dataset_name, used in the report
    path input_with_denovo_annotation
    path unfilter_cnv_qc_rmd  // <- This is your .Rmd file
    val output_name

    output:
    path "${output_name}" // Captures the actual output file name

    script:
    """

    Rscript -e "rmarkdown::render('$unfilter_cnv_qc_rmd', params=list(dataset_name='$dataset_name',
     dataset='$input_with_denovo_annotation', x_var_list='Confidence,Length,Num_Probes'),
        output_file='$output_name'
        )"
    """
}


process denovo_cnv_annotation{
    tag "denovo_cnv_annotation"

    input:
    path input_cnv
    path trio_file

    output:
    path "denovo_cnv_annotation.tsv"

    script:
    """
    cnv_trio_inheritance.py --cnv "$input_cnv" --pedigree $trio_file --type_col Type --overlap 0.5 --output denovo_cnv_annotation.tsv
    """
}

process merged_cnv_qc {
    tag "merged_cnv_qc"

    input:
    val dataset_name                      // dataset_name, used in the report
    path input_with_denovo_annotation
    path merged_cnv_qc_rmd  // <- This is your .Rmd file

    output:
    path "merged_cnv_qc.pdf" // Final rendered report

    script:
    """

    Rscript -e "rmarkdown::render('$merged_cnv_qc_rmd', params=list(dataset_name='$dataset_name',
            dataset='$input_with_denovo_annotation', x_var_list='Confidence_max,Two_Algorithm_Overlap'),
            output_file='merged_cnv_qc.pdf'
            )"
    """
}




workflow do_merged_cnv_qc {
    take:
    dataset_name
    merged_cnv
    trio_file

    main:
    merged_cnv_qc_rmd = file("${projectDir}/modules/qc_report_pdf/resources/bin/cnv_dataset_qc.Rmd")
    merged_cnv_denovo = denovo_cnv_annotation(merged_cnv, trio_file)
    merged_cnv_qc = merged_cnv_qc(dataset_name, merged_cnv_denovo, merged_cnv_qc_rmd)

    emit: merged_cnv_qc
}



workflow do_quantisnp_qc {
    take:
    dataset_name
    quantisnp_cnv
    trio_file

    main:
    unfilter_cnv_qc_rmd = file("${projectDir}/modules/qc_report_pdf/resources/bin/cnv_unfilter_qc_report.Rmd")
    quantisnp_denovo = denovo_cnv_annotation_type(quantisnp_cnv, trio_file)
    quantisnp_unfilter_cnv_qc = unfilter_cnv_qc(dataset_name, quantisnp_denovo, unfilter_cnv_qc_rmd,"quantisnp_unfilter_cnv_qc.pdf")

    emit: quantisnp_unfilter_cnv_qc
}

workflow do_penncnv_qc {
    take:
    dataset_name
    penncnv_cnv
    trio_file

    main:
    unfilter_cnv_qc_rmd = file("${projectDir}/modules/qc_report_pdf/resources/bin/cnv_unfilter_qc_report.Rmd")
    penncnv_denovo = denovo_cnv_annotation_type(penncnv_cnv, trio_file)
    penncnv_unfilter_cnv_qc = unfilter_cnv_qc(dataset_name, penncnv_denovo, unfilter_cnv_qc_rmd,"penncnv_unfilter_cnv_qc.pdf")

    emit: penncnv_unfilter_cnv_qc
}


workflow REPORT_PDF {
    take:
    dataset_name
    plink2samplemetadata_tsv
    trio_file
    penncnv_qc
    penncnv_cnv_raw
    quantisnp_cnv_raw
    merged_cnv


    main:
    sample_qc_rmd = file("${projectDir}/modules/qc_report_pdf/resources/bin/sample_qc_report.Rmd")

    sample_qc_report = sample_qc_report(dataset_name,
                    plink2samplemetadata_tsv,
                    penncnv_qc,
                    sample_qc_rmd
                    )

    penncnv_cnv = formating_penncnv(penncnv_cnv_raw)
    quantisnp_cnv = formating_quantisnp(quantisnp_cnv_raw)
    
    // Call the sub-workflows
    penncnv_unfilter_cnv_qc = do_penncnv_qc("${dataset_name} PennCNV Unfilter",penncnv_cnv, trio_file)
    quantisnp_unfilter_cnv_qc = do_quantisnp_qc("${dataset_name} QuantiSNP Unfilter",quantisnp_cnv, trio_file)
    merged_cnv_qc = do_merged_cnv_qc("${dataset_name} CNV merged dataset",merged_cnv, trio_file)

    emit:
    penncnv_unfilter_cnv_qc
    quantisnp_unfilter_cnv_qc
    merged_cnv_qc
    sample_qc_report
}
