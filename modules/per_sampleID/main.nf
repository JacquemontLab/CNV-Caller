#!/usr/bin/env nextflow

nextflow.enable.dsl=2



process CNV_calling {
    tag "CNV_calling: ${sample_id}"

    input:
    tuple val(sample_id), path(BAF_LRR_Probes), path(pfb_file), path(gcmodel_file), path(plink_metadata)


    output:
    tuple val(sample_id), path(BAF_LRR_Probes), path("${sample_id}.quantisnp.cnv"), path("${sample_id}.penncnv.cnv"), path("${sample_id}.penncnv.log")


    script:
    """
    # Default parameters can be set here
    gcdir="/usr/local/quantisnp/gc/b38/"
    minlength=1000
    minconf=15.0
    chr="1:23"
    hmm_file="/usr/local/bin/hhall.hmm"
    levels="/usr/local/quantisnp/bin/config/levels.dat"
    config="/usr/local/quantisnp/bin/config/params.dat"

    echo "Process Running: CNV_calling for ${sample_id}"

    set -euo pipefail

    echo "Extracting gender from sex file..."
    zgrep -P "${sample_id}\t" ${plink_metadata} | cut -f1,2 > sexfile.tsv
    gender=\$(cut -f2 sexfile.tsv)


    echo "Gender: \$gender"
    if [[ "\$gender" != "male" && "\$gender" != "female" ]]; then
        echo "ERROR: Invalid gender '\$gender' for sample ID ${sample_id}. Expected 'male' or 'female'."
        exit 1
    fi


    echo "Running PennCNV on autosomes..."
    /usr/bin/time -v perl /usr/local/bin/detect_cnv.pl --test \
        --pfbfile ${pfb_file} \
        --hmmfile \${hmm_file} \
        --logfile ${sample_id}.penncnv.log \
        --minlength \${minlength} \
        --minconf \${minconf} \
        --output ${sample_id}.penncnv.out \
        --gcmodelfile ${gcmodel_file} \
        ${BAF_LRR_Probes}
            
    echo "Running PennCNV on ChrX..."
    /usr/bin/time -v perl /usr/local/bin/detect_cnv.pl --test \
        --pfbfile ${pfb_file} \
        --hmmfile \${hmm_file} \
        --logfile ${sample_id}.penncnv.chrx.log \
        --minlength \${minlength} \
        --minconf \${minconf} \
        --output ${sample_id}.penncnv.chrx.out \
        --gcmodelfile ${gcmodel_file} \
        --sexfile sexfile.tsv \
        --chrx \
        ${BAF_LRR_Probes}


    echo "Running QuantiSNP..."
    export MCR_CACHE_ROOT="tmp/mcr_cache"
    mkdir -p "\$MCR_CACHE_ROOT"
    
    /usr/bin/time -v quantisnp --chr \${chr} \
        --outdir . \
        --sampleid ${sample_id} \
        --gender \${gender} \
        --config \${config} \
        --levels \${levels} \
        --gcdir \${gcdir} \
        --input-files ${BAF_LRR_Probes} \
        --doXcorrect --verbose &> ${sample_id}.quantisnp.log

    mv ${sample_id}.cnv ${sample_id}.quantisnp.cnv

    echo "Merging CNV results..."
    cat ${sample_id}.penncnv.chrx.out ${sample_id}.penncnv.out > ${sample_id}.penncnv.cnv
    """
}


process merge_cnv_quantisnp_penncnv {
    tag "merge_cnv_quantisnp_penncnv: ${sample_id}"

    input:
    tuple val(sample_id), file(BAF_LRR_Probes), file(quantisnp_file), file(penncnv_file), file(penncnv_logfile)

    output:
    path "${sample_id}.PennCNV_QC.tsv"
    path "${sample_id}.CNVs.tsv"

    script:
    """
    echo "Process Running: CNV_calling for ${sample_id}"

    merge_cnv_quantisnp_penncnv.sh "$quantisnp_file" "$penncnv_file" "$BAF_LRR_Probes" ${sample_id}.CNVs.tsv

    # You can add quality output generation here if needed
    extract_qc_penncnv.sh "$penncnv_logfile" ${sample_id}.PennCNV_QC.tsv
    """
}


params.list_path_to_BAF_LRR_Probes = "/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/per_sampleID/list.tsv"
params.BAF_LRR_Probes = "/home/flben/projects/rrg-jacquese/All_user_common_folder/RAW_DATA/Genetic/ALSPAC/BAF_LRR_Probes_by_sample/ALSPAC09894851.BAF_LRR_Probes.tsv"
params.pfb_file = "/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/penncnv_params/work/6d/c2e1e23412cba19b2d90889ac413c4/pfb.tsv"
params.gcmodel_file = "/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/penncnv_params/work/fb/1ed66c90f60a69ed6da1c17182c3b8/gcModel.tsv"
params.plink_metadata = "/lustre06/project/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/data_from_plink/work/27/ad6229a5f670de52fd1f7b8fabf8ac/plink_metadata.tsv"


workflow {
    // List of sample IDs you want to keep
    def selected_sample_ids = ["ALSPAC09897249", "ALSPAC09902309"]

    // 1. Input channel from TSV
    Channel
        .fromPath(params.list_path_to_BAF_LRR_Probes)
        .splitCsv(header: true, sep: '\t')
        .filter { row ->                                 
            row.SampleID in selected_sample_ids         // Keep only rows with desired SampleIDs
        }
        .map { row -> 
            [row.SampleID, file(row.Path_to_File_BAF_LRR_Probes)]
        }
        .set { sample_inputs }

    sample_inputs
        .view()

    // 2. Static file channels
    pfb_file_ch = file(params.pfb_file)
    gcmodel_file_ch = file(params.gcmodel_file)
    plink_metadata_ch = file(params.plink_metadata)

    cnv_inputs = sample_inputs
        .map { sample -> 
            tuple(sample[0], sample[1], pfb_file_ch, gcmodel_file_ch, plink_metadata_ch)
        }

    cnv_inputs | CNV_calling | merge_cnv_quantisnp_penncnv
}