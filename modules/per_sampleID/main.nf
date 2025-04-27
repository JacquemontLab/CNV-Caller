#!/usr/bin/env nextflow

nextflow.enable.dsl=2



process CNV_calling {
    tag "CNV_calling: ${sample_id}"

    input:
    tuple val(sample_id), path(BAF_LRR_Probes), path(pfb_file), path(gcmodel_file), path(plink_metadata), path(gcdir)


    output:
    tuple val(sample_id), path(BAF_LRR_Probes), path("${sample_id}.quantisnp.cnv"), path("${sample_id}.penncnv.cnv"), path("${sample_id}.penncnv.log")


    script:
    """
    # Default parameters can be set here
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
        --output ${sample_id}.penncnv.out \
        --gcmodelfile ${gcmodel_file} \
        ${BAF_LRR_Probes}
            
    echo "Running PennCNV on ChrX..."
    /usr/bin/time -v perl /usr/local/bin/detect_cnv.pl --test \
        --pfbfile ${pfb_file} \
        --hmmfile \${hmm_file} \
        --logfile ${sample_id}.penncnv.chrx.log \
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
        --gcdir ${gcdir} \
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
    tuple val(genome_version), file(regions_file)
    tuple val(sample_id), file(BAF_LRR_Probes), file(quantisnp_file), file(penncnv_file), file(penncnv_logfile)

    output:
    path "${sample_id}.PennCNV_QC.tsv"
    path "${sample_id}.CNVs.tsv"

    script:
    """
    echo "Process Running: CNV_calling for ${sample_id}"

    merge_cnv_quantisnp_penncnv.sh "$quantisnp_file" "$penncnv_file" "$BAF_LRR_Probes" "$regions_file" "$genome_version" ${sample_id}.CNVs.tsv

    extract_qc_penncnv.sh "$penncnv_logfile" ${sample_id}.PennCNV_QC.tsv
    """
}



params.list_path_to_BAF_LRR_Probes = "/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/penncnv_params/work/f4/7c4e30*/list_path_to_BAF_LRR_Probes.tsv"
params.pfb_file = "/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/penncnv_params/work/c2/cf2607*/pfb.tsv"
params.gcmodel_file = "/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/penncnv_params/work/a5/6ea189cd1a9453de83807dc408fdf0/gcModel.tsv"
params.gcDir = "/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/penncnv_params/resources/GRCh37_GCdir/"
params.plink_metadata = "/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/data_from_plink/work/76/e35544f2e885e3217ccf31b84e5ef2/plink_metadata.tsv"


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
    gcDir_ch = file(params.gcDir)

    cnv_inputs = sample_inputs
        .map { sample -> 
            tuple(sample[0], sample[1], pfb_file_ch, gcmodel_file_ch, plink_metadata_ch,gcDir_ch)
        }

    cnv_inputs | CNV_calling | merge_cnv_quantisnp_penncnv


    static_inputs = Channel
    .from(pfb_file_ch, gcmodel_file_ch, plink_metadata_ch, gcDir_ch)
    .map { files ->
        tuple(files[0], files[1], files[2], files[3])
    }

    // 4. Merge cnv_inputs and static_inputs and pass to the merge_cnv_quantisnp_penncnv process
    cnv_inputs.combine(static_inputs)
        .set { combined_inputs }

    combined_inputs | merge_cnv_quantisnp_penncnv
}