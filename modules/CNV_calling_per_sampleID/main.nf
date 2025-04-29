#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * This Nextflow pipeline performs Copy Number Variation (CNV) calling on selected samples using PennCNV and QuantiSNP, 
 * merges the results, and extracts quality control (QC) metrics. It filters input samples from a provided list, runs 
 * CNV detection on each sample, and organizes the outputs for downstream analysis.
 */

// NXF_OFFLINE=true nextflow run main.nf -resume

// Define the process to run CNV calling
process CNV_calling {
    tag "CNV_calling: ${sample_id}"

    input:
    tuple val(sample_id), path(BAF_LRR_Probes), path(pfb_file), path(gcmodel_file), path(from_plink_extracted_data), path(gcdir)

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
    zgrep -P "${sample_id}\t" ${from_plink_extracted_data} | cut -f1,2 > sexfile.tsv
    gender=\$(cut -f2 sexfile.tsv)
    if [[ -z "\$gender" ]]; then
        echo "ERROR: No gender information found for sample ID ${sample_id}."
        exit 1
    fi

    echo "Gender: \$gender"
    if [[ "\$gender" != "male" && "\$gender" != "female" ]]; then
        echo "ERROR: Invalid gender '\$gender' for sample ID ${sample_id}. Expected 'male' or 'female'."
        exit 1
    fi


    echo "Running PennCNV on autosomes..."
    /usr/bin/time -v perl /usr/local/bin/detect_cnv.pl --test \
        --conf \
        --pfbfile ${pfb_file} \
        --hmmfile \${hmm_file} \
        --logfile ${sample_id}.penncnv.log \
        --output ${sample_id}.penncnv.out \
        --gcmodelfile ${gcmodel_file} \
        ${BAF_LRR_Probes}
        
        
    echo "Running PennCNV on ChrX..."
    /usr/bin/time -v perl /usr/local/bin/detect_cnv.pl --test \
        --conf \
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


// Merge CNV callers results and extract QC metrics
process merge_cnv_callers_and_extract_qc {
    tag "merge_cnv_callers_and_extract_qc: ${sample_id}"

    input:
    tuple val(sample_id), file(BAF_LRR_Probes), file(quantisnp_file), file(penncnv_file), file(penncnv_logfile), val(genome_version), file(regions_file)

    output:
    path "${sample_id}.PennCNV_QC.tsv", emit: PennCNV_QC_tsv
    path "${sample_id}.CNVs.tsv", emit: CNVs_tsv

    script:
    """
    echo "Process Running: merge_cnv_callers_and_extract_qc for ${sample_id}"

    merge_cnv_quantisnp_penncnv.sh "$quantisnp_file" "$penncnv_file" "$BAF_LRR_Probes" "$regions_file" "$genome_version" ${sample_id}.CNVs.tsv

    extract_qc_penncnv.sh "$penncnv_logfile" ${sample_id}.PennCNV_QC.tsv
    """
}



workflow CNV_CALLING {
    take:
        cnv_inputs
        // genome_version
        // regions_file
        // sample_inputs
        // pfb_file
        // gcmodel_file
        // from_plink_extracted_data
        // gcDir
        // genome_version
        // regions_file

    main:
        // Map inputs to full tuple required by CNV_calling
        // cnv_inputs = sample_inputs
        //     .map { sample -> 
        //         tuple(sample[0], sample[1], pfb_file, gcmodel_file, from_plink_extracted_data, gcDir)
        //     }


        // Chain the two processes:
        cnv_inputs | CNV_calling
        //  | map { sample_id, BAF_LRR_Probes, quantisnp_file, penncnv_file, penncnv_logfile -> 
        //     tuple(sample_id, BAF_LRR_Probes, quantisnp_file, penncnv_file, penncnv_logfile, genome_version, regions_file)
        // } | merge_cnv_callers_and_extract_qc

    emit:
        output = CNV_calling.out
}




workflow {

    // Define all the input parameters manually
    list_path_to_BAF_LRR_Probes = "/lustre06/project/6008022/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/penncnv_params/work/dc/424e5e2256e1e965f0c3a95a91fb15/list_path_to_BAF_LRR_Probes.tsv"
    pfb_file = "/lustre06/project/6008022/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/penncnv_params/work/5d/2a4f9ea2fba78f07dd2385b1925729/pfb.tsv"
    gcmodel_file = "/lustre06/project/6008022/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/penncnv_params/work/2f/7c956bea5daf68b4ef3b8cb534b1a2/gcModel.tsv"
    from_plink_extracted_data = "/lustre06/project/6008022/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/data_from_plink/work/38/c40f4c897d9c2da4f9b98149b81a9a/from_plink_extracted_data.tsv"
    regions_file = "/lustre06/project/6008022/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/per_sampleID/resources/Genome_Regions_data.tsv"
    genome_version = "GRCh37"
    gcDir = "/lustre06/project/6008022/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/penncnv_params/resources/" + params.genome_version + "_GCdir/"

    // List of sample IDs you want to keep
    def selected_sample_ids = ["ALSPAC09897249", "ALSPAC09902309"]

    // 1. Input channel from TSV sample_inputs dictionnary : {SampleID:Path_to_File_BAF_LRR_Probes}
    Channel
        .fromPath(list_path_to_BAF_LRR_Probes)
        .splitCsv(header: true, sep: '\t')
        .filter { row ->                                 
            row.SampleID in selected_sample_ids         // Keep only rows with desired SampleIDs
        }
        .map { row -> 
            [row.SampleID, file(row.Path_to_File_BAF_LRR_Probes)]
        }
        .set { sample_inputs }

    // (Optional) view the inputs - useful for debugging
    sample_inputs.view()

    // Map inputs to full tuple required by CNV_calling
    cnv_inputs = sample_inputs
        .map { sample -> 
            tuple(sample[0], sample[1], file(pfb_file), file(gcmodel_file), file(from_plink_extracted_data), file(gcDir))
        }

    // Chain the two processes:
    cnv_inputs | CNV_calling | map { sample_id, BAF_LRR_Probes, quantisnp_file, penncnv_file, penncnv_logfile -> 
        tuple(sample_id, BAF_LRR_Probes, quantisnp_file, penncnv_file, penncnv_logfile, genome_version, file(regions_file))
    } | merge_cnv_callers_and_extract_qc

}