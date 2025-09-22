#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * This pipeline processes BAF+LRR probe files from a cohort of individuals
 * to generate the necessary input files for CNV calling: a PFB file and a GC model.
 * It includes the following steps:
 * 1. Index all BAF+LRR files.
 * 2. Identify the top 1000 samples by call rate.
 * 3. Generate a PFB file.
 * 4. Annotate SNPs with GC content using precomputed genomic windows.
 */


// Step 1: Generate a table of sample IDs and paths to BAF+LRR probe files
// process generate_list_of_path_to_BAF_LRR_Probes {
//     tag "generate_list_of_path_to_BAF_LRR_Probes"

//     input:
//     val directory_BAF_LRR_Probes_by_sample  // Input directory containing per-sample BAF+LRR files

//     output:
//     path "list_path_to_BAF_LRR_Probes.tsv" // TSV file mapping SampleID to full path

//     script:
//     """
//     echo "Process Running: generate_list_of_path_to_BAF_LRR_Probes"

//     # Create header, then find all files in the input directory and extract sample ID from filename
//     (echo -e "SampleID\tPath_to_File_BAF_LRR_Probes" &&
//     find "$directory_BAF_LRR_Probes_by_sample" -type f -name '*baf_lrr.tsv' | awk -F'/' '{
//         OFS = "\t";
//         sampleid = \$NF;
//         sub(/\\.baf_lrr\\.tsv/, "", sampleid);
//         print sampleid, \$0;
//     }' ) > list_path_to_BAF_LRR_Probes.tsv
//     """
// }




// Step 2: Select the 1000 samples with the highest call rate
process getBestSample {
    tag "get ${pfb_sample_size} best samples"

    input:
    path plink2samplemetadata_output    // File with sample ID, call rate, and imputed sex
    path list_path_to_BAF_LRR_Probes    // TSV file mapping SampleID to full path
    val pfb_sample_size     

    output:
    path "list_best_BAF_LRR_Probes.txt" // List of paths to top 1000 sample files

    script:
    """
    echo "Process Running: identify_1000_best_sampleid"
 
    # Extract top samples by call rate (assumed to be column 2), skipping header
    tail -n +2 "$plink2samplemetadata_output" | 
    awk 'NR==FNR { sample[\$1]; next } \$1 in sample' $list_path_to_BAF_LRR_Probes - | 
    sort -k2,2nr | cut -f1 | head -n ${pfb_sample_size} > SampleID_list
    
    # Filter original file list to keep only the above samples
    awk 'NR==FNR { sample[\$1]; next } \$1 in sample' SampleID_list $list_path_to_BAF_LRR_Probes | cut -f2 > list_best_BAF_LRR_Probes.txt

    # Compare the list of expected sample IDs with actual ones found (for QC/debugging)
    comm -23 <(cut -f1 SampleID_list | sort) <(cut -f1 $list_path_to_BAF_LRR_Probes | sort)
    """
}


// Step 3: Generate PFB (Population Frequency of B Allele) from selected samples
process generate_pfb {
    tag "generate_pfb"

    input:
    path list_best_BAF_LRR_Probes    // List of paths to top 1000 sample files

    output:
    path 'pfb.tsv'      // PFB file (Population Frequency of B Allele)

    script:
    """
    compile_pfb.py  ${list_best_BAF_LRR_Probes} pfb.tsv '${(task.memory.toMega() * 0.75) as int}MB'
    """
    }





// Step 4: Create a GC model file by mapping GC content to SNPs using genomic windows
process generate_gcmodel {
    tag "generate_gcmodel"


    input:
    path gc_content_windows // GC content by genomic window (e.g., from precomputed genome-wide scan)
    path pfb_file           // PFB file (Population Frequency of B Allele)

    output:
    path 'gcModel.tsv'      // Output file with GC content per SNP

    script:
    """
    echo "Process Running: generate_gcmodel"

    # Remove the first line of the input file and create a BED file with specific columns
    tail -n +2 "$pfb_file" | awk -v OFS="\\t" '{print \$2,\$3,\$3,\$1}' > SNP.bed &&

    # Create the header for the GC model output file 
    echo -e "Name\\tChr\\tPosition\\tGC_percent" > 'gcModel.tsv' &&

    # Intersect SNPs with genomic windows and compute GC content for each SNP position
    bedtools intersect -a SNP.bed -b "$gc_content_windows" -loj | awk -v OFS="\\t" '{print \$4,\$1,\$2,(\$8*100)}' >> 'gcModel.tsv'

    linecount=\$(wc -l < gcModel.tsv)

    if [ "\$linecount" -le 1 ]; then
        echo "ERROR: gcModel is an empty file. Halting execution." >&2
        exit 1
    fi
    """
}


// Workflow definition tying everything together
workflow PREPARE_PENNCNV_INPUTS {
    take:
        // Define inputs
        list_path_to_BAF_LRR
        plink2samplemetadata_output
        gc_content_windows
        pfb_sample_size

        
    main:

        // Step 1. Index all BAF+LRR files.
        // path_list = generate_list_of_path_to_BAF_LRR_Probes(directory_BAF_LRR_Probes_by_sample)

        // Step 2. Identify the top 1000 samples by call rate.
        getBestSample(
            plink2samplemetadata_output,
            list_path_to_BAF_LRR,
            pfb_sample_size
        )

        // Step 3. Generate a PFB file.
        pfb_file = generate_pfb(
            getBestSample.out
        )

        // Step 4. Annotate SNPs with GC content using precomputed genomic windows.
        gc_model = generate_gcmodel(
            gc_content_windows,
            pfb_file
        )

    emit:
        list_path_to_BAF_LRR
        pfb_file
        gc_model
}