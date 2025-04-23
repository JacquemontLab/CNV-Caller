#!/usr/bin/env nextflow

nextflow.enable.dsl=2
// NXF_OFFLINE=true nextflow run main.nf -resume -c nextflow.config -with-trace -with-report

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
process generate_list_of_path_to_BAF_LRR_Probes {
    tag "generate_list_of_path_to_BAF_LRR_Probes"

    executor "local"

    input:
    val directory_BAF_LRR_Probes_by_sample  // Input directory containing per-sample BAF+LRR files

    output:
    path "list_path_to_BAF_LRR_Probes.tsv" // TSV file mapping SampleID to full path

    script:
    """
    echo "Process Running: generate_list_of_path_to_BAF_LRR_Probes"

    # Create header, then find all files in the input directory and extract sample ID from filename
    (echo -e "SampleID\tPath_to_File_BAF_LRR_Probes" &&
    find "$directory_BAF_LRR_Probes_by_sample" -type f | awk -F'/' '{
        OFS = "\t";
        sampleid = \$NF;
        sub(/\\.BAF_LRR_Probes\\.tsv/, "", sampleid);
        print sampleid, \$0;
    }' ) > list_path_to_BAF_LRR_Probes.tsv
    """
}




// Step 2: Select the 1000 samples with the highest call rate
process identify_1000_best_sampleid {
    tag "identify_1000_best_sampleid"

    executor "local"

    input:
    path plink_metadata                  // File with sample ID, call rate, and imputed sex
    path list_path_to_BAF_LRR_Probes     // TSV file mapping SampleID to full path

    output:
    path "list_best_BAF_LRR_Probes.txt" // List of paths to top 1000 sample files

    script:
    """
    echo "Process Running: identify_1000_best_sampleid"
 
    # Extract top 1000 samples by call rate (assumed to be column 3), skipping header
    tail -n +2 "$plink_metadata" | sort -k3,3nr | cut -f1 | head -n 1000 > SampleID_list
    
    # Filter original file list to keep only those 1000 samples
    awk 'NR==FNR { sample[\$1]; next } \$1 in sample' SampleID_list $list_path_to_BAF_LRR_Probes | cut -f2 > list_best_BAF_LRR_Probes.txt

    # Compare the list of expected sample IDs with actual ones found (for QC/debugging)
    comm -23 <(cut -f1 SampleID_list | sort) <(cut -f1 $list_path_to_BAF_LRR_Probes | sort)
    """
}


// Step 3: Generate PFB (Population Frequency of B Allele) from selected samples
process generate_pfb {
    tag "generate_pfb"

    executor "local"

    input:
    path list_best_BAF_LRR_Probes    // List of paths to top 1000 sample files

    output:
    path 'pfb.tsv'      // PFB file (Population Frequency of B Allele)

    script:
    """
    echo "Process Running: generate_pfb"

    # Run Perl script to compile the PFB file from selected samples
    compile_pfb.py "$list_best_BAF_LRR_Probes" 'pfb.tsv' 64
    """
}


// Step 4: Create a GC model file by mapping GC content to SNPs using genomic windows
process generate_gcmodel {
    tag "generate_gcmodel"

    executor "local"

    input:
    path gc_content_windows // GC content by genomic window (e.g., from precomputed genome-wide scan)
    path pfb_file           // PFB file (Population Frequency of B Allele)

    output:
    path 'gcModel.tsv'      // Output file with GC content per SNP

    script:
    """
    echo "Process Running: generate_gcmodel"

    # Remove the first line of the input file and create a BED file with specific columns
    tail -n +2 "$pfb_file" | awk '{print \$2"\t"\$3"\t"\$3"\t"\$1}' > SNP.bed

    # Create the header for the GC model output file
    echo -e "Name\tChr\tPosition\tGC_percent" > 'gcModel.tsv'

    # Intersect SNPs with genomic windows and compute GC content for each SNP position
    bedtools intersect -a SNP.bed -b "$gc_content_windows" -loj | awk -F"\t" '{print \$4"\t"\$1"\t"\$2"\t"(\$8*100)}' >> 'gcModel.tsv'
    """
}



// Workflow definition tying everything together
workflow {
    // Input parameters (define in config or command line)
    directory_BAF_LRR_Probes_by_sample = file("/home/flben/projects/rrg-jacquese/All_user_common_folder/RAW_DATA/Genetic/ALSPAC/BAF_LRR_Probes_by_sample/")
    plink_metadata = file("/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/data_from_plink/work/76/e35544f2e885e3217ccf31b84e5ef2//plink_metadata.tsv")
    gc_content_windows = file("/home/flben/projects/rrg-jacquese/flben/cnv_annotation/scripts/workflow/CNV-Annotation-pipeline/modules/penncnv_params/resources/gc_content_1k_windows_GRCh37.bed")

    // Step 1
    generate_list_of_path_to_BAF_LRR_Probes(directory_BAF_LRR_Probes_by_sample)

    // Step 2
    identify_1000_best_sampleid(
        plink_metadata,
        generate_list_of_path_to_BAF_LRR_Probes.out
    )

    // Step 3
    generate_pfb(
        identify_1000_best_sampleid.out
    )

    // Step 4
    generate_gcmodel(
        gc_content_windows,
        generate_pfb.out
    )
}