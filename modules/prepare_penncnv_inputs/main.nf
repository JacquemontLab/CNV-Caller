#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * This pipeline processes BAF+LRR probe files from a cohort of individuals
 * to generate the necessary input files for CNV calling: a PFB file and a GC model.
 * Steps:
 * 1. Identify the top X samples by call rate.
 * 2. Generate a PFB file.
 * 3. Generate an HMM from top samples.
 * 4. Annotate SNPs with GC content using precomputed genomic windows.
 */


// Step 1: Select the X samples with the highest call rate
process getBestSample {
    tag "get ${pfb_max_sample_size} best samples"

    input:
    path plink2samplemetadata_output    // File with sample ID, call rate, and imputed sex
    path list_path_to_BAF_LRR_Probes    // TSV file mapping SampleID to full path
    val pfb_max_sample_size     

    output:
    path "list_best_BAF_LRR_Probes.txt" // List of paths to top X sample files

    script:
    """
    echo "Process Running: identify_${pfb_max_sample_size}_best_sampleid"
 
    # Extract top samples by call rate (assumed to be column 2), skipping header
    tail -n +2 "$plink2samplemetadata_output" | 
    awk 'NR==FNR { sample[\$1]; next } \$1 in sample' $list_path_to_BAF_LRR_Probes - | 
    sort -k2,2nr | cut -f1 | head -n ${pfb_max_sample_size} > SampleID_list
    
    # Filter original file list to keep only the above samples
    awk 'NR==FNR { sample[\$1]; next } \$1 in sample' SampleID_list $list_path_to_BAF_LRR_Probes | cut -f2 > list_best_BAF_LRR_Probes.txt

    # Compare the list of expected sample IDs with actual ones found (for QC/debugging)
    comm -23 <(cut -f1 SampleID_list | sort) <(cut -f1 $list_path_to_BAF_LRR_Probes | sort)
    """
}


// Step 2: Generate PFB (Population Frequency of B Allele) from selected samples
process generate_pfb {


    input:
    path list_best_BAF_LRR_Probes    // List of paths to top X sample files

    output:
    path 'pfb.tsv'      // PFB file (Population Frequency of B Allele)

    script:
    """
    echo ${list_best_BAF_LRR_Probes} | tr " " "\\n" | tr -d " " > batch_list.txt
    compile_pfb.py  batch_list.txt pfb.tsv '${(task.memory.toMega() * 0.75) as int}MB'
    """
    }




// Step 3: Generate HMM from the first 10 best samples
process generate_hmm {
    tag "building hmm model from default ${data_type} model"
    label "penncnv_quantisnp"

    input:
    path list_best_BAF_LRR_Probes    // List of paths to top X sample files
    path pfb_file                    // PFB file (Population Frequency of B Allele)
    val data_type

    output:
    path 'hmm_trained.hmm'      // HMM file

    script:
    """
    # take first 10 lines (sample paths) from the list
    # head -n 10 "$list_best_BAF_LRR_Probes" > list_baf_lrr.txt 
    echo ${list_best_BAF_LRR_Probes} | cut -d " " -f1-10 | tr " " "\\n" | tr -d " " > list_baf_lrr.txt

    # We start from a pre-existing HMM:
    #   - hhall.hmm  -> general high-density SNP arrays
    #   - hh550.hmm  -> Illumina HumanHap550-specific
    #   - wgs.hmm    -> whole-genome sequencing
    #
    # Each HMM differs in emission and transition probabilities because:
    #   - SNP arrays have fewer probes and noisier signals (hhall, hh550)
    #   - WGS has dense, uniform coverage (wgs)
    # Using the correct starting HMM improves CNV calling accuracy.

    if [[ "${data_type}" == "array"  ]];then
        model="hhall.hmm"
    elif [[ "${data_type}"  == "wgs" ]];then
        model="wgs.hmm"
    else
        echo "Invalid datatype. Use either 'array' or 'wgs'."
        exit 1 
    fi
    
    /usr/local/bin/timedev penncnv --train \
        --hmmfile "/usr/local/PennCNV-1.0.5/lib/\$model"  \
        --pfbfile "$pfb_file" \
        --listfile list_baf_lrr.txt \
        --output hmm_trained
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
        pfb_max_sample_size
        data_type

        
    main:
        // Step 1. Identify the top X samples by call rate.
        getBestSample(
            plink2samplemetadata_output,
            list_path_to_BAF_LRR,
            pfb_max_sample_size
        )

        best_sample_ch = getBestSample.out.splitCsv().collect()

        // Step 2. Generate a PFB file.
        pfb_file = generate_pfb(
            best_sample_ch 
        )

        // Step 3. Generate HMM from the top samples.
        hmm_file = generate_hmm(
            best_sample_ch ,
            pfb_file,
            data_type
        )

        // Step 4. Annotate SNPs with GC content using precomputed genomic windows.
        gc_model = generate_gcmodel(
            gc_content_windows,
            pfb_file
        )

    emit:
        list_path_to_BAF_LRR
        pfb_file
        hmm_file
        gc_model
}