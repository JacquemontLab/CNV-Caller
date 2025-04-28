#!/usr/bin/env python
 
# Importing necessary libraries
import pandas as pd  # For reading and manipulating data
import numpy as np  # For numerical operations (though not used here directly)
import os  # For interacting with the file system
import sys  # For handling system-specific parameters and inputs
from tqdm import tqdm  # For displaying a progress bar during iterations
from collections import defaultdict  # For creating dictionaries with default values
import multiprocessing as mp  # For parallel processing

# ------------------------------------------------------------------------------
# Script: SNP B Allele Frequency (BAF) Processor
# Author: Florian Bénitière
# Date: 17-04-2025
# Description: This script processes a list of SNP files, computes the sum and
#              count of B Allele Frequencies (BAF), and calculates the 
#              PFB (Proportional Frequency of B Allele) for each SNP across 
#              multiple files using parallel processing.
#              The results are merged, and the final output is written 
#              to a specified output file or printed to standard output.
# Warning:
# - This script processes large datasets and may require substantial memory.
#   Ensure your system has sufficient memory to handle large files, especially
#   when processing many files or files with a large number of SNPs. For memory
#   optimization, consider processing the data in smaller chunks or limiting
#   the number of parallel processes. You can track memory by using /usr/bin/time -v bash -c ''

#
# Input:
# - listfile: A text file containing a list of paths to SNP signal intensity files (one per line).
#   Each file should have columns like:
#     Name    Chr    Position    XXSampleIDXX.Log R Ratio    XXSampleIDXX.B Allele Freq
#   Example:
#     rs217682    14    62356724    -0.0639626951149    0.492056429742
#
# - output: Path to the output file where results will be saved. If not provided,
#           the output will be printed to the standard output.
# - num_cpu: Number of CPU cores to use for parallel processing.
#
# Example:
# python compile_pfb.py snp_file_list.txt output_results.tsv 4
# ------------------------------------------------------------------------------



# Function to read the SNP position file and return dictionaries for chromosomes and positions
def read_snpposfile(snpposfile):
    df = pd.read_csv(snpposfile, sep='\t', dtype=str, usecols=["Name", "Chr", "Position"])
    snp_chr = df.set_index("Name")["Chr"].to_dict()
    snp_pos = df.set_index("Name")["Position"].to_dict()
    return snp_chr, snp_pos

# Function to process an individual file and compute the sum and count of B Allele Frequencies (BAF)
def process_file(filepath):
    local_sum = defaultdict(float)
    local_count = defaultdict(int)
    try:
        df = pd.read_csv(filepath, sep='\t', dtype=str)
        baf_col = next(col for col in df.columns if "B Allele Freq" in col)
        df = df[["Name", baf_col]]
        df = df[df[baf_col].notnull() & ~df[baf_col].isin(["NA", "NaN"])]
        df[baf_col] = pd.to_numeric(df[baf_col], errors='coerce')
        df = df.dropna()
        for row in df.itertuples(index=False):
            name, baf = row
            local_sum[name] += baf
            local_count[name] += 1
    except Exception as e:
        print(f"Error processing {filepath}: {e}", file=sys.stderr)
    return local_sum, local_count

# Function to merge the partial results from each file
def merge_results(results):
    global_sum = defaultdict(float)
    global_count = defaultdict(int)

    for local_sum, local_count in tqdm(results, desc="Merging BAF", unit="file"):
        for name in local_sum:
            global_sum[name] += local_sum[name]
            global_count[name] += local_count[name]

    return global_sum, global_count

# Function to merge the partial results from each file
def main():# Read arguments
    listfile = sys.argv[1]  # First argument: listfile
    output = sys.argv[2]  # Second argument: output file
    num_cpu = int(sys.argv[3]) if len(sys.argv) > 3 else 1  # Set num_cpu to 1 by default if not provided

    print(f"Arg1: {listfile}")
    print(f"Arg2: {output}")
    print(f"Arg3: {num_cpu}")

    # Read the input files from the list file
    if listfile:
        with open(listfile) as f:
            inputfiles = [line.strip() for line in f if line.strip()]

    # Check if any input files were provided
    if not inputfiles:
        print("Error: No input files provided.", file=sys.stderr)
        sys.exit(1)
    
    # Read the SNP position file from the first input file
    snp_chr, snp_pos = read_snpposfile(inputfiles[0])

    # Use a multiprocessing pool to process the files in parallel
    with mp.Pool(processes=num_cpu) as pool:
        results = list(tqdm(pool.imap_unordered(process_file, inputfiles), total=len(inputfiles), desc="Collecting BAF per file"))

    # Combine all the partial sums and counts
    baf_sum, baf_count = merge_results(results)

    # Prepare the data for the final output
    names = sorted(list(snp_chr.keys()))
    pfb_data = {
        "Name": names,
        "Chr": [snp_chr.get(name, "NA") for name in names],
        "Position": [snp_pos.get(name, "NA") for name in names],
        "PFB": [round(baf_sum[name] / baf_count[name], 3) for name in names]
    }

    # Create a DataFrame with the computed data
    result_df = pd.DataFrame(pfb_data)

    # If SNP mapping is available, sort the results based on SNP names
    if snp_chr:
        result_df["Name"] = pd.Categorical(result_df["Name"], categories=list(snp_chr.keys()), ordered=True)
        result_df = result_df.sort_values("Name").reset_index(drop=True)
        
    # Output the final results to a file or stdout
    if output:
        result_df.to_csv(output, sep='\t', index=False)
    else:
        print("Error: No output file provided. Please specify an output file.", file=sys.stderr)
        sys.exit(1)


# Ensure the main function runs when the script is executed directly
if __name__ == "__main__":
    main()
