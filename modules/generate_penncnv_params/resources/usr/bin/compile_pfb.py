#!/usr/bin/env python

# ------------------------------------------------------------------------------
# Author: Benjamin Clark
# Date: 16-05-2025
# Description: This script processes a list of SNP files, computes the 
#              PFB (Proportional Frequency of B Allele) for each SNP across 
#              multiple files.
#
# Warning:
# - This script processes large datasets and may require substantial memory.
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
#
# Example:
# python compile_pfb_duckdb.py snp_file_list.txt output_results.tsv memory_limit
# ------------------------------------------------------------------------------

import polars as pl
import duckdb
import sys
import os
import psutil
import logging

def main():
    if len(sys.argv) < 2:
        print("Usage: python compile_pfb_duckdb.py <listfile> [output] [memory_limit]", file=sys.stderr)
        sys.exit(1)
        
    listfile = sys.argv[1]  # Required: listfile
    
    output = sys.argv[2] if len(sys.argv) >= 3 else "pfb.tsv"  # Default output file

    # Handle optional third argument for memory limit
    if len(sys.argv) >= 4:
        memory_limit = sys.argv[3]
    else:
        # Get 90% of total available system memory in bytes
        available_mem_bytes = psutil.virtual_memory().available ###!!broken!! 
        memory_limit_bytes = int(available_mem_bytes * 0.9)

        # Convert to gigabytes and format as a string for DuckDB
        memory_limit = f"{memory_limit_bytes // (1024**3)}GB"

    print(f"listfile: {listfile}")
    print(f"output: {output}")
    print(f"memory_limit: {memory_limit}")

    #turning input file into list
    if listfile:
        with open(listfile) as f:
            input_files = [line.strip() for line in f if line.strip()]

     # Check if any input files were provided
    if not input_files:
        print("Error: No input files provided.", file=sys.stderr)
        sys.exit(1)


    print(f"Number of input files to process: {len(input_files)}")
    
    #stream to parquet
    (pl.scan_csv(input_files, separator="\t",
                              has_header=False,
                              skip_rows=1,
                              new_columns=["Name","Chr","Position","LRR", "BAF"],
                              schema_overrides=[pl.String, pl.Categorical, pl.Int64, pl.Float64, pl.Float64])
                              .drop_nans()
                              .sink_parquet(pl.PartitionMaxSize(
                                                "./pfb_parquet/",
                                                max_size = 10_000_000
                                                ),
                                            mkdir = True,
                                            maintain_order = False, 
                                            compression = 'snappy',
                                            row_group_size = 2_500_000))
    logging.info("Parquet built")

    #Execute groupby, taking the average BAF per group of SNPs across samples
    duckdb.sql(
        f"""
        SET memory_limit = '{memory_limit}';
        SET preserve_insertion_order = false;
        SET temp_directory = './tmp_dir.tmp/';

            COPY (
                SELECT 
                    Name,
                    Chr,
                    Position,
                    ROUND(AVG(BAF), 3) AS PFB
                FROM read_parquet('pfb_parquet/*.parquet', hive_partitioning = true)
                GROUP BY Chr, Position, Name
            ) TO '{output}' (FORMAT 'csv', DELIMITER '\t', HEADER)
        """
    )

    #cleanup
    #os.remove("compiled_baf.parquet")

if __name__ == "__main__":
    main()
