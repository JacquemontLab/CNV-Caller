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


def main():
    listfile = sys.argv[1]  # First argument: listfile
    output = sys.argv[2]  # Second argument: output file
    memory_limit = sys.argv[3]  # Third argument: memory limit, example= 32GB

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


    #stream to parquet
    (pl.scan_csv(input_files, separator="\t",
                              has_header=False,
                              skip_rows=1,
                              new_columns=["Name","Chr","Position","LRR", "BAF"],
                              schema_overrides=[pl.String, pl.String, pl.Int32, pl.Float32, pl.Float32])
                              .drop_nans()
                              .sink_parquet("compiled_1k.parquet"))
    print("Built Parquet")

    #Execute groupby, taking the average BAF per group of SNPs across samples
    duckdb.sql(
        f"""
        SET memory_limit = '{memory_limit}';
        COPY (
            SELECT 
                Name,
                Chr,
                Position,
                ROUND(SUM(BAF) / COUNT(*), 3) AS PFB
            FROM read_parquet('compiled_1k.parquet')
            GROUP BY Name, Chr, Position
        ) TO '{output}' (FORMAT 'csv', DELIMITER '\t', HEADER)
        """
    )

    #cleanup
    os.remove("compiled_1k.parquet")

if __name__ == "__main__":
    main()
