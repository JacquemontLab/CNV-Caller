#!/usr/bin/env python

import polars as pl
import duckdb
import sys
import os


def main():
    listfile = sys.argv[1]  # First argument: listfile
    output = sys.argv[2]  # Second argument: output file

    print(f"Arg1: {listfile}")
    print(f"Arg2: {output}")

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
        SET memory_limit = '32GB';
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
#testing polars implementation, currently uses too much RAM
# (pl.scan_parquet("compiled_1k.parquet")
#                 .group_by("Name", "Chr", "Position")
#                 .agg(PFB = pl.col("BAF").mean())
#                 .sink_csv(output, separator="\t"))

if __name__ == "__main__":
    main()
