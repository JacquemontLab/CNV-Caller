#!/usr/bin/env python3
import argparse
import polars as pl
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Process CNV TSV file from QuantiSNP and classify CNVs as DUP, DEL, or MIX."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Path to the input TSV file"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Path to the output TSV file"
    )
    args = parser.parse_args()

    try:
        # Read the input TSV file
        df = pl.read_csv(
            args.input, 
            separator="\t", 
            infer_schema_length=1_000_000,
            schema_overrides={"Copy_Number": pl.Utf8,"SampleID": pl.Utf8}
        )

        # Process the Copy_Number column (parse as list of int)
        df = df.with_columns(
            (
                pl.col("Copy_Number")
                .str.split(",")
                .list.eval(pl.element().cast(pl.Int32))
            ).alias("Copy_Number")
        )

        # Classify CN type
        df = df.with_columns(
            pl.when(
                (pl.col("Copy_Number").list.max() >= 2) &
                (pl.col("Copy_Number").list.min() >= 2)
            )
            .then(pl.lit("DUP"))
            .when(
                (pl.col("Copy_Number").list.max() < 2) &
                (pl.col("Copy_Number").list.min() < 2)
            )
            .then(pl.lit("DEL"))
            .otherwise(pl.lit("MIX"))
            .alias("Type")
        )

        # Convert Copy_Number back to string list
        df = df.with_columns(
            pl.col("Copy_Number")
            .list.eval(pl.element().cast(pl.Utf8))
            .list.join(",")
            .alias("Copy_Number")
        )

        # Reorder columns: insert 'Type' at position 5
        cols = df.columns
        new_order = cols[:4] + ["Type"] + [col for col in cols[4:] if col != "Type" ]
        df = df.select(new_order)

        # Save to TSV
        df.write_csv(args.output, separator="\t")

    except Exception as e:
        print(f"[✘] Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
