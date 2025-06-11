import polars as pl

# Read the input TSV file
df = pl.read_csv("microarray_cnv.tsv", separator="\t")

# Process the Copy_Number column (parse it as list of int)
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

df = df.with_columns(
    pl.col("Copy_Number")
    .list.eval(pl.element().cast(pl.Utf8))  # cast elements to string inside the list
    .list.join(",")
    .alias("Copy_Number")
)

# Reorder columns: insert 'Type' at position 5
cols = df.columns
new_order = cols[:4] + ["Type"] + [col for col in cols[4:] if col != "Type" and col != "Confidence_max"]

df = df.select(new_order)

# Save to TSV
df.write_csv("penncnv_quantisnp_cnv_merged.tsv", separator="\t")

