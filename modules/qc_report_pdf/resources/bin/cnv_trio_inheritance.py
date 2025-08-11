#!/usr/bin/env python3
# =============================================================================
# CNV Inheritance Analysis Pipeline
# Author: Florian Bénitière and Benjamin Clark
# Date: 22/07/2025
# Implementation: Python (Polars)
#
# Description:
# This pipeline analyzes the inheritance of CNVs (Copy Number Variants) 
# in trio datasets. It filters CNVs based on pedigree information, 
# constructs parent–child mappings, performs overlap checks with bedtools, 
# and outputs an annotated child CNV table.
#
# Requirements:
# - Python packages: polars, pyarrow
# - External tools: zcat, awk, bedtools
#
# Input:
# - CNV file (.tsv.gz or .parquet file)
# - Pedigree file (tab-separated: SampleID, FatherID, MotherID)
#
# Output:
# - Annotated child CNV BED and TSV files
# =============================================================================


import polars as pl
import os
import argparse
from pathlib import Path
import subprocess
pl.enable_string_cache()

# ------------------- Parse Arguments -------------------

parser = argparse.ArgumentParser(
    description="CNV Inheritance Analysis Pipeline (Polars implementation)"
)
parser.add_argument(
    "-c", "--cnv", required=True,
    help="Path to merged CNV file (.tsv.gz OR .parquet) (SampleID, Chr, Start, End)"
)
parser.add_argument(
    "-p", "--pedigree", required=True,
    help="Path to pedigree file (.tsv) (SampleID, FatherID, MotherID)"
)
parser.add_argument(
    "-o", "--output", default="annotated_child_cnv.tsv",
    help="Path to the output child CNVs"
)
parser.add_argument(
    "-t", "--type_col", default=None,
    help="Column name of CNV type (DEL, DUP, MIX etc...) for per-type annotation; leave empty to skip"
)
parser.add_argument(
    "-f", "--overlap", default="0.5",
    help="Comma-separated list of reciprocal overlap fractions (e.g., 0.5,0.1)"
)

args = parser.parse_args()

# +
FILE_CNV = args.cnv
PEDIGREE_FILE = args.pedigree
OUTPUT = args.output
TYPE_COL = args.type_col
OVERLAPS = [float(x) for x in args.overlap.split(",")]

#Detect input type naively
tsv_mode = True
suffix = Path(FILE_CNV).suffixes
if '.parquet' in suffix:
    tsv_mode = False
elif ".tsv" not in suffix:
    print("Invalid input type or malformed extension")
    sys.exit(1)
# -

# ------------------- Filter CNV Calls -------------------

# +
family_info = pl.scan_csv(PEDIGREE_FILE, separator="\t",infer_schema_length=100000)

#load CNV lazily
if TYPE_COL:
    if tsv_mode:
        df_cnv = pl.scan_csv(FILE_CNV, 
                             separator = "\t", 
                             schema_overrides={f"{TYPE_COL}": pl.Categorical, 
                                           "Chr": pl.Categorical},infer_schema_length=100000)
        print(f"Loading {FILE_CNV} in tsv mode and using {TYPE_COL} as type column")
    #casting to categorical for performance 
    else:
        df_cnv = pl.scan_parquet(FILE_CNV)
        df_cnv = df_cnv.with_columns(pl.col(f"{TYPE_COL}").cast(pl.Categorical),
                                     pl.col("Chr").cast(pl.Categorical))
        print(f"Loading {FILE_CNV} in parquet mode and using {TYPE_COL} as type column")
else:
    if tsv_mode:
        print(f"Loading {FILE_CNV} in tsv mode without type column")
        df_cnv = pl.scan_csv(FILE_CNV, 
                             separator = "\t", 
                             schema_overrides={"Chr": pl.Categorical},infer_schema_length=100000)
    else:
        print(f"Loading {FILE_CNV} in parquet mode and without type column")
        df_cnv = pl.scan_parquet(FILE_CNV)
        df_cnv = df_cnv.with_columns(pl.col("Chr").cast(pl.Categorical))
        

print("Filtering for trios that are completely found in the cnv table...")

#keep only trios that have CNVs for every member 
sample_ids = df_cnv.select(pl.col("SampleID")).unique().collect().to_series()                         
family_trios =  family_info.filter(
    pl.col("SampleID").is_in(sample_ids.implode()) & #adding implodes as per https://github.com/pola-rs/polars/pull/22178
    pl.col("FatherID").is_in(sample_ids.implode()) &
    pl.col("MotherID").is_in(sample_ids.implode())
)



# +
print("Building bed files...")

#child cnvs subset from those available from pedigree file 
df_child = df_cnv.join(family_trios, on = "SampleID", how = "inner")
#father cnvs 
df_father = df_cnv.join(family_trios, left_on = "SampleID", 
                                     right_on = "FatherID", 
                                     suffix = ".child", 
                                     how = "inner")
#mother cnvs
df_mother = df_cnv.join(family_trios, left_on = "SampleID", 
                                     right_on = "MotherID", 
                                     suffix = ".child", 
                                     how = "inner")


#reorder to put child ID to front and drop parent IDs
standard_col = ["SampleID.child","Chr","Start","End"] 
f_order = standard_col + [col for col in df_father.collect_schema().names() if col not in standard_col]
m_order = standard_col + [col for col in df_mother.collect_schema().names() if col not in standard_col]
df_parents = pl.concat([df_father.select(f_order).drop("MotherID","SampleID"), 
                        df_mother.select(m_order).drop("FatherID", "SampleID")])

#build ID column for bedtools with concatentated groups (ie: Chr and type)
if TYPE_COL is not None:

    parents_out = df_parents.with_columns(pl.concat_str([pl.col("SampleID.child"),
                                                         pl.lit("_"),
                                                         pl.col("Chr"), 
                                                         pl.lit("_"), 
                                                         pl.col(f"{TYPE_COL}")])
                                          .alias("SampleID")).select("SampleID","Start","End")

    child_out = df_child.with_columns(pl.concat_str([pl.col("SampleID"),
                                                     pl.lit("_"),
                                                     pl.col("Chr"), 
                                                     pl.lit("_"), 
                                                     pl.col(f"{TYPE_COL}")])
                                          .alias("SampleID")).select("SampleID","Start","End")

else:

    parents_out = df_parents.with_columns(pl.concat_str([pl.col("SampleID.child"),
                                                         pl.lit("_"),
                                                         pl.col("Chr"), 
                                                         pl.lit("_")
                                                        ])
                                          .alias("SampleID")).select("SampleID","Start","End")

    child_out = df_child.with_columns(pl.concat_str([pl.col("SampleID"),
                                                     pl.lit("_"),
                                                     pl.col("Chr"), 
                                                     pl.lit("_") 
                                                     ])
                                      .alias("SampleID")).select("SampleID","Start","End")

#print out as bed files
child_bed = f"child_cnv.bed"
parent_bed = f"parents_cnv.bed"

child_out.sink_csv(child_bed, separator = "\t", include_header = False)
parents_out.sink_csv(parent_bed, separator = "\t", include_header = False)


# -

# ------------------- Run bedtools intersect -------------------

for ovlap in OVERLAPS:
    bash_cmd = f"bedtools intersect -a {child_bed} -b {parent_bed} -f {ovlap} -r -wa -u > intersect_ovlap{ovlap}.bed;\
bedtools intersect -a {child_bed} -b {parent_bed} -f {ovlap} -r -v > non_intersect_ovlap{ovlap}.bed"
    print("Running command:\n", bash_cmd)
    subprocess.run(bash_cmd, shell=True, check=True)

# +
print("Formatting intersections and combining results...")

#load non-intersects and intersects and add the overlap columns, all false for non-overlap and vice-versa
ovl_count = len(OVERLAPS)
count = 0
#for multiple overlap attempts we need to build a single table using multiple joins
while count < ovl_count:
    #for the first file no join needed 
    if(count == 0):
        non_intersects_file = pl.scan_csv(f"non_intersect_ovlap{OVERLAPS[count]}.bed", 
                                          separator = "\t", has_header = False,infer_schema_length=100000)
        #all false from non-interesects
        non_intersects = non_intersects_file.with_columns(pl.lit(False)
                                                          .alias(f"Observed_in_Parent_{OVERLAPS[count]}"))
        #all true from intersects
        intersects_file = pl.scan_csv(f"intersect_ovlap{OVERLAPS[count]}.bed", separator = "\t", has_header = False,infer_schema_length=100000)
        intersects = intersects_file.with_columns(pl.lit(True).alias(f"Observed_in_Parent_{OVERLAPS[count]}"))
        
        count += 1
    else:
        #non-intersects
        non_intersects_file = pl.scan_csv(f"non_intersect_ovlap{OVERLAPS[count]}.bed", 
                                          separator = "\t", has_header = False,infer_schema_length=100000)
        
        non_intersects_file = non_intersects_file.with_columns(pl.lit(False)
                                                          .alias(f"Observed_in_Parent_{OVERLAPS[count]}"))
        #join on all matching columns to add new intersect column
        non_intersects = non_intersects.join(non_intersects_file, on = ["column_1","column_2","column_3"])

        
        #intersects file
        intersects_file = pl.scan_csv(f"intersect_ovlap{OVERLAPS[count]}.bed", 
                                      separator = "\t", has_header = False,infer_schema_length=100000)
        
        intersects_file = intersects_file.with_columns(pl.lit(True).alias(f"Observed_in_Parent_{OVERLAPS[count]}"))
        intersects = intersects.join(intersects_file, on = ["column_1","column_2","column_3"])
        count += 1
        


#concat the two files
lookup = pl.concat([intersects, non_intersects])


lookup_id_fields  = ["SampleID", "Chr"]
new_lookup_cols = ["SampleID", "Chr", "Start", "End"]

if TYPE_COL:
    lookup_id_fields += [f"{TYPE_COL}"]
    new_lookup_cols += [f"{TYPE_COL}"]

#print(new_lookup_cols)
#reformat and deconstruct the id field
lookup = (lookup.with_columns(
            pl.col("column_1")
            .str.split_exact("_", 2)
            .struct.rename_fields(lookup_id_fields)
            .alias("ID_fields")
    ).unnest("ID_fields").drop("column_1")
)


lookup_dict = {"column_2" : "Start", "column_3" : "End"}
lookup = lookup.rename(lookup_dict)



if TYPE_COL:
    #cast back to make the next join work
    lookup = lookup.with_columns(pl.col("Chr").cast(pl.Categorical),
                                 pl.col(f"{TYPE_COL}").cast(pl.Categorical))
else:
    lookup = lookup.with_columns(pl.col("Chr").cast(pl.Categorical))

#Make sure the SampleID in the Child CNVs df is string to match with original child cnv dataframe
df_child  = df_child.with_columns(pl.col("SampleID").cast(pl.Utf8))


#perform join on child CNV columns 
result = df_child.join(lookup, on = new_lookup_cols,
                               how = "inner",
                               coalesce=True)


#reorder for final output 
standard_col = ["SampleID","Chr","Start","End"] 
f_order = standard_col + [col for col in result.collect_schema().names() if col not in standard_col]

#sink
# result.select(f_order).drop("FatherID", "MotherID").sink_csv(f"{OUTPUT}", separator = "\t")
# Collect pedigree file column names
pedigree_cols = family_info.collect_schema().names()

# Keep only SampleID from pedigree
cols_to_drop = [col for col in pedigree_cols if col != "SampleID"]
print(cols_to_drop)
# Drop them if they exist in the final DataFrame
lookup = result.select(f_order).drop(cols_to_drop, strict=False).sink_csv(f"{OUTPUT}", separator = "\t")


print("Pipeline completed successfully.")


print(f"💾 Results written to {OUTPUT}")


# --------------------- Cleanup Temporary Files ---------------------
files_to_remove = [
    parent_bed,
    child_bed
] + [f"intersect_ovlap{ov}.bed" for ov in OVERLAPS] + [f"non_intersect_ovlap{ov}.bed" for ov in OVERLAPS]

for f in files_to_remove:
    try:
        os.remove(f)
        print(f"🗑️ Removed {f}")
    except FileNotFoundError:
        print(f"⚠️ Skipped {f} (not found)")
