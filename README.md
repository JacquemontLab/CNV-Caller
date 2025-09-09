[![Jacquemont's Lab Header](labheader.png)](https://www.jacquemont-lab.org/)

[Git Repository CNV-Caller](https://github.com/JacquemontLab/CNV-Caller)

# Documentation of the CNV Calling Pipeline for Microarray
This repository contains a bioinformatics pipeline for the discovery of copy number variants (CNVs).
The workflow is implemented using Nextflow to ensure reproducibility and efficient execution.

## Pipeline Overview
This repository provides a Nextflow-based workflow for identifying and annotating CNVs on human the **GRCh37** or **GRCh38** human reference genome.


## Pipeline Versions / Modes

1. **Full Pipeline**

   * Runs **PennCNV** and **QuantiSNP** from raw BAF/LRR files.
   * Recommended when starting with raw genotype data and wanting end-to-end CNV detection.

2. **Partial Pipeline**

   * Starts from **PennCNV** and **QuantiSNP** merged outputs.
   * Skips raw CNV calling, useful if CNV files are already generated.

> Both modes share downstream processing steps such as CNV merging and reporting.


## Prerequisites to Run the Pipeline

* **dataset\_name**: Name of the dataset, used for directory and report naming.
* **genome\_version**: Either **GRCh38** or **GRCh37**.


### Specific Inputs to **Full Pipeline** (pipeline_mode=pipeline_full)

* **plink2samplemetadata\_tsv**: A TSV file containing:

```
SampleID  Call_Rate  Sex  FatherID  MotherID
```

> *FatherID* and *MotherID* are only required if `--report "true"` is used, to assess Mendelian precision as an indicator of data quality.

* **list\_sample\_baflrrpath**: A TSV file containing:

```
SampleID  path_to_BAF_LRR
```

> The BAF/LRR file name should follow the format: **{SampleID}.BAF\_LRR.tsv**.

* **list\_baflrr\_path**: A TXT file containing paths to BAF/LRR files (essentially the second column of the previous file, e.g., `cut -f2`).
* **batch\_size**: Number of batches to process.

⚠️ Be careful with the `.hmm` file used by PennCNV in the **callBatchCNVs** process.


#### BAF\_LRR File Format

**{SampleID}.BAF\_LRR.tsv** files should have the following columns:

```
Name  Chr  Position  {SampleID}.Log R Ratio  {SampleID}.B Allele Freq
```

* `Chr` should be formatted as `"1"`–`"22"`, `"X"`, or `"Y"`.

**Recommendations for SNP positions:**

* It is highly recommended to have only unique SNP positions (i.e., no more than one SNP per genomic coordinate).

From personal communication with the authors of QuantiSNP (Ioannis Ragoussis and Rui Li):

> "At the end of the day, the tool integrates signals from different positions, so we need only one per genomic coordinate. Sometimes multiple probes are designed for the same SNP just in case one probe may fail. It is okay to keep the better-performing probe (best call rate). Hope that helps."

From the PennCNV author Wang Kai:

> "You can keep just one SNP to avoid potential issues. One easy way is simply to modify the PFB file to remove duplicated SNPs."


### Specific Inputs to **Partial Pipeline** (pipeline_mode=pipeline_partial)

* **penncnv\_calls\_path**: Path to PennCNV raw CNVs file (.txt)
* **quantisnp\_calls\_path**: Path to QuantiSNP raw CNVs file (.txt)

**Optional:**

* **plink2samplemetadata\_tsv**: A TSV file containing `SampleID  Call_Rate  Sex  FatherID  MotherID`. *FatherID* and *MotherID* are only required if `--report "true"` is used.
* **penncnv\_qc\_path**: QC metrics from PennCNV (SampleID → LRR\_SD, BAF\_SD, WF).


## Running the CNV Calling Pipeline

Users on Compute Canada (CCDB, in the lab) are encouraged to refer directly to :

### **Partial Pipeline**

```bash
sbatch /path/to/CNV-Caller/setup/ccdb/run_CNV_caller_light.sh \
    --git_dir /path/to/CNV-Caller \
    --dataset_name SPARK_Array_GRCh37 \
    --genome_version GRCh37 \
    --plink2samplemetadata_tsv /path/to/SPARK/sample_metadata_from_plink.tsv \
    --penncnv_calls_path /path/to/SPARK/PennCNV_raw_calls.txt \
    --quantisnp_calls_path /path/to/SPARK/QuantiSNP_raw_calls.txt \
    --penncnv_qc_path /path/to/SPARK/PennCNV_QC.tsv \
    [--report true|false]
```

### **Full Pipeline**

```bash
   sbatch /path/to/CNV-Caller/setup/ccdb/run_CNV_caller_full.sh \
    --git_dir /path/to/CNV-Caller \
          --dataset_name SPARK_Array_GRCh37 \
          --genome_version GRCh37 \
          --plink2samplemetadata_tsv file.tsv \
          --list_sample_baflrrpath list.tsv \
          [--batch_size 20] \
          [--report true|false]
```


## CNV Calling in the Full Pipeline

Within the pipeline, the GC model and PFB files required by PennCNV are prepared.
CNVs are then called per individual on autosomes and chromosome X using the container `docker://flobenhsj/quantisnp_penncnv:v2.2` and the following tools:

### PennCNV

PennCNV ([Wang *et al.*, 2007](http://www.genome.org/cgi/doi/10.1101/gr.6861907))

PennCNV is run separately on autosomes and on chromosome X using the commands below.

* Default parameters available in the Docker image:

  * `hmm_file="/usr/local/PennCNV-1.0.5/lib/wgs.hmm"`

```bash
# Autosomes
perl detect_cnv.pl --test \
    --conf \
    --pfbfile ${pfb_file} \
    --hmmfile ${hmm_file} \
    --logfile ${sample_id}.penncnv.log \
    --output ${sample_id}.penncnv.out \
    --gcmodelfile ${gcmodel_file} \
    ${BAF_LRR_Probes}

# Chromosome X
perl detect_cnv.pl --test \
    --conf \
    --pfbfile ${pfb_file} \
    --hmmfile ${hmm_file} \
    --logfile ${sample_id}.penncnv.chrx.log \
    --output ${sample_id}.penncnv.chrx.out \
    --gcmodelfile ${gcmodel_file} \
    --sexfile sexfile.tsv \
    --chrx \
    ${BAF_LRR_Probes}
```

### QuantiSNP

QuantiSNP ([Colella *et al.*, 2007](https://doi.org/10.1093/nar/gkm076))

QuantiSNP is executed using default configuration files (`levels.dat`, `params.dat`) and GC content files provided for each chromosome.
Sex is inferred beforehand, and the appropriate `--gender` flag is passed to ensure proper calling on sex chromosomes.

* Default parameters available in the Docker image:

  * `chr="1:23"`
  * `gcdir=/usr/local/QuantiSNP-2.3/GC_correction/${genome_version}/GCdir/`
  * `levels="/usr/local/QuantiSNP-2.3/bin/config/levels.dat"`
  * `config="/usr/local/QuantiSNP-2.3/bin/config/params.dat"`

```bash
quantisnp --chr ${chr} \
    --outdir . \
    --sampleid ${sample_id} \
    --gender ${gender} \
    --config ${config} \
    --levels ${levels} \
    --gcdir ${gcdir} \
    --input-files ${BAF_LRR_Probes} \
    --doXcorrect --verbose
```


## Paragraph method, by Florian Bénitière, 04/08/2025 : **CNV Merging process**

Copy Number Variants (CNVs) were called using two programs: PennCNV ([Wang et al., 2007](http://www.genome.org/cgi/doi/10.1101/gr.6861907)) and QuantiSNP ([Colella et al., 2007](https://doi.org/10.1093/nar/gkm076)).

A pre-filtering step was applied to both CNV call sets to exclude low-confidence variants. Specifically, CNVs with confidence scores below 30 in PennCNV or below 15 in QuantiSNP were removed (based on results from the SPARK cohort; see SPARK_PennCNV_calls.pdf and SPARK_QuantiSNP_calls.pdf). We also required CNVs to be longer than 1 kilobase, a commonly accepted threshold for CNVs ([Eichler, 2008](https://www.nature.com/scitable/topicpage/copy-number-variation-and-human-disease-698/)). In addition, CNVs overlapping the pseudoautosomal regions (PARs or XTR) of the X chromosome with a copy number of 2 were excluded, as these regions normally contain two copies in males (one on chrX and one on chrY) and therefore do not represent true CNVs.

CNVs from both tools were then merged to generate a unified CNV set. For each CNV, we calculated the fraction of support across the two algorithms (**Two\_Algorithm\_Overlap**). The CNV **Type** was inferred based on the copy number values within each unified CNV region: regions with all values \≥2 were labeled duplications ("DUP"), all values <2 as deletions ("DEL"), and regions with mixed values as mixed ("MIX").

We next assessed the overlap of each CNV with **Problematic Regions** defined by the UCSC Genome Browser ([UCSC Track: Problematic Regions GRCh37](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=problematic) and [GRCh38](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=problematic)).


## Filtering Recommendations
### CNVs
For downstream analyses, we recommend filtering the final CNV set using the following criteria (based on results from the SPARK cohort; see `resources/docs/SPARK_CNV_dataset_report.pdf`):

* **Two\_Algorithm\_Overlap \≥ 0.5** (validated by at least 50% overlap between both algorithms)
* **ProblematicRegions\_Overlap < 0.5** (less than 50% overlap with UCSC-flagged problematic regions)
* **Type == DEL or DUP** (only CNVs with an unambiguous type)

The **Confidence\_max** can also be used as the **Copy\_Number**, but based on our results we do not recommend adding this complexity to the method.

### Individuals
According to Liu, J., Zhang, L., Xu, L. *et al.*
*Analysis of copy number variations in the sheep genome using 50K SNP BeadChip array.*
**BMC Genomics** 14, 229 (2013). [https://doi.org/10.1186/1471-2164-14-229](https://doi.org/10.1186/1471-2164-14-229)

The quality filters used were:

* **LRR standard deviation (SD) < 0.30**
* **BAF drift < 0.01**
* **Waviness factor between −0.05 and 0.05**

In addition, based on common practice in CNV detection studies (e.g., PennCNV and Illumina genotyping guidelines), we also recommend applying:

* **Call Rate \≥ 0.98**


## Content of the dataset **CNV_merged_dataset.tsv**

| **Column**      | **Description**                                                 |
| --------------- | --------------------------------------------------------------- |
| **SampleID**                    | Unique identifier of the sample.            |
| **Chr**                         | Chromosome where the CNV is located.            |
| **Start**                       | Start position of the CNV.            |
| **End**                         | End position of the CNV.            |
| **Type**                        | Type of CNV: `DEL` = deletion, `DUP` = duplication, `MIX` = mixed type (event includes both deletions and duplications, inferred from Copy\_Number states). |
| **Length**                      | Length of the CNV in base pairs.            |
| **Copy\_Number**                | Distinct copy number states observed among merged CNVs (e.g., 0, 1, 3, 4). Multiple values are separated by commas if differing between callers. |
| **Confidence\_max**             | Maximum confidence score among merged CNVs (reflecting the strongest supporting evidence).            |
| **Num\_Probes\_max**            | Maximum number of probes supporting the CNV across merged CNVs.            |
| **Num\_Merged\_CNVs**           | Number of CNVs merged into this event.            |
| **QuantiSNP\_Overlap**          | Fraction of the CNV region overlapping QuantiSNP calls.            |
| **PennCNV\_Overlap**            | Fraction of the CNV region overlapping PennCNV calls.            |
| **Two\_Algorithm\_Overlap**     | Fraction of the CNV region supported by both QuantiSNP and PennCNV.            |
| **ProblematicRegions\_Overlap** | Overlap with genomic regions flagged as problematic (e.g., low complexity, centromeres, telomeres, segmental duplications, all collected from UCSC Genome Browser resources).|


## Content of the dataset **sampleDB.tsv** if produced

| **Column**      | **Description**                                                 |
| --------------- | --------------------------------------------------------------- |
| **SampleID**   | Unique sample identifier from the PLINK dataset.                 |
| **LRR\_mean**   | Mean Log R Ratio (from PennCNV output).                         |
| **LRR\_median** | Median Log R Ratio (from PennCNV output).                       |
| **LRR\_SD**     | Standard deviation of Log R Ratio (from PennCNV output).        |
| **BAF\_mean**   | Mean B Allele Frequency (from PennCNV output).                  |
| **BAF\_median** | Median B Allele Frequency (from PennCNV output).                |
| **BAF\_SD**     | Standard deviation of B Allele Frequency (from PennCNV output). |
| **BAF\_DRIFT**  | B Allele Frequency drift metric (from PennCNV output).          |
| **WF**          | Waviness Factor (from PennCNV output).                          |
| **GCWF**        | GC-corrected Waviness Factor (from PennCNV output).             |
|*__COLUMNS OF `plink2samplemetadata_tsv`__*|                            |	


## Other output in results

- **{dataset_name}**
  - **calls_unfiltered**
    - **penncnv**
      - `PennCNV_CNV.tsv` : Formatted PennCNV called CNV calls (filtered for length >1bp) in TSV.
      - `PennCNV_QC.tsv` : Quality control summary for PennCNV calls.
      - `PennCNV_raw_calls.txt` : Original CNV call output from PennCNV.
    - **quantisnp**
      - `QuantiSNP_CNV.tsv` : Formatted QuantiSNP called CNV calls (filtered for length >1bp) in TSV.
      - `QuantiSNP_raw_calls.txt` : Original CNV call output from QuantiSNP.
  - `CNV_merged_dataset.tsv` : **Final dataset with merged CNVs** following the protocol above.
  - **docs**
    - `launch_report.txt` : Summary of the analysis run.
    - `merged_cnv_qc.pdf` : QC report for the merged CNV dataset.
    - `penncnv_unfilter_cnv_qc.pdf` : QC report for unfiltered PennCNV calls.
    - `quantisnp_unfilter_cnv_qc.pdf` : QC report for unfiltered QuantiSNP calls.
    - `sample_qc_report.pdf` : Sample-level QC metric report.
  - `sampleDB.tsv` : **Sample database** containing metadata and QC status for all genotyped individuals.


## Current Limitations of the pipeline

- For PennCNV, default files such as **hhall.hmm** or **wgs.hmm** are used (available in the Docker docker://flobenhsj/quantisnp_penncnv:v2.2):
[Git PennCNV](https://github.com/WGLab/PennCNV)

- For QuantiSNP, default files such as **levels.dat** and **params.dat** are used (available in the Docker docker://flobenhsj/quantisnp_penncnv:v2.2):
[Git QuantiSNP](https://github.com/cwcyau/quantisnp)

- Resource requirements of each step must be adjusted depending on the quantity of data analyzed.

- Currently, only a Nextflow workflow has been written.
