[![Jacquemont's Lab Header](labheader.png)](https://www.jacquemont-lab.org/)

[Git Repository CNV-Calling-Microarray](https://github.com/JacquemontLab/CNV-Calling-Microarray)

# Documentation of the CNV Calling and Annotation Pipeline for Microarray
This repository contains a bioinformatics pipeline for the discovery and annotation of copy number variants (CNVs).
The workflow is implemented using Nextflow to ensure reproducibility and efficient execution.

## Pipeline Overview
This repository provides a Nextflow-based workflow for identifying and annotating CNVs on human the **GRCh37** or **GRCh38** human reference genome.

### Workflow DAG
Below is a graphical representation of the workflow:

![Workflow DAG](dag.png)


### Prerequisite to run the Snakefile pipeline:
PLINK files: .bim, .fam, and .bed
Signal intensity files for each sample, named {SampleID}.BAF_LRR_Probes.tsv, with the following format:
Name\tChr\tPosition\t{SampleID}.Log R Ratio\t{SampleID}.B Allele Freq


### 1. Step TODO

TODO


## More Documentation
Here are some useful links about the plugins used in this pipeline, provided by VEP (Variant Effect Predictor).

### VEP

[Ensembl VEP Options](https://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html)\
[Ensembl VEP Plugins](https://useast.ensembl.org/info/docs/tools/vep/script/vep_plugins.html)\
[Ensembl VEP Consequences](https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html)

We used the VEP docker: ensembl-vep : 113.0


## Current Limitations of the pipeline

- For PennCNV, default files such as **hhall.hmm** are used (available in the Git repository):
[Git PennCNV](https://github.com/WGLab/PennCNV)

- For QuantiSNP, default files such as **levels.dat** and **params.dat** are used (available in the Git repository):
[Git QuantiSNP](https://github.com/cwcyau/quantisnp)

- Resource requirements of each step must be adjusted depending on the quantity of data analyzed.

- Currently, only the Nextflow workflow has been written.
