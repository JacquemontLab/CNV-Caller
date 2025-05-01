[Git Repository SNV-Annotation-SPARK
](https://github.com/JacquemontLab/CNV-Annotation-pipeline)

# Documentation of the CNV Calling and Annotation Pipeline

This repository contains a bioinformatics pipeline for discovering and annotating (CNVs).
The workflow is implemented with Nextflow to ensure reproducibility and efficient execution.

## Different Steps of the Pipeline or structured by modules
This repository provides a Nextflow-based workflow for identifying and annotating CNVs on human references genome version **GRCh37** or **GRCh38**.

### Workflow DAG
Below is a graphical representation of the workflow:

![Workflow DAG](dag.png)


### Prerequisite to run the Snakefile pipeline:
Preparing the PLINK files (.bim, .fam, .bed)
Having the Signal intensity files per SampleID stored in a directory such as named {SampleID}.BAF_LRR_Probes.tsv
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

- For PennCNV default files are used such as hmm profile **hhall.hmm** available on the git.
[Git PennCNV](https://github.com/WGLab/PennCNV)

- For QuantiSNP default files are used such as **levels.dat** and **params.dat** available on the git.
[Git QuantiSNP](https://github.com/cwcyau/quantisnp)

- Resource requirements of each step must be adjusted depending on the quantity of data analyzed.

- Currently, only the Nextflow workflow has been written.
