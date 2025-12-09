### Paragraph method, by Florian Bénitière, 04/08/2025

## CNV Calling

Copy Number Variants (CNVs) were called using two programs: PennCNV ([Wang et al., 2007](http://www.genome.org/cgi/doi/10.1101/gr.6861907)) and QuantiSNP ([Colella et al., 2007](https://doi.org/10.1093/nar/gkm076)).

A pre-filtering step was applied to both CNV call sets to exclude low-confidence variants. Specifically, CNVs with confidence scores below 30 in PennCNV or below 15 in QuantiSNP were removed (based on results from the SPARK cohort; see SPARK_PennCNV_calls.pdf and SPARK_QuantiSNP_calls.pdf). We also required CNVs to be longer than 1 kilobase, a commonly accepted threshold for CNVs ([Eichler, 2008](https://www.nature.com/scitable/topicpage/copy-number-variation-and-human-disease-698/)). In addition, CNVs overlapping the pseudoautosomal regions (PARs or XTR) of the X chromosome with a copy number of 2 were excluded, as these regions normally contain two copies in males (one on chrX and one on chrY) and therefore do not represent true CNVs.

CNVs from both tools were then merged to generate a unified CNV set. For each CNV, we calculated the fraction of support across the two algorithms. The CNV type was inferred based on the copy number values within each unified CNV region: regions with all values ≥2 were labeled duplications ("DUP"), all values <2 as deletions ("DEL"), and regions with mixed values as mixed ("MIX").

We next assessed the overlap of each CNV with **Problematic Regions** defined by the UCSC Genome Browser ([UCSC Track: Problematic Regions](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=problematic)).

For downstream analyses, we recommend filtering the final CNV set using the following criteria (based on results from the SPARK cohort; see SPARK_CNV_dataset_report.pdf):

* **Two\_Algorithm\_Overlap ≥ 0.7** (validated by at least 70% overlap between both algorithms)
* **ProblematicRegions\_Overlap < 0.5** (less than 50% overlap with UCSC-flagged problematic regions)
* **Type == DEL or DUP** (Only CNVs with an unambiguous type)






## Content of the dataset

| **Column**                      | **Description**            |
| ------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
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

