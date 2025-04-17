# BED files with columns: chrom, start, end, GC_fraction
- gc_content_windows_GRCh37.bed 
- gc_content_windows_GRCh38.bed


## Generated from
Reference genomes downloaded on 16/04/2025 on 
https://useast.ensembl.org/Homo_sapiens/Info/Index for GRCh38.p14
https://grch37.ensembl.org/Homo_sapiens/Info/Index for GRCh37.p13

## Then fasta files have been bgzip and indexed:
gunzip -c Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz | bgzip > Homo_sapiens.GRCh37.dna.primary_assembly.fa.bgz
module load samtools
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa.bgz

## Finally, the GC content was computed across fixed-width genomic windows (999 bp) using the command:
./compute_gc_content_windows.sh <genome_fasta.bgz> gc_content_windows_GRChXXXX.bed 999

./compute_gc_content_windows.sh Homo_sapiens.GRCh37.dna.primary_assembly.fa.bgz gc_content_windows_GRCh37.bed 999
./compute_gc_content_windows.sh Homo_sapiens.GRCh38.dna.primary_assembly.fa.bgz gc_content_windows_GRCh38.bed 999