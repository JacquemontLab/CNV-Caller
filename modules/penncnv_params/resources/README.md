# BED files with columns: chrom, start, end, GC_fraction
- gc_content_windows_GRCh37.bed 

- gc_content_windows_GRCh38.bed


## Generated from
Reference genomes downloaded on 16/04/2025 on 

https://useast.ensembl.org/Homo_sapiens/Info/Index for GRCh38.p14

wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz


https://grch37.ensembl.org/Homo_sapiens/Info/Index for GRCh37.p13

wget https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz


## Then fasta files have been bgzip and indexed:
module load tabix

gunzip -c Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz | bgzip > Homo_sapiens.GRCh37.dna.primary_assembly.fa.bgz

gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | bgzip > Homo_sapiens.GRCh38.dna.primary_assembly.fa.bgz


module load samtools

samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa.bgz

samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa.bgz


## Finally, the GC content was computed across fixed-width genomic windows (999 bp) using the command
## output (Chr\tStart\tEnd\tGC):
./compute_gc_content_windows.sh <genome_fasta.bgz> gc_content_windows_GRChXXXX.bed 999

./compute_gc_content_windows.sh Homo_sapiens.GRCh37.dna.primary_assembly.fa.bgz gc_content_1k_windows_GRCh37.bed 999

./compute_gc_content_windows.sh Homo_sapiens.GRCh38.dna.primary_assembly.fa.bgz gc_content_1k_windows_GRCh38.bed 999


## Create GCdir for QuantiSNP input
## outputs (Start\tEnd\tGC):
mkdir -p GRCh37_GCdir

awk -F'\t' '{
    OFS = "\t";
    if ($4 == "NA") {
        $4 = "";
    }
    {
        chr = $1;
        print $2, $3, $4 > "GRCh37_GCdir/" chr "_1k.txt";
    }
}' gc_content_1k_windows_GRCh37.bed

rm GRCh37_GCdir/GL*


mkdir -p GRCh38_GCdir

awk -F'\t' '{
    OFS = "\t";
    if ($4 == "NA") {
        $4 = "";
    }
    {
        chr = $1;
        print $2, $3, $4 > "GRCh38_GCdir/" chr "_1k.txt";
    }
}' gc_content_1k_windows_GRCh38.bed

rm GRCh38_GCdir/KI* GRCh38_GCdir/GL*