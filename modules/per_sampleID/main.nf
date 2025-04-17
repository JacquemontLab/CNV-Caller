









// Module process
process CNV_calling {
    tag "CNV_calling"

    input:
    File sample_to_path_vcf_gz
    File pfbfile_gz
    File gcmodelfile_gz
    File sexfile_gz
    String sample_id
    String workspace_id
    String? tsv_input_file = "IID_" + sample_id + "_vcf.tsv"
    String? hmmfile = "/usr/local/bin/hhall.hmm"
    String? levels = "/usr/local/quantisnp/bin/config/levels.dat"
    String? config = "/usr/local/quantisnp/bin/config/params.dat"
    String? gcdir = "/usr/local/quantisnp/gc/b38/"
    String? chr = "1:23"
    Int? minlength = 1000
    Float? minconf = 15.0


    output:
    path "{sample_id}.PennCNV_QC.tsv"
    path "{sample_id}.CNVs.tsv"

    script:
    """
    echo "Process Running: CNV_calling"

    # Print the number of cores
    nproc
    set -euo pipefail

    # Extract the sex of the sample
    zgrep -P "~{sample_id}\t" ~{sexfile_gz} | cut -f1,2 > sexfile.tsv
    gender=$(cut -f2 sexfile.tsv)

    echo "Gender: $gender."
    if [[ $gender != "male" && $gender != "female" ]]; then
        echo "ERROR: Invalid gender value for sample ID ~{sample_id}. Expected male or female."
        exit 1
    fi


    echo "Run PennCNV Autosome"
    /usr/bin/time -v bash -c '
        detectCNV="/usr/local/bin/detect_cnv.pl"
        perl $detectCNV --test \
            --pfbfile pfbfile.pfb \
            --hmmfile ~{hmmfile} \
            --logfile ~{sample_id}.penncnv.log \
            --minlength ~{minlength} \
            --minconf ~{minconf} \
            --output ~{sample_id}.penncnv.out \
            --gcmodelfile gcmodelfile.txt \
            ~{tsv_input_file}lite
    '
            
            
    echo "Run PennCNV ChrX"
    /usr/bin/time -v bash -c '
        detectCNV="/usr/local/bin/detect_cnv.pl"
        perl $detectCNV --test \
            --pfbfile pfbfile.pfb \
            --hmmfile ~{hmmfile} \
            --logfile ~{sample_id}.penncnv.chrx.log \
            --minlength ~{minlength} \
            --minconf ~{minconf} \
            --output ~{sample_id}.penncnv.chrx.out \
            --gcmodelfile gcmodelfile.txt \
            --sexfile sexfile.tsv \
            --chrx \
            ~{tsv_input_file}lite
    '
    
    echo "Run QuantiSNP"
    /usr/bin/time -v bash -c '                     
            quantisnp --chr ~{chr} \
                --outdir . \
                --sampleid ~{sample_id} \
                --gender $0 \
                --config ~{config} \
                --levels ~{levels} \
                --gcdir ~{gcdir} \
                --input-files ~{tsv_input_file}lite \
                --doXcorrect --verbose &> ~{sample_id}.quantisnp.log
            
            mv ~{sample_id}.cnv ~{sample_id}.quantisnp.cnv
    ' $gender
 
    cat ~{sample_id}.penncnv.chrx.out ~{sample_id}.penncnv.out > ~{sample_id}.penncnv.cnv

    merge_cnv_quantisnp_penncnv.sh ~{sample_id}.quantisnp.cnv ~{sample_id}.penncnv.cnv <probe_file> {sample_id}.CNVs.tsv

    """
}