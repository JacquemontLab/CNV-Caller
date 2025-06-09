

gsutil -m -u $workspace_id cp -r gs://fc-aou-datasets-controlled/v8/microarray/plink/ .

gsutil -m -u $workspace_id ls gs://fc-aou-datasets-controlled/


gsutil -m -u $workspace_id cp -r ${bucket_id}output_cnvcalling/cnv_annotation .



workspace_id="terra-vpc-sc-bf15287e"
bucket_id="gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/"

gsutil -m -u $workspace_id cp -r gs://fc-aou-datasets-controlled/v8/microarray/plink/ .
gsutil -m -u $workspace_id cp ${bucket_id}cnvcalling_inputs/from_pfb_gcModel.tsv.gz resources/



gsutil -m -u $workspace_id cp gs://fc-aou-datasets-controlled/v8/microarray/plink/ .



zcat quantisnp_cnv.tsv.gz | head -n 2000000  > test_quantisnp.tsv
zcat penncnv_cnv.tsv.gz | head -n 2000000  > test_penncnv.tsv


input_file_penncnv=test_penncnv.tsv
input_file_quantisnp=test_quantisnp.tsv



probe_file="resources/from_pfb_gcModel.tsv"
regions_file="resources/Genome_Regions_data.tsv"
genome_version="GRCh38"
output="microarray_cnv.tsv"


time ./scripts/merge_cnv_quantisnp_penncnv.sh ${input_file_quantisnp} ${input_file_penncnv} ${probe_file} ${regions_file} ${genome_version} ${output}



gsutil -m -u $workspace_id cp ${bucket_id}cnvcalling_inputs/manifest_tier_v8.tsv.gz .


./scripts/format_penncnv_cnv.sh inputs/merged.pcnv.cnv.gz penncnv_cnv.tsv
./scripts/format_quantisnp_cnv.sh inputs/merged.quanti.cnv.gz quantisnp_cnv.tsv


export TMPDIR=/home/jupyter/cnv_annotation/tmp

probe_file="resources/from_pfb_gcModel.tsv"
regions_file="resources/Genome_Regions_data.tsv"
genome_version="GRCh38"
output="microarray_cnv.tsv"

output_dir=out/

# Make sure output directory exists
mkdir -p "$output_dir"

num_cores=96
num_group=1000

# sample_list=list_sample.txt
sample_list=resources/list_sample.txt
# Split sample list
split -n l/$num_group "$sample_list" "$output_dir/sample_group_"


input_file_quantisnp=quantisnp_cnv.tsv
input_file_penncnv=penncnv_cnv.tsv

current_jobs=0
# Loop through each group file
for group_file in $(find "$output_dir" -name "sample_group_*" | sort); do
  group_id=$(basename "$group_file")

    echo $current_jobs
    (
    head -n 1 "$input_file_quantisnp" > "$output_dir/quantisnp_${group_id}.tsv"
    awk 'NR==FNR { samples[$1]; next } FNR > 1 && $1 in samples' "$group_file" "$input_file_quantisnp" >> "$output_dir/quantisnp_${group_id}.tsv"

    head -n 1 "$input_file_quantisnp" > "$output_dir/penncnv_${group_id}.tsv"
    awk 'NR==FNR { samples[$1]; next } FNR > 1 && $1 in samples' "$group_file" "$input_file_penncnv" >> "$output_dir/penncnv_${group_id}.tsv"

    # Run the script
    # time ./scripts/merge_cnv_quantisnp_penncnv.sh \
    #     "$output_dir/quantisnp_$group_id.tsv" \
    #     "$output_dir/penncnv_$group_id.tsv" \
    #     "$probe_file" "$regions_file" "$genome_version" \
    #     "$output_dir/output_$group_id.tsv"
    ) &

  current_jobs=$((current_jobs + 1))

  # If we've hit the limit, wait for one job to finish
  if [ "$current_jobs" -ge "$num_cores" ]; then
    wait -n
    current_jobs=$((current_jobs - 1))
  fi
done

# Wait for all remaining jobs
wait


output_dir=out/
merged_file="microarray_cnv.tsv"

# Get header from the first output file
head -n 1 "$(find "$output_dir" -name "output_*.tsv" | sort | head -n 1)" > "$merged_file"

# Append data (skip header) from all files
for file in $(find "$output_dir" -name "output_*.tsv" | sort); do
  tail -n +2 "$file" >> "$merged_file"
done








awk 'NR==FNR { samples[$1]; next } FNR > 1 && $1 in samples' "$group_file" "$input_file_quantisnp"















  awk "NR==FNR{samples[\$1]; next} \$1 in samples" "$group_file" '"$input_file_penncnv"' > "$output_dir/penncnv_$group_id.txt"


















zcat penncnv_cnv.tsv.gz | wc -l 









SCRIPT_DIR=./scripts/

qs_clean=qs_clean.tsv
pc_clean=pc_clean.tsv
probes_correct=probes_correct.tsv
head probes_correct.tsv

cp "$qs_clean"              ./qs_clean.tsv
cp "$pc_clean"              ./pc_clean.tsv
cp "$combined_bed"          ./combined_cnv.bed
cp "$merged_bed"            ./merged_raw.bed
cp "$merged_tsv"            ./merged.tsv
cp "$probes_bed"            ./probes.bed
cp "$probes_correct"        ./probes_correct.tsv
cp "$overlap_raw"           ./overlap_raw.tsv
cp "$overlap_raw_sorted"    ./overlap_raw_sorted.bed
cp "$qs_bed"                ./qs.bed
cp "$pc_bed"                ./pc.bed
cp "$overlap_bed"           ./overlap_len.bed
cp "$overlap_sum_bed"       ./overlap_sum_bed.bed
cp "$frac_bed"              ./frac_overlap.bed


./scripts/merge_cnv_quantisnp_penncnv.sh ${input_file_quantisnp} ${input_file_penncnv} ${probe_file} ${regions_file} ${genome_version} ${output}




gsutil -m -u $workspace_id cp ${bucket_id}cnvcalling_inputs/manifest_tier_v8.tsv.gz .


file=gs://fc-aou-datasets-controlled/pooled/microarray/vcf/v8_delta/1000000.21317004391.sorted.vcf.gz

gsutil -m -u $workspace_id cp $file .



gsutil -m -u $workspace_id cat gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cnvcalling_inputs/random_1000_samples.pfb.gz | zcat | head| wc -l



gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/cnv_annotation/inputs/


gsutil -m -u $workspace_id cp ${bucket_id}output_cnvcalling/cnv_annotation/inputs/merge.pcnv.cnv.gz .


zcat merge.pcnv.cnv.gz | head | awk '{gsub(/[[:space:]]+/, "\t"); print}' | cut -f1 | cut -d ':' -f1 | sort | uniq -c

zcat merge.pcnv.cnv.gz | awk '{gsub(/[[:space:]]+/, "\t"); print}' | cut -f5 | sort | uniq -c | wc -l




gsutil -m -u $workspace_id cat ${bucket_id}output_cnvcalling/penncnv/1000039.pcnv.out.gz | zcat | head

gsutil -m -u $workspace_id cat ${bucket_id}output_cnvcalling/quantisnp/1000000.quanti.cnv.gz | zcat | head




dir=fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/4bb9824a-dd6d-4286-ad88-8618620622d5/call-CNVcalling/shard-0/


gsutil -m -u $workspace_id cat gs://${dir}CNVcalling-0.log | wc -l
gs://fc-aou-datasets-controlled/pooled/microarray/vcf/v7_base/1000175.21037003001.sorted.vcf.gz,gs://fc-aou-datasets-controlled/pooled/microarray/vcf/v7_base/1000175.21037003001.sorted.vcf.gz.tbi

gsutil -m -u $workspace_id cp ${bucket_id}README.md .

job_id=e072a654-7e83-4bb3-bb5c-83efb9c76334

file=${bucket_id}cromwell-execution/RunPipeline/${job_id}/call-CNVcalling/shard-42/3365159.


# Transfer results to output_cnvcalling directory
workspace_id="terra-vpc-sc-a6da20a8"
bucket_id="gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/"
job_id=8d04138b-7912-492f-8918-9e9823f1451b

gsutil -m -u $workspace_id cp ${bucket_id}cromwell-execution/RunPipeline/${job_id}/*/*/*/*.gz ${bucket_id}output_cnvcalling/penncnv_quantisnp/

gsutil -m -u $workspace_id cp ${bucket_id}cromwell-execution/RunPipeline/*/call-CNVcalling/shard-*/*.gz ${bucket_id}output_cnvcalling/penncnv_quantisnp/
gsutil -m -u $workspace_id cp ${bucket_id}cromwell-execution/RunPipeline/${job_id}/call-CNVcalling/shard-*/glob*/*.gz ${bucket_id}output_cnvcalling/penncnv_quantisnp/


workspace_id="terra-vpc-sc-a6da20a8"
bucket_id="gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/"

gsutil -m -u $workspace_id cp ${bucket_id}cromwell-execution/RunPipeline/**/*.gz ${bucket_id}output_cnvcalling/penncnv_quantisnp/

gsutil -m -u $workspace_id cp ${bucket_id}cromwell-execution/RunPipeline/**.gz ${bucket_id}output_cnvcalling/penncnv_quantisnp/

gsutil -m -u $workspace_id ls ${bucket_id}cromwell-execution/RunPipeline/*/call-CNVcalling/shard-*/1098187.*
gsutil -m -u $workspace_id ls ${bucket_id}cromwell-execution/RunPipeline/*/call-CNVcalling/shard-*/1110799.*

1098187 1110799

gsutil -m -u $workspace_id cp ${file}* .

# Move .gz files to the new directory within the bucket
gsutil -m -u $workspace_id ls ${bucket_id}cromwell-execution/RunPipeline/${job_id}/call-CNVcalling/*/*.pcnv.chrx.log.gz | wc -l

gsutil -m -u $workspace_id ls ${bucket_id}cromwell-execution/RunPipeline/*/*/*/*.gz | head

gsutil -m -u $workspace_id cp ${bucket_id}cromwell-execution/RunPipeline/*/*/*/*.gz ${bucket_id}output_cnvcalling/penncnv_quantisnp/

gsutil -m -u $workspace_id cp ${bucket_id}cromwell-execution/RunPipeline/*/*/*/CNVcalling-*.log ${bucket_id}log_cnvcalling/penncnv_quantisnp/




gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv_quantisnp/ | head

gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv_quantisnp/ | wc -l

gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv_quantisnp/*.pcnv.chrx.out.gz | wc -l
gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv_quantisnp/*.pcnv.chrx.log.gz | wc -l
gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv_quantisnp/*.pcnv.log.gz | wc -l
gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv_quantisnp/*.pcnv.out.gz | wc -l
gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv_quantisnp/*.quanti.cnv.gz | wc -l
gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv_quantisnp/*.quanti.log.gz | wc -l
gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv_quantisnp/*.quanti.loh.gz | wc -l
gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv_quantisnp/*.quanti.qc.gz | wc -l


gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv_quantisnp/*.quanti.qc.gz | sed -E 's|.*/([^.]+)\.quanti.*|\1|' > sample_processed.txt

gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv_quantisnp/ | wc -l


comm -23 <(sort random_1000.txt) <(sort sample_processed.txt | sort) > new_sample_to_runv1.txt

zcat manifest.tsv.gz | tail -n +2 | awk '{print $1}' | shuf | sort | comm -23 - <(sort sample_processed.txt) | head -n 25000 > sample_ids_file_n25000.txt

gsutil -m -u $workspace_id cp sample_ids_file_n25000.txt ${bucket_id}cnvcalling_inputs/to_run/

zcat manifest.tsv.gz | awk '{print $1}' | head | shuf | sort | comm -12 <(sort sample_processed.txt) - > new_sample_to_run.txt

mkdir outputfiles


gsutil -m -u $workspace_id cp ${bucket_id}cnvcalling_inputs/sex_at_birth.tsv.gz .
zcat sex_at_birth.tsv.gz | head

gsutil -m -u $workspace_id cp ${bucket_id}cnvcalling_inputs/sex_file.tsv.gz .
zcat sex_file.tsv.gz | head

### Collect X samples output on workbench
list_sample=$(head -n 20000 sample_processed.txt)
for sample in $list_sample; do 
    echo $sample
    gsutil -q -m -u $workspace_id cp ${bucket_id}output_cnvcalling/penncnv_quantisnp/${sample}.quanti.cnv.gz /home/rstudio/outputfiles/
    gsutil -q -m -u $workspace_id cp ${bucket_id}output_cnvcalling/penncnv_quantisnp/${sample}.pcnv.out.gz /home/rstudio/outputfiles/
done


gsutil -q -m -u $workspace_id cp ${bucket_id}output_cnvcalling/penncnv_quantisnp/*.quanti.cnv.gz /home/rstudio/outputfiles/
gsutil -m -u $workspace_id cp ${bucket_id}output_cnvcalling/penncnv_quantisnp/*.pcnv.out.gz /home/rstudio/outputfiles/

ls /home/rstudio/outputfiles/*.quanti.cnv.gz | wc -l



job_id=8d04138b-7912-492f-8918-9e9823f1451b

number=3193
gsutil -m -u $workspace_id cat ${bucket_id}cromwell-execution/RunPipeline/${job_id}/call-CNVcalling/shard-${number}/CNVcalling-${number}.log | grep "Run Penn"


gsutil -m -u $workspace_id cat ${bucket_id}cromwell-execution/RunPipeline/${job_id}/call-CNVcalling/shard-${number}/stderr | grep "Run Penn"


number=8432
cat ${directory_name}/CNVcalling-${number}.log | grep "Run Penn"

1164543
gsutil -m -u $workspace_id  cp "gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cnvcalling_inputs/manifest.tsv.gz" .

gsutil -m -u $workspace_id cp gs://fc-aou-datasets-controlled/v7/microarray/vcf/manifest_deprecated.csv .
gsutil -m -u $workspace_id ls gs://fc-aou-datasets-controlled/v7/microarray/vcf/



comm -23 <(sort sample_processed.txt) <(cut -f1 manifest.tsv | sort)

touch output_dir_to_suppr
cat sample_deprecated.txt | while read sample 
do
    echo $sample
    gsutil -m -u $workspace_id ls ${bucket_id}cromwell-execution/RunPipeline/*/*/*/${sample}.pcnv.chrx.log.gz >> output_dir_to_suppr
done



workspace_id="terra-vpc-sc-a6da20a8"
bucket_id="gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/"


gsutil -m -u $workspace_id ls gs://fc-aou-datasets-controlled/v7/microarray/plink_v7.1/arrays.*

gsutil -m -u $workspace_id cp gs://fc-aou-datasets-controlled/v7/microarray/plink_v7.1/* plink/



### Check sex of the generated files of king
plink --bfile plink/arrays --check-sex --out check-sex_plink

gsutil -q -m -u $workspace_id cp check-sex_plink.sexcheck.tsv.gz ${bucket_id}




head sample_deprecated.txt

grep 1000608 logfiles/*/*



cat logfiles/samples5000_cpu1_mem3gb/CNVcalling-100.log

grep  "Failed to open file"  logfiles/*/* > output_suppr

sed 's/\.log.*//' output_suppr > output.txt

sed -z -i 's/logfiles\/samples50_cpu1_mem3gb\/CNVcalling-/4bb9824a-dd6d-4286-ad88-8618620622d5;/g' output.txt
sed -z -i 's/logfiles\/samples1000_cpu1_mem3gb\/CNVcalling-/6f0bff24-37a4-4fe0-8260-69afb5968ca5;/g' output.txt
sed -z -i 's/logfiles\/samples5000_cpu1_mem3gb\/CNVcalling-/7d6ff0f6-fee2-440a-81f1-1565fdc475c1;/g' output.txt
sed -z -i 's/logfiles\/samples8950_cpu1_mem3gb\/CNVcalling-/8d04138b-7912-492f-8918-9e9823f1451b;/g' output.txt
sed -z -i 's/logfiles\/samples10000_cpu1_mem3gb\/CNVcalling-/a8b52af5-5fb9-447a-a694-387073f27927;/g' output.txt


awk -F';' '{num=$2; print "gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/"$1"/call-CNVcalling/shard-"num"/CNVcalling-"num".log"}' output.txt > path_output.txt

gsutil -m -u $workspace_id cat gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/a8b52af5-5fb9-447a-a694-387073f27927/call-CNVcalling/shard-2474/CNVcalling-2474.log


cat path_output.txt | while read sample 
do
    echo $sample
    echo $(dirname $sample)
    # gsutil -m -u $workspace_id cat $sample | grep "Failed to open file"
    gsutil -m -u $workspace_id rm -r $(dirname $sample)
done




comm -23 <(sort random_1000.txt) <(sort sample_processed.txt | sort) > new_sample_to_runv1.txt


zgrep -E "male|female" sex_file.tsv.gz | awk '{print $1}' | shuf | sort | comm -12 <(sort new_sample_to_runv1.txt) - > new_sample_to_run.txt

gsutil -m -u $workspace_id mv ${bucket_id}cnvcalling_inputs/new_sample_to_run.txt ${bucket_id}cnvcalling_inputs/to_run/new_sample_to_run.txt

gsutil -m -u $workspace_id cat ${bucket_id}cnvcalling_inputs/to_run/header_new_sample_to_run.txt | wc -l

gsutil -m -u $workspace_id cp header_new_sample_to_run.txt ${bucket_id}cnvcalling_inputs/to_run/header_new_sample_to_run.txt



zgrep -E "male|female" sex_file.tsv.gz | awk '{print $1}' | shuf | sort | comm -23 <(sort random_1000.txt) - | head -n 8950 > sample_ids_file_n8950.txt


gsutil -m -u $workspace_id ls gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/cb1f0d9a-9575-4888-82fb-b453bc5087ba/call-CNVcalling/shard-1/

gsutil -m -u $workspace_id cat gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cnvcalling_inputs/to_run/sample_ids_file_n25000.txt


gsutil -m -u $workspace_id ls ${bucket_id}cnvcalling_inputs/to_run/header_new_sample_to_run.txt

iid=2561084

gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv_quantisnp/${iid}.pcnv.chrx.out.gz .
gsutil -m -u $workspace_id cp ${bucket_id}output_cnvcalling/penncnv_quantisnp/${iid}.pcnv.log.gz .
zcat ${iid}.pcnv.log.gz | wc -l

gsutil -m -u $workspace_id cp ${bucket_id}output_cnvcalling/penncnv_quantisnp/${iid}.pcnv.out.gz .
zcat ${iid}.pcnv.out.gz | wc -l

gsutil -m -u $workspace_id ls -la ${bucket_id}cromwell-execution/RunPipeline/4bb9824a-dd6d-4286-ad88-8618620622d5/call-CNVcalling/shard-0/${iid}.pcnv.out.gz
gsutil -m -u $workspace_id cp ${bucket_id}cromwell-execution/RunPipeline/4bb9824a-dd6d-4286-ad88-8618620622d5/call-CNVcalling/shard-0/${iid}.pcnv.out.gz .
zcat ${iid}.pcnv.out.gz | head

gsutil -m -u $workspace_id ls -la gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/741b4d60-09bb-4d37-852b-0bfb41cb6afa/call-CNVcalling/shard-0/glob-9f0876873ace9c952bdcdbb9c2697f09/${iid}.pcnv.out.gz
gsutil -m -u $workspace_id cp gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/741b4d60-09bb-4d37-852b-0bfb41cb6afa/call-CNVcalling/shard-0/glob-9f0876873ace9c952bdcdbb9c2697f09/${iid}.pcnv.out.gz .
zcat ${iid}.pcnv.out.gz | head


gsutil -m -u $workspace_id ls -la ${bucket_id}cromwell-execution/RunPipeline/4bb9824a-dd6d-4286-ad88-8618620622d5/call-CNVcalling/shard-0/${iid}.quanti.cnv.gz
gsutil -m -u $workspace_id cp ${bucket_id}cromwell-execution/RunPipeline/4bb9824a-dd6d-4286-ad88-8618620622d5/call-CNVcalling/shard-0/${iid}.quanti.cnv.gz 1.txt.gz
zcat 1.txt.gz | head

gsutil -m -u $workspace_id ls -la ${bucket_id}cromwell-execution/RunPipeline/741b4d60-09bb-4d37-852b-0bfb41cb6afa/call-CNVcalling/shard-0/glob-d19760cba9d2ba161d38e3b6659f9ba7/${iid}.quanti.cnv.gz
gsutil -m -u $workspace_id cp ${bucket_id}cromwell-execution/RunPipeline/741b4d60-09bb-4d37-852b-0bfb41cb6afa/call-CNVcalling/shard-0/glob-d19760cba9d2ba161d38e3b6659f9ba7/${iid}.quanti.cnv.gz 2.txt.gz
zcat 2.txt.gz | head




gsutil -m -u $workspace_id ls ${bucket_id}cromwell-execution/RunPipeline/4bb9824a-dd6d-4286-ad88-8618620622d5/call-CNVcalling/* | wc -l


gsutil -m -u $workspace_id ls gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cnvcalling_inputs/to_run
gsutil -m -u $workspace_id rm gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cnvcalling_inputs/to_run/header_new_sample_to_run.txt








mkdir pcnv_log
workspace_id="terra-vpc-sc-a6da20a8"
bucket_id="gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/"

gsutil -m -u $workspace_id cp gs://fc-secure-07cc78a0-0516-4955-a996-899e381a51f9/meta/GSA.sampleManifestAncestry.txt sex_file_from_autocallGender.tsv

awk 'NR==1 {print; next} {gsub(/\tM/, "\tmale"); gsub(/\tF/, "\tfemale"); print}' keep_column.txt > filter

cut -f2 sex_file_from_autocallGender.tsv | sort | uniq -c | sort -nr


comm -3 <(cut -f1 filter | sort) <(cut -f1 keep_column.txt | sort)
comm -3 <(cut -f1 sex_file_from_autocallGender.tsv | sort) <(zcat sex_file.tsv.gz | cut -f1 | sort)


gsutil -m -u $workspace_id cp sex_file_from_autocallGender.tsv.gz ${bucket_id}cnvcalling_inputs/sex_file_from_autocallGender.tsv.gz

gsutil -m -u $workspace_id cp ${bucket_id}output_cnvcalling/penncnv_quantisnp/*.pcnv.log.gz pcnv_log/


echo -e "IID\tLRR_mean\tLRR_median\tLRR_SD\tBAF_mean\tBAF_median\tBAF_SD\tBAF_DRIFT\tGCWF\tWF" > "QC_per_IID.tsv"
zcat pcnv_log/*.pcnv.log.gz | grep "quality summary" | sed -r ' s/NOTICE: quality summary for IID_//g; s/_vcf.tsvlite://g; s/LRR_mean=//g; s/LRR_median=//g; s/LRR_SD=//g; s/BAF_mean=//g; s/BAF_median=//g; s/BAF_SD=//g; s/BAF_DRIFT=//g;  s/GCWF=//g; s/WF=//g' | sed -r 's/\s+/\t/g' >> "QC_per_IID.tsv"
 

gsutil -m -u $workspace_id cp ${bucket_id}QC_per_IID.tsv .

gsutil -m -u $workspace_id mv ${bucket_id}check-sex_plink.sexcheck.tsv.gz ${bucket_id}cnvcalling_inputs/check-sex_plink.sexcheck.tsv.gz



job_id=3dcaa898-e00a-44dd-b9f9-c86c1a55388c


for job_id in 3dcaa898-e00a-44dd-b9f9-c86c1a55388c 4bb9824a-dd6d-4286-ad88-8618620622d5 6f0bff24-37a4-4fe0-8260-69afb5968ca5 7d6ff0f6-fee2-440a-81f1-1565fdc475c1 8d04138b-7912-492f-8918-9e9823f1451b a8b52af5-5fb9-447a-a694-387073f27927 e072a654-7e83-4bb3-bb5c-83efb9c76334/; do
echo $job_id
gsutil -m -u $workspace_id ls ${bucket_id}cromwell-execution/RunPipeline/${job_id}/call-CNVcalling | wc -l
done







${job_id}/*/*/*/*.gz


gsutil -m -u $workspace_id cp gs://fc-aou-datasets-controlled/v7/wgs/short_read/structural_variants/v7_offcycle/aux/aneuploidies/AoU_srWGS_SV.v7_offcycle.germline_allosome_aneuploidies.tsv .



jupyter@2e6735f95b9a:~$ gsutil -m -u $workspace_id ls ${bucket_id}cromwell-execution/RunPipeline/4bb9824a-dd6d-4286-ad88-8618620622d5/call-CNVcalling | wc -l
50
jupyter@2e6735f95b9a:~$ gsutil -m -u $workspace_id ls ${bucket_id}cromwell-execution/RunPipeline/6f0bff24-37a4-4fe0-8260-69afb5968ca5/call-CNVcalling | wc -l
1000
jupyter@2e6735f95b9a:~$ gsutil -m -u $workspace_id ls ${bucket_id}cromwell-execution/RunPipeline/7d6ff0f6-fee2-440a-81f1-1565fdc475c1/call-CNVcalling | wc -l
5000
jupyter@2e6735f95b9a:~$ gsutil -m -u $workspace_id ls ${bucket_id}cromwell-execution/RunPipeline/8d04138b-7912-492f-8918-9e9823f1451b/call-CNVcalling | wc -l
8950
jupyter@2e6735f95b9a:~$ gsutil -m -u $workspace_id ls ${bucket_id}cromwell-execution/RunPipeline/a8b52af5-5fb9-447a-a694-387073f27927/call-CNVcalling | wc -l
10000


workspace_id="terra-vpc-sc-a6da20a8"
bucket_id="gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/"
job_id=a8b52af5-5fb9-447a-a694-387073f27927


number=9015
gsutil -m -u $workspace_id ls ${bucket_id}cromwell-execution/RunPipeline/${job_id}/call-CNVcalling/shard-${number}/

gsutil -m -u $workspace_id rm -r ${bucket_id}cromwell-execution/RunPipeline/${job_id}/call-CNVcalling/shard-${number}/

for i in 1098187 1110799; do
    echo $i
    gsutil -m -u $workspace_id rm ${bucket_id}output_cnvcalling/penncnv_quantisnp/$i.*
done

zgrep -P "1006552\t" sex_file.tsv.gz > sexfile.tsv
gender=$(cut -f2 sexfile.tsv)

zgrep -P "~{sample_id}\t" ~{sexfile_gz} | cut -f1,2 > sexfile.tsv


for sample_id in $(zcat manifest.tsv.gz | cut -f1) ;do
    path_to_vcf=$(zgrep -P "${sample_id}\t" manifest.tsv.gz | cut -f2)
    zgrep -P "${sample_id}\t" sex_file.tsv.gz | cut -f1,2,3
    gsutil -m -u $workspace_id cat ${path_to_vcf} | head -n 5000 | zcat | grep Gender
done

# Output file
output_file="sex_info.tsv"
echo -e "Sample_ID\tGender_Source\tGender_Info" > $output_file

# Process each sample ID
zcat manifest.tsv.gz | tail -n +2 | cut -f1 | while read sample_id; do
    # Extract VCF path
    path_to_vcf=$(zgrep -P "${sample_id}\t" manifest.tsv.gz | cut -f2)
    
    # Extract sex information
    sex_info=$(zgrep -P "${sample_id}\t" sex_file.tsv.gz | cut -f1,2)
    
    # Extract gender-related lines from the VCF
    vcf_gender_info=$(gsutil -m -u $workspace_id cat ${path_to_vcf} | head -n 5000 | zcat | grep Gender)
    
    # Format and write output
    echo -e "${sex_info}\tplink" >> $output_file
    echo -e "$vcf_gender_info" | while read line; do
        echo -e "${sample_id}\tVCF\t${line}" >> $output_file
    done
done



gsutil -m -u "$workspace_id" ls ${bucket_id}cromwell-execution/RunPipeline/ | xargs -I {} sh -c 'echo "{}call-CNVcalling: $(gsutil -m -u terra-vpc-sc-a6da20a8 ls "{}call-CNVcalling" | wc -l)"'


gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/4bb9824a-dd6d-4286-ad88-8618620622d5/call-CNVcalling: 50
gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/6f0bff24-37a4-4fe0-8260-69afb5968ca5/call-CNVcalling: 996
gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/7d6ff0f6-fee2-440a-81f1-1565fdc475c1/call-CNVcalling: 4964
gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/7d9e54ed-929c-40ac-ac16-e40e98e398c1/call-CNVcalling: 25000
gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/8d04138b-7912-492f-8918-9e9823f1451b/call-CNVcalling: 8910
gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/a8b52af5-5fb9-447a-a694-387073f27927/call-CNVcalling: 9955
gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/e072a654-7e83-4bb3-bb5c-83efb9c76334/call-CNVcalling: 901



gsutil -m -u "$workspace_id" ls ${bucket_id}cromwell-execution/RunPipeline/7d9e54ed-929c-40ac-ac16-e40e98e398c1/call-CNVcalling/shard-5995/
gsutil -m -u "$workspace_id" ls "gs://${bucket_id}cromwell-execution/RunPipeline/7d9e54ed-929c-40ac-ac16-e40e98e398c1/call-CNVcalling/shard-*/**/*.gz"


workspace_id="terra-vpc-sc-a6da20a8"
bucket_id="gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/"
for jobid in 4bb9824a-dd6d-4286-ad88-8618620622d5 6f0bff24-37a4-4fe0-8260-69afb5968ca5 7d6ff0f6-fee2-440a-81f1-1565fdc475c1 7d9e54ed-929c-40ac-ac16-e40e98e398c1 8d04138b-7912-492f-8918-9e9823f1451b a8b52af5-5fb9-447a-a694-387073f27927 e072a654-7e83-4bb3-bb5c-83efb9c76334; do
    gsutil -m -u "$workspace_id" ls ${bucket_id}cromwell-execution/RunPipeline/${jobid}/call-CNVcalling/shard-*/**/*.gz > ${jobid}_gz_files_list.txt
done


for file in .pcnv.chrx.out.gz .pcnv.chrx.log.gz .pcnv.out.gz .pcnv.log.gz; do
    grep ${file} 7d9e54ed-929c-40ac-ac16-e40e98e398c1_gz_files_list.txt | wc -l
done


for file in .quanti.cnv.gz .quanti.loh.gz .quanti.log.gz .quanti.qc.gz ; do
    grep ${file} 7d9e54ed-929c-40ac-ac16-e40e98e398c1_gz_files_list.txt | wc -l
done

for file in .pcnv.chrx.out.gz .pcnv.chrx.log.gz .pcnv.out.gz .pcnv.log.gz; do 
    echo "${file}: $(grep "${file}" 7d9e54ed-929c-40ac-ac16-e40e98e398c1_gz_files_list.txt | wc -l)"
done


grep "${file}" /home/jupyter/7d9e54ed-929c-40ac-ac16-e40e98e398c1_gz_files_list.txt | grep -oP '(?<=/)[^/]+(?=\.gz)' | sort | uniq


grep -oP '(?<=/)[^/]+(?=\.gz)' /home/jupyter/7d9e54ed-929c-40ac-ac16-e40e98e398c1_gz_files_list.txt  | sort | uniq




gsutil -m -u "$workspace_id" ls ${bucket_id}cromwell-execution/RunPipeline/7d9e54ed-929c-40ac-ac16-e40e98e398c1/call-CNVcalling/shard-*/**/*.gz > 7d9e54ed-929c-40ac-ac16-e40e98e398c1_gz_files_list.txt


gsutil -m -u $workspace_id cp ${bucket_id}cromwell-execution/RunPipeline/*/call-CNVcalling/shard-*/**/*.gz

gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv_quantisnp/ | wc -l

gs://fc-secure-07cc78a0-0516-4955-a996-899e381a51f9/meta/GSA.sampleManifestAncestry.txt




workspace_id="terra-vpc-sc-a6da20a8"
bucket_id="gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/"
gsutil -m -u $workspace_id cp ${bucket_id}cnvcalling_inputs/sex_file.tsv.gz .
gsutil -m -u $workspace_id cp ${bucket_id}cnvcalling_inputs/sex_file_from_autocallGender.tsv.gz .





for sample_id in 1001355 1036867 1052179 1090122 1091623 1116947 1145516 1159666 1183030 1187687 1195452 1225778 1256163 1265358 1269444; do
    echo ${sample_id}
    directory=$(grep ${sample_id} *_gz_files_list.txt | head -n 1 | cut -d':' -f2,3  | sed 's|/[^/]*$||')
    echo $directory
    gsutil -m -u $workspace_id rm -r ${directory}
    gsutil -m -u $workspace_id rm ${bucket_id}output_cnvcalling/penncnv_quantisnp/${sample_id}.*
done




grep "7d9e54ed-929c-40ac-ac16-e40e98e398c1" merged_gz_files_list.txt | grep ".quanti.qc.gz" | grep -oP '(?<=/)[^/]+(?=\.gz)' | sort | uniq -d | wc -l

for file in $(grep ".quanti.qc.gz" merged_gz_files_list.txt | grep -oP '(?<=/)[^/]+(?=\.gz)' | sort | uniq -d);do
    echo $file
    # grep $file merged_gz_files_list.txt | wc -l 
    # grep $file merged_gz_files_list.txt | sed 's|/[^/]*$||'
    for file2 in $(grep $file merged_gz_files_list.txt); do 
     directory=$(echo $file2 | sed 's|/[^/]*$||' )
     gsutil -m -u $workspace_id cat ${directory}/CNVcalling-* | tail -n 2
    #  gsutil -m -u $workspace_id cat $file2 | zcat | wc -l
    done
done


path=gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/7d9e54ed-929c-40ac-ac16-e40e98e398c1/call-CNVcalling/shard-7565/
gsutil -m -u $workspace_id cp ${path}*.quanti.qc.gz .

gsutil -m -u $workspace_id ls ${path}
gsutil -m -u $workspace_id cat ${path}CNVcalling-7565.log

gsutil -m -u $workspace_id ls -la gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/7d9e54ed-929c-40ac-ac16-e40e98e398c1/call-CNVcalling/shard-7565/attempt-2/1216252.quanti.cnv.gz


grep -oP '(?<=/)[^/]+(?=\.gz)' /home/jupyter/7d9e54ed-929c-40ac-ac16-e40e98e398c1_gz_files_list.txt  | sort | uniq


grep ".quanti.qc.gz" merged_gz_files_list.txt | grep attempt | wc -l

path=gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cromwell-execution/RunPipeline/7d9e54ed-929c-40ac-ac16-e40e98e398c1/call-CNVcalling/shard-10462/
gsutil -m -u $workspace_id ls ${path}
gsutil -m -u $workspace_id cat ${path}CNVcalling-10462.log | tail


jobid="7d9e54ed-929c-40ac-ac16-e40e98e398c1"
### Remove those that don't have the "Done delocalization." tag at the end of the ${path}CNVcalling-*.log
gsutil -m -u "$workspace_id" ls ${bucket_id}cromwell-execution/RunPipeline/${jobid}/call-CNVcalling/shard-*/**/CNVcalling*.log > ${jobid}_log_file.txt


grep "01\/16\/2025" 6008022_lyon.txt | wc -l


grep "01\/16\/2025" *.txt | grep "\.sh" > test
grep "01\/16\/2025" *.txt | wc -l


grep "attempt" 7d9e54ed-929c-40ac-ac16-e40e98e398c1_log_file.txt |wc -l



for file in $(grep "." 7d9e54ed-929c-40ac-ac16-e40e98e398c1_log_file.txt);do
    echo "$file"
    if [ $(gsutil -m -u "$workspace_id" cat "$file" | tail -n 1 | grep "Done delocalization" | wc -l) -eq 1 ]; then
        echo "to not remove dir"
    else
        directory=$(echo $file | sed 's|/[^/]*$||')
        echo "to remove ${directory}"
        gsutil -m -u "$workspace_id" rm -r ${directory}
    fi
done




    # grep $file merged_gz_files_list.txt | wc -l 
    # grep $file merged_gz_files_list.txt | sed 's|/[^/]*$||'
    for file2 in $(grep $file merged_gz_files_list.txt); do 
     directory=$(echo $file2 | sed 's|/[^/]*$||' )
     gsutil -m -u $workspace_id cat ${directory}/CNVcalling-* | tail -n 2
    #  gsutil -m -u $workspace_id cat $file2 | zcat | wc -l
    done
done


workspace_id="terra-vpc-sc-a6da20a8"
bucket_id="gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/"

gsutil -m -u $workspace_id cp gs://fc-aou-datasets-controlled/v8/microarray/vcf/manifest.csv manifest_tier_v8.csv

awk -F',' '{OFS="\t"; $1=$1; print}' manifest_tier_v8.csv | cut -f1,2 > manifest_tier_v8.tsv


gsutil -m -u $workspace_id  cp "gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/cnvcalling_inputs/manifest.tsv.gz" manifest_tier_v7.tsv.gz



### Find those to redo:
# Identifier ceux qu'il faut importer : Même VCF

workspace_id="terra-vpc-sc-a6da20a8"
bucket_id="gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/"



gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv_quantisnp/*.quanti.qc.gz | sed -E 's|.*/([^.]+)\.quanti.*|\1|' > sample_processed.txt



gsutil -m -u $workspace_id cp ${bucket_id}cnvcalling_inputs/sex_file_from_autocallGender_tier_v7.tsv.gz .





# Output file
output_file="sex_info.tsv"
echo -e "Sample_ID\tGender_Source\tGender_Info" > $output_file

# Process each sample ID
for sample_id in $(cat person_id_to_find_autocallGender); do
    echo $sample_id
    # Extract VCF path
    path_to_vcf=$(zgrep -P "${sample_id}\t" manifest_tier_v8.tsv.gz | cut -f2)
    
    # Extract gender-related lines from the VCF
    vcf_gender_info=$(gsutil -m -u $workspace_id cat ${path_to_vcf} | head -n 5000 | zcat | grep Gender)
    
    # Format and write output
    echo -e "$vcf_gender_info" | while read line; do
        echo -e "${sample_id}\tVCF\t${line}" >> $output_file
    done
done



gsutil -m -u $workspace_id cat gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/output_cnvcalling/penncnv_quantisnp/1000004.pcnv.log.gz | zcat | head
1000004



gsutil -m -u "$workspace_id" ls ${bucket_id}cromwell-execution/RunPipeline/**/*.gz > gz_files_list.txt

awk -F'/' '{print $NF}' gz_files_list.txt > output.txt



gsutil -m -u $workspace_id cat gs://fc-secure-07cc78a0-0516-4955-a996-899e381a51f9/meta/GSA.sampleManifest.v8.Ancestry.txt | head




workspace_id="terra-vpc-sc-bf15287e"
bucket_id="gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/"
gsutil -m -u $workspace_id cp ${bucket_id}cnvcalling_inputs/manifest_tier_v8.tsv.gz .
gsutil -m -u $workspace_id cp ${bucket_id}cnvcalling_inputs/sex_file_v8.tsv.gz .


zcat manifest_tier_v8.tsv.gz | head 
zcat manifest_tier_v7.tsv.gz | head 


comm -12 <(zcat manifest_tier_v7.tsv.gz | cut -f1 | sort -u) <(zcat manifest_tier_v8.tsv.gz | cut -f1 | sort -u) > common_id_v7_v8.txt

grep -Fwf common_id_v7_v8.txt <(zcat manifest_tier_v7.tsv.gz) > v7_filtered.txt
grep -Fwf common_id_v7_v8.txt <(zcat manifest_tier_v8.tsv.gz) > v8_filtered.txt

sort -k1,1 v7_filtered.txt > sorted1.txt
sort -k1,1 v8_filtered.txt > sorted2.txt
join -t$'\t' -1 1 -2 1 sorted1.txt sorted2.txt > merged.txt

cut -f2 merged.txt > col2.txt
cut -f3 merged.txt > col3.txt
diff col2.txt col3.txt



zgrep -E "male|female" sex_file_v8.tsv.gz | awk '{print $1}' | shuf | sort | comm -23 - <(sort sample_processed.txt) | head -n 35000 > sample_ids_file_n35000.txt

gsutil -m -u $workspace_id cp sample_ids_file_n35000.txt ${bucket_id}cnvcalling_inputs/to_run/


gsutil -m -u $workspace_id cat gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cromwell-execution/RunPipeline/97430a36-e3cd-4b95-a5dc-8a0124b8b950/call-CNVcalling/shard-0/1000000.baf_lrr_perchr.tsv.gz | zcat





zgrep -E "male|female" sex_file_v8.tsv.gz | awk '{print $1}' | shuf | sort | comm -23 - <(sort sample_processed.txt) | comm -23 - <(sort sample_ids_file_n35000.txt) > all_sample_to_run.txt

split -l 35000 -d all_sample_to_run.txt chunk_

gsutil -m -u $workspace_id cat gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cnvcalling_inputs/random_1000.txt | head

gsutil -m -u $workspace_id cp gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cnvcalling_inputs/random_1000.txt .

<(zcat manifest_tier_v8.tsv.gz | awk '{print $1}')

comm -23 <(sort random_1000.txt) <(zcat manifest_tier_v7.tsv.gz | awk '{print $1}' | sort)

zcat manifest_tier_v7.tsv.gz > manifest_tier_v7.tsv
zcat manifest_tier_v8.tsv.gz > manifest_tier_v8.tsv
join -t$'\t' -1 1 -2 1 random_1000.txt manifest_tier_v7.tsv > merged.txt

awk 'NR==FNR {ids[$1]; next} $1 in ids' random_1000.txt manifest_tier_v7.tsv



gsutil -m -u $workspace_id cp random_1000_for_pfb.tsv ${bucket_id}cnvcalling_inputs/


cromshell status 886067bd-1386-4462-9003-3c2defaf913b

#############
#############
#############


#############
#############
#############

gsutil 

gsutil -m -u $workspace_id cat gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cnvcalling_inputs/to_run/chunk_00 | wc -l

gsutil -m -u $workspace_id cat gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cromwell-execution/RunPipeline/886067bd-1386-4462-9003-3c2defaf913b/call-CNVcalling/shard-30385/attempt-2/CNVcalling-30385.log


gsutil -m -u $workspace_id cat gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cromwell-execution/RunPipeline/886067bd-1386-4462-9003-3c2defaf913b/call-CNVcalling/shard-10001/stdout

ling/shard-12445/attempt-2/1336797.baf_lrr_perchr.tsv.gz...
Removing gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cromwell-execution/RunPipeline/886067bd-1386-4462-9003-3c2defaf913b/call-CNVcalling/shard-12445/attempt-2/1336797.pcnv.chrx.log.gz...
Removing gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cromwell-execution/RunPipeline/886067bd-1386-4462-9003-3c2defaf913b/call-CNVcalling/shard-12445/attempt-2/1336797.pcnv.chrx.out.gz...
Removing gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cromwell-execution/RunPipeline/886067bd-1386-4462-9003-3c2defaf913b/call-CNVcalling/shard-12445/attempt-2/1336797.pcnv.log.gz...
Removing gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cromwell-execution/RunPipeline/886067bd-1386-4462-9003-3c2defaf913b/call-CNVcalling/shard-12445/attempt-2/1336797.quanti.cnv.gz...
Removing gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cromwell-execution/RunPipeline/886067bd-1386-4462-9003-3c2defaf913b/call-CNVcalling/shard-12445/attempt-2/1336797.pcnv.out.gz...
Removing gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cromwell-execution/RunPipeline/886067bd-1386-4462-9003-3c2defaf913b/call-CNVcalling/shard-12445/attempt-2/1336797.quanti.qc.gz...
Removing gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cromwell-execution/RunPipeline/886067bd-1386-4462-9003-3c2defaf913b/call-CNVcalling/shard-12445/attempt-2/1336797.quanti.log.gz...
Removing gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187


## Reconstruct the directory path of those that are not complete
cat directories_status.tsv | awk -v bucket_id="$bucket_id" -v job_id="$job_id" '{
    split($1, a, /CNVcalling-|\.log/);
    print $0, "\t" bucket_id "cromwell-execution/RunPipeline/" job_id "/call-CNVcalling/shard-" a[2] "/" a[1]
}' | awk '{gsub(/no_attempt\//, ""); print}' > directories_all.tsv

comm -23 directories_with_gz.txt <(cut -f3 directories_all.tsv | sort)

gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/quantisnp/quanti/





gsutil -m -u $workspace_id ls gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/output_cnvcalling/penncnv_quantisnp/*.quanti.*.gz | head

gsutil -m -u $workspace_id cp gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/output_cnvcalling/penncnv_quantisnp/*.pcnv.*.gz ${bucket_id}output_cnvcalling/penncnv/


bucket_id="gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/"
gsutil -m -u $workspace_id cp gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/output_cnvcalling/penncnv_quantisnp/*.quanti.*.gz ${bucket_id}output_cnvcalling/quantisnp/




workspace_id="terra-vpc-sc-a6da20a8"

gsutil -m -u $workspace_id cp gs://fc-aou-datasets-controlled/v8/microarray/vcf/manifest.csv manifest.csv

gsutil -m -u $workspace_id cp gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cnvcalling_inputs/manifest_tier_v8.tsv.gz manifest_tier_v8.tsv.gz

comm -3 <(cat chunk_03 | sort) <(cut -d',' -f1 manifest.csv | sort)


comm -12 <(cat chunk_03 | sort) <(cut -d',' -f1 manifest.csv | sort) | wc -l


comm -12 <(zcat manifest_tier_v8.tsv.gz | cut -f1 | sort) <(cut -d',' -f1 manifest.csv | sort) | wc -l


for file in chunk*; do
    if [[ -f "$file" ]]; then
        echo "$file"
        wc -l $file
        #comm -12 <(zcat manifest_tier_v8.tsv.gz | cut -f1 | sort) <(cat $file | sort) > ${file}_v2
    fi
done



workspace_id="terra-vpc-sc-bf15287e"
bucket_id="gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/"

fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cnvcalling_inputs/to_run
gsutil -m -u $workspace_id cp chunk_*_v2 gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cnvcalling_inputs/to_run/


gsutil -m -u $workspace_id rm -r gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cromwell-execution/RunPipeline/fb6944a5-2e56-4828-84e7-8b4de512250d/





gsutil -m -u $workspace_id cp gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cnvcalling_inputs/* .



workspace_id="terra-vpc-sc-bf15287e"
gsutil -m -u $workspace_id cat gs://fc-secure-99703ba9-34f9-43ad-8b7d-287be45dbe12/README.md



workspace_id="terra-vpc-sc-xxxxxx"
bucket_id="gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/"
gsutil -m -u $workspace_id cp ${bucket_id}output_cnvcalling/cnv_annotation/README.md .

gsutil -m -u $workspace_id cat ${bucket_id}README.md
gsutil -m -u $workspace_id cat ${bucket_id}output_cnvcalling/cnv_annotation/inputs/penncnv_qc_score.tsv | wc -l
gsutil -m -u $workspace_id cat ${bucket_id}output_cnvcalling/cnv_annotation/inputs/penncnv_qc_score.tsv | head
gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/cn
gsutil -m -u $workspace_id rm -r ${bucket_id}output_cnvcalling/penncnv/x/



gsutil -m -u $workspace_id cp ${bucket_id}cromwell-execution/RunPipeline/${job_id}/**/*.quanti.*.gz ${bucket_id}output_cnvcalling/quantisnp/


## Check the number
timestamp=$(date "+%Y-%m-%d %H:%M:%S")
echo "Timestamp: $timestamp"

workspace_id="terra-vpc-sc-bf15287e"
bucket_id="gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/"

## Check the number
baf_lrr_count=$(gsutil -m -u "$workspace_id" ls "${bucket_id}output_cnvcalling/baf_lrr/" | wc -l)
echo "baf_lrr: $baf_lrr_count"

quantisnp_count=$(gsutil -m -u "$workspace_id" ls "${bucket_id}output_cnvcalling/quantisnp/" | wc -l)
quantisnp_divided=$((quantisnp_count / 4))
echo "quantisnp (divided by 4): $quantisnp_divided"

penncnv_count=$(gsutil -m -u "$workspace_id" ls "${bucket_id}output_cnvcalling/penncnv/" | wc -l)
penncnv_divided=$((penncnv_count / 4))
echo "penncnv (divided by 4): $penncnv_divided"


timestamp=$(date "+%Y-%m-%d %H:%M:%S")
echo "Timestamp: $timestamp"


OUTPUT_FILE="penncnv_qc_score.tsv" # Output TSV file

# Write the header to the output file
echo -e "IID\tLRR_mean\tLRR_median\tLRR_SD\tBAF_mean\tBAF_median\tBAF_SD\tBAF_DRIFT\tGCWF\tWF" > "$OUTPUT_FILE"

gsutil -m -u "$workspace_id" cat ${bucket_id}output_cnvcalling/penncnv/*.pcnv.log.gz | zcat \
    | grep "quality summary" \
    | sed -r 's/NOTICE: quality summary for IID_//g; \
              s/_vcf.tsvlite://g; \
              s/LRR_mean=//g; \
              s/LRR_median=//g; \
              s/LRR_SD=//g; \
              s/BAF_mean=//g; \
              s/BAF_median=//g; \
              s/BAF_SD=//g; \
              s/BAF_DRIFT=//g; \
              s/GCWF=//g; \
              s/WF=//g' \
    | sed -r 's/\s+/\t/g' >> "$OUTPUT_FILE"


# Extract relevant data from compressed log files, format it, and append to output file
zcat "$INPUT_DIR"/*.pcnv.log.gz \
    | grep "quality summary" \
    | sed -r 's/NOTICE: quality summary for IID_//g; \
              s/_vcf.tsvlite://g; \
              s/LRR_mean=//g; \
              s/LRR_median=//g; \
              s/LRR_SD=//g; \
              s/BAF_mean=//g; \
              s/BAF_median=//g; \
              s/BAF_SD=//g; \
              s/BAF_DRIFT=//g; \
              s/GCWF=//g; \
              s/WF=//g' \
    | sed -r 's/\s+/\t/g' >> "$OUTPUT_FILE"




mkdir penncnv
gsutil -m -u $workspace_id cp penncnv_qc_score.tsv ${bucket_id}output_cnvcalling/
gsutil -m -u $workspace_id cp plink_callrate.tsv ${bucket_id}output_cnvcalling/


workspace_id="terra-vpc-sc-bf15287e"
bucket_id="gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/"

gsutil -m -u $workspace_id cat ${bucket_id}output_cnvcalling/cnv_annotation/inputs | head


gsutil -m -u $workspace_id cp ${bucket_id}output_cnvcalling/plink_callrate.tsv .

workspace_id="terra-vpc-sc-bf15287e"
bucket_id="gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/"
gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/penncnv/ | head 


gsutil -m -u $workspace_id cat gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/output_cnvcalling/penncnv/*.pcnv.out.gz > merged.pcnv.out.gz

gsutil -m -u $workspace_id ls ${bucket_id}output_cnvcalling/



gsutil -m -u $workspace_id cp ${bucket_id}README.md .



gsutil -m -u $workspace_id cp gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cnvcalling_inputs/sex_file_v8.tsv.gz .


cut -f2 filtered_sexfile.tsv| sort | uniq -c



comm -12 <(zcat manifest_tier_v8.tsv.gz | cut -f1 | sort) <(zcat sex_file_v8.tsv.gz | cut -f1 | sort)

comm -12 <(zcat manifest_tier_v8.tsv.gz | cut -f1 | sort) <(zcat sex_file_v8.tsv.gz | cut -f1 | sort) | grep -Ff - <(zcat sex_file_v8.tsv.gz) > filtered sexfile.tsv
zcat sex_file_v8.tsv.gz | grep -Ff <(zcat manifest_tier_v8.tsv.gz | cut -f1 | sort -u) > filtered_sexfile.tsv

zgrep 1000061 sex_file_v8.tsv.gz
grep 1000061 filtered_sexfile.tsv



gsutil -m -u $workspace_id cat gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cnvcalling_inputs/sex_file_v8.tsv.gz | zcat | head


gsutil -m -u $workspace_id cp filtered_sexfile.tsv.gz gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cnvcalling_inputs/sex_file_v8.tsv.gz


cnv_annotation
inputs
scripts



jupyter@55e00d5e7dd6:~/cnv_annotation$ ls inputs/


gsutil -m -u $workspace_id cp -r cnv_annotation gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/output_cnvcalling/
gsutil -m -u $workspace_id ls gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/output_cnvcalling/



fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/cnvcalling_inputs

gsutil -m -u $workspace_id ls gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/output_cnvcalling/

sample=1000000


for sample in $(zcat sex_file_v8.tsv.gz | cut -f1);do
echo $sample
zgrep $sample sex_file_v8.tsv.gz | cut -f2
gsutil -m -u $workspace_id cat gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/output_cnvcalling/penncnv/${sample}.pcnv.chrx.log.gz | zcat  | sed -n '7p' | awk -F "is predicted | based" '{if (NF>1) print $2}'

#gsutil -m -u $workspace_id cat gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/output_cnvcalling/quantisnp/${sample}.quanti.qc.gz | zcat | sed -n '2p' | cut -f7
done
