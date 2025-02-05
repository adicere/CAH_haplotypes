#!/bin/bash

set -euo pipefail

: ' Script for merging vcf files and processing joint file to excel sheet
'

list=$1 # txt file with paths to vcf files
o_dir=$2 # output directory for a merged file

temp="${o_dir}temp"
res_dir="${o_dir}results"
qd_dir="${o_dir}QD_jsons"
mkdir $temp $res_dir $qd_dir

#filter out variants with super low quality and split multiallelic sites
files=$( cat "$list")
for f in $files
do
    filename=$(basename $f)
    basename=${filename%%.*}
    bcftools view -e 'QD < 1.5' $f | bcftools norm -m -snps --multi-overlaps . -a --atom-overlaps . -Oz -o $temp/"$basename"_splitted.vcf.gz
    #bcftools norm -m -snps -a --atom-overlaps . $f -Oz -o $temp/"$basename"_splitted.vcf.gz ## for server 
    tabix $temp/"$basename"_splitted.vcf.gz
done

#save paths to splitted files
find "$temp" -type f -name '*_splitted.vcf.gz' > "$temp"/splitted_files.txt

vcfs=$( cat "$temp"/splitted_files.txt)
for v in $vcfs
do
    python /home/katerine/Desktop/VDKN/russian_field_of_experiments/scripts/qual_calc.py $v $qd_dir
done

#merge all vcf files
bcftools merge --file-list $temp/"splitted_files.txt" -0 -m none -Oz -o $res_dir/merged.vcf.gz
tabix $res_dir/merged.vcf.gz
rm -rf "$temp"

#calculating average VAF for each sample - HANDLE PATHS
python /home/katerine/Desktop/VDKN/CAH_haplotypes/vaf.py $res_dir/merged.vcf.gz /media/katerine/XC2000/VDKN/results/

#removing unnecessary FORMAT fields except for GT
zcat $res_dir/merged.vcf.gz | bcftools annotate -x ^FORMAT/GT -Oz -o $res_dir/genotype_merged.vcf.gz
tabix $res_dir/genotype_merged.vcf.gz

#separate genotypes to chromosomes and annotation 
python /home/katerine/Desktop/VDKN/CAH_haplotypes/vcf_processing.py $res_dir/genotype_merged.vcf.gz \
/home/katerine/Desktop/VDKN/russian_field_of_experiments/scripts/variants.xlsx \
/media/katerine/XC2000/VDKN/results/avg_VAF.json \
/media/katerine/XC2000/VDKN/results/variant_VAF.json \
/home/katerine/Desktop/VDKN/russian_field_of_experiments/pseudo_variants.xlsx \
/media/katerine/XC2000/VDKN/results/



