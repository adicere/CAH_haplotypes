#!/bin/bash

set -euo pipefail 

: ' Script for merging vcf files and processing joint file to excel sheet 
'

list=$1 # txt file with paths to vcf files
o_dir=$2 # output directory for a merged file 

temp="${o_dir}temp"
res_dir="${o_dir}results"
mkdir $temp $res_dir

#split multiallelic sites
files=$( cat "$list")
for f in $files 
do
filename=$(basename $f)
basename=${filename%%.*}
bcftools norm -m -snps --multi-overlaps . -a --atom-overlaps . $f -Oz -o $temp/"$basename"_splitted.vcf.gz
tabix $temp/"$basename"_splitted.vcf.gz
done

#save paths to splitted files
find "$temp" -type f -name '*_splitted.vcf.gz' > $temp/splitted_files.txt 

#merge all vcf files 
bcftools merge --file-list $temp/"splitted_files.txt" -0 -m snps -Oz -o $res_dir/merged.vcf.gz 
tabix $res_dir/merged.vcf.gz
rm -rf "$temp"

#calculating average VAF for each sample - HANDLE PATHS
python /home/katerine/Desktop/VDKN/CAH_haplotypes/vaf.py $res_dir/merged.vcf.gz '/home/katerine/Desktop/VDKN/russian_field_of_experiments/results/'

#removing unnecessary FORMAT fields except for GT
zcat $res_dir/merged.vcf.gz | bcftools annotate -x ^FORMAT/GT -Oz -o $res_dir/genotype_merged.vcf.gz 
tabix $res_dir/genotype_merged.vcf.gz 

#separate genotypes to chromosomes 
python /home/katerine/Desktop/VDKN/CAH_haplotypes/vcf_processing.py $res_dir/genotype_merged.vcf.gz \
/home/katerine/Desktop/VDKN/russian_field_of_experiments/scripts/variants.xlsx \
/home/katerine/Desktop/VDKN/russian_field_of_experiments/results/VAF.json \
/home/katerine/Desktop/VDKN/russian_field_of_experiments/results/ 



