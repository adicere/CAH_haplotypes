#!/bin/bash

set -euo pipefail 

: ' Script for merging vcf files and processing joint file to excel sheet 
'

list=$1 # txt file with paths to vcf files
o_dir=$2 # output directory for a merged file 

#merge all vcf files 
bcftools merge --file-list "$list" -0 -m snps -Oz -o "$o_dir"/merged.vcf.gz
tabix "$o_dir"/merged.vcf.gz

#calculating average VAF for each sample - HANDLE PATHS
python /home/katerine/Desktop/VDKN/CAH_haplotypes/vaf.py "$o_dir"/merged.vcf.gz '/home/katerine/Desktop/VDKN/russian_field_of_experiments/'

#removing unnecessary FORMAT fields except for GT
zcat "$o_dir"/merged.vcf.gz | bcftools annotate -x ^FORMAT/GT -Oz -o "$o_dir"/genotype_merged.vcf.gz 
tabix "$o_dir"/genotype_merged.vcf.gz 

#separate genotypes to chromosomes 
python /home/katerine/Desktop/VDKN/CAH_haplotypes/vcf_processing.py "$o_dir"/genotype_merged.vcf.gz \
/home/katerine/Desktop/VDKN/russian_field_of_experiments/scripts/variants.xlsx \
/home/katerine/Desktop/VDKN/russian_field_of_experiments/VAF.json \
/home/katerine/Desktop/VDKN/russian_field_of_experiments/



