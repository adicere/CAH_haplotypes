#!/bin/bash

set -euo pipefail

while getopts 's:r:o:' flag; do
    case "${flag}" in
        s) 
            sample_list=${OPTARG} ;; # litt of VCF files to merge
        r)
            rename_list=${OPTARG} ;; # file with new sample names in the form of "old_name \t new_name"
        o) 
            o_dir=${OPTARG} ;; # directory for output files
    esac 
done

bcftools merge -l "${sample_list}" -0 -m none -Oz -o "${o_dir}merged_raw.vcf.gz"
bcftools index -f -t "${o_dir}merged_raw.vcf.gz"

bcftools reheader -s "${rename_list}" -o "${o_dir}merged.vcf.gz" "${o_dir}merged_raw.vcf.gz"
bcftools index -f -t "${o_dir}merged.vcf.gz"

bcftools annotate -x ^FORMAT/GT "${o_dir}merged.vcf.gz" | bcftools view -t chr6:32037643-32041345 | bcftools norm -m -snps -a --atom-overlaps . -Oz -o "${o_dir}splitted_genotype.vcf.gz"
bcftools index -f -t "${o_dir}splitted_genotype.vcf.gz"