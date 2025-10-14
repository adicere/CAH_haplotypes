#!/bin/bash

: ' 
Script combines recalibration of variant coordinates and add changes to the header to have right sample id 
'

set -euo pipefail

while getopts 'f:t:s:o:' flag; do
    case "${flag}" in 
        f) 
            file=${OPTARG} ;;
        t)
            type=${OPTARG} ;;
        s)
            sample=${OPTARG} ;;
        o)
            output_dir=${OPTARG} ;;
    esac
done

output_dir="${output_dir%/}/" 

splitted="${output_dir}splitted_vcfs/"
reheaded="${output_dir}reheaded_vcfs/"

mkdir -p $splitted $reheaded 

current_dir="$(dirname "$(realpath "$0")")"

#splittin multiallelic sites

bcftools norm -m -snps -a --atom-overlaps . -Oz "${file}" -o "$splitted${sample}_splitted.${type}.vcf.gz" 
bcftools index -t "$splitted${sample}_splitted.${type}.vcf.gz" 

#modifying VCF file
bcftools view -h "${file}" > ${reheaded}header.txt
sed -i "s/pool/$sample/" ${reheaded}header.txt
python "$current_dir/modify_vcf.py" "$file" ${reheaded}temp.txt
cat ${reheaded}header.txt ${reheaded}temp.txt | bgzip > "${reheaded}${sample}.${type}.rehead.vcf.gz" && tabix "${reheaded}${sample}.${type}.rehead.vcf.gz"

rm ${reheaded}header.txt ${reheaded}temp.txt 

echo "$splitted${sample}_splitted.${type}.vcf.gz"
echo "${reheaded}${sample}.${type}.rehead.vcf.gz"