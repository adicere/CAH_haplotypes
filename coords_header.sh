#!/bin/bash

: ' 
Script combines recalibration of variant coordinates and add changes to the header to have right sample id 
'

set -euo pipefail

while getopts 'h:d:s:o:' flag; do
    case "${flag}" in 
        h) 
            hc_file=${OPTARG} ;;
        d)
            dv_file=${OPTARG} ;;
        s)
            sample=${OPTARG} ;;
        o)
            output_dir=${OPTARG} ;;
    esac
done

current_dir="$(dirname "$(realpath "$0")")"

#modifying the HaplotypeCaller VCF file
bcftools view -h "${hc_file}" > hc_header.txt
sed -i "s/pool/$sample/" hc_header.txt
python "$current_dir/modify_vcf.py" "$hc_file" hc_temp.txt
cat hc_header.txt hc_temp.txt | bgzip > "$output_dir/${sample}.hc.rehead.vcf.gz" && tabix "$output_dir/${sample}.hc.rehead.vcf.gz"

#modifying the DeepVariant VCF file
bcftools view -h "${dv_file}" > dv_header.txt
sed -i "s/pool/$sample/" dv_header.txt
python "$current_dir/modify_vcf.py" "$dv_file" dv_temp.txt
cat dv_header.txt dv_temp.txt | bgzip > "$output_dir/${sample}.dv.rehead.vcf.gz" && tabix "$output_dir/${sample}.dv.rehead.vcf.gz"

rm hc_header.txt hc_temp.txt dv_header.txt dv_temp.txt 
