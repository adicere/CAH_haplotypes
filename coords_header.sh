#!/bin/bash

: ' 
Script combines recalibration of variant coordinates and add changes to the header to have right sample id 
'

set -euo pipefail

files=$1 # .txt list with paths to the samples
dir_out=$2 # output directory for processed vcf's

while IFS=$'\t' read -r hc dv; do
    id=`echo $hc | grep -o -E '[0-9]{8}'`
    bcftools view -h "$hc" > hc_header.txt
    sed -i "s/pool/$id/" hc_header.txt
    python /home/katerine/Desktop/VDKN/russian_field_of_experiments/scripts/Anton/process_vcf.py $hc hc_temp.txt
    cat hc_header.txt hc_temp.txt | bgzip > "$2${id}.haplotypecaller.phased.coord.rehead.vcf.gz" && tabix "$2${id}.haplotypecaller.phased.coord.rehead.vcf.gz"

    bcftools view -h "$dv" > dv_header.txt  
    sed -i "s/pool/$id/" dv_header.txt 
    python /home/katerine/Desktop/VDKN/russian_field_of_experiments/scripts/Anton/process_vcf.py $dv dv_temp.txt
    cat dv_header.txt dv_temp.txt | bgzip > "$2${id}.deepvariant.coord.rehead.vcf.gz" && tabix "$2${id}.deepvariant.coord.rehead.vcf.gz"

    rm hc_header.txt hc_temp.txt dv_header.txt dv_temp.txt 
done < "$files"

#creates a list for new .vcf.gz for the next step
find "$dir_out" -maxdepth 1 -type f -name "*.haplotypecaller.phased.coord.rehead.vcf.gz" > "$dir_out"/hc_gz_paths.txt
find "$dir_out" -maxdepth 1 -type f -name "*.deepvariant.coord.rehead.vcf.gz" > "$dir_out"/dv_gz_paths.txt  