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

# echo "$output_dir${sample}.${type}.rehead.vcf.gz" 
current_dir="$(dirname "$(realpath "$0")")"

#modifying VCF file
bcftools view -h "${file}" > header.txt
sed -i "s/pool/$sample/" header.txt
python "$current_dir/modify_vcf.py" "$file" temp.txt
cat header.txt temp.txt | bgzip > "$output_dir/${sample}.${type}.rehead.vcf.gz" && tabix "$output_dir/${sample}.${type}.rehead.vcf.gz"

rm header.txt temp.txt 