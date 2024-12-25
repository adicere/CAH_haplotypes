#!/bin/bash

: ' 
Script combines calibrating coordinates and changing the header to have sample id instead of pool
'

set -euo pipefail

files=$1 # .txt list with paths to the samples
dir_out=$2 # output directory for processed vcf's
var=$( cat $files) 
for f in $var
do
    id=`echo $f | grep -o -E '[0-9]{8}'` #extracts the first 8 numbers of each provided file in a list 
    bcftools view -h "$f" > header.txt 
    sed -i "s/pool/$id/" header.txt 
    #bcftools view -h ../vcfs/770924920801.haplotype.caller.phased.vcf | sed "s/pool/77092492/" | grep -v -E '^##FORMAT=<ID=(GT|AD|DP|GQ|PL|PS),' > header.txt
    python Anton/process_vcf.py $f temp.txt
    cat header.txt temp.txt | bgzip > "$2${id}.haplotype.caller.phased.coord.rehead.vcf.gz" && tabix "$2${id}.haplotype.caller.phased.coord.rehead.vcf.gz"
    rm header.txt temp.txt 
done

find "$(pwd)" -type f -wholename "$dir_out/*.haplotype.caller.phased.coord.rehead.vcf.gz" > gzip_paths.txt #creates a list for new .vcf.gz for the next step