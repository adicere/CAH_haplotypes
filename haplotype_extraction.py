import pandas as pd
import numpy as np
import argparse
import os 
from pathlib import Path
import subprocess

starter = argparse.ArgumentParser(prog='CYP21A2 haplotyping tool', 
                                  description='Haplotype extraction for single and family samples based on CYP21A2 amplicon sequencing')

input_group = starter.add_mutually_exclusive_group(required=True)
input_group.add_argument('-l', '--list', type=str, help='Provide .tsv/.csv list with full paths to .vcf.gz files to analyze and extract haplotypes')
input_group.add_argument('-d', '--dir', type=str, help='Provide a directory with .vcf.gz files to analyze and extract and extract haplotypes')
starter.add_argument('-o', '--output', type=str, help='Provide an output directory for results', required=True)

def modify_pos(file: str | Path) -> None:
    vcf = pd.read_csv(file, sep='\t', comment='#', header=None, index_col=False)  
    vcf[0] = 'chr6'
    vcf[1] = int(32037619) + vcf[1].astype(int)
    vcf.to_csv(f"{args.output}/tmp/body.txt", header=None, index=None, sep='\t')

def normalize_coordinates(vcf: str | Path, filename: str) -> None:
    header = f"{args.output}/tmp/header.txt"
    body = f"{args.output}/tmp/body.txt"
    out_vcf = f"{args.output}/recalibrated/{filename}.recalibrated.vcf.gz"

    subprocess.run(['bcftools', 'view', '-h', vcf, '-o', header], check=True)
    modify_pos(vcf)

    with open(out_vcf, 'wb') as vcf:
        p = subprocess.Popen(['bgzip'], stdin=subprocess.PIPE, stdout=vcf)

        with open(header, 'rb') as h, open(body, 'rb') as b:
            p.stdin.write(h.read())
            p.stdin.write(b.read())
            
        p.stdin.close() 
        p.wait()

        if p.returncode!=0:
            raise RuntimeError('bgzip failed')
        
    subprocess.run(['tabix', f"{o_dir}/recalibrated/{filename}.recalibrated.vcf.gz"], check=True)


args = starter.parse_args()

if args.dir:
    i_dir = os.path.abspath(args.dir)
o_dir = os.path.abspath(args.output)

folders = [f"{o_dir}/recalibrated", f"{o_dir}/results", f"{o_dir}/tmp"]
for dir in folders:
    os.makedirs(dir, exist_ok=True)

if args.list:
    with open(args.list) as f:
        files=pd.read_csv(f, sep='\t')
    
    for idx, row in files.iterrows():
        normalize_coordinates(row['vcf'], row['patient'])
    

if args.dir:
    for file in args.dir.glob('*.phased.vcf.gz'):
        if not file.is_file():
            continue
        
        name = os.path.basename(file)
        sample = name.split('.')[0]
        normalize_coordinates(file, sample)

def main() -> None:
    pass

if __name__ == '__main__':
    main()