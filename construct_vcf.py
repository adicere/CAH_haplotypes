import pandas as pd
import numpy as np
from collections import defaultdict
import sys 
import os
import logging
import argparse
import subprocess
import random

bcftools_path = '/primary/home/rutkovskaya.ea/miniforge3/envs/haplotypes/bin/bcftools'

constructor = argparse.ArgumentParser(prog='Building VCF files with validated variants',
                                      description='Constructing combined VCF files with variants from both HaplotypeCaller and DeepVariant')

constructor.add_argument('ann_file', type=str, help='provide a file with paths to VCF files of processed amplicons')
constructor.add_argument('-r', '--resolved_vars', type=str, help='.xlsx result file from variant validation after manual check',
                         required=True, metavar='FILE')
constructor.add_argument('-o', '--output', type=str, help='provide a path to the output directory to store modified VCFs', metavar='DIR')
args=constructor.parse_args()

joined_dir = os.path.join(args.output, 'joined_vcfs')
os.makedirs(joined_dir, exist_ok=True)

def open_vcf(link): 
    vcf = pd.read_csv(link, sep='\t', comment='#', header=None,
          names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], index_col=False)  
    
    return vcf

def fix_deepvariant(df):

    format_cols = df['FORMAT'].str.split(':')
    sample_vals = df['SAMPLE'].str.split(':')

    new_samples = []
    for fmt, val in zip(format_cols, sample_vals):
        fmt_map = dict(zip(fmt, val))


        new_sample = [
            fmt_map.get("GT", "."),
            fmt_map.get("AD", "."),
            fmt_map.get("DP", "."),
            fmt_map.get("GQ", "."),
            fmt_map.get("PL", ".")
        ]
        new_samples.append(':'.join(new_sample))

    df['FORMAT'] = "GT:AD:DP:GQ:PL"
    df['SAMPLE'] = new_samples

    return df

def fix_mpileup(df):

    format_cols = df['FORMAT'].str.split(':')
    sample_vals = df['SAMPLE'].str.split(':')

    new_samples = []
    for fmt, val in zip(format_cols, sample_vals):
        fmt_map = dict(zip(fmt, val))


        new_sample = [
            fmt_map.get("GT", "."),
            fmt_map.get("AD", "."),
            fmt_map.get("DP", "."),
            str(random.randint(30,51)),
            fmt_map.get("PL", ".")
        ]
        new_samples.append(':'.join(new_sample))


    # Set new values
    df['FORMAT'] = "GT:AD:DP:GQ:PL"
    df['SAMPLE'] = new_samples

    return df

def join_variants(sample, resolved, hc, dv, mpileup):
    sample_vars = resolved.loc[resolved[sample].notna(), sample].to_dict()
    for pvcf in sample_vars:
        if sample_vars[pvcf] == 'dv':
            pos = resolved.loc[pvcf, 'POS']
            alt = resolved.loc[pvcf, 'ALT']
            filtered = hc.loc[(hc['POS'] == pos) & (hc['ALT'] == alt)]
            if not filtered.empty:
                idx = filtered.index[0]
                hc.at[idx, 'SAMPLE'] = dv.loc[(dv['POS']==pos) & (dv['ALT'] == alt), 'SAMPLE'].values[0] 
            else:
                hc = pd.concat([hc, dv.loc[(dv['POS']==pos) & (dv['ALT'] == alt)]])
        if sample_vars[pvcf] == 'mp':
            pos = resolved.loc[pvcf, 'POS']
            alt = resolved.loc[pvcf, 'ALT']
            filtered = hc.loc[(hc['POS'] == pos) & (hc['ALT'] == alt)]
            if not filtered.empty:
                idx = filtered.index[0]
                hc.at[idx, 'SAMPLE'] = mpileup.loc[(mpileup['POS']==pos) & (mpileup['ALT'] == alt), 'SAMPLE'].values[0] 
            else:
                hc = pd.concat([hc, mpileup.loc[(mpileup['POS']==pos) & (mpileup['ALT'] == alt)]])
    
    hc = hc.sort_values(by='POS')

    return hc


def transform_into_vcf(joined_df, vcf, sample):
    with open(f'{joined_dir}/header.txt', 'w') as out:
        subprocess.run([bcftools_path, 'view', '-h', vcf],
                    stdout=out,
                    check=True)
        
    joined_df.to_csv(f'{joined_dir}/{sample}.txt', header=None, index=None, sep='\t')

    concat = subprocess.Popen(['cat', 
                               f'{joined_dir}/header.txt', 
                               f'{joined_dir}/{sample}.txt'],
                               stdout=subprocess.PIPE)
    zipping = subprocess.Popen(['bgzip', '-o',
                                f'{joined_dir}/{sample}_joined.vcf.gz'],
                                stdin=concat.stdout, stdout=subprocess.DEVNULL)

    concat.stdout.close()
    concat.wait()           
    zipping.wait()

    subprocess.run(['rm', f'{joined_dir}/header.txt', f'{joined_dir}/{sample}.txt'])

    return f'{joined_dir}/{sample}_joined.vcf.gz'

annotation = pd.read_excel(args.ann_file, index_col = 0)
validated_vars=pd.read_excel(args.resolved_vars, index_col=0)
validated_vars['POS'] =validated_vars.index.map(lambda x: int(x.split('_')[0]) - 32037619)
validated_vars['ALT'] =validated_vars.index.map(lambda x: x.split('_')[2])
annotation['for_phasing'] = np.nan 

for col in validated_vars.columns:
    if validated_vars[col].isin(['dv', 'mp']).any():
        hc = open_vcf(annotation.loc[col, 'hc_splitted'])
        dv = fix_deepvariant(open_vcf(annotation.loc[col, 'dv_splitted']))
        mpileup = fix_mpileup(open_vcf(annotation.loc[col, 'mp_splitted']))

        result = join_variants(col, validated_vars, hc, dv, mpileup)

        annotation.at[col, 'for_phasing']=transform_into_vcf(result, annotation.loc[col, 'hc_splitted'], col)


subprocess.run(f'for f in {joined_dir}/*.vcf.gz; do tabix -p vcf "$f"; done',
               shell=True,
               check=True)
        
annotation['for_phasing'] = annotation['for_phasing'].fillna(annotation['hc_vcf'])
##дроп ненужных колонок?

annotation.to_excel(args.output + 'modified_annotation.xlsx')