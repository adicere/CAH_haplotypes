import pandas as pd
import numpy as np 
import os 
import sys 
import argparse 
import subprocess
import glob

bcftools_path = '/home/rutkovskaya.ea/miniforge3/envs/haplotypes/bin/bcftools'
reference = '/home/rutkovskaya.ea/haplotypes/CYP21A2-amp_hg38.fasta'

updater = argparse.ArgumentParser(
    prog='Update Amplicon Annotation',
    description='Adds to annotation file a new column(s) with filepaths to corresponding output files after using naive calling (bcftools mpileup)'
)

updater.add_argument('ann_file', type=str, help='provide a path to the table with filepaths of processed amplicons')
updater.add_argument('-o', '--output', type=str, help='provide a path to the output directory to store the mpileup results', required = True, metavar = 'DIR')
updater.add_argument('-e', '--exclude', type=str, help='provide a txt file with samples that need to be excluded', metavar = 'FILE')

args=updater.parse_args()

# os.mkdir(args.output + 'mpileup')
ann_file = pd.read_excel(args.ann_file)
short=ann_file.iloc[:5, :]
res_dir=args.output + 'mpileup'

def exclude_samples(ann_file, sample_list):
    samples=sample_list.read().split('\n')
    ann_file['patient'] = ann_file['patient'].astype(str)
    filtered = ann_file[~ann_file['patient'].isin(samples)]
    
    return filtered

def naive_calling(ann_file):
    ann_file['mpileup'] = ''
    for idx, row in ann_file.iterrows():
        bam = row['bam']
        sample = str(row['orig_sample'])
        filename=f'{res_dir}/{sample}.pileup.vcf.gz'
        mpileup = subprocess.Popen([bcftools_path, 
                                    'mpileup',
                                    '-f', reference, 
                                    '-a', 'FORMAT/AD,FORMAT/DP', 
                                    '--no-BAQ', '-d', '3000',
                                    '-Ou', bam], stdout=subprocess.PIPE)
        calling =  subprocess.Popen([bcftools_path, 'call',
                                     '-mv', '-Oz', 
                                     '-o', filename], stdin=mpileup.stdout, stdout=subprocess.DEVNULL)
        mpileup.stdout.close()
        calling.wait()           
        mpileup.wait() 
        ann_file.at[idx, 'mpileup'] = str(filename)
    print(ann_file)    

if args.exclude is not None:
    to_exclude = open(args.exclude) 
    filtered_list=exclude_samples(ann_file, to_exclude)
else:
    res=naive_calling(short)
