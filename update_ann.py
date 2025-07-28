import pandas as pd
import numpy as np 
import os 
import sys 
import argparse 
import subprocess
import glob

bcftools_path = '/home/rutkovskaya.ea/miniforge3/envs/haplotypes/bin/bcftools'
reference = '/home/rutkovskaya.ea/haplotypes/CYP21A2-amp_hg38.fasta'
modificator='/home/rutkovskaya.ea/haplotypes/CAH_haplotypes/reheader.sh'


updater = argparse.ArgumentParser(
    prog='Update Amplicon Annotation',
    description=(
        """ Adds to annotation file a new column(s) with filepaths to corresponding output files after using naive calling (bcftools mpileup).
Results of naive calling would be stored at $OUTPUT_DIR/mpileup, all corrected by coordinates files at $OUTPUT_DIR/reheaded_vcf. 
        """
    ),
    formatter_class=argparse.RawTextHelpFormatter
)

updater.add_argument('ann_file', type=str, help='provide a path to the table with filepaths of processed amplicons')
updater.add_argument('-o', '--output', type=str, help='provide a path to the output directory to store the results', required = True, metavar = 'OUTPUT_DIR')
updater.add_argument('-e', '--exclude', type=str, help='provide a txt file with samples that need to be excluded', metavar = 'FILE')

args=updater.parse_args()

# amplicon_ann = pd.read_csv(args.ann_file, sep='\t', index_col=0)
amplicon_ann = pd.read_excel(args.ann_file, index_col=0)
# short=amplicon_ann.copy().iloc[:5, :]
mpileup_dir=os.path.join(args.output, 'mpileup')  
reheaded_dir = os.path.join(args.output, 'reheaded_vcfs')  
os.makedirs(mpileup_dir, exist_ok=True)   
os.makedirs(reheaded_dir, exist_ok=True)   

def exclude_samples(ann_file, sample_list):
    samples=sample_list.read().split('\n')
    updated = ann_file[~ann_file['patient'].isna()]
    updated['patient'] = updated['patient'].astype(int).astype(str)
    filtered = updated[~updated['patient'].isin(samples)]

    return filtered

def modify_vcf(ann_file):

    """
    Function recalibrates the coordinates of amplicon VCF files of corresponding caller,
    changes pool to sample_id in the header and writes paths to new files into
    original annotation
    """

    ann_file['hc_modified'] = ''
    ann_file['dv_modified'] = ''
    for idx, row in ann_file.iterrows():
        sample = str(row['orig_sample'])[:8]
        hc_vcf = row['hc_vcf']
        dv_vcf = row['deepvariant']
        subprocess.run(['bash', modificator, 
                        '-f', hc_vcf, 
                        '-t', 'hc', 
                        '-s', sample, 
                        '-o', reheaded_dir])
        subprocess.run(['bash', modificator, 
                        '-f', dv_vcf, 
                        '-t', 'dv', 
                        '-s', sample, 
                        '-o', reheaded_dir])
        ann_file.at[idx, 'hc_modified'] = f'{reheaded_dir}/{sample}.hc.rehead.vcf.gz'
        ann_file.at[idx, 'dv_modified'] = f'{reheaded_dir}/{sample}.dv.rehead.vcf.gz'
    
    return ann_file


def naive_calling(ann_file):

    """
    Function performs naive calling on bam files using bcftools mpileup+call, then
    modifies resulted VCF like in modify_vcf()
    """

    ann_file['mpileup'] = ''
    for idx, row in ann_file.iterrows():
        bam = row['bam']
        sample = str(row['orig_sample'])[:8]
        filename=f'{mpileup_dir}/{sample}.naive.vcf.gz'
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
        subprocess.run(['bash', modificator, 
                        '-f', filename, 
                        '-t', 'mpileup', 
                        '-s', sample, 
                        '-o', reheaded_dir])

        ann_file.at[idx, 'mpileup_modified'] = f'{reheaded_dir}/{sample}.mpileup.rehead.vcf.gz'
    
    return res
   

if args.exclude is not None:
    to_exclude = open(args.exclude) 
    filtered_list=exclude_samples(amplicon_ann, to_exclude)
else:
    res = modify_vcf(amplicon_ann)
    res=naive_calling(amplicon_ann)

res.to_excel(args.output+'modified_annotation.xlsx')