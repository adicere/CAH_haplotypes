import pandas as pd
import numpy as np 
import os 
import sys 
import argparse 
import subprocess

updater = argparse.ArgumentParser(
    prog='Update Amplicon Annotation',
    description='Adds to annotation file a new column(s) with filepaths to corresponding output files after using naive calling (bcftools mpileup)'
)

updater.add_argument('ann_file', type=str, help='provide a path to the table with filepaths of processed amplicons')
# updater.add_argument('-o', '--output', type=str, help='provide a path to the output directory to store the mpileup results', required = True, metavar = 'DIR')
updater.add_argument('-e', '--exclude', type=str, help='provide a txt file with samples that need to be excluded', metavar = 'FILE')

args=updater.parse_args()

ann_file = pd.read_excel(args.ann_file)


def exclude_samples(ann_file, sample_list):
    samples=sample_list.read().split('\n')
    ann_file['patient'] = ann_file['patient'].astype(str)
    filtered = ann_file[~ann_file['patient'].isin(samples)]
    
    return None

def naive_calling(ann_file):
    pass

if args.exclude is not None:
    to_exclude = open(args.exclude) 
    exclude_samples(ann_file, to_exclude)