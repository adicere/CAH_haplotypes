import pandas as pd
import argparse
import subprocess
import sys
import os
import numpy as np
import pysam 
import gzip
import re
import json

# parser = argparse.ArgumentParser(description='Transforms .vcf to .xlsx format')
# parser.add_argument('vcf_file', type=str, help='Path to .vcf file')
# parser.add_argument('-o', '--output_dir', type=str, help='Path to the output directory')
# args = parser.parse_args()

vcf_in = sys.argv[1] #merged vcf file
var_table = pd.read_excel(sys.argv[2]) #preprocessed excel file with variant's annotation (should contain POS, REF and ALT columns) for matching
json_file = sys.argv[3] #json with average VAF of each sample 
o_dir=str(sys.argv[4]) #output directory to store files

with open(json_file, 'rt') as f:
    vafs = json.load(f)

#preparing vcf file to further annotation
def vcf_to_table(vcf_file):
    #tranform vcf format into df - skipping header info and renaming columns
    num_header = 0
    with gzip.open(vcf_file, 'rt') as f:
        for line in f.readlines():
            if line.startswith("##"):
                num_header += 1
            else:
                break
    vcf = pd.read_csv(vcf_file, sep="\t", skiprows=num_header,  compression='gzip')
    vcf = vcf.rename({"#CHROM": "CHROM"}, axis=1)

    vcf_copy = vcf.copy()
    vcf_copy = vcf_copy.drop(vcf_copy.iloc[:, 9:], axis = 1)

    #splitting samples genotypes into chromosomes
    for col in vcf.iloc[:, 9:].columns:
        new_cols = vcf[col].str.split(r"[\/|]", expand=True)
        new_cols.columns = [f"{col}_1", f"{col}_2"]
        vcf_copy=vcf_copy.join(new_cols)
    
    vcf_copy=vcf_copy.drop(['ID', 'FILTER', 'INFO', 'FORMAT'], axis=1)
    for col in vcf_copy.columns:
        if col.startswith('770'):
            vcf_copy[col] = vcf_copy[col].replace('.', None).astype('Int64')

    return AF(vcf_copy)

#count allele frequency for each variant
def AF(df):
    subset=df.iloc[:,5:].T.groupby(lambda x: x.split('_')[0]).sum().T.replace(2,1)
    subset['AF'] = subset.sum(axis=1)
    df['AF'] = subset['AF']/(len(subset.columns)-1)
    return df

#matching result df with variant table and extracting corresponding columns (rsid, phenotype, ACMG classification) as table_annotated.xlsx
def variant_annotation(table, variants):
    variants['POS'] = variants['POS'].fillna(0).astype(int)
    matched = pd.merge(table,
    variants[['POS', 'REF', 'ALT', 'Effect', 'rs#', 'Allele associated phenotype', 'ACMG']],  on=['POS', 'REF', 'ALT'],  how='left') 
    matched.to_excel(os.path.join(o_dir+'table_annotated.xlsx'), index=False)
    return matched


#creating a separate excel file with added average VAF as VAF_annotated.xlsx
def add_VAF(df,vafs):
    VAF=[]
    for col in df.columns:
        if col.startswith('770'):
            res = re.search(r'([0-9]+)_', col)
            VAF.append(vafs[res.group(1)])
        else:
            VAF.append('')   
    df.loc[len(df)] = VAF
    df.to_excel(os.path.join(o_dir+'VAF_annotated.xlsx'), index=False)


phased = vcf_to_table(vcf_in)
annotated = variant_annotation(phased, var_table)
add_VAF(annotated, vafs)
