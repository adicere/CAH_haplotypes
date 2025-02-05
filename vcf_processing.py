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

vcf_in = sys.argv[1] #merged vcf file
var_table = pd.read_excel(sys.argv[2]) #preprocessed excel file with variant's annotation (should contain POS, REF and ALT columns) for matching

with open(sys.argv[3], 'rt') as f:
    avg_VAF = json.load(f)

with open(sys.argv[4], 'rt') as f:
    var_VAF = json.load(f)

diff_variants=pd.read_excel(sys.argv[5])
o_dir=str(sys.argv[6]) #output directory to store files

    
#preparing vcf file to further annotation
def vcf_to_table(vcf_file):
    #transform vcf format into df - skipping header info and renaming columns
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

    return cohort_AF(vcf_copy)

#count allele frequency for each variant
def cohort_AF(df):
    subset=df.iloc[:,5:].T.groupby(lambda x: x.split('_')[0]).sum().T.replace(2,1)
    subset['AF'] = subset.sum(axis=1)
    df['AF'] = subset['AF']/(len(subset.columns)-1)
    return df

#matching result df with variant table and extracting corresponding columns (rsid, phenotype, ACMG classification) as table_annotated.xlsx
def variant_annotation(table, variants):
    variants['POS'] = variants['POS'].fillna(0).astype(int)
    matched = pd.merge(table,
    variants[['POS', 'REF', 'ALT', 'Effect', 'rs#', 'Allele associated phenotype', 'ACMG']],  on=['POS', 'REF', 'ALT'],  how='left') 
    # matched.to_excel(os.path.join(o_dir+'table_annotated.xlsx'), index=False)
    return matched


#creating a separate excel file with added average VAF as VAF_annotated.xlsx
def VAF_table(df, avg, var):
    VAF=[]
    vafs=pd.DataFrame(var)
    result=pd.concat([df.iloc[:,0:4], vafs], axis=1)
    for col in result.columns:
        if col.startswith('770'):
            VAF.append(avg[col])
        else:
            VAF.append(np.nan)

    result.loc[len(result)]=VAF
    result.to_excel(os.path.join(o_dir+'cohort_VAFs.xlsx'), index=False)

#add chimeras assigning to the annotated file     
def chimeras(samples, pseudo): #function takes df with annotated variants and table with pVCF of differential positions that are found only in pseudogene 
    samples['pVCF']=samples['POS'].astype(str) + "-" + samples['REF'] + "-" + samples['ALT']
    results = {col: [] for col in samples.columns if col.startswith('770')}

    #assigning non-differential variants as 0 to make finding chimeric regions easier
    for index, row in samples.iterrows():
        if str(row['pVCF']) in pseudo['pVCF'].to_list():
            for col in results.keys():
                results[col].append(int(row[col])) if pd.isna(row[col]) == False else results[col].append(0)
        else:
            for col in results.keys():
                results[col].append(0)

    #if chromosome has 3 or more differential variants in a row it would be assigned as chimeric
    for key in results.keys():
        if np.any(np.convolve(results[key], np.ones(3, dtype=int), mode='valid') == 3):
            results[key] = 'chimera'
        else:
            results[key] =  np.nan

    new_row = {col: results.get(col, '') for col in samples.columns}
    samples.loc[len(samples)] = new_row
    samples.to_excel(os.path.join(o_dir+'annotated_variants.xlsx'), index=False)



phased = vcf_to_table(vcf_in)
annotated = variant_annotation(phased, var_table)
VAF_table(annotated, avg_VAF, var_VAF)
chimeras(annotated, diff_variants)