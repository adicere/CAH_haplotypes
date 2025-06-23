import pandas as pd
import numpy as np
from collections import defaultdict

def kinship_type(rel_i, rel_j):
    if rel_i == 'PARENT' and rel_j == 'PARENT':
        return 'spouse'
    if rel_i == 'PARENT' and rel_j in ('DIABETIC_ILL', 'PATIENT','SIBLING'):
        return 'child'
    if rel_i in ('DIABETIC_ILL','PATIENT','SIBLING') and rel_j == 'PARENT':
        return 'parent'
    if rel_i in ('DIABETIC_ILL', 'PATIENT','SIBLING') and rel_j in ('DIABETIC_ILL', 'PATIENT','SIBLING'):
        return 'sibling'
    return 'unknown'

def kinship_data(ann_file):
    kinship=defaultdict(dict)
    for fam, fam_df in ann_file.groupby('family'):
        records = fam_df.to_dict('records')
        for i in records:
            i_id, i_rel = i['patient'], i['aggr_relation']
            rel_dict = kinship[i_id]
            for j in records:
                j_id, j_rel = j['patient'], j['aggr_relation']
                if i_id == j_id:
                    continue
                rel_dict[j_id] = kinship_type(i_rel, j_rel)
    return kinship

def read_vcf(link, mode=None): 
    vcf = pd.read_csv(link, sep='\t', comment='#', header=None,
          names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], index_col=False) 
    vcf = vcf.assign(GT=vcf['SAMPLE'].map(lambda x: x.split(':')[0]))
    vcf['GT'] = vcf['GT'].str.replace('|', '/', regex=False)
    if mode == 'hc':
        vcf = vcf.assign(AD=vcf['SAMPLE'].map(lambda x: x.split(':')[1]))
        vcf['DP'] = vcf['AD'].map(lambda x: sum(map(int, x.split(','))) if pd.notna(x) else None)
        vcf['DP'] = vcf['DP'].astype(int)
    if mode == 'dv':
        vcf = vcf.assign(AD=vcf['SAMPLE'].map(lambda x: x.split(':')[3]))
        vcf['DP'] = vcf['AD'].map(lambda x: sum(map(int, x.split(','))) if pd.notna(x) else None)
        vcf['DP'] = vcf['DP'].astype(int)

    return vcf

def gt_vectors(vcf):
    gt_vector={}
    for idx, row in vcf.iterrows(): 
        if ',' in row['ALT']: 
            dp=int(row['AD'].split(',')[0])+int(row['AD'].split(',')[1]) + int(row['AD'].split(',')[2])
            if row['POS'] == 32038855:
                position = '32038855' + '_' + 'T' + "_" + 'TTTG'
                gt = 2
                vaf=(int(row['AD'].split(',')[1]) + int(row['AD'].split(',')[2]))/dp
                gt_vector[position] = [gt, vaf, dp]
            else:
                for num, alt in enumerate(row['ALT'].split(',')):
                    position = str(row['POS']) + '_' + row['REF'] + "_" + alt
                    gt = 1
                    vaf=int(row['AD'].split(',')[num+1])/dp
                    gt_vector[position] = [gt, vaf, row['DP']]
                    
        else: 
            position=str(row['POS']) + '_' + row['REF'] + "_" + row['ALT']
            gt=int(row['GT'].split('/')[0]) + int(row['GT'].split('/')[1])
            vaf=int(row['AD'].split(',')[1])/row['DP']
            gt_vector[position] = [gt, vaf, row['DP']]
    return gt_vector

def check_inheritance(patient, pos, kinship_data, gt_hc, gt_dv, gt_bam):
    closest_relatives=['child', 'sibling', 'parent']
    relatives = kinship_data.get(patient, {})
    if not relatives:
        return False
    inheritance=[]
    for rel, kinship_degree in relatives.items():
        if kinship_degree not in closest_relatives:
            continue
        rel_id = int(rel)

        hc_data = gt_hc.get(rel_id, {}).get(pos)
        dv_data = gt_dv.get(rel_id, {}).get(pos)
        bam_data = gt_bam.get(rel_id, {}).get(pos)

        if hc_data:
            genotype, vaf, dp = hc_data
            if genotype == 2 and vaf >= 0.85:
                msg = "homozygote" if dp >= 30 else "lowqual homozygote"
                print(f'found a {msg} {pos} in {rel_id} with VAF {vaf}, DP {dp}')
                inheritance.append(genotype)
            elif genotype == 1 and 0.2 <= vaf < 0.85:
                msg = "heterozygote" if dp >= 30 else "lowqual heterozygote"
                print(f'found a {msg} {pos} in {rel_id} with VAF {vaf}, DP {dp}')
                inheritance.append(genotype)
            else:
                print(f'found only low qual {pos} in {rel_id} with VAF {vaf}, DP {dp}, probably an artefact')

        elif dv_data:
            print(f'{pos} found only in DV data of {rel_id}')
            dv_genotype = dv_data[0]
            if bam_data:
                if bam_data[0] == dv_genotype:
                    print(f'DV variant in {rel_id} is supported by BAM')
                else:
                    print(f'DV genotype ({dv_genotype}) in {rel_id} differs from BAM data: {bam_data}')
                inheritance.append(bam_data[0])
            else:
                print(f'cannot support found DV variant in {rel_id} by BAM, probably an artefact')

        else:
            print(f'no {pos} found in {rel_id}')
        print(inheritance)
        return bool(inheritance)
    
def check_conflicts(ann_file, calling_region = None):
    kinship = kinship_data(ann_file)
    resolved = defaultdict(dict)
    gt_hc=defaultdict(dict)
    gt_dv=defaultdict(dict)
    gt_bam=defaultdict(dict)
    for idx, row in ann_file.iterrows():
        # dv
        vcf = read_vcf(row['dv'])
        vcf = vcf[vcf['FILTER'] == 'PASS']
        if calling_region:
            vcf = vcf[(vcf['POS'] >= calling_region[0]) & (vcf['POS'] <= calling_region[1])] 
        gt_dv[idx]=gt_vectors(vcf, mode = 'dv')
        # hc 
        vcf = read_vcf(row['hc'])
        if calling_region:
            vcf = vcf[(vcf['POS'] >= calling_region[0]) & (vcf['POS'] <= calling_region[1])]
        gt_hc[idx]=gt_vectors(vcf, mode= 'hc')
        #bam
        vcf = read_vcf(row['noBAQ'])
        if calling_region:
            vcf = vcf[(vcf['POS'] >= calling_region[0]) & (vcf['POS'] <= calling_region[1])]
        gt_bam[idx]=gt_vectors(vcf, mode= 'dv')
    
    for sample in ann_file.index:
        sample_hc = gt_hc.get(sample, {})
        sample_dv = gt_dv.get(sample, {})
        sample_bam = gt_bam.get(sample, {})

        print(f'checking variants for {idx}')
        for position in sample_hc.keys():
            hc, vaf, dp = sample_hc
            if position in sample_dv.keys():
                dv=sample_dv[position][0]
                if hc == dv:
                    print(f'no conflict for {position} in {idx} with {sample_hc[position]}')
                    resolved[idx][position] = 'hc'
                elif hc!=dv:
                    print(f'check {position}: HC - {sample_hc[position]}, DV - {dv}, bam - {sample_bam.get(position, 'Not found in bam')}')
                    if position in sample_bam.keys():
                            mp=sample_bam[position][0]
                            if mp == hc:
                                print(f'HC genotype is supported by bam for {position} in {idx};'
                                      f'HC - {sample_hc[position]}, DV - {dv}, bam - {sample_bam[position]}')
                                resolved[idx][position] = 'hc'
                            else:
                                print(f'DV genotype is supported by bam for {position} in {idx};'
                                      f'HC - {sample_hc[position]}, DV - {dv}, bam - {sample_bam[position]}')
                                resolved[idx][position] = 'dv'
                    else:
                        if (hc == 2 and vaf>0.85) or (hc==1 and vaf>= 0.2 and vaf <0.85):
                            if dp>=30:
                                print(f'choice in favor of HC genotype for {position} in {idx};'
                                      f'HC - {sample_hc[position]}, DV - {dv}')
                                resolved[idx][position] = 'hc'
                            else:
                                print(f'data is unclear; need to check inheritance for {position} in {idx}')
                                inheritance = check_inheritance(idx, position, kinship, gt_hc, gt_dv, gt_bam)
                                if inheritance: 
                                    print('saving this variant even with low DP')
                                    resolved[idx][position] = 'hc'
                        else:
                            print(f'data is unclear; need to check inheritance for {position} in {idx}')
                            inheritance = check_inheritance(idx, position, kinship, gt_hc, gt_dv, gt_bam)
                            if inheritance:
                                print('position found in relatives, but genotype needs manual inspection')
                                # ?
                            
                                

                    
annotation=pd.read_csv('/home/rutkovskaya.ea/haplotypes/text_files/amplicon_annotation_processed.tsv', sep='\t', index_col=0)
filtered = annotation.loc[~(annotation['family'].map(lambda x: int(x) < 0 if 'diab' not in x else False))]
files_list= pd.read_excel('/home/rutkovskaya.ea/haplotypes/text_files/callers_comparison_wo_singles.xlsx', index_col=0)
kinship = kinship_data(filtered)
calling_region = [32037643, 32041345]
