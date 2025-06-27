import pandas as pd
import numpy as np
from collections import defaultdict
import sys 
import os
import logging

logging.basicConfig(
    filename='check_results.log',
    filemode='w', 
    format='%(message)s',  
    level=logging.INFO
)

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

def gt_vectors(vcf, mode = None):
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

# def check_inheritance(patient, pos, kinship_data, gt_hc, gt_dv, gt_bam):
#     closest_relatives=['child', 'sibling', 'parent']
#     relatives = kinship_data.get(patient, {})
#     if not relatives:
#         return False
#     inheritance=[]
#     for rel, kinship_degree in relatives.items():
#         if kinship_degree not in closest_relatives:
#             continue
#         rel_id = str(rel)

#         hc_data = gt_hc.get(rel_id, {}).get(pos)
#         dv_data = gt_dv.get(rel_id, {}).get(pos)
#         bam_data = gt_bam.get(rel_id, {}).get(pos)

#         if hc_data:
#             genotype, vaf, dp = hc_data
#             if genotype == 2 and vaf >= 0.85:
#                 msg = "homozygote" if dp >= 30 else "lowqual homozygote"
#                 logging.info(f'found a {msg} {pos} in {rel_id} with VAF {vaf}, DP {dp}')
#                 inheritance.append(genotype)
#             elif genotype == 1 and 0.2 <= vaf < 0.85:
#                 msg = "heterozygote" if dp >= 30 else "lowqual heterozygote"
#                 logging.info(f'found a {msg} {pos} in {rel_id} with VAF {vaf}, DP {dp}')
#                 inheritance.append(genotype)
#             else:
#                 logging.info(f'found only low qual {pos} in {rel_id} with VAF {vaf}, DP {dp}, probably an artefact')

#         elif dv_data:
#             logging.info(f'{pos} found only in DV data of {rel_id}')
#             dv_genotype = dv_data[0]
#             if bam_data:
#                 if bam_data[0] == dv_genotype:
#                     logging.info(f'DV variant in {rel_id} is supported by BAM')
#                 else:
#                     logging.info(f'DV genotype ({dv_genotype}) in {rel_id} differs from BAM data: {bam_data}')
#                 inheritance.append(bam_data[0])
#             else:
#                 logging.info(f'cannot support found DV variant in {rel_id} by BAM, needs manual inspection')
#                 inheritance.append(dv_genotype)

#         elif bam_data:
#             logging.info(f'{pos} found only in BAM data of {rel_id}')
#             inheritance.append(bam_data[0])

#         else:
#             logging.info(f'no {pos} found in {rel_id}')
#         logging.info(inheritance)
#         return bool(inheritance)

def check_inheritance(patient, pos, kinship_data, gt_hc, gt_dv, gt_bam, ill_patients):
    closest_relatives=['child', 'sibling', 'parent', 'spouse']
    relatives = kinship_data.get(patient, {})
    if not relatives:
        return False
    inheritance=[]
    samples = [str(patient)] + [
        str(rel) for rel, degree in relatives.items()
        if degree in closest_relatives
    ]
    # print(samples)
    matrix = pd.DataFrame(index=samples, columns=['HC', 'DV', 'BAM'])

    for sample in samples:
        if sample in ill_patients:
            ill = sample
        sample_str = str(sample)
        matrix.loc[sample, 'HC'] = int(bool(gt_hc.get(sample_str, {}).get(pos)))
        matrix.loc[sample, 'DV'] = int(bool(gt_dv.get(sample_str, {}).get(pos)))
        matrix.loc[sample, 'BAM'] = int(bool(gt_bam.get(sample_str, {}).get(pos)))
    
    matrix['total'] = matrix.sum(axis=1)
    res = (matrix['total'] > 0).sum()
    # print(matrix)
    # print('*'*10)
    if res > 1:
        logging.info(f'more than one person in a family have {pos}, saving this variant for {patient}')
        if matrix.loc[ill, 'total'] == 0:
            logging.info(f'ill patient {ill} does not have parents variant, need manual inspection')
        return True
    else:
        log = (f'relatives do not have {pos}, patients data:',
                     f'HC - {gt_hc.get(str(patient), {}).get(pos)}',
                     f'DV - {gt_dv.get(str(patient), {}).get(pos)}',
                     f'BAM - {gt_bam.get(str(patient), {}).get(pos)}')
        logging.info(log)

        return False    

    
def check_conflicts(ann_file, kinship_data, ill_samples, calling_region = None):
    # kinship = kinship_data(ann_file)
    resolved = defaultdict(dict)
    gt_hc=defaultdict(dict)
    gt_dv=defaultdict(dict)
    gt_bam=defaultdict(dict)
    for idx, row in ann_file.iterrows():
        # dv
        vcf = read_vcf(row['dv'], mode= 'dv')
        vcf = vcf[vcf['FILTER'] == 'PASS']
        if calling_region:
            vcf = vcf[(vcf['POS'] >= calling_region[0]) & (vcf['POS'] <= calling_region[1])] 
        gt_dv[str(idx)]=gt_vectors(vcf)
        # hc 
        vcf = read_vcf(row['hc'], mode = 'hc')
        if calling_region:
            vcf = vcf[(vcf['POS'] >= calling_region[0]) & (vcf['POS'] <= calling_region[1])]
        gt_hc[str(idx)]=gt_vectors(vcf)
        #bam
        vcf = read_vcf(row['noBAQ'], mode = 'dv')
        if calling_region:
            vcf = vcf[(vcf['POS'] >= calling_region[0]) & (vcf['POS'] <= calling_region[1])]
        gt_bam[str(idx)]=gt_vectors(vcf)

    for sample in ann_file.index:
        sample_hc = gt_hc.get(str(sample), {})
        sample_dv = gt_dv.get(str(sample), {})
        sample_bam = gt_bam.get(str(sample), {})

        logging.info(f'checking variants for {sample}')
        for position in sample_hc.keys():
            hc, vaf, dp = sample_hc[position]
            if position in sample_dv.keys():
                dv=sample_dv[position][0]
                if hc == dv:
                    # logging.info(f'no conflict for {position} in {sample} with {sample_hc[position]}')
                    resolved[sample][position] = 'hc'
                elif hc!=dv:
                    logging.info(f'check {position}: HC - {sample_hc[position]}, DV - {dv}, bam - {sample_bam.get(position, 'Not found in bam')}')
                    if position in sample_bam.keys():
                            mp=sample_bam[position][0]
                            if mp == hc:
                                logging.info(f'HC genotype is supported by bam for {position} in {sample};'
                                      f'HC - {sample_hc[position]}, DV - {dv}, bam - {sample_bam[position]}')
                                resolved[sample][position] = 'hc'
                            else:
                                logging.info(f'DV genotype is supported by bam for {position} in {sample};'
                                      f'HC - {sample_hc[position]}, DV - {dv}, bam - {sample_bam[position]}')
                                resolved[sample][position] = 'dv'
                    else:
                        if (hc == 2 and vaf>0.85) or (hc==1 and vaf>= 0.2 and vaf <0.85):
                            if dp>=30:
                                logging.info(f'choice in favor of HC genotype for {position} in {sample};'
                                      f'HC - {sample_hc[position]}, DV - {dv}')
                                resolved[sample][position] = 'hc'
                            else:
                                logging.info(f'data is unclear; need to check inheritance for {position} in {sample}')
                                inheritance = check_inheritance(sample, position, kinship_data, gt_hc, gt_dv, gt_bam, ill_samples)
                                if inheritance: 
                                    logging.info('saving this variant even with low DP')
                                    resolved[sample][position] = 'hc'
                                else:
                                    logging.info('droping this variant, probably artefact')
                        else:
                            logging.info(f'data is unclear; need to check inheritance for {position} in {sample}')
                            inheritance = check_inheritance(sample, position, kinship_data, gt_hc, gt_dv, gt_bam, ill_samples)
                            if inheritance:
                                logging.info('position found in relatives, but genotype needs manual inspection')
                                # ?
                            else:
                                logging.info('droping this variant, probably artefact')
            else:
                logging.info(f'only called by HC for {position} in {sample}')
                if ((hc == 2 and vaf>0.85) or (hc==1 and vaf>= 0.2 and vaf <0.85)) and dp>=30:
                    logging.info(f'{position} in {sample} is OK')
                    resolved[sample][position] = 'hc'        
                else:
                    logging.info(f'unclear {position} varinant in {sample} with {sample_hc[position]}')
                    if position in sample_bam.keys():
                        mp=sample_bam[position][0]
                        if mp == hc:
                            logging.info(f'HC genotype is supported by bam for {position} in {sample}')
                            resolved[sample][position] = 'hc'
                        else:
                            logging.info(f'HC is conflicting with naive calling for {position} in {sample};' 
                                         f'HC - {sample_hc[position]}, BAM - {sample_bam[position]}')
                            resolved[sample][position] = 'mp'
                    else:
                        logging.info(f'{position} in {sample} with {sample_hc[position]} cannot be supported by BAM, checking inheritance')
                        inheritance = check_inheritance(sample, position, kinship_data, gt_hc, gt_dv, gt_bam, ill_samples)
                        if inheritance: 
                            logging.info('saving this variant, but genotype needs manual confirmation')
                            resolved[sample][position] = 'hc'
                        else:
                            logging.info('droping this variant, probably artefact')

        for position, dv_data in sample_dv.items():
            if position not in resolved[sample]:
                logging.info(f'{position} in {sample} found in DV, but not in HC, checking it')
                dv=dv_data[0]  
                bam = sample_bam.get(position)
                if bam:
                    mp=bam[0]
                    if dv == mp:
                        logging.info(f'DV {position} in {sample} is supported by BAM data {sample_bam[position]}')
                        resolved[sample][position] = 'dv'
                    else:
                        logging.info(f'DV and mpileup for {position} in {sample} are conflicted, need to check manually')
                else:
                    logging.info(f'DV {position} in {sample} cannot be supported by BAM data, need to check inheritance')
                    inheritance = check_inheritance(sample, position, kinship_data, gt_hc, gt_dv, gt_bam, ill_samples)
                    if inheritance == True:
                        logging.info('saving this variant due to inheritance')
                        resolved[sample][position] = 'dv'
                    else:
                        logging.info('droping this variant, probably artefact')         
        for position in sample_bam:
            if position not in resolved[sample]:
                logging.info(f'{position} in {sample} was called only by naive calling with {sample_bam[position][0]} genotype')
                resolved[sample][position] = 'mp'
       
    return resolved
                        
                            
                    
annotation=pd.read_csv('/home/rutkovskaya.ea/haplotypes/text_files/dev_files/amplicon_annotation_processed.tsv', sep='\t', index_col=0)
filtered = annotation[~annotation['patient'].isna()]
# filtered = annotation.loc[~(annotation['family'].map(lambda x: int(x) < 0 if 'diab' not in x else False))]
filtered['patient'] = filtered['patient'].astype(int)
files_list= pd.read_excel('/home/rutkovskaya.ea/haplotypes/text_files/callers_comparison.xlsx', index_col=0)
kinship = kinship_data(filtered)
filtered['patient'] = filtered['patient'].astype(str)
ill_patients = filtered.loc[filtered['group'] == 'ILL', 'patient'].tolist()
calling_region = [32037643, 32041345]
res = check_conflicts(files_list, kinship, ill_patients, calling_region)
