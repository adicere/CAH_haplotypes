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


def check_inheritance(patient, pos, kinship_data, gt_hc, gt_dv, gt_bam, ill_patients):
    closest_relatives=['child', 'sibling', 'parent', 'spouse']
    relatives = kinship_data.get(patient, {})
    if not relatives:
        return False, 'No data about relatives found for patient.'
    samples = [str(patient)] + [
        str(rel) for rel, degree in relatives.items()
        if degree in closest_relatives
    ]
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

    if res > 1:
        if matrix.loc[ill, 'total'] == 0:
            msg = '\n'.join([f"[INHERITANCE] {pos} was confirmed by family analysis - at least one more person has variant, saving this position, but genotype needs manual inspection",
                             f"[WARNING] Ill patient {ill} does not have parents variant {pos}, needs manual inspection"])
        else:
            msg = f'[INHERITANCE] {pos} was confirmed by family analysis - at least one more person has variant, saving this position, but genotype needs manual inspection'
    
        return True, msg
    else:
        msg = ', '.join([f"[INHERITANCE] Relatives do not have {pos}",
                     f"patient's data: HC - {gt_hc.get(str(patient), {}).get(pos)}",
                     f"DV - {gt_dv.get(str(patient), {}).get(pos)}",
                     f"BAM - {gt_bam.get(str(patient), {}).get(pos)}."])

        return False, msg  

def for_manual_inspection(sample, position, hc_data, dv_data, bam_data):
    bam_pos=int(position.split('_')[0]) - 32037619
    hc_raw = hc_data.get(sample, {}).get(position)
    dv_raw = dv_data.get(sample, {}).get(position)
    bam_raw = bam_data.get(sample, {}).get(position)

    hc = [str(x) for x in hc_raw] if hc_raw else ['N/A']
    dv = [str(x) for x in dv_raw] if dv_raw else ['N/A']
    bam = [str(x) for x in bam_raw] if bam_raw else ['N/A']

    sus_variant={
        'Sample': sample,
        'Variant': position,
        'BAM_position': bam_pos,
        'HC_call': '; '.join(hc), 
        'DV_call': '; '.join(dv), 
        'mpileup_call': '; '.join(bam)
    }
    
    return sus_variant


def check_conflicts(ann_file, kinship_data, ill_samples, calling_region = None):
    # kinship = kinship_data(ann_file)
    resolved = defaultdict(dict)
    gt_hc=defaultdict(dict)
    gt_dv=defaultdict(dict)
    gt_bam=defaultdict(dict)
    manual_check=[]
    all_logs = {}
    low_qual = defaultdict(lambda: {'Variant': [], 'Parameters': []})

    for idx, row in ann_file.iterrows():
        # dv
        vcf = read_vcf(row['dv_modified'], mode= 'dv')
        vcf = vcf[vcf['FILTER'] == 'PASS']
        if calling_region:
            vcf = vcf[(vcf['POS'] >= calling_region[0]) & (vcf['POS'] <= calling_region[1])] 
        gt_dv[str(idx)]=gt_vectors(vcf)
        # hc 
        vcf = read_vcf(row['hc_modified'], mode = 'hc')
        if calling_region:
            vcf = vcf[(vcf['POS'] >= calling_region[0]) & (vcf['POS'] <= calling_region[1])]
        gt_hc[str(idx)]=gt_vectors(vcf)
        #bam
        vcf = read_vcf(row['mpileup_modified'], mode = 'dv')
        if calling_region:
            vcf = vcf[(vcf['POS'] >= calling_region[0]) & (vcf['POS'] <= calling_region[1])]
        gt_bam[str(idx)]=gt_vectors(vcf)

    for sample in ann_file.index:
        sample_hc = gt_hc.get(str(sample), {})
        sample_dv = gt_dv.get(str(sample), {})
        sample_bam = gt_bam.get(str(sample), {})

        sample_checked = False
        log_buffer = []

        for position in sample_hc.keys():
            hc, vaf, dp = sample_hc[position]
            if position in sample_dv.keys():
                dv=sample_dv[position][0]
                if hc == dv:
                    resolved[sample][position] = 'hc'
                elif hc!=dv:
                    # log_buffer.append(f"check {position}: HC - {sample_hc[position]}, DV - {dv}, bam - {sample_bam.get(position, 'Not found in bam')}")
                    if position in sample_bam.keys():
                            mp=sample_bam[position][0]
                            if mp == hc:
                                sample_checked=True
                                log_buffer.append(f"[CONFLICT] Conflict at {position}; HC genotype was supported by mpileup data")
                                resolved[sample][position] = 'hc'
                            else:
                                sample_checked=True
                                log_buffer.append(f"[CONFLICT] Conflict at {position}; DV genotype was supported by mpileup data")
                                resolved[sample][position] = 'dv'
                    else:
                        if (hc == 2 and vaf>0.85) or (hc==1 and vaf>= 0.2 and vaf <0.85):
                            if dp>=30:
                                sample_checked=True
                                log_buffer.append(f"[CONFLICT] Conflict at {position};"
                                                  f"no mpileup data, but HC genotype is accepted due to correct VAF value")
                                resolved[sample][position] = 'hc'
                            else:
                                low_qual[str(sample)]['Variant'].append(position)
                                low_qual[str(sample)]['Parameters'].append(', '.join(str(p) for p in sample_hc[position]))
                                # log_buffer.append(f'data is unclear; need to check inheritance for {position} in {sample}')
                                inheritance, message = check_inheritance(sample, position, kinship_data, gt_hc, gt_dv, gt_bam, ill_samples)
                                if inheritance: 
                                    sample_checked=True
                                    log_buffer.append(message)
                                    resolved[sample][position] = 'hc'
                                else:
                                    sample_checked=True
                                    if message == 'No data about relatives found for patient.':
                                        manual_check.append(for_manual_inspection(str(sample), position, gt_hc, gt_dv, gt_bam))
                                        log_buffer.append(f"[INHERITANCE] Patient doesn't have relatives. {position} needs manual inspection.")
                                    else:
                                        log_buffer.append(message + ' Dropping this variant, probably an artefact')
                        else:
                            # log_buffer.append(f'data is unclear; need to check inheritance for {position} in {sample}')
                            inheritance, message = check_inheritance(sample, position, kinship_data, gt_hc, gt_dv, gt_bam, ill_samples)
                            if inheritance:
                                sample_checked=True
                                log_buffer.append(message)
                                resolved[sample][position] = 'hc'
                                manual_check.append(for_manual_inspection(str(sample), position, gt_hc, gt_dv, gt_bam))
                            else:
                                sample_checked=True
                                if message == 'No data about relatives found for patient.':
                                    manual_check.append(for_manual_inspection(str(sample), position, gt_hc, gt_dv, gt_bam))
                                    log_buffer.append(f"[INHERITANCE] Patient doesn't have relatives. {position} needs manual inspection.")
                                else:
                                    log_buffer.append(message + ' Dropping this variant, probably an artefact')
            else:
                # log_buffer.append(f'only called by HC for {position} in {sample}')
                if ((hc == 2 and vaf>0.85) or (hc==1 and vaf>= 0.2 and vaf <0.85)) and dp>=30:
                    # log_buffer.append(f'{position} in {sample} is OK')
                    resolved[sample][position] = 'hc'        
                else:
                    # log_buffer.append(f'unclear {position} variant in {sample} with {sample_hc[position]}')
                    if position in sample_bam.keys():
                        mp=sample_bam[position][0]
                        if mp == hc:
                            sample_checked=True
                            log_buffer.append(f'[CONFLICT] Conflict at {position}; HC genotype was supported by mpileup data')
                            resolved[sample][position] = 'hc'
                        else:
                            sample_checked=True
                            log_buffer.append(f'[CONFLICT] Conflict at {position}; Accepting mpileup genotype instead of HC' 
                                         f'HC - {sample_hc[position]}, BAM - {sample_bam[position]}')
                            resolved[sample][position] = 'mp'
                    else:
                        # log_buffer.append(f'{position} in {sample} with {sample_hc[position]} cannot be supported by BAM, checking inheritance')
                        inheritance, message = check_inheritance(sample, position, kinship_data, gt_hc, gt_dv, gt_bam, ill_samples)
                        if inheritance: 
                            sample_checked=True
                            log_buffer.append(message)
                            resolved[sample][position] = 'hc'
                            manual_check.append(for_manual_inspection(str(sample), position, gt_hc, gt_dv, gt_bam))
                        else:
                            sample_checked=True
                            if message == 'No data about relatives found for patient.':
                                manual_check.append(for_manual_inspection(str(sample), position, gt_hc, gt_dv, gt_bam))
                                log_buffer.append(f"[INHERITANCE] Patient doesn't have relatives. {position} needs manual inspection.")
                            else:
                                log_buffer.append(message + ' Dropping this variant, probably an artefact')

        for position, dv_data in sample_dv.items():
            if position not in resolved[sample]:
                # log_buffer.append(f'{position} in {sample} found in DV, but not in HC, checking it')
                dv=dv_data[0]  
                bam = sample_bam.get(position)
                if bam:
                    mp=bam[0]
                    if dv == mp:
                        # log_buffer.append(f'DV {position} in {sample} is supported by BAM data {sample_bam[position]}')
                        resolved[sample][position] = 'dv'
                    else:
                        sample_checked=True
                        log_buffer.append(f'[CONFLICT] Conflict at {position}; Accepting mpileup genotype instead of DV')
                        resolved[sample][position] = 'mp'
                else:
                    # log_buffer.append(f'DV {position} in {sample} cannot be supported by BAM data, need to check inheritance')
                    inheritance, message = check_inheritance(sample, position, kinship_data, gt_hc, gt_dv, gt_bam, ill_samples)
                    if inheritance == True:
                        sample_checked=True
                        log_buffer.append(message)
                        resolved[sample][position] = 'dv'
                        manual_check.append(for_manual_inspection(str(sample), position, gt_hc, gt_dv, gt_bam))
                    else:
                        sample_checked=True
                        if message == 'No data about relatives found for patient.':
                            manual_check.append(for_manual_inspection(str(sample), position, gt_hc, gt_dv, gt_bam))
                            log_buffer.append(f"[INHERITANCE] Patient doesn't have relatives. {position} needs manual inspection.")
                        else:
                            log_buffer.append(message + ' Dropping this variant, probably an artefact')         
        for position in sample_bam:
            if position not in resolved[sample]:
                mp, vaf, dp = sample_bam[position]
                if dp >= 30:
                    gt = 'heterozygous' if sample_bam[position][0] ==1 else 'homozygous'
                    sample_checked=True
                    log_buffer.append(f'[NAIVE VARIANT] {position} was called only by naive calling - {sample_bam[position]}')
                    resolved[sample][position] = 'mp'
                else:
                    low_qual[str(sample)]['Variant'].append(position)
                    low_qual[str(sample)]['Parameters'].append(', '.join(str(p) for p in sample_bam[position]))

        if sample_checked:
            all_logs[str(sample)] = log_buffer

    samples_with_warning = [s for s, logs in all_logs.items() if any("[WARNING]" in line for line in logs)]
    samples_without_warning = [s for s in all_logs if s not in samples_with_warning]

    for sample in samples_with_warning + samples_without_warning:
        logging.info(f"==== Results of variant check for sample {sample} ====")
        for line in all_logs[sample]:
            logging.info(line)
        logging.info("")

    return resolved, manual_check, low_qual
                        
                            
                    
annotation=pd.read_excel('/home/rutkovskaya.ea/haplotypes/text_files/samples_for_analysis.xlsx', index_col=0)
filtered = annotation[~annotation['patient'].isna()]
filtered['patient'] = filtered['patient'].astype(int)
kinship = kinship_data(filtered)
filtered['patient'] = filtered['patient'].astype(str)
ill_patients = filtered.loc[filtered['group'] == 'ILL', 'patient'].tolist()
calling_region = [32037643, 32041345]
res_var, manual_insp, lq_variants = check_conflicts(filtered, kinship, ill_patients, calling_region)

# pd.DataFrame(res_var).to_excel('/home/rutkovskaya.ea/haplotypes/text_files/resolved.xlsx')


# rows = []
# for sample, data in lq_variants.items():
#     for variant, param in zip(data['Variant'], data['Parameters']):
#         rows.append({'Sample': sample, 'Variant': variant, 'Parameters': param})

# low_qual_df = pd.DataFrame(rows)


# pd.DataFrame(manual_insp).to_excel('/home/rutkovskaya.ea/haplotypes/text_files/manual_inspection.xlsx', index = False)
# low_qual_df.to_excel('/home/rutkovskaya.ea/haplotypes/text_files/lq_variants.xlsx', index=False)