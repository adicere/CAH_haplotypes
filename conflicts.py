import pandas as pd 
from collections import defaultdict
import numpy as np

#collecting data about families
trio_info=pd.read_excel('/storage/haplotypes/excels/trio_info.xlsx')
trio_info['relatives'] = trio_info['relatives'].replace(np.nan, 'None', regex=False)
trio_info['relatives']=trio_info['relatives'].astype(str)

trios=defaultdict(dict)
for idx, row in trio_info.iterrows():
    if row['relatives']:
        relatives=[x for x in row['relatives'].split('; ')]
        for i in relatives:
            person=i.split(':')
            if len(person)>1:
                trios[row['sample']][person[1]] = int(person[0])

files_list= pd.read_excel('/storage/haplotypes/file_paths.xlsx', index_col=0)
calling_region = [32037643, 32041345]

def read_vcf(link, mode=None): 
    vcf = pd.read_csv(link, sep='\t', comment='#', header=None,
          names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], index_col=False) 
    vcf = vcf.assign(GT=vcf['SAMPLE'].map(lambda x: x.split(':')[0]))
    if mode == 'hc':
        vcf = vcf.assign(AD=vcf['SAMPLE'].map(lambda x: x.split(':')[1]))
        vcf['DP'] = vcf['AD'].map(lambda x: sum(map(int, x.split(','))) if pd.notna(x) else None)
        vcf['DP'] = vcf['DP'].astype(int)
        vcf['GT'] = vcf['GT'].str.replace('/', '|', regex=False)
    if mode == 'dv':
        vcf = vcf.assign(AD=vcf['SAMPLE'].map(lambda x: x.split(':')[3]))
        vcf['DP'] = vcf['AD'].map(lambda x: sum(map(int, x.split(','))) if pd.notna(x) else None)
        vcf['GT'] = vcf['GT'].str.replace('|', '/', regex=False)
        vcf['DP'] = vcf['DP'].astype(int)
    if mode == 'bam':
        vcf = vcf.assign(AD=vcf['SAMPLE'].map(lambda x: x.split(':')[5]))
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
                if mode == 'hc':
                    gt = '1|1'
                else:
                    gt='1/1'
                vaf=(int(row['AD'].split(',')[1]) + int(row['AD'].split(',')[2]))/dp
                gt_vector[position] = [gt, vaf, dp]
            else:
                for num, alt in enumerate(row['ALT'].split(',')):
                    position = str(row['POS']) + '_' + row['REF'] + "_" + alt
                    if mode == 'hc':
                        gt = '1|0'
                    else:
                        gt='1/0'
                    vaf=int(row['AD'].split(',')[num+1])/dp
                    gt_vector[position] = [gt, vaf, row['DP']]
                    
        else: 
            position=str(row['POS']) + '_' + row['REF'] + "_" + row['ALT']
            gt=row['GT']
            vaf=int(row['AD'].split(',')[1])/row['DP']
            gt_vector[position] = [gt, vaf, row['DP']]
    return gt_vector

def check_inheritance(patient, pos, trios, gt_hc, gt_dv, gt_bam=None):
    found_in_both_parents=0
    if trios[patient] != {}:
        print(trios[patient])
    else:
        print(f'{patient} does not have relatives')
    inheritance=[]
    for rel in trios[patient].keys():
        if rel == 'child' or rel == 'mother' or rel == 'father':
            rel_num=int(trios[patient][rel])
            try:
                print(f'found {pos} in {rel_num} HC data - {gt_hc[rel_num][pos]}')
            except KeyError:
                print(f'{pos} not found in gt_hc of {rel_num} - the {rel} of {patient}')
                try:
                    print(f'found {pos} in {rel_num} DV data - {gt_dv[rel_num][pos]}')
                except KeyError:
                    print(f'{pos} not found in gt_dv of {rel_num} - the {rel} of {patient}')

            if pos in gt_hc[rel_num]:
                print(f'{pos} found in HC data of {rel}')
                hc_genotype=gt_hc[rel_num][pos][0]
                hc=hc_genotype.split('|')
                hc=int(hc[0])+int(hc[1])
                if hc == 2 and gt_hc[rel_num][pos][1]>=0.85:
                    if gt_hc[rel_num][pos][2]>=30:
                        print(f'found a homozygote {pos} in {rel} with VAF {gt_hc[rel_num][pos][1]}')
                        inheritance.append(hc_genotype)
                        if rel == 'mother' or rel == 'father':
                            found_in_both_parents+=1
                    else:
                        print(f'found a lowqual homozygote {pos} in {rel} with right VAF {gt_hc[rel_num][pos][1]}, but low DP - {gt_hc[rel_num][pos][2]}')
                        inheritance.append(hc_genotype)
                        if rel == 'mother' or rel == 'father':
                            found_in_both_parents+=1
                if hc == 1 and (gt_hc[rel_num][pos][1]>=0.2 and gt_hc[rel_num][pos][1]<0.85):
                    if gt_hc[rel_num][pos][2]>=30:
                        print(f'found a heterozygote {pos} in {rel} with VAF {gt_hc[rel_num][pos][1]}')
                        inheritance.append(hc_genotype)
                        if rel == 'mother' or rel == 'father':
                            found_in_both_parents+=1
                    else:
                        print(f'found a lowqual heterozygote {pos} in {rel} with right VAF {gt_hc[rel_num][pos][1]}, but low DP - {gt_hc[rel_num][pos][2]}')
                        inheritance.append(hc_genotype)
                        if rel == 'mother' or rel == 'father':
                            found_in_both_parents+=1
                else:
                    print(f'found only low qual {pos} in {rel} with {gt_hc[rel_num][pos][1]} and {gt_hc[rel_num][pos][2]}, probably an artefact')
            if pos in gt_dv[rel_num] and pos not in gt_hc[rel_num]:
                print(f'{pos} found only in DV data of {rel}')
                dv_genotype=gt_dv[rel_num][pos][0]
                dv=dv_genotype.split('/')
                dv=int(dv[0])+int(dv[1])
                if pos in gt_bam[rel_num].keys():
                    pp_value=gt_bam[rel_num][pos][0]
                    pp=pp_value.split('/')
                    pp=int(pp[0])+int(pp[1])
                    if pp==2 and gt_bam[rel_num][pos][1]>=0.85:
                        print(f'found a homozygote DV variant in {rel_num}')
                        inheritance.append(pp_value)
                        if rel == 'mother' or rel == 'father':
                            found_in_both_parents+=1
                    if pp==1 and (gt_bam[rel_num][pos][1]>=0.2 and gt_bam[rel_num][pos][1]<0.85):
                        print(f'found a heterozygote DV variant in {rel_num}')
                        inheritance.append(pp_value)
                        if rel == 'mother' or rel == 'father':
                            found_in_both_parents+=1
                # if dv == 2 and gt_dv[rel_num][pos][1]>=0.85:
                #     if gt_dv[rel_num][pos][2]>=30:
                #         print(f'found a homozygote {pos} in {rel} with VAF {gt_dv[rel_num][pos][1]}')
                #         inheritance.append(dv_genotype)
                #     else:
                #         print(f'found a lowqual homozygote {pos} in {rel} with right VAF {gt_dv[rel_num][pos][1]}, but low DP - {gt_dv[rel_num][pos][2]}')
                #         inheritance.append(dv_genotype)
                # if dv == 1 and (gt_dv[rel_num][pos][1]>=0.2 and gt_dv[rel_num][pos][1]<0.85):
                #     if gt_dv[rel_num][pos][2]>=30:
                #         print(f'found a heterozygote {pos} in {rel} with VAF {gt_dv[rel_num][pos][1]}')
                #         inheritance.append(dv_genotype)
                #     else:
                #         print(f'found a lowqual heterozygote {pos} in {rel} with right VAF {gt_dv[rel_num][pos][1]}, but low DP - {gt_dv[rel_num][pos][2]}')
                #         inheritance.append(dv_genotype)
                else:
                    print(f'cannot support found DV variant by bam, probably an artefact')

            if pos not in gt_hc[rel_num] and pos not in gt_dv[rel_num]:
                print(f'no {pos} found in {rel}')
    
    print(inheritance)
    if len(inheritance)>0:
        return True
    else:
        return False
    
gt_hc=defaultdict(dict)
gt_dv=defaultdict(dict)
gt_bam=defaultdict(dict)
for idx, row in files_list.iterrows():
    # dv
    vcf = read_vcf(row['dv'], mode = 'dv')
    vcf = vcf[vcf['FILTER'] == 'PASS']
    vcf = vcf[(vcf['POS'] >= calling_region[0]) & (vcf['POS'] <= calling_region[1])]
    gt_dv[idx]=gt_vectors(vcf, mode = 'dv')
    # hc 
    vcf = read_vcf(row['hc'], mode='hc')
    vcf = vcf[(vcf['POS'] >= calling_region[0]) & (vcf['POS'] <= calling_region[1])]
    gt_hc[idx]=gt_vectors(vcf, mode= 'hc')
    #bam
    vcf = read_vcf(row['pileup'], mode='bam')
    vcf = vcf[(vcf['POS'] >= calling_region[0]) & (vcf['POS'] <= calling_region[1])]
    gt_bam[idx]=gt_vectors(vcf, mode= 'dv')

res = defaultdict(dict)
for idx, row in files_list.iterrows():
    # dv
    vcf = read_vcf(row['dv'], mode = 'dv')
    vcf = vcf[vcf['FILTER'] == 'PASS']
    vcf = vcf[(vcf['POS'] >= calling_region[0]) & (vcf['POS'] <= calling_region[1])]
    sample_dv=gt_vectors(vcf, mode = 'dv')
    # hc 
    vcf = read_vcf(row['hc'], mode='hc')
    vcf = vcf[(vcf['POS'] >= calling_region[0]) & (vcf['POS'] <= calling_region[1])]
    sample_hc=gt_vectors(vcf, mode= 'hc')

    #pileup
    vcf = read_vcf(row['pileup'], mode='bam')
    vcf = vcf[(vcf['POS'] >= calling_region[0]) & (vcf['POS'] <= calling_region[1])]
    sample_pp=gt_vectors(vcf, mode= 'dv')

    print(f'checking variants for {idx}')
    for position in sample_hc.keys():
        hc_value=sample_hc[position][0]
        hc=hc_value.split('|')
        hc=int(hc[0])+int(hc[1])
        if position in sample_dv.keys():
            dv_value=sample_dv[position][0]
            dv=dv_value.split('/')
            dv=int(dv[0])+int(dv[1])
            if hc == dv:
                # print(f'no conflict for {position} in {idx}, saving HC genotype')
                res[idx][position] = hc_value
            if hc!=dv:
                if hc == 2:
                    if sample_hc[position][1] > 0.85 and sample_hc[position][2]>=30:
                        res[idx][position] = hc_value
                    else:
                        print(f'check {position}: HC - {sample_hc[position]}, DV - {dv_value}, bam - {sample_pp.get(position, 'Not found in bam')}')
                        if position in sample_pp.keys():
                            pp_value=sample_pp[position][0]
                            pp=pp_value.split('/')
                            pp=int(pp[0])+int(pp[1])
                            if pp==2 and sample_pp[position][1]>=0.85:
                                res[idx][position]  = hc_value
                                print(f'HC genotype is supported by bam')
                            if pp==1 and (sample_pp[position][1]>=0.2 and sample_pp[position][1]<0.85):
                                print(f'HC genotype is not supported by bam, saving a new genotype')
                                print(sample_pp[position])
                                res[idx][position] = pp_value
                        else:
                            print('need to check inheritance')
                            inheritance = check_inheritance(idx, position, trios, gt_hc, gt_dv, gt_bam)
                            if inheritance == True:
                                print('saving this variant even with bad QUAL')
                                res[idx][position]  = hc_value
                            else:
                                print('droping this variant, probably artefact')
                                
                if hc == 1: 
                    if (sample_hc[position][1] >= 0.2 and sample_hc[position][1] < 0.85) and sample_hc[position][2]>=30:
                        res[idx][position]  = hc_value
                    else:
                        print(f'check {position}: HC - {sample_hc[position]}, DV - {dv_value}, {sample_pp.get(position, 'Not found in bam')}')
                        if position in sample_pp.keys():
                            pp_value=sample_pp[position][0]
                            pp=pp_value.split('/')
                            pp=int(pp[0])+int(pp[1])
                            if pp==1 and (sample_pp[position][1]>=0.2 and sample_pp[position][1]<0.85):
                                res[idx][position]  = hc_value
                                print(f'HC genotype is supported by bam')
                            if pp==2 and sample_pp[position][1]>=0.85:
                                print(f'HC genotype is not supported by bam, saving a new genotype')
                                print(sample_pp[position])
                                res[idx][position] = pp_value
                        else:
                            print('need to check inheritance')
                            inheritance = check_inheritance(idx, position, trios, gt_hc, gt_dv, gt_bam)
                            if inheritance == True:
                                print('saving this variant even with bad QUAL')
                                res[idx][position]  = hc_value
                            else:
                                print('droping this variant, probably artefact')
        else:
            hc_value=sample_hc[position][0]
            hc=hc_value.split('|')
            hc=int(hc[0])+int(hc[1])
            if hc == 2:
                if sample_hc[position][1] >= 0.85 and sample_hc[position][2]>=30:
                    res[idx][position]  = hc_value
                else:
                    print(idx)
                    print(f'check {position}: HC - {sample_hc[position]}, DV - not called, bam - {sample_pp.get(position, 'Not found in bam')}')
                    if position in sample_pp.keys():
                        pp_value=sample_pp[position][0]
                        pp=pp_value.split('/')
                        pp=int(pp[0])+int(pp[1])
                        if pp==2 and sample_pp[position][1]>=0.85:
                            res[idx][position]  = hc_value
                            print(f'genotype is supported by bam')
                        if pp==1 and (sample_pp[position][1]>=0.2 and sample_pp[position][1]<0.85):
                            print(f'HC genotype is not supported by bam, saving a new genotype')
                            res[idx][position]  = pp_value
                    else:
                        print('need to check inheritance')
                        inheritance = check_inheritance(idx, position, trios, gt_hc, gt_dv, gt_bam)
                        if inheritance == True:
                            print('saving this variant even with bad QUAL')
                            res[idx][position]  = hc_value
                        else:
                            print('droping this variant, probably artefact')
            if hc == 1:
                if (sample_hc[position][1] >= 0.2 and sample_hc[position][1] < 0.85):
                    res[idx][position]  = hc_value
                else:             
                    print(idx)
                    print(f'check {position}: HC - {sample_hc[position]}, DV - not called, {sample_pp.get(position, 'Not found in bam')}')
                    if position in sample_pp.keys():
                        pp_value=sample_pp[position][0]
                        pp=pp_value.split('/')
                        pp=int(pp[0])+int(pp[1])
                        if hc == pp and (sample_pp[position][1] >= 0.2 and sample_pp[position][1]<0.85):
                            res[idx][position]  = hc_value
                            print(f'HC genotype is supported by bam')
                        else:
                            print(f'HC genotype is not supported by bam, saving a new genotype')
                            print(sample_pp[position])
                            res[idx][position]  = pp_value
                    else:
                        print('need to check inheritance')
                        inheritance = check_inheritance(idx, position, trios, gt_hc, gt_dv, gt_bam)
                        if inheritance == True:
                            print('saving this variant even with bad QUAL')
                            res[idx][position]  = hc_value
                        else:
                            print('droping this variant, probably artefact')
    for position in sample_dv.keys():
        if position not in res[idx].keys():
            print(f'{position} in {idx} found in DV, but not in HC, checking it')
            dv_value=sample_dv[position][0]
            dv=dv_value.split('/')
            dv=int(dv[0])+int(dv[1])
            if dv == 2:
                if position in sample_pp.keys():
                    if sample_pp[position][1]> 0.85:
                        print(f'{position} not called by HC, saving it with {dv_value}, {sample_pp[position]}')
                        res[idx][position] = dv_value
                else:
                    print(f'low qual {position} in {idx}, checking inheritance; variant info: {dv_value}')
                    inheritance = check_inheritance(idx, position, trios, gt_hc, gt_dv, gt_bam)
                    if inheritance == True:
                        print('saving this variant even with bad QUAL')
                        res[idx][position] = dv_value
                    else:
                        print('droping this variant, probably artefact')
            if dv == 1:
                if position in sample_pp.keys():
                    if (sample_pp[position][1] >= 0.2 and sample_pp[position][1] < 0.85):
                        print(f'{position} not called by HC, saving it with {dv_value}, {sample_pp[position]}')
                        res[idx][position] = dv_value
                else:
                    print(f'low qual {position} in {idx}, checking inheritance; variant info: {dv_value}')
                    inheritance = check_inheritance(idx, position, trios, gt_hc, gt_dv, gt_bam)
                    if inheritance == True:
                        print('saving this variant even with bad QUAL')
                        res[idx][position] = dv_value
                    else:
                        print('droping this variant, probably artefact')