import pandas as pd
from collections import defaultdict
import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from tqdm import tqdm

#Dataframe with clinical metadata about patients - CAH clinical form, kinship between patients etc
clin_data=pd.read_csv('/home/rutkovskaya.ea/haplotypes/text_files/amplicon_annotation_processed.tsv', sep='\t', index_col=0)
clin_data['patient'] = clin_data['patient'].fillna(0).astype(int)

#Phased haplotypes in families
ann_variants=pd.read_excel('/home/rutkovskaya.ea/haplotypes/text_files/rs.xlsx')
ann_variants = ann_variants.fillna(0)

#Extracting all haplotypes 
haplotypes={} 
for col in ann_variants.iloc[:, 4:676].columns:
    haplotype = tuple([x for x in ann_variants[col] if x!=0])
    if haplotype not in haplotypes.keys():
        haplotypes[haplotype]=[col]
    else:
        haplotypes[haplotype].append(col)

#Creating a dataframe with haplotypes assigned with temporary name
TN_db=defaultdict(dict)
i=0
for key in haplotypes.keys():
    if key:
        i+=1
        TN_db['T'+ str(i)]['Variants'] = key
        TN_db['T'+ str(i)]['Samples'] = haplotypes[key]
        TN_db['T'+ str(i)]['Counts'] = len(haplotypes[key])
        TN_db['T'+ str(i)]['Freq'] = len(haplotypes[key])/684
    else:
        print(f'{key}, {haplotypes[key]}')
        wt = {key: haplotypes[key]}

TN_db_sorted = dict(
    sorted(TN_db.items(), key=lambda item: item[1]['Freq'], reverse=True))

#Gathering the statistics of discovered haplotypes in all groups of patients and their zygosity state
df = pd.DataFrame([{'Haplotype': h, 'Variants': ', '.join(d['Variants']), 'Samples': ', '.join(d['Samples']),
'Hap_count': d['Counts'], 'Cohort_freq': d['Freq']} for h, d in TN_db_sorted.items()])
df = df.assign(**{col: '' for col in ['CAH_hom', 'CAH_het', 'Healthy_hom', 'Healthy_het', 'Clinical_form']})

for idx, row in df.iterrows():
    samples=[]
    cah_hom=0
    cah_het=0
    h_hom=0
    h_het=0
    forms={'SW': 0, 'SV': 0, 'NC': 0, 'N/S': 0}
    alleles=list(row['Samples'].split(', '))
    for al in alleles:
        samples.append(al.split('_')[0])
    alleles = list(set(samples))

    for allele in alleles:
        if clin_data.loc[clin_data['patient']==int(allele), 'group'].values[0]=='ILL':
            if samples.count(allele) == 2:
                cah_hom+=1
            else:
                cah_het+=1
            if clin_data.loc[clin_data['patient']==int(allele), 'form'].values[0]=='SW':
                forms['SW'] +=1 
            if clin_data.loc[clin_data['patient']==int(allele), 'form'].values[0]=='SV':
                forms['SV'] +=1
            if clin_data.loc[clin_data['patient']==int(allele), 'form'].values[0]=='NC':
                forms['NC'] +=1
            else:
                forms['N/S'] +=1
        else:
            if samples.count(allele) ==2:
                h_hom+=1
            else:
                h_het+=1
           
    df.at[idx, 'CAH_hom'] = cah_hom
    df.at[idx, 'CAH_het'] = cah_het
    df.at[idx, 'Healthy_hom'] = h_hom
    df.at[idx, 'Healthy_het'] = h_het
    df.at[idx, 'Clinical_form'] = '; '.join(f"{k}: {v}" for k, v in forms.items() if v != 0)

df.to_excel('/home/rutkovskaya.ea/haplotypes/text_files/sorted_haplo.xlsx', index=False)

#Classifying variants into majors and minors and assigning haplotypes in accordance with developed nomenclature
arr = ['inframe', 'missense',
       'splicing', 'frameshift', 'stop']

#these minor variants are considered as exceptions due to theirs association with CAH in reports and articles
exceptions=['rs1246774295', 'rs6467_1', 'rs1463196531'] 


majors = defaultdict(list)
minors=[]
for k,v in TN_db_sorted.items():
    haplo = TN_db_sorted[k]['Variants']
    major_vars=[]
    for var in haplo:
        if (ann_variants.loc[ann_variants['rs#'] == var, 'Effect'].values[0] in arr) or var in exceptions:
            major_vars.append(var)
        
    if len(major_vars) > 0:
        major_vars = tuple(major_vars)
        if major_vars not in majors:
            majors[major_vars] = [k]
        else:
            majors[major_vars].append(k)
    else:
        minors.append(k)

assigned=defaultdict(dict)
assigned['CYP21A2#1']['Temporary_name'] = ''
assigned['CYP21A2#1']['Major variants'] = ' '
assigned['CYP21A2#1']['Phenotype accociated'] = ' '
assigned['CYP21A2#1']['Full variants'] = ' '
assigned['CYP21A2#1']['Samples'] = wt[()][0]

pats=['SW', 'SV', 'NC', 'CL']
for i, minor in enumerate(minors, start=1):
    haplo_name = f'CYP21A2#1.{i:03d}'
    assigned[haplo_name]['Temporary_name'] = minor
    assigned[haplo_name]['Major variants'] = ' '
    assigned[haplo_name]['Phenotype accociated'] = ' '
    assigned[haplo_name]['Full variants'] = df.loc[df['Haplotype'] == minor, 'Variants'].values[0]
    assigned[haplo_name]['Samples'] = df.loc[df['Haplotype'] == minor, 'Hap_count'].values[0]

for i, major in enumerate(majors, start=2):
    phen_acc=[]
    haplo_name = f'CYP21A2#{i}'
    assigned[haplo_name]['Temporary_name'] = ''
    assigned[haplo_name]['Major variants'] = ', '.join(major)
    assigned[haplo_name]['Phenotype accociated'] = ', '.join(phen_acc)
    assigned[haplo_name]['Full variants'] = ''
    assigned[haplo_name]['Samples'] = ''
    
    for rs in major:
        if ann_variants.loc[ann_variants['rs#'] == rs, 'Allele associated phenotype'].values[0] in pats:
            phen_acc.append(rs)
    for j, subtype in enumerate(majors[major], start=1):
        haplo_name = f'CYP21A2#{i}.{j:03d}'
        assigned[haplo_name]['Temporary_name'] = subtype
        assigned[haplo_name]['Major variants'] = ', '.join(major)
        assigned[haplo_name]['Phenotype accociated'] = ', '.join(phen_acc)
        assigned[haplo_name]['Full variants'] = df.loc[df['Haplotype'] == subtype, 'Variants'].values[0]
        assigned[haplo_name]['Samples'] = df.loc[df['Haplotype'] == subtype, 'Hap_count'].values[0]

nom_df=pd.DataFrame.from_dict(assigned)
nom_df = nom_df.transpose()

#Searching for chimeras in found haplotypes
nom_df['Chimera'] = ''

known_chimeras ={
        'CH-1': ['rs9378251', 'rs6467_1', 'rs387906510'],
        'CH-2': ['rs9378251', 'rs6467_1', 'rs387906510', 'rs6475'], 
        'CH-3': ['rs9378251', 'rs6467_1', 'rs387906510', 'rs6475', 'rs1554299737', 'rs12530380', 'rs6476', 'rs6471', 'rs267606756', 'rs7755898'], 
        'CH-5': ['rs9378251', 'rs6467_1', 'rs387906510', 'rs6475', 'rs1554299737', 'rs12530380', 'rs6476', 'rs267606756'], 
        'CH-6': ['rs9378251', 'rs6467_1'],
        'CH-7': ['rs9378251', 'rs6467_1', 'rs387906510', 'rs6475', 'rs1554299737', 'rs12530380', 'rs6476'], 
        'CH-8': ['rs9378251', 'rs6467_1', 'rs387906510', 'rs6475', 'rs1554299737', 'rs12530380', 'rs6476', 'rs6471', 'rs267606756', 'rs7755898', 'rs7769409'], 
        'Attenuated': ['rs9378251']
        }

deletion_30kb=set(['rs1246774295', 'rs909177624', 'rs573835051'])

for idx, row in nom_df.iterrows():
    best_match=''
    if deletion_30kb.issubset(set(row['Full variants'].split(', '))):
        for chimera in known_chimeras:
            if set(known_chimeras[chimera]).issubset(set(row['Full variants'].split(', '))):
                if not best_match:
                    best_match = chimera
                else:
                    if len(known_chimeras[chimera]) > len(known_chimeras[best_match]):
                        best_match=chimera
        if not best_match:
            print(row['Temporary_name'])
    nom_df.at[idx, 'Chimera'] = best_match
    
nom_df.to_excel('/home/rutkovskaya.ea/haplotypes/text_files/cyp21a2_assigned_v2.xlsx')


#Extracting diplotypes in patients and count them based on clinical group
patients_haplo=defaultdict(list)
for t_haplo in TN_db_sorted:
    for allele in TN_db_sorted[t_haplo]['Samples']:
        sample = allele.split('_')[0]
        patients_haplo[sample].append(t_haplo)

diplotypes=[]
for sample in patients_haplo:
    if len(patients_haplo[sample])>1:
        hap1 = re.search(r'#[0-9]', nom_df.loc[nom_df['Temporary_name'] == patients_haplo[sample][0]].index[0])[0]
        hap2 = re.search(r'#[0-9]', nom_df.loc[nom_df['Temporary_name'] == patients_haplo[sample][1]].index[0])[0]
    else:
        hap1 = '#1'
        hap2 = re.search(r'#[0-9]', nom_df.loc[nom_df['Temporary_name'] == patients_haplo[sample][0]].index[0])[0]
    diplotype = '/'.join(sorted([hap1, hap2]))
    group = clin_data.loc[clin_data['patient'] == int(sample), 'form'].values[0]
    diplotypes.append({
        'Diplotype': diplotype,
        'Patient': sample,
        'Phenotype': group
    })

diplo = pd.DataFrame(diplotypes)
diplo['Phenotype'] = diplo['Phenotype'].fillna('Healthy')
grouped = diplo.groupby(['Diplotype', 'Phenotype'])['Patient'].agg(['count', lambda x: list(x)]).reset_index()
grouped.columns = ['Diplotype', 'Phenotype', 'Count', 'Patients']

result = grouped.pivot(index='Diplotype', columns='Phenotype', values=['Patients', 'Count'])
result.columns = [f"{p}_{stat}" for stat, p in result.columns]  # flatten MultiIndex

result = result.fillna('')

for col in result.columns:
    if col.endswith('_Patients'):
        result[col] = result[col].apply(lambda x: ', '.join(map(str, x)) if isinstance(x, list) else x)

result.to_excel('/home/rutkovskaya.ea/haplotypes/text_files/diplotypes.xlsx')