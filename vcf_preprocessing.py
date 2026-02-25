import pandas as pd
import numpy as np
import argparse
import os 
from pathlib import Path
import subprocess
import glob
import gzip
import sqlite3

starter = argparse.ArgumentParser(prog='CYP21A2 haplotyping tool', 
                                  description='Haplotype extraction for single and family samples based on CYP21A2 amplicon sequencing')

input_group = starter.add_mutually_exclusive_group(required=True)
input_group.add_argument('-l', '--list', type=str, help='Provide .tsv/.csv list with full paths to .vcf.gz files to analyze and extract haplotypes')
input_group.add_argument('-d', '--dir', type=str, help='Provide a directory with .vcf.gz files to analyze and extract and extract haplotypes')
starter.add_argument('-o', '--output', type=str, help='Provide an output directory for results', required=True)

def modify_pos(file: str | Path) -> None:
    vcf = pd.read_csv(file, sep='\t', comment='#', header=None, index_col=False)  
    vcf[0] = 'chr6'
    vcf[1] = int(32037619) + vcf[1].astype(int)
    vcf.to_csv(f"{args.output}/tmp/body.txt", header=None, index=None, sep='\t')

def normalize_coordinates(vcf: str | Path, filename: str) -> None:
    header = f"{args.output}/tmp/header.txt"
    body = f"{args.output}/tmp/body.txt"
    out_vcf = f"{args.output}/recalibrated/{filename}.recalibrated.vcf.gz"

    subprocess.run(['bcftools', 'view', '-h', vcf, '-o', header], check=True)
    modify_pos(vcf)

    with open(out_vcf, 'wb') as vcf:
        p = subprocess.Popen(['bgzip'], stdin=subprocess.PIPE, stdout=vcf)

        with open(header, 'rb') as h, open(body, 'rb') as b:
            p.stdin.write(h.read())
            p.stdin.write(b.read())
            
        p.stdin.close() 
        p.wait()

        if p.returncode!=0:
            raise RuntimeError('bgzip failed')
        
    subprocess.run(['tabix', f"{o_dir}/recalibrated/{filename}.recalibrated.vcf.gz"], check=True)

def process_vevat(sqlite_file: str | Path) -> pd.DataFrame:
    connection = sqlite3.connect(sqlite_file)
    vevat = pd.read_sql_query('SELECT * FROM Variant', connection)
    # vevat =vevat[vevat['base__hugo']=='CYP21A2']
    vevat = vevat[['base__chrom', 'extra_vcf_info__pos', 'extra_vcf_info__ref', 'extra_vcf_info__alt', 'dbsnp__rsid', 
    'base__cchange', 'base__achange', 'vep_samovar_acmg__score', 'vep_samovar_acmg__verdict', 
    'vep_samovar_acmg__criteria', 'vevatacmg_postaggregator__gene_region_custom',
    'vevatacmg_postaggregator__variant_type_custom']]
    vevat['pVCF'] = (vevat[['base__chrom', 
                    'extra_vcf_info__pos', 'extra_vcf_info__ref', 'extra_vcf_info__alt']].astype(str).agg('-'.join, axis=1))
    
    vevat = vevat.rename(columns={'base__chrom': 'CHROM', 'extra_vcf_info__pos': 'POS', 'extra_vcf_info__ref': 'REF', 
    'extra_vcf_info__alt': 'ALT','dbsnp__rsid': 'rsID', 'base__cchange': 'cDNA', 'base__achange': 'Protein', 
    'vep_samovar_acmg__score': 'SamoVarACMG Score', 'vep_samovar_acmg__verdict': 'SamoVarACMG Verdict', 
    'vep_samovar_acmg__criteria': 'SamoVarACMG Criteria', 'vevatacmg_postaggregator__gene_region_custom': 'GeneRegion_v2', 
    'vevatacmg_postaggregator__variant_type_custom': 'VariantType_v2'})

    vevat['rsID'] = vevat['rsID'].fillna(vevat['pVCF'])

    return vevat
    
def transform_vcf(df: pd.DataFrame, out_vcf: str) -> pd.DataFrame:
    required = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']

    pre_vcf=df.copy()
    pre_vcf = pre_vcf.filter(items=required) 
    pre_vcf['POS'] = pre_vcf['POS'].astype(int)
    pre_vcf = pre_vcf.sort_values(by='POS')
    pre_vcf['FILTER'] = 'PASS'
    pre_vcf['INFO'] = '.'
    if 'ID' in pre_vcf.columns:
        pre_vcf['ID'] = pre_vcf['ID'].fillna('.')
    else:
        pre_vcf['ID'] = '.'
    pre_vcf['FORMAT'] = 'GT:GQ:DP:AD:VAF:PL'
    pre_vcf['SAMPLE'] = '1/1:26:24:0,24:1.00:34,27,0'
    pre_vcf['QUAL'] = np.round(np.random.uniform(300, 1000, size=len(pre_vcf)),1)
    
    pre_vcf = pre_vcf[required]

    with open(out_vcf, "w", encoding="utf-8") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##source=Custom_CYP21A2_Study\n")
        f.write("##reference=GRCh38\n")
        f.write(f"##contig=<ID={pre_vcf['CHROM'][0]},length=138394717,assembly=GRCh38>\n")
        f.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">\n')
        f.write('##INFO=<ID=HGVS,Number=1,Type=String,Description="HGVS notation">\n')
        f.write('##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance from ClinVar">\n')
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
        f.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">\n')
        f.write('##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant Allele Fraction">\n')
        f.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">\n')
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")

        for idx, row in pre_vcf.iterrows():
            chrom = str(row["CHROM"])
            pos   = int(row["POS"])
            vid   = str(row["ID"])
            ref   = str(row["REF"])
            alt   = str(row["ALT"])
            qual  = int(row['QUAL'])
            filt  = str(row["FILTER"])
            info  = str(row["INFO"])
            fmt   = str(row["FORMAT"])
            sample_val = row.get("SAMPLE", "1/1:26:24:0,24:1.00:34,27,0")

            f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\t{fmt}\t{sample_val}\n")

    return pre_vcf


def merge_samples(path_dir: str | Path, annotation: pd.DataFrame) -> None:
    vcfs = glob.glob(f"{path_dir}/*.recalibrated.vcf.gz")
    vcf_dict = {Path(vcf).name.replace(".recalibrated.vcf.gz", ""): vcf for vcf in vcfs}
    order = annotation["patient"].astype(str)
    sorted_vcfs = [vcf_dict[p] for p in order if p in vcf_dict]
    missing = set(order) - set(vcf_dict.keys())
    if missing:
        raise ValueError(f"Missing VCFs for patients: {missing}")

    with open(f"{o_dir}/tmp/sample_list.txt", 'w') as f:
        f.write('\n'.join(sorted_vcfs))

    renamer = annotation[['sample_name', 'patient']]
    renamer.to_csv(f"{o_dir}/tmp/rename_merged.txt", sep='\t', header=False, index=False)

    subprocess.run(['./merging.sh', '-s', f"{o_dir}/tmp/sample_list.txt", 
                    '-r', f"{o_dir}/tmp/rename_merged.txt",
                    '-o', f"{o_dir}/tmp/"], check=True)

    

def annotate_vars(genotypes: str | Path, manual: pd.DataFrame) -> pd.DataFrame:
    num_header = 0
    with gzip.open(genotypes, 'rt') as f:
        for line in f.readlines():
            if line.startswith("##"):
                num_header += 1
            else:
                break
    vcf = pd.read_csv(genotypes, sep="\t", skiprows=num_header,  compression='gzip')
    vcf = vcf.rename({"#CHROM": "CHROM"}, axis=1)

    sm_ngp = '/primary/home/rutkovskaya.ea/sm-ngp/'
    os.makedirs(f"{o_dir}/tmp/vevat", exist_ok=True)
    transform_vcf(vcf, f"{o_dir}/tmp/cohort.vcf")
    
    subprocess.run(["conda", "run", "-n", "sm-ngp",
        'python', f"{sm_ngp}ngp.py", 'run', '-S', 'vep_annotation',  '-i', 
        f"{o_dir}/tmp/", '-o', f"{o_dir}/tmp/vevat", '--af-popmax-threshold', '0.99', '-w', '-A', "cohort"
    ], cwd=sm_ngp, check=True)

    vevat_ann = process_vevat(f"{o_dir}/tmp/vevat/cohort.vevat_annotated.sqlite")
    vevat_ann.to_csv(f"{o_dir}/results/cohort_variants_samovar.tsv", sep='\t', index=False)

    vcf_copy = vcf.copy()
    vcf_copy = vcf_copy.drop(vcf_copy.iloc[:, 9:], axis = 1)

    #splitting samples genotypes into chromosomes
    for col in vcf.iloc[:, 9:].columns:
        new_cols = vcf[col].str.split(r"[\/|]", expand=True)
        new_cols.columns = [f"{col}_1", f"{col}_2"]
        vcf_copy=vcf_copy.join(new_cols)
    
    vcf_copy=vcf_copy.drop(['ID', 'FILTER', 'INFO', 'FORMAT', 'QUAL'], axis=1)
    for col in vcf_copy.columns:
        if col.startswith('770'):
            vcf_copy[col] = vcf_copy[col].replace('.', None).astype('Int64')
    
    merged = pd.merge(vcf_copy, manual_ann[['CHROM', 'POS', 'REF', 'ALT', 'cDNA', 'Protein', 'rsID', 'Allele associated phenotype', 'ACMG']],
                      how='left', on=['CHROM', 'POS', 'REF', 'ALT'])
    merged = pd.merge(merged, vevat_ann[['CHROM', 'POS', 'REF', 'ALT', 'rsID', 'cDNA', 'Protein',
                                 'SamoVarACMG Verdict', 'GeneRegion_v2', 'VariantType_v2']], how='left',
                                 on=['CHROM', 'POS', 'REF', 'ALT'], suffixes=('_manual','_samovar'))
    
    merged['cDNA_manual'] = merged['cDNA_manual'].fillna(merged['cDNA_samovar'])
    merged['Protein_manual'] = merged['Protein_manual'].fillna(merged['Protein_samovar'])
    merged['rsID_manual'] = merged['rsID_manual'].fillna(merged['rsID_samovar'])
    merged['ACMG'] = merged['ACMG'].fillna(merged['SamoVarACMG Verdict'])
    merged_filtered = merged.drop(merged.columns[merged.columns.str.contains('_samovar|SamoVar')], axis=1)
    merged_filtered = merged_filtered.rename(columns={'cDNA_manual': 'cDNA', 'Protein_manual': 'Protein', 
'rsID_manual': 'rsID'})
    
    merged_filtered.to_csv(f"{o_dir}/results/genotypes_annotated.tsv", sep='\t', index=False)


args = starter.parse_args()

if args.dir:
    i_dir = os.path.abspath(args.dir)
o_dir = os.path.abspath(args.output)

folders = [f"{o_dir}/recalibrated", f"{o_dir}/results", f"{o_dir}/tmp"]
for dir in folders:
    os.makedirs(dir, exist_ok=True)

manual_ann = pd.read_csv('/primary/home/rutkovskaya.ea/haplotypes/CAH_haplotypes/references/CYP21A2_all_var_normalized.tsv', sep='\t')

if args.list:
    with open(args.list) as f:
        ann=pd.read_csv(f, sep='\t')
    
    # for idx, row in ann.iterrows():
    #     normalize_coordinates(row['vcf'], row['patient'])
    
    merge_samples(f"{o_dir}/recalibrated", ann)
    annotate_vars(f"{o_dir}/tmp/splitted_genotype.vcf.gz", manual_ann)
    # vevat_ann = process_vevat(f"{o_dir}/tmp/vevat/splitted_genotype.vevat_annotated.sqlite")
    # vevat_ann.to_csv(f"{o_dir}/results/cohort_variants_samovar", sep='\t', index=False)


if args.dir:
    for file in args.dir.glob('*.phased.vcf.gz'):
        if not file.is_file():
            continue
        
        name = os.path.basename(file)
        sample = name.split('.')[0]
        normalize_coordinates(file, sample)

def main() -> None:
    pass

if __name__ == '__main__':
    main()