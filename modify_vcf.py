import pandas as pd
import sys

link = sys.argv[1]

def read_vcf(link): 
    vcf = pd.read_csv(link, sep='\t', comment='#', header=None,
          names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], index_col=False)  
    # process vcf
    contig_start = int(32037620)
    vcf = vcf.assign(CHROM='chr6',
                    POS=contig_start + vcf['POS'] - 1)

    return vcf

read_vcf(link).to_csv(sys.argv[2], header=None, index=None, sep='\t')
