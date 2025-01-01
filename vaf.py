import pysam
import json 
import sys 
import os 
import numpy as np

vcf_file = sys.argv[1]
o_dir = sys.argv[2]

def VAF_calc(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    VAF = {}

    for record in vcf:
        for sample_name, sample_data in record.samples.items():
            if sample_name not in VAF:
                VAF[sample_name] = []
            ad = sample_data.get('AD')
            dp = sample_data.get('DP')

            if ad and dp:
                ref, alt = ad[0], ad[1]
                vaf = alt/dp if dp > 0 else 0
                VAF[sample_name].append(vaf)
            else:
                pass 

    for key in VAF.keys():
        VAF[key] = round(np.mean(VAF.get(key)), 3)

    with open(os.path.join(o_dir+'VAF.json'), 'w') as file:
        json.dump(VAF, file)

    return VAF


VAF_calc(vcf_file)