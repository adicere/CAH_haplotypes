import pysam
import json 
import sys 
import os 
import numpy as np

vcf_file = sys.argv[1] #path to merged vcf file 
o_dir = str(sys.argv[2]) #path to output directory

def VAF_calc(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    VAF = {}

    #calculating VAF for each variant of a sample and writing them in dictionary
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
    
    #calculating average VAF across samples 
    for key in VAF.keys():
        VAF[key] = round(np.mean(VAF.get(key)), 3)

    #writing results to json file in provided output directory
    with open(os.path.join(o_dir+'VAF.json'), 'w') as file:
        json.dump(VAF, file)

    return VAF


VAF_calc(vcf_file)