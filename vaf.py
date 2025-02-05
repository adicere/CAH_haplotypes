import pysam
import json 
import sys 
import os 
import numpy as np

vcf_file = sys.argv[1] #path to merged vcf file 
o_dir = str(sys.argv[2]) #path to output directory

def VAF_calc(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    variant_VAF = {}
    avg_VAF={}
    positions=list(range(32037500,32041346,1))
    #calculating VAF for each variant in CYP21A2 borders for a sample and writing them in dictionary
    for record in vcf:
        if record.pos in positions:
            for sample_name, sample_data in record.samples.items():
                if sample_name not in variant_VAF:
                    variant_VAF[sample_name] = []
                ad = sample_data.get('AD')

                if len(ad)!=1:
                    ref, alt = ad[0], ad[1]
                    dp = ref+alt
                    vaf = alt/dp if dp > 0 else 0
                    variant_VAF[sample_name].append(round(vaf, 3))
                else:
                    variant_VAF[sample_name].append(np.nan)
        else:
            for sample_name, sample_data in record.samples.items():
                variant_VAF[sample_name].append(np.nan)

                    
    
    #calculating average VAF across samples 
    # for key in avg_VAF.keys():
    #     avg_VAF[key] = round(np.mean(avg_VAF.get(key)), 3)
    for key in variant_VAF.keys():
        avg_VAF[key] = round(np.nanmean(variant_VAF.get(key)),3)
    #writing results to json file in provided output directory
    with open(os.path.join(o_dir+'variant_VAF.json'), 'w') as file:
        json.dump(variant_VAF, file)
    
    with open(os.path.join(o_dir+'avg_VAF.json'), 'w') as file:
        json.dump(avg_VAF, file)

    return variant_VAF


VAF_calc(vcf_file)