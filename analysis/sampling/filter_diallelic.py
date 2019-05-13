'''
filter_diallelic.py - remove nondiallelic SNPs from VCFs already filtered w/ bcftools

usage:
python3.5 filter_diallelic.py [input vcf] [output vcf]

'''

import sys
from cyvcf2 import VCF, Writer
from tqdm import tqdm

fname = sys.argv[-2]
outname = sys.argv[-1]

vcf_in = VCF(fname)
vcf_out = Writer(outname, vcf_in)

kept_count = 0
total_count = 0

print('Reading from {}'.format(fname))
print('Writing to {}'.format(outname))

for record in tqdm(vcf_in):
    total_count += 1
    if len(record.ALT) == 1:
        vcf_out.write_record(record)
        kept_count += 1
    else:
        continue

print('{0} records kept out of {1} total.'.format(kept_count, total_count))
print('{} records removed.'.format(total_count - kept_count))
