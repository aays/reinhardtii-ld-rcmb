'''
annotation_tracts.py - create dataset of rho + polymorphism
using 'tracts' of annotations
'''

import argparse
from tqdm import tqdm
import antr
from cyvcf2 import VCF
import csv

def args():
    parser = argparse.ArgumentParser(
        description='create dataset of rho + polymorphism', 
        usage='python3.5 annotation_tracts.py [options]')

    parser.add_argument('-f', '--filename', required=True,
                        type=str, help='Annotation lookup filename')
    parser.add_argument('-t', '--table', required=True,
                        type=str, help='Annotation table w/ rho (.txt.gz)')
    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='VCF (.vcf.gz)')
    parser.add_argument('-c', '--chrom', required=True,
                        type=str, help='Chromosome')
    parser.add_argument('-o', '--out', required=True,
                        type=str, help='File to write to')

    args = parser.parse_args()

    return args.filename, args.table, args.vcf, args.chrom, args.out


def get_tract_snps(vcf, chrom, start, end):
    ''' (str, str, int, int) -> int
    helper function that returns # of variants in given window
    '''
    v = VCF(vcf)
    region = '{chrom}:{start}-{end}'.format(chrom=chrom, start=start, end=end)
    snp_count = len([record for record in v.__call__(region)])
    return snp_count

def get_tract_rho(table, chrom, start, end):
    ''' (str, str, int, int) -> (float, int)
    helper function that returns cumulative sum of rho + count of sites
    in input window
    '''
    p = antr.Reader(table)
    rho_vals, rho_count = 0.0, 0
    for record in p.fetch(chrom, start, end):
        rho_vals += record.ld_rho
        rho_count += 1
    return rho_vals, rho_count

def create_output_line(table, vcf, chrom, tract_start, tract_end,
        current_annotation):
    ''' (str, str, str, int, int) -> dict
    helper function that creates dict for csv.DictWriter in parse_tracts
    '''
    rho_vals, rho_count = get_tract_rho(table, chrom, tract_start, tract_end)
    snp_count = get_tract_snps(vcf, chrom, tract_start, tract_end)
    tract_length = tract_end - tract_start
    snp_density = snp_count / tract_length
    out_dict = {'chrom': chrom, 'start': tract_start, 'end': tract_end,
            'tract_length': tract_length, 'rho_sum': rho_vals, 'rho_count': rho_count, 
            'rho_tract': rho_vals / rho_count, 'snp_count': snp_count,
            'snp_density': snp_count / tract_length, 'is_intergenic': 0,
            'is_utr5': 0, 'is_in_CDS': 0, 'is_intronic': 0, 'is_utr3': 0}
    # assign annotation
    lookup = {'i': 'is_intergenic', '5': 'is_utr5', 'c': 'is_in_CDS',
              'n': 'is_intronic', '3': 'is_utr3'}
    for k in lookup:
        if current_annotation == k:
            out_dict[lookup[k]] = 1
    if current_annotation != 'i':
        out_dict['is_genic'] = 1
    elif current_annotation == 'i':
        out_dict['is_genic'] = 0
    return out_dict


def parse_tracts(filename, table, vcf, chrom, out):
    ''' (str, str, str, str, str) -> None
    iterates through lookup string and parses out continuous tracts of the same
    annotation, before passing to create_output_line to get snp counts etc
    '''
    fieldnames = ['chrom', 'start', 'end', 'tract_length', 'rho_sum',
        'rho_count', 'rho_tract', 'snp_count', 'snp_density', 'is_intergenic',
        'is_genic', 'is_utr5', 'is_in_CDS', 'is_intronic', 'is_utr3']

    with open(filename, 'r') as f:
        lookup_string = f.read().rstrip()
    print('There are {n} sites in the lookup string.'.format(n=len(lookup_string)))
    current_annotation = 'N'
    in_tract = False
    with open(out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i, annotation in tqdm(enumerate(lookup_string)):
            if annotation == 'N':
                continue
            elif not annotation == current_annotation: # new tract starts
                tract_start = i
                current_annotation = annotation
                continue
            elif annotation == current_annotation:
                try: # check for upcoming new tract
                    next_site = lookup_string[i + 1]
                except IndexError: # end of chrom/string
                    next_site = 'N'
                if current_annotation != next_site: # tract ends at next site
                    tract_end = i
                    out_dict = create_output_line(table, vcf, chrom, tract_start,
                            tract_end, current_annotation)
                    writer.writerow(out_dict)
                elif current_annotation == next_site:
                    continue
                
def main():
    filename, table, vcf, chrom, out = args()
    print('Parsing tracts for {chrom}'.format(chrom=chrom))
    parse_tracts(filename, table, vcf, chrom, out)
    print('Done!')

if __name__ == '__main__':
    main()

        

