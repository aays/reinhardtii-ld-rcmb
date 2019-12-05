'''
snp_density.py - get SNP density in windows and append to LDhelmet file
'''

import argparse
from tqdm import tqdm
from cyvcf2 import VCF
import csv
from copy import deepcopy
import antr

def args():
    parser = argparse.ArgumentParser(
        description='get SNP density in windows', 
        usage='python3.5 snp_density.py [options]')

    parser.add_argument('-f', '--filename', required=True,
                        type=str, help='LDhelmet file (summarised in 2kb windows)')
    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='VCF (.vcf.gz)')
    parser.add_argument('-c', '--chrom', required=True,
                        type=str, help='Chromosome')
    parser.add_argument('-t', '--table', required=True,
                        type=str, help='Annotation table')
    parser.add_argument('-l', '--lookup', required=False,
                        type=str, help='Lookup string (if exists)')
    parser.add_argument('-o', '--outfile', required=True,
                        type=str, help='File to write to')

    args = parser.parse_args()

    return args.filename, args.vcf, args.chrom, args.table, args.lookup, args.outfile

def create_lookup(table, chrom):
    ''' (str, str) -> str
    uses annotation table to create a lookup string for all sites

    modified from rcmb_correlates.py

    i: intergenic
    c: cds
    n: intron
    5: utr5
    3: utr3
    N: unknown
    '''

    p = antr.Reader(table)
    lookup_string = ''
    start = next(p.fetch(chrom)).pos
    if start > 1:
        for i in range(0, start):
            lookup_string += 'N'

    for rec in tqdm(p.fetch(chrom)):
        if rec.is_intergenic:
            lookup_string += 'i'
        elif rec.is_in_CDS:
            lookup_string += 'c'
        elif rec.is_intronic:
            lookup_string += 'n'
        elif rec.is_utr5:
            lookup_string += '5'
        elif rec.is_utr3:
            lookup_string += '3'
        else:
            lookup_string += 'N'

    with open(chrom + '_temp_lookup', 'w') as f:
        f.write(lookup_string)

    return lookup_string

def count_annotations(lookup_string, start, end):
    ''' (str, int, int) -> dict
    uses lookup string to return annotation counts in window
    '''
    snippet = lookup_string[start:end]
    ant_dict = dict.fromkeys(['i', 'c', 'n', '5', '3', 'N'], 0)
    for ant in ant_dict:
        ant_dict[ant] = snippet.count(ant)

    return ant_dict
    

def count_window(vcf, chrom, start, end):
    ''' (str, str, int, int) -> int
    helper function for get_density that counts SNPs
    in a given window
    '''
    v = VCF(vcf)
    region = '{chrom}:{start}-{end}'.format(chrom=chrom, start=start, end=end)
    snp_count = len([record for record in v.__call__(region)])
    return snp_count


def get_density(filename, vcf, chrom, lookup_string, outfile):
    ''' (str, str, str, str) -> None
    iterates through summarised LDhelmet infile and appends
    columns containing SNP count and SNP density

    assumes VCF is pre-filtered and only contains SNPs - should
    be using the same filtered SNP-only VCF that was used for
    LDhelmet input anyways
    '''
    with open(filename, 'r', newline='') as f_in:
        reader = csv.DictReader(f_in)
        fieldnames = reader.fieldnames

        with open(outfile, 'w', newline='') as f_out:
            fieldnames.extend(['snp_count', 'snp_density'])
            fieldnames.extend([ant + '_count' for ant in ['i', 'c', 'n', '5',
                '3', 'N']])
            writer = csv.DictWriter(f_out, fieldnames=fieldnames)
            writer.writeheader()

            for line in tqdm(reader):
                start, end = int(line['block_start']), int(line['block_end'])
                window_size = end - start
                snp_count = count_window(vcf, chrom, start, end)
                snp_density = snp_count / window_size
                ant_dict = count_annotations(lookup_string, start, end)

                line_out = deepcopy(line)
                line_out['snp_count'] = snp_count
                line_out['snp_density'] = snp_density
                for ant in ant_dict:
                    line_out[ant + '_count'] = str(ant_dict[ant])

                writer.writerow(line_out)


def main():
    filename, vcf, chrom, table, lookup, outfile = args()
    if not lookup:
        print('Generating lookup string...')
        lookup_string = create_lookup(table, chrom)
    elif lookup:
        with open(lookup, 'r') as f:
            lookup_string = f.read().strip()
    print('SNP density and annotation counts for {chrom}'.format(chrom=chrom))
    # https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-cheaply-in-python
    line_count = sum(1 for line in open(filename))
    print('There are {l} lines in the file.'.format(l=line_count))
    get_density(filename, vcf, chrom, lookup_string, outfile)
    print('Done.')

if __name__ == '__main__':
    main()

        

