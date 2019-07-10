'''
hotspot_annotation.py - get # of sites from each annotation that are in hotspots
'''

import sys
import argparse
from tqdm import tqdm
import antr
import csv

def args():
    parser = argparse.ArgumentParser(
        description='get # of sites from each annotation that are in hotspots', 
        usage='python3.5 hotspot_annotation.py [options]')

    parser.add_argument('-f', '--filename', required=True,
                        type=str, help='2 kb recombination estimate csv file.')
    parser.add_argument('-t', '--table', required=True,
                        type=str, help='Annotation table.')
    parser.add_argument('-c', '--chrom', required=True,
                        type=str, help='Chromosome')
    parser.add_argument('-o', '--out', required=True,
                        type=str, help='File to write to.')

    args = parser.parse_args()

    return args.filename, args.table, args.chrom, args.out


def create_full_lookup(table, context_size, chrom):
    '''(str, int, str) -> str
    uses annotation table to create a lookup string
    for all sites (modified from rcmb_correlates.py)
    c: CDS
    i: intronic
    f: utr5
    t: utr3
    0: intergenic (non-gene-proximate)
    1: site upstream of gene
    2: site downstream of gene
    3: site upstream and downstream of gene
    N: unknown

    in the lookup string, string[x] represents position x + 1
    '''

    p = antr.Reader(table)
    initial_lookup = ''
    genic_lookup = ''
    proximity_lookup = ''
    start = next(p.fetch(chrom)).pos
    if start > 1:
        for i in range(0, start):
            initial_lookup += 'N'
            genic_lookup += 'N'
            proximity_lookup += 'N'
    
    lookup_codes = {
        'is_in_CDS': 'c',
        'is_intronic': 'i',
        'is_utr5': 'f',
        'is_utr3': 't',
        'is_intergenic': '0'
    }
    genic_codes = ['c', 'i', 'f', 't']

    print('Initial run through...')
    for rec in tqdm(p.fetch(chrom)):
        if rec.is_genic:
            genic_lookup += 'g'
        elif rec.is_intergenic:
            genic_lookup += 'i'
        for annotation in lookup_codes.keys():
            if getattr(rec, annotation):
                initial_lookup += lookup_codes[annotation]
                break
    
    print('Done.')
    print('Generating full string...')
    start = len(proximity_lookup)
    for pos in tqdm(range(start, len(genic_lookup))):
        if genic_lookup[pos] == 'g':
            proximity_lookup += initial_lookup[pos]
        elif genic_lookup[pos] == 'i':
            upstream, downstream = False, False

            if 'g' in genic_lookup[pos:pos+context_size]:
                upstream = True
            if 'g' in genic_lookup[pos-context_size:pos]:
                downstream = True

            if upstream and not downstream:
                proximity_lookup += '1'
                continue
            elif downstream and not upstream:
                proximity_lookup += '2'
                continue
            elif upstream and downstream:
                proximity_lookup += '3'
                continue
            else:
                proximity_lookup += '0'

    return genic_lookup, proximity_lookup
    

def parse_annotations(filename, table, chrom, out):
    ''' (str, str, str, str)
    generates lookups w/ create_full_lookup and creates two counters -
    one keeping track of hotspot annotations, and the other total annotations
    '''

    correlates = [
        'intergenic', 'intronic', 'utr5', 'utr3', 
        'CDS', 'upstream', 'downstream', 'both'
    ]

    hotspot_counter = dict.fromkeys(correlates, 0)
    total_counter = dict.fromkeys(correlates, 0)

    print('Generating lookup string for {chrom}...'.format(chrom=chrom))

    genic_lookup, proximity_lookup = create_full_lookup(table, 2000, chrom)

    print('Done.')
    print('Parsing annotations...')

    reverse_codes = {
        'c': 'CDS', 'i': 'intronic', 'f': 'utr5', 't': 'utr3', '0': 'intergenic', 
        '1': 'upstream', '2': 'downstream', '3': 'both', 'N': 'N'
    }
    
    with open(filename, 'r') as f:
        reader = csv.DictReader(f, delimiter=',')
        for record in tqdm(reader):
            if record['rate_ratio'] == 'NA':
                continue
            elif float(record['rate_ratio']) < 5.0:
                start, end = int(record['block_start']), int(record['block_end'])
                if end > len(proximity_lookup):
                    end = len(proximity_lookup)
                for pos in range(start, end):
                    annotation = reverse_codes[proximity_lookup[pos]]
                    if annotation != 'N':
                        total_counter[annotation] += 1
            elif float(record['rate_ratio']) >= 5.0: # hotspots only
                start, end = int(record['block_start']), int(record['block_end'])
                if end > len(proximity_lookup):
                    end = len(proximity_lookup)
                for pos in range(start, end):
                    annotation = reverse_codes[proximity_lookup[pos]]
                    if annotation != 'N':
                        hotspot_counter[annotation] += 1
                        total_counter[annotation] += 1

    with open(out, 'w') as f_out:
        f_out.write('correlate hotspot_count total_count\n')
        for correlate in correlates:
            line_out = ' '.join([
                correlate, str(hotspot_counter[correlate]), str(total_counter[correlate])
            ])
            f_out.write(line_out + '\n')

def main():
    filename, table, chrom, out = args()
    parse_annotations(filename, table, chrom, out)

if __name__ == '__main__':
    main()

