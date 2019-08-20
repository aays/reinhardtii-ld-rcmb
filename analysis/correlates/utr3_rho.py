'''
utr3_rho.py - iterate through file containing intergenic tracts
and three_prime_UTR lines in GFF

will return list of UTRs, their rho vals/counts, and the nearest
intergenic tracts - could then left join with the actual intergenic
tract data to compare rho values

need to first create filtered GFF containing just 3' UTRs, starts,
ends, and info column

ie grep 'three_prime_UTR' [gff] | cut -f 1,4,5,9 > [bed]
'''

import argparse
from tqdm import tqdm
import antr
import csv
from copy import deepcopy
import numpy as np
import sys

def args():
    parser = argparse.ArgumentParser(
        description="get rho for 3' UTRs", 
        usage='python3.5 utr3_rho.py [options]')

    parser.add_argument('-b', '--bed', required=True,
                        type=str, help='Filtered bed file with utr3s')
    parser.add_argument('-f', '--fname', required=True,
                        type=str, help='File listing intergenic tracts')
    parser.add_argument('-t', '--table', required=True,
                        type=str, help='Annotation table with rho values')
    parser.add_argument('-o', '--out', required=True,
                        type=str, help='File to write to')

    args = parser.parse_args()

    return args.bed, args.fname, args.table, args.out

def create_intergenic_lookups(fname):
    ''' (str) -> dict, dict
    reads in intergenic tract tsv from intergenic_tract_lengths.py
    to create nested dict lookup table

    first layer of dict keys are chr names and second layer
    has starts as keys and ends as values

    ie lookup['chromosome_1'][673] returns 18765 (end of tract)

    second dict is just np array of starts for quick lookup
    '''
    with open(fname, 'r', newline='') as f:
        records = [{k: record[k] for k in ['chrom', 'start', 'end']}
                    for record in csv.DictReader(f, delimiter='\t')]
        lookup = {}
        print('Creating first lookup...')
        for record in tqdm(records):
            if record['chrom'] not in lookup.keys():
                lookup[record['chrom']] = {}
            if record['chrom'] in lookup.keys():
                lookup[record['chrom']][int(record['start'])] = int(record['end'])
        print('Creating second lookup...')
        start_arrays = dict.fromkeys(lookup.keys())
        for chrom in tqdm(start_arrays):
            start_arrays[chrom] = np.array(sorted(list(lookup[chrom].keys())))
    return lookup, start_arrays

def parse_utrs(bed, lookup, start_arrays, table, out):
    ''' (str, dict, dict, str, str) -> None
    uses lookups made above + annotation to both get rho in
    3' UTRs and 'assign' them their nearest intergenic tracts
    '''
    with open(out, 'w', newline='') as f_out:
        fieldnames_out = [
            'chrom', 'utr3_start', 'utr3_end', 'start', 'end',
            'utr3_rho_vals', 'utr3_rho_count', 'utr3_rho_window'
            ]
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames_out)
        writer.writeheader()
        with open(bed, 'r', newline='') as f_in:
            fieldnames = ['chrom', 'start', 'end', 'info']
            reader = csv.DictReader(f_in, delimiter='\t', fieldnames=fieldnames)
            print('Parsing UTRs...')
            for utr in tqdm(reader):
                chrom = utr['chrom']
                start = int(utr['start'])
                end = int(utr['end'])
                possible_starts = start_arrays[chrom][start_arrays[chrom] > end]
                if possible_starts.size:
                    tract_start = possible_starts.min()
                else:
                    continue # UTR is next to 'end of chromosome' tract
                utr3_rho_vals = 0.0
                utr3_rho_count = 0
                p = antr.Reader(table)
                for record in p.fetch(chrom, start, end):
                    try:
                        assert record.is_utr3
                    except:
                        print('Error - not UTR?')
                        print(record.chrom, record.pos)
                        sys.exit()
                    if record.ld_rho != 'NA':
                        utr3_rho_vals += record.ld_rho
                        utr3_rho_count += 1
                try:
                    utr3_rho_window = utr3_rho_vals / utr3_rho_count
                except ZeroDivisionError:
                    utr3_rho_window = 0
                out_dict = {
                    'chrom': chrom, 'utr3_start': start, 'utr3_end': end,
                    'start': tract_start, 'end': lookup[chrom][tract_start],
                    'utr3_rho_vals': utr3_rho_vals, 'utr3_rho_count': utr3_rho_count,
                    'utr3_rho_window': utr3_rho_window}
                writer.writerow(out_dict)


def main():
    bed, fname, table, out = args()
    lookup, start_arrays = create_intergenic_lookups(fname)
    parse_utrs(bed, lookup, start_arrays, table, out)
    print('Done!')
    
if __name__ == '__main__':
    main()


