'''
get_COs_rho.py - quick and dirty script to get mean rho at each Liu CO
'''

import argparse
import antr
import csv
from tqdm import tqdm
from copy import deepcopy

def args():
    parser = argparse.ArgumentParser(
        description='get mean rho at each Liu CO', 
        usage='python3.5 get_COs_rho.py [options]')

    parser.add_argument('-f', '--filename', required=True,
                        type=str, help='Crossovers file')
    parser.add_argument('-t', '--table', required=True,
                        type=str, help='Annotation table w/ rho values.')
    parser.add_argument('-o', '--out', required=True,
                        type=str, help='File to write to.')

    args = parser.parse_args()

    return args.filename, args.table, args.out

def parse_crossovers(filename, table, out):
    with open(filename, 'r') as f:
        crossovers = [line for line in csv.DictReader(f, delimiter='\t')]
    with open(out, 'w') as f_out:
        fieldnames = ['cross', 'tetrad', 'individual', 'chromosome', 'left_bound', 'right_bound', 
            'mid_point', 'length', 'rho_total', 'rho_count', 'rho_window']
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()

        for co in tqdm(crossovers):
            chrom = str(co['chromosome'])
            start, end = int(co['left_bound']), int(co['right_bound'])
            p = antr.Reader(table)
            rho_vals = [record.ld_rho for record in p.fetch(chrom, start, end)]

            out_dict = deepcopy(co)
            out_dict['rho_total'] = sum(rho_vals)
            out_dict['rho_count'] = len(rho_vals)
            out_dict['rho_window'] = sum(rho_vals) / len(rho_vals)
            
            writer.writerow(out_dict)



def main():
    filename, table, out = args()
    parse_crossovers(filename, table, out)

if __name__ == '__main__':
    main()

        

