'''
utr_rho.py - iterate through file containing intergenic tracts
and return rho for adjacent 5' UTRs
'''

import argparse
from tqdm import tqdm
import antr
import csv

def args():
    parser = argparse.ArgumentParser(
        description="get rho for 5' UTRs", 
        usage='python3.5 script.py [options]')

    parser.add_argument('-f', '--fname', required=True,
                        type=str, help='File listing intergenic tracts')
    parser.add_argument('-t', '--table', required=True,
                        type=str, help='Annotation table with rho values')
    parser.add_argument('-o', '--outname', required=True,
                        type=str, help='File to write to')

    args = parser.parse_args()

    return args.fname, args.table, args.outname


def parse_tracts(fname, table, outname):
    with open(outname, 'w', newline='') as f_out:
        fieldnames = [
            'chrom', 'start', 'end', 'tract_size',
            'rho_vals', 'rho_count', 'rho_window',
            'utr_start', 'utr_end', 'utr_rho_vals',
            'utr_rho_count', 'utr_rho_window'
        ]
        f_out.write('\t'.join(fieldnames) + '\n')
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        lengths = antr.chlamy_lengths()
        with open(fname, 'r', newline='') as f_in:
            reader = csv.DictReader(f_in, delimiter='\t')
            for tract in tqdm(reader):
                utr_rho_vals = 0.0
                utr_rho_count = 0
                if int(tract['start']) > int(tract['end']):
                    continue
                else: 
                    chrom, utr_start = tract['chrom'], int(tract['end']) + 1
                    chrom_length = lengths[chrom]
                    p = antr.Reader(table)
                    first_iter = True
                    for record in p.fetch(chrom, utr_start, chrom_length):
                        if first_iter and not record.is_utr5 and not record.is_utr3:
                            print('wtf')
                            print(record.chrom, record.pos)
                            break
                        else:
                            first_iter = False
                        if record.is_utr5 or record.is_utr3:
                            utr_rho_vals += record.ld_rho
                            utr_rho_count += 1
                        elif not record.is_utr5 and not record.is_utr3 and not first_iter:
                            utr_end = record.pos - 1
                            out_dict = tract
                            out_dict['utr_start'] = utr_start
                            out_dict['utr_end'] = utr_end
                            out_dict['utr_rho_vals'] = utr_rho_vals
                            out_dict['utr_rho_count'] = utr_rho_count
                            out_dict['utr_rho_window'] = utr_rho_vals / utr_rho_count
                            writer.writerow(out_dict)
                            break

def main():
    fname, table, outname = args()
    parse_tracts(fname, table, outname)

if __name__ == '__main__':
    main()

        

