'''
intergenic_tract_lengths.py - use annotation table to get length of intergenic tracts
and rho vals + counts in them
'''

import argparse
from tqdm import tqdm
import antr

def args():
    parser = argparse.ArgumentParser(
        description='get length of intergenic tracts', 
        usage='python3.5 intergenic_tract_length.py [options]')

    parser.add_argument('-t', '--table', required=True,
                        type=str, help='Annotation table')
    parser.add_argument('-o', '--out', required=True,
                        type=str, help='File to write to')

    args = parser.parse_args()

    return args.table, args.out

def get_lengths(table, out):
    with open(out, 'w') as f:
        colnames = [
            'chrom', 'start', 'end', 'rho_vals',
            'rho_count', 'tract_size', 'rho_window'
        ]
        f.write('\t'.join(colnames) + '\n')
        p = antr.Reader(table)
        in_tract = False
        rho_vals = 0.0
        rho_count = 0
        for record in tqdm(p):
            if record.is_intergenic:
                if not in_tract:
                    in_tract = True
                    current_chrom = record.chrom
                    start = record.pos
                    if record.ld_rho != 'NA':
                        rho_vals += record.ld_rho
                        rho_count += 1
                elif in_tract:
                    if record.ld_rho != 'NA' and record.chrom == current_chrom:
                        rho_vals += record.ld_rho
                        rho_count += 1
                    elif record.chrom != current_chrom: # hit end of chrom
                        in_tract = False
                        end = record.pos
                        out = [
                            record.chrom, start, end, rho_vals, 
                            rho_count, end - start - 1, rho_vals / rho_count
                        ]
                        out = [str(item) for item in out]
                        f.write('\t'.join(out) + '\n')
                        # reset
                        rho_vals = 0.0
                        rho_count = 0
            elif record.is_genic:
                if in_tract:
                    in_tract = False
                    end = record.pos - 1
                    if start == end:
                        tract_size = 0
                    else:
                        tract_size = end - start - 1
                    out = [
                        record.chrom, start, end, rho_vals, 
                        rho_count, tract_size, rho_vals / rho_count
                    ]
                    out = [str(item) for item in out]
                    f.write('\t'.join(out) + '\n')
                    # reset
                    rho_vals = 0.0
                    rho_count = 0
                elif not in_tract:
                    continue


def main():
    table, out = args()
    get_lengths(table, out)

if __name__ == '__main__':
    main()

        

