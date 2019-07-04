'''
add_ld_rho.py - takes in LDhelmet files + annotation table
to create new annotation table w/ LDhelmet rho values appended
'''

import argparse
import ant
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = 'add rho vals to annotation table',
                                     usage = 'python3.5 add_ld_rho.py [options]')

    parser.add_argument('-t', '--table', required = True,
                        type = str, help = 'Annotation table file')
    parser.add_argument('-d', '--ldhelmet_dir', required = True,
                        type = str, help = 'Directory containing chromosomal LDhelmet outputs')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'File to write to')

    args = parser.parse_args()

    return args.table, args.ldhelmet_dir, args.outfile

def prep_header(table, outfile):
    with open(outfile, 'w') as f_out:
        p = ant.Reader(table)
        for item in p.header:
            f_out.write(item.strip() + '\n')
        f_out.write('##ld_rho=LDhelmet recombination rate on Quebec dataset. p/bp [FLOAT]' + '\n')

        # add column names
        p.cols.append('ld_rho')
        p.cols = [item.strip() for item in p.cols]
        col_line = '#' + '\t'.join(p.cols) + '\n'
        f_out.write(col_line)

def write_lines(table, ldhelmet_dir, outfile):
    for i in range(1, 18):
        with open(outfile, 'a') as f_out:
            with open(ldhelmet_dir + 'chromosome_{}.txt'.format(i)) as f:
                for line in tqdm(f):
                    if line.startswith(('#', 'ver')):
                        continue
                    else:
                        split = line.split(' ')
                        start, end, rho = int(split[0]), int(split[1]), float(split[2])
                        p = ant.Reader(table).fetch('chromosome_{}'.format(i), start - 1, end - 1, raw=True)
                        for record in p:
                            record = record + '\t' + str(rho) + '\n'
                            f_out.write(record)

def main():
    table, ldhelmet_dir, outfile = args()
    prep_header(table, outfile)
    write_lines(table, ldhelmet_dir, outfile)

if __name__ == '__main__':
    main()


