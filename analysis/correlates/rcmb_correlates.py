'''
rcmb_correlates.py - use annotation table to correlate RR with annotations

hardcoded annotations are utr5, CDS, intronic, utr3, and intergenic
'''

import argparse
import antr
from tqdm import tqdm
from collections import OrderedDict

def args():
    parser = argparse.ArgumentParser(
        description = 'use annotation table to correlate RR with annotations', 
        usage='python3.5 rcmb_correlates.py [options]')

    parser.add_argument('-t', '--table', required=True,
                        type=str, help='Annotation table (.txt.gz)')
    parser.add_argument('-w', '--windowsize', required=True,
                        type=int, help='Window size')
    parser.add_argument('-x', '--gene_context', required=True,
                        type=int, help='Report rho for regions of x kb up/downstream of genes')
    parser.add_argument('-c', '--chrom', required=True,
                        type=str, help='Compute GC content [optional]')
    parser.add_argument('-o', '--out', required=True,
                        type=str, help='File to write to.')

    args = parser.parse_args()

    return args.table, args.windowsize, args.gene_context, args.chrom, args.out

# chromosome lengths - hardcoded for chlamy
lengths = {'chromosome_1': 8033585,
'chromosome_2': 9223677,
'chromosome_3': 9219486,
'chromosome_4': 4091191,
'chromosome_5': 3500558,
'chromosome_6': 9023763,
'chromosome_7': 6421821,
'chromosome_8': 5033832,
'chromosome_9': 7956127,
'chromosome_10': 6576019,
'chromosome_11': 3826814,
'chromosome_12': 9730733,
'chromosome_13': 5206065,
'chromosome_14': 4157777,
'chromosome_15': 1922860,
'chromosome_16': 7783580,
'chromosome_17': 7188315}

def create_lookup(table, context_size, chrom):
    '''(str, int, str) -> str
    uses annotation table to create a lookup string
    for all intergenic sites

    g: genic
    0: intergenic
    1: site upstream of gene
    2: site downstream of gene
    3: site upstream and downstream of gene
    N: unknown
    '''

    p = antr.Reader(table)
    genic_lookup = ''
    proximity_lookup = ''
    start = next(p.fetch(chrom)).pos
    if start > 1:
        for i in range(0, start):
            genic_lookup += 'N'
            proximity_lookup += 'N'

    for rec in p.fetch(chrom):
        if rec.is_genic:
            genic_lookup += 'g'
            continue
        elif rec.is_intergenic:
            genic_lookup += 'i'
            continue

    start = len(proximity_lookup)
    for pos in tqdm(range(start, len(genic_lookup))):
        if genic_lookup[pos] == 'g':
            proximity_lookup += 'g'
            continue
        elif genic_lookup[pos] == 'i':
            upstream, downstream = False, False

            if 'g' in genic_lookup[pos:pos+context_size]: upstream = True
            if 'g' in genic_lookup[pos-context_size:pos]: downstream = True

            if upstream and not downstream:
                proximity_lookup += '1'
                continue
            elif downstream and not upstream:
                proximity_lookup += '2'
                continue
            if upstream and downstream:
                proximity_lookup += '3'
                continue
            else:
                proximity_lookup += '0'

    return genic_lookup, proximity_lookup

def rho_annotations(table, windowsize, gene_context, chrom, out):
    '''(str, int, int, str, str) -> None
    iterates through annotation table and collects rho values for each annotation

    will first create a lookup string to speed up upstream/downstream/both calc
    '''

    print('Chromosome {chrom} selected.'.format(chrom=chrom))
    print('Creating intergenic lookup...')
    genic_lookup, proximity_lookup = create_lookup(table, gene_context, chrom)
    print('Done.')
    correlates = ['is_intergenic', 'is_utr5', 'is_in_CDS', 'is_intronic', 
                  'is_utr3', 'upstream', 'downstream', 'both']

    windows = list(range(0, lengths[chrom], windowsize)) + [lengths[chrom]]
    p = antr.Reader(table)

    print('Starting windowed correlate calc...')
    with open(out, 'w') as f:

        # prep header
        title1 = ' '.join([item + '_total' for item in correlates])
        title2 = ' '.join([item + '_count' for item in correlates])
        header = ' '.join(['chrom', 'start', 'end', title1, title2, 'total_count'])
        f.write(header + '\n')

        # iterate through chromosome
        for i in range(len(windows) - 1):
            rho = OrderedDict.fromkeys(correlates, 0.0)
            count = OrderedDict.fromkeys(correlates, 0)
            total_count = 0

            for record in tqdm(p.fetch(chrom, windows[i], windows[i+1])):
                for key in rho.keys():
                    if key in ['upstream', 'downstream', 'both']:
                        continue

                    elif getattr(record, key) and not record.ld_rho == 'NA':
                        if key == 'is_intergenic' and getattr(record, 'is_intergenic'):
                            intergenic_type = proximity_lookup[record.pos]
                            if intergenic_type == '1':
                                rho['upstream'] += record.ld_rho
                                count['upstream'] += 1
                            elif intergenic_type == '2':
                                rho['downstream'] += record.ld_rho
                                count['downstream'] += 1
                            elif intergenic_type == '3':
                                rho['both'] += record.ld_rho
                                count['both'] += 1
                            elif intergenic_type == '0':
                                rho['is_intergenic'] += record.ld_rho
                                count['is_intergenic'] += 1
                        else:
                            rho[key] += record.ld_rho
                            count[key] += 1
                            total_count += 1

            # prep output
            rhovals = list(rho.values())
            countvals = list(count.values())
            total_count = str(total_count)
            totals = ' '.join([str(v) for v in rhovals])
            counts = ' '.join([str(v) for v in countvals])

            line_out = ' '.join([chrom, str(windows[i]), str(windows[i+1]), totals, counts, total_count])
            f.write(line_out + '\n')
    print('Complete.')
    print('File written to {out}.'.format(out=out))
    print('Good job!')

def main():
    rho_annotations(*args())

if __name__ == '__main__':
    main()

        

