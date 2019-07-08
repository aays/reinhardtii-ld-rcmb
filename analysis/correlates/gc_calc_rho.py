'''
gc_calc_rho.py - calculate chromosomal GC content (total and at 4D sites) and rho from FASTA files
'''

import argparse
import antr
import sys
from tqdm import tqdm
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import OrderedDict

def args():
    parser = argparse.ArgumentParser(
        description='GC content and rho',
        usage='python3.5 gc_calc_rho.py [options]')

    parser.add_argument('-f', '--filename', required=True,
                        type=str, help='Input FASTA file')
    parser.add_argument('-a', '--annotation', required=True,
                        type=str, help='Annotation table')
    parser.add_argument('-w', '--windowsize', required=True,
                        type=int, help='Windowsize')
    parser.add_argument('-r', '--chrom', required=True,
                        type=str, help='Chromosome to calculate GC for')
    parser.add_argument('-o', '--outfile', required=True,
                        type=str, help='File to write to')

    args = parser.parse_args()

    return args.filename, args.annotation, args.windowsize, \
           args.chrom, args.outfile

def get_consensus(filename):
    ''' (str) -> str
    reads in input aligned fasta and generates a consensus
    consensus uses AlignIO's default threshold of 70%
    '''
    # https://www.biostars.org/p/284637/
    alignment = AlignIO.read(filename, 'fasta')
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus(threshold=0.7, ambiguous='N')
    return str(consensus)

def gc_content_calc(consensus, annotation, windowsize, chrom, outfile):
    ''' (str, str, int, str, str) -> None
    uses the consensus sequence generated above to calculate:
    
    1) total GC content
    2) GC content at selectively unconstrained sites (intronic, intergenic, 4D).
    3) cumulative rho + # of sites with rho estimates in window
    
    writes count of GC nucleotides + Ns + all non-N sites to specified outfile.
    '''
    print('Selected chromosome {chrom}'.format(chrom=chrom))
    lengths = antr.chlamy_lengths()
    chrom_length = lengths[chrom]
    with open(outfile, 'w') as f:
        f.write('start end GC GC4 N 4D_sites total_sites rho_total rho_count\n')
        seq_index = 0
        for window in tqdm(range(0, chrom_length, windowsize)):
            counter = OrderedDict.fromkeys(
                ['GC_count', 'GC4_count', 'N_count', '4D_sites', 'total_sites'], 0)
            rho = OrderedDict.fromkeys(
                ['rho_total', 'rho_count'], 0.0
            )
            window_start = window
            if window + windowsize > chrom_length:
                window_end = chrom_length
            else:
                window_end = window + windowsize
            p = antr.Reader(annotation)
            for record in p.fetch(chromosome, window_start, window_end):
                if not record.ld_rho == 'NA':
                    rho['rho_total'] += record.ld_rho
                    rho['rho_count'] += 1
                if record.is_fold4: 
                    if consensus[seq_index] in ['G', 'C']:
                        counter['GC_count'] += 1
                        counter['GC4_count'] += 1
                        counter['4D_sites'] += 1
                        counter['total_sites'] += 1
                    elif consensus[seq_index] in ['A', 'T']:
                        counter['4D_sites'] += 1
                        counter['total_sites'] += 1
                    elif consensus[seq_index] == 'N':
                        counter['N_count'] += 1
                elif not record.is_fold4: 
                    if consensus[seq_index] in ['G', 'C']:
                        counter['GC_count'] += 1
                        counter['total_sites'] += 1
                    elif consensus[seq_index] in ['A', 'T']:
                        counter['total_sites'] += 1
                    elif consensus[seq_index] == 'N':
                        counter['N_count'] += 1
                        
            window_out = ' '.join([str(num) for num in [window_start, window_end]])
            line_out_counts = ' '.join([str(i) for i in list(counter.values())])
            line_out_counter += ' '.join([str(i) for i in list(rho.values())])
            f.write(window_out + ' ' + line_out_counts + '\n')


def main():
    filename, annotation, windowsize, region, outfile = args()
    print('Obtaining consensus...')
    consensus = get_consensus(filename)
    print('Done.')
    print('Calculating GC content...')
    gc_content_calc(consensus, annotation, windowsize, region, outfile)
    print('Done.')
    print('Hooray!')

if __name__ == '__main__':
    main()

        
