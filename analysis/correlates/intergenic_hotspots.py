'''
intergenic_hotspots.py - get hotspot stats for each 
intergenic tract
'''

import argparse
from tqdm import tqdm
import antr
import csv
import sys
import numpy as np
from copy import deepcopy

def args():
    parser = argparse.ArgumentParser(
        description='hotspot stats for intergenic tracts', 
        usage='python3.5 intergenic_hotspots.py [options]')

    parser.add_argument('-f', '--fname', required=True,
                        type=str, help='2 kb recombination estimate csv file')
    parser.add_argument('-i', '--intergenic', required=True,
                        type=str, help='tsv file containing intergenic tracts')
    parser.add_argument('-o', '--out', required=True,
                        type=str, help='File to write to')

    args = parser.parse_args()

    return args.fname, args.intergenic, args.out


def create_hotspot_lookup(fname):
    """ (str) -> dict, dict
    iterates through full dataset of rho values summarised in 2 kb windows to
    create two dicts to be used as lookup by parse_tracts() below

    modelled on create_intergenic_lookups from utr3_rho.py

    second dict returned contains numpy arrays for quick indexing
    """
    with open(fname, 'r', newline='') as f:
        hotspots = {}
        reader = csv.DictReader(f)
        hotspot_count = 0
        counter = 0
        for record in tqdm(reader):
            counter += 1
            if record['chr'] not in hotspots.keys():
                hotspots[record['chr']] = {}
            if record['rate_ratio'] == 'NA':
                continue
            elif float(record['rate_ratio']) >= 5.0:
                start = int(record['block_start'])
                hotspots[record['chr']][start] = {
                    'chrom': record['chr'],
                    'end': int(record['block_end']),
                    'flank_rate': float(record['flank_rate']),
                    'block_rate': float(record['block_rate']),
                    'rate_ratio': float(record['rate_ratio'])
                    }
                hotspot_count += 1
        print('Hotspot lookup generated.')
        print('{h} of {t} records retained.'.format(h=hotspot_count,
            t=counter))
        assert len(hotspots.keys()) == 17 # hardcoding for chlamy
        hotspot_arrays = dict.fromkeys(hotspots.keys())
        for chrom in tqdm(hotspot_arrays):
            hotspot_arrays[chrom] = np.array(sorted(list(hotspots[chrom].keys())))
        return hotspots, hotspot_arrays


def parse_tracts(hotspots, hotspot_arrays, intergenic, out):
    """ (dict, dict, str, str) -> None
    uses lookup tables generated above to iterate through intergenic tract file
    and get hotspot stats for each tract

    takes into account three possible cases -
    1. hotspot starts before tract but ends within it
    2. hotspot is contained entirely within tract
    3. hotspot starts within tract, but ends outside it

    for cases 1 and 3, there can only be one hotspot for each by definition - 
    more than one found will break the fxn
    """
    with open(out, 'w', newline='') as f_out:
        fieldnames = [
            'chrom', 'start', 'end', 'rho_vals', 'rho_count',
            'tract_size', 'rho_window', 'flank_sum', 'block_sum',
            'rate_sum', 'n_hotspots', 'sites_in_hotspot'
            ]
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        with open(intergenic, 'r', newline='') as f_in:
            reader = csv.DictReader(f_in, delimiter='\t')
            for tract in tqdm(reader):
                tract_hotspots = np.array([]) # populate with any matching starts
                chrom = tract['chrom']
                flank_sum, block_sum, rate_sum = 0.0, 0.0, 0.0
                n_hotspots, sites_in_hotspot = 0, 0
                tract_start = int(tract['start'])
                tract_end = int(tract['end'])
                # hotspots that are fully contained in tract
                fully_within_tract = (hotspot_arrays[chrom] >= tract_start) & \
                        (hotspot_arrays[chrom] <= tract_end) & \
                        (hotspot_arrays[chrom] + 2000 <= tract_end)
                tract_hotspots = np.append(tract_hotspots, 
                        hotspot_arrays[chrom][fully_within_tract])
                sites_in_hotspot += 2000 * tract_hotspots.size
                # hotspots that start before tract but end in it
                end_within_tract = (hotspot_arrays[chrom] + 2000 >= tract_start) & \
                        (hotspot_arrays[chrom] + 2000 <= tract_start + 2000) # has to end within first 2 kb
                end_within = hotspot_arrays[chrom][end_within_tract]
                tract_hotspots = np.append(tract_hotspots, end_within)
                try:
                    assert end_within.size <= 1
                except:
                    print(tract, end_outside)
                    sys.exit()
                if end_within.size:
                    sites_in_hotspot += ((end_within + 2000) - tract_start)[0]
                # hotspots that start in tract but end outside
                end_outside_tract = (hotspot_arrays[chrom] <= tract_end) & \
                        (hotspot_arrays[chrom] + 2000 >= tract_end)
                end_outside = hotspot_arrays[chrom][end_outside_tract]
                tract_hotspots = np.append(tract_hotspots, end_outside)
                try:
                    assert end_outside.size <= 1
                except:
                    print(tract, end_outside)
                    sys.exit()
                if end_outside.size:
                    sites_in_hotspot += (tract_end - end_outside)[0]
                # populate out dict
                if tract_hotspots.size != 0:
                    for start in tract_hotspots:
                        full_dict = hotspots[chrom][start]
                        flank_sum += full_dict['flank_rate']
                        block_sum += full_dict['block_rate']
                        rate_sum += full_dict['rate_ratio']
                out_dict = deepcopy(tract)
                out_dict['flank_sum'] = flank_sum
                out_dict['block_sum'] = block_sum
                out_dict['rate_sum'] = rate_sum
                out_dict['n_hotspots'] = tract_hotspots.size
                out_dict['sites_in_hotspot'] = sites_in_hotspot
                writer.writerow(out_dict)


def main():
    fname, intergenic, out = args()
    print('Creating lookups...')
    hotspots, hotspot_arrays = create_hotspot_lookup(fname)
    print('Parsing tracts...')
    parse_tracts(hotspots, hotspot_arrays, intergenic, out)
    print('Done.')

if __name__ == '__main__':
    main()

        

