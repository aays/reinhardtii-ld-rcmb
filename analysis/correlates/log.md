
to do:
- calculate recombination rates over all annotations in the genome
- calculate which annotations are enriched for hotspots

## 2/7/2019

can port in earlier scripts for this

first - need to create a version of the annotation table w/ LDhelmet rho values tacked on

```bash
# getting files/scripts in order

cd data/correlates
ln -sv /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/annotation/concatenated_GFF/annotation_table.txt.gz* .

pro
cd analysis/correlates
cp -v ../../../analysis/fiftyshadesofgreen/annotation_parser/add_ld_rho.py .
cp -v ../../../analysis/fiftyshadesofgreen/annotation_parser/ant.py .
cp -v ../../../analysis/fiftyshadesofgreen/annotation_parser/antr.py .
```

running revamped `add_ld_rho.py` to create new annotation table:

```bash
time python3.5 analysis/correlates/add_ld_rho.py \
--table data/correlates/annotation_table.txt.gz \
--ldhelmet_dir data/ldhelmet/block_100 \
--outfile data/correlates/annotation_table_rho.txt
```

## 3/7/2019

that took about 11 hours - now to bgzip and tabix:

```bash
bgzip annotation_table_rho.txt
tabix -p vcf annotation_table_rho.txt.gz
```

porting annotation table parser and correlate estimation script:

```bash
pro
cd analysis/correlates
cp -v ../../../analysis/fiftyshadesofgreen/annotation_parser/antr_correlate_dict_chrom.py .
```

test run after renaming:

```bash
mkdir -p data/correlates/chrom

time python3.5 analysis/correlates/ant_correlate_chrom.py \
--table data/correlates/annotation_table_rho.txt.gz \
--windowsize 100000 \
--correlates utr5 CDS intronic utr3 intergenic \
--gene_context 2000 \
--chromosome chromosome_1 > data/correlates/chrom/chromosome_1.txt
```

this is taking forever - trudging through the genome at 6 records
per second. a better way to do this would be to write a separate
script that creates a lookup string for each site in a chrom indicating whether
it's near a gene (i.e. 0 for no genes, 1 for upstream of a gene,
2 for downstream of a gene, 3 for 'both')

after that's done, the script could query the string positionally
to immediately assess whether the site is intergenic/upstream/downstream/both

overall - looks like this correlation script needs a bit of a rewrite today

update - here goes:

```bash
time python3.5 analysis/correlates/rcmb_correlates.py \
--table data/correlates/annotation_table_rho.txt.gz \
--windowsize 100000 \
--gene_context 2000 \
--chrom chromosome_15 \
--out data/correlates/chrom/chromosome_15.txt
```

## 4/7/2019

so after some debugging, the script is good to go and runs super fast!
chromosome 15 was done with in about 3 min and 30 seconds! 

let's queue this over the remaining chromosomes in parallel:

```bash
parallel -j 16 -i sh \
-c 'time python3.5 analysis/correlates/rcmb_correlates.py \
--table data/correlates/annotation_table_rho.txt.gz \
--windowsize 1000000 --gene_context 2000 \
--chrom chromosome_{} --out data/correlates/chrom/chromosome_{}.txt' -- {1..14} 16 17
```

and now to analyse these in `correlates_analysis.Rmd`

other correlate-related things to do once this is done:

1. hotspot enrichment by annotation
2. RR and GC at fine and broad scales

## 5/7/2019

wait, why did I not use 2 kb windows here?

```bash
parallel -j 17 -i sh \
-c 'time python3.5 analysis/correlates/rcmb_correlates.py \
--table data/correlates/annotation_table_rho.txt.gz \
--windowsize 2000 --gene_context 2000 \
--chrom chromosome_{} --out data/correlates/chrom/chromosome_{}.txt' -- {1..17}
```

alright, annotation analysis (in `correlates_analysis.Rmd`) is now done

to do:

1. GC content at broad and fine scales
2. hotspot enrichment by annotation

## 8/7/2019

today: GC content at broad and fine scales

can modify GC content script from mt locus paper (`gc_calc.py`) for this - 
now in a new script called `gc_calc_rho.py` that parses the annotation table
with an LD rho column

first pass:

```bash
mkdir -p data/correlates/gc_content

time python3.5 analysis/correlates/gc_calc_rho.py \
--filename data/references/fasta/filtered/chromosome_5.fasta \
--annotation data/correlates/annotation_table_rho.txt.gz \
--windowsize 2000 \
--chrom chromosome_5 \
--outfile data/correlates/gc_content/chromosome_5.txt
```

looks good - now for the rest, in parallel:

```bash
parallel -j 16 -i sh -c \
'time python3.5 analysis/correlates/gc_calc_rho.py \
--filename data/references/fasta/filtered/chromosome_{}.fasta \
--annotation data/correlates/annotation_table_rho.txt.gz \
--windowsize 2000 \
--chrom chromosome_{} \
--outfile data/correlates/gc_content/chromosome_{}.txt' -- {1..4} {6..17}
```

## 9/7/2019

all done - will do analyses in `gc_analysis.Rmd`

## 10/7/2019

today: hotspot enrichment

annotation counts in and out of hotspots - test run:

```bash
mkdir -p data/correlates/hotspot_enrichment

time python3.5 analysis/correlates/hotspot_annotation.py \
--filename data/ldhelmet/block_5/chromosome_15_summarised.txt \
--table data/correlates/annotation_table_rho.txt.gz \
--chrom chromosome_15 \
--out data/correlates/hotspot_enrichment/chromosome_15.txt
```

full run post-debugging:

```bash
parallel -j 17 -i sh -c \
'time python3.5 analysis/correlates/hotspot_annotation.py \
--filename data/ldhelmet/block_5/chromosome_{}_summarised.txt \
--table data/correlates/annotation_table_rho.txt.gz \
--chrom chromosome_{} \
--out data/correlates/hotspot_enrichment/chromosome_{}.txt' -- {1..17}
```

downstream analysis is in `hotspot_enrichment_analysis.Rmd`

## 6/8/2019

what is the distribution of intergenic tract lengths? how does recombination rate
vary by tract length? (and does it explain the 'both' pattern)

getting genes only:

```bash
grep 'gene' data/references/final.strict.GFF3 | cut -f 1,4,5 > data/references/genes.bed
```

getting intergenic regions w/ python:

```python
>>> from tqdm import tqdm
>>> import csv
>>> fname = 'data/references/genes.bed'
>>> outname = 'data/references/intergenic.bed'
>>> with open(fname, 'r') as f:
  2     lines = [line for line in csv.reader(f, delimiter='\t')]
>>> with open(outname, 'w', newline='') as f:
  2     f.write('\t'.join(['chrom', 'start', 'end', 'next_start', 'next_end', 'intergen_start', 'intergen_end', 'dist'
    ]))
  3     f.write('\n')
  4     w = csv.writer(f, delimiter='\t')
  5     for i in tqdm(range(len(lines) - 1)):
  6         chrom, start, end = lines[i]
  7         start, end = int(start), int(end)
  8         next_chrom, next_start, next_end = lines[i + 1]
  9         next_start, next_end = int(next_start), int(next_end)
 10         if start == next_start: # two of same record
 11             continue
 12         if chrom == next_chrom:
 13             dist = next_start - end - 1
 14             out_start = end + 1
 15             out_end = next_start - 1
 16             w.writerow([chrom, start, end, next_start, next_end, out_start, out_end, dist])
 17         else:
 18             continue
100%|███████████████████████████████████████████████████████████████████████| 21853/21853 [00:00<00:00, 190694.61it/s]
```

getting rho values and counts for each of these:

```python
>>> with open(outname, 'w', newline='') as f_out:
  2     with open(fname, 'r', newline='') as f:
  3         lines = [line for line in csv.DictReader(f, delimiter='\t')]
  4         fieldnames = list(lines[0].keys())
  5         fieldnames.extend(['rho_vals', 'rho_count', 'rho_window'])
  6         w = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
  7         for line in tqdm(lines):
  8             p = antr.Reader(table)
  9             rho_vals, rho_count = 0.0, 0
 10             for record in p.fetch(line['chrom'], int(line['intergen_start']), int(line['intergen_end'])):
 11                 if not record.ld_rho == 'NA':
 12                     rho_vals += record.ld_rho
 13                     rho_count += 1
 14             out_dict = line
 15             out_dict['rho_vals'] = rho_vals
 16             out_dict['rho_count'] = rho_count
 17             try:
 18                 out_dict['rho_window'] = rho_vals / rho_count
 19             except ZeroDivisionError:
 20                 out_dict['rho_window'] = 0.0
 21             w.writerow(out_dict)
```

nvm - need to account for overlapping records in the GFF

trying a script:

```bash
time python3.5 analysis/correlates/make_intergenic_bed.py \
--fname data/references/genes.bed \
--out data/references/intergenic.bed
```

that didn't work either - let's use the annotation table instead:

```bash
time python3.5 analysis/correlates/intergenic_tract_lengths.py \
--table data/correlates/annotation_table_rho.txt.gz \
--out data/correlates/intergenic_tract_rho.tsv
```

## 7/8/2019

script above took 1.5 hrs - but there was a bug in the script not accounting
for tracts 'spanning chromosomes' (ie when parser hits end of chrom)


## 9/8/2019

schematic:
let L = size of intergenic tract
- if L < 2 kb:
    - both = L
- if 2 kb < L < 4 kb:
    - D = L - 2 kb (ie size over 2)
        - eg 2.2 kb
    - upstream = downstream = D
        - e.g. 0.2 kb, since 0-0.2 is downstream of gene_left but >2 kb out from gene_right, and 2.0-2.2 is upstream of gene_right but >2kb out from gene_left
    - both = L - 2D
        - eg 0-0.2 downstream (gene_left), 2-2.2 upstream (gene_right), both = 1.8 kb
    - intergenic = 0
- if L > 4 kb
    - upstream = downstream = 2 kb
    - intergenic = L - 4 kb
        - ie remainder of tract that is not upstream/downstream of list
    - both = 0

## 10/8/2019

how does RR vary in UTRs by the sizes of adjacent intergenic tracts?

writing a new script that takes in the output of `intergenic_tract_lengths.py` 
and appends rho values for adjacent UTRs

## 11/8/2019

script is done:

```bash
time python3.5 analysis/correlates/utr_rho.py \
--fname data/correlates/intergenic_tract_rho.tsv \
--table data/correlates/annotation_table_rho.txt.gz \
--outname data/correlates/utr_tract_rho.tsv
```

looks good - there do seem to be some odd genes that don't have UTRs at the
starts in the GFF - looking at them manually just now it seems they just
don't have UTRs listed

```
2573it [01:24, 28.02it/s]wtf
chromosome_2 8000026
4662it [02:34, 28.30it/s]wtf
chromosome_4 2783437
6104it [03:22, 30.30it/s]wtf
chromosome_6 4293777
6701it [03:42, 28.48it/s]wtf
chromosome_6 8000002
8651it [04:47, 29.00it/s]wtf
chromosome_8 5000003
14619it [08:08, 39.42it/s]wtf
chromosome_15 1945
17279it [09:42, 29.64it/s]
```

most UTRs were parsed correctly though - let's have a look at the rho values
back in `intergenic_tract_analysis.Rmd`

## 15/8/2019

to do:
- repeat above analysis for 3' UTRs
- look at hotspots in 5' and 3' UTRs - how do they vary by tract size?

amending UTR script to also include 3' UTR values is going to be
tough if iterating through intergenic tracts - need to take start of tract
and iterate backwards (?) 

what is the longest 3' UTR in the genome?

```R
> cols <- c('chrom', 'source', 'type', 'start', 'end', 'na', 'strand', 'frame', 'id')
> d <- read_tsv('data/references/phytozome.gff', skip = 1, col_names = cols)
> utrs <- d %>% filter(type == 'three_prime_UTR')
> utrs %>% mutate(length = end - start) %>% select(length) %>% summary()
     length
 Min.   :   0.0
 1st Qu.: 400.0
 Median : 648.0
 Mean   : 771.7
 3rd Qu.: 995.0
 Max.   :2999.0
```

I don't like this solution... but I guess the site iteration could work off 
`p.fetch(chrom, utr_end - 3000, utr_end)` for now 

wait - why did the first script also include `if record.is_utr3` if the
iteration would always start at the end of a tract?? let's rerun that just
to make sure it didn't affect the above results in any way

```bash
time python3.5 analysis/correlates/utr_rho.py \
--fname data/correlates/intergenic_tract_rho.tsv \
--table data/correlates/annotation_table_rho.txt.gz \
--outname data/correlates/utr_tract_rho2.tsv
```

alright, we're good:

```bash
$ diff data/correlates/utr_tract_rho.tsv data/correlates/utr_tract_rho2.tsv | wc -l
0
```

and now for the utr3 script:

```bash
time python3.5 analysis/correlates/utr3_rho.py \
--fname data/correlates/intergenic_tract_rho.tsv \
--table data/correlates/annotation_table_rho.txt.gz \
--outname data/correlates/utr3_tract_rho.tsv
```

## 19/8/2019

`utr3_rho.py` doesn't account for there being multiple 3' UTRs for some genes...

need to make a filtered GFF with just `three_prime_UTR` lines, and then
match with intergenic tracts by proximity

could read in intergenic tract file and create lookup table

calculate `rho_vals` and `rho_count` for all UTRs, and then for each row
also have columns listing start/end for nearest intergenic tract

getting UTRs out:

```bash
grep 'three_prime_UTR' data/references/final.strict.GFF3 | cut -f 1,4,5,9 > data/references/utr3_all.bed
```

also bringing in the 9th (info) column - could use gene names to
make sure that intergenic tracts are matched correctly (ie shouldn't
have two tracts attached to the same 'chain' of successive UTRs)

first pass attempt:

```bash
time python3.5 analysis/correlates/utr3_rho.py \
--bed data/references/utr3_all.bed \
--fname data/correlates/intergenic_tract_rho.tsv \
--table data/correlates/annotation_table_rho.txt.gz \
--out data/correlates/utr3_tract_rho.tsv
```

I screwed up - got to remove scaffolds from the subsetted GFF above

```bash
grep 'three_prime_UTR' data/references/final.strict.GFF3 | \
cut -f 1,4,5,9 | \
grep -v 'scaffold' > data/references/utr3_all.bed
```

the script above still did all the chromosomes successfully though,
and took ~16 min

back to the Rmd file I go to work w/ this dataset

## 20/8/2019

so it looks like the pattern holds in 3' UTRs for the most part

I've added a new bin, 'extended', for intergenic tracts > 8 kb (since
the intergenic annotation outnumbers upstream/downstreams sites in
these tracts)

in extended tracts, rho in both UTR types is approx. equivalent;
in long tracts, five prime rho is higher than three prime,
but in medium tracts three prime rho is higher than five prime

in short tracts, both UTRs are approx equivalent, but three prime slightly higher

to do:
- compare hotspots in UTRs/tracts by tract size
- validate 2 kb definition of upstream/downstream - could it be that rho 'peaks' < 2 kb from genes

the data from the hotspot enrichment analysis are insufficient for joining w/ tract data -
will need to write a new script that parses hotspots and lists counts/other info with intergenic tracts

need to create a 'combined' 2 kb summary file for this:

```bash
cd data/ldhelmet/block_5
touch genome_summarised.txt
head -n 1 chromosome_1_summarised.txt >> genome_summarised.txt
for i in {1..17}; do
    awk 'NR%2==0' chromosome_${i}_summarised.txt >> genome_summarised.txt;
done
```

the `awk` command here prints every second line, effectively ensuring
non-overlapping 2 kb windows starting from position 0

using the block 5 dataset here since that is specifically for hotspot detection

this looks messy and will probably break, but here's a naive run:

```bash
time python3.5 analysis/correlates/intergenic_hotspots.py \
--fname data/ldhelmet/block_5/genome_summarised.txt \
--intergenic data/correlates/intergenic_tract_rho.tsv \
--out data/correlates/intergenic_tract_hotspots.tsv
```

that worked?? and only took 1.5 seconds???? oh man I love numpy

back to the Rmd I go with this

## 9/9/2019

to do:
- write script that pulls rho for 'edges' of intergenic tracts
    - script should take in a 'gene proximate size' as input

here goes, with 2 kb:

```bash
time python3.5 analysis/correlates/intergenic_tract_proximal.py \
--fname data/correlates/intergenic_tract_rho.tsv \
--table data/correlates/annotation_table_rho.txt.gz \
--windowsize 2000 \
--outfile data/correlates/intergenic_flanks_2kb.tsv
```

script is working - projected to take ~10 min to complete, which isn't awful

also running for 500 bp and 1 kb

```bash
time python3.5 analysis/correlates/intergenic_tract_proximal.py \
--fname data/correlates/intergenic_tract_rho.tsv \
--table data/correlates/annotation_table_rho.txt.gz \
--windowsize 1000 \
--outfile data/correlates/intergenic_flanks_1kb.tsv

time python3.5 analysis/correlates/intergenic_tract_proximal.py \
--fname data/correlates/intergenic_tract_rho.tsv \
--table data/correlates/annotation_table_rho.txt.gz \
--windowsize 500 \
--outfile data/correlates/intergenic_flanks_500.tsv
```



