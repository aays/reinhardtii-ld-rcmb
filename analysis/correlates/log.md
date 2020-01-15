
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

## 12/9/2019

100 bp flanks:

```bash
time python3.5 analysis/correlates/intergenic_tract_proximal.py \
--fname data/correlates/intergenic_tract_rho.tsv \
--table data/correlates/annotation_table_rho.txt.gz \
--windowsize 100 \
--outfile data/correlates/intergenic_flanks_100.tsv
```

## 29/11/2019

post-review - need to account for SNP density in comparing
rho between annotations (ie CDS vs intron, genic vs intergenic)

re-orienting myself - got chromosomal VCFs in `data/references/vcf/filtered`
and LDhelmet 2 kb window files in `data/ldhelmet/block_100`

need to write a script that takes in a vcf + summarised ldhelmet outfile
as input and computes SNP density for each window before adding it on as a separate column

first pass:

```bash
time python3.5 analysis/correlates/snp_density.py \
--filename data/ldhelmet/block_100/chromosome_15_summarised.txt \
--vcf data/references/vcf/filtered/chromosome_15.vcf.gz \
--chrom chromosome_15 \
--outfile data/ldhelmet/block_100/chromosome_15_density.txt
```

looks good - now over the genome:

```bash
for i in {1..14} 16 17; do
    time python3.5 analysis/correlates/snp_density.py \
    --filename data/ldhelmet/block_100/chromosome_${i}_summarised.txt \
    --vcf data/references/vcf/filtered/chromosome_${i}.vcf.gz \
    --chrom chromosome_${i} \
    --outfile data/ldhelmet/block_100/chromosome_${i}_density.txt;
done
```

## 4/12/2019


alright, now to fit a model with annotation and polymorphism (ie `rcmb_rate ~
percent_CDS + polymorphism`)

though I think this will be most easily done if the above script is further
modified to add in # intergenic/genic/CDS/intron sites as separate columns

the count column lookup is as follows:

```
i: intergenic
c: cds
n: intron
5: utr5
3: utr3
N: unknown
```

running the script on chr 15 as a test:

```bash
time python3.5 analysis/correlates/snp_density.py \
--filename data/ldhelmet/block_100/chromosome_15_summarised.txt \
--vcf data/references/vcf/filtered/chromosome_15.vcf.gz \
--chrom chromosome_15 \
--table data/correlates/annotation_table_rho.txt.gz \
--outfile data/ldhelmet/block_100/chromosome_15_density.txt
```

this script also now saves the lookup as `chromosome_i_temp_lookup` so that it
can be provided on subsequent runs and doesn't have to be recreated every single time

running on all other chromosomes (running in parallel bugged out for some reason...)

```bash
time for i in {1..14} 16 17; do
    time python3.5 analysis/correlates/snp_density.py \
    --filename data/ldhelmet/block_100/chromosome_${i}_summarised.txt \
    --vcf data/references/vcf/filtered/chromosome_${i}.vcf.gz \
    --chrom chromosome_${i} \
    --table data/correlates/annotation_table_rho.txt.gz \
    --outfile data/ldhelmet/block_100/chromosome_${i}_density.txt;
done
```

## 5/12/2019

loop above took 1h 30min, but looks like it was successful

having a look at this in `correlates_analysis.Rmd`

update - I don't think this approach is the way to go - it offers some insight
into genic vs intergenic per 2 kb window, but we're working with rho averages/window

a better thing to do would be to use the lookup strings in `data/annotation_lookup`
and create a dataset that shows 'tracts' of annotations, looking something like this:

```
chr start end rho_sum rho_count rho_tract snp_count is_intergenic is_CDS is_intronic ...
chromosome_1 1000 4500 7 3500 0.002 60 1 0 0 ...
```

bc then that way we can fit `rho_window ~ snp_density + is_CDS` with `is_CDS` encoded as
a categorical variable, and then see how much `is_CDS:1` affects `rho_window`

## 8/12/2019

renaming all the lookups that were created by `snp_density.py`

```bash
for i in {1..17}; do
    mv -v data/annotation_lookups/chromosome_${i}_temp_lookup data/annotation_lookups/chromosome_${i}.txt;
done
```

alright, first pass at script is done - here goes:

```bash
time python3.5 analysis/correlates/annotation_tracts.py \
--filename data/annotation_lookups/chromosome_15.txt \
--table data/correlates/annotation_table_rho.txt.gz \
--vcf data/references/vcf/filtered/chromosome_15.vcf.gz \
--chrom chromosome_15 \
--out data/correlates/annotation_tracts/chromosome_15.txt
```

worked cleanly and was done in two minutes! sometimes I write good code! 

```bash
for i in {1..14} 16 17; do
    time python3.5 analysis/correlates/annotation_tracts.py \
    --filename data/annotation_lookups/chromosome_${i}.txt \
    --table data/correlates/annotation_table_rho.txt.gz \
    --vcf data/references/vcf/filtered/chromosome_${i}.vcf.gz \
    --chrom chromosome_${i} \
    --out data/correlates/annotation_tracts/chromosome_${i}.txt;
    sleep 1;
done
```

## 13/12/2019

redoing the bootstrapping stuff on the server (since it's going to take a while)

creating a new hardcoded script called `intergenic_tract_boot.R` to do
this which lifts from `intergenic_tract_analysis.Rmd`

## 16/12/2019

ended up doing this in the console bc the script broke
due to some bizarre rounding error 

took a while:

```R
> utr3_cis <- c('short', 'long') %>%
+   map_dfr(~ utr_boot_fxn(
+     rename(
+       utr3_tracts, utr_rho_vals = utr3_rho_vals, utr_rho_count = utr3_rho_coun+     20000, .)) %>%
+   mutate(utr = 'three_prime')
Starting boot
2019-12-15 17:58:35
Generating CIs for short
Completed sampling for short
Done
2019-12-15 20:57:07
Starting boot
2019-12-15 20:57:07
Generating CIs for long
Completed sampling for long
Done
2019-12-15 21:00:05
> utr_cis <- bind_rows(utr5_cis, utr3_cis)
> utr_cis
     conf.low   conf.high   bin         utr
1 0.002758457 0.003667423 short  five_prime
2 0.002645785 0.004273968  long  five_prime
3 0.003035521 0.004493206 short three_prime
4 0.002821646 0.004234891  long three_prime
> head(utr_bars_all)
# A tibble: 4 x 9
  name  bin   total_rho total_count sd_rho     n mean_rho   se_rho utr
  <chr> <chr>     <dbl>       <dbl>  <dbl> <int>    <dbl>    <dbl> <chr>
1 1     long      2880.      894459 0.0177  1474  0.00322 0.000462 five_prime
2 1     short    22384.     7292922 0.0130 14792  0.00307 0.000107 five_prime
3 2     long      5317.     1576256 0.0215  1920  0.00337 0.000491 three_prime
4 2     short    45869.    13370248 0.0227 17282  0.00343 0.000173 three_prime
> utr_bars_cis <- utr_bars_all %>% left_join(utr_cis, by = c('utr', 'bin'))
> write_csv(utr_bars_cis, here('data/correlates/utr_bars_cis_corrected.csv'))
> write_csv(utr_cis, here('data/correlates/utr_cis_only.csv'))
```

the utr cis are now in `data/correlates`

repeating this for intergenic tracts

## 24/12/2019

need to redo the intergenic tract analysis, but with the same
multiple regression approach used for the other correlate stuff

first need to make an all chromosomes VCF (since the existing
script works on the whole genome)

```python
# in data/references/vcf/filtered
from cyvcf2 import VCF, Writer
from tqdm import tqdm

fnames = ['chromosome_{}.vcf.gz'.format(i) for i in range(1, 18)]

vcf_in = VCF(fnames[0])
vcf_out = Writer('all_chrom.vcf', vcf_in)

for fname in fnames:
    print(fname)
    for record in tqdm(VCF(fname)):
        vcf_out.write_record(record)
```

actually - what if I just used python to directly
update the outfile from `intergenic_tract_proximal.py`?

```python
>>> tract_fname = 'data/correlates/intergenic_flanks_2kb.tsv'
>>> from tqdm import tqdm
>>> from cyvcf2 import VCF
>>> import csv
>>> with open(tract_fname, 'r', newline='') as f:
  2     tract_lines = [line for line in csv.DictReader(f, delimiter='\t')]
>>> tract_lines[0].keys()
dict_keys(['left_window', 'start', 'right_vals', 'right_count', 
'tract_size', 'chrom', 'windowsize', 'end', 'right_window', 
'left_count', 'left_vals'])
>>> from copy import deepcopy
  2 with open('data/correlates/intergenic_flanks_2kb_snps.tsv', 'w') as f:
  3     fieldnames = list(tract_lines[0].keys())
  4     fieldnames.extend(['snp_count', 'snp_density'])
  5     writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter = '\t')
  6     writer.writeheader()
  7     vcf_fname = 'data/references/vcf/filtered/all_chrom.vcf.gz'
  8     v = VCF(vcf_fname)
  9     for line in tqdm(tract_lines):
 10         chrom, start, end = line['chrom'], int(line['start']), int(line['
    end'])
 11         region = '{c}:{s}-{e}'.format(c=chrom, s=start, e=end)
 12         snp_count = len([r for r in v.__call__(region)])
 13         out_dict = deepcopy(line)
 14         out_dict['snp_count'] = snp_count
 15         out_dict['snp_density'] = snp_count / (end - start)
 16         writer.writerow(out_dict)
100%|█████████████████████████████████| 16192/16192 [00:31<00:00, 510.29it/s]
```

in Rmd - set bin size as a binary categorical predictor
and regress RR ~ bin + diversity

## 26/12/2019

so this column is useful for full tract comparisons, but there also needs to be
a separate column for SNP density specifically in the 'left window' and 'right window'
of longer (>2kb) tracts

for shorter tracts, report SNP density for full tract - same as before

```python
>>> tract_fname = 'data/correlates/intergenic_flanks_2kb_snps.tsv'
>>> from tqdm import tqdm
>>> from cyvcf2 import VCF
>>> import csv
>>> with open(tract_fname, 'r', newline='') as f:
  2     tract_lines = [line for line in csv.DictReader(f, delimiter='\t')]
>>> fieldnames = ['chrom', 'start', 'end', 'tract_size', 'windowsize', 'left_vals',
  2 'left_count', 'left_window', 'right_vals', 'right_count', 'right_window', 'snp_count', 'snp_density']
>>> from copy import deepcopy
>>> with open('data/correlates/intergenic_flanks_2kb_tract_snps.tsv', 'w') as f:
  2     fieldnames.extend(['left_snp_count', 'left_snp_density', 'right_snp_count', 'right_snp_density'])
  3     writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
  4     writer.writeheader()
  5     vcf_fname = 'data/references/vcf/filtered/all_chrom.vcf.gz'
  6     v = VCF(vcf_fname)
  7     for line in tqdm(tract_lines):
  8         chrom, start, end, tract_size = line['chrom'], int(line['start']), int(line['end']), int(line['tract_size
    '])
  9         if tract_size > 4000:
 10             left_region = '{}:{}-{}'.format(chrom, start, start + 2000)
 11             right_region = '{}:{}-{}'.format(chrom, end - 2000, end)
 12         elif tract_size <= 4000:
 13             half_tract_size = round(tract_size / 2)
 14             left_region = '{}:{}-{}'.format(chrom, start, start + half_tract_size)
 15             right_region = '{}:{}-{}'.format(chrom, start + half_tract_size, end)
 16         left_snp_count = len([r for r in v.__call__(left_region)])
 17         right_snp_count = len([r for r in v.__call__(right_region)])
 18         out_dict = deepcopy(line)
 19         out_dict['left_snp_count'], out_dict['right_snp_count'] = left_snp_count, right_snp_count
 20         out_dict['left_snp_density'] = left_snp_count / 2000
 21         out_dict['right_snp_density'] = right_snp_count / 2000
 22         writer.writerow(out_dict)
100%|█████████████████████████████████████████████████████████████████████████| 16192/16192 [00:56<00:00, 285.38it/s]
```

## 28/12/2019

updating the intergenic proximal script to get full 2 kb to the 'right' in < 4 kb tracts
instead of splitting tracts in half

```bash
time python3.5 analysis/correlates/intergenic_tract_proximal.py \
--fname data/correlates/intergenic_tract_rho.tsv \
--table data/correlates/annotation_table_rho.txt.gz \
--windowsize 2000 \
--split right \
--outfile data/correlates/intergenic_flanks_2kb_right.tsv
```

done in 10 min like before

adding SNPs using adapted python code above:

(n.b. lines 20 and 21 above are wrong - count should only be divided by 2000 if
tract length > 4 kb - not going to be using that output file anyways though)

```python
>>> tract_fname = 'data/correlates/intergenic_flanks_2kb_right.tsv'
>>> from tqdm import tqdm; from cyvcf2 import VCF; import csv
>>> from copy import deepcopy
>>> with open(tract_fname, 'r', newline='') as f:
  2     tract_lines = [line for line in csv.DictReader(f, delimiter='\t')]
>>> fieldnames = ['chrom', 'start', 'end', 'tract_size', 'windowsize', 'left_vals', 'left_count',
  2 'left_window', 'right_vals', 'right_count', 'right_window']
>>> with open('data/correlates/intergenic_flanks_2kb_right_snps.tsv', 'w') as f:
  2     fieldnames.extend(['snp_count', 'left_snp_count', 'left_snp_density',
  3                        'right_snp_count', 'right_snp_density'])
  4     writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
  5     writer.writeheader()
  6     vcf_fname = 'data/references/vcf/filtered/all_chrom.vcf.gz'
  7     v = VCF(vcf_fname)
  8     for line in tqdm(tract_lines):
  9         chrom, start, end = line['chrom'], int(line['start']), int(line['end'])
 10         tract_size = int(line['tract_size'])
 11         region_str = '{}:{}-{}'
 12         if tract_size > 4000:
 13             left_len, right_len = 2000, 2000
 14             left_region = region_str.format(chrom, start, start + 2000)
 15             right_region = region_str.format(chrom, end - 2000, end)
 16         elif 2000 < tract_size <= 4000:
 17             left_len = (end - 2000) - start
 18             right_len = 2000
 19             right_region = region_str.format(chrom, end - 2000, end)
 20             left_region = region_str.format(chrom, start, end - 2000)
 21         elif tract_size <= 2000:
 22             half_tract_size = round(tract_size / 2)
 23             left_len, right_len = half_tract_size, half_tract_size
 24             left_region = region_str.format(chrom, start, start + half_tract_size)
 25             right_region = region_str.format(chrom, start + half_tract_size, end)
 26         overall_region = region_str.format(chrom, start, end)
 27         total_snp_count = len([r for r in v.__call__(overall_region)])
 28         left_snp_count = len([r for r in v.__call__(left_region)])
 29         right_snp_count = len([r for r in v.__call__(right_region)])
 30         out_dict = deepcopy(line)
 31         out_dict['snp_count'] = total_snp_count
 32         out_dict['left_snp_count'], out_dict['right_snp_count'] = left_snp_count, right_snp_count
 33         out_dict['left_snp_density'] = left_snp_count / left_len
 34         out_dict['right_snp_density'] = right_snp_count / right_len
 35         writer.writerow(out_dict)
100%|█████████████████████████████████████████████████████████████████████████| 16192/16192 [01:27<00:00, 185.63it/s]
```

now to have a look at this in the Rmd

update - this looks good - going to repeat for 'left' 2 kb tracts (ie downstream of genes)
for completeness

```bash
time python3.5 analysis/correlates/intergenic_tract_proximal.py \
--fname data/correlates/intergenic_tract_rho.tsv \
--table data/correlates/annotation_table_rho.txt.gz \
--windowsize 2000 \
--split left \
--outfile data/correlates/intergenic_flanks_2kb_left.tsv
```

adding SNPs, as before (I honestly should have turned this into a script long ago...)

```python
>>> tract_fname = 'data/correlates/intergenic_flanks_2kb_left.tsv'
>>> from tqdm import tqdm; from cyvcf2 import VCF; import csv; from copy import deepcopy
>>> with open(tract_fname, 'r', newline='') as f:
  2     tract_lines = [line for line in csv.DictReader(f, delimiter='\t')]
>>> fieldnames = ['chrom', 'start', 'end', 'tract_size', 'windowsize', 'left_vals', 'left_count',
  2 'left_window', 'right_vals', 'right_count', 'right_window']
  3 with open('data/correlates/intergenic_flanks_2kb_left_snps.tsv', 'w') as f:
  4     fieldnames.extend(['snp_count', 'left_snp_count', 'left_snp_density',
  5                        'right_snp_count', 'right_snp_density'])
  6     writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
  7     writer.writeheader()
  8     vcf_fname = 'data/references/vcf/filtered/all_chrom.vcf.gz'
  9     v = VCF(vcf_fname)
 10     for line in tqdm(tract_lines):
 11         chrom, start, end = line['chrom'], int(line['start']), int(line['end'])
 12         tract_size = int(line['tract_size'])
 13         region_str = '{}:{}-{}'
 14         if tract_size > 4000:
 15             left_len, right_len = 2000, 2000
 16             left_region = region_str.format(chrom, start, start + 2000)
 17             right_region = region_str.format(chrom, end - 2000, end)
 18         elif 2000 < tract_size <= 4000:
 19             left_len = 2000
 20             right_len = tract_size - 2000
 21             left_region = region_str.format(chrom, start, start + 2000)
 22             right_region = region_str.format(chrom, start + 2000, end)
 23         elif tract_size <= 2000:
 24             half_tract_size = round(tract_size / 2)
 25             left_len, right_len = half_tract_size, half_tract_size
 26             left_region = region_str.format(chrom, start, start + half_tract_size)
 27             right_region = region_str.format(chrom, start + half_tract_size, end)
 28         overall_region = region_str.format(chrom, start, end)
 29         total_snp_count = len([r for r in v.__call__(overall_region)])
 30         left_snp_count = len([r for r in v.__call__(left_region)])
 31         right_snp_count = len([r for r in v.__call__(right_region)])
 32         out_dict = deepcopy(line)
 33         out_dict['snp_count'] = total_snp_count
 34         out_dict['left_snp_count'], out_dict['right_snp_count'] = left_snp_count, right_snp_count
 35         out_dict['left_snp_density'] = left_snp_count / left_len
 36         out_dict['right_snp_density'] = right_snp_count / right_len
 37         writer.writerow(out_dict)
100%|█████████████████████████████████████████████████████████████████████████| 16192/16192 [01:21<00:00, 199.17it/s]
```

## 15/1/2020

redoing hotspot enrichment w/ the new hotspot dataset

moved the older files into `data/correlates/hotspot_enrichment/old`

```bash
parallel -j 17 -i sh -c \
'time python3.5 analysis/correlates/hotspot_annotation.py \
--filename data/ldhelmet/block_10/chromosome_{}_summarised.txt \
--table data/correlates/annotation_table_rho.txt.gz \
--chrom chromosome_{} \
--out data/correlates/hotspot_enrichment/chromosome_{}.txt' -- {1..17}
```

took 8.5 min - nice - off we go to rerun `hotspot_enrichment_analysis.Rmd` with
these new files


































