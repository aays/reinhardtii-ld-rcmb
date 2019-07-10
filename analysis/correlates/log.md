
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














