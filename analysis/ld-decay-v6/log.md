
log for LD decay analysis, with new VCF

## 09/01/2021

grabbing VCF:

```bash
cd data/ld-decay-2021
ln -sv /scratch/research/projects/chlamydomonas/Cincerta_deNovo/analysis/assembly_V3/chlamy_v6_project/QC17_CC4532/SNP_filtering/QC17_CC4532.all_sites.all_isolates.vcf.gz* .
mv -v QC17_CC4532.all_sites.all_isolates.vcf.gz v6.vcf.gz
mv -v QC17_CC4532.all_sites.all_isolates.vcf.gz.tbi v6.vcf.gz.tbi
mkdir vcf
```

splitting into chromosomes:

```python
from cyvcf2 import VCF
from cyvcf2 import Writer
from tqdm import tqdm

for i in tqdm(range(1, 18)):
    if i < 10:
        chrom = 'chromosome_0{}'.format(i)
    elif i >= 10:
        chrom = 'chromosome_{}'.format(i)
    reader = VCF('v6.vcf.gz')
    writer = Writer('vcf/' + chrom + '.vcf', reader)
    print('starting {}'.format(chrom))
    writer.write_header()
    for record in tqdm(reader(chrom)):
        writer.write_record(record)
    print('done {}'.format(chrom))
```

trying out plink on the first chr:

```bash
mkdir -p data/ld-decay-v6/chromosome_1
mkdir -p data/ld-decay-v6/finals
time ./bin/plink \
--vcf data/ld-decay-v6/vcf/chromosome_01.vcf \
--make-bed \
--out data/ld-decay-v6/chromosome_1/chromosome_1 \
--allow-extra-chr

# r2 calc
time ./bin plink \
--bfile data/ld-decay-v6/chromosome_1/chromosome_1 \
--r2 --out data/ld-decay-v6/finals/chromosome_1 \
--ld-window-r2 0 --ld-window-kb 100 \
--allow-extra-chr
```

looks good - scaling up:

```bash
time for i in {2..9}; do
    mkdir -p data/ld-decay-v6/chromosome_${i}
    time ./bin/plink \
    --vcf data/ld-decay-v6/vcf/chromosome_0${i}.vcf \
    --make-bed \
    --out data/ld-decay-v6/chromosome_${i}/chromosome_${i} \
    --allow-extra-chr
done

time for i in {2..9}; do
    time ./bin/plink \
    --bfile data/ld-decay-v6/chromosome_${i}/chromosome_${i} \
    --r2 --out data/ld-decay-v6/finals/chromosome_${i} \
    --ld-window-r2 0 \
    --ld-window-kb 100 \
    --allow-extra-chr;
done

# remaining 7
time for i in {10..17}; do
    mkdir -p data/ld-decay-v6/chromosome_${i}
    time ./bin/plink \
    --vcf data/ld-decay-v6/vcf/chromosome_${i}.vcf \
    --make-bed \
    --out data/ld-decay-v6/chromosome_${i}/chromosome_${i} \
    --allow-extra-chr
done

time for i in {10..17}; do
    time ./bin/plink \
    --bfile data/ld-decay-v6/chromosome_${i}/chromosome_${i} \
    --r2 --out data/ld-decay-v6/finals/chromosome_${i} \
    --ld-window-r2 0 \
    --ld-window-kb 100 \
    --allow-extra-chr;
done

```
removing unneeded files:

```bash
rm -v data/ld-decay-v6/finals/*nosex
rm -v data/ld-decay-v6/finals/*log
```


## 12/1/2021

today: actually fitting the models:

```R
# in data/ld-decay-v6/finals dir

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(fs)
library(stringr)
library(janitor)

fnames <- dir_ls(regexp = 'chromosome_[0-9]{1,2}\\.ld')

write_ld_file <- function(fname) {
    chrname <- fname %>%
        str_extract('chromosome_[0-9]{1,2}')
    outname <- paste0(chrname, '.csv')
    d <- read_delim(fname, delim = ' ', col_types = cols()) %>%
        clean_names() %>%
        select(bp_a, bp_b, r2) %>%
        mutate_at(
            vars(contains('bp')),
            .funs = list(num = ~ str_extract(., '[0-9]+'))
        ) %>%
        select(-bp_a, -bp_b) %>%
        select(pos1 = bp_a_num, pos2 = bp_b_num, r2) %>%
        mutate_all(as.numeric) %>%
        mutate(chrom = chrname) %>%
        select(chrom, everything())
    write_csv(d, outname)
}

fnames %>%
    walk(~ write_ld_file(.))
```

looks good - now for fits:

```R
library(readr)
library(dplyr)
library(stringr)
library(magrittr)
library(fs)
library(purrr)

fnames <- dir_ls(regexp = 'chromosome_[0-9]{1,2}\\.csv')

weir_hill <- function(fname) {

    chrname <- str_extract(fname, 'chromosome_[0-9]{1,2}')
    outname <- paste0(chrname, '_fit.csv')
    equation <- '((10 + p*d)/(22 + (13*p*d) + (p*d)^2))*(1 + (((3 + (p*d))/(24*(22 + (13*p*d) +\
                   (p*d)^2))) * (12 + (12*p*d) + (p*d)^2)))' %>%
                   str_replace(., '\\n', '') %>%
                   str_replace(' {2,}', '')

    d <- read_csv(fname, col_types = cols()) %>%
        mutate(d = abs(pos2 - pos1)) %>%
        select(d, r2)
    predicted_decay <- nls(
            paste('r2', '~', equation),
            data = d, control = list(maxiter = 500), start = list(p = 0.5)
        ) %>%
        predict() %>%
        as_tibble()
    colnames(predicted_decay) <- 'r2'
    predicted_decay$d <- d$d
    
    predicted_decay %<>%
        group_by(d) %>%
        summarise(r2 = mean(r2)) %>%
        arrange(d) %>%
        mutate(chrom = chrname) %>%
        select(chrom, d, r2)

    write_csv(predicted_decay, outname)
}

fnames %>%
    walk(~ weir_hill(.))
```

and finally, getting the fitted rho values:

```R
library(readr)
library(dplyr)
library(stringr)
library(magrittr)
library(fs)
library(purrr)
library(broom)

fnames <- dir_ls(regexp = 'chromosome_[0-9]{1,2}\\.csv')

weir_hill_rho <- function(fname) {

    chrname <- str_extract(fname, 'chromosome_[0-9]{1,2}')
    outname <- paste0(chrname, '_fit.csv')
    equation <- '((10 + p*d)/(22 + (13*p*d) + (p*d)^2))*(1 + (((3 + (p*d))/(24*(22 + (13*p*d) +\
                   (p*d)^2))) * (12 + (12*p*d) + (p*d)^2)))' %>%
                   str_replace(., '\\n', '') %>%
                   str_replace(' {2,}', '')

    d <- read_csv(fname, col_types = cols()) %>%
        mutate(d = abs(pos2 - pos1)) %>%
        select(d, r2)

    fit <- nls(
            paste('r2', '~', equation),
            data = d, control = list(maxiter = 500), start = list(p = 0.5)) %>%
            tidy() %>%
            mutate(chrom = chrname) %>%
            select(chrom, everything())

    message(chrname)
    return(fit)

}

d_final <- fnames %>%
    map_dfr(~ weir_hill_rho(.))

write_csv(d_final, 'fitted_rho.csv')
```

next up - analyse + plot decay stats over chrs - doing this in
an Rmd file (`decay_analysis.Rmd`)






















