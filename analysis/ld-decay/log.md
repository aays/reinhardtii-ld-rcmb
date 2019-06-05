
log for LD decay analysis

## 4/6/2019

getting plink 1.9 into the proj dir:

```bash
ln -sv ~/apps/plink_linux_x86_64/plink bin/
```

prepping a sample chromosome for starters:

```bash
mkdir -p data/ld-decay/chromosome_1
mkdir -p data/ld-decay/finals
time ./bin/plink \
--vcf data/references/vcf/filtered/chromosome_1.vcf.gz \
--make-bed \
--out data/ld-decay/chromosome_1/chromosome_1 \
--allow-extra-chr
```

r2 calc:

```bash
time ./bin/plink \
--bfile data/ld-decay/chromosome_1/chromosome_1 \
--r2 \
--out data/ld-decay/finals/chromosome_1 \
--ld-window-r2 0 \
--ld-window-kb 100 \
--allow-extra-chr
```

where `ld-window-r2` means that anything with r2 < 0.2 isn't filtered out
(which is the default setting in plink)

looks good - though will need a script to parse that output - it doesn't
seem to use regular delimiters

running over other chromosomes:

```bash
# file prep
for i in {2..17}; do
    mkdir -p data/ld-decay/chromosome_${i};
    time ./bin/plink \
    --vcf data/references/vcf/filtered/chromosome_${i}.vcf.gz \
    --make-bed \
    --out data/ld-decay/chromosome_${i}/chromosome_${i} \
    --allow-extra-chr ;
done

# r2 calc
for i in {2..17}; do
    time ./bin/plink \
    --bfile data/ld-decay/chromosome_${i}/chromosome_${i} \
    --r2 \
    --out data/ld-decay/finals/chromosome_${i} \
    --ld-window-r2 0 \
    --ld-window-kb 100 \
    --allow-extra-chr;
done
```

removing unneeded files:

```bash
rm -v data/ld-decay/finals/*nosex
rm -v data/ld-decay/finals/*log
```

cleaning up the plink outputs in R:

```R
library(readr)
library(dplyr)
library(purrr)
library(fs)
library(stringr)

fnames <- dir_ls(regexp = 'chromosome_[0-9]{1,2}\\.ld')

write_ld_file <- function(fname) {
    outname <- str_extract(fname, 'chromosome_[0-9]{1,2}') %>%
        paste0(., '.txt')
    d <- read_table(fname, col_types = cols()) %>%
             select(chr = CHR_A, snp1 = SNP_A, snp2 = SNP_B, r2 = R2)
    write_csv(d, outname)
}

fnames %>%
    walk(~ write_ld_file(.))
```

## 5/6/2019

today: fit W-H equation to LD estimates from plink

wait - the script from yesterday didn't work - doing this again

```R
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

```
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























