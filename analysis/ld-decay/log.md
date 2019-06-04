
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

fnames <- dir_ls(pattern = '*ld')

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




















