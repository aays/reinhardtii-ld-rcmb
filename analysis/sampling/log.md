
## 12/5/2019

to do: get variant/sequencing stats

1. initial # of SNPs
2. mean coverage
3. SNPs left after filtering 
    - GQ > 30
    - MAF > 0.1 to excl. singletons
    - diallelic

first, porting over references from prev folder:

```bash
# vcfs
cd data/references/vcf
mkdir raw filtered
cd raw
ln -sv ../../../../../data/vcfs/no_clones/*vcf* .

# fastas - back to root
cd data/references/fasta
ln -sv ../../../../data/fastas/no_clones/* .
```

the `vcf_filter.py` script from pyVCF should
allow for most of the filtering needed

counting initial SNPs:

```python
>>> from cyvcf2 import VCF
>>> fname = 'raw/chromosome_1.vcf.gz'
>>> v = VCF(fname)
>>> counter = 0
>>> from tqdm import tqdm
  2 for record in tqdm(v):
  3     if record.is_snp:
  4         counter += 1
672564it [00:09, 74396.75it/s]
>>> counter
480144
```

doing this for the remaining chromosomes:

```python
>>> fname = 'raw/chromosome_{}.vcf.gz'
>>> for i in tqdm(range(2, 18)):
  2     current_chrom = fname.format(i)
  3     for record in VCF(current_chrom):
  4         if record.is_snp:
  5             counter += 1
100%|█████████████████████████████████████████| 16/16 [01:48<00:00,  6.69s/it]
>>> counter
6497950
```

6,497,950 initial SNPs prior to filtering

getting the `vcf_filter` script (in `analysis/sampling`)

```
wget https://raw.githubusercontent.com/jamescasbon/PyVCF/master/scripts/vcf_filter.py
```

first pass run on `chromosome_1.vcf.gz`:

```bash
time python3.5 analysis/sampling/vcf_filter.py \
data/references/vcf/raw/chromosome_1.vcf.gz
--mgq 30 \
--snp-only > \
data/references/vcf/filtered/chromosome_1.vcf
```

alright - this isn't working and this script's
documentation is confusing me so let's try `bcftools` instead

## 13/5/2019

getting the latest version (1.9):

```bash
cd ~/apps
git clone git://github.com/samtools/htslib.git
git clone git://github.com/samtools/bcftools.git
cd bcftools
make
```

here's a sample bcftools filter command:

```bash
./bin/bcftools filter -i 'AVG(GQ)>30 && TYPE="snp" && MAF>0.1' data/references/vcf/raw/chromosome_1.vcf.gz \
--output test.vcf
```

this already removes many variants:

```bash
zgrep -v '#' data/references/vcf/raw/chromosome_1.vcf.gz | wc -l
# 672564

grep -v '#' test.vcf | wc -l
# 398053
```

running this over all chromosomes (although a python script will be needed to 
keep just the diallelic SNPs)

```bash
mkdir -p data/references/vcf/filtered/temp
for i in {1..17}; do
    echo "currently on chromosome ${i}" ;
    ./bin/bcftools filter -i 'AVG(GQ)>30 && TYPE="snp" && MAF>0.1' \
    data/references/vcf/raw/chromosome_${i}.vcf.gz \
    --output data/references/vcf/filtered/temp/chromosome_${i}.vcf ;
done
```

quick count comparison eye test

```
for i in {1..17}; do
    echo "chromosome_${i}" ;
    zgrep -v '#' data/references/vcf/raw/chromosome_${i}.vcf.gz | wc -l ;
    grep -v '#' data/references/vcf/filtered/chromosome_${i}.vcf | wc -l ;
    printf "\n" ;
done
```

all of these have had about ~200k variants removed each - no anomalies

python script to just keep diallelic SNPs:

```python
import sys
from cyvcf2 import VCF, Writer
from tqdm import tqdm

fname = sys.argv[-2]
outname = sys.argv[-1]

vcf_in = VCF(fname)
vcf_out = Writer(outname, vcf_in)

kept_count = 0
total_count = 0

print('Reading from {}'.format(fname))
print('Writing to {}'.format(outname))

for record in tqdm(vcf_in):
    total_count += 1
    if len(record.ALT) == 1:
        vcf_out.write_record(record)
        kept_count += 1
    else:
        continue

print('{0} records kept out of {1} total.'.format(kept_count, total_count))
```
    
quick test on chromosome 1:

```bash
time python3.5 analysis/sampling/filter_diallelic.py \
data/references/vcf/filtered/temp/chromosome_1.vcf \
data/references/vcf/filtered/chromosome_1.vcf
```

looks good - doing this across the remainder:

```bash
for i in {2..17}; do
    time python3.5 analysis/sampling/filter_diallelic.py \
    data/references/vcf/filtered/temp/chromosome_${i}.vcf \
    data/references/vcf/filtered/chromosome_${i}.vcf ;
done
```

how many SNPs does this leave us with in total?

```python
>>> from cyvcf2 import VCF
>>> from tqdm import tqdm
>>> fname = 'chromosome_{}.vcf'
>>> counter = 0
  2 for i in tqdm(range(1, 18)):
  3     current_chrom = fname.format(i)
  4     for record in VCF(current_chrom):
  5         assert record.is_snp and len(record.ALT) == 1
  6         counter += 1
100%|█████████████████████████████████████████████████████████████████████████████████| 17/17 [00:49<00:00,  2.83s/it]
>>> counter
4736814
```

finally, bgzip and tabix references:

```bash
for i in {1..17}; do
    echo "chromosome_${i}";
    bgzip chromosome_${i}.vcf ;
    tabix chromosome_${i}.vcf.gz ;
done
```






