log for LDhelmet analysis

## 31/5/2019

goal - two LDhelmet runs over C. reinhardtii genome (block = 5 and block = 100)

first, need to generate fastas out of filtered VCFs

creating script symlinks + symlink to reference fasta:

```bash
ln -sv /scratch/research/repos/vcf2fasta/vcf2fasta.py bin/vcf2fasta.py

cd data/references/fasta
ln -sv /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/chlamy.5.3.w_organelles_mtMinus.fa .

cd ../../
```

python script with subprocess to generate fastas:

```python
import subprocess
from tqdm import tqdm
import time

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

vcf2fasta_cmd = 'time ./bin/vcf2fasta.py -v data/references/vcf/filtered/{chrom}.vcf.gz \
-r data/references/fasta/chlamy.5.3.w_organelles_mtMinus.fasta \
-i {chrom}:1-{length} \
-s CC2935 CC2936 CC2937 CC2938 CC3059 CC3060 CC3061 CC3062 CC3063 CC3064 \
CC3065 CC3068 CC3071 CC3073 CC3075 CC3076 CC3079 CC3084 CC3086 \
--min_GQ 30 > data/references/fasta/filtered/{chrom}.fasta'

for chrom in tqdm(sorted(lengths.keys())):
    chrom_length = lengths[chrom]
    subprocess.call(vcf2fasta_cmd.format(chrom=chrom, length=chrom_length), shell=True)
    print('done {chrom}'.format(chrom=chrom))
    time.sleep(3)

```

all done - took about 2 hours (should really look into speeding up this script if possible
using the long format fasta conversion)


## 1/6/2019

LDhelmet runs! 

prepping dirs

```bash
mkdir -p data/ldhelmet
mkdir -p data/ldhelmet/block_5
mkdir -p data/ldhelmet/block_100
cp -v ../ldhelmet-sims/data/mut_mat data/ldhelmet/ # get mutation matrix
```

LDhelmet run script (`run_ldhelmet.sh`)

```bash
block=$1
outdir=block_$1

for fname in data/references/fasta/filtered/*fasta; do
    base=$(basename $fname .fasta)
    echo "Currently on ${base}"

    time ./bin/ldhelmet find_confs \
    --num_threads 10 \
    --window_size 50 \
    --output_file data/ldhelmet/${outdir}/${base}.conf ${fname}

    sleep 1

    echo "table_gen for ${base}"

    time ./bin/ldhelmet table_gen \
    --num_threads 10 \
    --conf_file data/ldhelmet/${outdir}/${base}.conf \
    --theta 0.03 \
    --rhos 0.0 0.1 10.0 1.0 100.0 \
    --output_file data/ldhelmet/${outdir}/${base}.lk > table_gen 2> table_gen2
    
    sleep 1

    rm table_gen*

    sleep 1

    time ./bin/ldhelmet pade \
    --num_threads 10 \
    --conf_file data/ldhelmet/${outdir}/${base}.conf \
    --theta 0.03 \
    --output_file data/ldhelmet/${outdir}/${base}.pade

    sleep 1

    time ./bin/ldhelmet rjmcmc \
    --num_threads 30 \
    --window_size 50 \
    --seq_file ${fname} \
    --lk_file data/ldhelmet/${outdir}/${base}.lk \
    --pade_file data/ldhelmet/${outdir}/${base}.pade \
    --num_iter 1000000 \
    --burn_in 100000 \
    --block_penalty ${block} \
    --mut_mat_file data/ldhelmet/mut_mat \
    --output_file data/ldhelmet/${outdir}/${base}.post

    sleep 1

    time ./bin/ldhelmet post_to_text \
    --mean \
    --perc 0.025 \
    --perc 0.50 \
    --perc 0.975 \
    --output_file data/ldhelmet/${outdir}/${base}.txt \
    data/ldhelmet/${outdir}/${base}.post

    sleep 3

    echo "Removing temp files..."
    rm data/ldhelmet/${outdir}/${base}.conf
    rm data/ldhelmet/${outdir}/${base}.lk
    rm data/ldhelmet/${outdir}/${base}.pade
    rm data/ldhelmet/${outdir}/${base}.post

done 
```

running this over block = 5:

```bash
time bash run_ldhelmet.sh 5
```





























