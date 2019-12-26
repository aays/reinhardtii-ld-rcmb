#!/bin/bash

# usage:
# bash variant_calling.sh [dirname] [reference fasta]
#
# e.g. bash variant_calling data/fastq/ data/references/chlamy.5.3.fasta
# this script only contains non-Flowers fastq data - see Flowers 2015 for download/assembly of those samples

dirname=$1
ref_fasta=$2

mkdir -p data/sam
mkdir -p data/bam
mkdir -p data/vcf

# download paired-end fastqs from ENA before running this script
# https://www.ebi.ac.uk/ena/data/view/PRJEB27323

for strain in {82..96}; do
    fname=ERR26411${strain}

    echo "bwa aln ${fname}"
    time bwa mem ${ref_fasta} ${dirname}/${fname}_1.fastq.gz ${dirname}/${fname}_2.fastq.gz > data/sam/${fname}.sam;

    echo "sam to bam ${fname}"
    time samtools view -S -b data/sam/${fname}.sam > data/bam/${fname}.bam;

    echo "sort bam ${fname}"
    samtools sort -T ${ref_fasta} -o data/bam/${fname}.sorted.bam data/bam/${fname}.bam;

    echo "index bam ${fname}"
    samtools index data/bam/${fname}.sorted.bam;

    echo "create gvcf ${fname}"
    java -Xmx16g -jar /opt/gatk/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R chlamy.5.3.w_organelles_mtMinus.fasta \
        --ploidy=1 \
        --includeNonVariantSites=true \
        --heterozygosity=0.02,
        --indel_heterozygosity=0.002 \
        -I data/bam/${fname}.sorted.bam \
        -o data/vcf/${fname}.g.vcf;
done

# run genotype gvcfs
for i in {1..17}; do
    java -jar /opt/gatk/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    -R ../data/chlamy.5.3.w_organelles_mtMinus.fasta \
    --variant ERR2641182.g.vcf \
    --variant ERR2641183.g.vcf \
    --variant ERR2641184.g.vcf \
    --variant ERR2641185.g.vcf \
    --variant ERR2641186.g.vcf \
    --variant ERR2641187.g.vcf \
    --variant ERR2641188.g.vcf \
    --variant ERR2641189.g.vcf \
    --variant ERR2641190.g.vcf \
    --variant ERR2641191.g.vcf \
    --variant ERR2641192.g.vcf \
    --variant ERR2641194.g.vcf \
    --variant ERR2641195.g.vcf \
    --variant ERR2641196.g.vcf \
    -L chromosome_${i} \
    -o data/vcf/chromosome_${i}.vcf;
done

# bgzip + tabix
for i in {1..17}; do
    bgzip data/vcf/chromosome_${i}.vcf
    sleep 1
    tabix -p vcf data/vcf/chromosome_${i}.vcf.gz
done 

# filtering

## quality filtering + MAF filtering

mkdir -p data/vcf/filtered
mkdir -p data/vcf/filtered/temp

for i in {1..17}; do
    echo "chromosome_${i}";
    ./bin/bcftools filter -i 'AVG(GQ)>30 && TYPE="snp" && MAF>0.1' \
    data/vcf/chromosome_${i}.vcf.gz \
    data/vcf/filtered/temp/chromosome_${i}.vcf;
done

## filtering for diallelic sites
for i in {1..17}; do
    time python3.5 analysis/sampling/filter_diallelic.py \
    data/vcf/filtered/temp/chromosome_${i}.vcf \
    data/vcf/filtered/chromosome_${i}.vcf

## bgzip and tabix filtered vcfs
for i in {1..17}; do
    echo "chromosome_${i}";
    bgzip data/vcf/filtered/chromosome_${i}.vcf ;
    tabix -p vcf data/vcf/filtered/chromosome_${i}.vcf.gz;
done
