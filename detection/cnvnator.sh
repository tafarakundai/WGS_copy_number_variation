#!/bin/bash

# CNVnator CNV calling pipeline
# WGS Copy Number Variation Analysis

# Choose the appropriate bin size based on your genome coverage. For this we used 100bp bin sizes 

# Load modules
module add chpc/BIOMODULES
module add CNVnator

conda activate cnvnator

# Change to your working directory
cd /mnt/tafara/AsianZebu/Final_Data2/CNV

echo "Running CNVnator analysis"

echo "Step 1: Generate root file"
# Chromosome naming was not chr1 so had to use the specific chromosome names 

for i in *.bam; do
    cnvnator \
    -root "${i%.bam}_100.root" \
    -chrom \
    NC_091760.1 NC_091761.1 NC_091762.1 NC_091763.1 \
    NC_091764.1 NC_091765.1 NC_091766.1 NC_091767.1 \
    NC_091768.1 NC_091769.1 NC_091770.1 NC_091771.1 \
    NC_091772.1 NC_091773.1 NC_091774.1 NC_091775.1 \
    NC_091776.1 NC_091777.1 NC_091778.1 NC_091779.1 \
    NC_091780.1 NC_091781.1 NC_091782.1 NC_091783.1 \
    NC_091784.1 NC_091785.1 NC_091786.1 NC_091787.1 \
    NC_091788.1 NC_091789.1 \
    -tree "$i"
done

echo "Step 2: Calculate statistics"
for i in *_100.root; do cnvnator -root "$i" -his 100 -d ./ref_chromosomes/; done
for i in *_100.root; do cnvnator -root "$i" -stat 100; done
for i in *_100.root; do cnvnator -root "$i" -eval 100; done
for i in *_100.root; do cnvnator -root "$i" -partition 100; done


echo "Step 3: Calling CNVs"
for i in *_100.root; do cnvnator -root "$i" -call 100 > "${i%_100.root}_100_cnvs.csv"; done

echo "Step 3: Genotyping"
awk '{ print $2 } END { print "exit" }' IBR1_dedup_100_cnvs.csv | cnvnator -root IBR1_dedup_100.root -genotype 100

echo "CNVnator analysis completed"





