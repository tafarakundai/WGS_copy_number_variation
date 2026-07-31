#!/bin/bash

# First create a mapping file (chrom_map.txt)
echo -e "NC_091760.1\t1" > chrom_map.txt
echo -e "NC_091761.1\t2" >> chrom_map.txt
echo -e "NC_091762.1\t3" >> chrom_map.txt
echo -e "NC_091763.1\t4" >> chrom_map.txt
echo -e "NC_091764.1\t5" >> chrom_map.txt
echo -e "NC_091765.1\t6" >> chrom_map.txt
echo -e "NC_091766.1\t7" >> chrom_map.txt
echo -e "NC_091767.1\t8" >> chrom_map.txt
echo -e "NC_091768.1\t9" >> chrom_map.txt
echo -e "NC_091769.1\t10" >> chrom_map.txt
echo -e "NC_091770.1\t11" >> chrom_map.txt
echo -e "NC_091771.1\t12" >> chrom_map.txt
echo -e "NC_091772.1\t13" >> chrom_map.txt
echo -e "NC_091773.1\t14" >> chrom_map.txt
echo -e "NC_091774.1\t15" >> chrom_map.txt
echo -e "NC_091775.1\t16" >> chrom_map.txt
echo -e "NC_091776.1\t17" >> chrom_map.txt
echo -e "NC_091777.1\t18" >> chrom_map.txt
echo -e "NC_091778.1\t19" >> chrom_map.txt
echo -e "NC_091779.1\t20" >> chrom_map.txt
echo -e "NC_091780.1\t21" >> chrom_map.txt
echo -e "NC_091781.1\t22" >> chrom_map.txt
echo -e "NC_091782.1\t23" >> chrom_map.txt
echo -e "NC_091783.1\t24" >> chrom_map.txt
echo -e "NC_091784.1\t25" >> chrom_map.txt
echo -e "NC_091785.1\t26" >> chrom_map.txt
echo -e "NC_091786.1\t27" >> chrom_map.txt
echo -e "NC_091787.1\t28" >> chrom_map.txt
echo -e "NC_091788.1\t29" >> chrom_map.txt
echo -e "NC_091789.1\t30" >> chrom_map.txt

# Rename chromosomes using bcftools and output a compressed VCF
bcftools annotate --rename-chrs chrom_map.txt IBR_snps.vcf.gz -Oz -o IBR_snps_1.vcf.gz

# Index the renamed VCF (may be needed for some operations)
bcftools index IBR_snps_1.vcf.gz

