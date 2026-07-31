# 1. Load modules
module add chpc/BIOMODULES
module add samtools
module add perl
module add cnvcaller

# 2. Activate conda environment
eval "$(conda shell.bash hook)"
conda activate cnvcaller

# 3. Create directories for results
mkdir -p results

# 4. Generate reference database with 800bp windows
perl CNVReferenceDB.pl -f /mnt/lustre/users/mmalima/All_projects/AsianZebu/Final_Data2/CNV/CNVcaller/NIAB_ARS_BosIndicus_Tharparkar_1.0.fa -w 800 -o results/Tharparkar_refDB

#CREATE THE LINK FILE ON YOUR OWN 

#Create the reference.fa.sa file:
sawriter NIAB_ARS_BosIndicus_Tharparkar_1.0.fa.sa NIAB_ARS_BosIndicus_Tharparkar_1.0.fa

#ALIGN THE KMER TO THE FASTA REFERENCE GENOME
blasr kmer.fa NIAB_ARS_BosIndicus_Tharparkar_1.0.fa \
  --sa NIAB_ARS_BosIndicus_Tharparkar_1.0.fa.sa \
  --out kmer.aln -m 5 \
  --noSplitSubreads --minMatch 15 --maxMatch 20 \
  --advanceHalf --advanceExactMatches 10 \
  --fastMaxInterval --fastSDP --aggressiveIntervalCut --bestn 10

#Generate duplicated window record file
module add python/3.9.6

python 0.2.Kmer_Link.py kmer.aln 800 Tharparkar.link

# Run Individual.Process.sh on example samples

  bash /home/apps/chpc/bio/cnvcaller/Individual.Process.sh \
    -b /mnt/tafara/CNVcaller/IBR1_dedup.bam \
    -h IBR1 \
    -d Tharparkar.link \
    -s NC_091789.1
    
  bash /home/apps/chpc/bio/cnvcaller/Individual.Process.sh \
    -b /mnt/tafara/CNVcaller/IBR2_dedup.bam \
    -h IBR2 \
    -d Tharparkar.link \
    -s NC_091789.1
 
  bash /home/apps/chpc/bio/cnvcaller/Individual.Process.sh \
     -b /mnt/tafara/CNVcaller/IBR3_dedup.bam \
     -h IBR3 \
     -d Tharparkar.link \
    -s NC_091789.1

  bash /home/apps/chpc/bio/cnvcaller/Individual.Process.sh \
     -b /mnt/tafara/CNVcaller/IBR4_dedup.bam \
     -h IBR4 \
     -d Tharparkar.link \
    -s NC_091789.1
    
   bash /home/apps/chpc/bio/cnvcaller/Individual.Process.sh \
      -b /mnt/tafara/CNVcaller/IBR5_dedup.bam \
      -h IBR5 \
      -d Tharparkar.link \
    -s NC_091789.1

# CNV DETECTION
# create a list of the samples in each breed, this is because CNVcaller first run individual analysis then merge later

ls RD_normalized/IBR* > single_list.txt

#create an exclude list as well
touch exclude_list.txt

# 9. Run CNV discovery
bash CNV.Discovery.sh \
  -l single_list.txt \
  -e exclude_list.txt \
  -f 0.1 \
  -h 3 \
  -r 0.5 \
  -p primaryCNVR_IBR.txt \
  -m mergedCNVR_IBR.txt

# If you have more breeds then you can also run code below to run for all samples at once 

ls RD_normalized/* > main_list.txt

touch exclude_list.txt


bash CNV.Discovery.sh \
  -l main_list.txt \
  -e exclude_list.txt \
  -f 0.1 \
  -h 3 \
  -r 0.5 \
  -p primaryCNVR_AA.txt \
  -m mergedCNVR_AA.txt
  
  bash CNV.Discovery.sh \
    -l main_listB.txt \
    -e exclude_list.txt \
    -f 0.1 \
    -h 3 \
    -r 0.5 \
    -p primaryCNVR_Africa.txt \
    -m mergedCNVR_Africa.txt
  
# GENOTYPING

python Genotype.py --cnvfile mergedCNVR_AA.txt --outprefix Genotype

