#!/usr/bin/bash
#PBS -l select=1:ncpus=24:mpiprocs=24
#PBS -l walltime=96:00:00
#PBS -q smp
#PBS -P ######
#PBS -o /home/tafara/out_IBR.txt
#PBS -e /home/tafara/err_IBR.txt
#PBS -N Bargur
#PBS -M tafarakmavunga@gmail.com
#PBS -m abe

# Load necessary modules
module load chpc/BIOMODULES
module add samtools
module add bwa
module add bcftools
module add vcftools
module add fastp
module add htslib
module add picard/2.20.3


# Paths to data
DATA_DIR=/home/tafara/Bargur
REFERENCE=/home/tafara/Final_Data2/NIAB_ARS_BosIndicus_Tharparkar_1.0.fa.gz
REFERENCE_UNCOMPRESSED=${REFERENCE%.gz}  # Remove .gz extension

# Ensure reference exists
if [ ! -f "$REFERENCE" ]; then
    echo "Reference file not found: $REFERENCE"
    exit 1
fi

# Uncompress reference if required (only once)
if [ ! -f "$REFERENCE_UNCOMPRESSED" ]; then
    echo "Uncompressing reference genome..."
    gunzip -c "$REFERENCE" > "$REFERENCE_UNCOMPRESSED"
fi

# Index the reference genome with BWA if not already indexed
if [ ! -f "${REFERENCE}.bwt" ]; then
    echo "Indexing reference genome with BWA..."
    bwa index "$REFERENCE"
fi

# Create fai index if not already present
if [ ! -f "${REFERENCE_UNCOMPRESSED}.fai" ]; then
    echo "Creating FASTA index..."
    samtools faidx "$REFERENCE_UNCOMPRESSED"
fi

# Define tools
BWA="bwa"
SAMTOOLS="samtools"
VCFTOOLS="vcftools"
BCFTOOLS="bcftools"
FASTP="fastp"
PICARD="/apps/chpc/bio/picard/2.20.3/picard"

# Output directories
TRIMMED_DIR=$DATA_DIR/trimmed
ALIGN_DIR=$DATA_DIR/aligned
VARIANT_DIR=$DATA_DIR/variants
STATS_DIR=$DATA_DIR/statistics
mkdir -p $TRIMMED_DIR $ALIGN_DIR $VARIANT_DIR $STATS_DIR

# Array of sample prefixes
SAMPLES=("IBR1" "IBR2" "IBR3" "IBR4" "IBR5")

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    READ1=$DATA_DIR/${SAMPLE}_1.fastq.gz
    READ2=$DATA_DIR/${SAMPLE}_2.fastq.gz

    if [ ! -f "$READ1" ] || [ ! -f "$READ2" ]; then
        echo "Input files not found for sample $SAMPLE"
        continue
    fi

    echo "Processing $SAMPLE..."
    # Quality trimming
    $FASTP -i $READ1 -I $READ2 \
        -o $TRIMMED_DIR/${SAMPLE}_trimmed_1.fq.gz \
        -O $TRIMMED_DIR/${SAMPLE}_trimmed_2.fq.gz \
        --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 50 \
        --thread 12 --cut_right --overrepresentation_analysis \
        --html $STATS_DIR/${SAMPLE}_fastp.html --json $STATS_DIR/${SAMPLE}_fastp.json

    # Alignment with BWA-MEM
    $BWA mem -t 20 "$REFERENCE" $TRIMMED_DIR/${SAMPLE}_trimmed_1.fq.gz $TRIMMED_DIR/${SAMPLE}_trimmed_2.fq.gz | \
    $SAMTOOLS view -@ 20 -bS - | \
    $SAMTOOLS sort -@ 20 -o $ALIGN_DIR/${SAMPLE}_aligned_sorted.bam

    $SAMTOOLS index $ALIGN_DIR/${SAMPLE}_aligned_sorted.bam

    # Add read groups and mark duplicates
    $PICARD AddOrReplaceReadGroups I=$ALIGN_DIR/${SAMPLE}_aligned_sorted.bam O=$ALIGN_DIR/${SAMPLE}_rg.bam \
        RGID=${SAMPLE} RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=${SAMPLE} CREATE_INDEX=true

    $PICARD MarkDuplicates I=$ALIGN_DIR/${SAMPLE}_rg.bam O=$ALIGN_DIR/${SAMPLE}_dedup.bam \
        M=$STATS_DIR/${SAMPLE}_marked_dup_metrics.txt CREATE_INDEX=true
