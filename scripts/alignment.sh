#!/usr/bin/bash
# ============================================================
# SIMPLIFIED BAM GENERATION PIPELINE (for teaching purposes)
# ============================================================

# 1. Define your files 
REFERENCE=reference.fa          # reference genome (already uncompressed & indexed)
READ1=sample_1.fastq.gz         # forward reads
READ2=sample_2.fastq.gz         # reverse reads
SAMPLE=sample1                  # sample name, used for read groups & output naming

OUTDIR=results
mkdir -p $OUTDIR

# 2. Trim adapters & low-quality bases with fastp
fastp \
    -i $READ1 -I $READ2 \
    -o $OUTDIR/${SAMPLE}_trimmed_1.fq.gz \
    -O $OUTDIR/${SAMPLE}_trimmed_2.fq.gz \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --html $OUTDIR/${SAMPLE}_fastp.html \
    --json $OUTDIR/${SAMPLE}_fastp.json

# 3. Align reads to the reference with BWA-MEM 
bwa mem $REFERENCE \
    $OUTDIR/${SAMPLE}_trimmed_1.fq.gz \
    $OUTDIR/${SAMPLE}_trimmed_2.fq.gz \
    > $OUTDIR/${SAMPLE}.sam

#  4. Convert SAM -> sorted BAM 
samtools sort -o $OUTDIR/${SAMPLE}_sorted.bam $OUTDIR/${SAMPLE}.sam
samtools index $OUTDIR/${SAMPLE}_sorted.bam

# 5. Add read groups (required by many downstream tools, e.g. GATK) 
picard AddOrReplaceReadGroups \
    I=$OUTDIR/${SAMPLE}_sorted.bam \
    O=$OUTDIR/${SAMPLE}_rg.bam \
    RGID=${SAMPLE} RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=${SAMPLE} \
    CREATE_INDEX=true

# 6. Mark PCR/optical duplicates 
picard MarkDuplicates \
    I=$OUTDIR/${SAMPLE}_rg.bam \
    O=$OUTDIR/${SAMPLE}_dedup.bam \
    M=$OUTDIR/${SAMPLE}_dup_metrics.txt \
    CREATE_INDEX=true

echo "Done. Final BAM: $OUTDIR/${SAMPLE}_dedup.bam"
