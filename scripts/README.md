# CNV Analysis Scripts

This directory contains scripts used for whole genome sequencing (WGS) copy number variation (CNV) analysis.

## Pipeline Steps

### 1. Read Alignment

`alignment.sh`

Performs alignment of sequencing reads to the reference genome and prepares BAM files for downstream CNV analysis.

---

### 2. CNV Detection

`cnvnator.sh`

Identifies copy number variations using read-depth based CNV detection with CNVnator.

`cnvcaller.sh`

Identifies CNVs using CNVcaller.

`lumpy.sh`

Identifies CNVs using CNVcaller.

---

### 3. Structural Variant Detection

`lumpyexpress.config`

Configuration file for running LUMPY Express structural variant detection.

---

### 4. CNV Integration

`merge.R`

Combines CNV calls from different callers.

---

### 5. CNV Annotation

`annotation.R`

Annotates identified CNVs with genomic information.

---

### 6. Visualization
check postprocessing folder which has various scripts for bar plots etc

`figures.R`

Generates figures and summary plots.
