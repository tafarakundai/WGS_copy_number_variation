#!/bin/bash

# Step 1: Save the header from one of the files
head -n 1 CNVnator.tsv > consensus_merged.tsv

# Step 2: Process each tool's calls
bedtools intersect -a strict_consensus_regions.bed -b CNVnator.tsv -wa -wb | cut -f1-3,7- > cnvnator_hits.tsv
bedtools intersect -a strict_consensus_regions.bed -b CNVcaller.tsv -wa -wb | cut -f1-3,7- > cnvcaller_hits.tsv
bedtools intersect -a strict_consensus_regions.bed -b Lumpy.txt -wa -wb | cut -f1-3,7- > lumpy_hits.tsv

# Step 3: Merge calls from all three tools
awk '
BEGIN { OFS="\t" }
# Read CNVnator calls
FILENAME == "cnvnator_hits.tsv" {
    key = $1 "\t" $2 "\t" $3
    for (i=4; i<=NF; i++) {
        cnvnator[key][i-3] = $i
    }
    next
}
# Read CNVcaller calls
FILENAME == "cnvcaller_hits.tsv" {
    key = $1 "\t" $2 "\t" $3
    for (i=4; i<=NF; i++) {
        cnvcaller[key][i-3] = $i
    }
    next
}
# Read Lumpy calls
FILENAME == "lumpy_hits.tsv" {
    key = $1 "\t" $2 "\t" $3
    for (i=4; i<=NF; i++) {
        lumpy[key][i-3] = $i
    }
    next
}
# Process consensus regions
FILENAME == "strict_consensus_regions.bed" {
    key = $1 "\t" $2 "\t" $3
    if (key in cnvnator && key in cnvcaller && key in lumpy) {
        # Print coordinates
        printf "%s\t%s\t%s", $1, $2, $3
        
        # Calculate length (if needed)
        len = $3 - $2
        printf "\t%s", len
        
        # Process each sample (assuming 60 samples)
        for (i=1; i<=60; i++) {
            sum = cnvnator[key][i] + cnvcaller[key][i] + lumpy[key][i]
            # Majority vote: -1 if =2 dels, 1 if =2 dups, else 0
            if (sum <= -2) call = -1
            else if (sum >= 2) call = 1
            else call = 0
            printf "\t%s", call
        }
        
        # Calculate Total column
        total = 0
        for (i=1; i<=60; i++) {
            sum = cnvnator[key][i] + cnvcaller[key][i] + lumpy[key][i]
            if (sum <= -2 || sum >= 2) total++
        }
        printf "\t%s\n", total
    }
}
' cnvnator_hits.tsv cnvcaller_hits.tsv lumpy_hits.tsv strict_consensus_regions.bed >> consensus_merged.tsv

# Cleanup temporary files
rm -f cnvnator_hits.tsv cnvcaller_hits.tsv lumpy_hits.tsv