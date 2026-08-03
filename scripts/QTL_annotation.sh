#!/bin/bash

#Annotate candidate CNV regions to QTL 

#Load modules if working on cluster
module add chpc/BIOMODULES
module add bedtools/2.27.1

# Convert CNV file to bed 
tail -n +2 regions_renamed.txt | awk -F'\t' 'BEGIN{OFS="\t"} {print "chr"$1, $2, $3, $4}' > CNVRs_regions.bed

# Download QTL database from AnimalQTLdb
# Clean QTL file 
awk -F'\t' 'BEGIN{OFS="\t"}
    $1 ~ /^Chr\./ && $2 < $3 {
        gsub(/^Chr\./, "chr", $1);
        print $1, $2, $3, $4
    }' QTLdb_cattleARS_UCD1.bed > QTLdb_cleaned.bed


#Run the intersection (can copy the code and make a small bash file to run)
#!/bin/bash

# Define input/output files
CNVR_FILE="CNVRs_regions.bed"
QTL_FILE="QTLdb_cleaned.bed"
OUTPUT_FILE="CNVRs_with_QTLinfo.bed"

# Clean and sort both files
sed 's/[[:space:]]*$//' "$CNVR_FILE" | sort -k1,1 -k2,2n > CNVRs.sorted.bed
sed 's/[[:space:]]*$//' "$QTL_FILE" | sort -k1,1 -k2,2n > QTLdb.sorted.bed

# Intersect to find CNVRs overlapping with QTLs
bedtools intersect -a CNVRs.sorted.bed -b QTLdb.sorted.bed -wa -wb > "$OUTPUT_FILE"

# Preview result
echo "Output written to $OUTPUT_FILE"
echo "Example lines:"
head "$OUTPUT_FILE"


# FILTER UNIQUE QTLS (copy the code into a python file)

import pandas as pd
import re

# Load file
df = pd.read_csv("QTL_all_regions.txt", sep="\t")

# Extract base QTL name (without the ID in parentheses)
df['QTL_base_name'] = df['QTL_name'].apply(lambda x: re.sub(r'\s+\(\d+\)', '', x))

# Group by base name and aggregate
aggregated_df = df.groupby('QTL_base_name').agg({
    'Chromosome': 'first',
    'start': 'first',
    'end': 'first',
    'Type': 'first',
    'QTL_chr': 'first',
    'start.1': 'min',
    'end.1': 'max',
    'QTL_name': 'first'  # keep one full name with ID
}).reset_index()

# Save to file
aggregated_df.to_csv("aggregated_qtls_by_trait2.txt", sep="\t", index=False)

# Optional: print result
print(aggregated_df)


#LOGIC FOR FILTERING 
Same trait in the same CNVR region- merged into one row (removes redundancy)
Same trait in different CNVR regions - kept as separate rows
Same trait on different chromosomes - kept as separate rows
