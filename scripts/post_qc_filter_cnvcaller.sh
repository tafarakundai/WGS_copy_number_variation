#!/bin/bash

# CNV Filtering Script
# Filters CNVs based on silhouette score and length criteria

input_file="Genotype.tsv"
output_file="filtered_CNVs2.tsv"

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found!"
    exit 1
fi

echo "Filtering CNVs from $input_file..."
echo "Criteria:"
echo "- Silhouette score > 0.6"
echo "- Length = 50kb for deletions and both types"
echo "- Length < 500kb for duplications"
echo ""

# Create output file with header
head -n 1 "$input_file" > "$output_file"

# Filter the data with debugging
awk -F'\t' '
BEGIN {OFS="\t"}
NR == 1 {
    # Find the silhouette_score column index
    for (i = 1; i <= NF; i++) {
        if ($i == "silhouette_score") {
            sil_col = i
            break
        }
    }
    if (sil_col == 0) {
        print "Error: silhouette_score column not found!" > "/dev/stderr"
        exit 1
    }
    print "Silhouette score column found at position:", sil_col > "/dev/stderr"
    next
}
{
    # Calculate CNV length (end - start)
    cnv_length = $3 - $2
    
    # Get silhouette score using the found column index
    silhouette = $(sil_col)
    
    # Debug: Print first few rows to check values
    if (NR <= 5) {
        print "Row", NR, "Silhouette:", silhouette, "Length:", cnv_length > "/dev/stderr"
    }
    
    # Convert silhouette to number and check if it is valid
    silhouette_num = silhouette + 0
    
    # Filter based on silhouette score first
    if (silhouette_num > 0.55) {
        
        # Get genotype counts from the last columns
        # Find column indices for genotype columns
        dd_col = NF - 6   # dd
        ad_col = NF - 5   # Ad  
        aa_col = NF - 4   # AA
        ab_col = NF - 3   # AB
        bb_col = NF - 2   # BB
        bc_col = NF - 1   # BC
        m_col = NF        # M
        
        dd_count = $(dd_col) + 0
        ad_count = $(ad_col) + 0
        aa_count = $(aa_col) + 0
        ab_count = $(ab_col) + 0
        bb_count = $(bb_col) + 0
        bc_count = $(bc_col) + 0
        m_count = $(m_col) + 0
        
        # Determine CNV type based on genotypes
        deletion_count = dd_count + ad_count
        duplication_count = ab_count + bb_count + bc_count
        both_count = m_count
        
        # Apply length filters based on CNV type
        passed_filter = 0
        
        # Check for deletions (Ad or dd): length = 50kb
        if (deletion_count > 0 && cnv_length <= 50000) {
            passed_filter = 1
        }
        
        # Check for duplications (AB, BB, BC): length < 500kb
        if (duplication_count > 0 && cnv_length < 500000) {
            passed_filter = 1
        }
        
        # Check for both type (M): length = 50kb
        if (both_count > 0 && cnv_length <= 50000) {
            passed_filter = 1
        }
        
        # Print if passed any filter
        if (passed_filter == 1) {
            print $0
        }
    }
}' "$input_file" >> "$output_file"

# Count results
original_count=$(tail -n +2 "$input_file" | wc -l)
filtered_count=$(tail -n +2 "$output_file" | wc -l)

echo "Filtering completed!"
echo "Original CNVs: $original_count"
echo "Filtered CNVs: $filtered_count"
echo "Output saved to: $output_file"

# Display first few lines of output for verification
echo ""
echo "First 5 rows of filtered data:"
head -n 6 "$output_file"