#!/bin/bash

# Navigate to the new directory containing CNV files

# Create log file with timestamp
log_file="cnv_filtering_log_$(date +%Y%m%d_%H%M%S).txt"

# Filter each CSV file based on the criteria:
# q0 < 0.5 (column 9)
# CNV_size > 1000 (column 3, assuming size is in bp)
# e-val1 < 0.01 (column 5)

for file in *.csv; do
    if [[ -f "$file" ]]; then
        # Extract base sample name (remove path and any suffixes)
        sample_name=$(basename "$file" | cut -d'_' -f1)
        
        # Output to both terminal and log file
        echo "Processing $file (sample: $sample_name)..." | tee -a "$log_file"
        
        # Create filtered filename with just sample name
        filtered_file="${sample_name}_filtered.csv"
        
        # Extract header and save to filtered file
        head -n 1 "$file" > "$filtered_file"
        
        # Filter data based on criteria and append to filtered file
        # Skip header (NR>1), then apply filters:
        # $3 > 1000 (CNV_size > 1kb)
        # $5 < 0.01 (e-val1 < 0.01)
        # $9 < 0.5 (q0 < 0.5)
        awk -F'\t' 'NR>1 && $3 > 1000 && $5 < 0.01 && $9 < 0.5' "$file" >> "$filtered_file"
        
        # Count original and filtered lines
        orig_lines=$(($(wc -l < "$file") - 1))  # Subtract header
        filt_lines=$(($(wc -l < "$filtered_file") - 1))  # Subtract header
        
        echo "  Original CNVs: $orig_lines" | tee -a "$log_file"
        echo "  Filtered CNVs: $filt_lines" | tee -a "$log_file"
        echo "  Filtered file: $filtered_file" | tee -a "$log_file"
        echo "" | tee -a "$log_file"
    fi
done

echo "Filtering complete!" | tee -a "$log_file"
echo "Filtered files saved with '_filtered.csv' suffix" | tee -a "$log_file"
echo "Log file saved as: $log_file" | tee -a "$log_file"