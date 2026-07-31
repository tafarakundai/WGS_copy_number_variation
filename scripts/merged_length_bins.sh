#!/bin/bash
# Script to analyze CNVRs by size category and type (Del/Dup/Both)

# Output files
SIZE_DATA="cnvr_size_data.tsv"
TEMP_BED="temp_cnvr.bed"

# Define size ranges in kb
SIZE_RANGES=("1_2" "2_5" "5_10" "10_20" "20_50" "50_100" "100+")
SIZE_MIN=(1000 2000 5000 10000 20000 50000 100000)
SIZE_MAX=(2000 5000 10000 20000 50000 100000 999999999)

# Initialize output file with header
echo -e "Size\tDel\tDup\tBoth" > "$SIZE_DATA"

# Initialize counters for each size range
for i in "${!SIZE_RANGES[@]}"; do
    del_counts[$i]=0
    dup_counts[$i]=0
    both_counts[$i]=0
done

# Step 1: Extract header to identify columns
header=$(head -n1 filtered_CNVs.tsv)
header_arr=($header)

# Step 2: Process the genotype file
awk 'BEGIN {FS=OFS="\t"}
NR==1 {print; next}  # Skip header
{
    # Calculate length
    len = $3 - $2
    
    # Count types of calls
    del_count = 0
    dup_count = 0
    
    # Start from column 9 (first sample)
    for(i=9; i<=NF; i++) {
        if($i == "dd" || $i == "Ad") del_count++
        else if($i == "AB" || $i == "BB" || $i == "BC") dup_count++
        else if($i == "M") { del_count++; dup_count++ }
    }
    
    # Only keep CNVR with at least 2 samples
    total_count = 0
    for(i=9; i<=NF; i++) {
        if($i != "AA" && $i != "NA") total_count++
    }
    
    if(total_count >= 2) {
        # Determine type (Del/Dup/Both)
        type = "NA"
        if(del_count >= 2 && dup_count >= 2) type = "Both"
        else if(del_count >= 2) type = "Del"
        else if(dup_count >= 2) type = "Dup"
        
        if(type != "NA") {
            print $1, $2, $3, len, type
        }
    }
}' filtered_CNVs2.tsv | cut -f1-5 > "$TEMP_BED"

# Step 3: Count CNVRs by size and type
while IFS=$'\t' read -r chr start end length type; do
    # Skip if not numeric (header)
    if ! [[ "$length" =~ ^[0-9]+$ ]]; then
        continue
    fi
    
    # Find appropriate size bin
    size_bin=-1
    for i in "${!SIZE_MIN[@]}"; do
        if [ "$length" -ge "${SIZE_MIN[$i]}" ] && [ "$length" -lt "${SIZE_MAX[$i]}" ]; then
            size_bin=$i
            break
        fi
    done
    
    # If valid size bin found
    if [ $size_bin -ne -1 ]; then
        if [ "$type" == "Del" ]; then
            del_counts[$size_bin]=$((del_counts[$size_bin] + 1))
        elif [ "$type" == "Dup" ]; then
            dup_counts[$size_bin]=$((dup_counts[$size_bin] + 1))
        elif [ "$type" == "Both" ]; then
            both_counts[$size_bin]=$((both_counts[$size_bin] + 1))
        fi
    fi
done < "$TEMP_BED"

# Step 4: Write size counts to file
for i in "${!SIZE_RANGES[@]}"; do
    echo -e "${SIZE_RANGES[$i]}\t${del_counts[$i]}\t${dup_counts[$i]}\t${both_counts[$i]}" >> "$SIZE_DATA"
done

echo "Analysis complete. Results saved to $SIZE_DATA"

# Clean up
rm -f "$TEMP_BED"