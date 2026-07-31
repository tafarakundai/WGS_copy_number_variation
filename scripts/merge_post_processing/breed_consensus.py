import pandas as pd
from collections import defaultdict

# Define all 12 breeds
breeds = {
    "IBR": ["IBR1", "IBR2", "IBR3", "IBR4", "IBR5"],
    "IDN": ["IDN1", "IDN2", "IDN3", "IDN4", "IDN5"],
    "IHL": ["IHL1", "IHL2", "IHL3", "IHL4", "IHL5"],
    "IKA": ["IKA1", "IKA2", "IKA3", "IKA4", "IKA5"],
    "KDR": ["KDR1", "KDR2", "KDR3", "KDR4", "KDR5"],
    "LWC": ["LWC1", "LWC2", "LWC3", "LWC4", "LWC5"],
    "PDJ": ["PDJ3", "PDJ4", "PDJ8", "PDJ9", "PDJ10"],
    "PDN": ["PDN4", "PDN5", "PDN6", "PDN7", "PDN8"],
    "PRS": ["PRS1", "PRS2", "PRS3", "PRS4", "PRS5"],
    "PSH": ["PSH2", "PSH4", "PSH6", "PSH7", "PSH8"],
    "YPZ": ["YPZ1", "YPZ2", "YPZ3", "YPZ4", "YPZ5"],
    "YSN": ["YSN1", "YSN2", "YSN3", "YSN4", "YSN5"]
}

# Load consensus regions
consensus = pd.read_csv("strict_consensus_regions.bed", sep="\t", header=None, names=["chr", "start", "end"])

# Initialize results
results = []

# Process each filtered caller file
for caller in ["CNVnator", "CNVcaller", "Lumpy"]:
    filename = f"{caller}_consensus.tsv"
    df = pd.read_csv(filename, sep="\t")
    
    # For each consensus region
    for _, region in consensus.iterrows():
        chr_, start, end = region["chr"], region["start"], region["end"]
        
        # Find matching CNVs in this caller
        matches = df[(df["chr"] == chr_) & (df["start"] == start) & (df["end"] == end)]
        
        if not matches.empty:
            row = {
                "chr": chr_,
                "start": start,
                "end": end,
                "caller": caller,
            }
            
            # Count per breed
            for breed, samples in breeds.items():
                # Get samples present in this dataframe
                breed_samples = [s for s in samples if s in df.columns]
                
                if breed_samples:  # Only proceed if samples exist
                    # Count DEL (-1) and DUP (1)
                    del_count = (matches[breed_samples] == -1).sum().sum()
                    dup_count = (matches[breed_samples] == 1).sum().sum()
                    
                    row.update({
                        f"{breed}_DEL": del_count,
                        f"{breed}_DUP": dup_count,
                        f"{breed}_TOTAL": del_count + dup_count,
                    })
            
            results.append(row)

# Save detailed stats
detailed_df = pd.DataFrame(results)
detailed_df.to_csv("breed_stats_detailed.tsv", sep="\t", index=False)

# Create summary across all callers
summary = detailed_df.groupby(["chr", "start", "end"]).sum().reset_index()
summary.to_csv("breed_stats_summary.tsv", sep="\t", index=False)

# Create breed totals
breed_totals = []
for breed in breeds:
    breed_totals.append({
        "breed": breed,
        "DEL": summary[f"{breed}_DEL"].sum(),
        "DUP": summary[f"{breed}_DUP"].sum(),
        "TOTAL_CNVs": summary[f"{breed}_TOTAL"].sum(),
        "N_samples": len(breeds[breed])
    })

pd.DataFrame(breed_totals).to_csv("breed_totals.tsv", sep="\t", index=False)