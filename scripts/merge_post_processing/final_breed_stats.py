import pandas as pd

# Load data
summary = pd.read_csv("breed_stats_summary.tsv", sep="\t")
consensus = pd.read_csv("strict_consensus_regions.bed", sep="\t", header=None, names=["chr", "start", "end"])

# Add length to each region
consensus["length"] = consensus["end"] - consensus["start"] + 1

# Merge with summary stats
merged = pd.merge(summary, consensus, on=["chr", "start", "end"])

# Define breeds (from your list)
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

# Calculate statistics
results = []
for breed in breeds:
    breed_del = merged[merged[f"{breed}_DEL"] > 0]
    breed_dup = merged[merged[f"{breed}_DUP"] > 0]
    breed_both = merged[(merged[f"{breed}_DEL"] > 0) & (merged[f"{breed}_DUP"] > 0)]
    breed_any = merged[(merged[f"{breed}_DEL"] > 0) | (merged[f"{breed}_DUP"] > 0)]
    
    results.append({
        "breed": breed,
        # Region counts
        "DEL_only": breed_del.shape[0] - breed_both.shape[0],
        "DUP_only": breed_dup.shape[0] - breed_both.shape[0],
        "BOTH": breed_both.shape[0],
        "Total_regions": breed_any.shape[0],
        # Length sums (bp)
        "DEL_length": breed_del["length"].sum(),
        "DUP_length": breed_dup["length"].sum(),
        "BOTH_length": breed_both["length"].sum(),
        "Total_length": breed_any["length"].sum(),
        "N_samples": len(breeds[breed])
    })

# Create final dataframe
stats_df = pd.DataFrame(results)

# Reorder columns
cols = [
    "breed", "N_samples",
    "Total_regions", "Total_length",
    "DEL_only", "DEL_length",
    "DUP_only", "DUP_length", 
    "BOTH", "BOTH_length"
]
stats_df = stats_df[cols]

# Save results
stats_df.to_csv("enhanced_breed_stats.tsv", sep="\t", index=False)