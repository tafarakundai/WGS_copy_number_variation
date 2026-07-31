import pandas as pd

# Define all breeds (update with your actual breed-sample mappings)
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
consensus = pd.read_csv("strict_consensus_regions.bed", sep="\t", header=None, 
                       names=["chr", "start", "end"])

# Calculate length and classify size
consensus["length"] = consensus["end"] - consensus["start"] + 1

def classify_size(length):
    if length < 1000:
        return "<1k"
    elif length < 2000:
        return "1_2k"
    elif length < 5000:
        return "2_5k"
    elif length < 10000:
        return "5_10k"
    elif length < 20000:
        return "10_20k"
    elif length < 50000:
        return "20_50k"
    elif length < 100000:
        return "50_100k"
    else:
        return ">100k"

consensus["size_class"] = consensus["length"].apply(classify_size)

# Load breed stats to determine DEL/DUP/BOTH
summary = pd.read_csv("breed_stats_summary.tsv", sep="\t")
merged = pd.merge(consensus, summary, on=["chr", "start", "end"])

# Classify CNV type per region (DEL/DUP/BOTH)
def classify_type(row):
    has_del = any(row[f"{breed}_DEL"] > 0 for breed in breeds)
    has_dup = any(row[f"{breed}_DUP"] > 0 for breed in breeds)
    if has_del and has_dup:
        return "BOTH"
    elif has_del:
        return "DEL"
    elif has_dup:
        return "DUP"
    else:
        return "UNK"

merged["cnv_type"] = merged.apply(classify_type, axis=1)

# Count size classes by type
size_counts = pd.crosstab(
    index=merged["size_class"],
    columns=merged["cnv_type"],
    margins=True,
    margins_name="Total"
)

# Calculate length sums by type and size class
length_sums = merged.groupby(["size_class", "cnv_type"])["length"].sum().unstack()
length_sums["Total"] = length_sums.sum(axis=1)

# Combine counts and lengths
result = pd.concat([
    size_counts.rename(columns=lambda x: f"{x}_count"),
    length_sums.rename(columns=lambda x: f"{x}_bp")
], axis=1)

# Reorder columns
cols = []
for t in ["DEL", "DUP", "BOTH", "Total"]:
    cols.extend([f"{t}_count", f"{t}_bp"])
result = result[cols]

# Save results
result.to_csv("cnv_size_distribution.tsv", sep="\t")
consensus.to_csv("consensus_with_lengths.tsv", sep="\t", index=False)

print("Results saved to:")
print("- cnv_size_distribution.tsv (size class counts)")
print("- consensus_with_lengths.tsv (all CNVs with lengths)")