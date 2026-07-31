import pandas as pd

# Load your detailed CNV file
df = pd.read_csv('lumpy_region_analysis_detailed.csv')

# Map your region_type values to simple CNV types
def classify_region_type(rt):
    rt = rt.lower()
    if 'deletion' in rt or 'del' in rt:
        return 'DEL'
    elif 'duplication' in rt or 'dup' in rt:
        return 'DUP'
    elif 'mixed' in rt or 'both' in rt:
        return 'BOTH'
    else:
        return 'UNK'

df['CNVtype'] = df['region_type'].apply(classify_region_type)

# Save output with CNVtype column
output = df[['chr', 'start', 'end', 'CNVtype']]
output.to_csv('combined_consensus2.bed', sep='\t', index=False, header=False)

# Print summary
print(f"Total regions processed: {len(output)}")
print("CNV type distribution:")
print(output['CNVtype'].value_counts())
print("\nFirst few lines:")
print(output.head())
