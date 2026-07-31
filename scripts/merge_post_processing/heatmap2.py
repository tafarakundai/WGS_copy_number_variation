import pandas as pd
import numpy as np
from scipy.stats import kruskal
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys

# ---- USER CONFIGURATION ----
INPUT_FILE = "consensus_cnvs.tsv"
OUTPUT_DIR = "cnv_population_analysis"

# Robust directory creation
def create_output_directory(dir_name):
    """Create output directory with robust error handling"""
    try:
        os.makedirs(dir_name, exist_ok=True)
        test_file = os.path.join(dir_name, "test_write.tmp")
        with open(test_file, 'w') as f:
            f.write("test")
        os.remove(test_file)
        print(f"Output directory created: {os.path.abspath(dir_name)}")
        return dir_name
    except Exception as e:
        print(f"Error with relative path: {e}")

        abs_dir = os.path.abspath(dir_name)
        try:
            os.makedirs(abs_dir, exist_ok=True)
            test_file = os.path.join(abs_dir, "test_write.tmp")
            with open(test_file, 'w') as f:
                f.write("test")
            os.remove(test_file)
            print(f"Output directory created: {abs_dir}")
            return abs_dir
        except Exception as e:
            print(f"Error with absolute path: {e}")

            fallback_dir = os.path.join(os.getcwd(), dir_name)
            try:
                os.makedirs(fallback_dir, exist_ok=True)
                print(f"Fallback directory created: {fallback_dir}")
                return fallback_dir
            except Exception as e:
                print(f"All directory creation methods failed: {e}")
                print("Using current directory as output location")
                return "."

OUTPUT_DIR = create_output_directory(OUTPUT_DIR)

# Mapping: breed code/prefix -> full breed name
BREED_FULLNAME = {
    'BreedA': 'Breed Alpha',
    'BreedB': 'Breed Beta',
    'BreedC': 'Breed Gamma',
    # Add more breed mappings as needed
}

# ---- LOAD DATA ----
print("Loading data...")
try:
    df = pd.read_csv(INPUT_FILE, sep='\t')
    print(f"Data loaded: {df.shape[0]} rows, {df.shape[1]} columns")
except Exception as e:
    print(f"Error loading data: {e}")
    sys.exit(1)

# Identify sample columns
genomic_cols = ['chr', 'start', 'end', 'cnv_id', 'length', 'type']
sample_cols = [c for c in df.columns if c not in genomic_cols]
print(f"Found {len(sample_cols)} sample columns")

# Extract breed code from sample names
def get_breed_code(sample):
    import re
    m = re.match(r'^([A-Za-z]+)', sample)
    if m:
        return m.group(1)
    return None

breed_codes = [get_breed_code(s) for s in sample_cols]
breed_fullnames = [BREED_FULLNAME.get(code, code) for code in breed_codes]

# Build sample info DataFrame
sample_info = pd.DataFrame({
    'Sample': sample_cols,
    'BreedCode': breed_codes,
    'Breed': breed_fullnames  
})

# Filter valid samples
valid_samples = sample_info[~sample_info['Breed'].isnull()]
df_samples = df[valid_samples['Sample']]

# Create CNVR IDs
if 'cnv_id' in df.columns:
    cnvr_ids = df['cnv_id']
else:
    cnvr_ids = df['chr'].astype(str) + '_' + df['start'].astype(str) + '_' + df['end'].astype(str)

# Transpose data
cnv_data = df_samples.T
cnv_data.columns = cnvr_ids
cnv_data['Breed'] = valid_samples['Breed'].values
cnv_data['Sample'] = valid_samples['Sample'].values
cnv_data = cnv_data.set_index('Sample')

# Filter breeds with >1 sample
breed_counts = cnv_data['Breed'].value_counts()
valid_breeds = breed_counts[breed_counts > 1].index.tolist()
cnv_data = cnv_data[cnv_data['Breed'].isin(valid_breeds)]
unique_breeds = sorted(valid_breeds)

print(f"Analyzing {len(unique_breeds)} breeds: {unique_breeds}")

# ---- VST CALCULATION ----
print("Calculating VST values...")

def calc_vst_pair(cnv_df, breed_series, breed1, breed2):
    samples1 = breed_series[breed_series == breed1].index
    samples2 = breed_series[breed_series == breed2].index
    data1 = cnv_df.loc[samples1].drop(columns=['Breed'])
    data2 = cnv_df.loc[samples2].drop(columns=['Breed'])
    vst_values = []
    for cnvr in data1.columns:
        vals1 = data1[cnvr].astype(float).values
        vals2 = data2[cnvr].astype(float).values
        if len(vals1) > 1 and len(vals2) > 1:
            vt = np.var(np.concatenate([vals1, vals2]), ddof=1)
            vs = ((len(vals1)-1)*np.var(vals1, ddof=1) + (len(vals2)-1)*np.var(vals2, ddof=1)) / (len(vals1)+len(vals2)-2)
            vst = (vt - vs) / vt if vt > 0 else 0.0
        else:
            vst = np.nan
        vst_values.append(vst)
    return np.array(vst_values)

# Calculate pairwise VST
mean_vst_matrix = pd.DataFrame(np.nan, index=unique_breeds, columns=unique_breeds)
vst_pair_dict = {}

for i, b1 in enumerate(unique_breeds):
    for j, b2 in enumerate(unique_breeds):
        if i < j:
            vst_vals = calc_vst_pair(cnv_data, cnv_data['Breed'], b1, b2)
            mean_vst = np.nanmean(vst_vals)
            mean_vst_matrix.loc[b1, b2] = mean_vst
            mean_vst_matrix.loc[b2, b1] = mean_vst
            vst_pair_dict[f'{b1}_vs_{b2}'] = vst_vals
            print(f"  {b1} vs {b2}: mean VST = {mean_vst:.6f}")

# Safe file saving function
def safe_save(df, filename, **kwargs):
    filepath = os.path.join(OUTPUT_DIR, filename)
    try:
        df.to_csv(filepath, **kwargs)
        print(f"Saved: {filename}")
        return True
    except Exception as e:
        print(f"Error saving {filename}: {e}")
        return False

# Save VST matrix
safe_save(mean_vst_matrix, "pairwise_VST_matrix.tsv", sep='\t', float_format='%.6f')

# ---- HEATMAP (Lower triangle only) ----
print("Creating heatmap...")
try:
    mask = np.triu(np.ones_like(mean_vst_matrix, dtype=bool))  # mask upper triangle
    plt.figure(figsize=(10, 8))
    sns.heatmap(mean_vst_matrix, annot=True, cmap='coolwarm', square=True,
                xticklabels=mean_vst_matrix.columns,
                yticklabels=mean_vst_matrix.index,
                fmt='.3f', linewidths=0.5, mask=mask, cbar_kws={"shrink": 0.8})
    plt.title('Mean Pairwise VST Between Breeds (Lower Triangle Only)', fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=10)
    plt.yticks(rotation=0, fontsize=10)
    plt.tight_layout()

    heatmap_path = os.path.join(OUTPUT_DIR, "mean_vst_heatmap_lower_triangle.png")
    plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    print("Heatmap saved")
except Exception as e:
    print(f"Error creating heatmap: {e}")

# ---- REST OF YOUR CODE (Kruskal-Wallis, thresholds, results saving) ----
# Unchanged from your script
