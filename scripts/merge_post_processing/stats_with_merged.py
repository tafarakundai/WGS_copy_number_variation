import pandas as pd
from collections import defaultdict
import os

# Configuration - Update these paths if needed
INPUT_FILES = {
    'consensus': "strict_consensus_regions.bed",
    'CNVnator': "CNVnator.tsv",
    'CNVcaller': "CNVcaller.tsv", 
    'Lumpy': "Lumpy.txt"
}

OUTPUT_FILES = {
    'typed_consensus': "consensus_regions_with_types.bed",
    'stats': "consensus_type_statistics.txt"
}

def load_caller_file(filename):
    """Flexible loader for CNV caller files"""
    try:
        df = pd.read_csv(filename, sep='\t')
        print(f"  Loaded {filename} ({len(df)} rows)")
        
        # Standardize column names
        col_map = {}
        for col in df.columns:
            col_lower = col.lower()
            if 'type' in col_lower or 'svtype' in col_lower:
                col_map[col] = 'CNVtype'
            elif 'chrom' in col_lower or 'chr' in col_lower:
                col_map[col] = 'chr'
            elif 'start' in col_lower or 'begin' in col_lower:
                col_map[col] = 'start'
            elif 'end' in col_lower or 'stop' in col_lower:
                col_map[col] = 'end'
        
        if col_map:
            df = df.rename(columns=col_map)
        
        return df
    except Exception as e:
        print(f"  Error loading {filename}: {str(e)}")
        return None

def main():
    print("Starting CNV type classification...")
    
    # 1. Load consensus regions
    print(f"\nLoading consensus regions from {INPUT_FILES['consensus']}")
    consensus = pd.read_csv(INPUT_FILES['consensus'], sep='\t', header=None,
                          names=["chr", "start", "end"])
    print(f"  Found {len(consensus)} consensus regions")

    # 2. Load all caller files
    print("\nLoading caller files:")
    callers = {}
    for name, path in INPUT_FILES.items():
        if name != 'consensus':
            callers[name] = load_caller_file(path)
    
    # 3. Classify each consensus region
    print("\nClassifying consensus regions...")
    type_counts = defaultdict(int)
    
    def classify_region(row):
        chr_, start, end = row['chr'], row['start'], row['end']
        types = set()
        
        for caller_name, df in callers.items():
            if df is None:
                continue
                
            # Find overlapping calls
            matches = df[(df['chr'] == chr_) & 
                        (df['start'] <= end) & 
                        (df['end'] >= start)]
            
            # Record types
            for _, match in matches.iterrows():
                cnv_type = str(match.get('CNVtype', '')).upper()
                if 'DEL' in cnv_type or '0' in cnv_type or '-1' in cnv_type:
                    types.add('DEL')
                if 'DUP' in cnv_type or '1' in cnv_type:
                    types.add('DUP')
        
        # Determine classification
        if not types:
            classification = 'UNK'
        elif 'DEL' in types and 'DUP' in types:
            classification = 'BOTH'
        elif 'DEL' in types:
            classification = 'DEL'
        elif 'DUP' in types:
            classification = 'DUP'
        else:
            classification = 'UNK'
        
        type_counts[classification] += 1
        return classification

    consensus['CNVtype'] = consensus.apply(classify_region, axis=1)
    
    # 4. Save outputs
    print("\nSaving results:")
    print(f"  Typed consensus to {OUTPUT_FILES['typed_consensus']}")
    consensus.to_csv(OUTPUT_FILES['typed_consensus'], sep='\t', header=False, index=False)
    
    print(f"  Statistics to {OUTPUT_FILES['stats']}")
    with open(OUTPUT_FILES['stats'], 'w') as f:
        f.write("Consensus CNV Type Classification Report\n")
        f.write("="*50 + "\n")
        f.write(f"Total consensus regions: {len(consensus)}\n")
        f.write(f"Deletions (DEL): {type_counts['DEL']}\n")
        f.write(f"Duplications (DUP): {type_counts['DUP']}\n")
        f.write(f"Mixed (BOTH): {type_counts['BOTH']}\n")
        f.write(f"Unknown: {type_counts.get('UNK', 0)}\n")
    
    # 5. Print summary
    print("\nClassification Complete!")
    print(f"• Deletions (DEL): {type_counts['DEL']}")
    print(f"• Duplications (DUP): {type_counts['DUP']}")
    print(f"• Mixed (BOTH): {type_counts['BOTH']}")
    if type_counts.get('UNK', 0) > 0:
        print(f"• Unknown: {type_counts['UNK']} (check input files)")

if __name__ == "__main__":
    main()