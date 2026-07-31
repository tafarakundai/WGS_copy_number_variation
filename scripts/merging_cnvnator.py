#!/usr/bin/env python3

import os
import sys
import subprocess
import pandas as pd
from collections import defaultdict
import re

def parse_coordinates(coord_str):
    """Parse coordinate string like 'NC_091760.1:296501-318900'"""
    chrom, positions = coord_str.split(':')
    start, end = map(int, positions.split('-'))
    return chrom, start, end

def calculate_overlap(start1, end1, start2, end2):
    """Calculate reciprocal overlap between two intervals"""
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    
    if overlap_start >= overlap_end:
        return 0.0
    
    overlap_length = overlap_end - overlap_start
    length1 = end1 - start1
    length2 = end2 - start2
    
    # Calculate reciprocal overlap (minimum of both directions)
    overlap1 = overlap_length / length1
    overlap2 = overlap_length / length2
    
    return min(overlap1, overlap2)

def read_cnv_files(file_pattern="*filtered.csv"):
    """Read all CNV files and return consolidated data"""
    cnv_data = []
    
    # Get all filtered.csv files
    import glob
    files = glob.glob(file_pattern)
    
    if not files:
        print(f"No files found matching pattern: {file_pattern}")
        return []
    
    print(f"Found {len(files)} CNV files")
    
    for file_path in files:
        # Extract sample name from filename
        sample_name = os.path.basename(file_path).replace('_filtered.csv', '').replace('.filtered.csv', '')
        
        try:
            # Read CSV file
            df = pd.read_csv(file_path, sep='\t')
            
            for _, row in df.iterrows():
                chrom, start, end = parse_coordinates(row['coordinates'])
                cnv_data.append({
                    'sample': sample_name,
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'cnv_type': row['CNV_type'],
                    'size': row['CNV_size']
                })
                
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            continue
    
    return cnv_data

def find_population_cnvs(cnv_data, min_individuals=3):
    """Find CNVs present in at least min_individuals, regardless of type"""
    # Group CNVs by chromosome only (not by type)
    chrom_cnvs = defaultdict(list)
    
    for cnv in cnv_data:
        chrom_cnvs[cnv['chrom']].append(cnv)
    
    population_cnvs = []
    
    for chrom, cnvs in chrom_cnvs.items():
        # Sort by start position
        cnvs.sort(key=lambda x: x['start'])
        
        processed = set()
        
        for i, cnv1 in enumerate(cnvs):
            if i in processed:
                continue
                
            overlapping_cnvs = [cnv1]
            processed.add(i)
            
            for j, cnv2 in enumerate(cnvs[i+1:], i+1):
                if j in processed:
                    continue
                    
                # Check for overlap (regardless of CNV type)
                overlap = calculate_overlap(cnv1['start'], cnv1['end'], 
                                         cnv2['start'], cnv2['end'])
                
                if overlap >= 0.5:  # 50% reciprocal overlap
                    overlapping_cnvs.append(cnv2)
                    processed.add(j)
            
            # Check if we have enough individuals (regardless of CNV type)
            unique_samples = set(cnv['sample'] for cnv in overlapping_cnvs)
            
            if len(unique_samples) >= min_individuals:
                # Calculate merged coordinates
                min_start = min(cnv['start'] for cnv in overlapping_cnvs)
                max_end = max(cnv['end'] for cnv in overlapping_cnvs)
                
                # Create sample -> cnv_type mapping
                sample_types = {}
                for cnv in overlapping_cnvs:
                    sample_types[cnv['sample']] = cnv['cnv_type']
                
                population_cnvs.append({
                    'chrom': chrom,
                    'start': min_start,
                    'end': max_end,
                    'sample_types': sample_types,  # Changed from 'samples' and 'cnv_type'
                    'count': len(unique_samples)
                })
    
    return population_cnvs

def merge_overlapping_cnvrs(population_cnvs):
    """Merge overlapping CNVRs while preserving individual sample CNV types"""
    # Group by chromosome
    chrom_cnvs = defaultdict(list)
    for cnv in population_cnvs:
        chrom_cnvs[cnv['chrom']].append(cnv)
    
    merged_cnvrs = []
    
    for chrom, cnvs in chrom_cnvs.items():
        # Sort by start position
        cnvs.sort(key=lambda x: x['start'])
        
        # Merge overlapping regions regardless of type
        i = 0
        while i < len(cnvs):
            current = cnvs[i]
            merged_start = current['start']
            merged_end = current['end']
            merged_sample_types = current['sample_types'].copy()
            
            # Check for overlaps with subsequent CNVs
            j = i + 1
            while j < len(cnvs):
                next_cnv = cnvs[j]
                
                # Check if there's any overlap
                if next_cnv['start'] <= merged_end:
                    # Calculate overlap
                    overlap = calculate_overlap(merged_start, merged_end,
                                             next_cnv['start'], next_cnv['end'])
                    
                    if overlap >= 0.5:  # 50% reciprocal overlap
                        # Merge coordinates
                        merged_end = max(merged_end, next_cnv['end'])
                        
                        # Merge sample types - each sample keeps its own CNV type
                        for sample, cnv_type in next_cnv['sample_types'].items():
                            if sample not in merged_sample_types:
                                merged_sample_types[sample] = cnv_type
                            else:
                                # Same sample with different CNV types in overlapping regions
                                # This could happen if a sample has both a deletion and duplication
                                # that overlap. Keep the first one encountered, but warn user.
                                if merged_sample_types[sample] != cnv_type:
                                    print(f"Warning: Sample {sample} has conflicting CNV types in overlapping region {chrom}:{merged_start}-{merged_end}")
                                    print(f"  Using: {merged_sample_types[sample]}, ignoring: {cnv_type}")
                        
                        # Remove merged CNV
                        cnvs.pop(j)
                    else:
                        j += 1
                else:
                    break
            
            # Add merged CNVR
            merged_cnvrs.append({
                'chrom': chrom,
                'start': merged_start,
                'end': merged_end,
                'sample_types': merged_sample_types,
                'cnvr_id': f"{chrom}:{merged_start}-{merged_end}"
            })
            
            i += 1
    
    return merged_cnvrs

def create_cnvr_matrix(merged_cnvrs, all_samples):
    """Create matrix with CNVRs and sample genotypes"""
    
    # Sort CNVRs by chromosome and position
    merged_cnvrs.sort(key=lambda x: (x['chrom'], x['start']))
    
    # Create matrix
    matrix_data = []
    
    for cnvr in merged_cnvrs:
        row = {
            'chr': cnvr['chrom'],
            'start': cnvr['start'],
            'end': cnvr['end'],
            'cnvr_id': cnvr['cnvr_id']
        }
        
        # Add sample columns - each sample gets its actual CNV type
        for sample in all_samples:
            if sample in cnvr['sample_types']:
                cnv_type = cnvr['sample_types'][sample]
                if cnv_type == 'duplication':
                    row[sample] = 1
                elif cnv_type == 'deletion':
                    row[sample] = -1
                else:
                    row[sample] = 0  # Unknown type
            else:
                row[sample] = 0  # Sample absent from this CNVR
        
        matrix_data.append(row)
    
    return pd.DataFrame(matrix_data)

def filter_by_population(cnv_data, populations):
    """Filter CNVs to include only samples from specified populations"""
    population_prefixes = populations
    filtered_data = []
    
    for cnv in cnv_data:
        sample_prefix = cnv['sample'][:3]  # Assuming first 3 characters indicate population
        if sample_prefix in population_prefixes:
            filtered_data.append(cnv)
    
    return filtered_data

def main():
    # Define populations
    populations = ['BGO', 'BZP', 'ZAG', 'MGZ']
    
    print("Step 1: Reading CNV files...")
    cnv_data = read_cnv_files("*filtered.csv")
    
    if not cnv_data:
        print("No CNV data found. Please check your file pattern.")
        return
    
    print(f"Loaded {len(cnv_data)} CNV records")
    
    # Get all samples
    all_samples = sorted(set(cnv['sample'] for cnv in cnv_data))
    print(f"Found {len(all_samples)} samples")
    
    # Filter by populations if needed
    if populations:
        print("Step 2: Filtering by populations...")
        cnv_data = filter_by_population(cnv_data, populations)
        all_samples = sorted(set(cnv['sample'] for cnv in cnv_data))
        print(f"After population filtering: {len(all_samples)} samples")
    
    print("Step 3: Finding population-based CNVs (=3 individuals)...")
    population_cnvs = find_population_cnvs(cnv_data, min_individuals=3)
    print(f"Found {len(population_cnvs)} population-based CNVs")
    
    print("Step 4: Merging overlapping CNVRs (50% reciprocal overlap)...")
    merged_cnvrs = merge_overlapping_cnvrs(population_cnvs)
    print(f"Created {len(merged_cnvrs)} CNVRs after merging")
    
    print("Step 5: Creating CNVR matrix...")
    cnvr_matrix = create_cnvr_matrix(merged_cnvrs, all_samples)
    
    # Save results
    print("Step 6: Saving results...")
    
    # Save CNVR matrix
    cnvr_matrix.to_csv('cnvr_matrix.csv', index=False)
    print("CNVR matrix saved to: cnvr_matrix.csv")
    
    # Save CNVR summary
    cnvr_summary = cnvr_matrix[['chr', 'start', 'end', 'cnvr_id']].copy()
    cnvr_summary['length'] = cnvr_summary['end'] - cnvr_summary['start']
    cnvr_summary['n_samples'] = cnvr_matrix[all_samples].abs().sum(axis=1)
    cnvr_summary.to_csv('cnvr_summary.csv', index=False)
    print("CNVR summary saved to: cnvr_summary.csv")
    
    # Print statistics
    print(f"\nFinal Statistics:")
    print(f"- Total CNVRs: {len(cnvr_matrix)}")
    print(f"- Total samples: {len(all_samples)}")
    print(f"- Average CNVR length: {cnvr_summary['length'].mean():.0f} bp")
    print(f"- Average samples per CNVR: {cnvr_summary['n_samples'].mean():.1f}")
    
    # Print population distribution
    pop_counts = defaultdict(int)
    for sample in all_samples:
        pop = sample[:3]  # First 3 characters
        pop_counts[pop] += 1
    
    print(f"\nPopulation distribution:")
    for pop, count in sorted(pop_counts.items()):
        print(f"- {pop}: {count} samples")

if __name__ == "__main__":
    main()