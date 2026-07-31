#!/usr/bin/env python3
"""
VCF Structural Variant Processor
Processes LUMPY/SVtyper VCF files to extract structural variant information
"""

import re
import os
import argparse
from typing import List, Dict, Tuple, Optional

def parse_info_field(info_string: str) -> Dict[str, str]:
    """Parse the INFO field from VCF and return a dictionary of key-value pairs"""
    info_dict = {}
    
    # Split by semicolon and process each key-value pair
    for item in info_string.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            # Flag fields (no value)
            info_dict[item] = 'True'
    
    return info_dict

def extract_sv_info(vcf_line: str) -> Optional[Dict[str, str]]:
    """Extract structural variant information from a VCF line"""
    fields = vcf_line.strip().split('\t')
    
    if len(fields) < 8:
        return None
    
    chrom = fields[0]
    pos = fields[1]
    ref = fields[3]
    alt = fields[4]
    info = fields[7]
    
    # Parse INFO field
    info_dict = parse_info_field(info)
    
    # Extract SVTYPE
    svtype = info_dict.get('SVTYPE', 'UNKNOWN')
    
    # Extract size information
    svlen = info_dict.get('SVLEN', '0')
    
    # Handle multiple SVLEN values (take the first one)
    if ',' in svlen:
        svlen = svlen.split(',')[0]
    
    # Convert SVLEN to absolute value for size
    try:
        size = abs(int(svlen))
    except ValueError:
        size = 0
    
    # Extract END position
    end_pos = info_dict.get('END', pos)
    
    # For BND (breakend) variants, calculate approximate size differently
    if svtype == 'BND':
        # Try to extract partner position from ALT field
        bnd_match = re.search(r'NC_091760\.1:(\d+)', alt)
        if bnd_match:
            partner_pos = int(bnd_match.group(1))
            current_pos = int(pos)
            size = abs(partner_pos - current_pos)
            end_pos = str(max(current_pos, partner_pos))
        else:
            size = 0
            end_pos = pos
    
    return {
        'SVTYPE': svtype,
        'CHROMOSOME': chrom,
        'START': pos,
        'END': end_pos,
        'SIZE': str(size)
    }

def process_vcf_file(input_file: str, output_file: str, sv_types: List[str] = None) -> None:
    """
    Process VCF file and extract structural variant information
    
    Args:
        input_file: Path to input VCF file
        output_file: Path to output text file
        sv_types: List of SV types to filter (None = all types)
    """
    
    variants = []
    
    try:
        with open(input_file, 'r') as f:
            for line in f:
                # Skip header lines
                if line.startswith('#'):
                    continue
                
                # Process variant line
                sv_info = extract_sv_info(line)
                if sv_info:
                    # Filter by SV type if specified
                    if sv_types is None or sv_info['SVTYPE'] in sv_types:
                        variants.append(sv_info)
    
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        return
    except Exception as e:
        print(f"Error reading VCF file: {e}")
        return
    
    # Write results to output file
    try:
        with open(output_file, 'w') as f:
            # Write header
            f.write("SVTYPE\tCHROMOSOME\tSTART\tEND\tSIZE\n")
            
            # Write variant data
            for variant in variants:
                f.write(f"{variant['SVTYPE']}\t{variant['CHROMOSOME']}\t"
                       f"{variant['START']}\t{variant['END']}\t{variant['SIZE']}\n")
        
        print(f"Successfully processed {len(variants)} variants.")
        print(f"Output written to: {output_file}")
        
        # Print summary statistics
        sv_counts = {}
        for variant in variants:
            svtype = variant['SVTYPE']
            sv_counts[svtype] = sv_counts.get(svtype, 0) + 1
        
        print("\nSummary by SV type:")
        for svtype, count in sorted(sv_counts.items()):
            print(f"  {svtype}: {count}")
            
    except Exception as e:
        print(f"Error writing output file: {e}")

def main():
    """Main function with command line interface"""
    parser = argparse.ArgumentParser(
        description="Process VCF files from LUMPY/SVtyper to extract structural variant information",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all SV types
  python vcf_processor.py input.vcf output.txt
  
  # Process only deletions and duplications
  python vcf_processor.py input.vcf output.txt --sv-types DEL DUP
  
  # Process only large variants (>1kb) - you can modify the script for size filtering
  python vcf_processor.py input.vcf output.txt --sv-types DEL DUP INV
        """
    )
    
    parser.add_argument('input_vcf', help='Input VCF file')
    parser.add_argument('output_file', help='Output text file')
    parser.add_argument('--sv-types', nargs='*', 
                       help='SV types to include (DEL, DUP, INV, BND, INS, CNV). Default: all types')
    
    args = parser.parse_args()
    
    # Validate SV types if provided
    valid_sv_types = {'DEL', 'DUP', 'INV', 'BND', 'INS', 'CNV'}
    if args.sv_types:
        invalid_types = set(args.sv_types) - valid_sv_types
        if invalid_types:
            print(f"Warning: Invalid SV types ignored: {invalid_types}")
            args.sv_types = [t for t in args.sv_types if t in valid_sv_types]
    
    # Process the VCF file
    process_vcf_file(args.input_vcf, args.output_file, args.sv_types)

# Function for direct use without command line
def process_vcf_simple(input_file: str, output_file: str = None, sv_types: List[str] = None):
    """
    Simplified function for direct use in scripts
    
    Args:
        input_file: Path to input VCF file
        output_file: Path to output file (default: input_file + '_processed.txt')
        sv_types: List of SV types to filter (default: all)
    """
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_processed.txt"
    
    process_vcf_file(input_file, output_file, sv_types)
    return output_file

if __name__ == "__main__":
    main()