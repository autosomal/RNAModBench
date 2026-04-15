#!/usr/bin/env python3
"""
Post-process NanoNm results to generate standardized output format.
Converts NanoNm predictions to BED-like format with filtering.
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path
import glob

def process_nanonm_results(input_dir, output_file, ratio_threshold=0.1, coverage_threshold=20):
    """
    Process NanoNm results and convert to standardized format.
    
    Parameters:
    -----------
    input_dir : str
        Path to NanoNm results directory
    output_file : str
        Path to processed output file
    ratio_threshold : float
        Modification ratio threshold for filtering (default: 0.1)
    coverage_threshold : int
        Minimum coverage threshold for filtering (default: 20)
    """
    
    # Find NanoNm result files
    result_files = glob.glob(os.path.join(input_dir, "*.txt")) + \
                   glob.glob(os.path.join(input_dir, "*.tsv")) + \
                   glob.glob(os.path.join(input_dir, "*.csv"))
    
    if not result_files:
        print(f"No result files found in {input_dir}", file=sys.stderr)
        sys.exit(1)
    
    all_results = []
    
    for result_file in result_files:
        try:
            # Try different separators
            for sep in ['\t', ',', ' ']:
                try:
                    df = pd.read_csv(result_file, sep=sep, header=0)
                    break
                except:
                    continue
            else:
                continue
            
            # Standardize format
            standardized = pd.DataFrame()
            
            # Map chromosome/contig column
            if 'chr' in df.columns:
                standardized['Chr'] = df['chr']
            elif 'chromosome' in df.columns:
                standardized['Chr'] = df['chromosome']
            elif 'contig' in df.columns:
                standardized['Chr'] = df['contig']
            elif 'transcript_id' in df.columns:
                standardized['Chr'] = df['transcript_id']
            else:
                standardized['Chr'] = 'unknown'
            
            # Map position column
            if 'pos' in df.columns:
                standardized['Start'] = df['pos']
            elif 'position' in df.columns:
                standardized['Start'] = df['position']
            elif 'start' in df.columns:
                standardized['Start'] = df['start']
            else:
                standardized['Start'] = 0
            
            standardized['End'] = standardized['Start']
            standardized['Status'] = 'Mod'
            
            # Get modification ratio/probability
            ratio_col = None
            prob_col = None
            coverage_col = None
            
            for col in df.columns:
                col_lower = col.lower()
                if any(term in col_lower for term in ['ratio', 'mod_ratio', 'modification_ratio']):
                    ratio_col = col
                elif any(term in col_lower for term in ['prob', 'score', 'probability', 'confidence']):
                    prob_col = col
                elif any(term in col_lower for term in ['coverage', 'depth', 'support', 'count']):
                    coverage_col = col
            
            # Set probability
            if prob_col:
                standardized['Prob'] = df[prob_col]
            elif ratio_col:
                standardized['Prob'] = df[ratio_col]
            else:
                standardized['Prob'] = 1.0
            
            # Set strand
            if 'strand' in df.columns:
                standardized['Strand'] = df['strand']
            else:
                standardized['Strand'] = '*'
            
            # Set modification ratio
            if ratio_col:
                standardized['mod_ratio'] = df[ratio_col]
            else:
                standardized['mod_ratio'] = standardized['Prob']
            
            # Set coverage
            if coverage_col:
                standardized['Coverage'] = df[coverage_col]
            else:
                standardized['Coverage'] = 0
            
            standardized['Modification_Type'] = 'Nm'
            
            # Apply filtering
            filtered = standardized[
                (standardized['mod_ratio'] > ratio_threshold) &
                (standardized['Coverage'] >= coverage_threshold)
            ].copy()
            
            all_results.append(filtered)
            
        except Exception as e:
            print(f"Error processing {result_file}: {e}", file=sys.stderr)
    
    if not all_results:
        print(f"No valid NanoNm results found in {input_dir} after filtering", file=sys.stderr)
        # Create empty output file with headers
        empty_df = pd.DataFrame(columns=['Chr', 'Start', 'End', 'Status', 'Prob', 'Strand', 'mod_ratio', 'Coverage', 'Modification_Type'])
        empty_df.to_csv(output_file, sep='\t', index=False, header=True)
        print(f"Empty results saved to {output_file}")
        return
    
    # Combine all results
    final_results = pd.concat(all_results, ignore_index=True)
    
    # Ensure integer positions
    final_results['Start'] = final_results['Start'].astype(int)
    final_results['End'] = final_results['End'].astype(int)
    
    # Reorder columns
    final_results = final_results[['Chr', 'Start', 'End', 'Status', 'Prob', 'Strand', 'mod_ratio', 'Coverage', 'Modification_Type']]
    
    # Save results
    final_results.to_csv(output_file, sep='\t', index=False, header=True)
    
    print(f"NanoNm processing complete. {len(final_results)} modifications detected.")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python postprocess_nanonm.py <input_dir> <output_file>")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    ratio_threshold = 0.1
    coverage_threshold = 20
    
    process_nanonm_results(input_dir, output_file, ratio_threshold, coverage_threshold)
