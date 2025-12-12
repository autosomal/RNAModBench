#!/usr/bin/env python3
"""
Post-process NanoSPA results to generate standardized output format.
Converts NanoSPA predictions to BED-like format with filtering.
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path
import glob

def process_nanospa_results(input_dir, output_file, prob_threshold=0.5):
    """
    Process NanoSPA results and convert to standardized format.
    
    Parameters:
    -----------
    input_dir : str
        Path to NanoSPA results directory
    output_file : str
        Path to processed output file
    prob_threshold : float
        Probability threshold for filtering (default: 0.5)
    """
    
    # Find NanoSPA result files
    result_files = glob.glob(os.path.join(input_dir, "*.txt")) + \
                   glob.glob(os.path.join(input_dir, "*.tsv")) + \
                   glob.glob(os.path.join(input_dir, "*.csv"))
    
    # Look for m6A and psU results
    m6a_files = [f for f in result_files if 'm6a' in f.lower() or 'm6A' in f]
    psu_files = [f for f in result_files if 'psu' in f.lower() or 'psU' in f]
    
    all_results = []
    
    # Process m6A results
    for result_file in m6a_files:
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
            
            # Map columns
            if 'chr' in df.columns:
                standardized['Chr'] = df['chr']
            elif 'chromosome' in df.columns:
                standardized['Chr'] = df['chromosome']
            elif 'contig' in df.columns:
                standardized['Chr'] = df['contig']
            else:
                standardized['Chr'] = 'unknown'
            
            if 'pos' in df.columns:
                standardized['Start'] = df['pos']
            elif 'position' in df.columns:
                standardized['Start'] = df['position']
            else:
                standardized['Start'] = 0
            
            standardized['End'] = standardized['Start']
            standardized['Status'] = 'Mod'
            
            # Get probability/score
            prob_col = None
            for col in df.columns:
                if any(term in col.lower() for term in ['prob', 'score', 'pvalue', 'confidence']):
                    prob_col = col
                    break
            
            if prob_col:
                standardized['Prob'] = df[prob_col]
            else:
                standardized['Prob'] = 1.0
            
            standardized['Strand'] = df.get('strand', '*')
            standardized['mod_ratio'] = df.get('mod_ratio', standardized['Prob'])
            standardized['Modification_Type'] = 'm6A'
            
            # Apply filtering
            filtered = standardized[
                (standardized['Prob'] > prob_threshold) &
                (standardized['mod_ratio'] > 0.1)
            ].copy()
            
            all_results.append(filtered)
            
        except Exception as e:
            print(f"Error processing {result_file}: {e}", file=sys.stderr)
    
    # Process psU results
    for result_file in psu_files:
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
            
            # Standardize format (similar to m6A)
            standardized = pd.DataFrame()
            
            if 'chr' in df.columns:
                standardized['Chr'] = df['chr']
            elif 'chromosome' in df.columns:
                standardized['Chr'] = df['chromosome']
            else:
                standardized['Chr'] = 'unknown'
            
            if 'pos' in df.columns:
                standardized['Start'] = df['pos']
            elif 'position' in df.columns:
                standardized['Start'] = df['position']
            else:
                standardized['Start'] = 0
            
            standardized['End'] = standardized['Start']
            standardized['Status'] = 'Mod'
            
            # Get probability/score
            prob_col = None
            for col in df.columns:
                if any(term in col.lower() for term in ['prob', 'score', 'pvalue', 'confidence']):
                    prob_col = col
                    break
            
            if prob_col:
                standardized['Prob'] = df[prob_col]
            else:
                standardized['Prob'] = 1.0
            
            standardized['Strand'] = df.get('strand', '*')
            standardized['mod_ratio'] = df.get('mod_ratio', standardized['Prob'])
            standardized['Modification_Type'] = 'psU'
            
            # Apply filtering
            filtered = standardized[
                (standardized['Prob'] > prob_threshold) &
                (standardized['mod_ratio'] > 0.1)
            ].copy()
            
            all_results.append(filtered)
            
        except Exception as e:
            print(f"Error processing {result_file}: {e}", file=sys.stderr)
    
    if not all_results:
        print(f"No valid NanoSPA results found in {input_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Combine all results
    final_results = pd.concat(all_results, ignore_index=True)
    
    # Ensure integer positions
    final_results['Start'] = final_results['Start'].astype(int)
    final_results['End'] = final_results['End'].astype(int)
    
    # Reorder columns
    if 'Modification_Type' in final_results.columns:
        final_results = final_results[['Chr', 'Start', 'End', 'Status', 'Prob', 'Strand', 'mod_ratio', 'Modification_Type']]
    else:
        final_results = final_results[['Chr', 'Start', 'End', 'Status', 'Prob', 'Strand', 'mod_ratio']]
    
    # Save results
    final_results.to_csv(output_file, sep='\t', index=False, header=True)
    
    print(f"NanoSPA processing complete. {len(final_results)} modifications detected.")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python postprocess_nanospa.py <input_dir> <output_file>")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    prob_threshold = 0.5  # Default threshold
    
    process_nanospa_results(input_dir, output_file, prob_threshold)