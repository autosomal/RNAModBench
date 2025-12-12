#!/usr/bin/env python3
"""
Post-process xPore results to generate standardized output format.
Converts xPore predictions to BED-like format with filtering.
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path
import glob

def process_xpore_results(input_dir, output_file, prob_threshold=0.5):
    """
    Process xPore results and convert to standardized format.
    
    Parameters:
    -----------
    input_dir : str
        Path to xPore diffmod output directory
    output_file : str
        Path to processed output file
    prob_threshold : float
        Probability threshold for filtering (default: 0.5)
    """
    
    # Find the results file (xPore uses different naming conventions)
    result_files = glob.glob(os.path.join(input_dir, "*.csv")) + \
                   glob.glob(os.path.join(input_dir, "*.tsv")) + \
                   glob.glob(os.path.join(input_dir, "results", "*.csv"))
    
    if not result_files:
        print(f"Could not find xPore results in {input_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Use the most likely result file
    result_file = result_files[0]
    for f in result_files:
        if 'diffmod' in f.lower() or 'results' in f.lower():
            result_file = f
            break
    
    # Read xPore results
    try:
        df = pd.read_csv(result_file, sep=',', header=0)
    except Exception as e:
        print(f"Error reading input file {result_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Look for probability column (xPore uses different naming)
    prob_col = None
    for col in df.columns:
        if any(term in col.lower() for term in ['prob', 'pvalue', 'p_value', 'significance']):
            prob_col = col
            break
    
    if prob_col is None:
        print(f"Could not find probability column in {result_file}", file=sys.stderr)
        sys.exit(1)
    
    # Convert to standardized format
    # xPore typically has: transcript_id, position, reference_kmer, probability
    standardized = pd.DataFrame()
    
    # Map common column patterns
    if 'transcript_id' in df.columns:
        standardized['Chr'] = df['transcript_id']
    elif 'contig' in df.columns:
        standardized['Chr'] = df['contig']
    elif 'chr' in df.columns:
        standardized['Chr'] = df['chr']
    else:
        standardized['Chr'] = 'unknown'
    
    if 'position' in df.columns:
        standardized['Start'] = df['position']
    elif 'pos' in df.columns:
        standardized['Start'] = df['pos']
    else:
        standardized['Start'] = 0
    
    standardized['End'] = standardized['Start']
    standardized['Status'] = 'Mod'
    standardized['Prob'] = df[prob_col]
    standardized['Strand'] = df.get('strand', '*')
    standardized['mod_ratio'] = df.get('mod_ratio', df[prob_col])
    
    # Apply filtering
    filtered = standardized[
        (standardized['Prob'] > prob_threshold) &
        (standardized['mod_ratio'] > 0.1)
    ].copy()
    
    # Ensure integer positions
    filtered['Start'] = filtered['Start'].astype(int)
    filtered['End'] = filtered['End'].astype(int)
    
    # Reorder columns
    filtered = filtered[['Chr', 'Start', 'End', 'Status', 'Prob', 'Strand', 'mod_ratio']]
    
    # Save results
    filtered.to_csv(output_file, sep='\t', index=False, header=True)
    
    print(f"xPore processing complete. {len(filtered)} modifications detected.")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python postprocess_xpore.py <input_dir> <output_file>")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    prob_threshold = 0.5  # Default threshold
    
    process_xpore_results(input_dir, output_file, prob_threshold)