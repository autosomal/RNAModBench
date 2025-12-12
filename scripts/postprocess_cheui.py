#!/usr/bin/env python3
"""
Post-process CHEUI results to generate standardized output format.
Converts CHEUI predictions to BED-like format with filtering.
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path

def process_cheui_results(input_file, output_file, prob_threshold=0.999, ratio_threshold=0.1):
    """
    Process CHEUI results and convert to standardized format.
    
    Parameters:
    -----------
    input_file : str
        Path to CHEUI raw output file
    output_file : str
        Path to processed output file
    prob_threshold : float
        Probability threshold for filtering (default: 0.999)
    ratio_threshold : float
        Modification ratio threshold for filtering (default: 0.1)
    """
    
    # Read CHEUI results
    try:
        df = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        print(f"Error reading input file {input_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Check required columns
    required_cols = ['contig', 'position', 'probability', 'stoichiometry']
    if not all(col in df.columns for col in required_cols):
        print(f"Missing required columns in {input_file}. Expected: {required_cols}", file=sys.stderr)
        sys.exit(1)
    
    # Create standardized output format
    result_df = pd.DataFrame({
        'Chr': df['contig'],
        'Start': df['position'] + 4,  # Adjust position
        'End': df['position'] + 5,    # BED format
        'Status': 'Mod',
        'Prob': df['probability'],
        'Strand': '*',
        'mod_ratio': df['stoichiometry']
    })
    
    # Apply filtering criteria
    filtered_df = result_df[
        (result_df['Prob'] > prob_threshold) & 
        (result_df['mod_ratio'] > ratio_threshold)
    ].copy()
    
    # Ensure integer positions
    filtered_df['Start'] = filtered_df['Start'].astype(int)
    filtered_df['End'] = filtered_df['End'].astype(int)
    
    # Save results
    filtered_df.to_csv(
        output_file,
        sep='\t',
        index=False,
        columns=['Chr', 'Start', 'End', 'Status', 'Prob', 'Strand', 'mod_ratio']
    )
    
    print(f"CHEUI processing complete. {len(filtered_df)} modifications detected.")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python postprocess_cheui.py <input_file> <output_file> <prob_threshold> <ratio_threshold>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    prob_threshold = float(sys.argv[3])
    ratio_threshold = float(sys.argv[4])
    
    process_cheui_results(input_file, output_file, prob_threshold, ratio_threshold)