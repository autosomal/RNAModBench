#!/usr/bin/env python3
"""
Post-process yanocomp results to generate standardized output format.
Converts yanocomp predictions to BED-like format with filtering.
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path

def process_yanocomp_results(input_file, output_file, pvalue_threshold=0.05):
    """
    Process yanocomp results and convert to standardized format.
    
    Parameters:
    -----------
    input_file : str
        Path to yanocomp raw output file (BED format)
    output_file : str
        Path to processed output file
    pvalue_threshold : float
        P-value threshold for filtering (default: 0.05)
    """
    
    # Read yanocomp results (BED format)
    try:
        df = pd.read_csv(input_file, sep='\t', header=None)
    except Exception as e:
        print(f"Error reading input file {input_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Check minimum required columns (BED format has at least 3)
    if df.shape[1] < 3:
        print(f"Insufficient columns in {input_file}. Expected at least 3 columns for BED format.", file=sys.stderr)
        sys.exit(1)
    
    # Set column names based on BED format
    bed_cols = ['Chr', 'Start', 'End']
    for i, col in enumerate(bed_cols):
        if i < df.shape[1]:
            df.columns.values[i] = col
    
    # Add additional columns if present
    if df.shape[1] > 3:
        df.columns.values[3] = 'Name' if df.shape[1] > 3 else 'Score'
    if df.shape[1] > 4:
        df.columns.values[4] = 'Score'
    if df.shape[1] > 5:
        df.columns.values[5] = 'Strand'
    
    # Extract p-value from the file if available (yanocomp may include it in additional columns)
    # If not available, we'll use a default filtering based on the input being already filtered
    pvalue_col = None
    for col in df.columns:
        if any(term in str(col).lower() for term in ['pvalue', 'p_value', 'pval', 'significance']):
            pvalue_col = col
            break
    
    # Create standardized output format
    standardized = pd.DataFrame()
    standardized['Chr'] = df['Chr']
    standardized['Start'] = df['Start']
    standardized['End'] = df['End']
    standardized['Status'] = 'Mod'  # yanocomp output is already filtered
    
    # Try to extract probability/score
    if 'Score' in df.columns:
        standardized['Prob'] = df['Score']
    elif pvalue_col is not None:
        standardized['Prob'] = 1 - df[pvalue_col]  # Convert p-value to probability
    else:
        standardized['Prob'] = 1.0  # Default if no score available
    
    standardized['Strand'] = df.get('Strand', '*')
    standardized['mod_ratio'] = df.get('Score', standardized['Prob'])  # Use score as proxy
    
    # Apply additional filtering if p-value is available
    if pvalue_col is not None:
        filtered = standardized[df[pvalue_col] < pvalue_threshold].copy()
    else:
        filtered = standardized.copy()
    
    # Ensure integer positions
    filtered['Start'] = filtered['Start'].astype(int)
    filtered['End'] = filtered['End'].astype(int)
    
    # Reorder columns
    filtered = filtered[['Chr', 'Start', 'End', 'Status', 'Prob', 'Strand', 'mod_ratio']]
    
    # Save results
    filtered.to_csv(output_file, sep='\t', index=False, header=True)
    
    print(f"yanocomp processing complete. {len(filtered)} modifications detected.")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python postprocess_yanocomp.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    pvalue_threshold = 0.05  # Default threshold
    
    process_yanocomp_results(input_file, output_file, pvalue_threshold)