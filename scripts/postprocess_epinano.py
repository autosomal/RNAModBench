#!/usr/bin/env python3
"""
Post-process Epinano results to generate standardized output format.
Converts Epinano predictions to BED-like format with filtering.
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path

def process_epinano_results(input_file, output_file, delta_threshold=0.1):
    """
    Process Epinano results and convert to standardized format.
    
    Parameters:
    -----------
    input_file : str
        Path to Epinano raw output file (CSV format)
    output_file : str
        Path to processed output file
    delta_threshold : float
        Delta sum error threshold for filtering (default: 0.1)
    """
    
    # Read Epinano results
    try:
        data = pd.read_csv(input_file, header=0, sep=",")
    except Exception as e:
        print(f"Error reading input file {input_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Check required columns
    required_cols = ['chr_pos', 'delta_sum_err', 'z_score_prediction']
    if not all(col in data.columns for col in required_cols):
        print(f"Missing required columns in {input_file}. Expected: {required_cols}", file=sys.stderr)
        sys.exit(1)
    
    # Split chr_pos column
    split_cols = data['chr_pos'].str.split(expand=True)
    data['Chr'] = split_cols[0]
    data['Pos'] = pd.to_numeric(split_cols[1], errors='coerce')
    data['Base'] = split_cols[2]
    data['Strand'] = split_cols[3]
    
    # Convert delta_sum_err to numeric
    data['delta_sum_err'] = pd.to_numeric(data['delta_sum_err'], errors='coerce')
    
    # Filter valid data
    valid_data = data[
        data['Pos'].notna() & 
        (data['Pos'] > 0) & 
        data['Strand'].isin(['+', '-']) & 
        data['z_score_prediction'].isin(['mod', 'unm']) &
        data['delta_sum_err'].notna()
    ].copy()
    
    # Apply filtering criteria
    condition = (
        (valid_data['z_score_prediction'] == "mod") & 
        (valid_data['delta_sum_err'] > delta_threshold) & 
        (
            ((valid_data['Strand'] == "+") & (valid_data['Base'] == "A")) | 
            ((valid_data['Strand'] == "-") & (valid_data['Base'] == "T"))
        )
    )
    
    # Create final result
    result = pd.DataFrame({
        "Chr": valid_data['Chr'],
        "Start": valid_data['Pos'].astype(int),
        "End": valid_data['Pos'].astype(int) + 1,  # BED format
        "Status": valid_data['z_score_prediction'],
        "delta_sum_err": valid_data['delta_sum_err'],
        "Strand": valid_data['Strand']
    }).loc[condition]
    
    # Reorder columns
    result = result[["Chr", "Start", "End", "Status", "delta_sum_err", "Strand"]]
    
    # Save results
    result.to_csv(output_file, sep='\t', index=False, header=True)
    
    print(f"Epinano processing complete. {len(result)} modifications detected.")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python postprocess_epinano.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    delta_threshold = 0.1  # Default threshold
    
    process_epinano_results(input_file, output_file, delta_threshold)