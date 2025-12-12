#!/usr/bin/env python3
"""
Post-process DENA results to generate standardized output format.
Converts DENA predictions to BED-like format with filtering.
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path

def process_dena_results(input_file, output_file, ratio_threshold=0.1, coverage_threshold=20):
    """
    Process DENA results and convert to standardized format.
    
    Parameters:
    -----------
    input_file : str
        Path to DENA raw output file (TSV format)
    output_file : str
        Path to processed output file
    ratio_threshold : float
        Modification ratio threshold for filtering (default: 0.1)
    coverage_threshold : int
        Minimum coverage threshold (default: 20)
    """
    
    # Read DENA results
    try:
        data_dena = pd.read_csv(input_file, header=None, sep='\t')
    except Exception as e:
        print(f"Error reading input file {input_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Check minimum required columns
    if data_dena.shape[1] < 6:
        print(f"Insufficient columns in {input_file}. Expected at least 6 columns.", file=sys.stderr)
        sys.exit(1)
    
    # Create standardized output format
    dena = pd.DataFrame({
        "Chr": data_dena[0],
        "Start": data_dena[1],
        "End": data_dena[1],
        "Strand": "*",
        "Status": data_dena[5].apply(lambda x: "mod" if x > ratio_threshold else "unm"),
        "mod_ratio": data_dena[5]
    })
    
    # Calculate coverage
    data_dena['coverage'] = data_dena[3] + data_dena[4]
    
    # Filter for coverage and modifications
    dena = dena[data_dena['coverage'] > coverage_threshold]
    dena = dena[dena['Status'] == "mod"]
    
    # Select and reorder columns
    dena = dena[['Chr', 'Start', 'End', 'Status', 'mod_ratio', 'Strand']]
    
    # Adjust start position for BED format
    dena['Start'] = dena['Start'] - 1
    
    # Save results
    dena.to_csv(output_file, sep='\t', index=False, header=True)
    
    print(f"DENA processing complete. {len(dena)} modifications detected.")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python postprocess_dena.py <input_file> <output_file> <ratio_threshold> <coverage_threshold>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    ratio_threshold = float(sys.argv[3])
    coverage_threshold = int(sys.argv[4])
    
    process_dena_results(input_file, output_file, ratio_threshold, coverage_threshold)