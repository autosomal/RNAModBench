#!/usr/bin/env python3
"""
Post-process ELIGOS2 results to generate standardized output format.
Converts ELIGOS2 predictions to BED-like format with filtering.
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path

def process_eligos2_results(input_file, output_file, padj_threshold=0.0001, oddr_threshold=1.2):
    """
    Process ELIGOS2 results and convert to standardized format.
    
    Parameters:
    -----------
    input_file : str
        Path to ELIGOS2 raw output file
    output_file : str
        Path to processed output file
    padj_threshold : float
        Adjusted p-value threshold for filtering (default: 0.0001)
    oddr_threshold : float
        Odds ratio threshold for filtering (default: 1.2)
    """
    
    # Read ELIGOS2 results
    try:
        data_eligos = pd.read_csv(input_file, sep='\t', header=0)
    except Exception as e:
        print(f"Error reading input file {input_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Check required columns
    required_cols = ['ref', 'total_reads', 'chrom', 'start_loc', 'end_loc', 'strand', 'adjPval', 'oddR', 'ESB_test']
    if not all(col in data_eligos.columns for col in required_cols):
        print(f"Missing required columns in {input_file}. Expected: {required_cols}", file=sys.stderr)
        sys.exit(1)
    
    # Filter for adenosine modifications
    data_eligos = data_eligos[data_eligos['ref'] == "A"]
    
    # Filter for sufficient coverage
    data_eligos = data_eligos[data_eligos['total_reads'] > 20]
    
    # Select required columns
    eligos = data_eligos[['chrom', 'start_loc', 'end_loc', 'strand', 'adjPval', 'oddR', 'ESB_test']].copy()
    
    # Adjust positions
    eligos['start_loc'] += 1
    
    # Set modification status
    eligos['ESB_test'] = np.where(
        (~np.isnan(eligos['adjPval'])) & (~pd.isna(eligos['adjPval'])) &
        (~np.isnan(eligos['oddR'])) & (~pd.isna(eligos['oddR'])) &
        (eligos['adjPval'] < padj_threshold) & (eligos['oddR'] > oddr_threshold),
        "Mod", "Unmod"
    )
    
    # Filter for modifications
    eligos = eligos[eligos['ESB_test'] == "Mod"]
    
    # Select final columns
    eligos = eligos[['chrom', 'start_loc', 'end_loc', 'strand', 'ESB_test', 'adjPval']]
    
    # Rename columns
    eligos.columns = ['Chr', 'Start', 'End', 'Strand', 'Status', 'Padj']
    
    # Reorder columns
    eligos = eligos[['Chr', 'Start', 'End', 'Status', 'Padj', 'Strand']]
    
    # Adjust start position for BED format
    eligos['Start'] = eligos['Start'] - 1
    
    # Save results
    eligos.to_csv(output_file, sep='\t', index=False, header=True)
    
    print(f"ELIGOS2 processing complete. {len(eligos)} modifications detected.")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python postprocess_eligos2.py <input_file> <output_file> <padj_threshold> <oddr_threshold>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    padj_threshold = float(sys.argv[3])
    oddr_threshold = float(sys.argv[4])
    
    process_eligos2_results(input_file, output_file, padj_threshold, oddr_threshold)