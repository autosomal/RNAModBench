#!/usr/bin/env python3
"""
Post-process DRUMMER results to generate standardized output format.
Converts DRUMMER predictions to BED-like format with filtering.
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path

def process_drummer_results(input_file, output_file, pvalue_threshold=0.05, frac_diff_threshold=0.1):
    """
    Process DRUMMER results and convert to standardized format.
    
    Parameters:
    -----------
    input_file : str
        Path to DRUMMER summary.txt file
    output_file : str
        Path to processed output file
    pvalue_threshold : float
        P-value threshold for filtering (default: 0.05)
    frac_diff_threshold : float
        Fraction difference threshold (default: 0.1)
    """
    
    # Read DRUMMER summary results
    try:
        table = pd.read_csv(input_file, header=0, sep="\t")
    except Exception as e:
        print(f"Error reading input file {input_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Filter for adenosine modifications (reference base = "A")
    data_drummer = table[table['reference_base'] == "A"]
    
    # Filter for sufficient coverage
    data_drummer = data_drummer[data_drummer['depth_ctrl'] > 20]
    
    # Select required columns
    drummer = data_drummer[['transcript_id', 'transcript_pos', 'OR_padj', 'G_padj', 'frac_diff']]
    
    # Set modification status based on p-value
    drummer['G_padj'] = np.where(drummer['OR_padj'] < pvalue_threshold, "Mod", "Unmod")
    
    # Filter for significant fraction difference
    drummer = drummer[np.abs(drummer['frac_diff']) > frac_diff_threshold]
    
    # Add strand information
    drummer['strand'] = "*"
    drummer['End'] = drummer['transcript_pos']
    
    # Select and rename columns
    drummer = drummer[['transcript_id', 'transcript_pos', 'End', 'strand', 'G_padj', 'OR_padj', 'frac_diff']]
    drummer.columns = ['Chr', 'Start', 'End', 'Strand', 'Status', 'Pvalue', 'FracDiff']
    
    # Reorder columns
    drummer = drummer[['Chr', 'Start', 'End', 'Status', 'Pvalue', 'Strand', 'FracDiff']]
    
    # Adjust start position for BED format
    drummer['Start'] = drummer['Start'] - 1
    
    # Save results
    drummer.to_csv(output_file, sep='\t', index=False, header=True)
    
    print(f"DRUMMER processing complete. {len(drummer)} modifications detected.")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python postprocess_drummer.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_drummer_results(input_file, output_file)