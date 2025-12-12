#!/usr/bin/env python3
"""
Post-process Nanocompore results to generate standardized output format.
Converts Nanocompore predictions to BED-like format with filtering.
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path
import glob

def process_nanocompore_results(input_dir, output_file, pvalue_threshold=0.05, lor_threshold=0.5):
    """
    Process Nanocompore results and convert to standardized format.
    
    Parameters:
    -----------
    input_dir : str
        Path to Nanocompore results directory
    output_file : str
        Path to processed output file
    pvalue_threshold : float
        P-value threshold for filtering (default: 0.05)
    lor_threshold : float
        Log odds ratio threshold for filtering (default: 0.5)
    """
    
    # Find the nanocompore_results.tsv file
    result_file = os.path.join(input_dir, "outnanocompore_results.tsv")
    if not os.path.exists(result_file):
        # Try alternative names
        alt_files = ["nanocompore_results.tsv", "results.tsv", "out.tsv"]
        for alt_file in alt_files:
            test_file = os.path.join(input_dir, alt_file)
            if os.path.exists(test_file):
                result_file = test_file
                break
        else:
            # Try to find any tsv file
            tsv_files = glob.glob(os.path.join(input_dir, "*.tsv"))
            if not tsv_files:
                print(f"Could not find Nanocompore results in {input_dir}", file=sys.stderr)
                sys.exit(1)
            result_file = tsv_files[0]
    
    # Read Nanocompore results
    try:
        data_nanocompore = pd.read_csv(result_file, header=0, sep='\t')
    except Exception as e:
        print(f"Error reading input file {result_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Check required columns
    required_cols = ['ref_kmer', 'ref_id', 'pos', 'Logit_LOR', 'GMM_logit_pvalue', 'strand']
    if not all(col in data_nanocompore.columns for col in required_cols):
        print(f"Missing required columns in {result_file}. Expected: {required_cols}", file=sys.stderr)
        sys.exit(1)
    
    # Filter for adenosine modifications (A in the middle of 5-mer)
    data_nanocompore = data_nanocompore[data_nanocompore['ref_kmer'].str[2] == 'A']
    
    # Select required columns
    nanocompore = data_nanocompore[['ref_id', 'pos', 'Logit_LOR', 'GMM_logit_pvalue', 'strand']].copy()
    
    # Adjust position
    nanocompore['pos'] += 3
    nanocompore['pos_end'] = nanocompore['pos']
    
    # Handle Logit_LOR (replace 'NC' with NaN)
    nanocompore['Logit_LOR'] = nanocompore['Logit_LOR'].replace('NC', np.nan)
    
    # Convert to numeric
    nanocompore['Logit_LOR'] = pd.to_numeric(nanocompore['Logit_LOR'], errors='coerce')
    
    # Set modification status
    nanocompore['Status'] = np.where(
        (~np.isnan(nanocompore['GMM_logit_pvalue'])) & (~pd.isna(nanocompore['GMM_logit_pvalue'])) &
        (~np.isnan(nanocompore['Logit_LOR'])) & (~pd.isna(nanocompore['Logit_LOR'])) &
        (nanocompore['GMM_logit_pvalue'] < pvalue_threshold) &
        (np.abs(nanocompore['Logit_LOR']) > lor_threshold),
        "Mod", "Unmod"
    )
    
    # Filter for modifications
    nanocompore = nanocompore[nanocompore['Status'] == "Mod"]
    
    # Select final columns
    nanocompore = nanocompore[['ref_id', 'pos', 'pos_end', 'strand', 'Status', 'GMM_logit_pvalue']]
    
    # Rename columns
    nanocompore.columns = ['Chr', 'Start', 'End', 'Strand', 'Status', 'Pvalue']
    
    # Reorder columns
    nanocompore = nanocompore[['Chr', 'Start', 'End', 'Status', 'Pvalue', 'Strand']]
    
    # Adjust start position for BED format
    nanocompore['Start'] = nanocompore['Start'] - 1
    
    # Save results
    nanocompore.to_csv(output_file, sep='\t', index=False, header=True)
    
    print(f"Nanocompore processing complete. {len(nanocompore)} modifications detected.")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python postprocess_nanocompore.py <input_dir> <output_file> <pvalue_threshold> <lor_threshold>")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    pvalue_threshold = float(sys.argv[3])
    lor_threshold = float(sys.argv[4])
    
    process_nanocompore_results(input_dir, output_file, pvalue_threshold, lor_threshold)