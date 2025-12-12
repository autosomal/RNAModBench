#!/usr/bin/env python3
"""
Post-process m6Anet results to generate standardized output format.
Converts m6Anet predictions to BED-like format with filtering.
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path
import glob

def process_m6anet_results(input_dir, output_file, prob_threshold=0.5, ratio_threshold=0.1):
    """
    Process m6Anet results and convert to standardized format.
    
    Parameters:
    -----------
    input_dir : str
        Path to m6Anet inference output directory
    output_file : str
        Path to processed output file
    prob_threshold : float
        Probability threshold for filtering (default: 0.5)
    ratio_threshold : float
        Modification ratio threshold for filtering (default: 0.1)
    """
    
    # Find the data.site_proba.csv file
    proba_file = os.path.join(input_dir, "data.site_proba.csv")
    if not os.path.exists(proba_file):
        # Try to find it in subdirectories
        proba_files = glob.glob(os.path.join(input_dir, "**", "data.site_proba.csv"), recursive=True)
        if not proba_files:
            print(f"Could not find data.site_proba.csv in {input_dir}", file=sys.stderr)
            sys.exit(1)
        proba_file = proba_files[0]
    
    # Read m6Anet results
    try:
        data_m6anet = pd.read_csv(proba_file, header=0, sep=",")
    except Exception as e:
        print(f"Error reading input file {proba_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Check minimum required columns
    if data_m6anet.shape[1] < 6:
        print(f"Insufficient columns in {proba_file}. Expected at least 6 columns.", file=sys.stderr)
        sys.exit(1)
    
    # Create standardized output format
    m6anet = pd.DataFrame({
        "TranscriptID": data_m6anet.iloc[:, 0],
        "Start": data_m6anet.iloc[:, 1] + 1,
        "End": data_m6anet.iloc[:, 1] + 1,
        "Status": data_m6anet.iloc[:, 3],
        "Prob": data_m6anet.iloc[:, 3],
        "mod_ratio": data_m6anet.iloc[:, 5]
    })
    
    # Set modification status
    m6anet['Status'] = np.where(
        (~np.isnan(m6anet['Status'])) & (~pd.isna(m6anet['Status'])) & 
        (m6anet['Status'] > prob_threshold),
        "Mod", "Unmod"
    )
    
    # Filter for modifications
    m6anet = m6anet[m6anet['Status'] == "Mod"]
    
    # Filter for minimum modification ratio
    m6anet = m6anet[m6anet['mod_ratio'] > ratio_threshold]
    
    # Add strand information
    m6anet['strand'] = "*"
    
    # Select final columns
    df_m6anet_final = m6anet[['TranscriptID', 'Start', 'End', 'strand', 'Status', 'Prob', 'mod_ratio']]
    
    # Rename columns
    df_m6anet_final.columns = ['Chr', 'Start', 'End', 'Strand', 'Status', 'Prob', 'mod_ratio']
    
    # Reorder columns
    df_m6anet_final = df_m6anet_final[['Chr', 'Start', 'End', 'Status', 'Prob', 'Strand', 'mod_ratio']]
    
    # Adjust start position for BED format
    df_m6anet_final['Start'] = df_m6anet_final['Start'] - 1
    
    # Save results
    df_m6anet_final.to_csv(output_file, sep='\t', index=False, header=True)
    
    print(f"m6Anet processing complete. {len(df_m6anet_final)} modifications detected.")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python postprocess_m6anet.py <input_dir> <output_file> <prob_threshold> <ratio_threshold>")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    prob_threshold = float(sys.argv[3])
    ratio_threshold = float(sys.argv[4])
    
    process_m6anet_results(input_dir, output_file, prob_threshold, ratio_threshold)