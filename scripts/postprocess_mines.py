#!/usr/bin/env python3
"""
Post-process MINES results to generate standardized output format.
Converts MINES predictions to BED-like format with filtering.
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path

def process_mines_results(input_file, output_file, coverage_threshold=20, ratio_threshold=0.1):
    """
    Process MINES results and convert to standardized format.
    
    Parameters:
    -----------
    input_file : str
        Path to MINES raw output file (BED format)
    output_file : str
        Path to processed output file
    coverage_threshold : int
        Minimum coverage threshold (default: 20)
    ratio_threshold : float
        Modification ratio threshold for filtering (default: 0.1)
    """
    
    # Read MINES results
    try:
        # Read BED file assuming no header
        data_mines = pd.read_csv(input_file, sep="\t", header=None)
    except Exception as e:
        print(f"Error reading input file {input_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Check minimum required columns
    if data_mines.shape[1] < 8:
        print(f"Insufficient columns in {input_file}. Expected at least 8 columns.", file=sys.stderr)
        sys.exit(1)
    
    # Set column names
    data_mines.columns = [
        "Chr", "Start", "End", "5-mer", "unique_key", "Strand", "mod_ratio", "coverage"
    ] + [f"extra_{i}" for i in range(data_mines.shape[1] - 8)]
    
    # Convert mod_ratio and coverage to numeric
    data_mines["mod_ratio"] = pd.to_numeric(data_mines["mod_ratio"], errors='coerce')
    data_mines["coverage"] = pd.to_numeric(data_mines["coverage"], errors='coerce')
    
    # Drop rows with conversion failures
    data_mines = data_mines.dropna(subset=["mod_ratio", "coverage"])
    
    # Apply filtering criteria
    data_mines = data_mines[data_mines["coverage"] >= coverage_threshold]
    data_mines = data_mines[data_mines["mod_ratio"] > ratio_threshold]
    
    # Adjust End position
    data_mines["End"] = data_mines["Start"]
    
    # Add Status column
    data_mines["Status"] = "Mod"
    
    # Select final columns
    mines = data_mines[["Chr", "Start", "End", "Status", "mod_ratio", "Strand"]]
    
    # Adjust start position for BED format
    mines["Start"] = mines["Start"] - 1
    
    # Save results
    mines.to_csv(output_file, sep='\t', index=False, header=True, quoting=3)
    
    print(f"MINES processing complete. {len(mines)} modifications detected.")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python postprocess_mines.py <input_file> <output_file> <coverage_threshold> <ratio_threshold>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    coverage_threshold = int(sys.argv[3])
    ratio_threshold = float(sys.argv[4])
    
    process_mines_results(input_file, output_file, coverage_threshold, ratio_threshold)