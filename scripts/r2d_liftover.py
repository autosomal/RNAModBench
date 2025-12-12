#!/usr/bin/env python3
"""
R2Dtool liftover wrapper script.
This script coordinates liftover operations using R2Dtool.
"""

import subprocess
import sys
import os
import argparse
from pathlib import Path

def run_r2d_liftover(input_file, gtf_file, output_file, r2d_path="r2d"):
    """
    Run R2Dtool liftover command.
    
    Parameters:
    -----------
    input_file : str
        Path to input BED/TXT file
    gtf_file : str
        Path to reference GTF file
    output_file : str
        Path to output file
    r2d_path : str
        Path to R2Dtool executable (default: "r2d")
    """
    
    # Check if R2Dtool is available
    try:
        subprocess.run([r2d_path, "--help"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(f"R2Dtool not found at {r2d_path}", file=sys.stderr)
        sys.exit(1)
    
    # Run liftover command
    cmd = [
        r2d_path, "liftover",
        "-H",  # Include header
        "-g", gtf_file,
        "-i", input_file
    ]
    
    try:
        with open(output_file, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True, check=True)
        
        print(f"Liftover completed successfully.")
        print(f"Input: {input_file}")
        print(f"Output: {output_file}")
        
        # Check if output file was created
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            line_count = sum(1 for line in open(output_file) if line.strip())
            print(f"Lifted over {line_count} entries")
        else:
            print(f"Warning: Output file is empty or not created")
            
    except subprocess.CalledProcessError as e:
        print(f"R2Dtool failed with error: {e}", file=sys.stderr)
        print(f"Error output: {e.stderr}", file=sys.stderr)
        sys.exit(1)

def batch_liftover(input_dir, gtf_file, output_dir, r2d_path="r2d"):
    """
    Run liftover for multiple files in a directory.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing input files
    gtf_file : str
        Path to reference GTF file
    output_dir : str
        Directory for output files
    r2d_path : str
        Path to R2Dtool executable
    """
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all .txt and .bed files
    input_files = []
    for ext in ['.txt', '.bed', '.tsv']:
        input_files.extend(Path(input_dir).glob(f"*{ext}"))
    
    if not input_files:
        print(f"No input files found in {input_dir}")
        return
    
    print(f"Found {len(input_files)} files to process")
    
    # Process each file
    for input_file in input_files:
        output_file = os.path.join(output_dir, f"{input_file.stem}_liftover{input_file.suffix}")
        print(f"\nProcessing {input_file.name}...")
        
        try:
            run_r2d_liftover(str(input_file), gtf_file, output_file, r2d_path)
        except Exception as e:
            print(f"Failed to process {input_file.name}: {e}")
            continue

def main():
    parser = argparse.ArgumentParser(description='R2Dtool liftover wrapper')
    parser.add_argument('-i', '--input', required=True, 
                       help='Input file or directory')
    parser.add_argument('-g', '--gtf', required=True,
                       help='Reference GTF file')
    parser.add_argument('-o', '--output', required=True,
                       help='Output file or directory')
    parser.add_argument('--r2d-path', default='r2d',
                       help='Path to R2Dtool executable (default: r2d)')
    parser.add_argument('--batch', action='store_true',
                       help='Process all files in input directory')
    
    args = parser.parse_args()
    
    if args.batch:
        # Batch processing
        batch_liftover(args.input, args.gtf, args.output, args.r2d_path)
    else:
        # Single file processing
        run_r2d_liftover(args.input, args.gtf, args.output, args.r2d_path)

if __name__ == "__main__":
    main()