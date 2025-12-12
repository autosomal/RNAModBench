#!/usr/bin/env python3
"""
Extract 5-mer context sequences from genomic positions.
This script adds 5-mer sequence context to modification calls.
"""

import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq


def extract_5mer_context(chr_name, position, chr_sequences, strand='+'):
    """
    Extract 5-mer sequence context around a genomic position.
    
    Parameters:
    -----------
    chr_name : str
        Chromosome name
    position : int
        Genomic position (1-based)
    chr_sequences : dict
        Dictionary of chromosome sequences
    strand : str
        Strand information ('+' or '-')
    
    Returns:
    --------
    str
        5-mer sequence context
    """
    
    # Handle chromosome name variations
    seq_key = chr_name
    if seq_key not in chr_sequences:
        # Try different prefix combinations
        candidates = [seq_key, f"chr{seq_key}", seq_key.lstrip('0')]
        for candidate in candidates:
            if candidate in chr_sequences:
                seq_key = candidate
                break
        else:
            return 'NNNNN'
    
    # Get chromosome sequence
    sequence = chr_sequences[seq_key]
    if not sequence:
        return 'NNNNN'
    
    # Calculate target region (center ± 2)
    target_start = position - 2
    target_end = position + 2
    seq_len = len(sequence)
    
    # Handle boundary conditions
    start = max(target_start, 1)
    end = min(target_end, seq_len)
    
    if start > end:
        return 'NNNNN'
    
    # Extract sequence
    extracted = sequence[start-1:end]
    
    # Pad to 5 bp and convert to uppercase
    extracted = extracted.ljust(5, 'N')[:5].upper()
    
    # Reverse complement if on negative strand
    if strand == '-':
        seq_obj = Seq(extracted)
        extracted = str(seq_obj.reverse_complement())
    
    return extracted

    Extract 5-mer sequence context around a genomic position.
    
    Parameters:
    -----------
    chr_name : str
        Chromosome name
    position : int
        Genomic position (1-based)
    chr_sequences : dict
        Dictionary of chromosome sequences
    strand : str
        Strand information ('+' or '-')
    
    Returns:
    --------
    str
        5-mer sequence context
    """
    
    # Handle chromosome name variations
    seq_key = chr_name
    if seq_key not in chr_sequences:
        # Try different prefix combinations
        candidates = [seq_key, f"chr{seq_key}", seq_key.lstrip('0')]
        for candidate in candidates:
            if candidate in chr_sequences:
                seq_key = candidate
                break
        else:
            return 'NNNNN'
    
    # Get chromosome sequence
    sequence = chr_sequences[seq_key]
    if not sequence:
        return 'NNNNN'
    
    # Calculate target region (center ± 2)
    target_start = position - 2
    target_end = position + 2
    seq_len = len(sequence)
    
    # Handle boundary conditions
    start = max(target_start, 1)
    end = min(target_end, seq_len)
    
    if start > end:
        return 'NNNNN'
    
    # Extract sequence
    extracted = sequence[start-1:end]
    
    # Pad to 5 bp and convert to uppercase
    extracted = extracted.ljust(5, 'N')[:5].upper()
    
    # Reverse complement if on negative strand
    if strand == '-':
        seq_obj = Seq(extracted)
        extracted = str(seq_obj.reverse_complement())
    
    return extracted

def extract_5mer_context(chr_name, position, chr_sequences, strand='+'):
    """
    Extract 5-mer sequence context around a genomic position.
    
    Parameters:
    -----------
    chr_name : str
        Chromosome name
    position : int
        Genomic position (1-based)
    chr_sequences : dict
        Dictionary of chromosome sequences
    strand : str
        Strand information ('+' or '-')
    
    Returns:
    --------
    str
        5-mer sequence context
    """
    
    # Handle chromosome name variations
    seq_key = chr_name
    if seq_key not in chr_sequences:
        # Try different prefix combinations
        candidates = [seq_key, f"chr{seq_key}", seq_key.lstrip('0')]
        for candidate in candidates:
            if candidate in chr_sequences:
                seq_key = candidate
                break
        else:
            return 'NNNNN'
    
    # Get chromosome sequence
    sequence = chr_sequences[seq_key]
    if not sequence:
        return 'NNNNN'
    
    # Calculate target region (center ± 2)
    target_start = position - 2
    target_end = position + 2
    seq_len = len(sequence)
    
    # Handle boundary conditions
    start = max(target_start, 1)
    end = min(target_end, seq_len)
    
    if start > end:
        return 'NNNNN'
    
    # Extract sequence
    extracted = sequence[start-1:end]
    
    # Pad to 5 bp and convert to uppercase
    extracted = extracted.ljust(5, 'N')[:5].upper()
    
    # Reverse complement if on negative strand
    if strand == '-':
        seq_obj = Seq(extracted)
        extracted = str(seq_obj.reverse_complement())
    
    return extracted

def main(input_file, fasta_file, output_file=None):
    """
    Main function to process input file and add 5-mer context.
    
    Parameters:
    -----------
    input_file : str
        Path to input BED/TSV file
    fasta_file : str
        Path to reference FASTA file
    output_file : str
        Path to output file (default: stdout)
    """
    
    # Read FASTA file and build chromosome sequence dictionary
    chr_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Remove possible 'chr' prefix and convert to uppercase
        chr_key = record.id.lstrip('chr').upper()
        chr_sequences[chr_key] = str(record.seq)
    
    # Determine output destination
    if output_file and output_file != 'stdout':
        out_f = open(output_file, 'w')
    else:
        out_f = sys.stdout
    
    try:
        with open(input_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            
            # Add 5mer column to header
            new_header = header + ['5mer']
            print('\t'.join(new_header), file=out_f)
            
            for row in reader:
                # Parse input row
                chr_name = row[0].strip().upper()
                pos = int(row[2]) if len(row) > 2 else int(row[1])  # Handle different formats
                
                # Extract 5-mer context
                strand = row[5] if len(row) > 5 else '*'
                fivemer = extract_5mer_context(chr_name, pos, chr_sequences, strand)
                
                # Output result
                new_row = row + [fivemer]
                print('\t'.join(new_row), file=out_f)
    
    finally:
        if output_file and output_file != 'stdout':
            out_f.close()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_5mer.py <input_file> <reference.fasta> [output_file]")
        print("If output_file is not specified, results are printed to stdout")
        sys.exit(1)
    
    input_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    main(input_file, fasta_file, output_file)
