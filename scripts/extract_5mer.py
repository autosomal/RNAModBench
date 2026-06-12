#!/usr/bin/env python3
"""
Extract 5-mer context sequences from genomic/transcriptomic positions.
Given a bed-like file of modification positions, append 5-mer context
as a new column (center base ± 2 neighbors) from a reference FASTA.

[Input Format]
    Any bed-like / tsv file with at least the first three columns:
      col 0: Chromosome / transcript ID
      col 1: start (if file has >=3 columns, col 1 is used instead of col 2 as position)
      col 2: end (this script actually reads col 2 when present, otherwise col 1)
      col 5: strand (optional, defaults to '*'; reverse-complement 5-mer on minus strand)
    Separator: \t, with header line
[Output Format]
    Same columns as input, with an additional `5mer` column appended at the end.
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq


def extract_5mer_context(chr_name, position, chr_sequences, strand='+'):
    """
    Extract 5-mer sequence context around a genomic/transcriptomic position.

    Parameters:
    -----------
    chr_name : str
        Chromosome/transcript ID. Matched case-insensitive after stripping 'chr'.
    position : int
        1-based position on the chromosome/transcript.
    chr_sequences : dict
        Pre-parsed dict: {<normalized_chrom_id>: <sequence string>}.
    strand : str
        '+' / '-' / '*'. On '-' strand the reverse complement is returned so the 5-mer
        always corresponds to the 5'->3' order of the reads bearing the modification.

    Returns:
    --------
    str
        5-bp 5-mer context (uppercase). Out-of-range positions return 'NNNNN'.
    """

    # Handle chromosome name variations ('chr1' <-> '1' <-> 'CHR1')
    seq_key = chr_name
    if seq_key not in chr_sequences:
        candidates = [seq_key, f"chr{seq_key}", seq_key.lstrip('0')]
        for candidate in candidates:
            if candidate in chr_sequences:
                seq_key = candidate
                break
        else:
            return 'NNNNN'

    sequence = chr_sequences[seq_key]
    if not sequence:
        return 'NNNNN'

    # center +/- 2; position is 1-based
    target_start = position - 2
    target_end = position + 2
    seq_len = len(sequence)

    # Handle boundary conditions - pad with 'N' rather than truncating directly
    start = max(target_start, 1)
    end = min(target_end, seq_len)

    if start > end:
        return 'NNNNN'

    # Python slicing is 0-based [start, end)
    extracted = sequence[start - 1:end]

    # Left-pad with 'N' to ensure the k-mer is 5 bp
    if target_start < 1:
        extracted = 'N' * (1 - target_start) + extracted
    if target_end > seq_len:
        extracted = extracted + 'N' * (target_end - seq_len)

    extracted = extracted.upper()[:5]

    # Minus strand -> reverse complement: ensure the returned k-mer matches the order actually observed in reads
    if strand == '-':
        extracted = str(Seq(extracted).reverse_complement())

    return extracted


def main(input_file, fasta_file, output_file=None):
    """
    Read bed-like input file and append the 5-mer column.

    Parameters:
    -----------
    input_file : str
        Path to input BED/TSV file (first line must be header).
    fasta_file : str
        Path to reference FASTA file (genomic or transcriptomic).
    output_file : str or None
        Path to output file. If None or "stdout", prints to stdout.
    """

    # Build the sequence dictionary in one pass; store both the original ID and
    # the uppercased ID with any 'chr' prefix stripped, to avoid mismatches such
    # as the reference using 'chr1' while the input uses '1'.
    chr_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_str = str(record.seq)
        # original ID
        chr_sequences[record.id] = seq_str
        # uppercase ID with the 'chr' prefix removed (used as a fallback key)
        key_no_chr = record.id.replace('chr', '', 1).upper() if record.id.lower().startswith('chr') else record.id.upper()
        if key_no_chr != record.id:
            chr_sequences.setdefault(key_no_chr, seq_str)

    if output_file and output_file != 'stdout':
        out_f = open(output_file, 'w')
    else:
        out_f = sys.stdout

    try:
        with open(input_file, 'r') as f:
            header = next(f).rstrip('\n\r')
            out_f.write(header + '\t5mer\n')

            for line in f:
                row = line.rstrip('\n\r').split('\t')
                if len(row) < 2:
                    continue

                chr_name = row[0].strip().upper()
                # Prefer col 2 (end) as position; otherwise use col 1 (start)
                pos_col = row[2].strip() if len(row) > 2 and row[2].strip() else row[1].strip()
                try:
                    pos = int(float(pos_col))
                except ValueError:
                    fivemer = 'NNNNN'
                else:
                    strand = row[5] if len(row) > 5 else '*'
                    fivemer = extract_5mer_context(chr_name, pos, chr_sequences, strand)

                out_f.write('\t'.join(row) + '\t' + fivemer + '\n')

    finally:
        if out_f is not sys.stdout:
            out_f.close()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_5mer.py <input_file> <reference.fasta> [output_file]")
        print("  input_file:       bed-like / tsv with at least Chr/Start columns, with header")
        print("  reference.fasta:  genome or transcriptome FASTA")
        print("  output_file:      optional; defaults to stdout when omitted")
        sys.exit(1)

    input_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else None

    main(input_file, fasta_file, output_file)
