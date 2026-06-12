#!/usr/bin/env python3
"""
Extract 5-mer context sequences from genomic/transcriptomic positions.
Given a bed-like file of modification positions, append 5-mer context
as a new column (center base ± 2 neighbors) from a reference FASTA.

[输入格式]
    任意 bed-like / tsv 文件，含至少前三列：
      col 0: 染色体/转录本 ID
      col 1: start（若文件 ≥3 列则用 col 1 替代 col 2 做位置）
      col 2: end（本脚本实际读取 col 2，如果存在，否则用 col 1）
      col 5: strand（可选，缺省用 '*'，负链时取反向互补 5-mer）
    分隔符: \t，含表头
[输出格式]
    与输入相同的列，在末尾追加一列 `5mer`。
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
        '+' / '-' / '*'. '-' 时返回反向互补序列，使 5-mer 始终
        对应修饰所在 reads 的 5'→3' 顺序。

    Returns:
    --------
    str
        5-bp 5-mer context (uppercase)。边界外位置返回 'NNNNN'。
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

    # center ± 2；position 是 1-based
    target_start = position - 2
    target_end = position + 2
    seq_len = len(sequence)

    # Handle boundary conditions — pad with 'N' 而非直接截断
    start = max(target_start, 1)
    end = min(target_end, seq_len)

    if start > end:
        return 'NNNNN'

    # Python 切片是 0-based [start, end)
    extracted = sequence[start - 1:end]

    # Left-pad with 'N' to ensure the k-mer is 5 bp
    if target_start < 1:
        extracted = 'N' * (1 - target_start) + extracted
    if target_end > seq_len:
        extracted = extracted + 'N' * (target_end - seq_len)

    extracted = extracted.upper()[:5]

    # 负链 → 反向互补：保证返回的 k-mer 是 reads 实际观察到的顺序
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

    # 一次性构建序列字典；同时存储原始 ID 和大写无 'chr' 的 ID，避免
    # 参考文件用 'chr1' 而输入文件用 '1' 这类不一致问题。
    chr_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_str = str(record.seq)
        # 原始 ID
        chr_sequences[record.id] = seq_str
        # 去掉 chr 前缀后的大写（作为 fallback key）
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
                # 优先使用 col 2 (end) 作为位置；否则用 col 1 (start)
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
        print("  input_file:       bed-like / tsv 至少有 Chr/Start 两列，含表头")
        print("  reference.fasta:  genome 或 transcriptome FASTA")
        print("  output_file:      可选；缺省则打印到 stdout")
        sys.exit(1)

    input_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else None

    main(input_file, fasta_file, output_file)
