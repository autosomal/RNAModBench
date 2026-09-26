"""Transcript-strand imputation for callsets that carry no strand column.

Background (2026-09-18)
-----------------------
Several tools report a position without a strand -- Nanom6A writes ``'*'`` for
**all** of its rows (25 callsets, ~286 k rows).  ``common/annotate.py`` treats
anything that is not ``'-'`` as ``'+'``, so the minus-strand half of those rows
(the modified adenosine shows up as a genomic ``T``) were annotated with a
spurious ``dist_center`` (-5..+5) and a ``pos_center`` 2-5 bp away from the real
site.  The calls themselves are correct: reverse-complementing their 5mer lifts
the DRACH share from 0.000 to 0.734 (Arabidopsis) / 0.998 (HeLa), exactly the
share the plus-strand rows show.

This module recovers the strand from two independent evidence layers:

1. the tool's own read-alignment BED (``config.READ_STRAND_BED``) -- the strand
   of the reads that produced the calls.  Decisive and sample-specific; where
   two genes overlap in antisense the *majority* strand of the covering reads is
   taken.  A wrong majority is self-limiting: the centre-base check then rejects
   the row, which is what would have happened to an unresolved row anyway, so
   taking the majority can only recover calls, never fake them.
2. the species exon annotation (``config.GTF_EXON``, the one the candidate
   universe is built from) -- used for the rows the read layer cannot settle,
   and marked ``ambiguous`` where exons of both strands overlap.

Per-position stranded depth
---------------------------
Merging intervals would destroy resolution (in a bacterium every block merges
into one genome-spanning segment).  Instead the depth of a position ``p`` is
counted with two binary searches per strand::

    depth(p) = #{start <= p} - #{end < p}
             = searchsorted(starts, p, 'right') - searchsorted(ends, p, 'left')

which is exact, vectorised and needs no interval algebra.

Coordinates are 0-based; the exon files are 1-based inclusive and the BED files
0-based half-open, each converted at the single point where it is read.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

from .config import GTF_EXON
from .match import fix_chromosome

#: imputation status codes
NO_EXON = 0
PLUS = 1
MINUS = 2
BOTH = 3

STATUS_LABEL = {NO_EXON: "no_exon", PLUS: "+", MINUS: "-", BOTH: "ambiguous"}

#: chrom -> (plus_starts, plus_ends, minus_starts, minus_ends); ends inclusive.
StrandIndex = dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]

_EXON_CACHE: dict[str, StrandIndex] = {}


def _iter_exons(gtf_path: Path):
    """Yield ``(contig, start_1based, end_1based, strand)`` for every exon."""
    p = Path(gtf_path)
    if p.suffix == ".exon" or str(p).endswith(".gtf.exon"):
        with open(p) as fh:
            rdr = csv.reader(fh)
            try:
                next(rdr)  # header
            except StopIteration:
                return
            for row in rdr:
                if len(row) < 8 or row[2] != "exon":
                    continue
                yield row[0], int(row[3]), int(row[4]), row[6]
        return
    with open(p) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] != "exon":
                continue
            yield f[0], int(f[3]), int(f[4]), f[6]


def _build(by_chrom: dict[str, list]) -> StrandIndex:
    out: StrandIndex = {}
    for chrom, iv in by_chrom.items():
        ps = np.array(sorted(s for s, e, b in iv if b == PLUS), dtype=np.int64)
        pe = np.array(sorted(e for s, e, b in iv if b == PLUS), dtype=np.int64)
        ms = np.array(sorted(s for s, e, b in iv if b == MINUS), dtype=np.int64)
        me = np.array(sorted(e for s, e, b in iv if b == MINUS), dtype=np.int64)
        out[chrom] = (ps, pe, ms, me)
    return out


def load_exon_index(species: str) -> StrandIndex:
    """Exon intervals (0-based, inclusive ends) per chromosome and strand."""
    if species in _EXON_CACHE:
        return _EXON_CACHE[species]
    gtf = GTF_EXON.get(species)
    by_chrom: dict[str, list] = defaultdict(list)
    if gtf is not None and Path(gtf).exists():
        for contig, start, end, strand in _iter_exons(gtf):
            bit = PLUS if strand == "+" else (MINUS if strand == "-" else 0)
            if bit == 0:
                continue
            by_chrom[fix_chromosome(contig)].append((start - 1, end - 1, bit))
    _EXON_CACHE[species] = _build(by_chrom)
    return _EXON_CACHE[species]


def load_alignment_index(path) -> StrandIndex:
    """Read-alignment intervals from a tool's own BED (e.g. ``extract.bed12``).

    Columns ``chrom, start(0-based), end(exclusive), ..., strand``.  BED12
    blocks (columns 11/12) are irrelevant: every block of a read carries the
    read's strand.
    """
    p = Path(path)
    by_chrom: dict[str, list] = defaultdict(list)
    if p.exists():
        with open(p) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 6 or f[0].startswith(("#", "track", "browser")):
                    continue
                bit = PLUS if f[5] == "+" else (MINUS if f[5] == "-" else 0)
                if bit == 0:
                    continue
                try:
                    start, end = int(f[1]), int(f[2]) - 1
                except ValueError:
                    continue
                if end < start:
                    continue
                by_chrom[fix_chromosome(f[0])].append((start, end, bit))
    return _build(by_chrom)


def depths(index: StrandIndex, chrom: np.ndarray, pos: np.ndarray
           ) -> tuple[np.ndarray, np.ndarray]:
    """``(n_plus, n_minus)`` covering intervals per position (exact, vectorised)."""
    chrom = np.asarray(chrom).astype(str)
    pos = np.asarray(pos, dtype=np.int64)
    nplus = np.zeros(len(pos), dtype=np.int64)
    nminus = np.zeros(len(pos), dtype=np.int64)
    if not index:
        return nplus, nminus
    for c in np.unique(chrom):
        seg = index.get(fix_chromosome(c))
        if seg is None:
            continue
        ps, pe, ms, me = seg
        sel = np.flatnonzero(chrom == c)
        p = pos[sel]
        if len(ps):
            nplus[sel] = (np.searchsorted(ps, p, side="right")
                          - np.searchsorted(pe, p, side="left"))
        if len(ms):
            nminus[sel] = (np.searchsorted(ms, p, side="right")
                           - np.searchsorted(me, p, side="left"))
    return nplus, nminus


def codes_from_depths(nplus: np.ndarray, nminus: np.ndarray,
                      majority: bool = False) -> np.ndarray:
    """Status per row from the two depths.

    ``majority=False`` (exon layer): a position covered by exons of both strands
    is reported as ``BOTH`` -- gene models have no meaningful "count".
    ``majority=True`` (read layer): the strand with more covering reads wins;
    an exact tie stays ``BOTH``.
    """
    code = np.zeros(len(nplus), dtype=np.int8)
    if majority:
        code[(nplus > nminus)] = PLUS
        code[(nminus > nplus)] = MINUS
        code[(nplus == nminus) & (nplus > 0)] = BOTH
        return code
    only_plus = (nplus > 0) & (nminus == 0)
    only_minus = (nminus > 0) & (nplus == 0)
    code[only_plus] = PLUS
    code[only_minus] = MINUS
    code[(nplus > 0) & (nminus > 0)] = BOTH
    return code


def strand_codes(chrom: np.ndarray, pos: np.ndarray, species: str) -> np.ndarray:
    """Exon-annotation status codes for positions (0/1/2/3)."""
    nplus, nminus = depths(load_exon_index(species), chrom, pos)
    return codes_from_depths(nplus, nminus, majority=False)


def impute_frame(df, species: str, strand_col: str = "strand",
                 alignment_index: StrandIndex | None = None
                 ) -> tuple[np.ndarray, np.ndarray, dict]:
    """``(strand_values, status, stats)`` for a callset frame.

    ``strand_values`` keeps the original value for rows that already carry
    ``+``/``-``; the rest get ``'+'``/``'-'`` when a layer is decisive and stay
    ``'*'`` otherwise.  ``status`` is the code of the layer that decided
    (``-1`` = untouched) and ``stats`` counts the outcome per layer.
    """
    cur = (df[strand_col].astype(str).to_numpy() if strand_col in df.columns
           else np.array([""] * len(df)))
    known = np.isin(cur, ["+", "-"])
    out_strand = cur.copy()
    status = np.full(len(df), -1, dtype=np.int8)
    stats = {"n_read_plus": 0, "n_read_minus": 0, "n_read_both": 0,
             "n_exon_plus": 0, "n_exon_minus": 0,
             "n_ambiguous": 0, "n_no_exon": 0}
    todo = ~known
    if todo.any():
        chrom = df["chrom"].astype(str).to_numpy()[todo]
        pos = df["pos_raw"].to_numpy()[todo]
        code = np.zeros(todo.sum(), dtype=np.int8)  # 0 = undecided so far
        if alignment_index:
            rp, rm = depths(alignment_index, chrom, pos)
            rcode = codes_from_depths(rp, rm, majority=True)
            stats["n_read_plus"] = int((rcode == PLUS).sum())
            stats["n_read_minus"] = int((rcode == MINUS).sum())
            stats["n_read_both"] = int((rcode == BOTH).sum())
            code = rcode
        undecided = code == 0
        if undecided.any():
            ep, em = depths(load_exon_index(species), chrom[undecided], pos[undecided])
            ecode = codes_from_depths(ep, em, majority=False)
            stats["n_exon_plus"] = int((ecode == PLUS).sum())
            stats["n_exon_minus"] = int((ecode == MINUS).sum())
            stats["n_ambiguous"] = int((ecode == BOTH).sum())
            stats["n_no_exon"] = int((ecode == NO_EXON).sum())
            code[undecided] = ecode
        new = np.array(["*"] * todo.sum(), dtype=object)
        new[code == PLUS] = "+"
        new[code == MINUS] = "-"
        out_strand[todo] = new
        status[todo] = code
    return out_strand, status, stats
