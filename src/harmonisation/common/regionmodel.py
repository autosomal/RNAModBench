"""Transcript region model for metagene (GUITAR-style) analysis.

Segments every transcript into ``5'UTR / CDS / 3'UTR`` (protein-coding) or a
single ``ncRNA body``, and adds the ``1 kb`` windows upstream of the TSS and
downstream of the TES, all stored as genomic 0-based half-open intervals.  A
site is then placed on the **transcript-direction normalised coordinate** of the
segment it falls in -- the x-axis of a GUITAR metagene plot
(``1kb | 5'UTR | CDS | 3'UTR | 1kb``).

Differences from the published R workflow (``code/guitar/*.r``), all deliberate:

* the legacy script forced ``strand <- "+"`` for every BED line, so minus-strand
  sites were read on the plus-strand coordinate; here the site's own strand is
  honoured (``strand_aware=False`` reproduces the legacy behaviour);
* a site overlapping several transcripts is assigned to the **longest** one
  (one record per site, deterministic) rather than counted once per overlapping
  transcript -- ``tx_rule="all"`` gives the multiplicity-weighted variant.
"""

from __future__ import annotations

# --- RNAModBench path bootstrap (added when this file was deposited) ----------
import os as _rb_os, pathlib as _rb_pl


def _rb_find(start):
    for p in (start, *start.parents):
        if (p / "RNAMOD_BENCH_ROOT").exists():
            return p
    return start


_RB = _rb_pl.Path(_rb_os.environ.get("RNAMODBENCH_ROOT") or _rb_find(_rb_pl.Path(__file__).resolve().parent))
_XB = _rb_pl.Path(_rb_os.environ.get("RNAMODBENCH_LOCAL") or (_RB / "_local"))
# --------------------------------------------------------------------------- #
import pickle
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .config import SITES_ROOT
from .match import fix_chromosome

MODEL_DIR = (_XB / "reference/regionmodels")

UTR5, CDS, UTR3, F5, F3, NCR = range(6)
KINDS = {UTR5: "five_prime_UTR", CDS: "CDS", UTR3: "three_prime_UTR",
         F5: "five_prime_flank", F3: "three_prime_flank", NCR: "ncRNA_body"}
#: plot order of the metagene x-axis (transcript direction)
SEGMENT_ORDER = [F5, UTR5, CDS, UTR3, F3]
NCR_ORDER = [F5, NCR, F3]
FLANK = 1000


@dataclass
class RegionIndex:
    """Interval tables keyed by ``(chrom, strand)``."""

    tables: dict[tuple[str, str], dict[str, np.ndarray]]
    n_tx: int = 0

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("wb") as fh:
            pickle.dump(self, fh, protocol=4)

    @staticmethod
    def load(path: Path) -> "RegionIndex":
        with path.open("rb") as fh:
            return pickle.load(fh)

    def n_segments(self) -> int:
        return sum(len(t["start"]) for t in self.tables.values())


# --------------------------------------------------------------------------- #
# GTF parsing
# --------------------------------------------------------------------------- #
def _attrs(field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for kv in field.split(";"):
        kv = kv.strip()
        if not kv:
            continue
        if " " in kv:
            k, v = kv.split(" ", 1)
            out[k] = v.strip().strip('"')
        elif "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def parse_gtf(gtf: Path, *, tx_class: str = "mrna",
              flank: int = FLANK) -> RegionIndex:
    """Build the region index for one transcript class.

    ``tx_class``: ``"mrna"`` (protein-coding) or ``"ncrna"`` (everything else).
    """
    exons: dict[str, list[tuple[int, int]]] = {}
    cdss: dict[str, list[tuple[int, int]]] = {}
    meta: dict[str, tuple[str, str]] = {}

    with gtf.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS"):
                continue
            a = _attrs(f[8])
            tid = a.get("transcript_id") or a.get("Parent", "").split(":")[-1]
            if not tid:
                continue
            coding = any(a.get(k) == "protein_coding" for k in
                         ("transcript_biotype", "biotype",
                          "transcript_type", "gene_type"))
            if coding != (tx_class == "mrna"):
                continue
            meta.setdefault(tid, (fix_chromosome(f[0]), f[6]))
            if f[2] == "exon":
                exons.setdefault(tid, []).append((int(f[3]) - 1, int(f[4])))
            else:
                cdss.setdefault(tid, []).append((int(f[3]) - 1, int(f[4])))

    rows: dict[tuple[str, str], list[tuple[int, int, int, int, int, int]]] = {}
    for tid, (chrom, strand) in meta.items():
        blocks = sorted(exons.get(tid, []), key=lambda b: b[0],
                        reverse=(strand == "-"))
        blocks = [(s, e) for s, e in blocks if e > s]
        if not blocks:
            continue
        tx_len = sum(e - s for s, e in blocks)
        cds_tx = _cds_tx_range(blocks, strand, cdss.get(tid, []))

        segs: list[tuple[int, int, int, int, int]] = []   # gs, ge, kind, tx_a, len
        cum = 0
        for s, e in blocks:
            blen = e - s
            for ta, tb, kind in _classify(cum, blen, cds_tx):
                segs.append((_to_genomic(s, e, strand, cum, ta),
                             _to_genomic(s, e, strand, cum, tb),
                             kind, ta, tb - ta))
            cum += blen
        tss, tes = blocks[0][0], blocks[-1][1]
        if strand == "+":
            segs.append((max(0, tss - flank), tss, F5, 0, flank))
            segs.append((tes, tes + flank, F3, 0, flank))
        else:
            segs.append((tes, tes + flank, F3, 0, flank))
            segs.append((max(0, tss - flank), tss, F5, 0, flank))

        for gs, ge, kind, ta, slen in segs:
            if ge > gs and slen > 0:
                rows.setdefault((chrom, strand), []).append(
                    (gs, ge, kind, ta, slen, tx_len))

    tables: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    for key, lst in rows.items():
        arr = np.asarray(lst, dtype=np.int64)
        arr = arr[np.lexsort((arr[:, 1], arr[:, 0]))]
        tables[key] = {
            "start": np.ascontiguousarray(arr[:, 0]),
            "end": np.ascontiguousarray(arr[:, 1]),
            "kind": np.ascontiguousarray(arr[:, 2], dtype=np.int8),
            "tx_start": np.ascontiguousarray(arr[:, 3]),
            "seg_len": np.ascontiguousarray(arr[:, 4]),
            "tx_len": np.ascontiguousarray(arr[:, 5]),
        }
    return RegionIndex(tables, len(meta))


def _to_genomic(bs: int, be: int, strand: str, cum: int, tx: int) -> int:
    """Transcript coordinate -> genomic coordinate inside one exon block."""
    return bs + (tx - cum) if strand == "+" else be - (tx - cum)


def _cds_tx_range(blocks: list[tuple[int, int]], strand: str,
                  cds: list[tuple[int, int]]) -> tuple[int, int] | None:
    """Transcript-coordinate span ``[start, end)`` of the CDS, or None."""
    lo = hi = None
    for cs, ce in cds:
        cum = 0
        for bs, be in blocks:
            a, b = max(cs, bs), min(ce, be)
            if b > a:
                ta = (a - bs) if strand == "+" else (be - a)
                tb = (b - bs) if strand == "+" else (be - b)
                t0, t1 = cum + min(ta, tb), cum + max(ta, tb)
                lo = t0 if lo is None else min(lo, t0)
                hi = t1 if hi is None else max(hi, t1)
            cum += be - bs
    return None if lo is None else (lo, hi)


def _classify(cum: int, blen: int,
              cds_tx: tuple[int, int] | None) -> list[tuple[int, int, int]]:
    """Split one exon block (transcript coords ``[cum, cum+blen)``) into pieces."""
    a, b = cum, cum + blen
    if cds_tx is None:
        return [(a, b, NCR)]
    cs, ce = cds_tx
    out: list[tuple[int, int, int]] = []
    if a < cs:
        out.append((a, min(b, cs), UTR5))
    if min(b, ce) > max(a, cs):
        out.append((max(a, cs), min(b, ce), CDS))
    if b > ce:
        out.append((max(a, ce), b, UTR3))
    return out


# --------------------------------------------------------------------------- #
# site assignment
# --------------------------------------------------------------------------- #
def assign_sites(idx: RegionIndex, chrom: str, strand: str,
                 positions: np.ndarray, *, tx_rule: str = "longest",
                 max_overlap: int = 300):
    """Place ``positions`` (0-based) on the region model of one chrom/strand.

    Returns ``(site_index, kind, norm_pos, tx_len)``; ``site_index`` indexes
    ``positions`` for the records that fell inside some segment.
    """
    tab = idx.tables.get((chrom, strand))
    empty = (np.zeros(0, np.int64), np.zeros(0, np.int8),
             np.zeros(0, np.float64), np.zeros(0, np.int64))
    pos = np.asarray(positions, dtype=np.int64)
    if tab is None or pos.size == 0:
        return empty
    starts, ends = tab["start"], tab["end"]
    si, ki, nz, tl = [], [], [], []
    for i in np.argsort(pos, kind="stable"):
        p = int(pos[i])
        hi = int(np.searchsorted(starts, p, side="right"))
        lo = max(0, hi - max_overlap)
        sel = np.flatnonzero(ends[lo:hi] > p)
        if sel.size == 0:
            continue
        sel += lo
        if tx_rule == "longest":
            sel = sel[np.argmax(tab["tx_len"][sel]):sel.size][:1]
        for j in sel:
            j = int(j)
            si.append(i)
            ki.append(tab["kind"][j])
            local = p - starts[j] if strand == "+" else ends[j] - 1 - p
            nz.append((tab["tx_start"][j] + local + 0.5) / tab["seg_len"][j])
            tl.append(tab["tx_len"][j])
    if not si:
        return empty
    return (np.asarray(si, np.int64), np.asarray(ki, np.int8),
            np.asarray(nz, np.float64), np.asarray(tl, np.int64))
