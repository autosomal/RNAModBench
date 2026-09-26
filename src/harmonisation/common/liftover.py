"""Transcript -> genome liftover, implemented from the GTF (no external binary).

Why this exists
---------------
The legacy pipeline used ``r2d liftover``.  For Human / Mouse / E. coli that
reproduces the archived ``*_liftover.txt`` files exactly, but for **Arabidopsis
it does not**: ~46 % of rows come back at a *different* exon (checked on the
Nanocompore input: 137/255 rows differ, and the reference base under the r2d
coordinate is not the modified base any more).  A plain exon-arithmetic mapper
written here reproduces the archived files **row-for-row for every species**
(see ``02b_validate_liftover.py``), and it is 100 % deterministic - so the
rebuild uses this instead of shelling out to r2d.

Conventions (verified against the legacy ``*_liftover.txt`` files)
-----------------------------------------------------------------
* input: transcript id + transcript position, 0-based, exon-aware
* output: ``chromosome`` (GTF chromosome), ``start`` = 0-based genomic position,
  ``end`` = ``start + 1``; ``strand`` taken from the GTF
* transcripts missing from the model are dropped (counted and reported)

The model is cached next to the other derived references, so the GTF is parsed
once per species.
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
import pandas as pd

from .config import SITES_ROOT

MODEL_DIR = (_RB / "data/_refs/txmodels")


@dataclass(frozen=True)
class Transcript:
    chrom: str
    strand: str
    #: exon blocks in **transcript order** (0-based, half-open), so the first
    #: block starts at transcript position 0
    blocks: tuple[tuple[int, int], ...]
    length: int


def _parse_gtf(gtf: Path) -> dict[str, Transcript]:
    """exon features -> {transcript_id: Transcript} (blocks in transcript order).

    Both the bare ``transcript_id`` and, when the GTF carries a separate
    ``transcript_version`` (Ensembl/GENCODE style), the versioned
    ``transcript_id.version`` are registered, because the tool outputs use the
    versioned spelling while the GTF attributes are unversioned.
    """
    exons: dict[str, list[tuple[int, int]]] = {}
    meta: dict[str, tuple[str, str]] = {}
    with gtf.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            attrs: dict[str, str] = {}
            for kv in f[8].split(";"):
                kv = kv.strip()
                if " " in kv:
                    k, v = kv.split(" ", 1)
                    attrs[k] = v.strip().strip('"')
                elif "=" in kv:                      # GFF3-style attributes
                    k, v = kv.split("=", 1)
                    attrs[k.strip()] = v.strip()
            tid = attrs.get("transcript_id") or attrs.get("Parent", "").split(":")[-1]
            if not tid:
                continue
            keys = [tid]
            ver = attrs.get("transcript_version")
            if ver:
                keys.append(f"{tid}.{ver}")
            for key in keys:
                exons.setdefault(key, []).append((int(f[3]) - 1, int(f[4])))
                meta.setdefault(key, (f[0], f[6]))
    model: dict[str, Transcript] = {}
    for tid, blocks in exons.items():
        chrom, strand = meta[tid]
        ordered = sorted(blocks, key=lambda b: b[0], reverse=(strand == "-"))
        length = sum(e - s for s, e in ordered)
        model[tid] = Transcript(chrom, strand, tuple(ordered), length)
    return model


def load_model(gtf: Path, *, rebuild: bool = False) -> dict[str, Transcript]:
    """Cached transcript model for one GTF."""
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    cache = MODEL_DIR / (gtf.name + ".txmodel.pkl")
    if cache.exists() and not rebuild and cache.stat().st_mtime >= gtf.stat().st_mtime:
        with cache.open("rb") as fh:
            return pickle.load(fh)
    model = _parse_gtf(gtf)
    with cache.open("wb") as fh:
        pickle.dump(model, fh, protocol=4)
    return model


def map_positions(model: dict[str, Transcript], tx_ids: pd.Series,
                  positions: pd.Series) -> pd.DataFrame:
    """Map transcript positions to genomic coordinates.

    ``positions`` are 0-based transcript coordinates (the convention used by the
    legacy inputs).  Returns a frame with ``chromosome``, ``start`` (0-based),
    ``end`` (= start + 1), ``strand`` and ``mapped`` (False when the transcript
    is absent from the model or the position is outside it).
    """
    tx = tx_ids.astype(str).to_numpy()
    pos = pd.to_numeric(positions, errors="coerce").to_numpy()
    n = len(tx)
    chrom = np.empty(n, dtype=object)
    strand = np.empty(n, dtype=object)
    gstart = np.full(n, -1, dtype=np.int64)
    mapped = np.zeros(n, dtype=bool)
    for i in range(n):
        p = pos[i]
        t = model.get(tx[i])
        if t is None and "." in tx[i]:        # versioned id -> bare id fallback
            t = model.get(tx[i].split(".")[0])
        if t is None or not np.isfinite(p):
            continue
        p = int(p)
        if p < 0 or p >= t.length:
            continue
        offset = p
        for s, e in t.blocks:
            size = e - s
            if offset < size:
                chrom[i] = t.chrom
                strand[i] = t.strand
                if t.strand == "+":
                    gstart[i] = s + offset
                else:
                    gstart[i] = e - 1 - offset
                mapped[i] = True
                break
            offset -= size
    return pd.DataFrame({
        "chromosome": chrom, "start": gstart, "end": gstart + 1,
        "strand": strand, "mapped": mapped,
    })


def liftover_like_r2d(df: pd.DataFrame, model: dict[str, Transcript],
                      tx_col: str = "Chr", pos_col: str = "Start") -> pd.DataFrame:
    """Reproduce the legacy ``*_liftover.txt`` layout for a legacy input frame.

    The result keeps the original columns and prepends the genomic block, exactly
    like ``r2d liftover -H`` does, so downstream parsers are unchanged.
    """
    mapped = map_positions(model, df[tx_col], df[pos_col])
    keep = mapped["mapped"].to_numpy()
    head = pd.DataFrame({
        "chromosome": mapped.loc[keep, "chromosome"].to_numpy(),
        "start": mapped.loc[keep, "start"].to_numpy(),
        "end": mapped.loc[keep, "end"].to_numpy(),
        # r2d writes two empty columns (name / score) before the strand
        "name": "",
        "score": "",
        "strand": mapped.loc[keep, "strand"].to_numpy(),
    })
    return pd.concat([head, df[keep].reset_index(drop=True)], axis=1)
