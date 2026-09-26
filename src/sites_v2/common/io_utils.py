"""I/O helpers for the sites_v2 rebuild.

Design notes
------------
* Everything is UTF-8 tab-separated text; float formatting is left to pandas.
* FASTA access is cached per process (pyfaidx random access).  For whole-
  chromosome scans the sequence is fetched once per chromosome and sliced.
* Hashing policy: sha256 for files < 512 MB, ``size+mtime`` for anything larger
  (CephFS: hashing multi-GB files is slow and can stall on cold data).
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
import hashlib
import os
import re
import time
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np
import pandas as pd
from pyfaidx import Fasta

from .config import RESULT_RNA002, RESULT_RNA004

HASH_SIZE_LIMIT = 512 * 1024 * 1024  # 512 MB


# --------------------------------------------------------------------------- #
# small file metadata
# --------------------------------------------------------------------------- #
def file_fingerprint(path: Path) -> str:
    """sha256 for small files, ``size:mtime`` for large ones."""
    try:
        st = path.stat()
    except FileNotFoundError:
        return "missing"
    if st.st_size >= HASH_SIZE_LIMIT:
        return f"size={st.st_size};mtime={int(st.st_mtime)}"
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return "sha256=" + h.hexdigest()


def rel(path: Path) -> str:
    """Path relative to the benchmark project root (for compact manifests)."""
    for root in (RESULT_RNA002, RESULT_RNA004):
        try:
            return str(Path("result") / path.relative_to(root)) if root.name == "result" else str(
                Path("result_RNA004") / path.relative_to(root))
        except ValueError:
            continue
    try:
        return str(path.relative_to(Path(str(_RB))))
    except ValueError:
        return str(path)


# --------------------------------------------------------------------------- #
# TSV
# --------------------------------------------------------------------------- #
def read_table(path: Path | str, **kw) -> pd.DataFrame:
    """Read a TSV; empty files return an empty frame instead of raising."""
    path = Path(path)
    if path.stat().st_size == 0:
        return pd.DataFrame()
    kw.setdefault("sep", "\t")
    kw.setdefault("dtype", str)
    kw.setdefault("keep_default_na", False)
    return pd.read_csv(path, **kw)


def write_table(df: pd.DataFrame, path: Path, *, index: bool = False,
                float_format: str = "%.6g") -> int:
    """Write a TSV (creating parents) and return the row count."""
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=index, float_format=float_format)
    return len(df)


# --------------------------------------------------------------------------- #
# FASTA
# --------------------------------------------------------------------------- #
class FastaCache:
    """Per-process cache of pyfaidx Fasta objects and whole-chromosome strings."""

    def __init__(self) -> None:
        self._fasta: dict[Path, Fasta] = {}
        self._chrom: dict[tuple[Path, str], str] = {}

    def fasta(self, path: Path | str) -> Fasta:
        path = Path(path)
        if path not in self._fasta:
            self._fasta[path] = Fasta(str(path), as_raw=False,
                                      sequence_always_upper=True)
        return self._fasta[path]

    def keys(self, path: Path | str) -> list[str]:
        return list(self.fasta(path).keys())

    def chrom(self, path: Path | str, chrom: str) -> str | None:
        """Whole-chromosome sequence (careful: human chr1 ≈ 250 MB)."""
        key = (Path(path), chrom)
        if key not in self._chrom:
            fa = self.fasta(path)
            if chrom not in fa:
                return None
            self._chrom[key] = str(fa[chrom][:])
        return self._chrom[key]

    def fetch(self, path: Path | str, chrom: str, start0: int, end0: int) -> str | None:
        """1-bp-resolution slice ``[start0, end0)`` (0-based, half-open)."""
        s = self.chrom(path, chrom)
        if s is None:
            return None
        if start0 < 0 or end0 > len(s):
            return None
        return s[start0:end0]


FASTA = FastaCache()


#: Organelle labels are spelled differently by every source in this project:
#: the callsets say ``chrm`` / ``chrpt``, TAIR10 says ``Mt`` / ``Pt``, the mouse
#: and human primary assemblies say ``MT``.  Without these, mitochondrial calls
#: silently fail to resolve (annotate / universe / evaluation all go through this
#: one function, so they all lose the same rows).
_CHROM_ALIASES: dict[str, tuple[str, ...]] = {
    "m": ("MT", "Mt", "mito", "chrM", "M"),
    "mt": ("Mt", "mito", "M"),
    "pt": ("Pt", "cp", "chrP", "P"),
    "c": ("Pt", "CP", "chrC"),
    "chrm": ("MT", "Mt", "M"),
    "chrpt": ("Pt", "P"),
}


def _scaffold_form(chrom: str) -> str | None:
    """``chrun_jh584304v1`` -> ``JH584304.1`` (GenBank unplaced-scaffold style)."""
    core = chrom[3:] if chrom.lower().startswith("chr") else chrom
    m = re.match(r"^un_([a-z]+\d+)v(\d+)$", core.lower())
    return f"{m.group(1).upper()}.{m.group(2)}" if m else None


def normalize_chrom_for_fasta(fa_keys: Iterable[str], chrom: str) -> str | None:
    """Map an arbitrary chromosome label onto one present in a FASTA index.

    Handles the ``1`` vs ``chr1`` vs ``Chromosome`` families used by the
    project's references and bed files, the organelle spellings in
    :data:`_CHROM_ALIASES`, and the ``un_<acc>v<n>`` / ``<ACC>.<n>`` scaffold
    styles.  Every candidate is checked against the actual index, so an alias can
    never invent a chromosome the assembly does not have.
    """
    keys = list(fa_keys)
    keyset = set(keys)
    low = {k.lower(): k for k in keys}
    stripped = chrom[3:] if chrom.lower().startswith("chr") else chrom
    cands = [chrom, f"chr{chrom}", stripped, chrom.upper(), chrom.lower(),
             stripped.upper(), stripped.lower()]
    for c in (*cands, *_CHROM_ALIASES.get(str(chrom).lower(), ()),
              *_CHROM_ALIASES.get(str(stripped).lower(), ()),
              (_scaffold_form(chrom) or "")):
        if not c:
            continue
        if c in keyset:
            return c
        if c.lower() in low:
            return low[c.lower()]
    # Ensembl bacteria names the single replicon plainly: "chromosome"
    if len(keys) == 1 and str(chrom).lower() in {"chromosome", "chr", "1", "genome"}:
        return keys[0]
    return None


# --------------------------------------------------------------------------- #
# misc
# --------------------------------------------------------------------------- #
def ensure_dirs(*paths: Path) -> None:
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)


def chunk_ranges(length: int, chunk: int = 2_000_000) -> Iterator[tuple[int, int]]:
    for start in range(0, length, chunk):
        yield start, min(start + chunk, length)


def to_int(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").astype("Int64")


def now_stamp() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")
