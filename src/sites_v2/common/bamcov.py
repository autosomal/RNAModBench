"""Per-position coverage from a genomic BAM.

Primary path: ``samtools depth -a -b positions.bed`` (single C pass over the
BAM; ~4 s for 40k positions on a 1M-read sample).  ``pysam.count_coverage`` is
kept as a fallback but is 30-100x slower on the same data.

Both paths count *every* aligned read (mapping/base quality 0) and exclude the
default flag set (UNMAP/SECONDARY/QCFAIL/DUP), so coverage is comparable across
samples and tools.
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
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from .config import MANIFEST_DIR
from .match import fix_chromosome

#: candidate samtools binaries (benchmark-revision has no samtools of its own)
SAMTOOLS_CANDIDATES = [
    "samtools",
    str(_XB / "miniconda3/envs/asaif/bin/samtools"),
    str(_XB / "miniconda3/envs/nanom6A/bin/samtools"),
    str(_XB / "miniconda3/envs/chipseq/bin/samtools"),
]


def find_samtools() -> str | None:
    for cand in SAMTOOLS_CANDIDATES:
        p = shutil.which(cand) if "/" not in cand else (cand if Path(cand).is_file() else None)
        if p:
            return p
    return None


def _bam_reference_map(bam_path: Path) -> dict[str, str]:
    import pysam

    with pysam.AlignmentFile(str(bam_path), "rb") as af:
        return {fix_chromosome(r): r for r in af.references}


def coverage_at_positions(bam_path: Path | str, chrom_positions: dict[str, np.ndarray],
                          chunk: int = 2_000_000) -> dict[str, np.ndarray]:
    """Coverage for every (chrom, pos) in ``chrom_positions``.

    Parameters
    ----------
    bam_path
        Coordinate-sorted, indexed BAM.
    chrom_positions
        {chromosome: int64 array of 0-based positions}; labels are normalised
        internally so ``1``/``chr1``/``Chromosome`` all resolve.

    Returns
    -------
    {chromosome: float array aligned with the input position array} (NaN where
    the chromosome is absent from the BAM).
    """
    bam_path = Path(bam_path)
    out: dict[str, np.ndarray] = {c: np.full(np.asarray(p).shape, np.nan, dtype=float)
                                  for c, p in chrom_positions.items()}
    ref_map = _bam_reference_map(bam_path)
    samtools = find_samtools()
    if samtools is not None:
        _coverage_samtools(bam_path, samtools, chrom_positions, ref_map, out)
    else:
        _coverage_pysam(bam_path, chrom_positions, ref_map, out, chunk)
    return out


def _coverage_samtools(bam_path: Path, samtools: str,
                       chrom_positions: dict[str, np.ndarray],
                       ref_map: dict[str, str],
                       out: dict[str, np.ndarray]) -> None:
    rows = []
    for chrom, pos in chrom_positions.items():
        ref = ref_map.get(fix_chromosome(chrom))
        if ref is None or len(pos) == 0:
            continue
        for p in np.asarray(pos, dtype=np.int64):
            rows.append((ref, int(p), int(p) + 1))
    if not rows:
        return
    rows.sort()
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", suffix=".bed", dir=str(MANIFEST_DIR),
                                     delete=False) as fh:
        bed = Path(fh.name)
        for ref, s, e in rows:
            fh.write(f"{ref}\t{s}\t{e}\n")
    try:
        cmd = [samtools, "depth", "-a", "-q", "0", "-Q", "0", "-d", "0", "-b", str(bed),
               str(bam_path)]
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        depth: dict[tuple[str, int], int] = {}
        for line in res.stdout.splitlines():
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 3:
                continue
            depth[(fix_chromosome(f[0]), int(f[1]) - 1)] = int(f[2])  # depth is 1-based
        for chrom, pos in chrom_positions.items():
            arr = out[chrom]
            if arr.size == 0:
                continue
            key_chrom = fix_chromosome(chrom)
            arr[:] = [depth.get((key_chrom, int(p)), 0) for p in np.asarray(pos, dtype=np.int64)]
    finally:
        bed.unlink(missing_ok=True)


def _coverage_pysam(bam_path: Path, chrom_positions: dict[str, np.ndarray],
                    ref_map: dict[str, str], out: dict[str, np.ndarray],
                    chunk: int) -> None:
    import pysam

    with pysam.AlignmentFile(str(bam_path), "rb") as af:
        lengths = {fix_chromosome(r): af.get_reference_length(r) for r in af.references}
        for chrom, pos in chrom_positions.items():
            ref = ref_map.get(fix_chromosome(chrom))
            values = out[chrom]
            if ref is None or values.size == 0:
                continue
            L = lengths[fix_chromosome(chrom)]
            pos = np.asarray(pos, dtype=np.int64)
            order = np.argsort(pos, kind="mergesort")
            pos_sorted = pos[order]
            for start in range(0, L, chunk):
                end = min(start + chunk, L)
                lo = np.searchsorted(pos_sorted, start, side="left")
                hi = np.searchsorted(pos_sorted, end, side="left")
                if hi <= lo:
                    continue
                cov = af.count_coverage(ref, start, end, quality_threshold=0)
                total = np.asarray(cov[0], dtype=np.int32)
                for arr in cov[1:]:
                    total = total + np.asarray(arr, dtype=np.int32)
                idx = pos_sorted[lo:hi] - start
                values[order[lo:hi]] = total[idx]
