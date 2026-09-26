"""Site annotation: strand-aware 5mer, DRACH context, centre correction, GLORI
distance.

Conventions (biology / coordinates)
-----------------------------------
* ``pos_raw`` is the 0-based start of a 1-bp interval; BED length is
  ``end - start = 1`` (0-based half-open, never ``+1``).
* The *plus-strand* 5mer around a site ``p`` is ``genome[p-2 : p+3]``.  For a
  site reported on the minus strand the transcript-orientation 5mer is the
  reverse complement of that slice.  All downstream motif logic uses the
  transcript-orientation 5mer.
* ``pos_center`` is the nearest position (within +/- ``CENTER_SEARCH_MAX``)
  whose *genome* base equals the expected base in transcript orientation
  (m6A -> A on '+', T on '-'; m5C -> C on '+', G on '-'; ...).  Offsets are
  reported, never applied to ``pos_raw``.
* ``dist_drach_a`` is the signed offset to the nearest DRACH-modified A
  (motif ``[AGT][AG]AC[ACT]`` with the A at index 2), again strand-aware.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import CENTER_SEARCH_MAX, DRACH_REGEX, MOD_REF_BASE
from .io_utils import FASTA, normalize_chrom_for_fasta

# byte codes
_A, _C, _G, _T = 65, 67, 71, 84
_COMPLEMENT = {_A: _T, _T: _A, _C: _G, _G: _C}
DRACH_ALLOWED = ((_A, _G, _T), (_A, _G), (_A,), (_C,), (_A, _C, _T))

ANNOTATION_COLUMNS = [
    "ref_base", "five_mer_raw", "five_mer_center", "center_base_expected",
    "is_drach_raw", "is_drach_center", "pos_center", "dist_center",
    "dist_drach_a", "coverage", "cov_source", "dist_glori", "glori_exact",
    "offset_flag", "center_status",
]

#: ``center_status`` values -- the strict, chain-aware verdict per row
#: (see ``common.center``).  A row whose strand is unknown cannot be judged.
CENTER_STATUS_LABEL = {0: "unknown_strand", 1: "ok", 2: "off_base",
                       3: "no_expectation"}


def _seq_array(seq: str) -> np.ndarray:
    return np.frombuffer(seq.encode("ascii"), dtype=np.uint8)


def _revcomp_bytes(mat: np.ndarray) -> np.ndarray:
    """Reverse complement a (n, 5) uint8 matrix."""
    out = mat[:, ::-1].copy()
    mapped = np.full(out.shape, _A, dtype=np.uint8)
    for src, dst in _COMPLEMENT.items():
        mapped[out == src] = dst
    return mapped


def _gather(arr: np.ndarray, pos: np.ndarray, offsets: list[int]) -> tuple[np.ndarray, np.ndarray]:
    """arr[pos[:, None] + offsets] with clipping; also returns the validity mask."""
    L = arr.size
    idx = pos[:, None] + np.asarray(offsets, dtype=np.int64)[None, :]
    valid = (idx >= 0) & (idx < L)
    return arr[np.clip(idx, 0, L - 1)], valid


def five_mer_matrix(arr: np.ndarray, pos: np.ndarray, strand: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(n,5) transcript-orientation 5mer matrix + validity mask."""
    mat, valid = _gather(arr, pos, [-2, -1, 0, 1, 2])
    row_valid = valid.all(axis=1)
    minus = strand == "-"
    if minus.any():
        mat = mat.copy()
        mat[minus] = _revcomp_bytes(mat[minus])
    return mat, row_valid


def _as_str(mat: np.ndarray, row_valid: np.ndarray) -> np.ndarray:
    out = np.full(mat.shape[0], "NNNNN", dtype=object)
    if row_valid.any():
        sub = mat[row_valid]
        out[row_valid] = [bytes(r).decode("ascii", "replace") for r in sub]
    return out


def _drach_mask(mat: np.ndarray, row_valid: np.ndarray) -> np.ndarray:
    ok = row_valid.copy()
    for j, allowed in enumerate(DRACH_ALLOWED):
        ok &= np.isin(mat[:, j], allowed)
    return ok


def _nearest_offset(expected: np.ndarray, kmax: int) -> np.ndarray:
    """First True offset in 0, +1, -1, +2, -2, ... order; NaN if none."""
    n = expected.shape[1]
    order = [0] + [o for k in range(1, kmax + 1) for o in (k, -k)]
    out = np.full(expected.shape[0], np.nan)
    filled = np.zeros(expected.shape[0], dtype=bool)
    cols = list(range(n))
    for off in order:
        col = cols[off + kmax]
        hit = expected[:, col] & ~filled
        out[hit] = off
        filled |= hit
        if filled.all():
            break
    return out


def annotate_frame(df: pd.DataFrame, *, genome: Path, mod_type: str,
                   coverage: np.ndarray | None = None,
                   cov_source: str = "", glori: dict[str, np.ndarray] | None = None,
                   search_max: int = CENTER_SEARCH_MAX) -> pd.DataFrame:
    """Add annotation columns to a callset frame (returns a new frame)."""
    out = df.drop(columns=[c for c in ANNOTATION_COLUMNS if c in df.columns]).copy()
    out["ref_base"] = ""
    out["five_mer_raw"] = "NNNNN"
    out["five_mer_center"] = "NNNNN"
    out["center_base_expected"] = False
    out["is_drach_raw"] = False
    out["is_drach_center"] = False
    out["pos_center"] = pd.NA
    out["dist_center"] = np.nan
    out["dist_drach_a"] = np.nan
    out["coverage"] = np.nan
    out["cov_source"] = cov_source
    out["dist_glori"] = np.nan
    out["glori_exact"] = False
    out["offset_flag"] = ""
    out["center_status"] = CENTER_STATUS_LABEL[0]

    pos = out["pos_raw"].to_numpy(dtype=np.int64)
    strand = out["strand"].astype(str).to_numpy()
    chrom = out["chrom"].astype(str).to_numpy()
    #: a strand that is neither '+' nor '-' (Nanom6A writes '*') makes every
    #: transcript-orientation statement -- 5mer, DRACH, centre -- undefined.  It
    #: used to be silently treated as '+', which moved ~half of those sites
    #: 2-5 bp through ``pos_center``; now those rows are simply not judged
    #: (``32_impute_strand`` fills the strand from the reads / exons first).
    known = np.isin(strand, ["+", "-"])

    if coverage is not None and len(coverage) == len(out):
        out["coverage"] = coverage

    fa = FASTA.fasta(genome)
    fa_keys = list(fa.keys())
    target = MOD_REF_BASE.get(mod_type)
    genome_targets: dict[str, int] = {}
    if target is not None:
        t = ord(target)
        genome_targets = {"+": t, "-": {_A: _T, _T: _A, _C: _G, _G: _C}[t], "*": t}

    for c in np.unique(chrom):
        keys = np.flatnonzero(chrom == c)
        if keys.size == 0:
            continue
        fa_name = normalize_chrom_for_fasta(fa_keys, c)
        if fa_name is None:
            continue
        seq = FASTA.chrom(genome, fa_name)
        arr = _seq_array(seq)
        p = pos[keys]
        s = strand[keys]

        ref = arr[np.clip(p, 0, arr.size - 1)]
        out.loc[out.index[keys], "ref_base"] = [chr(b) for b in ref]

        # ``five_mer_raw`` = GENOME-orientation 5mer at ``pos_raw``
        # (``genome[p-2:p+3]``), exactly as the legacy ``extract_5mer.py`` sliced
        # it.  It is deliberately NOT strand-aware so it reconciles with the old
        # ``output/`` pipeline.  ALL motif logic (DRACH / centre) uses the
        # strand-aware ``five_mer_center`` below instead.
        mat_plus, valid_plus = _gather(arr, p, [-2, -1, 0, 1, 2])
        out.loc[out.index[keys], "five_mer_raw"] = _as_str(mat_plus, valid_plus.all(axis=1))

        # transcript-orientation 5mer at pos_raw (reverse-complemented on '-')
        mat, row_valid = five_mer_matrix(arr, p, s)
        out.loc[out.index[keys], "is_drach_raw"] = _drach_mask(mat, row_valid) & known[keys]

        # --- centre correction -------------------------------------------------
        if target is not None:
            offs = list(range(-search_max, search_max + 1))
            gathered, valid = _gather(arr, p, offs)
            exp = np.zeros(gathered.shape, dtype=bool)
            k = known[keys]
            for j, off in enumerate(offs):
                want = np.where(s == "-", genome_targets["-"], genome_targets["+"])
                exp[:, j] = (gathered[:, j] == want) & valid[:, j] & k
            dist_center = _nearest_offset(exp, search_max)
            pos_center = p + np.nan_to_num(dist_center, nan=0).astype(np.int64)
            base_ok = ~np.isnan(dist_center)
            out.loc[out.index[keys], "dist_center"] = dist_center
            # store as nullable Int64 so the TSV shows "50240" not "50240.0"
            ser = pd.Series(np.where(base_ok, pos_center, np.nan),
                            index=out.index[keys]).astype("Int64")
            out.loc[out.index[keys], "pos_center"] = ser
            out.loc[out.index[keys], "center_base_expected"] = base_ok

            # 5mer / DRACH at the *centred* position (strand-aware)
            mat_c, valid_c = five_mer_matrix(arr, pos_center, s)
            out.loc[out.index[keys], "five_mer_center"] = _as_str(mat_c, valid_c & base_ok)
            out.loc[out.index[keys], "is_drach_center"] = _drach_mask(mat_c, valid_c & base_ok)

        # --- nearest DRACH A ---------------------------------------------------
        kmax = search_max + 2  # motif centre may sit up to 2 bp off centre
        offs = [o for k in range(0, kmax + 1) for o in ((0,) if k == 0 else (k, -k))]
        exp = np.zeros((p.size, len(offs)), dtype=bool)
        for j, off in enumerate(offs):
            mat_o, valid_o = five_mer_matrix(arr, p + off, s)
            exp[:, j] = _drach_mask(mat_o, valid_o)
        dist_drach = np.full(p.size, np.nan)
        filled = np.zeros(p.size, dtype=bool)
        for j, off in enumerate(offs):
            hit = exp[:, j] & ~filled
            dist_drach[hit] = off
            filled |= hit
        dist_drach[~known[keys]] = np.nan
        out.loc[out.index[keys], "dist_drach_a"] = dist_drach

    # --- GLORI distance -------------------------------------------------------
    if glori:
        from .match import exact_hit_mask, nearest_distance
        dist = np.full(len(out), np.nan)
        exact = np.zeros(len(out), dtype=bool)
        for c in np.unique(chrom):
            gpos = glori.get(c)
            if gpos is None or gpos.size == 0:
                continue
            keys = np.flatnonzero(chrom == c)
            dist[keys] = nearest_distance(pos[keys], gpos)
            exact[keys] = exact_hit_mask(pos[keys], gpos)
        out["dist_glori"] = dist
        out["glori_exact"] = exact

    # --- strict centre-base verdict ------------------------------------------
    # ``ok`` / ``off_base`` / ``unknown_strand`` / ``no_expectation`` -- the
    # only statement about a call's centre base that can be checked against the
    # reference.  (``center_base_expected`` searches +/- 5 bp and is therefore a
    # diagnostic, never a verdict.)
    from .center import strict_status as _strict_status

    code = _strict_status(out["strand"].astype(str).to_numpy(),
                          out["ref_base"].astype(str).to_numpy(), mod_type)
    out["center_status"] = [CENTER_STATUS_LABEL[int(c)] for c in code]
    return out


def load_glori(path: Path) -> dict[str, np.ndarray]:
    """GLORI bed (0-based, end-start==1) -> {chrom: sorted positions}."""
    from .match import fix_chromosome

    out: dict[str, list[int]] = {}
    with open(path, "r") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                start = int(f[1])
            except ValueError:
                continue
            out.setdefault(fix_chromosome(f[0]), []).append(start)
    return {c: np.array(sorted(set(v)), dtype=np.int64) for c, v in out.items()}
