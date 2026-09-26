"""Vectorised site matching helpers.

The window/exact/recall functions are copied verbatim from
``code/revision/common/match.py`` (verified there to reproduce the published
hit rates) so that the harmonisation evaluation and the NA-1/NA-9 revision analyses
share one implementation.  Two additions are made here:

* :func:`signed_nearest_distance` - signed distance ``tool - reference`` used by
  the per-tool offset audit (positive = tool site is downstream).
* :func:`nearest_reference_position` - the reference coordinate that minimises
  the distance, so callers can report *which* site was matched.

All coordinates are 0-based positions of a single reference nucleotide.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# chromosome normalisation  (verbatim from code/revision/common/match.py)
# --------------------------------------------------------------------------- #
_SPECIAL = {"m": "chrm", "mt": "chrm", "mito": "chrm", "x": "chrx", "y": "chry",
            "23": "chrx", "24": "chry", "25": "chrm"}


def fix_chromosome(chrom) -> str:
    """Normalise any chromosome label to lowercase ``chr*`` form."""
    if chrom is None or (isinstance(chrom, float) and np.isnan(chrom)):
        return "unknown"
    s = str(chrom).strip().lower()
    if s.startswith("chr"):
        base = s[3:]
    else:
        base = s
    if base in _SPECIAL:
        return _SPECIAL[base]
    return f"chr{base}"


# --------------------------------------------------------------------------- #
# grouping helpers (verbatim)
# --------------------------------------------------------------------------- #
def group_sorted_positions(df: pd.DataFrame, chrom_col: str = "chrom",
                           pos_col: str = "pos") -> dict[str, np.ndarray]:
    """Return {chromosome: sorted unique 1-D int64 array of positions}."""
    out: dict[str, np.ndarray] = {}
    for chrom, sub in df.groupby(chrom_col, sort=False):
        out[chrom] = np.unique(sub[pos_col].to_numpy(dtype=np.int64))
    return out


def group_sorted_sites(df: pd.DataFrame, chrom_col: str = "chrom",
                       pos_col: str = "pos") -> dict[str, dict[str, np.ndarray]]:
    """Return {chromosome: {'pos': sorted positions, 'idx': original row index}}."""
    out: dict[str, dict[str, np.ndarray]] = {}
    for chrom, sub in df.groupby(chrom_col, sort=False):
        order = np.argsort(sub[pos_col].to_numpy(dtype=np.int64), kind="mergesort")
        pos = sub[pos_col].to_numpy(dtype=np.int64)[order]
        idx = sub.index.to_numpy()[order]
        out[chrom] = {"pos": pos, "idx": idx}
    return out


# --------------------------------------------------------------------------- #
# core matching (verbatim)
# --------------------------------------------------------------------------- #
def window_hit_mask(tool_pos: np.ndarray, glori_pos: np.ndarray, w: int) -> np.ndarray:
    """Boolean mask: does each tool site have a reference site within +/- w nt?"""
    if glori_pos.size == 0 or tool_pos.size == 0:
        return np.zeros(tool_pos.shape, dtype=bool)
    glori_pos = np.asarray(glori_pos, dtype=np.int64)
    left = np.searchsorted(glori_pos, tool_pos - w, side="left")
    right = np.searchsorted(glori_pos, tool_pos + w, side="right")
    return right > left


def exact_hit_mask(tool_pos: np.ndarray, glori_sorted: np.ndarray) -> np.ndarray:
    """Boolean mask: is the tool site an exact single-nucleotide match?"""
    if glori_sorted.size == 0 or tool_pos.size == 0:
        return np.zeros(tool_pos.shape, dtype=bool)
    i = np.searchsorted(glori_sorted, tool_pos, side="left")
    i_clipped = np.clip(i, 0, glori_sorted.size - 1)
    return glori_sorted[i_clipped] == tool_pos


def nearest_distance(tool_pos: np.ndarray, glori_sorted: np.ndarray) -> np.ndarray:
    """Distance from each tool site to the closest reference site (nt)."""
    if tool_pos.size == 0:
        return np.zeros(0, dtype=float)
    if glori_sorted.size == 0:
        return np.full(tool_pos.shape, np.nan, dtype=float)
    n = glori_sorted.size
    if n == 1:
        return np.abs(tool_pos - glori_sorted[0]).astype(float)
    i = np.clip(np.searchsorted(glori_sorted, tool_pos, side="left"), 1, n - 1)
    left_val = glori_sorted[i - 1]
    right_val = glori_sorted[i]
    return np.minimum(np.abs(tool_pos - left_val), np.abs(tool_pos - right_val)).astype(float)


def glori_recall_mask(glori_pos: np.ndarray, tool_sorted: np.ndarray, w: int) -> np.ndarray:
    """Boolean mask over reference sites: is a tool site within +/- w nt?"""
    if tool_sorted.size == 0 or glori_pos.size == 0:
        return np.zeros(glori_pos.shape, dtype=bool)
    left = np.searchsorted(tool_sorted, glori_pos - w, side="left")
    right = np.searchsorted(tool_sorted, glori_pos + w, side="right")
    return right > left


# --------------------------------------------------------------------------- #
# additions for the offset audit
# --------------------------------------------------------------------------- #
def nearest_reference_position(tool_pos: np.ndarray,
                               reference_sorted: np.ndarray) -> np.ndarray:
    """Reference coordinate minimising |tool - reference| (NaN if no refs)."""
    if tool_pos.size == 0:
        return np.zeros(0, dtype=float)
    if reference_sorted.size == 0:
        return np.full(tool_pos.shape, np.nan, dtype=float)
    n = reference_sorted.size
    if n == 1:
        return np.full(tool_pos.shape, float(reference_sorted[0]))
    i = np.clip(np.searchsorted(reference_sorted, tool_pos, side="left"), 1, n - 1)
    left_val = reference_sorted[i - 1]
    right_val = reference_sorted[i]
    take_left = (tool_pos - left_val) <= (right_val - tool_pos)
    out = np.where(take_left, left_val, right_val)
    return out.astype(float)


def signed_nearest_distance(tool_pos: np.ndarray,
                            reference_sorted: np.ndarray) -> np.ndarray:
    """Signed distance ``tool_pos - nearest_reference_pos`` (float, NaN padded)."""
    ref = nearest_reference_position(tool_pos, reference_sorted)
    return tool_pos.astype(float) - ref


def match_sites(tool_df: pd.DataFrame, glori: dict[str, np.ndarray], window: int,
                chrom_col: str = "chrom", pos_col: str = "pos"
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (window_hit, exact_hit, nearest_dist) aligned with ``tool_df``."""
    n = len(tool_df)
    win_hit = np.zeros(n, dtype=bool)
    exact_hit = np.zeros(n, dtype=bool)
    dist = np.full(n, np.nan, dtype=float)
    grouped = group_sorted_sites(tool_df, chrom_col=chrom_col, pos_col=pos_col)
    for chrom, grp in grouped.items():
        gpos = glori.get(chrom)
        if gpos is None or gpos.size == 0:
            continue
        pos, idx = grp["pos"], grp["idx"]
        win_hit[idx] = window_hit_mask(pos, gpos, window)
        exact_hit[idx] = exact_hit_mask(pos, gpos)
        dist[idx] = nearest_distance(pos, gpos)
    return win_hit, exact_hit, dist


def iter_chrom_pairs(tool_df: pd.DataFrame, glori: dict[str, np.ndarray],
                     chrom_col: str = "chrom", pos_col: str = "pos"
                     ) -> Iterable[tuple]:
    """Yield ``(chrom, tool_positions, glori_positions)`` for shared chromosomes."""
    grouped = group_sorted_positions(tool_df, chrom_col=chrom_col, pos_col=pos_col)
    for chrom, pos in grouped.items():
        g = glori.get(chrom)
        if g is not None and g.size:
            yield chrom, pos, g
