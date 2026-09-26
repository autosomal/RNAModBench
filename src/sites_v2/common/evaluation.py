"""Shared evaluation helpers: universe loading and window-wise confusion counts.

All functions work on 0-based single-nucleotide positions grouped per
chromosome as sorted int64 arrays.  Chromosome labels are normalised with
``fix_chromosome`` on both sides so ``1``/``chr1``/``Chromosome`` all match.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import MOD_REF_BASE
from .io_utils import read_table
from .match import fix_chromosome
from .metrics import summarize, window_hit_mask

#: genome-orientation bases that can carry each modification (strand aware):
#: m6A sits on A ('+') or T ('-'); m5C on C ('+') or G ('-');
#: Psi/m1Psi on T ('+') or A ('-'); inosine on A/T; Nm any base.
MOD_GENOME_BASES: dict[str, set[str]] = {
    "m6A": {"A", "T"}, "inosine": {"A", "T"},
    "m5C": {"C", "G"},
    "Psi": {"T", "A"}, "m1Psi": {"T", "A"},
    "Nm": {"A", "C", "G", "T"},
}


def load_universe(path: Path, mod_type: str, min_cov: int = 10
                  ) -> dict[str, np.ndarray]:
    """Universe positions for one modification type / coverage threshold.

    ``path`` may be plain or gzipped TSV with columns
    ``chrom pos base coverage drach in_glori``.
    """
    usecols = ["chrom", "pos", "base", "coverage"]
    df = pd.read_csv(path, sep="\t", usecols=usecols, dtype={"chrom": str, "base": str})
    want = MOD_GENOME_BASES.get(mod_type, {"A", "C", "G", "T"})
    df = df[(df["coverage"] >= min_cov) & (df["base"].isin(want))]
    df["chrom"] = [fix_chromosome(c) for c in df["chrom"]]
    out: dict[str, np.ndarray] = {}
    for chrom, sub in df.groupby("chrom", sort=False):
        out[chrom] = np.unique(sub["pos"].to_numpy(dtype=np.int64))
    return out


def positions_in_universe(df: pd.DataFrame, universe: dict[str, np.ndarray]
                          ) -> tuple[dict[str, np.ndarray], int]:
    """Restrict a callset frame to universe positions.

    Returns ({chrom: sorted call positions}, n_calls_out_of_universe).
    """
    out: dict[str, np.ndarray] = {}
    n_out = 0
    for chrom, sub in df.groupby("chrom", sort=False):
        u = universe.get(chrom)
        pos = sub["pos_raw"].to_numpy(dtype=np.int64)
        if u is None or u.size == 0:
            n_out += pos.size
            continue
        idx = np.searchsorted(u, pos)
        idx_clip = np.clip(idx, 0, u.size - 1)
        inside = u[idx_clip] == pos
        n_out += int((~inside).sum())
        if inside.any():
            out[chrom] = np.unique(pos[inside])
    return out, n_out


def reference_in_universe(reference: dict[str, np.ndarray],
                          universe: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Reference sites that are inside the universe (testable reference sites)."""
    out: dict[str, np.ndarray] = {}
    for chrom, ref in reference.items():
        u = universe.get(chrom)
        if u is None or u.size == 0:
            continue
        idx = np.searchsorted(u, ref)
        idx_clip = np.clip(idx, 0, u.size - 1)
        keep = u[idx_clip] == ref
        if keep.any():
            out[chrom] = ref[keep]
    return out


def confusion_windows(universe: dict[str, np.ndarray],
                      reference_u: dict[str, np.ndarray],
                      calls_u: dict[str, np.ndarray],
                      windows: list[int]) -> dict[int, dict]:
    """TP/FP/FN/TN per window (all inputs already restricted to the universe)."""
    n_u = int(sum(v.size for v in universe.values()))
    counts = {w: [0, 0, 0] for w in windows}  # tp, fp, fn
    for chrom, u in universe.items():
        calls = calls_u.get(chrom, np.zeros(0, dtype=np.int64))
        refs = reference_u.get(chrom, np.zeros(0, dtype=np.int64))
        for w in windows:
            hit_call = window_hit_mask(calls, refs, w) if calls.size else np.zeros(0, bool)
            tp = int(hit_call.sum())
            fp = int(calls.size - tp)
            hit_ref = window_hit_mask(refs, calls, w) if refs.size else np.zeros(0, bool)
            fn = int(refs.size - int(hit_ref.sum()))
            counts[w][0] += tp
            counts[w][1] += fp
            counts[w][2] += fn
    out: dict[int, dict] = {}
    for w, (tp, fp, fn) in counts.items():
        tn = n_u - tp - fp - fn
        rec = {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "universe": n_u,
               "n_calls_in_universe": tp + fp, "n_reference_in_universe": tp + fn}
        rec.update(summarize(tp, fp, fn, tn))
        out[w] = rec
    return out


def strand_base_consistency(df: pd.DataFrame, mod_type: str, genome_fa: Path
                            ) -> tuple[int, int]:
    """(#calls whose +/- strand matches the reference base, #calls checked).

    A call is inconsistent when it is reported on '+' but the genome base is the
    minus-strand version of the target (or vice versa).  Unknown strand ('*') is
    always counted as consistent (it is checked, but cannot be contradicted).
    """
    from .io_utils import FASTA, normalize_chrom_for_fasta

    want = MOD_REF_BASE.get(mod_type)
    if want is None:
        return 0, 0
    plus = want
    minus = {"A": "T", "T": "A", "C": "G", "G": "C"}[want]
    fa_keys = FASTA.keys(genome_fa)
    ok = bad = 0
    for chrom, sub in df.groupby("chrom", sort=False):
        fa_name = normalize_chrom_for_fasta(fa_keys, chrom)
        if fa_name is None:
            continue
        seq = FASTA.chrom(genome_fa, fa_name)
        arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
        pos = sub["pos_raw"].to_numpy(dtype=np.int64)
        keep = (pos >= 0) & (pos < arr.size)
        if not keep.any():
            continue
        bases = np.array([chr(b) for b in arr[pos[keep]]])
        strands = sub["strand"].astype(str).to_numpy()[keep]
        expect_plus = strands == "+"
        expect_minus = strands == "-"
        good = (expect_plus & (bases == plus)) | (expect_minus & (bases == minus)) \
            | (~expect_plus & ~expect_minus)
        ok += int(good.sum())
        bad += int((~good).sum())
    return ok, ok + bad


def load_reference(path: Path) -> dict[str, np.ndarray]:
    """Load a reference bed (GLORI / Curlcake truth) as {chrom: sorted positions}."""
    out: dict[str, list[int]] = {}
    with open(path) as fh:
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
