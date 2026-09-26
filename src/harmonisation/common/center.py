"""Chain-aware centre-base judgement (the "is this call chemically possible?" rule).

Why this module exists
----------------------
``common/annotate.py`` annotates a callset but deliberately *reports* offsets
instead of judging them: ``center_base_expected`` is True whenever the expected
base is found anywhere within +/- ``CENTER_SEARCH_MAX`` (5 bp), which is almost
always, so it can never serve as a pass/fail criterion.  The only strict
statement is the transcript-orientation centre base at ``pos_raw`` itself:

    m6A     : transcript 'A'  -> genome 'A' on '+', genome 'T' on '-'
    m5C     : transcript 'C'  -> genome 'C' / 'G'
    Psi     : transcript 'U'  -> genome 'T' / 'A'   ({U -> T} DNA alphabet)
    m1Psi   : transcript 'U'  -> genome 'T' / 'A'
    inosine : transcript 'A'  -> genome 'A' / 'T'   (edited A)
    Nm      : any base (no single-base expectation)

A call whose centre base does not match is *chemically impossible* for that
modification: either the coordinate is wrong, the reference frame is wrong, or
the row is noise.  ``29_anchor_audit`` used the loose, strand-blind version
(``ref_base in {expected, complement}``) which lets a half-wrong callset pass;
this module is the strict version shared by the audit (29), the strand
imputation (32) and the hard filter (33) so the three can never disagree.

All coordinates are 0-based single nucleotides; the module never converts
between 0/1-based systems.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd

from .config import MOD_REF_BASE

#: complement of the expected base (genome orientation on the minus strand).
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

#: DRACH in transcript orientation: [AGT][AG]AC[ACT] with the A at index 2.
DRACH_ALLOWED = (frozenset("AGT"), frozenset("AG"), frozenset("A"),
                 frozenset("C"), frozenset("ACT"))

_RC = str.maketrans("ACGTN", "TGCAN")

#: status codes returned by :func:`strict_status`.
UNKNOWN_STRAND = 0
OK = 1
OFF_BASE = 2
NO_EXPECTATION = 3


def expected_base_for(mod_type: str, strand: str) -> str | None:
    """Genome base the modification must sit on, given the transcript strand.

    ``None`` means "no single-base expectation" -- either the modification can
    sit on any base (Nm) or the strand is unknown, in which case the caller has
    to fall back to the loose check (or treat the row as unverifiable).
    """
    base = MOD_REF_BASE.get(str(mod_type))
    if base is None:
        return None
    if strand == "+":
        return base
    if strand == "-":
        return COMPLEMENT[base]
    return None


def expected_base_pair(mod_type: str) -> tuple[str, ...]:
    """``(expected, complement)`` -- the strand-blind (loose) expectation."""
    base = MOD_REF_BASE.get(str(mod_type))
    if base is None:
        return ()
    return (base, COMPLEMENT[base])


def strict_status(strand: np.ndarray, ref_base: np.ndarray,
                  mod_type: str) -> np.ndarray:
    """uint8 status per row: 0 unknown strand, 1 ok, 2 off-base, 3 no expectation.

    ``strand``/``ref_base`` are string arrays (anything that is not exactly
    ``'+'``/``'-'`` counts as unknown -- notably Nanom6A's ``'*'``).
    """
    strand = np.asarray(strand).astype(str)
    ref = np.asarray(ref_base).astype(str)
    base = MOD_REF_BASE.get(mod_type)
    if base is None:
        return np.full(len(ref), NO_EXPECTATION, dtype=np.uint8)
    comp = COMPLEMENT[base]
    plus = strand == "+"
    minus = strand == "-"
    ok = (plus & (ref == base)) | (minus & (ref == comp))
    known = plus | minus
    out = np.full(len(ref), UNKNOWN_STRAND, dtype=np.uint8)
    out[known & ok] = OK
    out[known & ~ok] = OFF_BASE
    return out


def drach_mask(five_mer: Iterable[str], revcomp: bool = False) -> np.ndarray:
    """DRACH hit mask for genome-orientation 5mers.

    ``revcomp=True`` first reverse-complements each 5mer, i.e. it evaluates the
    5mer as seen from the minus strand.  A row set that goes from ~0 to a
    realistic DRACH share under ``revcomp=True`` is proof that the coordinates
    are fine and only the *strand* was missing -- the diagnostic used to tell
    "no strand info" apart from "coordinate bug".
    """
    arr = np.asarray(list(five_mer)).astype(str)
    if revcomp:
        arr = np.array([s.translate(_RC)[::-1] if len(s) == 5 else "NNNNN"
                        for s in arr], dtype=object)
    ok = np.ones(len(arr), dtype=bool)
    for i, allowed in enumerate(DRACH_ALLOWED):
        col = np.array([s[i] if len(s) == 5 else "N" for s in arr])
        ok &= np.isin(col, list(allowed))
    return ok


def drach_rate(five_mer: Iterable[str], revcomp: bool = False) -> float:
    m = drach_mask(five_mer, revcomp=revcomp)
    return float(m.mean()) if len(m) else float("nan")


def centre_base_report(df: pd.DataFrame, mod_type: str) -> dict:
    """Per-callset summary used by the audit and the hard filter.

    Returns the strict / loose shares, the unknown-strand share and the two
    DRACH views (as-is vs reverse-complemented) for the rows that fail the
    strict check -- the pair that separates "missing strand" from "wrong
    coordinate".
    """
    n = len(df)
    out = {"n_calls": n, "mod_type": mod_type,
           "frac_strict": np.nan, "frac_loose": np.nan, "frac_unknown_strand": np.nan,
           "frac_off_base": np.nan, "drach_failed": np.nan, "drach_failed_revcomp": np.nan,
           "dist_center_hist": ""}
    if n == 0:
        return out
    strand = df["strand"].astype(str) if "strand" in df else pd.Series([""] * n)
    ref = df["ref_base"].astype(str) if "ref_base" in df else pd.Series([""] * n)
    st = strict_status(strand.to_numpy(), ref.to_numpy(), mod_type)
    out["frac_strict"] = float((st == OK).mean())
    out["frac_unknown_strand"] = float((st == UNKNOWN_STRAND).mean())
    out["frac_off_base"] = float((st == OFF_BASE).mean())
    pair = expected_base_pair(mod_type)
    out["frac_loose"] = float(ref.isin(pair).mean()) if pair else np.nan
    fail = st == OFF_BASE
    if fail.any() and "five_mer_raw" in df.columns:
        sub = df.loc[fail, "five_mer_raw"].astype(str)
        out["drach_failed"] = drach_rate(sub)
        out["drach_failed_revcomp"] = drach_rate(sub, revcomp=True)
    if "dist_center" in df.columns:
        d = pd.to_numeric(df["dist_center"], errors="coerce")
        d = d[fail].dropna().astype(int)
        if len(d):
            counts = d.value_counts().sort_index().head(12)
            out["dist_center_hist"] = ",".join(f"{k}:{v}" for k, v in counts.items())
    return out
