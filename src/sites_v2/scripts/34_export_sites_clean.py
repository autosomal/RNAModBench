#!/usr/bin/env python3
"""34 -- export a CLEAN site collection from ``sites_v2/callsets``.

Why (2026-09-18, user request)
------------------------------
One tidy table per callset that keeps **every informative callset column**
with a clear per-tool meaning, dropping only true redundancy:

* ``mod_ratio`` is dropped for tools where it is numerically identical to
  ``score`` (ratio-type tools) or empty -- it survives only where it carries
  a *second* value (CHEUI: Prob vs stoichiometry; m6Anet: probability vs
  mod_ratio);
* ``frac_diff`` exists only for xPore (KO-vs-WT differential mod rate),
  ``mod_label``/``src_coverage`` only for Dorado/modkit -- those columns
  disappear automatically in files where they are empty;
* the score semantics differ per tool, so ``score_type`` travels with
  ``score`` (Prob | mod_ratio | m6a_ratio | percent_modified/100 | FDR |
  adjPval | Pvalue | delta_sum_err).

Not exported (constant per file or guaranteed by the pipeline, documented in
``sites_clean/README.md`` instead): identity columns (path), strand_src,
cov_source semantics, dist_center/pos_center (all 0 after the 33 hard
filter), and the parser ``note`` bookkeeping.

Usage
-----
conda activate benchmark-revision
python .../34_export_sites_clean.py
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
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import CALLSET_ROOT, SITES_ROOT  # noqa: E402

OUT_DIR = (_RB / "data/sites_clean")

#: fixed base order; conditional columns (mod_ratio / frac_diff / mod_label /
#: src_coverage / offset_flag) are inserted after ``score`` / ``coverage``
#: respectively when they carry information for that callset.
COLUMNS = ["chrom", "start", "end", "strand", "score", "score_type",
           "mod_ratio", "frac_diff", "mod_label", "src_coverage", "coverage",
           "ref_base", "five_mer_raw", "center_status",
           "offset_flag", "drach"]

READ_COLS = {"chrom", "pos_raw", "strand", "score", "score_type", "mod_ratio",
             "frac_diff", "mod_label", "src_coverage", "coverage", "cov_source",
             "ref_base", "five_mer_raw", "center_status", "offset_flag",
             "is_drach_raw"}


def main() -> None:
    frames = []
    for f in sorted(CALLSET_ROOT.rglob("*.tsv")):
        parts = f.relative_to(CALLSET_ROOT).parts
        if len(parts) != 6:
            continue
        platform, species, group, mod_type, tool, fname = parts
        try:
            df = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False,
                             usecols=lambda c: c in READ_COLS)
        except pd.errors.EmptyDataError:
            df = pd.DataFrame(columns=sorted(READ_COLS))
        base_cols = [c for c in COLUMNS
                     if c not in ("mod_ratio", "frac_diff", "mod_label",
                                  "src_coverage", "offset_flag", "drach")]
        if mod_type == "m6A":
            base_cols.append("drach")
        if df.empty:
            # empty callset -> header-only placeholder so the export mirrors
            # callsets one-to-one (user request 2026-09-18)
            dest = OUT_DIR / f.relative_to(CALLSET_ROOT)
            dest.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(columns=base_cols).to_csv(dest, sep="\t", index=False)
            continue

        def num(c):
            return pd.to_numeric(df[c], errors="coerce") if c in df.columns \
                else pd.Series(np.nan, index=df.index, dtype="Float64")

        start = num("pos_raw").astype("Int64")
        chrom = df["chrom"].astype(str)
        if species == "Curlcake":
            # fix_chromosome's generic "chr*" normalisation produced the ugly
            # internal name chrcurlcake1 -- the native contig (FASTA headers,
            # exon annotation) is Curlcake1..4
            chrom = chrom.str.replace("^chrcurlcake", "Curlcake", regex=True)
        out = pd.DataFrame({
            "chrom": chrom,
            "start": start,
            "end": start + 1,
            "strand": df["strand"].astype(str),
            "score": df["score"].astype(str) if "score" in df.columns else "",
            "score_type": df["score_type"].astype(str) if "score_type" in df.columns else "",
            "mod_ratio": num("mod_ratio"),
            "frac_diff": num("frac_diff"),
            "mod_label": df["mod_label"].astype(str) if "mod_label" in df.columns else "",
            "src_coverage": num("src_coverage").astype("Int64"),
            "coverage": num("coverage").astype("Int64"),
            "cov_source": df["cov_source"].astype(str) if "cov_source" in df.columns else "",
            "ref_base": df["ref_base"].astype(str) if "ref_base" in df.columns else "",
            "five_mer_raw": df["five_mer_raw"].astype(str) if "five_mer_raw" in df.columns else "",
            "center_status": df["center_status"].astype(str) if "center_status" in df.columns else "",
            "offset_flag": df["offset_flag"].astype(str) if "offset_flag" in df.columns else "",
        })
        cols = list(COLUMNS)
        cols.remove("drach")
        if mod_type == "m6A":
            out["drach"] = (df["is_drach_raw"].astype(str).str.lower().eq("true")
                            if "is_drach_raw" in df.columns else pd.NA)
            cols.append("drach")

        # --- de-redundancy ------------------------------------------------
        s = pd.to_numeric(out["score"], errors="coerce")
        # mod_ratio only where it is a *second* value (CHEUI, m6Anet)
        mr = pd.to_numeric(out["mod_ratio"], errors="coerce")
        if mr.notna().any():
            same = (s.round(4) == mr.round(4)) | (mr.isna() & s.isna())
            if same.all():
                cols.remove("mod_ratio")
        else:
            cols.remove("mod_ratio")
        # drop conditional columns that are empty in this callset
        # ("" is not NaN for object dtype -- test string emptiness explicitly)
        for c in ("frac_diff", "src_coverage"):
            if c in cols and not pd.to_numeric(out[c], errors="coerce").notna().any():
                cols.remove(c)
        for c in ("mod_label", "offset_flag"):
            if c in cols and not out[c].astype(str).str.len().gt(0).any():
                cols.remove(c)
        # center_status: keep, but shorten the label for the export
        out["center_status"] = out["center_status"].replace(
            {"no_reference_base": "no_ref_base"})
        out = out.dropna(subset=["start"])
        frames.append((f.relative_to(CALLSET_ROOT), out[cols]))

    total = 0
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for rel, sub in frames:
        dest = OUT_DIR / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        sub.to_csv(dest, sep="\t", index=False)
        total += len(sub)
    print(f"sites_clean: {len(frames)} files, {total:,} rows -> {OUT_DIR}")


if __name__ == "__main__":
    main()
