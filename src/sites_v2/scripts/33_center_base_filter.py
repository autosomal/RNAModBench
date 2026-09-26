#!/usr/bin/env python3
"""33 -- hard-filter every callset down to chemically possible centre bases.

Rule (strict, chain-aware)
--------------------------
A call is kept only when the reference base at ``pos_raw`` is the base the
modification can sit on, *in the orientation of its transcript*:

    m6A     -> genome 'A' on '+', 'T' on '-'
    m5C     -> 'C' / 'G'
    Psi, m1Psi -> 'T' / 'A'      ({U -> T} DNA alphabet)
    inosine -> 'A' / 'T'
    Nm      -> any base (kept, ``no_reference_base``)

Everything else is deleted from the callset: ``off_base`` (the reference simply
does not carry that base) and ``unknown_strand`` (the verdict is impossible to
make -- after ``32_impute_strand`` those are ~5 rows out of 4.48 M).

This is *not* the ``center_base_expected`` column: that one searched +/- 5 bp
and was therefore ~100 % everywhere.  It is also stricter than
``29_anchor_audit``'s loose check, which accepted the complement of the expected
base regardless of strand and would have passed a callset that is wrong on half
of its minus-strand rows.

Order matters: ``32_impute_strand`` -> ``03_annotate_callsets`` -> **33**.

Outputs
-------
sites_v2/evaluation/tables/center_base_strict_audit.tsv  per-callset audit
sites_v2/evaluation/tables/center_filter_removed.csv     every deleted row
sites_v2/manifest/center_base_filter.csv                 per-callset counts

Usage
-----
conda activate benchmark-revision
python .../33_center_base_filter.py [--apply] [--sample S] [--tool T]

Without ``--apply`` nothing is written to the callsets (audit + removal list
only), so the damage can be inspected before it is done.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.center import (OFF_BASE, OK, UNKNOWN_STRAND, centre_base_report,  # noqa: E402
                           strict_status)
from common.config import CALLSET_ROOT, MANIFEST_DIR, TABLE_DIR  # noqa: E402
from common.io_utils import read_table, write_table  # noqa: E402
from common.manifest import Inventory, log_time, setup_logger  # noqa: E402

#: columns copied into the removal list so a deleted row can be re-identified
_KEEP_COLS = ("chrom", "pos_raw", "strand", "ref_base", "five_mer_raw",
              "mod_ratio", "score", "coverage")


def iter_callsets():
    for f in sorted(CALLSET_ROOT.rglob("*.tsv")):
        parts = f.relative_to(CALLSET_ROOT).parts
        if len(parts) != 6:
            continue
        platform, species, group, mod_type, tool, fname = parts
        yield f, platform, species, group, mod_type, tool, fname[:-4]


def keep_mask(df: pd.DataFrame, mod_type: str) -> pd.Series:
    """Boolean mask of rows to keep.

    Keep ``ok`` **and** ``no_expectation``: a modification without a single-base
    expectation (Nm) has nothing to be wrong against -- 33's first run wrongly
    deleted all 18,698 NanoNm rows by accepting only ``ok``; Nm is exempt.
    """
    if "center_status" in df.columns:
        return df["center_status"].astype(str).isin(["ok", "no_expectation"])
    code = strict_status(df["strand"].astype(str).to_numpy(),
                         df["ref_base"].astype(str).to_numpy(), mod_type)
    return pd.Series((code == OK) | (code == 3), index=df.index)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true",
                    help="rewrite the callsets (default: audit + removal list only)")
    ap.add_argument("--sample", action="append", default=None)
    ap.add_argument("--tool", action="append", default=None)
    args = ap.parse_args()

    logger = setup_logger("33_center_base_filter")
    inv = Inventory("33_center_base_filter")
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    audit_rows, removal_frames = [], []
    with log_time(logger, "centre-base filter" + (" (apply)" if args.apply else " (dry-run)")):
        for path, platform, species, group, mod_type, tool, sample in iter_callsets():
            if args.sample and sample not in args.sample:
                continue
            if args.tool and tool not in args.tool:
                continue
            df = read_table(path)
            rec = {"platform": platform, "species": species, "dataset_group": group,
                   "mod_type": mod_type, "tool": tool, "sample": sample,
                   "n_before": len(df), "n_kept": 0, "n_off_base": 0,
                   "n_unknown_strand": 0, "n_no_expectation": 0, "verdict": "empty"}
            if df.empty:
                audit_rows.append(rec)
                continue

            rep = centre_base_report(df, mod_type)
            keep = keep_mask(df, mod_type)
            code = strict_status(df["strand"].astype(str).to_numpy(),
                                 df["ref_base"].astype(str).to_numpy(), mod_type)
            rec["n_kept"] = int(keep.sum())
            rec["n_off_base"] = int((code == OFF_BASE).sum())
            rec["n_unknown_strand"] = int((code == UNKNOWN_STRAND).sum())
            rec["n_no_expectation"] = int((code == 3).sum())
            rec["frac_strict"] = rep["frac_strict"]
            rec["frac_loose"] = rep["frac_loose"]
            rec["frac_off_base"] = rep["frac_off_base"]
            rec["frac_unknown_strand"] = rep["frac_unknown_strand"]
            rec["drach_failed"] = rep["drach_failed"]
            rec["drach_failed_revcomp"] = rep["drach_failed_revcomp"]
            rec["dist_center_hist_failed"] = rep["dist_center_hist"]
            rec["verdict"] = ("no_reference_base" if rec["n_no_expectation"] == len(df)
                              else "ok" if rec["n_kept"] == len(df)
                              else "filtered")
            rec["pct_removed"] = round(100 * (len(df) - rec["n_kept"]) / len(df), 4)

            dropped = df.loc[~keep]
            if len(dropped):
                cols = [c for c in _KEEP_COLS if c in dropped.columns]
                rem = dropped[cols].copy()
                rem.insert(0, "sample", sample)
                rem.insert(0, "tool", tool)
                rem.insert(0, "mod_type", mod_type)
                rem.insert(0, "dataset_group", group)
                rem.insert(0, "species", species)
                rem.insert(0, "platform", platform)
                rem["reason"] = np.where(
                    dropped["strand"].astype(str).isin(["+", "-"]),
                    "off_base", "unknown_strand")
                removal_frames.append(rem)

            if args.apply and len(dropped):
                write_table(df.loc[keep].reset_index(drop=True), path)
            if rec["verdict"] == "filtered":
                logger.info("  %-22s %-30s -%-7d (%.2f%%) off=%d unk=%d",
                            sample, tool, len(dropped), rec["pct_removed"],
                            rec["n_off_base"], rec["n_unknown_strand"])
            audit_rows.append(rec)

    adf = pd.DataFrame(audit_rows)
    write_table(adf, TABLE_DIR / "center_base_strict_audit.tsv")
    inv.record(TABLE_DIR / "center_base_strict_audit.tsv", n_rows=len(adf))

    if removal_frames:
        rdf = pd.concat(removal_frames, ignore_index=True)
        write_table(rdf, TABLE_DIR / "center_filter_removed.csv")
        inv.record(TABLE_DIR / "center_filter_removed.csv", n_rows=len(rdf))
    else:
        rdf = pd.DataFrame()
        write_table(rdf, TABLE_DIR / "center_filter_removed.csv")
        inv.record(TABLE_DIR / "center_filter_removed.csv", n_rows=0)

    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    cpath = MANIFEST_DIR / "center_base_filter.csv"
    write_table(adf, cpath)
    inv.record(cpath, n_rows=len(adf))
    inv.flush()

    if len(adf):
        n_before = int(adf["n_before"].sum())
        n_kept = int(adf["n_kept"].sum())
        logger.info("rows %d -> %d (removed %d, %.3f%%)",
                    n_before, n_kept, n_before - n_kept,
                    100 * (n_before - n_kept) / max(n_before, 1))
        logger.info("verdicts:\n%s", adf["verdict"].value_counts().to_string())
        worst = adf[adf["n_kept"] < adf["n_before"]].sort_values(
            "pct_removed", ascending=False).head(20)
        if len(worst):
            logger.warning("hardest hit:\n%s",
                           worst[["sample", "tool", "mod_type", "n_before",
                                  "n_kept", "n_off_base", "n_unknown_strand",
                                  "pct_removed"]].to_string(index=False))
    logger.info("-> %s", TABLE_DIR / "center_base_strict_audit.tsv")
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
