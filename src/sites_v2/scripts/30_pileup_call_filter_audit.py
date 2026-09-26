#!/usr/bin/env python3
"""30 -- pileup call-filter audit report-only + hard assertion : no no-call row may reach a
Dorado callset, and every source must declare the same call filter.

Motivation (2026-09-18)
-----------------------
``modkit pileup`` emits one row per *covered position* x *modification code*,
including rows where nothing was called (``percent_modified = 0``).  The HeLa
RNA004 family was filtered at source (``dorado_model_split/``: every row has
``valid_coverage >= 20`` and ``percent_modified > 0``), but the Curlcake family
(``RNA004_result/dorado_model/*_pileup.bed``) was not: there the no-call rows are
~60-99 % of every file, so the "callsets" collapsed onto the construct's base
composition -- m6A / m5C / Psi all looked like the expected base only ~58 % of the
time, which is exactly a random position set.  Once the declared call filter
(``config.RNA004_SOURCES`` 5th element -> ``parse_dorado_pileup(min_cov, min_pct)``)
is applied, the same channels sit at 91-97 % (Psi 63-79 %, in line with HeLa's own
pseU channels at 72-90 %).

This script audits the contract:

* every declared Dorado source: rows, no-call rows (``percent_modified == 0``),
  rows below the declared coverage floor, the declared filter, and the retained
  share;
* every Dorado callset under ``callsets/RNA004/``: rows with ``mod_ratio == 0``
  (must be **0**) and rows below the coverage floor of their source family;
* exit code 1 when a callset still carries no-call rows, so the pipeline can
  never silently regress to the 2026-09-18 state.

Outputs
-------
sites_v2/evaluation/tables/pileup_call_filter_audit.tsv

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/sites_v2/scripts/30_pileup_call_filter_audit.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import CALLSET_ROOT, SAMPLES_BY_NAME, TABLE_DIR  # noqa: E402
from common.io_utils import read_table, write_table  # noqa: E402
from common.manifest import Inventory, log_time, setup_logger  # noqa: E402
from common.registry import resolve_rna004_sources  # noqa: E402

RNA004_CALLSETS = CALLSET_ROOT / "RNA004"
SAMPLES = ("HeLa_RNA004_WT", "HeLa_RNA004_IVT", "Curlcake_RNA004_IVT")


def _cols(df: pd.DataFrame) -> tuple[str | None, str | None, str | None]:
    """(coverage, percent-modified, name) column names as written by modkit."""
    def find(cands):
        for c in cands:
            if c in df.columns:
                return c
        return None
    return (find(["valid_coverage", "coverage"]),
            find(["percent_modified", "pct_modified"]),
            find(["name", "mod"]))


def audit_source(path: Path, declared: dict) -> dict:
    df = read_table(path)
    c_cov, c_pct, c_name = _cols(df)
    n = len(df)
    stats = {"kind": "source", "sample_group": "", "tool": path.name,
             "mod_type": "", "n_rows": n, "n_no_call": "", "n_below_floor": "",
             "declared_filter": json.dumps(declared, sort_keys=True) if declared else "",
             "retained_frac": "", "violation": ""}
    if n == 0:
        return stats
    if c_pct:
        stats["n_no_call"] = int((pd.to_numeric(df[c_pct], errors="coerce")
                                  .fillna(0) <= 0).sum())
    min_cov = declared.get("min_cov")
    if c_cov and min_cov is not None:
        cov = pd.to_numeric(df[c_cov], errors="coerce").fillna(0)
        stats["n_below_floor"] = int((cov < min_cov).sum())
    if c_name:
        codes = sorted(set(df[c_name].astype(str).str.split("#").str[0]))
        stats["mod_type"] = "|".join(codes)
    # retained share under the declared filter (what the extractor will keep)
    keep = pd.Series(True, index=df.index)
    if c_cov and min_cov is not None:
        keep &= pd.to_numeric(df[c_cov], errors="coerce").fillna(0) >= min_cov
    if c_pct and declared.get("min_pct") is not None:
        keep &= pd.to_numeric(df[c_pct], errors="coerce").fillna(0) > declared["min_pct"]
    stats["retained_frac"] = round(float(keep.mean()), 4)
    return stats


def audit_callset(path: Path, group: str) -> dict:
    df = read_table(path)
    n = len(df)
    rec = {"kind": "callset", "sample_group": group, "tool": path.parent.name,
           "mod_type": path.parent.parent.name, "n_rows": n,
           "n_no_call": "", "n_below_floor": "", "declared_filter": "",
           "retained_frac": "", "n_off_base": "", "violation": ""}
    if n == 0:
        return rec
    ratio_col = "mod_ratio" if "mod_ratio" in df.columns else ("score" if "score" in df.columns else None)
    if ratio_col:
        ratio = pd.to_numeric(df[ratio_col], errors="coerce").fillna(0)
        rec["n_no_call"] = int((ratio <= 0).sum())
        if rec["n_no_call"]:
            rec["violation"] = f"{rec['n_no_call']} rows with {ratio_col} == 0 (no-call)"
    if "src_coverage" in df.columns:
        cov = pd.to_numeric(df["src_coverage"], errors="coerce").fillna(0)
        rec["n_below_floor"] = int((cov < 20).sum())
        if rec["n_below_floor"]:
            rec["violation"] = (rec["violation"] + "; " if rec["violation"] else "") + \
                f"{rec['n_below_floor']} rows with valid_coverage < 20"
    # reference-base contract (2026-09-18): modkit only calls a modification on
    # the base it belongs to, so a Dorado row that is not on that base is either
    # a coordinate defect or a no-call that slipped through.  ``33`` deletes
    # them; this assertion makes sure it did.
    if "center_status" in df.columns:
        status = df["center_status"].astype(str)
        off = int((~status.isin(["ok", "no_expectation"])).sum())
        rec["n_off_base"] = off
        if off:
            rec["violation"] = (rec["violation"] + "; " if rec["violation"] else "") + \
                f"{off} rows not on the modification's reference base"
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-fail", action="store_true",
                    help="report only, never exit non-zero")
    args = ap.parse_args()

    logger = setup_logger("30_pileup_call_filter_audit")
    inv = Inventory("30_pileup_call_filter_audit")
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    with log_time(logger, "pileup call-filter audit"):
        for sample_name in SAMPLES:
            sample = SAMPLES_BY_NAME.get(sample_name)
            if sample is None:
                continue
            for src in resolve_rna004_sources(sample):
                if "dorado_pileup" != src.parser or src.path is None:
                    continue
                declared = json.loads(src.parser_arg) if src.parser_arg else {}
                rec = audit_source(src.path, declared)
                rec["sample_group"] = sample_name
                rec["tool"] = src.tool
                rec["mod_type"] = src.mod_type
                rows.append(rec)

        for f in sorted(RNA004_CALLSETS.rglob("*.tsv")):
            parts = f.relative_to(RNA004_CALLSETS).parts
            if len(parts) != 5 or not parts[3].startswith("Dorado"):
                continue
            rows.append(audit_callset(f, f"{parts[0]}/{parts[1]}"))

    df = pd.DataFrame(rows)
    out_path = TABLE_DIR / "pileup_call_filter_audit.tsv"
    write_table(df, out_path)
    inv.record(out_path, n_rows=len(df))
    inv.flush()

    if len(df):
        src = df[df["kind"] == "source"]
        cal = df[df["kind"] == "callset"]
        if len(src):
            logger.info("sources: %d | no-call rows %s / %s total rows",
                        len(src),
                        int(pd.to_numeric(src["n_no_call"], errors="coerce").fillna(0).sum()),
                        int(pd.to_numeric(src["n_rows"], errors="coerce").fillna(0).sum()))
        if len(cal):
            logger.info("callsets audited: %d", len(cal))
        bad = df[df["violation"].astype(str) != ""]
        if len(bad):
            logger.error("VIOLATIONS in %d entries:\n%s", len(bad),
                         bad[["kind", "sample_group", "tool", "n_rows", "violation"]]
                         .to_string(index=False))
        else:
            logger.info("no violation: no callset carries no-call rows")
    logger.info("-> %s", out_path)
    logger.info("log: %s", logger.log_path)

    bad = df[df["violation"].astype(str) != ""] if len(df) else df
    sys.exit(0 if args.no_fail or len(bad) == 0 else 1)


if __name__ == "__main__":
    main()
