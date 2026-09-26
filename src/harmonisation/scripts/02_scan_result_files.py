#!/usr/bin/env python3
"""02 -- scan result/ + result_RNA004/ for per-replicate files and flag anything
the extraction patterns would miss.

Rationale: the legacy per-tool naming is inconsistent (``_remove_chr.txt``,
``_liftover.txt``, ``_result1_``, ``result1/``, ``*_5mer.txt`` ...).  A tool x
sample pair whose file exists only under a variant spelling would silently drop
out of the rebuild, which is exactly the "missing biological replicates" problem
raised by the reviewer.  This script cross-checks every file against the
configured patterns and writes:

harmonisation/manifest/result_file_inventory.csv   every candidate file + matched pattern
harmonisation/manifest/unmatched_files.csv         files no pattern matched (review list)
harmonisation/manifest/replicate_matrix.csv        species x tool x replicate availability

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/02_scan_result_files.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import (MANIFEST_DIR, RESULT_RNA002, RNA002_TOOLS,
                           SAMPLES_BY_NAME)
from common.io_utils import write_table
from common.manifest import Inventory, setup_logger
from common.registry import canonical_from_dir, curlcake_sample_from_stem

#: file-name markers that identify an extractable callsets (any tool)
CALLSET_MARKERS = ("_remove_chr.txt", "_5mer.txt", ".txt", ".tsv", ".bed")


def main() -> None:
    logger = setup_logger("02_scan_result_files")
    inv = Inventory("02_scan_result_files")
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)

    inv_rows, unmatched = [], []
    for tool in RNA002_TOOLS:
        tool_dir = RESULT_RNA002 / tool.result_subdir
        if not tool_dir.is_dir():
            continue
        for d in sorted(p for p in tool_dir.iterdir() if p.is_dir()):
            sample, amb = canonical_from_dir(d.name)
            files = [f for f in d.rglob("*") if f.is_file() and f.stat().st_size > 0]
            files = [f for f in files if any(m in f.name for m in CALLSET_MARKERS)]
            if not files:
                continue
            for f in files:
                rel_name = str(f.relative_to(d))
                matched = None
                for pat in tool.patterns:
                    if f.match(pat) or Path(rel_name).match(pat):
                        matched = pat
                        break
                if matched is None and sample is None and tool.tool.startswith(("CHEUI", "ELIGOS")):
                    pass
                if matched is None:
                    # stem-based Curlcake resolution
                    cc = curlcake_sample_from_stem(f.stem)
                    if cc is not None:
                        matched = "(curlcake stem)"
                row = {"tool": tool.tool, "dir": d.name, "sample": sample or "",
                       "file": rel_name, "size": f.stat().st_size,
                       "matched_pattern": matched or ""}
                inv_rows.append(row)
                if matched is None:
                    unmatched.append(row)

    inv_df = pd.DataFrame(inv_rows)
    write_table(inv_df, MANIFEST_DIR / "result_file_inventory.csv")
    inv.record(MANIFEST_DIR / "result_file_inventory.csv", n_rows=len(inv_df))
    logger.info("files scanned: %d (in %d tool dirs)", len(inv_df),
                inv_df["tool"].nunique() if len(inv_df) else 0)

    if unmatched:
        udf = pd.DataFrame(unmatched)
        write_table(udf, MANIFEST_DIR / "unmatched_files.csv")
        inv.record(MANIFEST_DIR / "unmatched_files.csv", n_rows=len(udf))
        logger.warning("unmatched files: %d", len(udf))
        for _, r in udf.iterrows():
            logger.warning("  %-14s %-34s %s", r["tool"], r["dir"], r["file"])

    # ---------- per-replicate availability matrix ---------------------------
    rows = []
    if len(inv_df):
        for _, r in inv_df.iterrows():
            sample = r["sample"]
            if not sample:
                continue
            rows.append({"species": SAMPLES_BY_NAME[sample].species,
                         "dataset_group": SAMPLES_BY_NAME[sample].dataset_group,
                         "sample": sample,
                         "replicate_tag": SAMPLES_BY_NAME[sample].replicate_tag,
                         "tool": r["tool"], "dir": r["dir"], "file": r["file"]})
    if rows:
        mat = pd.DataFrame(rows)
        piv = (mat.pivot_table(index=["tool"], columns="sample", values="file",
                               aggfunc="count").fillna(0).astype(int))
        write_table(piv.reset_index(), MANIFEST_DIR / "replicate_matrix.csv")
        inv.record(MANIFEST_DIR / "replicate_matrix.csv", n_rows=len(piv))
        logger.info("replicate availability (raw result/ files):\n%s%s",
                    "\n", piv.to_string())

    inv.flush()
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
