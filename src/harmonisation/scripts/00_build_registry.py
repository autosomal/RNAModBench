#!/usr/bin/env python3
"""00 -- build the harmonisation sample / tool registry and classify gaps.

Outputs
-------
harmonisation/manifest/sample_registry.csv      one row per canonical sample
harmonisation/manifest/sample_tool_registry.csv one row per (sample, tool, mod_type)
harmonisation/manifest/pending.csv              expected-but-missing pairs (legacy cross-check)

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/00_build_registry.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))  # code/harmonisation

from common.config import (CALLSET_ROOT, EVALUATION_ROOT, LOG_DIR, MANIFEST_DIR,
                           SITES_ROOT, UNIVERSE_ROOT)
from common.io_utils import ensure_dirs, write_table
from common.manifest import Inventory, log_time, setup_logger
from common.registry import build_registry


def main() -> None:
    ensure_dirs(SITES_ROOT, CALLSET_ROOT, UNIVERSE_ROOT, EVALUATION_ROOT,
                MANIFEST_DIR, LOG_DIR, EVALUATION_ROOT / "tables")
    logger = setup_logger("00_build_registry")
    inv = Inventory("00_build_registry")

    with log_time(logger, "registry build"):
        samples, pairs, pending = build_registry()

    write_table(samples, MANIFEST_DIR / "sample_registry.csv")
    inv.record(MANIFEST_DIR / "sample_registry.csv", n_rows=len(samples))
    write_table(pairs, MANIFEST_DIR / "sample_tool_registry.csv")
    inv.record(MANIFEST_DIR / "sample_tool_registry.csv", n_rows=len(pairs))
    write_table(pending, MANIFEST_DIR / "pending.csv")
    inv.record(MANIFEST_DIR / "pending.csv", n_rows=len(pending))

    logger.info("samples: %d (platforms: %s)", len(samples),
                samples["platform"].value_counts().to_dict())
    logger.info("(sample, tool) pairs: %d", len(pairs))
    logger.info("status counts:\n%s", pairs["status"].value_counts().to_string())
    ok = pairs[pairs["status"] == "ok"]
    logger.info("ok pairs by species x platform:\n%s",
                ok.groupby(["platform", "species"]).size().to_string())
    logger.info("ok pairs by tool (RNA002):\n%s",
                ok[ok["platform"] == "RNA002"].groupby("tool").size().to_string())
    logger.info("ok pairs by tool (RNA004):\n%s",
                ok[ok["platform"] == "RNA004"].groupby("tool").size().to_string())

    amb = pairs[pairs["ambiguous_dir"] == 1]
    if not amb.empty:
        logger.warning("ambiguous directory mappings (%d rows):", len(amb))
        for _, r in amb.drop_duplicates(["sample", "tool"]).iterrows():
            logger.warning("  %s x %s -> %s", r["sample"], r["tool"], r["note"])

    missing = pairs[pairs["status"] == "missing"]
    logger.warning("missing pairs: %d (see sample_tool_registry.csv)", len(missing))
    if not pending.empty:
        logger.warning("pending vs legacy copy map: %d rows -> %s",
                       len(pending), MANIFEST_DIR / "pending.csv")
        for _, r in pending.head(60).iterrows():
            logger.warning("  PENDING %-22s %-16s %s", r["sample"], r["tool"], r["reason"])

    inv.flush()
    logger.info("registry written to %s", MANIFEST_DIR)
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
