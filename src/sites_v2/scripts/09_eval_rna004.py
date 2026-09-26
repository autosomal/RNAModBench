#!/usr/bin/env python3
"""09 -- RNA004 evaluation: Dorado threshold scan, IVT false positives, DRACH
specificity (R3-8, R3-10, E8).

RNA004 call sets are labelled ``Dorado_<hac|sup>@<version>_<model>[_drachGuitar|
_otherMod]``; the Dorado confidence is ``percent_modified/100`` in the ``score``
column, so every model can be re-thresholded without re-reading the pileups.

Panels
------
* ``rna004_dorado_scan.tsv`` - for each Dorado model x threshold in
  {5,10,20,50}%: calls on the HeLa WT / HeLa IVT / Curlcake IVT samples, FP per
  10^6 candidate sites on the IVT controls, and precision vs GLORI on HeLa WT
  (windows 0 and 2).
* ``rna004_tool_eval.tsv`` - the non-Dorado RNA004 tools (m6Anet, NanoSPA,
  NanoPsu, ELIGOS2_solo/_diff, DRUMMER): calls, IVT FP rates, WT precision and
  the share of calls in a DRACH context (``is_drach_center``).
* ``rna004_drach_specificity.tsv`` - DRACH share per tool x sample.

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/sites_v2/scripts/09_eval_rna004.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import (CALLSET_ROOT, GLORI, PRIMARY_WINDOW, SAMPLES_BY_NAME,
                           TABLE_DIR, UNIVERSE_ROOT, WINDOWS)
from common.evaluation import (confusion_windows, load_reference, load_universe,
                               positions_in_universe, reference_in_universe)
from common.io_utils import read_table, write_table
from common.manifest import Inventory, log_time, setup_logger

DORADO_THRESHOLDS = [5, 10, 20, 50]  # percent_modified


def universe_path(sample) -> Path:
    plain = UNIVERSE_ROOT / sample.platform / sample.species / f"{sample.canonical}__universe.tsv"
    gz = plain.with_suffix(".tsv.gz")
    return plain if plain.exists() else gz


def iter_rna004_callsets():
    base = CALLSET_ROOT / "RNA004"
    for f in sorted(base.rglob("*.tsv")):
        species, group, mod_type, tool, fname = f.relative_to(base).parts
        yield fname[:-4], species, group, mod_type, tool, f


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-cov", type=int, default=10)
    args = ap.parse_args()

    logger = setup_logger("09_eval_rna004")
    inv = Inventory("09_eval_rna004")
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    scan_rows, tool_rows, drach_rows = [], [], []
    cc_truth = None

    with log_time(logger, "RNA004 evaluation"):
        for sample_name, species, group, mod_type, tool, path in iter_rna004_callsets():
            sample = SAMPLES_BY_NAME.get(sample_name)
            if sample is None:
                continue
            up = universe_path(sample)
            if not up.exists():
                logger.warning("[%s] no universe, skipped", sample_name)
                continue
            df = read_table(path)
            if df.empty:
                continue
            df["pos_raw"] = pd.to_numeric(df["pos_raw"]).astype(np.int64)
            df["score"] = pd.to_numeric(df.get("score", np.nan), errors="coerce")
            universe = load_universe(up, mod_type, args.min_cov)
            n_u = int(sum(v.size for v in universe.values()))

            is_ivt = sample.condition_class == "IVT"
            ref = None
            if GLORI.get(sample.species) and mod_type == "m6A":
                ref = reference_in_universe(load_reference(GLORI[sample.species]), universe)

            thresholds = ([t / 100.0 for t in DORADO_THRESHOLDS]
                          if (df["score_type"].astype(str) == "percent_modified/100").all()
                          and df["score"].notna().any() else [None])

            for thr in thresholds:
                sub = df if thr is None else df[df["score"] >= thr]
                calls_u, n_out = positions_in_universe(sub, universe)
                n_in = int(sum(v.size for v in calls_u.values()))
                rec = {
                    "sample": sample_name, "species": species,
                    "condition_class": sample.condition_class,
                    "dataset_group": group, "mod_type": mod_type, "tool": tool,
                    "threshold_pct": None if thr is None else int(round(thr * 100)),
                    "n_calls": len(sub), "n_calls_in_universe": n_in,
                    "n_calls_out_of_universe": n_out, "n_universe": n_u,
                    "calls_per_1e6_candidates": 1e6 * n_in / n_u if n_u else np.nan,
                }
                if ref is not None:
                    conf = confusion_windows(universe, ref, calls_u, [0, PRIMARY_WINDOW])
                    for w, c in conf.items():
                        rec[f"precision_w{w}"] = c["precision"]
                        rec[f"recall_w{w}"] = c["recall"]
                        rec[f"tp_w{w}"] = c["tp"]
                        rec[f"fp_w{w}"] = c["fp"]
                        rec[f"fn_w{w}"] = c["fn"]
                else:
                    rec["fp_per_1e6_candidates"] = 1e6 * n_in / n_u if n_u else np.nan
                if thr is not None:
                    scan_rows.append(rec)
                else:
                    tool_rows.append(rec)

            drach_rows.append({
                "sample": sample_name, "condition_class": sample.condition_class,
                "mod_type": mod_type, "tool": tool, "n_calls": len(df),
                "pct_drach_raw": 100 * df["is_drach_raw"].astype(bool).mean(),
                "pct_drach_center": 100 * df["is_drach_center"].astype(bool).mean(),
                "pct_expected_base_at_raw": 100 * df["center_base_expected"].astype(bool).mean(),
            })

    for name, data in (("rna004_dorado_scan.tsv", scan_rows),
                       ("rna004_tool_eval.tsv", tool_rows),
                       ("rna004_drach_specificity.tsv", drach_rows)):
        if data:
            write_table(pd.DataFrame(data), TABLE_DIR / name)
            inv.record(TABLE_DIR / name, n_rows=len(data))
            logger.info("%s: %d rows", name, len(data))
    inv.flush()
    logger.info("tables -> %s", TABLE_DIR)
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
