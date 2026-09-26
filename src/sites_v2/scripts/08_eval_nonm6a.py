#!/usr/bin/env python3
"""08 -- synthetic Curlcake ground truth and the HeLa non-m6A panel.

The non-m6A tools (m5C / Psi / m1Psi / Nm) were only ever run on HeLa and the
Curkcake constructs.  This script therefore evaluates them where a reference
exists and reports them as "negative-control + reproducibility" elsewhere:

* **Curkcake fully m6A-modified constructs** (Curlcake_m6A_rep2, Curlcake_m6A_rep1):
  ground truth = all DRACH A sites of ``cc.fasta`` (both strands).  TP/FP/FN/TN
  are computed inside the sample's universe exactly like the GLORI analysis.
* **HeLa non-m6A (WT and IVT)**: no reference exists; the table reports the
  call burden per 10^6 candidate sites, the WT/IVT ratio on the shared universe
  and the cross-replicate reproducibility (Jaccard) for m5C/Psi/m1Psi/Nm tools.
* **Curkcake non-m6A positives** are not available (the modified constructs
  carry m6A only) - recorded as a gap.

Outputs
-------
sites_v2/evaluation/tables/curlcake_truth.tsv
sites_v2/evaluation/tables/hela_nonm6a.tsv

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/sites_v2/scripts/08_eval_nonm6a.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import (CALLSET_ROOT, CURLCAKE_FASTA, PRIMARY_WINDOW,
                           SAMPLES_BY_NAME, TABLE_DIR, UNIVERSE_ROOT, WINDOWS)
from common.evaluation import (confusion_windows, load_universe,
                               positions_in_universe, reference_in_universe)
from common.io_utils import read_table, write_table
from common.manifest import Inventory, log_time, setup_logger
from common.refs import curlcake_drach_a

NONM6A_GROUPS = {"HeLa_WT", "HeLa_IVT"}


def universe_path(sample) -> Path:
    plain = UNIVERSE_ROOT / sample.platform / sample.species / f"{sample.canonical}__universe.tsv"
    gz = plain.with_suffix(".tsv.gz")
    return plain if plain.exists() else gz


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-cov", type=int, default=10)
    args = ap.parse_args()

    logger = setup_logger("08_eval_nonm6a")
    inv = Inventory("08_eval_nonm6a")
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    truth_rows, hela_rows = [], []

    with log_time(logger, "non-m6A evaluation"):
        # ---------------------------------------------------------- Curlcake --
        cc_truth = curlcake_drach_a(CURLCAKE_FASTA)
        for sample_name in ("Curlcake_m6A_rep2", "Curlcake_m6A_rep1"):
            sample = SAMPLES_BY_NAME.get(sample_name)
            if sample is None:
                continue
            up = universe_path(sample)
            if not up.exists():
                logger.warning("[%s] no universe, skipped", sample_name)
                continue
            for tool_dir in sorted((CALLSET_ROOT / sample.platform / sample.species /
                                    sample.dataset_group).glob("*/*")):
                mod_type = tool_dir.parent.name
                tool = tool_dir.name
                path = tool_dir / f"{sample_name}.tsv"
                if not path.exists():
                    continue
                universe = load_universe(up, mod_type, args.min_cov)
                ref = reference_in_universe(cc_truth, universe)
                df = read_table(path)
                df["pos_raw"] = pd.to_numeric(df["pos_raw"]).astype(np.int64)
                calls_u, n_out = positions_in_universe(df, universe)
                conf = confusion_windows(universe, ref, calls_u, WINDOWS)
                for w, c in conf.items():
                    truth_rows.append({
                        "sample": sample_name, "dataset_group": sample.dataset_group,
                        "mod_type": mod_type, "tool": tool, "window": w,
                        "n_calls": len(df), "n_calls_in_universe": c["n_calls_in_universe"],
                        "n_calls_out_of_universe": n_out,
                        "n_truth_in_universe": c["n_reference_in_universe"],
                        "tp": c["tp"], "fp": c["fp"], "fn": c["fn"], "tn": c["tn"],
                        "precision": c["precision"], "recall": c["recall"],
                        "f1": c["f1"], "mcc": c["mcc"],
                        "universe": c["universe"], "min_cov": args.min_cov,
                    })
                logger.info("  curlcake %-20s %-24s calls=%d", sample_name, tool, len(df))
        if truth_rows:
            df_truth = pd.DataFrame(truth_rows)
            write_table(df_truth, TABLE_DIR / "curlcake_truth.tsv")
            inv.record(TABLE_DIR / "curlcake_truth.tsv", n_rows=len(df_truth))

        # ------------------------------------------------------------- HeLa ----
        for group in sorted(NONM6A_GROUPS):
            samples = [s for s in SAMPLES_BY_NAME.values()
                       if s.dataset_group == group and s.species == "Human"]
            for s in samples:
                for tool_dir in sorted((CALLSET_ROOT / s.platform / s.species /
                                        group).glob("*/*")):
                    mod_type = tool_dir.parent.name
                    if mod_type == "m6A":
                        continue
                    tool = tool_dir.name
                    path = tool_dir / f"{s.canonical}.tsv"
                    if not path.exists():
                        continue
                    up = universe_path(s)
                    universe = load_universe(up, mod_type, args.min_cov) if up.exists() else {}
                    n_u = int(sum(v.size for v in universe.values()))
                    df = read_table(path)
                    df["pos_raw"] = pd.to_numeric(df["pos_raw"]).astype(np.int64)
                    calls_u, n_out = positions_in_universe(df, universe) if universe else ({}, len(df))
                    n_in = int(sum(v.size for v in calls_u.values()))
                    hela_rows.append({
                        "sample": s.canonical, "dataset_group": group,
                        "condition_class": s.condition_class, "replicate_tag": s.replicate_tag,
                        "mod_type": mod_type, "tool": tool,
                        "n_calls": len(df), "n_calls_in_universe": n_in,
                        "n_calls_out_of_universe": n_out, "n_universe": n_u,
                        "calls_per_1e6_candidates": 1e6 * n_in / n_u if n_u else np.nan,
                        "pct_glori_undefined": np.nan,
                        "note": "no reference available for this modification type",
                    })
        if hela_rows:
            hela = pd.DataFrame(hela_rows)
            # cross-replicate reproducibility (mean pairwise Jaccard) per tool
            jac_rows = []
            for (mod_type, tool), sub in hela.groupby(["mod_type", "tool"]):
                wt = sub[sub["condition_class"] == "WT"]
                ivt = sub[sub["condition_class"] == "IVT"]
                for label, part in (("WT", wt), ("IVT", ivt)):
                    if len(part) < 2:
                        continue
                    sets = []
                    for _, r in part.iterrows():
                        s = SAMPLES_BY_NAME[r["sample"]]
                        up = universe_path(s)
                        u = load_universe(up, mod_type, args.min_cov) if up.exists() else {}
                        p = CALLSET_ROOT / s.platform / s.species / s.dataset_group / \
                            mod_type / tool / f"{s.canonical}.tsv"
                        d = read_table(p)
                        d["pos_raw"] = pd.to_numeric(d["pos_raw"]).astype(np.int64)
                        c, _ = positions_in_universe(d, u) if u else ({}, 0)
                        sets.append(c)
                    vals = []
                    for i in range(len(sets)):
                        for j in range(i + 1, len(sets)):
                            a, b = sets[i], sets[j]
                            ch = sorted(set(a) | set(b))
                            inter = sum(np.intersect1d(a.get(c, np.zeros(0, np.int64)),
                                                       b.get(c, np.zeros(0, np.int64))).size for c in ch)
                            union = sum(np.union1d(a.get(c, np.zeros(0, np.int64)),
                                                   b.get(c, np.zeros(0, np.int64))).size for c in ch)
                            if union:
                                vals.append(inter / union)
                    jac_rows.append({"mod_type": mod_type, "tool": tool,
                                     "panel": label, "n_replicates": len(sets),
                                     "mean_pairwise_jaccard": float(np.mean(vals)) if vals else np.nan})
            if jac_rows:
                hela = hela.merge(pd.DataFrame(jac_rows),
                                  on=["mod_type", "tool"], how="left")
            write_table(hela, TABLE_DIR / "hela_nonm6a.tsv")
            inv.record(TABLE_DIR / "hela_nonm6a.tsv", n_rows=len(hela))
            logger.info("HeLa non-m6A panel: %d rows", len(hela))

    inv.flush()
    logger.info("tables -> %s", TABLE_DIR)
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
