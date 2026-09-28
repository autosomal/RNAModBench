#!/usr/bin/env python
"""Figure 5 legacy-vs-revision contrast (for the response letter, not a panel).

The published panel A/B/C numbers came from ``Total_new1/NGS.ipynb`` and
``Total_new1/mod_ratio_hit_rate.ipynb``; the notebook's own WT table
(``Total_new1/comprehensive_metrics_data.csv``, 13 tools x 3 species) is the
authoritative record of those legacy values -- it is read here, never recomputed,
because the legacy recipe has no per-replicate structure to recompute from:

    TP = interval-overlap hits, FP = Total_Detected - TP,
    FN = n_GLORI - TP (every GLORI site assumed measurable),
    TN = genome_total - n_GLORI - FP (genome_total a hard-coded constant),
    one legacy aggregate per species (mouse = one study only).

The table pairs every legacy value with the revision value (mean over the drawn
independent units at the primary window / coverage, from
``fig5a_mean_rank_by_group_tool.tsv``) so the response letter can state exactly
how much of the change is the metric definition rather than the data.

Outputs -> figures/figure5/tables/fig5_legacy_vs_revision.tsv

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/figures/figure5/src/38_legacy_fig5_compare.py
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

HERE = _RB / "src/harmonisation"
sys.path.insert(0, str(HERE))

from common import config as C                                        # noqa: E402
from common.io_utils import write_table                              # noqa: E402
from common.manifest import setup_logger                             # noqa: E402

LEGACY_CSV = ((_XB / "code/next_postprocessing/Total_new1/comprehensive_metrics_data.csv"))
OUT = (_RB / "figures/figure5")
TAB = OUT / "tables"

#: the legacy ``Tools.txt`` labels vs the canonical tool ids
LEGACY_TOOL_ALIAS = {"ELIGOS_diff": "ELIGOS2_diff", "ELIGOS_solo": "ELIGOS2_solo",
                     "Yanocomp": "yanocomp"}
LEGACY_SPECIES = {"Human (HeLa)": "Human"}
METRIC_MAP = {"Precision_WT": "precision", "Recall_WT": "recall",
              "F1_WT": "f1", "MCC_WT": "mcc"}


def main() -> None:
    logger = setup_logger("38_legacy_fig5_compare")
    TAB.mkdir(parents=True, exist_ok=True)

    legacy = pd.read_csv(LEGACY_CSV)
    rows = []
    for _, r in legacy.iterrows():
        species = LEGACY_SPECIES.get(str(r["Species"]), str(r["Species"]))
        tool = LEGACY_TOOL_ALIAS.get(str(r["Tools"]), str(r["Tools"]))
        for col, metric in METRIC_MAP.items():
            rows.append({"species": species, "tool": tool, "metric": metric,
                         "legacy_value": float(r[col]),
                         "legacy_total_detected": int(r["Total_Detected_WT"]),
                         "legacy_hits": int(r["Hit_Count_WT"])})
    lg = pd.DataFrame(rows)

    new = pd.read_csv(TAB / "fig5a_mean_rank_by_group_tool.tsv", sep="\t")
    new = new[new["metric"].isin(METRIC_MAP.values())]
    new = new[["species", "tool", "metric", "n_units", "metric_mean", "metric_sd",
               "metric_min", "metric_max"]].rename(
        columns={"metric_mean": "revision_mean", "metric_sd": "revision_sd",
                 "metric_min": "revision_min", "metric_max": "revision_max"})

    cmp = lg.merge(new, on=["species", "tool", "metric"], how="outer")
    cmp["delta_revision_minus_legacy"] = cmp["revision_mean"] - cmp["legacy_value"]
    cmp["ratio_revision_over_legacy"] = cmp["revision_mean"] / cmp["legacy_value"]
    cmp["legacy_source"] = str(LEGACY_CSV.relative_to(C.PROJECT))
    cmp["revision_source"] = ("harmonisation/evaluation/tables/m6a_glori_confusion.tsv "
                              f"(window={C.PRIMARY_WINDOW}, min_cov={C.C_MIN_DEFAULT}, "
                              "explicit universe, per independent unit)")
    cmp = cmp.sort_values(["species", "tool", "metric"])
    write_table(cmp, TAB / "fig5_legacy_vs_revision.tsv")

    unmatched = cmp[cmp["legacy_value"].isna() | cmp["revision_mean"].isna()]
    if len(unmatched):
        logger.warning("rows without both sides: %s",
                       "; ".join(f"{r.species}/{r.tool}/{r.metric}"
                                 for r in unmatched.itertuples()))
    logger.info("legacy-vs-revision rows: %d (tools matched: %d)",
                len(cmp), cmp["tool"].nunique())
    for metric in ("precision", "recall", "f1", "mcc"):
        s = cmp[cmp["metric"] == metric].dropna(subset=["legacy_value", "revision_mean"])
        if len(s):
            logger.info("%-9s legacy mean=%.4f  revision mean=%.4f  ratio=%.2f",
                        metric, s["legacy_value"].mean(), s["revision_mean"].mean(),
                        (s["revision_mean"].mean() / s["legacy_value"].mean()))
    logger.info("output: %s", TAB / "fig5_legacy_vs_revision.tsv")


if __name__ == "__main__":
    main()
