#!/usr/bin/env python3
"""Export a figure-ready per-replicate table (R3-2 / E6 error-bar decisions).

Reads the evaluation tables and writes one row per (dataset_group, tool, window)
carrying what a plot actually needs: the individual replicate values, the mean,
the between-replicate SD, the site-level bootstrap CI of each replicate, the
cross-replicate Jaccard, and an explicit verdict on whether an error bar may be
drawn at all.

The bar rule is encoded from the independence model in ``common/config.py``:

  mean +/- SD + the individual points   replicate_class == "replicate" and n >= 2
  two separate points, NO bar           replicate_class == "cross_study"
  single value, NO bar                  n == 1

``precision`` here is the manuscript's "hit rate" = TP/(TP+FP), and the two error
sources stay separate on purpose: ``*_ci_lo/hi`` is sequence-sampling error inside
one replicate, ``*_sd`` is biological/run-to-run spread across replicates.  They
answer different reviewer questions and must not be merged into one bar.

Usage
-----
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/export_figure_ready.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

CODE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE))

from common.config import (PRIMARY_WINDOW, TABLE_DIR, WINDOWS, independent_units,
                           replicate_class, SAMPLES)
from common.manifest import setup_logger

GROUP_OF = {s.dataset_group for s in SAMPLES}


def members(group: str) -> list[str]:
    return [s.canonical for s in SAMPLES if s.dataset_group == group]


def independent_n(group: str) -> int:
    ms = [s for s in SAMPLES if s.dataset_group == group]
    return len(independent_units(ms))


def main() -> None:
    logger = setup_logger("export_figure_ready")
    conf = pd.read_csv(TABLE_DIR / "m6a_glori_confusion.tsv", sep="\t", engine="python")
    repro = pd.read_csv(TABLE_DIR / "reproducibility.tsv", sep="\t", engine="python")
    anchor = pd.read_csv(TABLE_DIR / "legacy_published_check.tsv", sep="\t", engine="python")

    #: per-window species anchor recomputed by NA-1 (NOT the manuscript's
    #: single-sample number -- that one only exists for one window, see
    #: PUBLISHED_LEGACY below).  "Human (HeLa)" -> "Human".
    pub = {(str(r.species).split(" (")[0], int(r.window)): float(r.published_or_NA1)
           for _, r in anchor.iterrows()}
    #: the values reviewer 3 quotes (R3-1): Arabidopsis 30.33 %, mouse 22.67 %,
    #: HeLa 22.09 % -- legacy, one sample per species, mean over tools.
    PUBLISHED_LEGACY = {"Arabidopsis": 0.303294, "Mouse": 0.226700, "Human": 0.220900}
    summ = repro[repro.support == "group_summary"].set_index(["dataset_group", "tool"])

    rows = []
    for w in (0, PRIMARY_WINDOW):
        c = conf[conf.window == w].copy()
        for (group, tool), g in c.groupby(["dataset_group", "tool"]):
            g = g.sort_values("sample")
            klass = replicate_class(group)
            n_indep = independent_n(group)
            pre = g.precision.to_numpy(dtype=float)
            rec = g.recall.to_numpy(dtype=float)
            nsites = g.n_calls_in_universe.to_numpy(dtype=float)
            key = (group, tool)
            jac = mean_jac = np.nan
            if key in summ.index:
                s = summ.loc[key]
                jac = float(s.mean_pairwise_jaccard) if pd.notna(s.mean_pairwise_jaccard) else np.nan
                mean_jac = float(s.min_pairwise_jaccard) if pd.notna(s.min_pairwise_jaccard) else np.nan
            if klass == "cross_study":
                bar, why = "none", "different studies -> concordance, not replication"
            elif n_indep < 2:
                bar, why = "none", "n = 1 independent sequencing unit"
            else:
                bar, why = "mean±SD + points", f"{n_indep} independent replicates"
            rows.append({
                "platform": g.platform.iloc[0], "species": g.species.iloc[0],
                "dataset_group": group, "tool": tool, "window": w, "min_cov": g.min_cov.iloc[0],
                "replicates": "|".join(g["sample"]), "n_rows_in_table": len(g),
                "n_independent": n_indep, "replicate_class": klass,
                "error_bar": bar, "error_bar_reason": why,
                "precision_each": "|".join(f"{v:.4f}" for v in pre),
                "precision_mean": float(np.mean(pre)),
                "precision_sd": float(np.std(pre, ddof=1)) if len(pre) > 1 else np.nan,
                "precision_min": float(np.min(pre)), "precision_max": float(np.max(pre)),
                "precision_site_ci_lo": "|".join(f"{v:.4f}" for v in g.precision_ci_lo),
                "precision_site_ci_hi": "|".join(f"{v:.4f}" for v in g.precision_ci_hi),
                "recall_each": "|".join(f"{v:.4f}" for v in rec),
                "recall_mean": float(np.mean(rec)),
                "recall_sd": float(np.std(rec, ddof=1)) if len(rec) > 1 else np.nan,
                "sites_called_each": "|".join(f"{int(v)}" for v in nsites),
                "sites_called_mean": float(np.mean(nsites)),
                "tp_sum": int(g.tp.sum()), "fp_sum": int(g.fp.sum()), "fn_sum": int(g.fn.sum()),
                "mean_pairwise_jaccard": jac, "min_pairwise_jaccard": mean_jac,
                "na1_species_anchor_this_window": pub.get((g.species.iloc[0], w), np.nan),
                "manuscript_legacy_hit_rate": PUBLISHED_LEGACY.get(g.species.iloc[0], np.nan),
            })

    out = pd.DataFrame(rows).sort_values(["window", "species", "dataset_group", "tool"])
    path = TABLE_DIR / "figure_ready_replicates.tsv"
    out.to_csv(path, sep="\t", index=False)

    logger.info("rows=%d -> %s", len(out), path)
    for grp in sorted(out.dataset_group.unique()):
        s = out[(out.dataset_group == grp) & (out.window == PRIMARY_WINDOW)]
        if s.empty:
            continue
        logger.info("%-18s %-22s n_indep=%d =%d", grp,
                    s.error_bar.iloc[0], s.n_independent.iloc[0], len(s))
    logger.info("precision_site_ci_* precision_sd ")


if __name__ == "__main__":
    main()
