#!/usr/bin/env python
"""NA-4b -- threshold-sensitivity curves for representative m6A tools of each class.

The response letter (R1-5 / E4) states that "threshold-sensitivity curves for
representative tools of each class show that the main conclusions are robust
across a range of calling thresholds".  This script produces that evidence: for
m6Anet (probabilistic model), DENA (random forest / ratio), Nanom6A (XGBoost,
mod-ratio cutoff), ELIGOS2 (both the solo and the diff variant) and Nanocompore
(differential test, GMM p-value) the per-unit confusion matrix against GLORI is
recomputed on the manuscript (w = 2 bp, cov >= 10 universe, one row per
independent sequencing unit) as the calling threshold is swept.

Threshold semantics (from the callset ``score`` column, never re-derived):

=============  =============  ============================================
tool           score_type     rule
=============  =============  ============================================
m6Anet         Prob           keep ``score >= t``   (default 0.5)
DENA           m6a_ratio      keep ``score >= t``   (manuscript cutoff 0.1)
Nanom6A        mod_ratio      keep ``score >= t``   (manuscript cutoff 0.5)
Nanocompore    Pvalue         keep ``score <= t``   (GMM p-value, default 0.05)
ELIGOS2_solo   adjPval/Padj   keep ``score <= t``   (FDR-adjusted p-value)
ELIGOS2_diff   adjPval/Padj   keep ``score <= t``
=============  =============  ============================================

The score endpoints of every callset are tabulated as well: for m6Anet the file
starts exactly at 0.5 and for DENA at 0.1, i.e. the distributed callsets are
already cut at the documented threshold and the sweep can only be tightened
there; the Nanom6A file spans 0.1-1.0, so its sweep brackets the documented 0.5
cutoff; the ELIGOS2 files contain adjusted p-values capped at 1e-4, and the
Nanocompore export is filtered at the tool default p < 0.05 (with |logit LOR| >
0.5) upstream, so its sweep tightens from 0.05 as well.

2026-10-10 (R1-5 shipping round): Nanocompore is added as the second
differential-testing representative, the five variants become six, and the
``no_filter`` baseline still anchors every curve to the published callset.

Outputs
-------
tables/figS11_sweep_per_unit.tsv        unit x tool x threshold metrics
tables/figS11_sweep_summary.tsv         species x tool x threshold (mean over units)
tables/figS11_default_summary.tsv       the documented operating point of each tool
tables/figS11_score_endpoints.tsv       score distribution of every callset
tables/figS11_tool_rank_by_fraction.tsv PPV / F1 ordering of the six variants
                                      across retained-call fractions

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/figures/figureS11/src/figS11_sweep.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
R21_ANALYSIS = HERE.parents[2] / "analysis/coverage_rank_stability/analysis"
sys.path.insert(0, str(R21_ANALYSIS))
sys.path.insert(0, str(HERE))

from r21_common import (C, METRICS, PANELS, POOLED, WINDOW,                     # noqa: E402
                        callset_path, calls_by_stratum, load_banded_universe,
                        load_reference, read_table, setup_logger, strata_confusion,
                        universe_file, write_table)

OUT = HERE.parent
TAB = OUT / "tables"
LOG = OUT / "logs"

#: tool -> (score column rule, grid, documented operating point, class)
TOOLS: dict[str, dict] = {
    "m6Anet": {"rule": "ge", "grid": [0.5, 0.6, 0.7, 0.8, 0.9, 0.95],
               "default": 0.5, "class": "probabilistic model"},
    "DENA": {"rule": "ge", "grid": [0.1, 0.15, 0.2, 0.25, 0.3, 0.4],
             "default": 0.1, "class": "random forest"},
    "Nanom6A": {"rule": "ge", "grid": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9],
                "default": 0.5, "class": "XGBoost"},
    #: 2026-10-10 (R1-5 shipping round): the exported callset is already filtered at the
    #: tool default p < 0.05 (plus |logit LOR| > 0.5), so the grid tightens from there and
    #: 0.05 is the documented operating point the Figure S11 open marker sits on
    "Nanocompore": {"rule": "le", "grid": [0.05, 0.01, 1e-3, 1e-4, 1e-6, 1e-8],
                    "default": 0.05, "class": "differential test"},
    "ELIGOS2_solo": {"rule": "le", "grid": [1e-4, 1e-5, 1e-6, 1e-8, 1e-12, 1e-20],
                     "default": 1e-4, "class": "differential test"},
    "ELIGOS2_diff": {"rule": "le", "grid": [1e-4, 1e-5, 1e-6, 1e-8, 1e-12, 1e-20],
                     "default": 1e-4, "class": "differential test"},
}
NONE_LABEL = "no_filter"


def drawn_units(logger=None) -> list:
    """The Fig. 5A drawn units (Arabidopsis 3, HeLa 3, mouse mES_WT)."""
    out = []
    for species, group, units in PANELS:
        for s in C.SAMPLES:
            if (s.platform == C.PLATFORM_RNA002 and s.species == species
                    and s.dataset_group == group and s.canonical in units):
                out.append(s)
    if logger:
        logger.info("drawn units: %s", ", ".join(s.canonical for s in out))
    return out


def score_endpoints(df: pd.DataFrame) -> dict:
    s = pd.to_numeric(df["score"], errors="coerce").dropna()
    if not s.size:
        return {}
    q = np.percentile(s, [0, 1, 5, 25, 50, 75, 95, 100])
    return {"n_rows": int(len(df)), "n_scored": int(s.size),
            "score_min": float(q[0]), "score_p1": float(q[1]), "score_p5": float(q[2]),
            "score_p25": float(q[3]), "score_median": float(q[4]),
            "score_p75": float(q[5]), "score_p95": float(q[6]), "score_max": float(q[7])}


def main() -> None:
    logger = setup_logger("figS11_sweep", log_dir=LOG)
    TAB.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rows: list[dict] = []
    eps: list[dict] = []

    for sample in drawn_units(logger):
        uni = load_banded_universe(universe_file(sample))
        ref = load_reference(C.GLORI[sample.species])
        logger.info("[%s] universe pooled=%d", sample.canonical, uni.n_pooled)
        for tool, spec in TOOLS.items():
            path = callset_path(sample, tool)
            if not path.exists():
                logger.warning("[%s] callset missing: %s", sample.canonical, tool)
                continue
            df = read_table(path)
            if "score" not in getattr(df, "columns", []):
                logger.warning("[%s] %s: no score column, skipped",
                               sample.canonical, tool)
                continue
            df = df.copy()
            df["pos_raw"] = pd.to_numeric(df["pos_raw"]).astype(np.int64)
            df["chrom"] = df["chrom"].astype(str)
            score_type = str(df.get("score_type", pd.Series(["?"])).iloc[0])
            ep = score_endpoints(df)
            eps.append({"species": sample.species, "dataset_group": sample.dataset_group,
                        "sample": sample.canonical, "tool": tool,
                        "score_type": score_type, "rule": spec["rule"],
                        "class": spec["class"],
                        "documented_default": spec["default"], **ep})

            # baseline = full callset (the manuscript operating point), then the sweep
            points: list[tuple[str, float | None]] = [(NONE_LABEL, None)]
            lo = ep.get("score_min", np.nan)
            hi = ep.get("score_max", np.nan)
            for t in spec["grid"]:
                if spec["rule"] == "ge":
                    keep_any = np.isfinite(hi) and hi >= t
                else:
                    keep_any = np.isfinite(lo) and lo <= t
                if keep_any:
                    points.append((f"{t:g}", float(t)))

            n_pos_full = None
            for label, t in points:
                sub = df if t is None else (
                    df[pd.to_numeric(df["score"], errors="coerce") >= t]
                    if spec["rule"] == "ge"
                    else df[pd.to_numeric(df["score"], errors="coerce") <= t])
                calls, n_out, n_in = calls_by_stratum(sub, uni, [POOLED])
                conf = strata_confusion(uni, ref, calls, [POOLED])[POOLED]
                if n_pos_full is None:
                    n_pos_full = conf["n_calls_in_universe"]
                mcc = conf.get("mcc", np.nan)
                rows.append({
                    "species": sample.species, "dataset_group": sample.dataset_group,
                    "sample": sample.canonical, "tool": tool, "class": spec["class"],
                    "score_type": score_type, "rule": spec["rule"],
                    "threshold": label, "threshold_value": (np.nan if t is None else t),
                    "is_default": bool(t is not None and abs(t - spec["default"]) < 1e-12),
                    "is_full_callset": bool(t is None),
                    "n_rows": int(len(sub)),
                    "n_calls_in_universe": conf["n_calls_in_universe"],
                    "n_calls_out_of_universe": int(n_out),
                    "n_reference_in_universe": conf["n_reference_in_universe"],
                    "n_glori": conf["n_glori"], "universe": conf["universe"],
                    "retained_fraction": (conf["n_calls_in_universe"] / n_pos_full
                                          if n_pos_full else np.nan),
                    "tp": conf["tp"], "fp": conf["fp"], "fn": conf["fn"],
                    "tn": conf["tn"], **{m: conf[m] for m in METRICS},
                })
            logger.info("[%s] %-13s score %s %.4g..%.4g  points=%d",
                        sample.canonical, tool, score_type,
                        ep.get("score_min", np.nan), ep.get("score_max", np.nan),
                        len(points))

    per_unit = pd.DataFrame(rows)
    write_table(per_unit, TAB / "figS11_sweep_per_unit.tsv")
    write_table(pd.DataFrame(eps), TAB / "figS11_score_endpoints.tsv")

    # ---- summary over units ------------------------------------------------- #
    agg = (per_unit.groupby(["species", "dataset_group", "tool", "class", "score_type",
                             "rule", "threshold", "threshold_value", "is_default"],
                            dropna=False, as_index=False)
           .agg(n_units=("sample", "nunique"),
                retained_fraction_mean=("retained_fraction", "mean"),
                n_calls_mean=("n_calls_in_universe", "mean"),
                **{f"{m}_mean": (m, "mean") for m in METRICS},
                **{f"{m}_min": (m, "min") for m in METRICS},
                **{f"{m}_max": (m, "max") for m in METRICS}))
    write_table(agg, TAB / "figS11_sweep_summary.tsv")
    write_table(agg[agg["is_default"]], TAB / "figS11_default_summary.tsv")

    # ---- ordering of the five variants across retained fractions ------------ #
    frac_rows = []
    for (species, group), sub in per_unit.groupby(["species", "dataset_group"], sort=False):
        for tool in TOOLS:
            t = sub[sub["tool"] == tool]
            if t.empty:
                continue
            for m in ["precision", "f1"]:
                # per retained-fraction bucket: mean metric of this tool
                for name, lo, hi in [("full", -0.01, 1e-9), ("f>=0.75", 0.75, 1.01),
                                     ("0.5-0.75", 0.5, 0.75), ("0.25-0.5", 0.25, 0.5),
                                     ("<0.25", -0.01, 0.25)]:
                    if name == "full":
                        sel = t[t["is_full_callset"]]
                    else:
                        sel = t[(~t["is_full_callset"])
                                & (t["retained_fraction"] >= lo)
                                & (t["retained_fraction"] < hi)]
                    if sel.empty:
                        continue
                    frac_rows.append({"species": species, "dataset_group": group,
                                      "bucket": name, "tool": tool, "metric": m,
                                      "value": float(pd.to_numeric(sel[m]).mean())})
    fr = pd.DataFrame(frac_rows)
    if not fr.empty:
        fr["rank"] = fr.groupby(["species", "bucket", "metric"])["value"].rank(
            method="min", ascending=False)
    write_table(fr, TAB / "figS11_tool_rank_by_fraction.tsv")

    logger.info("per-unit rows=%d  tools=%d  score endpoints=%d  (%.1f min)",
                len(per_unit), per_unit["tool"].nunique(), len(eps),
                (time.time() - t0) / 60)
    # the two extremes of every curve (full callset -> strictest grid point)
    for (species, tool), sub in per_unit.groupby(["species", "tool"], sort=False):
        sub = sub.sort_values("retained_fraction", ascending=False)
        first, last = sub.iloc[0], sub.iloc[-1]
        logger.info("%-12s %-13s full: PPV=%.3f F1=%.3f (retained %.2f) -> "
                    "strictest t=%s: PPV=%.3f F1=%.3f (retained %.2f)",
                    species, tool, first["precision"], first["f1"],
                    first["retained_fraction"], last["threshold"],
                    last["precision"], last["f1"], last["retained_fraction"])
    logger.info("output: figS11_sweep_per_unit.tsv, figS11_sweep_summary.tsv, "
                "figS11_default_summary.tsv, figS11_score_endpoints.tsv, "
                "figS11_tool_rank_by_fraction.tsv")
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
