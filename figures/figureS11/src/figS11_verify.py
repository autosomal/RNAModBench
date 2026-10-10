#!/usr/bin/env python
"""NA-4b gate -- the curves must start from the manuscript's own numbers.

Checks (all hard unless marked *info*):

1. anchor: every ``no_filter`` row (the full callset = the published operating
   point) equals ``harmonisation/evaluation/tables/m6a_glori_confusion.tsv`` at
   w = 2 / min_cov = 10 -- all unit x tool rows, counts and metrics, max_diff = 0;
2. monotone tightening: within each unit x tool, ``n_calls_in_universe`` never
   increases as the threshold becomes stricter;
3. ``retained_fraction`` equals ``n_calls_in_universe / full-callset calls`` and
   lies in (0, 1];
4. contingency integrity of every row (tp+fp == calls, tn == U-tp-fp-fn,
   fn <= n_glori);
5. exactly one row per unit x tool carries ``is_default``;
6. *info*: is the distributed callset already cut at the documented threshold?
   (score_min >= default for ``ge`` tools, score_max <= default for ``le`` tools)

Exit code != 0 if any hard check fails.

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/figures/figureS11/src/figS11_verify.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "analysis/coverage_rank_stability/analysis"))
sys.path.insert(0, str(HERE))

from figS11_sweep import TOOLS                                              # noqa: E402
from r21_common import C, METRICS, setup_logger, write_table              # noqa: E402

OUT = HERE.parent
PER_UNIT = OUT / "tables" / "figS11_sweep_per_unit.tsv"
ENDPOINTS = OUT / "tables" / "figS11_score_endpoints.tsv"
SANITY = OUT / "tables" / "figS11_sanity_check.tsv"
CONFUSION = C.TABLE_DIR / "m6a_glori_confusion.tsv"
KEY = ["species", "dataset_group", "sample", "tool"]
COUNT_COLS = ["tp", "fp", "fn", "tn", "universe", "n_calls_in_universe",
              "n_reference_in_universe"]


def _row(check, scope, n_total, n_pass, max_diff, status, note=""):
    return {"check": check, "scope": scope, "n_total": n_total, "n_pass": n_pass,
            "max_abs_diff": max_diff, "status": status, "note": note}


def check_anchor(mine: pd.DataFrame) -> dict:
    conf = pd.read_csv(CONFUSION, sep="\t")
    conf = conf[(conf["window"] == 2) & (conf["min_cov"] == 10)]
    full = mine[mine["is_full_callset"]]
    cols = COUNT_COLS + METRICS
    both = full[KEY + cols].merge(conf[KEY + cols], on=KEY, suffixes=("_ours", "_conf"))
    max_diff = 0.0
    n_pass = 0
    for _, r in both.iterrows():
        d = max(abs(float(r[f"{c}_ours"]) - float(r[f"{c}_conf"])) for c in cols
                if pd.notna(r[f"{c}_ours"]) and pd.notna(r[f"{c}_conf"]))
        max_diff = max(max_diff, d)
        n_pass += int(d < 1e-9)
    missing = len(full) - len(both)
    return _row("anchor_full_callset_vs_confusion", "unit x tool", len(both), n_pass,
                max_diff, "pass" if n_pass == len(both) and missing == 0 else "FAIL",
                f"{CONFUSION.name} (w=2, min_cov=10); rows missing: {missing} "
                "(the full callset IS the published operating point)")


def check_monotone(mine: pd.DataFrame) -> dict:
    n_total = n_pass = 0
    max_diff = 0.0
    for (_s, tool), g in mine.groupby(["sample", "tool"], sort=False):
        rule = TOOLS[tool]["rule"]
        g = g[~g["is_full_callset"]].copy()
        v = pd.to_numeric(g["threshold_value"], errors="coerce").to_numpy(dtype=float)
        c = pd.to_numeric(g["n_calls_in_universe"], errors="coerce").to_numpy(dtype=float)
        order = np.argsort(v if rule == "ge" else -v)          # loose -> strict
        d = np.diff(c[order])
        n_total += 1
        bad = float(max(d.max(), 0)) if d.size else 0.0
        max_diff = max(max_diff, bad)
        n_pass += int(bad == 0)
    return _row("threshold_tightening_monotone", "unit x tool", n_total, n_pass,
                max_diff, "pass" if n_pass == n_total else "FAIL",
                "a stricter threshold can never retain more calls")


#: tables are written with ``float_format="%.6g"`` -> 6 significant digits, so the
#: stored fraction can differ from the recomputed one by ~5e-7
FRACTION_TOL = 5e-6


def check_fractions(mine: pd.DataFrame) -> dict:
    full = mine[mine["is_full_callset"]].set_index(["sample", "tool"])["n_calls_in_universe"]
    n_total = n_pass = 0
    max_diff = 0.0
    for _, r in mine.iterrows():
        base = full.get((r["sample"], r["tool"]), np.nan)
        exp = r["n_calls_in_universe"] / base if base else np.nan
        got = r["retained_fraction"]
        n_total += 1
        d = abs(float(got) - float(exp)) if np.isfinite(exp) and np.isfinite(got) else 9.9
        max_diff = max(max_diff, d if d < 9 else 0.0)
        ok_in_range = np.isfinite(got) and 0 <= got <= 1 + 1e-12
        n_pass += int(d < FRACTION_TOL and ok_in_range)
    return _row("retained_fraction_consistent", "every row", n_total, n_pass,
                max_diff, "pass" if n_pass == n_total else "FAIL",
                f"retained_fraction == calls / full-callset calls and lies in [0, 1] "
                f"(tolerance {FRACTION_TOL:g} for the %.6g table format)")


def check_integrity(mine: pd.DataFrame) -> dict:
    tp = pd.to_numeric(mine["tp"], errors="coerce").astype(float)
    fp = pd.to_numeric(mine["fp"], errors="coerce").astype(float)
    fn = pd.to_numeric(mine["fn"], errors="coerce").astype(float)
    calls = pd.to_numeric(mine["n_calls_in_universe"], errors="coerce").astype(float)
    glori = pd.to_numeric(mine["n_glori"], errors="coerce").astype(float)
    ok = ((tp + fp - calls).abs() < 1e-9) & (fn <= glori) & (tp <= calls)
    ok &= ((pd.to_numeric(mine["tn"])
            - (pd.to_numeric(mine["universe"]) - tp - fp - fn)).abs() < 1e-9)
    max_diff = float(max((tp + fp - calls).abs().max(),
                         (pd.to_numeric(mine["tn"])
                          - (pd.to_numeric(mine["universe"]) - tp - fp - fn)).abs().max()))
    return _row("contingency_integrity", "every row", len(mine), int(ok.sum()),
                max_diff, "pass" if ok.all() else "FAIL",
                "tp+fp == calls; tn == U-tp-fp-fn; fn <= n_glori; tp <= calls")


def check_default_flag(mine: pd.DataFrame) -> dict:
    g = mine.groupby(["sample", "tool"])["is_default"].sum()
    ok = g == 1
    return _row("documented_default_present", "unit x tool", len(g), int(ok.sum()),
                0.0, "pass" if ok.all() else "FAIL",
                "exactly one grid point per unit x tool equals the documented "
                "operating point: " + ", ".join(
                    f"{t}={TOOLS[t]['default']:g}" for t in TOOLS))


def check_precut(endpoints: pd.DataFrame) -> dict:
    """info: does the distributed callset already sit at the documented cut?"""
    rows = []
    for _, r in endpoints.iterrows():
        spec = TOOLS[r["tool"]]
        if spec["rule"] == "ge":
            pre = float(r["score_min"]) >= float(spec["default"]) - 1e-9
        else:
            pre = float(r["score_max"]) <= float(spec["default"]) * 1.000001
        rows.append(pre)
    n_pre = int(np.sum(rows))
    detail = "; ".join(
        f"{r['tool']} on {r['sample']}: {r['score_min']:.3g}-{r['score_max']:.3g} "
        f"(default {TOOLS[r['tool']]['default']:g})"
        for _, r in endpoints.iterrows())
    return _row("callset_precut_at_default", "unit x tool", len(rows), n_pre, 0.0,
                "info", "callsets already cut at the documented threshold (so the "
                f"sweep can only tighten): {n_pre}/{len(rows)} | {detail}")


def main() -> None:
    logger = setup_logger("figS11_verify", log_dir=OUT / "logs")
    mine = pd.read_csv(PER_UNIT, sep="\t")
    endpoints = pd.read_csv(ENDPOINTS, sep="\t")
    checks = [check_anchor(mine), check_monotone(mine), check_fractions(mine),
              check_integrity(mine), check_default_flag(mine), check_precut(endpoints)]
    out = pd.DataFrame(checks)
    write_table(out, SANITY)
    for _, r in out.iterrows():
        logger.info("%-34s %-6s n=%d/%d max_diff=%s  %s", r["check"], r["status"],
                    r["n_pass"], r["n_total"],
                    f"{r['max_abs_diff']:.3g}" if pd.notna(r["max_abs_diff"]) else "-",
                    r["note"][:150])
    if (out["status"] == "FAIL").any():
        logger.error("HARD CHECKS FAILED: %s",
                     ", ".join(out[out["status"] == "FAIL"]["check"]))
        raise SystemExit(1)
    logger.info("all hard checks passed; sanity table -> %s", SANITY)
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
