#!/usr/bin/env python
"""Figure 5B, revision rebuild (R1-2 / R2-1 / R3-2).

The published panel B (``mod_ratio_hit_rate.ipynb``) plotted, for four tools
that report a stoichiometry (MINES, m6Anet, DENA, Nanom6A), the GLORI hit rate
inside ten bins of the tool's own modification ratio.  It read the legacy
``output/<group>/...`` aggregates and matched with ``contains_glori`` on a
closed interval, i.e. one legacy sample per species (no replicate structure).

This script recomputes the same quantity from the reconstructable per-replicate
call sets:

  ratio      the tool's per-site modification ratio: ``mod_ratio`` when the
             callset carries it, otherwise ``score`` when ``score_type`` denotes
             a ratio (same rule as ``15_mod_ratio_replicate_agreement.py``;
             values clipped to [0, 1]);
  hit        ``dist_glori <= 2`` -- the primary matching window of the revision
             (the published panel used an interval overlap, which is the w = 0
             equivalent and only credits an exact anchor).  The y axis is
             labelled ``PPV vs. GLORI (2 bp)``: the quantity is precision
             (positive predictive value) against the reference, not "hit rate";
  bins       [0, 0.1, ..., 1.0], labels ``0.0-0.1`` ... identical to the panel's;
  replicates per bin, then mean across the independent units of the group, with
             the unit range drawn as a band -- never a merged call set.

Replicate structure and the drawn units are identical to Figure 5A
(``35_fig5a_metric_ranks.py``): Arabidopsis 3 biological replicates, HeLa 3,
Mouse the single study ``mES_WT``.

Outputs -> 04_revision_analysis/fig5_revision/{tables,figures}

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/src/sites_v2/scripts/36_fig5b_modratio_hitrate.py
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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from common import config as C                                        # noqa: E402
from common.figstyle import apply as apply_style                     # noqa: E402
from common.figstyle import save                                     # noqa: E402
from common.io_utils import write_table                              # noqa: E402
from common.manifest import setup_logger                             # noqa: E402

OUT = (_RB / "figures/figure5")
TAB, FIG = OUT / "tables", OUT / "figures"

TOOLS = ["MINES", "m6Anet", "DENA", "Nanom6A"]
TOOL_COLOR = {"MINES": "#4C72B0", "m6Anet": "#DD8452",
              "DENA": "#55A868", "Nanom6A": "#C44E52"}

BINS = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
BIN_LABELS = [f"{BINS[i]:.1f}-{BINS[i + 1]:.1f}" for i in range(len(BINS) - 1)]

WINDOW = C.PRIMARY_WINDOW
MIN_SITES_PER_BIN = 5          # below this a bin is not drawn for that unit
RATIO_SCORE_TYPES = {"mod_ratio", "m6a_ratio", "ratio", "percent_modified/100"}

PANELS: list[tuple[str, str, list[str]]] = [
    ("Arabidopsis", "Arabidopsis_WT",
     ["Arabidopsis_WT_rep1", "Arabidopsis_WT_rep2", "Arabidopsis_WT_rep3"]),
    ("Mouse", "Mouse_WT", ["mES_WT"]),
    ("Human", "HeLa_WT", ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"]),
]


# --------------------------------------------------------------------------- #
def callset_path(species: str, group: str, tool: str, sample: str) -> Path:
    return (C.CALLSET_ROOT / C.PLATFORM_RNA002 / species / group / "m6A"
            / tool / f"{sample}.tsv")


def extract_ratio(df: pd.DataFrame, logger) -> pd.Series:
    """Per-site modification ratio: ``mod_ratio`` if populated, else ratio-score."""
    mr = (pd.to_numeric(df["mod_ratio"], errors="coerce")
          if "mod_ratio" in df.columns else pd.Series(np.nan, index=df.index))
    if mr.notna().any():
        return mr.clip(0, 1)
    st = ""
    if "score_type" in df.columns:
        seen = [str(v).strip().lower() for v in df["score_type"].dropna().unique() if str(v).strip()]
        st = seen[0] if seen else ""
    if st in RATIO_SCORE_TYPES or "ratio" in st:
        return pd.to_numeric(df["score"], errors="coerce").clip(0, 1)
    logger.warning("no ratio column (score_type=%r) in %s rows", st, len(df))
    return pd.Series(np.nan, index=df.index)


def per_unit_bins(logger) -> pd.DataFrame:
    """One row per (species, group, sample, tool, bin) with sites/hits/rate."""
    rows = []
    for species, group, units in PANELS:
        for sample in units:
            for tool in TOOLS:
                path = callset_path(species, group, tool, sample)
                if not path.exists():
                    logger.warning("callset missing: %s", path)
                    continue
                df = pd.read_csv(path, sep="\t", low_memory=False)
                if df.empty:
                    continue
                ratio = extract_ratio(df, logger)
                dist = pd.to_numeric(df.get("dist_glori"), errors="coerce")
                d = pd.DataFrame({"ratio": ratio, "dist": dist}).dropna()
                if d.empty:
                    logger.warning("%s/%s/%s: no ratio+dist rows", species, tool, sample)
                    continue
                d = d.assign(bin=pd.cut(d["ratio"], bins=BINS,
                                        labels=BIN_LABELS, include_lowest=True),
                             hit=d["dist"] <= WINDOW)
                for b in BIN_LABELS:
                    sub = d[d["bin"] == b]
                    n = int(sub.shape[0])
                    if n == 0:
                        continue
                    rows.append({
                        "species": species, "dataset_group": group,
                        "sample": sample, "tool": tool, "mod_ratio_bin": b,
                        "n_sites": n, "n_hits": int(sub["hit"].sum()),
                        "hit_rate": float(sub["hit"].mean()),
                        "usable": n >= MIN_SITES_PER_BIN,
                    })
    return pd.DataFrame(rows)


def summarise(per_unit: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for species, group, units in PANELS:
        for tool in TOOLS:
            for b in BIN_LABELS:
                sub = per_unit[(per_unit["species"] == species)
                               & (per_unit["tool"] == tool)
                               & (per_unit["mod_ratio_bin"] == b)]
                used = sub[sub["usable"]]
                rows.append({
                    "species": species, "dataset_group": group, "tool": tool,
                    "mod_ratio_bin": b,
                    "n_units_total": len(units),
                    "n_units_used": int(used["sample"].nunique()),
                    "n_sites_per_unit": "|".join(str(v) for v in sub["n_sites"]),
                    "n_sites_mean": float(sub["n_sites"].mean()) if len(sub) else np.nan,
                    "hit_rate_mean": float(used["hit_rate"].mean()) if len(used) else np.nan,
                    "hit_rate_min": float(used["hit_rate"].min()) if len(used) else np.nan,
                    "hit_rate_max": float(used["hit_rate"].max()) if len(used) else np.nan,
                    "hit_rate_each": "|".join(f"{v:.4f}" for v in used["hit_rate"]),
                })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
def figure(summary: pd.DataFrame) -> None:
    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.4), constrained_layout=True)
    x = np.arange(len(BIN_LABELS))

    for ax, (species, group, units) in zip(axes, PANELS):
        for tool in TOOLS:
            s = (summary[(summary["species"] == species) & (summary["tool"] == tool)]
                 .set_index("mod_ratio_bin").reindex(BIN_LABELS))
            y = s["hit_rate_mean"].to_numpy(dtype=float)
            lo = s["hit_rate_min"].to_numpy(dtype=float)
            hi = s["hit_rate_max"].to_numpy(dtype=float)
            ok = np.isfinite(y)
            if not ok.any():
                continue
            col = TOOL_COLOR[tool]
            if np.isfinite(lo).sum() > 1:      # unit range, no error bars
                ax.fill_between(x[ok], lo[ok], hi[ok], color=col, alpha=.15, lw=0)
            ax.plot(x[ok], y[ok], "-o", color=col, ms=5, lw=1.8, label=tool)
        ax.set_xticks(x)
        ax.set_xticklabels(BIN_LABELS, rotation=45, ha="right", fontsize=8.5)
        n = len(units)
        ax.set_title(f"{species} ({'n = 1 study' if n == 1 else f'n = {n}'})",
                     fontweight="bold", fontsize=15, pad=10)
        ax.set_xlabel("Modification ratio reported by the tool")
        ax.set_ylim(0, 1)
        ax.set_xlim(-0.4, len(BIN_LABELS) - 0.6)
    # naming (user decision 2026-09-19): the quantity is precision/PPV against the
    # GLORI reference, so the axis says PPV rather than the old "GLORI hit rate"
    axes[0].set_ylabel(f"PPV vs. GLORI ({WINDOW} bp)")
    axes[0].legend(frameon=False, fontsize=10, loc="upper left")

    FIG.mkdir(parents=True, exist_ok=True)
    stem = FIG / "Fig5B_modratio_hitrate"
    save(fig, str(stem))
    print("wrote", stem.with_suffix(".pdf").name, "+ .png")


# --------------------------------------------------------------------------- #
def main() -> None:
    logger = setup_logger("36_fig5b_modratio_hitrate")
    TAB.mkdir(parents=True, exist_ok=True)

    per_unit = per_unit_bins(logger)
    write_table(per_unit, TAB / "fig5b_modratio_bins_per_replicate.tsv")
    summary = summarise(per_unit)
    write_table(summary, TAB / "fig5b_bin_summary.tsv")

    for species, group, units in PANELS:
        for tool in TOOLS:
            s = summary[(summary["species"] == species) & (summary["tool"] == tool)]
            drawn = int(np.isfinite(s["hit_rate_mean"]).sum())
            logger.info("%-12s %-8s bins drawn=%2d/%d  n_units_used=%s",
                        species, tool, drawn, len(BIN_LABELS),
                        sorted(s.loc[np.isfinite(s["hit_rate_mean"]), "n_units_used"]
                               .unique().tolist()))
    figure(summary)
    logger.info("output: %s", OUT)


if __name__ == "__main__":
    main()
