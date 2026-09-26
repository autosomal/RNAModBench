#!/usr/bin/env python
"""R1-4 / R3-8 / R3-9: what the negative controls actually show.

The reviewers objected that (i) KO/KD libraries are not clean negatives and
(ii) the "near-perfect IVT false-positive control" claim was overstated.  Both
are answered with the replicate-aware evaluation tables, reported per
*sequencing unit* so a run that was split in two cannot be counted twice.

Figures in ``04_revision_analysis/negative_controls/figures``:

  FigR11_ivt_null          false-positive density of every unmodified control,
                           per tool, chemistry and modification type
  FigR12_partial_negatives KO/KD depletion on the shared measurable universe,
                           with the purified-site count that depends on it

Tables: ivt_null_per_unit.tsv, ko_kd_paired.tsv.
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

from common.config import SITES_ROOT                    # noqa: E402
from common.figstyle import apply as apply_style        # noqa: E402
from common.figstyle import save                        # noqa: E402

EVAL = (_RB / "data/evaluation/tables")
OUT = (_RB / "analysis/negative_controls")
TAB, FIG = OUT / "tables", OUT / "figures"
apply_style()

WT_COLOR, NULL_COLOR = "#1f5c8b", "#a33f3f"


# --------------------------------------------------------------------------- #
def fig_ivt_null() -> None:
    d = pd.read_csv(EVAL / "controls_ivt_fpr.tsv", sep="\t")
    # one observation per (tool, mod, system, sequencing unit)
    per_unit = (d.groupby(["species", "dataset_group", "mod_type", "tool",
                           "sequencing_unit"], as_index=False)
                 .agg(calls=("n_calls", "mean"),
                      fp_per_1e6=("fp_per_1e6_candidates", "mean"),
                      universe=("n_universe", "mean")))
    agg = (per_unit.groupby(["species", "dataset_group", "mod_type", "tool"],
                            as_index=False)
                .agg(n_units=("sequencing_unit", "nunique"),
                     fp_per_1e6=("fp_per_1e6", "median"),
                     fp_lo=("fp_per_1e6", "min"),
                     fp_hi=("fp_per_1e6", "max")))
    agg.to_csv(TAB / "ivt_null_per_unit.tsv", sep="\t", index=False)

    systems = (agg.groupby(["species", "dataset_group", "mod_type"]).size()
               .reset_index()[["species", "dataset_group", "mod_type"]]
               .values.tolist())
    systems.sort(key=lambda s: (s[2] != "m6A", s[0], s[1]))
    n = len(systems)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.7 * ncols, 4.3 * nrows),
                             constrained_layout=True, squeeze=False)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    for ax, (sp, grp, mod) in zip(axes.ravel(), systems):
        sub = agg[(agg.species == sp) & (agg.dataset_group == grp)
                  & (agg.mod_type == mod)].sort_values("fp_per_1e6")
        y = np.arange(len(sub))
        ax.errorbar(sub.fp_per_1e6, y,
                    xerr=[sub.fp_per_1e6 - sub.fp_lo, sub.fp_hi - sub.fp_per_1e6],
                    fmt="o", color=WT_COLOR, ms=5, lw=1.2, capsize=2)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r.tool} (n={r.n_units})" for r in sub.itertuples()],
                           fontsize=10)
        ax.set_xscale("symlog", linthresh=1)
        ax.set_xticks([0, 1, 10, 100, 1000, 10000])
        ax.xaxis.set_major_formatter(plt.ScalarFormatter())
        ax.minorticks_off()
        ax.set_xlabel("FP per 10$^6$ candidate sites")
        ax.set_title(f"{grp} — {mod}", fontweight="bold", fontsize=11)
        ax.axvline(1, color="0.6", lw=0.8, ls=(0, (2, 2)))
    stem = FIG / "FigR11_ivt_null"
    save(fig, str(stem))
    print("wrote", stem.name)


# --------------------------------------------------------------------------- #
def fig_partial_negatives() -> None:
    ko = pd.read_csv(EVAL / "ko_kd_metrics.tsv", sep="\t")
    pur = pd.read_csv(EVAL / "purified_sites.tsv", sep="\t")
    ko = ko.merge(pur[["species", "sample_wt", "sample_ko", "tool", "n_purified_sites"]],
                  on=["species", "sample_wt", "sample_ko", "tool"], how="left")
    ko.to_csv(TAB / "ko_kd_paired.tsv", sep="\t", index=False)

    groups = [("Arabidopsis", "fip37 knockdown"), ("Mouse", "Mettl3 knockout")]
    fig, axes = plt.subplots(2, 2, figsize=(6.0 * 2, 5.0 * 2),
                             constrained_layout=True)
    for k, (sp, lab) in enumerate(groups):
        sub = ko[ko.species == sp]
        if sub.empty:
            for i in (0, 1):
                axes[k, i].axis("off")
            continue
        ax = axes[k, 0]
        tools = sorted(sub.tool.unique())
        for i, t in enumerate(tools):
            s = sub[sub.tool == t]
            ax.plot([i, i], [s.wt_fraction.min(), s.wt_fraction.max()],
                    color="0.6", lw=1.0, zorder=1)
            ax.plot([i - .12] * len(s), s.wt_fraction, "o", color=WT_COLOR, ms=6,
                    zorder=3, label="WT" if i == 0 else None)
            ax.plot([i + .12] * len(s), s.ko_fraction, "o", color=NULL_COLOR, ms=6,
                    zorder=3, label=lab.split()[0] if i == 0 else None)
        ax.set_xticks(range(len(tools)))
        ax.set_xticklabels(tools, rotation=55, ha="right", fontsize=10)
        ax.set_ylabel("detected sites / shared measurable universe")
        ax.set_ylim(0, None)
        ax.set_title(f"{sp}: depletion on the shared universe", fontweight="bold")
        ax.legend(frameon=False, fontsize=10)

        ax = axes[k, 1]
        rows = (sub.groupby("tool")
                    .agg(ratio=("ko_wt_ratio", "median"),
                         ratio_lo=("ko_wt_ratio", "min"),
                         ratio_hi=("ko_wt_ratio", "max"),
                         purified=("n_purified_sites", "median"),
                         ko_only=("ko_only", "median")))
        rows = rows.sort_values("ratio")
        y = np.arange(len(rows))
        ax.errorbar(rows.ratio, y,
                    xerr=[rows.ratio - rows.ratio_lo, rows.ratio_hi - rows.ratio],
                    fmt="o", color="#4a4a4a", ms=5, lw=1.2, capsize=2)
        ax.axvline(1, color="0.5", lw=1.0, ls=(0, (3, 3)))
        ax.axvline(0, color=NULL_COLOR, lw=1.2)
        ax.set_yticks(y)
        ax.set_yticklabels(rows.index, fontsize=10)
        ax.set_xscale("symlog", linthresh=0.1)
        ax.set_xlabel("KO/KD : WT detection ratio on shared universe")
        ax.set_title(f"{sp}: residual signal (0 = clean negative)",
                     fontweight="bold")
    stem = FIG / "FigR12_partial_negatives"
    save(fig, str(stem))
    print("wrote", stem.name)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    fig_ivt_null()
    fig_partial_negatives()


if __name__ == "__main__":
    main()
