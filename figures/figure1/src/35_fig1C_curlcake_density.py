#!/usr/bin/env python
"""Figure 1C, replicate-aware rebuild (revision request 2026-09-20).

The published panel plotted one line per tool: the total number of Curlcake
calls per dataset (Curlcake_m6A, Curlcake_IVT), computed from replicate-union
aggregates -- one number per dataset, no replicate information (R3-2).

This rebuild keeps the published layout (two datasets, tools on x ordered by
the m6A mean) but every plotted number is the mean over the dataset's
independent sequencing units (Curlcake_IVT: rep1 + rep3, rep2_partial belongs
to the same run; Curlcake_m6A: rep1 + rep2), with one dot per unit.  Tools that
were never run on a dataset get no marker (no fake zeros).

Outputs
-------
figures/figure1/figures/panels/Figure1_rev_C_curlcake_counts.{pdf,png}
figures/figure1/tables/fig1C_tool_counts.tsv

Usage
-----
conda run -n benchmark-revision --no-capture-output python scripts/35_fig1C_curlcake_density.py
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

HERE = _RB / "src/harmonisation"
sys.path.insert(0, str(HERE))

from common import config as C                    # noqa: E402
from common.figstyle import apply as apply_style  # noqa: E402
from common.figstyle import save                  # noqa: E402
from common.io_utils import write_table           # noqa: E402
from common.manifest import setup_logger          # noqa: E402

OUT = (_RB / "figures/figure1")
PANEL_DIR = (_RB / "figures/figure1/figures/panels")
TAB_DIR = (_RB / "figures/figure1/tables")

PLATFORM = "RNA002"
MOD = "m6A"
COLORS = {"Curlcake_m6A": "#F5B264", "Curlcake_IVT": "#3778A0"}
CONDS = ("Curlcake_m6A", "Curlcake_IVT")
DISPLAY = {"yanocomp": "Yanocomp"}


def unit_members(group: str) -> dict[str, list[str]]:
    """{sequencing-unit key: [representative first, other members]}, as in 28."""
    members = [s for s in C.SAMPLES
               if s.dataset_group == group and s.platform == PLATFORM]
    names = {s.canonical for s in members}
    out = {}
    for unit, mem in C.independent_units(members).items():
        rep = unit if unit in names else mem[0]
        out[unit] = [rep] + [m for m in mem if m != rep]
    return out


def unit_counts(group: str, tool: str, logger) -> dict[str, float]:
    """{unit key: distinct (chrom, pos_raw) sites}, one number per unit.

    Same dedup key and representative-fallback rule as 28_fig_tool_counts.
    """
    rows = []
    for unit, members in unit_members(group).items():
        n = None
        used = None
        for name in members:                    # representative first
            f = (C.CALLSET_ROOT / PLATFORM / "Curlcake" / group / MOD
                 / tool / f"{name}.tsv")
            if f.exists():
                df = pd.read_csv(f, sep="\t", usecols=["chrom", "pos_raw"],
                                 dtype=str, keep_default_na=False)
                n = int(df.drop_duplicates().shape[0])
                used = name
                break
        if n is None:
            logger.info("no callset for %s/%s (unit %s)", group, tool, unit)
            continue
        rows.append({"dataset_group": group, "tool": tool, "unit": unit,
                     "sample": used, "n_sites": n})
    return rows


def main() -> None:
    logger = setup_logger("35_fig1C_curlcake_density")
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    PANEL_DIR.mkdir(parents=True, exist_ok=True)

    inputs: list[dict] = []
    means: dict[str, dict[str, float]] = {c: {} for c in CONDS}
    for cond in CONDS:
        for tool in sorted(C.article_tool_scope(PLATFORM, "Curlcake", cond, MOD)):
            rows = unit_counts(cond, tool, logger)
            inputs.extend(rows)
            if rows:
                means[cond][tool] = float(np.mean([r["n_sites"] for r in rows]))

    tab = pd.DataFrame(inputs)
    write_table(tab, (_RB / "figures/figure1/tables/fig1C_tool_counts.tsv"))
    logger.info("units tabulated: %d", len(tab))

    # tools that were run on either condition; ordered by the m6A mean, as the
    # published panel did
    order = sorted({t for c in CONDS for t in means[c]},
                   key=lambda t: -means["Curlcake_m6A"].get(t, 0.0))

    apply_style()
    fig, ax = plt.subplots(figsize=(5.6, 2.54), constrained_layout=True)
    x = np.arange(len(order))
    for cond in CONDS:
        col = COLORS[cond]
        ys = [means[cond].get(t, np.nan) for t in order]
        ax.plot(x, ys, "-", marker="o" if cond == "Curlcake_m6A" else "s",
                ms=6.5, lw=1.8, color=col, zorder=3,
                label=f"{cond} (mean of {len(unit_members(cond))})")
        # one dot per independent sequencing unit
        for i, t in enumerate(order):
            for r in tab[(tab.dataset_group == cond) & (tab.tool == t)].itertuples():
                ax.plot(i, r.n_sites, "o", ms=4.2, mfc=col, mec="white",
                        mew=.7, alpha=.95, zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels([DISPLAY.get(t, t) for t in order], rotation=45,
                       ha="right", fontsize=10)
    ax.set_ylabel("Counts")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=10, loc="upper right")
    FIG = (_RB / "figures/figure1/figures/panels/Figure1_rev_C_curlcake_counts")
    save(fig, str(FIG))
    logger.info("wrote %s", FIG.with_suffix(".pdf").name)


if __name__ == "__main__":
    main()
