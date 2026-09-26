#!/usr/bin/env python
"""Figure 1D, replicate-aware rebuild (user request 2026-09-20).

Published panel: grouped horizontal bars of the number of detected Curlcake
sites that fall inside the canonical RRACH motif, one bar pair per tool,
computed from replicate-union aggregates (no replicate information, R3-2).

Rebuild: same layout and colours, but each bar is the mean over the dataset's
independent sequencing units (Curlcake_IVT: rep1 + rep3, rep2_partial belongs
to the same run; Curlcake_m6A: rep1 + rep2) and one dot per unit sits at the
bar end.  Legacy conventions kept: RRACH = ``[AG][AG]AC[ACT]`` on the 5-mer
centred at the site, and the per-construct ``Start >= 37`` trim (0-based
``pos_raw >= 36``); no centre-base filter is applied (that filter post-dates
the published panel).  Tools that were never run on a dataset get no bar (no
fake zeros -- the published panel drew 0 bars for the two IVT tools that were
never run, DRUMMER / ELIGOS2_diff).

Outputs
-------
figures/figure1/figures/panels/Figure1_rev_D_rrach_counts.{pdf,png}
figures/figure1/tables/fig1D_rrach_counts.tsv

Usage
-----
conda run -n benchmark-revision --no-capture-output python scripts/36_fig1D_rrach_counts.py
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
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
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
#: canonical RRACH on the DNA alphabet, as the legacy candidate generator
#: (DENA --motif 'RRACH') enumerated it
RRACH = re.compile(r"[AG][AG]AC[ACT]")
#: legacy per-construct trim (1-based Start >= 37); pos_raw is 0-based
MIN_POS_1BASED = 37


def unit_members(group: str) -> dict[str, list[str]]:
    members = [s for s in C.SAMPLES
               if s.dataset_group == group and s.platform == PLATFORM]
    names = {s.canonical for s in members}
    out = {}
    for unit, mem in C.independent_units(members).items():
        rep = unit if unit in names else mem[0]
        out[unit] = [rep] + [m for m in mem if m != rep]
    return out


def rrach_counts(group: str, tool: str, logger) -> list[dict]:
    """One row per independent unit: RRACH-inside site count (deduplicated)."""
    rows = []
    for unit, members in unit_members(group).items():
        n = None
        used = None
        for name in members:
            f = (C.CALLSET_ROOT / PLATFORM / "Curlcake" / group / MOD
                 / tool / f"{name}.tsv")
            if not f.exists():
                continue
            df = pd.read_csv(f, sep="\t",
                             usecols=["chrom", "pos_raw", "five_mer_raw"],
                             dtype=str, keep_default_na=False)
            ok = df.five_mer_raw.str.upper().str.match(RRACH)
            df = df[ok].copy()
            pos = pd.to_numeric(df.pos_raw)
            keep = pos + 1 >= MIN_POS_1BASED
            n = int(df[keep].drop_duplicates(["chrom", "pos_raw"]).shape[0])
            used = name
            break
        if n is None:
            logger.info("no callset for %s/%s", group, tool)
            continue
        rows.append({"dataset_group": group, "tool": tool, "unit": unit,
                     "sample": used, "n_rrach": n})
    return rows


def main() -> None:
    logger = setup_logger("36_fig1D_rrach_counts")
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    PANEL_DIR.mkdir(parents=True, exist_ok=True)

    tab = pd.DataFrame()
    means: dict[str, dict[str, float]] = {c: {} for c in CONDS}
    for cond in CONDS:
        rows_all: list[dict] = []
        for tool in sorted(C.article_tool_scope(PLATFORM, "Curlcake", cond, MOD)):
            rows = rrach_counts(cond, tool, logger)
            rows_all.extend(rows)
            if rows:
                means[cond][tool] = float(np.mean([r["n_rrach"] for r in rows]))
        tab = pd.concat([tab, pd.DataFrame(rows_all)], ignore_index=True)
    write_table(tab, (_RB / "figures/figure1/tables/fig1D_rrach_counts.tsv"))

    order = sorted({t for c in CONDS for t in means[c]},
                   key=lambda t: -means["Curlcake_m6A"].get(t, 0.0))

    apply_style()
    fig, ax = plt.subplots(figsize=(5.6, 2.54), constrained_layout=True)
    y = np.arange(len(order))          # barh: index 0 at the bottom, so the
    h, off = 0.34, 0.19                # m6A-largest tool ends up at the
    for i, t in enumerate(order):      # bottom, exactly like the published panel
        for cond, dy in (("Curlcake_m6A", +off), ("Curlcake_IVT", -off)):
            if t not in means[cond]:
                continue               # never run on this dataset: no bar
            col = COLORS[cond]
            m = means[cond][t]
            ax.barh(i + dy, m, height=h, color=col, edgecolor="black",
                    linewidth=.8, zorder=3,
                    label=None)
            units = tab[(tab.dataset_group == cond) & (tab.tool == t)]
            ax.plot(units.n_rrach, np.full(len(units), i + dy), "o", ms=4.2,
                    mfc=col, mec="black", mew=.6, zorder=4)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=COLORS[c], ec="black", lw=.8)
               for c in CONDS]
    labels = [f"{c} (mean of {len(unit_members(c))})" for c in CONDS]
    ax.set_yticks(y)
    ax.set_yticklabels([DISPLAY.get(t, t) for t in order], fontsize=10)
    ax.set_xlabel("Counts (RRACH)")
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.legend(handles, labels, frameon=False, fontsize=10, loc="upper right")
    FIG = (_RB / "figures/figure1/figures/panels/Figure1_rev_D_rrach_counts")
    save(fig, str(FIG))
    logger.info("wrote %s", FIG.with_suffix(".pdf").name)


if __name__ == "__main__":
    main()
