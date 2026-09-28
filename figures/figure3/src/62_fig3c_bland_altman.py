#!/usr/bin/env python3
"""Panel C (layout v5): binned Bland-Altman agreement of the WT-treated pairs.

One small panel per tool x species (4 x 3), each with its own axes:

  x = mean of the two modification ratios of a site,
  y = difference (WT - treated),
  solid line = bias (mean difference), dashed = 95% limits of agreement
  (bias +/- 1.96 SD), dotted = zero.

The numbers themselves are never printed on the figure (house rule) - they are
in `tables/fig3b_bland_altman.tsv`; the points come from
`tables/fig3b_ba_points.tsv.gz`, written by `57_fig3b_modratio_wt_treatment.py`
on exactly the pairing rule used by every other table of this figure.
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

from common import config as C                                        # noqa: E402
from common.figstyle import apply as apply_style                      # noqa: E402
from common.manifest import setup_logger                              # noqa: E402

OUT = (_RB / "figures/figure3")
TAB, FIG, LOG = OUT / "tables", OUT / "figures", OUT / "logs"
W, H = 6.66, 1.62
FS = {"title": 8.6, "label": 7.6, "tick": 7.2}
SPECIES = ["Arabidopsis", "Mouse", "Human"]
TOOLS = ["m6Anet", "MINES", "Nanom6A", "DENA"]
EFF_COLOR = "#4a4a4a"
TOOL_COLORS = {"m6Anet": "#2166AC", "MINES": "#B2182B",
               "Nanom6A": "#E69F00", "DENA": "#4D9221"}


def main() -> None:
    logger = setup_logger("62_fig3c_bland_altman", log_dir=LOG)
    apply_style()
    pts = pd.read_csv(TAB / "fig3b_ba_points.tsv.gz", sep="\t")

    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(1, len(SPECIES), left=0.075, right=0.995, top=0.86,
                          bottom=0.16, wspace=0.30)
    edges = np.arange(0.05, 1.0, 0.05)
    for ci, species in enumerate(SPECIES):
        ax = fig.add_subplot(gs[0, ci])
        for tool in TOOLS:
            sub = pts[(pts.species == species) & (pts.tool == tool)]
            if sub.empty:
                continue
            b = pd.cut(sub["mean"], edges)
            g = sub.groupby(b, observed=True)["diff"]
            x = np.array([iv.mid for iv in g.groups.keys()])
            m = g.mean().to_numpy()
            sd = g.std(ddof=1).to_numpy()
            n = g.count().to_numpy()
            keep = n >= 200
            ci95 = 1.96 * sd / np.sqrt(np.maximum(n, 1))
            ax.errorbar(x[keep], m[keep], yerr=ci95[keep], color=TOOL_COLORS[tool],
                        lw=1.0, marker="o", ms=2.6, capsize=0, elinewidth=0.8,
                        zorder=3)
        ax.axhline(0, color="0.55", lw=0.7, ls=(0, (2.5, 2.5)), zorder=1)
        ax.set_xlim(0.10, 0.90)
        ax.set_ylim(-0.35, 0.35)
        ax.set_xticks([0.2, 0.4, 0.6, 0.8])
        ax.set_yticks([-0.2, 0.0, 0.2])
        ax.tick_params(labelsize=FS["tick"], length=2.5, width=0.7,
                       labelleft=(ci == 0))
        ax.set_title(species, fontsize=FS["title"], fontweight="bold", pad=3)
        if ci == 0:
            ax.set_ylabel("\u0394 ratio (WT \u2212 treated)", fontsize=FS["label"],
                          labelpad=2)
        if ci == 1:
            ax.set_xlabel("mean modification ratio", fontsize=FS["label"], labelpad=2)
    fig.text(0.004, 0.97, "C", fontsize=FS["title"], fontweight="bold",
             ha="left", va="top")
    FIG.mkdir(parents=True, exist_ok=True)
    stem = FIG / "Fig3C_bland_altman"
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=300)
    plt.close(fig)
    logger.info("wrote %s.{pdf,png} (%.2f x %.2f in)", stem.name, W, H)


if __name__ == "__main__":
    main()
