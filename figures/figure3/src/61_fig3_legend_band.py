#!/usr/bin/env python3
"""Figure 3, layout v2: the single legend band under the four rows.

The band carries every key the four panels need and nothing else:

  group 1  13 m6A tools, marker + colour exactly as in panel A (read from
           ``tables/fig3a_tool_style.tsv``, written by `56_...`), two columns;
  group 2  condition colours of panels B/C/D (wild type / treated);
  group 3  panel D linetypes (DENA / m6Anet / Nanom6A) and the panel C/D line
           widths (majority consensus thick / single unit thin).

No axes, no frame, Arial, minimum text 7.2 pt (house floor).  Output:
``figures/Fig3_legend_band.{pdf,png}`` at the figure width (6.66 x 0.42 in).
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

import matplotlib as mpl
import pandas as pd

sys.path.insert(0, str(_RB / "src/harmonisation"))
import common.config as C                                          # noqa: E402
from common.figstyle import apply as apply_style                    # noqa: E402
from common.manifest import setup_logger                            # noqa: E402
from matplotlib.lines import Line2D                                # noqa: E402
import matplotlib.pyplot as plt                                    # noqa: E402

OUT = (_RB / "figures/figure3")
TAB, FIG, LOG = OUT / "tables", OUT / "figures", OUT / "logs"
W, H = 6.66, 0.34
FS = 7.2
WT_COLOR, TRT_COLOR = "#3778A0", "#F5B264"


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    logger = setup_logger("61_fig3_legend_band", log_dir=LOG)
    apply_style()
    style = pd.read_csv(TAB / "fig3a_tool_style.tsv", sep="\t")

    fig = plt.figure(figsize=(W, H))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

    
    #: manuscript use "WT"; one wording for the condition everywhere.
    conditions = [
        Line2D([], [], linestyle="none", marker="o", markersize=4.2,
               markerfacecolor=WT_COLOR, markeredgecolor="none", label="WT"),
        Line2D([], [], linestyle="none", marker="o", markersize=4.2,
               markerfacecolor=TRT_COLOR, markeredgecolor="none", label="treated"),
        Line2D([], [], color="0.35", lw=0.9, label="majority consensus"),
        Line2D([], [], color="0.35", lw=0.35, label="individual unit"),
    ]
    leg2 = fig.legend(handles=conditions, loc="center left",
                      bbox_to_anchor=(0.10, 0.5), frameon=False, fontsize=FS,
                      ncol=2, handletextpad=0.42, labelspacing=0.30,
                      borderaxespad=0.0)
    fig.add_artist(leg2)

    lines = [Line2D([], [], color="0.35", lw=0.9, ls=ls, label=t)
             for t, ls in (("DENA", "solid"), ("m6Anet", "dashed"),
                           ("Nanom6A", "dotted"))]
    leg3 = fig.legend(handles=lines, loc="center left",
                      bbox_to_anchor=(0.55, 0.5), frameon=False, fontsize=FS,
                      ncol=2, handletextpad=0.42, labelspacing=0.30,
                      borderaxespad=0.0)
    fig.add_artist(leg3)

    stem = FIG / "Fig3_legend_band"
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=300)
    plt.close(fig)
    logger.info("wrote %s.{pdf,png} (%.2f x %.2f in)", stem.name, W, H)


if __name__ == "__main__":
    main()
