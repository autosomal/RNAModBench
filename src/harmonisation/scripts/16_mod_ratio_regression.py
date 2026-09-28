#!/usr/bin/env python
"""
Revision Fig. S3C (redo, per-replicate) — mod-ratio regression panel.

Same VISUAL as the manuscript's sup3 panel C (x = GLORI modification ratio,
y = tool predicted modification ratio, one regression line per tool, dashed
identity line, r/n in the legend) but computed PER BIOLOGICAL REPLICATE and
with agreement metrics (Lin's CCC) added, to answer R3-2 (replicates) and
R3-5 (association vs agreement).

Reads the cached tables produced by 15_mod_ratio_replicate_agreement.py
(sources callsets).  Outputs a square, Arial, no-gridline, closed-box panel.
"""

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
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

OUT = str(_RB / "analysis/mod_ratio_replicates")
TAB = f"{OUT}/tables"
FIG = f"{OUT}/figures"

rcParams.update({
    "font.family": "Arial", "font.size": 16, "axes.titlesize": 20,
    "axes.labelsize": 18, "axes.linewidth": 2.0, "xtick.labelsize": 14,
    "ytick.labelsize": 14, "legend.fontsize": 11.5, "axes.grid": False,
    "figure.dpi": 120, "savefig.dpi": 300,
})

PANELS = [("Arabidopsis", "Arabidopsis_WT"), ("Mouse", "Mouse_WT"), ("Human", "HeLa_WT")]
TOOLS = ["m6Anet", "Nanom6A", "DENA", "MINES"]
TOOL_COLOR = {"m6Anet": "#3aae9a", "Nanom6A": "#e07b1a", "DENA": "#a6d3ea", "MINES": "#3b6fb0"}
REP_ALPHA = [0.55, 0.55, 0.55, 0.55]


def box(ax):
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(2.0); s.set_color("black")
    ax.grid(False)


def fit(xs, ys):
    if len(xs) < 2 or np.allclose(np.std(xs), 0):
        return None
    b, a = np.polyfit(xs, ys, 1)
    return a, b


def main():
    sites = pd.read_csv(f"{TAB}/mod_ratio_matched_sites.tsv", sep="\t")
    summ = pd.read_csv(f"{TAB}/mod_ratio_summary_by_group.tsv", sep="\t")

    fig, axes = plt.subplots(1, 3, figsize=(7.2 * 3, 7.2))
    for ax, (sp, grp) in zip(axes, PANELS):
        box(ax)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_box_aspect(1)
        ax.plot([0, 1], [0, 1], ls="--", lw=1.2, color="black", zorder=1)
        legend_handles = []
        for tool in TOOLS:
            sub = sites[(sites.species == sp) & (sites.group == grp) & (sites.tool == tool)]
            col = TOOL_COLOR[tool]
            drew = False
            for rep in sorted(sub["replicate_tag"].dropna().unique()):
                rs = sub[sub.replicate_tag == rep]
                g = rs["glori_ratio"].values; t = rs["tool_ratio"].values
                f = fit(g, t)
                if f is None:
                    continue
                a, b = f
                xx = np.linspace(0, 1, 2)
                ax.plot(xx, np.clip(a + b * xx, 0, 1), color=col, lw=1.3,
                        alpha=0.6, zorder=2)
                drew = True
            # bold pooled-fit line (the "headline" line, like the original panel)
            f = fit(sub["glori_ratio"].values, sub["tool_ratio"].values)
            if f is not None:
                a, b = f
                xx = np.linspace(0, 1, 2)
                ax.plot(xx, np.clip(a + b * xx, 0, 1), color=col, lw=2.6, zorder=3)
            # legend text from the per-replicate summary
            row = summ[(summ.species == sp) & (summ.group == grp) & (summ.tool == tool)]
            if len(row):
                r = row["pearson_r_mean"].values[0]; sd = row["pearson_r_sd"].values[0]
                ccc = row["ccc_mean"].values[0]; n = row["n_overlap_mean"].values[0]
                sd_txt = f"±{sd:.2f}" if pd.notna(sd) else ""
                lab = f"{tool}  r={r:.2f}{sd_txt}, CCC={ccc:.2f}, n≈{n/1000:.1f}k"
            else:
                lab = tool
            h = plt.Line2D([], [], color=col, lw=2.6, label=lab)
            legend_handles.append(h)
        ax.set_title(sp)
        ax.set_xlabel("GLORI Modification Ratio")
        if sp == "Arabidopsis":
            ax.set_ylabel("Tool Predicted Modification Ratio")
        ax.legend(handles=legend_handles, loc="upper left", frameon=True,
                  edgecolor="black", fontsize=11.5)
    fig.suptitle("Tool-predicted vs GLORI modification ratio — regression per "
                 "biological replicate (bold = pooled fit)", y=1.02)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(f"{FIG}/mod_ratio_regression_S3C_style.png", bbox_inches="tight")
    fig.savefig(f"{FIG}/mod_ratio_regression_S3C_style.pdf", bbox_inches="tight")
    print("[done]", f"{FIG}/mod_ratio_regression_S3C_style.png")


if __name__ == "__main__":
    main()
