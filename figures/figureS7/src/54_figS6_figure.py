#!/usr/bin/env python3
"""54 -- Figure S6 (rebuilt 2026-09-21): non-m6A tools against unmodified controls.

Panel contract (approved 2026-09-21).  The questions of R3-9 come first and the
HeLa call numbers are compressed into one panel:

A  unmodified Curlcake controls -- in-universe calls per construct (filled = the
   two independent constructs, open = the depth-matched subset, bar = mean of
   the independent constructs), for the five tools that were actually run on the
   synthetic constructs.  CHEUI-m5C has no column here at all (it was never run
   on them); the caption states that, which is what R1-6 asks for;
B  overlap within a condition against overlap between conditions, one point per
   replicate pair (nine WT x IVT pairs on the shared candidate universe);
C  where each unit's calls sit inside the tool's own 0-1 score: one point (that
   unit's median) and one whisker (its interquartile range) per independent
   sequencing unit, WT against unmodified IVT.  This is the mechanism behind the
   R3-9 failures (saturation, a narrow band, overlap, inversion) rather than a
   separation statistic -- the pooled AUC is Fig. 7F and the per-pair AUCs are
   tabulated; the reproducibility numbers (k-of-n) are quoted in the caption;
D  HeLa calls per unit plus the unmodified-IVT-to-WT ratio: the ratio strip
   carries the mean-of-counts ratio with its unit-level bootstrap interval, and
   the coverage-matched and raw-union values are tabulated in Table S5.  The
   union-ratio diamond was retired on 2026-09-28 with the rest of the pooled
   reporting units.

Layout rules (2026-09-21, third revision -- " fonts too small relative to the plot / tidier /
D may be slanted / compress the canvas "):

* **compact canvas**: 170 x 240 mm (481.9 x 680.3 pt), i.e. the SI print width,
  instead of A4.  At 10 pt ticks that is a font/page ratio of 2.1 % versus
  1.6 % on A4, which is what made the A4 draft read as "text too small for the
  figure";
* **full-width rows**: B, C and D share the six-column tool grid, so their
  columns line up top to bottom; A has its own five-column axis because it
  covers only the five tools that were run on the synthetic constructs (the
  user asked for the CHEUI-m5C column to disappear from A on 2026-09-21, and
  accepted that A's columns therefore no longer sit under B's).  Every row
  carries the tool names at 45 degrees;
* **no titles anywhere** (house rule of 2026-09-19): panels are identified by
  their bold letters in the page margin plus the caption, tools by their axis
  labels.  Panel C carries the tool name as the facet *x-axis label*;
* **all tool labels are drawn at 45 degrees** (the only orientation that fits
  the 0.92 in column pitch: a horizontal label is 0.96 in wide, a 90 deg label
  hangs 0.96 in below the axis, 45 deg needs 0.77 in);
* legends live *inside* the panels in the empty corner of each drawing area and
  a post-draw audit converts every legend box to data coordinates and asserts
  that no data point falls inside it (`legend_data_violations`); the y ranges of
  A, B and D carry the headroom that keeps those corners empty;
* panel C is a 0.55 in **score-location strip** since the 2026-09-21 fifth
  revision: one median point + one q25-q75 whisker per independent unit, WT
  (filled) against unmodified IVT (open), each tool on its own native 0-1 score.
  It replaces both the six density facets of the third version (they mixed three
  score kinds and needed a rescaling per facet to be legible) and the AUC strip
  of the fourth version (thin -- the pooled per-tool AUC is Fig. 7F -- and it
  read as "the benchmark cannot separate the conditions" instead of as the
  mechanism the reviewer asks about).  The densities stay in
  ``tables/s6_score_density_per_unit.tsv``, the AUC evidence in
  ``tables/s6_score_separation.tsv`` and ``tables/s6_score_location_per_unit.tsv``
  carries exactly the seven numbers this panel draws;
* log axes carry named ticks only -- no tick is ever labelled "0";
* the finished page goes through ``pagelayout.assert_page_clean``.

Usage:
    conda run -n benchmark-revision --no-capture-output python \
        figures/figureS7/src/54_figS6_figure.py
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
import json
import sys
import time
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve()
PROJECT = _RB
COMMON = (_RB / "src/harmonisation/common")
sys.path.insert(0, str(COMMON))

import figstyle                                     # noqa: E402
import pagelayout as pl                             # noqa: E402

OUT = (_RB / "figures/figureS7")
TAB = (_RB / "figures/figureS7/tables")
FIG = (_RB / "figures/figureS7/figures")
LOG = (_RB / "figures/figureS7/logs")

# --------------------------------------------------------------------------- #
# palette, tool order, style
# --------------------------------------------------------------------------- #
C_WT, C_IVT = "#4B81B8", "#E8A76B"        # HeLa WT / unmodified IVT
C_CROSS = "0.35"                          # between-condition comparison
C_ACC = "#B3261E"                         # control burden (Currlake / subset)
C_DARK, C_GREY = "0.15", "0.55"

CLASSES: list[tuple[str, list[str]]] = [
    ("FP-dominated", ["CHEUI_m5C", "NanoMUD_psi", "NanoMUD_m1psi"]),
    ("Intermediate", ["NanoNm"]),
    ("Sparse", ["NanoPsu", "NanoSPA_psU"]),
]
TOOL_ORDER: list[str] = [t for _c, ts in CLASSES for t in ts]
DISPLAY = {
    "CHEUI_m5C": "CHEUI-m5C",
    "NanoMUD_psi": "NanoMUD-\u03a8",
    "NanoMUD_m1psi": "NanoMUD-m1\u03a8",
    "NanoNm": "NanoNm",
    "NanoPsu": "NanoPsu",
    "NanoSPA_psU": "NanoSPA-\u03a8",
}
WT_UNITS = ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"]
IVT_UNITS = ["HeLa_IVT_rep1", "HeLa_IVT_rep2", "HeLa_IVT_rep3"]
#: CHEUI-m5C has no synthetic-construct run in this study
NO_CURLCAKE = ["CHEUI_m5C"]

#: one notch up from the A4 draft: on the smaller canvas these read as a normal
#: journal figure, and the smallest (facet ticks) stays above the 7 pt floor
STYLE = dict(letter=18.0, tick=10.0, tooltick=9.5, label=11.5, legend=10.5)

# --------------------------------------------------------------------------- #
# page geometry -- 170 x 240 mm canvas, four full-width rows, one tool grid
# --------------------------------------------------------------------------- #
PAGE_W, PAGE_H = 170.0 / 25.4, 240.0 / 25.4        # 6.693 x 9.449 in
LEFT, RIGHT, TOP, BOTTOM = 0.85, 0.12, 0.30, 0.28
CONTENT_L = LEFT
CONTENT_W = PAGE_W - LEFT - RIGHT                  # 5.723 in
NCOLS = len(TOOL_ORDER)
XLIM = (-0.6, NCOLS - 0.4)                         # data x = tool index
#: panel A covers only the five tools that were run on the synthetic constructs
#: (CHEUI-m5C was not), so it gets its own five-column axis; B, C and D keep the
#: shared six-column grid and therefore line up with each other only
A_TOOLS: list[str] = [t for t in TOOL_ORDER if t not in NO_CURLCAKE]
XLIM_A = (-0.6, len(A_TOOLS) - 0.4)

#: C is the AUC separation strip (0.55 in), so PLOT_H carries the three full
#: plotting rows only; the height that C released went to A, B and D

#: names, and the three panels should read as one block".  B and C therefore gave
#: up their 0.82 in label bands; the released height went into the plots, and C's
#: strip (which was cut off at 0.55 in) more than doubled.
PLOT_H = {"A": 1.73, "B": 1.73, "D": 1.54}
C_STRIP_H = 1.15
STRIP_H, STRIP_GAP = 0.50, 0.10                    # panel D ratio strip
#: a 45 deg tool label is measured, not estimated: at 9.5 pt the bbox of
#: "NanoMUD-m1\u03a8" is a 0.77 in square, i.e. it hangs 0.77 in below the axis
#: (a 90 deg label would need 0.96 in, a horizontal one is 0.96 in wide and
#: does not fit the 0.92 in column pitch -- 45 deg is the only orientation that
#: fits every row)
XLAB_H = 0.82
ROW_GAP = 0.10

_LAYOUT: list[dict] = []
#: every log-scaled axis on the page, with its tick labels (the "no fake 0"
#: contract: a log axis may never carry a tick labelled ``0``)
_LOG_AXES: list[dict] = []
#: legend boxes with data points inside them (the "legend must not cover data"
#: contract of 2026-09-21); empty is the only acceptable value
_LEGEND_HITS: list[dict] = []
#: legends that stick out of their own panel (a legend that overflows is
#: silently clipped by the canvas -- the page gate cannot see that)
_LEGEND_OUT: list[dict] = []
_SLACK = 0.0


def apply_style() -> None:
    figstyle.apply()
    mpl.rcParams.update({
        "font.size": STYLE["tick"],
        "axes.labelsize": STYLE["label"],
        "xtick.labelsize": STYLE["tick"],
        "ytick.labelsize": STYLE["tick"],
        "legend.fontsize": STYLE["legend"],
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.9,
        "xtick.major.width": 0.9,
        "ytick.major.width": 0.9,
        "xtick.major.size": 3.2,
        "ytick.major.size": 3.2,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "figure.dpi": 120,
    })


def rect(top_in: float, height_in: float, *, left_in: float | None = None,
         width_in: float | None = None) -> list[float]:
    """Axes rectangle in figure fractions from inches measured off the page top."""
    left_in = CONTENT_L if left_in is None else left_in
    width_in = CONTENT_W if width_in is None else width_in
    return [left_in / PAGE_W, (PAGE_H - top_in - height_in) / PAGE_H,
            width_in / PAGE_W, height_in / PAGE_H]


def x_of(i: float, *, xlim: tuple[float, float] | None = None,
         width_in: float | None = None) -> float:
    """Page x (inches) of data coordinate ``i`` on a panel's tool grid."""
    xlim = XLIM if xlim is None else xlim
    width_in = CONTENT_W if width_in is None else width_in
    return CONTENT_L + (i - xlim[0]) / (xlim[1] - xlim[0]) * width_in


def jitter(n: int, width: float) -> np.ndarray:
    return np.linspace(-width, width, n) if n > 1 else np.zeros(1)


def tool_axis(ax, *, labels: bool = True, tools: list[str] | None = None,
              xlim: tuple[float, float] | None = None) -> None:
    """Tool index axis; 45 deg labels are the only orientation that fits."""
    tools = TOOL_ORDER if tools is None else tools
    xlim = XLIM if xlim is None else xlim
    ax.set_xlim(*xlim)
    ax.set_xticks(range(len(tools)))
    if labels:
        ax.set_xticklabels([DISPLAY[t] for t in tools], rotation=45,
                           ha="right", rotation_mode="anchor",
                           fontsize=STYLE["tooltick"])
    else:
        ax.set_xticklabels([])


def _legend(ax, handles: list, loc: str, ncol: int = 1,
            fontsize: float | None = None,
            columnspacing: float = 1.3) -> None:
    ax.legend(handles=handles, loc=loc, ncol=ncol, frameon=False,
              fontsize=fontsize or STYLE["legend"], handletextpad=0.40,
              columnspacing=columnspacing, labelspacing=0.30, borderpad=0.0,
              borderaxespad=0.25)


# --------------------------------------------------------------------------- #
# panel A -- unmodified Curlcake controls
# --------------------------------------------------------------------------- #
def panel_a(fig, top_in: float, height_in: float, cc: pd.DataFrame) -> plt.Axes:
    ax = fig.add_axes(rect(top_in, height_in))
    tool_axis(ax, tools=A_TOOLS, xlim=XLIM_A)
    # 160 leaves the top of the panel free for the key (all data are <= 97)
    ax.set_ylim(0.0, 160.0)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_ylabel("Calls in the\ncandidate-site set", y=0.52, labelpad=7.0)
    for i, tool in enumerate(A_TOOLS):
        x = float(i)
        sub = cc[(cc.tool == tool) & (cc.row_type == "construct")]
        ind = sub[sub.construct_role == "independent"].sort_values("construct")
        deep = sub[sub.construct_role != "independent"]
        vals = ind["n_calls_in_universe"].to_numpy(float)
        xs = x - 0.08 + jitter(len(vals), 0.07)
        ax.plot(xs, vals, marker="o", ls="none", ms=4.8, mfc=C_ACC, mec=C_ACC,
                zorder=4)
        if len(deep):
            dv = float(deep["n_calls_in_universe"].iloc[0])
            ax.plot([x + 0.15], [dv], marker="o", ls="none", ms=4.8,
                    mfc="white", mec=C_ACC, mew=1.0, zorder=4)
        mean = float(ind["n_calls_in_universe"].mean())
        ax.plot([x - 0.18, x + 0.18], [mean, mean], color=C_DARK, lw=1.7,
                zorder=5, solid_capstyle="butt")
    # three entries in one row: this panel only carries the five tools that were
    # run on the constructs, so there is no "not run" key to explain
    _legend(ax, [
        Line2D([], [], marker="o", ls="none", ms=4.8, mfc=C_ACC, mec=C_ACC,
               label="Independent"),
        Line2D([], [], marker="o", ls="none", ms=4.8, mfc="white", mec=C_ACC,
               mew=1.0, label="Depth-matched"),
        Line2D([], [], color=C_DARK, lw=1.7, label="Construct mean"),
    ], loc="upper right", ncol=3)

    _LAYOUT.append({"panel": "A", "n_tools": len(A_TOOLS), "gap_marker": False,
                    "grid": "own five-column axis (B/C/D keep the six-column grid)",
                    "n_positive_points": int(
                        (cc[(cc.row_type == "construct")][
                            "n_calls_in_universe"] > 0).sum()),
                    "n_zero_points": int(
                        (cc[(cc.row_type == "construct")][
                            "n_calls_in_universe"] == 0).sum())})
    return ax


# --------------------------------------------------------------------------- #
# panel B -- overlap within a condition and between conditions
# --------------------------------------------------------------------------- #
def panel_b(fig, top_in: float, height_in: float, pairs: pd.DataFrame) -> plt.Axes:
    ax = fig.add_axes(rect(top_in, height_in))
    tool_axis(ax, labels=False)      # names are printed under panel D only
    ax.set_yscale("log")
    # 1.6 leaves the top of the panel free for the single-row key while keeping
    # every pair point (max 0.552) below it
    ax.set_ylim(0.006, 1.6)
    ticks = [0.01, 0.02, 0.05, 0.1, 0.2, 0.4]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{t:g}" for t in ticks])
    ax.minorticks_off()
    _LOG_AXES.append({"panel": "B", "yscale": ax.get_yscale(),
                      "ytick_labels": [t.get_text() for t in ax.get_yticklabels()]})
    ax.set_ylabel("Jaccard index", y=0.52, labelpad=7.0)
    floor = 0.0072                                     # only if a pair is 0
    n_zero = 0
    for i, tool in enumerate(TOOL_ORDER):
        x = float(i)
        d = pairs[pairs.tool == tool]
        for comp, off, col in (("within_WT", -0.22, C_WT),
                               ("cross", 0.0, C_CROSS),
                               ("within_IVT", 0.22, C_IVT)):
            sub = d[d.comparison == comp]
            vals = sub["jaccard_uni"].to_numpy(float)
            zero = vals <= 0.0
            n_zero += int(zero.sum())
            xs = x + off + jitter(len(vals), 0.075 if len(vals) > 3 else 0.05)
            if (~zero).any():
                ax.plot(xs[~zero], vals[~zero], marker="o", ls="none", ms=3.4,
                        mfc=col, mec="none", alpha=0.8, zorder=3)
            if zero.any():
                ax.plot(xs[zero], np.full(int(zero.sum()), floor),
                        marker="o", ls="none", ms=3.4, mfc="white", mec=col,
                        mew=1.0, zorder=4)
            pos = vals[~zero]
            if pos.size:
                ax.plot([x + off - 0.09, x + off + 0.09], [pos.mean()] * 2,
                        color=C_DARK, lw=1.5, zorder=5, solid_capstyle="butt")
    _legend(ax, [
        Line2D([], [], marker="o", ls="none", ms=4.2, mfc=C_WT, mec="none",
               label="Within WT"),
        Line2D([], [], marker="o", ls="none", ms=4.2, mfc=C_IVT, mec="none",
               label="Within IVT"),
        Line2D([], [], marker="o", ls="none", ms=4.2, mfc=C_CROSS, mec="none",
               label="WT \u00d7 IVT"),
        Line2D([], [], color=C_DARK, lw=1.5, label="Pair mean"),
    ], loc="upper right", ncol=4)
    _LAYOUT.append({"panel": "B", "n_pair_points": int(len(pairs)),
                    "n_zero_pairs": n_zero})
    return ax


# --------------------------------------------------------------------------- #
# panel C -- score location per unit (where the calls sit in the tool's range)
# --------------------------------------------------------------------------- #
def panel_c(fig, top_in: float, height_in: float, loc: pd.DataFrame) -> plt.Axes:
    """Where each unit's calls sit inside the tool's own 0-1 score.

    Fifth version of this panel (2026-09-21).  The six density facets of the
    third version plotted three different score kinds side by side and only
    became legible after a per-facet rescaling; the AUC strip of the fourth
    version was thin (the pooled per-tool AUC is Fig. 7F) and read as "the
    benchmark cannot separate the conditions" rather than as the mechanism R3-9
    asks about.  This panel shows the mechanism directly: the point is one
    unit's median score and the whisker its interquartile range, WT (filled)
    against unmodified IVT (open), so all four signal-level reasons are visible
    in one strip -- the NanoMUD pair pinned at 1.000 with a zero-width range,
    NanoPsu and NanoSPA-Psi in a 0.955-0.980 band in *both* conditions, NanoNm
    low and overlapping, and CHEUI-m5C scoring the unmodified library at or
    above the wild type.

    Each tool is drawn on its own native 0-1 score (the modification ratio for
    CHEUI-m5C and NanoNm, the reported probability for the other four); the
    caption states that, and it is the reason this panel is not read as a
    separation statistic.
    """
    ax = fig.add_axes(rect(top_in, height_in))
    tool_axis(ax, labels=False)      # names are printed under panel D only
    # 1.05 keeps the markers of the saturated tools (score exactly 1.000) inside
    # the axes instead of clipping their upper half
    ax.set_ylim(0.0, 1.05)
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_yticklabels(["0", "0.5", "1"])
    ax.set_ylabel("Score", y=0.52, labelpad=7.0)

    #: drawing area -- C had no key at all (filled/open were caption-only)
    _legend(ax, [
        Line2D([], [], marker="o", ls="none", ms=4.2, mfc=C_WT, mec=C_WT,
               label="Human WT"),
        Line2D([], [], marker="o", ls="none", ms=4.2, mfc="white", mec=C_IVT,
               label="unmodified IVT"),
    ], loc="center right")
    n_marks = 0
    saturated: list[str] = []
    medians: dict[str, dict[str, float]] = {}
    for i, tool in enumerate(TOOL_ORDER):
        x = float(i)
        d = loc[loc.tool == tool]
        medians[tool] = {}
        for cond, off, col, filled in (("WT", -0.16, C_WT, True),
                                       ("IVT", 0.16, C_IVT, False)):
            sub = d[d.condition == cond].sort_values("sample")
            offs = jitter(len(sub), 0.05)
            for u, r in enumerate(sub.itertuples()):
                mid, lo, hi = float(r.median), float(r.q25), float(r.q75)
                ax.errorbar([x + off + offs[u]], [mid],
                            yerr=[[mid - lo], [hi - mid]],
                            fmt="o", ms=4.0, mfc=col if filled else "white",
                            mec=col, mew=0.9, ecolor=col, elinewidth=1.0,
                            capsize=2.2, zorder=4)
                n_marks += 1
                medians[tool][str(r.sample)] = mid
        if float(d["frac_ge_0p9"].min()) >= 1.0:
            saturated.append(tool)
    _LAYOUT.append({"panel": "C", "kind": "score_location_strip",
                    "n_marks": n_marks, "per_tool": len(TOOL_ORDER),
                    "marks_per_tool": 6,
                    "y_basis": "each tool's own 0-1 score",
                    "y_limits": list(ax.get_ylim()),
                    "saturated_tools": saturated,
                    "unit_medians": medians})
    return ax


# --------------------------------------------------------------------------- #
# panel D -- HeLa calls per unit + unmodified-control ratio
# --------------------------------------------------------------------------- #
def panel_d(fig, top_in: float, counts_h: float, strip_gap: float,
            strip_h: float, rep: pd.DataFrame, summ: pd.DataFrame,
            ratio: pd.DataFrame) -> tuple[plt.Axes, plt.Axes]:
    counts = fig.add_axes(rect(top_in, counts_h))
    tool_axis(counts, labels=False)
    counts.set_yscale("log")
    counts.set_ylim(80.0, 5.0e5)
    counts.set_yticks([1e2, 1e3, 1e4, 1e5])
    counts.set_yticklabels(["10$^{2}$", "10$^{3}$", "10$^{4}$", "10$^{5}$"])
    counts.minorticks_off()
    _LOG_AXES.append({"panel": "D", "yscale": counts.get_yscale(),
                      "ytick_labels": [t.get_text()
                                       for t in counts.get_yticklabels()]})
    counts.set_ylabel("Calls per unit", y=0.45, labelpad=7.0)
    for i, tool in enumerate(TOOL_ORDER):
        x = float(i)
        for cond, off, col, filled in (("WT", -0.12, C_WT, True),
                                       ("IVT", 0.12, C_IVT, False)):
            sub = rep[(rep.tool == tool) & (rep.condition == cond)]
            vals = sub.sort_values("sample")["n_calls_in_universe"].to_numpy(float)
            counts.plot(x + off + jitter(len(vals), 0.045), vals, marker="o",
                        ls="none", ms=4.2, mfc=col if filled else "white",
                        mec=col, mew=0.9, zorder=4)
        s = summ[summ.tool == tool]
        for cond, col, a, b in (("WT", C_WT, -0.34, -0.02),
                                ("IVT", C_IVT, 0.02, 0.34)):
            u = float(s.loc[s.condition == cond, "union_in_universe"].iloc[0])
            counts.plot([x + a, x + b], [u, u], color=col, lw=1.3,
                        ls=(0, (3, 2)), zorder=3)

    #: (filled WT, open IVT, dashed union per condition) instead of pointing at
    #: panel B's key; the free upper-right corner takes it without touching data
    _legend(counts, [
        Line2D([], [], marker="o", ls="none", ms=4.2, mfc=C_WT, mec=C_WT,
               label="Human WT"),
        Line2D([], [], marker="o", ls="none", ms=4.2, mfc="white", mec=C_IVT,
               label="unmodified IVT"),
        Line2D([], [], color=C_DARK, lw=1.3, ls=(0, (3, 2)), label="union"),
    ], loc="upper right")

    strip = fig.add_axes(rect(top_in + counts_h + strip_gap, strip_h))
    tool_axis(strip)
    strip.set_ylim(0.35, 3.35)
    strip.set_yticks([1.0, 2.0, 3.0])
    strip.set_yticklabels(["1", "2", "3"])
    strip.axhline(1.0, color=C_GREY, lw=1.0, ls=(0, (2, 2)), zorder=1)
    # the label sits left of the strip's own tick labels *and* clear of the
    # counts axis above it (the renderer re-places axis labels on the second
    # draw, so the clearance has to be explicit rather than incidental)
    strip.set_ylabel("Ratio", y=0.5, labelpad=12.0)
    #: 2026-09-28 (union drop): the compression strip carries only the
    #: mean-of-counts ratio (open circle + whiskers, its unit-level bootstrap
    #: CI); the union-ratio diamond was retired with the rest of the pooled
    #: reporting units, so the key has a single entry
    _legend(strip, [
        Line2D([], [], marker="o", ls="none", ms=4.2, mfc="white", mec=C_DARK,
               mew=1.0, label="mean-of-counts ratio"),
    ], loc="upper left", ncol=1)
    for i, tool in enumerate(TOOL_ORDER):
        r = ratio[ratio.tool == tool].iloc[0]
        mean = float(r.ratio_mean_counts)
        strip.errorbar([float(i)], [mean],
                       yerr=[[mean - float(r.ci_lo)], [float(r.ci_hi) - mean]],
                       fmt="o", ms=4.2, mfc="white", mec=C_DARK, mew=1.0,
                       ecolor=C_CROSS, elinewidth=1.0, capsize=2.2, zorder=4)
    _LAYOUT.append({"panel": "D",
                    "counts_ylim": list(counts.get_ylim()),
                    "ratio_ylim": list(strip.get_ylim())})
    return counts, strip


# --------------------------------------------------------------------------- #
# data
# --------------------------------------------------------------------------- #
def load() -> dict[str, pd.DataFrame]:
    return {
        "cc": pd.read_csv((_RB / "figures/figureS7/tables/s6_curlcake_per_construct.tsv"), sep="\t"),
        "pairs": pd.read_csv((_RB / "figures/figureS7/tables/s6_jaccard_pairs.tsv"), sep="\t"),
        "loc": pd.read_csv((_RB / "figures/figureS7/tables/s6_score_location_per_unit.tsv"), sep="\t"),
        "rep": pd.read_csv((_RB / "figures/figureS7/tables/s6_counts_per_replicate.tsv"), sep="\t"),
        "summ": pd.read_csv((_RB / "figures/figureS7/tables/s6_counts_summary.tsv"), sep="\t"),
        "ratio": pd.read_csv((_RB / "figures/figureS7/tables/s6_ratio_ci.tsv"), sep="\t"),
    }


# --------------------------------------------------------------------------- #
# legend-vs-data audit (a legend may never cover a data point)
# --------------------------------------------------------------------------- #
def _data_points(ax) -> np.ndarray:
    pts = []
    for ln in ax.lines:
        xy = np.asarray(ln.get_xydata(), dtype=float)
        if xy.size:
            pts.append(xy)
    for coll in ax.collections:
        segs = getattr(coll, "get_segments", lambda: [])()
        arr = [np.asarray(s, dtype=float) for s in segs if len(s)]
        if arr:
            pts.append(np.vstack(arr))
    return np.vstack(pts) if pts else np.zeros((0, 2))


def legend_bounds_audit(fig) -> list[dict]:
    """Every legend box must sit inside its own axes.

    A legend wider than its panel is clipped by the canvas edge without any
    error, so it is measured explicitly (the A key overflowed by 0.4 in on the
    first compressed render and lost its last entry).
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    out: list[dict] = []
    for ax in fig.axes:
        leg = ax.get_legend()
        if leg is None or not ax.axison:
            continue
        lb = leg.get_window_extent(renderer)
        ab = ax.get_window_extent(renderer)
        slack = 1.0                                    # px
        if (lb.x0 < ab.x0 - slack or lb.x1 > ab.x1 + slack
                or lb.y0 < ab.y0 - slack or lb.y1 > ab.y1 + slack):
            out.append({"panel": f"axes{fig.axes.index(ax)}",
                        "legend_w_in": round(lb.width / fig.dpi, 3),
                        "panel_w_in": round(ab.width / fig.dpi, 3),
                        "overflow_in": round(
                            max(lb.x1 - ab.x1, ab.x0 - lb.x0,
                                lb.y1 - ab.y1, ab.y0 - lb.y0) / fig.dpi, 3)})
    return out


def legend_data_audit(fig) -> list[dict]:
    """Data points inside a legend box, tested in data coordinates.

    Bbox intersection is *not* usable here (a long dashed union line or a
    density curve would produce false positives); the user's rule of
    2026-09-21 is that the marks themselves must be clear of the legend.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    hits: list[dict] = []
    for ax in fig.axes:
        leg = ax.get_legend()
        if leg is None or not ax.axison:
            continue
        box = leg.get_window_extent(renderer)
        inv = ax.transData.inverted()
        corners = inv.transform([(box.x0, box.y0), (box.x1, box.y1)])
        x0, x1 = sorted(corners[:, 0])
        y0, y1 = sorted(corners[:, 1])
        pts = _data_points(ax)
        if not len(pts):
            continue
        inside = ((pts[:, 0] >= x0) & (pts[:, 0] <= x1)
                  & (pts[:, 1] >= y0) & (pts[:, 1] <= y1))
        if inside.any():
            hits.append({"panel": ax.get_title() or f"axes{fig.axes.index(ax)}",
                         "legend": [t.get_text() for t in leg.get_texts()][:2],
                         "n_points_inside": int(inside.sum())})
    return hits


# --------------------------------------------------------------------------- #
def build(data: dict[str, pd.DataFrame]) -> plt.Figure:
    global _SLACK
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))

    # only the bottom-most axis of the B/C/D block carries the tool names
    heights = {"A": PLOT_H["A"] + XLAB_H, "B": PLOT_H["B"] + 0.06,
               "C": C_STRIP_H + 0.06,
               "D": PLOT_H["D"] + STRIP_GAP + STRIP_H + XLAB_H}
    tops: dict[str, float] = {}
    cur = TOP
    for letter in "ABCD":
        tops[letter] = cur
        cur += heights[letter] + ROW_GAP
    _SLACK = (PAGE_H - BOTTOM) - (cur - ROW_GAP)
    if _SLACK < 0:
        raise SystemExit(f"layout overflows the page by {-_SLACK:.3f} in")

    for letter in "ABCD":
        pl.margin_letter(fig, letter, y=1.0 - tops[letter] / PAGE_H,
                         fontsize=STYLE["letter"])

    panel_a(fig, tops["A"], PLOT_H["A"], data["cc"])
    panel_b(fig, tops["B"], PLOT_H["B"], data["pairs"])
    panel_c(fig, tops["C"], C_STRIP_H, data["loc"])
    panel_d(fig, tops["D"], PLOT_H["D"], STRIP_GAP, STRIP_H,
            data["rep"], data["summ"], data["ratio"])
    return fig


# --------------------------------------------------------------------------- #
def print_preview(src: Path, dst: Path, width_mm: float = 169.0) -> None:
    """Down-scale to the printed width for a human check (300 dpi)."""
    from PIL import Image
    im = Image.open(src)
    target = int(round(width_mm / 25.4 * 300))
    im.resize((target, int(round(im.height * target / im.width))),
              Image.LANCZOS).save(dst, dpi=(300, 300))
    print(f"[print] {dst} ({width_mm:.0f} mm wide)", flush=True)


def main() -> None:
    t0 = time.time()
    for d in (FIG, LOG):
        d.mkdir(parents=True, exist_ok=True)
    apply_style()
    data = load()
    missing = [t for t in TOOL_ORDER
               if t not in set(data["ratio"].tool) | set(data["loc"].tool)]
    if missing:
        raise SystemExit(f"tables are missing tools: {missing}")
    fig = build(data)

    pl.assert_page_clean(fig, min_pt=7.0)

    _LEGEND_HITS.extend(legend_data_audit(fig))
    for hit in _LEGEND_HITS:
        print(f"[legend] data under legend: {hit}", flush=True)
    _LEGEND_OUT.extend(legend_bounds_audit(fig))
    for hit in _LEGEND_OUT:
        print(f"[legend] legend overflows its panel: {hit}", flush=True)

    fig.savefig((_RB / "figures/figureS7/figures/FigureS6_rev.pdf"))
    fig.savefig((_RB / "figures/figureS7/figures/FigureS6_rev.png"), dpi=300)
    fig.savefig((_RB / "figures/figureS7/figures/FigureS6_rev_600dpi.pdf"), dpi=600)
    fig.savefig((_RB / "figures/figureS7/figures/FigureS6_rev_600dpi.png"), dpi=600)
    print_preview((_RB / "figures/figureS7/figures/FigureS6_rev.png"), (_RB / "figures/figureS7/figures/FigureS6_print_preview.png"))

    texts = pl._text_artists(fig)
    n_titles = sum(1 for ax in fig.axes if ax.get_title().strip())
    report = {"canvas_mm": [round(PAGE_W * 25.4, 1), round(PAGE_H * 25.4, 1)],
              "page_in": [PAGE_W, PAGE_H],
              "page_pt": [round(PAGE_W * 72, 3), round(PAGE_H * 72, 3)],
              "style_pt": STYLE,
              "font_page_ratio_pct": round(STYLE["tick"] / (PAGE_W * 72) * 100, 2),
              "panels": _LAYOUT,
              "log_axes": _LOG_AXES,
              "legend_data_violations": _LEGEND_HITS,
              "legend_out_of_panel": _LEGEND_OUT,
              "n_titles": n_titles,
              "slack_in": round(_SLACK, 3),
              "n_text": len(texts),
              "min_font_pt": round(min(t.get_fontsize() for _o, _a, t in texts), 2)}
    ((_RB / "figures/figureS7/logs/54_layout.json")).write_text(json.dumps(report, indent=2))
    plt.close(fig)
    print(f"[done] FigureS6_rev.pdf/png in {time.time() - t0:.1f} s", flush=True)
    print(f"[log] {LOG / '54_layout.json'}", flush=True)
    if _LEGEND_HITS:
        raise SystemExit("a legend covers data -- move it before shipping")
    if _LEGEND_OUT:
        raise SystemExit("a legend overflows its panel -- shorten it first")


if __name__ == "__main__":
    main()
