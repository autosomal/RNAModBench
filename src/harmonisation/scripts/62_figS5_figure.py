#!/usr/bin/env python
"""62 -- Figure S5 (revised): the detail the main Figure 6 does not carry.

The rebuilt main Figure 6 (``65_fig6_combination_page.py``) already shows the
selection criterion, the selected trade-off, the per-unit trajectories, the
configurations themselves, the coverage of the site sets and their control
burden.  This supplementary figure is therefore *complementary by
construction*: it draws only what the main figure cannot hold, and no quantity
appears twice.

A  the full greedy path 1 -> 13 tools (the main figure stops at k = 5): union
   recall with the best-single-tool reference, union PPV, and the intersection
   -- both its PPV and the recall that intersecting costs (0.1-2.2% at k = 5,
   exactly 0 once six or more tools are combined).  The shaded band marks the
   k range the main figure covers;
B  per-tool site quality: coverage (log) of each configuration's own calls;
C  per-tool site quality: DRACH fraction of the same calls;
D  the control burden per tool: FP/10 kb of every configuration on Curlcake IVT
   (filled) and HeLa IVT (open), one dot per sample.  The main figure shows the
   controls only for the five selected *sets*; this is the per-configuration
   view, and because controls are not species-specific it carries no group
   dimension;
E  the DRACH composition of the site sets behind the trade-off (the best single
   tool, what the union adds, both intersections and the GLORI reference).

Rows B, C and D are three full-width tiers *sharing one x axis*: the 13 tool
configurations, whose names are printed once, **slanted 45 deg** under row D
(see ``TOOL_ROT``): neighbouring names are parallel lines 0.47 in apart, so no
name can collide with its neighbour however long it is.  Inside every position
the four independent groups sit as separate clusters (Arabidopsis | mouse study
A | mouse study B | HeLa) -- dots per unit plus the mean within the group, never
a mean across groups and never a mean across the two mouse studies.

The selection criterion itself is deliberately *not* re-drawn: it is the
subject of main-text Fig. 6D and of the frozen table ``figS5_search_space.tsv``.
The paste-ready legend states the two facts a reader needs there -- PPV >= p0
held for every enumerated combination of every k (p0 is a formal guardrail, so
the criterion is "maximum mean recall at each k"), and the greedy forward
selection reproduces every exhaustive optimum.

Layout (2026-09-21, user decisions): a compact **240 x 175 mm page** (1:1, no
rescaling) so the text reads large relative to the figure; row A keeps three
species columns (Arabidopsis | Mouse, both studies in one axes, solid study A /
dashed study B | HeLa) because it plots curves against k, while rows B-E are
single full-width panels.  The only slanted text in the figure is the block of
13 tool names under row D (45 deg, parallel to one another); the k values of row
A and the five set names of row E stay horizontal.  Panel letters A-E are
printed once per row.

House rules: drawn at the final printed size, >= 9 pt printed everywhere,
Arial, no gridlines, no in-panel annotation text, bold panel letters, vector PDF
+ 300 dpi PNG.

Outputs -> figures/figureS6/figures/FigureS5_rev
(the frozen evidence tables stay in ``figures/figure6/tables``)

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    src/harmonisation/scripts/62_figS5_figure.py
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
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.transforms import Bbox

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from common import config as C                       # noqa: E402
from common.figstyle import apply as figstyle_apply  # noqa: E402
from common.manifest import setup_logger             # noqa: E402
#: rotated text must be judged by its quad, not its (much larger) axis-aligned
#: bounding box -- the shared gate works on the ink rectangle of the rotated box
from common.pagelayout import (quad_gap, quads_intersect,  # noqa: E402
                               rect_quad, text_quad)

#: S5 lives next to the other supplementary revisions; its evidence tables stay
#: in the Fig. 6 family directory because ``40_fig6_combination.py`` and
#: ``61_figS5_sitequality.py`` produce them and the rebuilt main Figure 6 reads
#: the same frozen tables.
OUT = (_RB / "figures/figureS6")
FROZEN = (_RB / "figures/figure6")
TAB, FIG, LOG = FROZEN / "tables", OUT / "figures", OUT / "logs"

# ------------------------------------------------------------------ canvas --- #
#: 240 x 175 mm (9.4488 x 6.8898 in), printed 1:1.  The page is deliberately
#: smaller than A4: the text has to read large *relative to the figure*, and the
#: 13 tool names are staggered over two lines so the width can come down too.
#: 2026-09-23: 9.4488 in (240 mm) was wider than the SI text column (178 mm),
#: so the page was placed at 0.74 and the 9 pt labels printed at 6.7 pt.  The
#: canvas is now 178 mm wide and prints 1:1.
CANVAS_W, CANVAS_H = 7.0079, 6.8898
MIN_PT = 9.0
#: (letter, panel height, gap below the panel) in inches, top down.  Row A is
#: the three curve columns; rows B-D share one tool axis (only D carries the
#: names, hence its larger gap); row E carries the five set names.
BANDS = [("A", 1.14, 0.40), ("B", 0.76, 0.15), ("C", 0.58, 0.15),
         ("D", 0.74, 0.90), ("E", 0.74, 0.38)]
TOP, HEADER, LEGEND_H = 0.02, 0.46, 0.42
#: the left margin carries the panel letters *and* the rotated y-axis labels of
#: four rows, now at 10.5 pt with 9.5 pt tick labels ("0.0001" in row D): the
#: widest tick text decides how far left the label sits, measured 0.92 in keeps
#: every label clear of the 12 pt panel letters (0.88 left them 1.6 px overlapped)
L_LEFT, L_GAP, L_RIGHT = 0.92, 0.26, 0.06
FS = {"tick": 9.5, "tool": 9.5, "label": 10.5, "title": 10.5, "legend": 9.0,
      "letter": 12.0}
#: marker for one independent unit
UNIT_S = 14.0
#: the 13 tool names are slanted 45 deg under row D: neighbouring names sit
#: 0.66 in apart along x, i.e. ~0.47 in apart perpendicular to their own
#: baseline, far more than the 0.13 in cap height -- so no name can collide with
#: its neighbour however long it is, and every name begins at its own tick
TOOL_ROT = 45.0

#: row A columns: the mouse column carries BOTH independent studies in one axes
#: (two shades, solid/dashed, filled/open), never averaged; single-group columns
#: have a one-element ``groups`` list so every panel can loop uniformly
COLS = [{"sp": "Arabidopsis", "groups": ["Arabidopsis_WT"], "colour": "#1E888B",
         "title": "Arabidopsis"},
        {"sp": "Mouse", "groups": ["studyA", "studyB"], "colour": "#C07A16",
         "title": "Mouse"},
        {"sp": "Human", "groups": ["HeLa_WT"], "colour": "#3778A0",
         "title": "HeLa"}]
#: per-study drawing style: the two mouse shades plus the visual channels the
#: panels use to keep the studies apart (line style, dot fill)
GROUP_STYLE = {"Arabidopsis_WT": {"colour": "#1E888B", "ls": "-", "fill": True},
               "studyA": {"colour": "#E8A33D", "ls": "-", "fill": True},
               "studyB": {"colour": "#9C5B12", "ls": (0, (4, 2)), "fill": False},
               "HeLa_WT": {"colour": "#3778A0", "ls": "-", "fill": True}}
#: the four independent groups of rows B / C / E, in the order of row A, each
#: group drawn as its own cluster inside one x position
GROUPS = ("Arabidopsis_WT", "studyA", "studyB", "HeLa_WT")
GROUP_LABEL = {"Arabidopsis_WT": "Arabidopsis", "studyA": "mouse study A",
               "studyB": "mouse study B", "HeLa_WT": "HeLa"}
#: cluster geometry per row type: x offsets of the four clusters inside one
#: position, the jitter of the dots within a cluster and the half-width of the
#: within-group mean bar (all in x data units, i.e. one tool / one set apart)
CLUSTER_TOOL = {"offset": (-0.27, -0.09, 0.09, 0.27), "spread": 0.040,
                "half": 0.055}
CLUSTER_SET = {"offset": (-0.22, -0.075, 0.075, 0.22), "spread": 0.045,
               "half": 0.085}
UNION_C, ISECT_C = "#1f5c8b", "#e08a2e"
CONTROLS = ("Curlcake IVT", "HeLa IVT")
#: the five sets of panel E: the best single tool, what combining adds, the two
#: intersections and the reference.  The five-tool union itself is deliberately
#: left out -- it is the *sum* of the single-tool set and the marginal set, and
#: the trade-off is read from its two components.
SETS = ["single (k=1)", "union marginal", "intersection (k=2)",
        "intersection (k=5)", "GLORI reference"]
#: horizontal, two lines, under panel E: the wide page leaves 2.2 in per set, so
#: the names never have to be rotated
SET_LABELS = ["single tool\n$k$ = 1", "union\nmarginal",
              "intersection\n$k$ = 2", "intersection\n$k$ = 5",
              "GLORI\nreference"]


def apply_style() -> None:
    figstyle_apply()
    mpl.rcParams.update({
        "font.size": FS["tick"],
        "axes.labelsize": FS["label"],
        "xtick.labelsize": FS["tick"],
        "ytick.labelsize": FS["tick"],
        "legend.fontsize": FS["legend"],
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 2.4, "ytick.major.size": 2.4,
        "legend.frameon": False,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "mathtext.fontset": "custom", "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
        "figure.dpi": 120,
    })


# ------------------------------------------------------------------- data ---- #
def load() -> dict[str, pd.DataFrame]:
    return {
        "space": pd.read_csv(TAB / "figS5_search_space.tsv", sep="\t"),
        "greedy": pd.read_csv(TAB / "figS5_greedy_1to13.tsv", sep="\t"),
        "members": pd.read_csv(TAB / "fig6_selected_members.tsv", sep="\t"),
        "qual": pd.read_csv(TAB / "figS5_site_quality.tsv", sep="\t"),
        "toolq": pd.read_csv(TAB / "figS5_tool_quality.tsv", sep="\t"),
        "ctltool": pd.read_csv(TAB / "fig6_negative_control_fp_bytool.tsv",
                               sep="\t"),
    }


def tool_order(members: pd.DataFrame, space: pd.DataFrame) -> list[str]:
    """Display order for the per-tool rows B / C / D.

    Earliest k at which a tool is selected in *any* study, then the best
    single-tool recall it reaches in *any* study (never a mean across the two
    independent mouse studies), then the name.
    """
    first: dict[str, int] = {}
    for r in members[members.in_combination].itertuples():
        first[r.tool] = min(first.get(r.tool, 99), int(r.k))
    k1 = space[space.k == 1]
    recall = k1.assign(tool=k1.combination).groupby("tool").union_recall_mean.max()
    return sorted(set(members.tool),
                  key=lambda t: (first.get(t, 99), -recall.get(t, 0.0), t))


def _offsets(n: int, half: float = 0.16) -> np.ndarray:
    return np.linspace(-half, half, n) if n > 1 else np.zeros(1)


def _decades(lo: float, hi: float, max_n: int = 5) -> list[float]:
    ticks = [10.0 ** e for e in range(-6, 7) if lo <= 10.0 ** e <= hi]
    while len(ticks) > max_n:
        ticks = ticks[::2]
    return ticks


# ----------------------------------------------------------------- layout ---- #
def build_stack() -> tuple[float, dict[str, list[float]], dict[str, float]]:
    """Bands from the top of an **exact A4 landscape page** (print contract).

    Any slack stays between the last panel and the legend, so the page size is
    always 841.9 x 595.3 pt and nothing is rescaled on the way into the SI.
    """
    fig_h = CANVAS_H
    used = TOP + HEADER + LEGEND_H + sum(h + g for _, h, g in BANDS)
    assert used <= fig_h + 1e-6, \
        f"bands need {used:.2f} in but the A4 landscape canvas is {fig_h:.2f} in"
    rects: dict[str, list[float]] = {}
    letters: dict[str, float] = {}
    y = fig_h - TOP - HEADER
    for letter, h, gap in BANDS:
        letters[letter] = y
        rects[letter] = [y - h, y]
        y -= h + gap
    return fig_h, rects, letters


def _rect(x_in: float, y_in: float, w_in: float, h_in: float,
          fig_w: float, fig_h: float) -> list[float]:
    return [x_in / fig_w, y_in / fig_h, w_in / fig_w, h_in / fig_h]


def col_axes(fig, rects, letter, fig_w, fig_h) -> list:
    """Row A: one axes per species column (the curves run against k)."""
    y0, y1 = rects[letter]
    w = (CANVAS_W - L_LEFT - L_RIGHT - (len(COLS) - 1) * L_GAP) / len(COLS)
    axes = []
    for i, _ in enumerate(COLS):
        ax = fig.add_axes(_rect(L_LEFT + i * (w + L_GAP), y0, w, y1 - y0,
                                fig_w, fig_h))
        ax.tick_params(length=2.2, pad=1.4)
        axes.append(ax)
    return axes


def full_axes(fig, rects, letter, fig_h):
    """Rows B-E: one panel across the full text width.

    The tool axis of B / C / D is shared (same x limits, labels only under D),
    and row D has no group dimension at all (controls are not species-specific),
    so repeating anything in columns would suggest measurements that do not
    exist.
    """
    y0, y1 = rects[letter]
    ax = fig.add_axes(_rect(L_LEFT, y0, CANVAS_W - L_LEFT - L_RIGHT, y1 - y0,
                            CANVAS_W, fig_h))
    ax.tick_params(length=2.2, pad=1.4)
    return ax


# ----------------------------------------------------------------- panels ---- #
def panel_a(ax, col: dict, greedy: pd.DataFrame) -> None:
    """The full greedy path 1 -> 13 tools: the union and the cost of intersecting.

    Nothing here repeats the main figure: it stops at k = 5 (the shaded band),
    while the path beyond the five selected tools is only visible here.  Mouse:
    one curve per study (solid = A, dashed = B), never averaged.
    """
    ax.axvspan(0.55, 5.0, color="0.93", zorder=0)
    for g in col["groups"]:
        d = greedy[(greedy.species == col["sp"])
                   & (greedy.group == g)].sort_values("k")
        ls = GROUP_STYLE[g]["ls"]
        ax.plot(d.k, 100 * d.union_recall_mean, color=UNION_C, lw=1.3, ls=ls,
                marker="o", ms=2.8, zorder=4)
        ax.plot(d.k, 100 * d.union_precision_mean, color=UNION_C, lw=1.1,
                ls=ls, marker="o", ms=2.8, mfc="white", mec=UNION_C, zorder=3)
        ax.plot(d.k, 100 * d.isect_precision_mean, color=ISECT_C, lw=1.1,
                ls=ls, marker="s", ms=2.6, mfc="white", mec=ISECT_C, zorder=3)
        ax.plot(d.k, 100 * d.isect_recall_mean, color=ISECT_C, lw=1.1, ls=ls,
                marker="^", ms=2.8, mfc="white", mec=ISECT_C, zorder=3)
        ref = 100 * float(d[d.k == 1].union_recall_mean.iloc[0])
        ax.axhline(ref, color="0.55", lw=0.9, ls=(0, (1.6, 1.6)), zorder=1)
    ax.set_xticks([1, 3, 5, 7, 9, 11, 13])
    ax.set_xlim(0.55, 13.45)
    ax.set_ylim(-2, 105)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlabel("tools combined", fontsize=FS["label"])
    # the column headers sit above row A; HEADER in build_stack() reserves them
    ax.set_title(col["title"], fontsize=FS["title"], fontweight="bold", pad=3)


def _per_unit_dots(ax, sub: pd.DataFrame, x: float, values, colour: str, *,
                   half: float = 0.34, filled: bool = True,
                   spread: float = 0.13) -> None:
    """Dots (one per unit) plus that group's mean bar, at one x position.

    ``filled`` keeps the two mouse studies apart (filled = study A, open =
    study B); the bar is always the mean *within* one group.
    """
    for off, (_, r) in zip(_offsets(len(sub), spread), sub.iterrows()):
        v = values(r)
        if np.isfinite(v):
            ax.scatter([x + off], [v], s=UNIT_S, color=colour,
                       facecolors=colour if filled else "white",
                       edgecolors=colour, linewidths=0.0 if filled else 0.7,
                       zorder=3)
    vals = np.asarray([values(r) for _, r in sub.iterrows()], dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size:
        ax.plot([x - half, x + half], [vals.mean()] * 2, color="0.25", lw=1.7,
                zorder=4, solid_capstyle="butt")


def _group_cluster(ax, x: float, per_group: dict[str, pd.DataFrame], values,
                   geom: dict) -> None:
    """The four independent groups as four clusters inside one x position."""
    for g, off in zip(GROUPS, geom["offset"]):
        sub = per_group.get(g)
        if sub is None or not len(sub):
            continue
        _per_unit_dots(ax, sub, x + off, values, GROUP_STYLE[g]["colour"],
                       half=geom["half"], filled=GROUP_STYLE[g]["fill"],
                       spread=geom["spread"])


def panel_b(ax, toolq: pd.DataFrame, tools: list[str],
            xpos: dict[str, float]) -> None:
    """Coverage of every configuration's own calls (log), four groups per tool."""
    for t in tools:
        sub = toolq[toolq.tool == t]
        _group_cluster(ax, xpos[t],
                       {g: sub[sub.group == g] for g in GROUPS},
                       lambda r: r.coverage_median, CLUSTER_TOOL)
    ax.set_yscale("log")
    vals = pd.to_numeric(toolq.coverage_median, errors="coerce").dropna()
    assert len(vals), "per-tool coverage table is empty"
    lo = max(1.0, 10 ** np.floor(np.log10(max(float(vals.min()), 1))))
    hi = 10 ** np.ceil(np.log10(float(vals.max()) + 1e-9))
    ax.set_ylim(lo, hi)
    ticks = _decades(lo, hi, 4)
    ax.set_yticks(ticks)
    # 10^n, not 0.0001: row D's axis is the one the reader compares to the
    # tabulated FP-per-10 kb values
    ax.set_yticklabels([f"$10^{{{int(round(np.log10(t)))}}}$" for t in ticks])
    ax.minorticks_off()
    # two lines: at 10.5 pt the single-line label is longer than the panel is
    # tall and would run into the DRACH label of row C below it
    ax.set_ylabel("coverage\n(reads)", fontsize=FS["label"])


def panel_c(ax, toolq: pd.DataFrame, tools: list[str],
            xpos: dict[str, float]) -> None:
    """DRACH fraction of every configuration's own calls, four groups per tool."""
    for t in tools:
        sub = toolq[toolq.tool == t]
        _group_cluster(ax, xpos[t],
                       {g: sub[sub.group == g] for g in GROUPS},
                       lambda r: 100 * r.frac_drach, CLUSTER_TOOL)
    ax.set_ylim(0, 105)
    ax.set_yticks([0, 50, 100])
    ax.set_ylabel("DRACH (%)", fontsize=FS["label"])


def panel_d(ax, ctltool: pd.DataFrame, tools: list[str],
            xpos: dict[str, float]) -> None:
    """Control burden per tool: Curlcake IVT (filled) and HeLa IVT (open).

    One dot per sample; the Curlcake values exist for 10 of the 13
    configurations (the three HeLa-only ones have no Curlcake run), which the
    paste-ready legend states.  This row carries the horizontal tool names of
    the whole B-D block.
    """
    for i, cname in enumerate(CONTROLS):
        s = ctltool[ctltool.control == cname]
        for t in tools:
            st = s[s.tool == t]
            if not len(st):
                continue
            base = xpos[t] + (-0.15 if i == 0 else 0.15)
            _per_unit_dots(ax, st, base, lambda r: r.fp_per_10kb,
                           "0.45", half=0.11, filled=i == 0, spread=0.06)
    vals = ctltool.fp_per_10kb.dropna()
    vals = vals[vals > 0]
    lo = max(1e-5, 10 ** np.floor(np.log10(vals.min()) - 0.05))
    hi = 10 ** np.ceil(np.log10(vals.max() + 1e-9) + 0.05)
    ax.set_yscale("log")
    ax.set_ylim(lo, hi)
    ticks = _decades(lo, hi, 4)
    ax.set_yticks(ticks)
    # 10^n, not 0.0001: row D's axis is the one the reader compares to the
    # tabulated FP-per-10 kb values
    ax.set_yticklabels([f"$10^{{{int(round(np.log10(t)))}}}$" for t in ticks])
    ax.minorticks_off()
    ax.set_ylabel("FP / 10 kb", fontsize=FS["label"])


def panel_e(ax, qual: pd.DataFrame) -> None:
    """DRACH fraction of the five site sets, four groups per set."""
    for i, name in enumerate(SETS):
        sub = qual[qual.set == name]
        _group_cluster(ax, float(i), {g: sub[sub.group == g] for g in GROUPS},
                       lambda r: 100 * r.frac_drach, CLUSTER_SET)
    ax.set_xlim(-0.55, len(SETS) - 0.45)
    ax.set_ylim(0, 105)
    ax.set_yticks([0, 50, 100])
    ax.set_ylabel("DRACH (%)", fontsize=FS["label"])


# ----------------------------------------------------------------- legend ---- #
#: which entries belong to which key line.  All fourteen used to sit in one
#: five-column strip at the foot of the page, so the reader had to jump the whole
#: canvas to decode row A.
#: "shaded = main-text Fig. 6" left the figure for the caption: a 178 mm canvas
#: with the column titles on top of row A fits five entries in one line
KEY_A = ["union recall", "union PPV", "intersection recall", "intersection PPV",
         "best single tool"]
KEY_REST = ["one dot = one unit", "group mean", "Arabidopsis", "mouse study A",
            "mouse study B", "HeLa", "Curlcake IVT", "HeLa IVT"]


def handles() -> list:
    """Every legend entry of the page (``KEY_A`` / ``KEY_REST`` pick the rows)."""
    return [
        Line2D([], [], color=UNION_C, marker="o", ms=3.6, lw=1.2,
               label="union recall"),
        Line2D([], [], color=UNION_C, marker="o", ms=3.6, mfc="white", mew=0.9,
               lw=1.1, label="union PPV"),
        Line2D([], [], color=ISECT_C, marker="^", ms=3.6, mfc="white", mew=0.9,
               lw=1.1, label="intersection recall"),
        Line2D([], [], color=ISECT_C, marker="s", ms=3.4, mfc="white", mew=0.9,
               lw=1.1, label="intersection PPV"),
        Line2D([], [], color="0.55", lw=0.9, ls=(0, (1.6, 1.6)),
               label="best single tool"),
        Line2D([], [], color="0.93", lw=7.0,
               label="shaded = main-text Fig. 6"),
        Line2D([], [], color="0.25", marker="o", ms=3.0, ls="none",
               label="one dot = one unit"),
        Line2D([], [], color="0.25", lw=1.7, label="group mean"),
        Line2D([], [], color=GROUP_STYLE["Arabidopsis_WT"]["colour"],
               marker="o", ms=3.4, ls="none", label="Arabidopsis"),
        Line2D([], [], color=GROUP_STYLE["studyA"]["colour"], marker="o",
               ms=3.4, ls="none", label="mouse study A"),
        Line2D([], [], color=GROUP_STYLE["studyB"]["colour"], marker="o",
               ms=3.4, ls="none", mfc="white", mew=0.8,
               label="mouse study B"),
        Line2D([], [], color=GROUP_STYLE["HeLa_WT"]["colour"], marker="o",
               ms=3.4, ls="none", label="HeLa"),
        Line2D([], [], color="0.45", marker="o", ms=3.4, ls="none", mfc="0.45",
               mec="0.45", label="Curlcake IVT"),
        Line2D([], [], color="0.45", marker="o", ms=3.4, ls="none", mfc="white",
               mec="0.45", label="HeLa IVT"),
    ]


# ------------------------------------------------------------------ checks --- #
def layout_report(fig) -> dict:
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    items = []
    for t in fig.findobj(Text):
        if not t.get_visible() or not str(t.get_text()).strip():
            continue
        try:
            bb = t.get_window_extent(renderer=rend)
            quad = text_quad(t, rend)
        except Exception:                                      # pragma: no cover
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        if bb.width < 2.0 and bb.height < 2.0:
            continue          # detached tick labels sit at the figure origin
        xs = [p[0] for p in quad]
        ys = [p[1] for p in quad]
        items.append({"s": str(t.get_text())[:32], "fs": float(t.get_fontsize()),
                      "rot": float(t.get_rotation()),
                      "bb": [min(xs), min(ys), max(xs), max(ys)],
                      "quad": quad})
    # quad against quad, on the ink rectangle: a 45 deg label has a bounding box
    # several times its own area, so bbox tests would report collisions that are
    # never printed (and hide the ones that are)
    overlaps = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if quads_intersect(items[i]["quad"], items[j]["quad"], pad=2.0):
                overlaps.append((items[i]["s"], items[j]["s"]))
    fonts: dict[float, int] = {}
    for it in items:
        fonts[it["fs"]] = fonts.get(it["fs"], 0) + 1
    w_px = float(fig.get_size_inches()[0]) * float(fig.dpi)
    h_px = float(fig.get_size_inches()[1]) * float(fig.dpi)
    clipped = [(it["s"], round(it["bb"][0], 1), round(it["bb"][2], 1))
               for it in items
               if it["bb"][0] < -0.5 or it["bb"][2] > w_px + 0.5]
    # y-axis labels are vertical by definition; what must never be rotated are
    # the *category* labels on the x axes -- the tick labels (k values, set
    # names) and the staggered tool names, which are Text artists under row D
    category_texts = [t for ax in fig.axes for t in ax.get_xticklabels()]
    category_texts += [t for ax in fig.axes for t in ax.texts]
    rotated = [(str(t.get_text())[:24], float(t.get_rotation()))
               for t in category_texts
               if t.get_visible() and str(t.get_text()).strip()
               and abs((float(t.get_rotation()) + 180) % 180) > 1e-6]
    return {"n_text": len(items),
            "min_fontsize": min((it["fs"] for it in items), default=float("nan")),
            "fonts": dict(sorted(fonts.items())), "overlaps": overlaps[:20],
            # the full list: the 13 slanted tool names are checked by count
            "rotated": rotated, "clipped": clipped,
            "canvas_px": [w_px, h_px]}


def slant_clearance(fig, *, pad: float = 6.0) -> list[tuple[str, float]]:
    """Clearance (px) between every slanted tool name and row E underneath.

    The 45 deg names hang below row D; if the D->E gap is too small they dangle
    into row E's axes -- no text touches another text, so the collision gate
    stays silent while the two panels visually merge (found in the 2026-09-21
    review: the longest name reached 0.09 in *below* row E's top edge).
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    w = (CANVAS_W - L_LEFT - L_RIGHT) / CANVAS_W
    rows = sorted([a for a in fig.axes
                   if abs(a.get_position().width - w) < 1e-3],
                  key=lambda a: -a.get_position().y0)      # B, C, D, E
    d, e = rows[2], rows[3]
    e_rect = rect_quad(e.get_window_extent(rend))
    out = []
    for tl in d.get_xticklabels():
        if not tl.get_visible() or not str(tl.get_text()).strip():
            continue
        out.append((str(tl.get_text()),
                    float(quad_gap(text_quad(tl, rend), e_rect))))
    assert out, "row D carries no slanted tool names"
    return out


def check_tables(t: dict[str, pd.DataFrame], tools: list[str], logger) -> None:
    """Anchors: the figure may only draw what the frozen tables state."""
    assert len(tools) == 13, f"expected 13 configurations, got {len(tools)}"
    space = t["space"]
    # the paste-ready legend states that p0 is a formal guardrail; if that ever
    # stops holding, the sentence has to change with it
    assert bool(space.feasible.all()), \
        "some enumerated combination falls below p0 (the legend states otherwise)"
    for col in COLS:
        for g in col["groups"]:
            gp = t["greedy"][t["greedy"].group == g].sort_values("k")
            assert sorted(gp.k) == list(range(1, 14)), \
                f"{g}: greedy path incomplete"
            # the legend claims the greedy path reproduces every exhaustive
            # optimum -- compare the member *sets* (the two tables order the
            # members differently, the values agree to 1e-9)
            sel = space[(space.group == g) & space.selected]
            for k in range(1, 6):
                want = frozenset(sel[sel.k == k].combination.iloc[0].split("+"))
                got = frozenset(gp[gp.k == k].combination.iloc[0].split("+"))
                assert want == got, \
                    f"{g}: greedy k={k} is not the exhaustive optimum"
            q = t["toolq"][t["toolq"].group == g]
            assert set(q.tool) == set(tools), f"{g}: per-tool quality incomplete"
            for name in SETS:
                assert (t["qual"][t["qual"].group == g].set == name).sum() >= 1, \
                    f"{g}: site-quality set {name!r} missing"
    for cname in CONTROLS:
        s = t["ctltool"][t["ctltool"].control == cname]
        assert len(s) and set(s.tool) <= set(tools), f"{cname}: no per-tool rows"
    logger.info("table anchors OK: %d independent groups (the two mouse studies "
                "are kept separate and never averaged), 13-step greedy path per "
                "study (identical to the enumerated optima at k = 1..5), 13 "
                "tools and %d site sets per group, both controls",
                len(GROUPS), len(SETS))


# -------------------------------------------------------------------- main --- #
def build_figure(t: dict[str, pd.DataFrame],
                 tools: list[str]) -> tuple[plt.Figure, dict]:
    apply_style()
    xpos = {tool: float(i) for i, tool in enumerate(tools)}
    fig_h, rects, letters = build_stack()
    assert fig_h <= CANVAS_H, f"canvas {fig_h:.2f} in exceeds A4 ({CANVAS_H})"
    fig = plt.figure(figsize=(CANVAS_W, fig_h))

    ax = {letter: full_axes(fig, rects, letter, fig_h) for letter in "BCDE"}
    ax["A"] = col_axes(fig, rects, "A", CANVAS_W, fig_h)

    for j, col in enumerate(COLS):
        panel_a(ax["A"][j], col, t["greedy"])
    panel_b(ax["B"], t["toolq"], tools, xpos)
    panel_c(ax["C"], t["toolq"], tools, xpos)
    panel_d(ax["D"], t["ctltool"], tools, xpos)
    panel_e(ax["E"], t["qual"])

    # "coverage" only: at 10.5 pt the single-line label with its unit is taller
    # than the panel and would reach the DRACH label of row C (the unit is in
    # the caption, which states "coverage (reads, logarithmic)")
    ylabels = {"A": "recall / PPV (%)", "B": "coverage",
               "C": "DRACH (%)", "D": "FP / 10 kb", "E": "DRACH (%)"}
    ax["A"][0].set_ylabel(ylabels["A"], fontsize=FS["label"])
    for j in range(1, len(COLS)):
        ax["A"][j].set_ylabel("")
    for key in "BCDE":
        ax[key].set_ylabel(ylabels[key], fontsize=FS["label"])

    # rows B-D share one tool axis: ticks everywhere, the 13 names printed once,
    # slanted 45 deg under row D (no axis title: the names *are* the axis text);
    # row E carries its own horizontal set names
    tool_ticks = [xpos[t] for t in tools]
    for key in ("B", "C", "D"):
        ax[key].set_xlim(-0.75, len(tools) - 0.25)
        ax[key].set_xticks(tool_ticks)
        ax[key].set_xticklabels([])
    ax["D"].set_xticklabels(tools, rotation=TOOL_ROT, ha="right",
                            fontsize=FS["tool"])
    ax["D"].tick_params(pad=4.0)
    ax["E"].set_xticks(range(len(SETS)))
    ax["E"].set_xticklabels(SET_LABELS, rotation=0, ha="center",
                            fontsize=FS["tool"])

    for letter, y_in in letters.items():
        fig.text(0.010, (y_in - 0.005) / fig_h, letter, fontsize=FS["letter"],
                 fontweight="bold", va="top", ha="left")

    all_h = {h.get_label(): h for h in handles()}
    fig.legend(handles=[all_h[l] for l in KEY_A], loc="upper center",
               bbox_to_anchor=(0.5, 1 - 0.005 / fig_h), ncol=5, frameon=False,
               handlelength=1.4, columnspacing=1.2, handletextpad=0.4,
               labelspacing=0.30, borderpad=0.0)          # travels with row A
    fig.legend(handles=[all_h[l] for l in KEY_REST], loc="lower center",
               bbox_to_anchor=(0.5, 0.004), ncol=4, frameon=False,
               handlelength=1.4, columnspacing=1.2, handletextpad=0.4,
               labelspacing=0.30)
    return fig, {"fig_h": fig_h, "tools": tools}


def main() -> None:
    logger = setup_logger("62_figS5_figure", log_dir=LOG)
    apply_style()
    t = load()
    tools = tool_order(t["members"], t["space"])
    check_tables(t, tools, logger)
    fig, meta = build_figure(t, tools)
    rep = layout_report(fig)
    logger.info("canvas %.3f x %.3f in (%.0f x %.0f mm, printed 1:1)", CANVAS_W,
                meta["fig_h"], CANVAS_W * 25.4, meta["fig_h"] * 25.4)
    logger.info("tool order: %s", " > ".join(meta["tools"]))
    logger.info("text objects %d; font sizes %s; smallest %.2f pt; overlaps %d; "
                "rotated %d", rep["n_text"], rep["fonts"], rep["min_fontsize"],
                len(rep["overlaps"]), len(rep["rotated"]))
    for a, b in rep["overlaps"]:
        logger.warning("text overlap: %r / %r", a, b)
    for c in rep["clipped"]:
        logger.warning("text outside the canvas: %r", c)
    for s, rot in rep["rotated"]:
        if abs(rot - TOOL_ROT) > 1e-6 or s not in set(meta["tools"]):
            logger.warning("unexpected rotated text: %r at %.1f deg", s, rot)
    assert rep["min_fontsize"] >= MIN_PT, \
        f"smallest text {rep['min_fontsize']:.2f} pt < {MIN_PT} pt"
    assert not rep["overlaps"], f"{len(rep['overlaps'])} overlapping text pairs"
    assert not rep["clipped"], f"text outside the canvas: {rep['clipped'][:6]}"
    # slanted text is allowed for exactly one thing: the 13 tool names at
    # TOOL_ROT degrees.  Anything else rotated (or a name at another angle) is
    # a regression -- the k values of row A and the set names of row E must
    # stay horizontal.
    bad_rot = [(s, a) for s, a in rep["rotated"]
               if not (abs(a - TOOL_ROT) < 1e-6 and s in set(meta["tools"]))]
    assert len(rep["rotated"]) == len(meta["tools"]) and not bad_rot, \
        f"unexpected rotated text: {bad_rot[:4]} (only the {len(meta['tools'])} " \
        f"tool names at {TOOL_ROT:g} deg are allowed)"
    # the slanted names must stay clear of row E below them
    clear = slant_clearance(fig)
    worst = min(clear, key=lambda kv: kv[1])
    logger.info("slanted-name clearance to row E: %.1f px (tightest: %r)",
                worst[1], worst[0])
    assert worst[1] >= 6.0, \
        f"slanted tool name {worst[0]!r} dangles into row E ({worst[1]:.1f} px)"
    legs = [lg.get_window_extent(renderer=fig.canvas.get_renderer())
            for lg in fig.legends]
    for i, bb in enumerate(legs):
        assert float(bb.x0) >= -1.0 and float(bb.x1) <= rep["canvas_px"][0] + 1.0, \
            f"legend {i} does not fit the A4 width"
    logger.info("legends inside the canvas: %s",
                [[round(float(v), 1) for v in (b.x0, b.x1)] for b in legs])

    FIG.mkdir(parents=True, exist_ok=True)
    stem = FIG / "FigureS5_rev"
    # exact A4 landscape page (no tight bbox): the print size must equal the canvas
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches=None, pad_inches=0)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches=None,
                pad_inches=0)
    plt.close(fig)
    logger.info("wrote %s.pdf + %s.png", stem, stem)


if __name__ == "__main__":
    main()
