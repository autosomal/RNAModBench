#!/usr/bin/env python
"""65 -- main-text Figure 6, rebuilt at the printed size (2026-09-21).

Why this script exists
----------------------
The 2026-09-20 rebuild (``40_fig6_combination.py::fig_main``) was drawn on a
``figsize=(5.0*3, 3.4*4)`` canvas = 380 x 367 mm, while the manuscript places the
figure with ``\\includegraphics[width=0.95\\textwidth]{Figure6.pdf}`` = 479.5 pt =
169 mm, i.e. at a **0.445 scale**: its 8-14 pt text printed at 3.7-6.2 pt against
9-10 pt body text.  This script draws the evidence on a canvas whose width *is*
the printed width (0.95 \\textwidth = 6.66 in), so every font size written here is
the printed one, and each row answers one reviewer question.

Rows (three species columns: Arabidopsis | mouse study A/B | HeLa; the two mouse
studies stay separate series in every row and are never combined):

  A  the trade-off: union recall inside the measurable universe (filled), the
     same recall against **every** GLORI site (thin dashed -- the absolute
     coverage reviewer R3-4 asks about), union PPV (open) and intersection PPV
     (open squares) against k; error bars, SD across units (single-unit mouse
     studies carry none); thin dotted grey, best single tool;
  B  the same trade-off per independent unit over k = 1..5 -- union and
     intersection in one panel, the evidence that it is not an averaging
     artefact;
  C  the tools themselves, in two blocks sharing the 13 rows: **left**, what each
     configuration achieves on its own (recall and PPV vs. GLORI (2 bp), one dot
     per unit plus the mean; the 13 tool names are printed once, in the left
     margin of column 0 -- the three columns share those rows); **right**, the
     main effect of each configuration at fixed k (cell colour: change in mean
     recall when the tool is present, from the full enumeration) with the
     configurations the criterion actually selects marked on top (mouse study A
     light circles, study B dark squares);
  D  the selection criterion itself (reviewer R3-m2): every enumerated k = 1..5
     combination of the group in the (mean union recall, mean union PPV) plane,
     coloured by k, with the chance level p0 and the selected optimum of each k
     (accent-coloured rings with a white halo, joined by a line) -- "maximise
     mean recall subject to mean PPV >= p0" is readable straight off the panel;
  E  coverage CDF of every site set -- single tool, the sites the union adds, the
     union, the selected intersection and the GLORI reference (thin lines, one
     per unit; thick line, mean CDF): the added sites are as well covered as the
     best single tool's calls, the reference itself is the worst covered set;
  F  the price over the whole path: FP / 10 kb of the selected k = 1..5 sets on
     the unmodified controls (filled = Curlcake IVT, open = HeLa IVT, one dot per
     unit, line = mean) -- reviewer R3-4's sensitivity-precision trade-off.

Every number is read from the frozen tables of ``40_fig6_combination.py`` and
``61_figS6_sitequality.py`` -- nothing is recomputed here.  House rules: printed
size 1:1, >= 7 pt everywhere, Arial, no gridlines, no in-panel annotation text,
bold panel letters, vector PDF + 300 dpi PNG.

2026-09-30 revision 
------------------------------------------------------
The two grey reference lines of row A no longer borrow the study A/B line style:
each one carries a single fixed style on the whole page -- **thin dashed** for
the recall against every GLORI site, **thin dotted** for the best single tool
(the same dotted line the supplementary companion, Figure S7A, uses for the same
quantity).  As drawn before, the key said "all GLORI (dashed)" while the line was
solid in two of the three columns and solid for mouse study A, and the two greys
were told apart only by dash density.  The two mouse studies are now separated by
a lighter/darker grey instead, and the row keys state the study A/B line
convention.  Every row carries the species titles (2026-09-30), and the y range of row D always contains the chance level p0, so the
dashed p0 line is visible in all three columns (it used to be clipped in
Arabidopsis and mouse).

Layout, 2026-09-30 (the legend moved below every facet)
-----------------------------------------------------------------
Each row key sits **below the panels of its own row** (``KEY_BELOW`` under the
panel's lower edge, clear of the x ticks and the x axis label), not in the title
band above them; the title band is ROW_TITLE high everywhere and seats the row
letter plus the three species titles.

Outputs -> ``figures/figure6/figures/Figure6_rev.{pdf,png}``

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    figures/figure6/src/65_fig6_combination_page.py
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
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.transforms import Bbox

HERE = _RB / "src/harmonisation"
sys.path.insert(0, str(HERE))

from common import config as C                       # noqa: E402
from common.figstyle import apply as figstyle_apply  # noqa: E402
from common.manifest import setup_logger             # noqa: E402

OUT = (_RB / "figures/figure6")
TAB, FIG, LOG = OUT / "tables", OUT / "figures", OUT / "logs"

# ------------------------------------------------------------------ canvas --- #
#: 0.95 * \textwidth of USG.cls (paper 210 mm, page margins 16 mm -> 178 mm);
#: the exported PDF page width must equal this or the print size is a lie.
CANVAS_W = 6.66
#: The figure shares its page with the caption; the page is US Letter with 16 mm
#: margins (text block 9.74 in) and the Figure 6 caption of the revision measures
#: 1.63 in on page 21 (pdftotext -bbox), which leaves 8.11 in for the artwork.

#: species titles now sit on all six rows, which costs 0.03 in of band on each of
#: rows C, E and F; the canvas grew by 0.06 in (8.01 -> 8.07) after the row-C gap
#: and the bottom margin gave back 0.03 in, so the figure still shares page 21
#: with its caption.
MAX_H = 8.10
MIN_PT = 7.0
#: (row letter, [(panel height, gap below the panel) ...]) in inches
ROWS = [("A", [(0.66, 0.48)]),      # trade-off across k
        ("B", [(0.66, 0.48)]),      # trade-off per independent unit
        ("C", [(1.16, 0.32)]),      # 13 tools: own performance | effect + members
        ("D", [(0.66, 0.48)]),      # selection criterion (search space)
        ("E", [(0.52, 0.48), (0.52, 0.52)])]   # E site sets, F control cost

#: every row key moved out of the title band into the space *below* the panels it
#: belongs to; the gaps seat one 7.2 pt key line under the x ticks and the x axis
#: label, and the last row's gap seats the last key too (nothing is left below it
#: but the page margin).
#: Earlier rounds: 2026-09-23 split the single ten-entry foot strip into one key
#: per row (so the reader would not have to travel to the foot of the page for
#: every panel) and brought the page down from 8.46 in so that the figure plus
#: its caption fits the text block.
TOP, ROW_TITLE, SUB_TITLE, LEGEND_H = 0.02, 0.18, 0.18, 0.03
#: Distance from a panel's bottom edge to the top edge of its row key: the x tick
#: labels and the x axis label reach 0.32 in below a panel, and ``layout_report``
#: asserts that no legend text touches either of them.
KEY_BELOW = 0.35
#: left margin: the membership matrix carries the 13 tool names as y tick labels
#: ("NanoSPA_m6A" = 0.69 in at 7.2 pt Arial), plus the tick pad
L_LEFT, L_GAP, L_RIGHT = 0.80, 0.30, 0.06
FS = {"tick": 7.5, "tool": 7.2, "label": 8.0, "title": 8.5, "legend": 7.2,
      "letter": 9.5}

COLS = [{"sp": "Arabidopsis", "groups": ("Arabidopsis_WT",), "colour": "#1E888B",
         "title": "Arabidopsis"},
        {"sp": "Mouse", "groups": ("studyA", "studyB"), "colour": "#E8A33D",
         "colour_b": "#9C5B12", "title": "Mouse"},
        {"sp": "Human", "groups": ("HeLa_WT",), "colour": "#3778A0",
         "title": "Human"}]
UNION_C, ISECT_C = "#1f5c8b", "#e08a2e"
COV_C = "#4a7d3f"
K_COLOURS = ["#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c"]
#: the selected optimum of every k (row D here, panel A of Figure S6): an accent colour
#: no other artist uses, a thick ring and a white halo, so the five points still
#: read on top of the dense k-coloured search cloud (2026-09-21)
HILITE = "#b02418"
TIERS = (1, 2, 5)
TIER_LABEL = {1: "$k$ = 1", 2: "$k$ = 2", 5: "$k$ = 5"}
CONTROLS = ("Curlcake IVT", "HeLa IVT")
#: E sets of the site-quality CDF (the union marginal sites are what combining
#: adds; the intersections are the precision cores read off the same layer)
QUAL_SETS = ["single (k=1)", "union marginal", "union (k=5)",
             "intersection (k=2)", "GLORI reference"]

#: the k of the two composite sets left the key and now lives in the caption, so
#: the printed names are the bare set names -- ``union marginal`` stays distinct
#: from ``union`` because it is the set the union *adds* to the best single tool
QUAL_LABELS = ["single\ntool", "union\nmarginal", "union", "intersection",
               "GLORI\nreference"]
QUAL_COLOURS = ["#8c8c8c", "#1f5c8b", "#6baed6", "#e08a2e", "0.15"]

#: 2026-09-30: the two grey reference lines of row A, each with ONE fixed style
#: page-wide (see the module docstring).  ``GLORI_LS`` also keys the mouse study
#: B series of rows A and B; ``BEST_LS`` is the dotted line of Figure S7A.
GLORI_LS = (0, (4, 2))
BEST_LS = (0, (1.6, 1.6))

#: union recall with the denominator changed from the measurable universe to
#: *every* GLORI site of the group, so the key names it as the share of that pool
GLORI_LABEL = "share of all GLORI sites"
#: grey level of those two lines; a single-study column keeps the mid grey, the
#: mouse column's two independent studies take the lighter/darker pair the rest
#: of the page uses for study A and study B
GREY_SINGLE = (0.35, 0.55)
GREY_STUDY_A = (0.55, 0.72)
GREY_STUDY_B = (0.18, 0.42)
#: the neutral swatch of the "mouse study A / B" style keys
STUDY_KEY_GREY = 0.30


def apply_style() -> None:
    figstyle_apply()
    mpl.rcParams.update({
        "font.size": FS["tick"],
        "axes.titlesize": FS["title"],
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
    """Frozen evidence tables only -- this script never recomputes a number."""
    return {
        "sel": pd.read_csv(TAB / "fig6_combination_selected.tsv", sep="\t"),
        "allk": pd.read_csv(TAB / "fig6_per_unit_allk.tsv", sep="\t"),
        "members": pd.read_csv(TAB / "fig6_selected_members.tsv", sep="\t"),
        "space": pd.read_csv(TAB / "figS6_search_space.tsv", sep="\t"),
        "qual": pd.read_csv(TAB / "figS6_site_quality.tsv", sep="\t"),
        "hist": pd.read_csv(TAB / "figS6_coverage_hist.tsv", sep="\t"),
        "toolq": pd.read_csv(TAB / "figS6_tool_quality.tsv", sep="\t"),
        "single": pd.read_csv(TAB / "fig6_single_tool_metrics.tsv", sep="\t"),
        "effects": pd.read_csv(TAB / "fig6_tool_effects.tsv", sep="\t"),
        "ctl": pd.read_csv(TAB / "fig6_negative_control_fp.tsv", sep="\t"),
    }


def tool_order(members: pd.DataFrame, space: pd.DataFrame) -> list[str]:
    """Tools ordered by the k at which any group first selects them.

    Ties (and tools that are never selected) follow the mean single-tool recall,
    so the membership matrix reads as a staircase and the rows below run from the
    strongest configuration to the weakest.
    """
    first: dict[str, int] = {}
    for r in members[members.in_combination].itertuples():
        first[r.tool] = min(first.get(r.tool, 99), int(r.k))
    k1 = space[space.k == 1]
    recall = k1.assign(tool=k1.combination).groupby("tool").union_recall_mean.mean()
    return sorted(set(members.tool),
                  key=lambda t: (first.get(t, 99), -recall.get(t, 0.0), t))


def _offsets(n: int, half: float = 0.16) -> np.ndarray:
    return np.linspace(-half, half, n) if n > 1 else np.zeros(1)


def _col_sel(df: pd.DataFrame, col: dict, key: str = "group") -> pd.DataFrame:
    return df[df[key].isin(col["groups"])]


def _gcolour(col: dict, gi: int) -> str:
    """Group colour inside a column: mouse study B uses the darker shade."""
    return col["colour"] if gi == 0 else col.get("colour_b", col["colour"])


def _mean(df: pd.DataFrame, column: str) -> float:
    v = pd.to_numeric(df[column], errors="coerce").dropna()
    return float(v.mean()) if len(v) else np.nan


def _sci_label(t: float) -> str:
    """10^n label for a log-axis decade (0.01 -> 10$^{-2}$)."""
    return f"$10^{{{int(round(np.log10(t)))}}}$"


def _decades(lo: float, hi: float, max_n: int = 4) -> list[float]:
    """Decade tick values inside ``[lo, hi]``, thinned until ``max_n`` remain.

    A 0.5 in panel cannot show six 7.5 pt tick labels without them touching.
    """
    ticks = [10.0 ** e for e in range(-3, 6) if lo <= 10.0 ** e <= hi]
    while len(ticks) > max_n:
        ticks = ticks[::2]
    return ticks


# ----------------------------------------------------------------- layout ---- #
def build_stack() -> tuple[float, dict[str, list[list[float]]], dict[str, float]]:
    """Canvas height, per-row list of ``[y0, y1]`` panel rectangles, letters.

    Every row now carries the species titles, so every band is ROW_TITLE high;
    the second panel of a row (row F) takes SUB_TITLE, which is the same value.
    """
    fig_h = (TOP + LEGEND_H
             + sum((ROW_TITLE if i == 0 else SUB_TITLE) + h + g
                   for _, panels in ROWS
                   for i, (h, g) in enumerate(panels)))
    rects: dict[str, list[list[float]]] = {}
    letters: dict[str, float] = {}
    pool, idx = "ABCDEFGH", 0
    y = fig_h - TOP
    for row, panels in ROWS:
        rects[row] = []
        for i, (h, gap) in enumerate(panels):
            if idx < len(pool):            # one letter per panel, not per row
                letters[pool[idx]] = y - 0.01
                idx += 1
            y -= ROW_TITLE if i == 0 else SUB_TITLE
            rects[row].append([y - h, y])
            y -= h + gap
    return fig_h, rects, letters


def _rect(x_in: float, y_in: float, w_in: float, h_in: float,
          fig_w: float, fig_h: float) -> list[float]:
    return [x_in / fig_w, y_in / fig_h, w_in / fig_w, h_in / fig_h]


def panel_axes(fig, rects, row, index, fig_w, fig_h) -> list:
    """One axes per species column for panel ``index`` of row ``row``."""
    y0, y1 = rects[row][index]
    w = (CANVAS_W - L_LEFT - L_RIGHT - 2 * L_GAP) / 3
    axes = []
    for i, _ in enumerate(COLS):
        ax = fig.add_axes(_rect(L_LEFT + i * (w + L_GAP), y0, w, y1 - y0,
                                fig_w, fig_h))
        ax.tick_params(length=2.2, pad=1.4)
        axes.append(ax)
    return axes


# ----------------------------------------------------------------- panels ---- #
def _ref_greys(col: dict, gi: int) -> tuple[float, float]:
    """Grey level of the two reference lines, ``(all GLORI, best single tool)``.

    A single-study column (Arabidopsis, human) keeps the mid greys; in the mouse
    column the two independent studies take the lighter/darker pair, which is how
    the page tells study A from study B wherever a line style is not available.
    """
    if len(col["groups"]) == 1:
        return GREY_SINGLE
    return GREY_STUDY_A if gi == 0 else GREY_STUDY_B


def panel_a(ax, col: dict, sel: pd.DataFrame) -> dict:
    """Mean trade-off with the SD across units (mouse: one line per study).

    Returns the two grey reference series so ``build_figure`` can check that the
    key's swatches are the styles actually drawn (2026-09-30).
    """
    d = _col_sel(sel, col).sort_values(["group", "k"])
    refs: dict[str, list] = {"all_glori": [], "best_single": [], "series": []}
    for gi, (_, sub) in enumerate(d.groupby("group", sort=True)):
        ls = "-" if gi == 0 else GLORI_LS
        grey_all, grey_best = _ref_greys(col, gi)
        for column, colour, mk, mfc, lw, ms in (
                ("union_recall", UNION_C, "o", None, 1.3, 3.0),
                ("union_precision", UNION_C, "o", "white", 1.2, 3.0),
                ("isect_precision", ISECT_C, "s", "white", 1.1, 2.8)):
            series = ax.errorbar(sub.k, 100 * sub[f"{column}_mean"],
                                 yerr=100 * sub[f"{column}_sd"].fillna(0.0).to_numpy(),
                                 color=colour, marker=mk, ms=ms, mfc=mfc or colour,
                                 mew=0.9 if mfc else 0.0, lw=lw, ls=ls, capsize=1.3,
                                 elinewidth=0.6,
                                 zorder=4 if column != "union_precision" else 3)
            refs["series"].append((gi, series.lines[0]))
        # recall against **every** GLORI site: one fixed style in every column
        # (it used to inherit the study A/B style, which made the key entry
        # "all GLORI (dashed)" false in two of the three columns)
        cont = ax.errorbar(sub.k, 100 * sub.union_recall_all_mean,
                           yerr=100 * sub.union_recall_all_sd.fillna(0.0).to_numpy(),
                           color=(grey_all, grey_all, grey_all), marker="",
                           lw=0.9, ls=GLORI_LS, capsize=1.3, elinewidth=0.6,
                           zorder=2)
        refs["all_glori"].append(cont.lines[0])
        # best single tool: fixed thin dotted, one level per study
        ref = 100 * float(sub[sub.k == 1].union_recall_mean.iloc[0])
        refs["best_single"].append(
            ax.axhline(ref, color=(grey_best, grey_best, grey_best), lw=0.9,
                       ls=BEST_LS, zorder=1))
    ax.set_xticks(range(1, 6))
    ax.set_ylim(0, 105)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlabel("tools combined ($k$)", fontsize=FS["label"])
    ax.set_ylabel("recall / PPV (%)", fontsize=FS["label"])
    return refs


def panel_b(ax, col: dict, allk: pd.DataFrame) -> dict:
    """Per-unit trajectories of the same trade-off over k = 1..5.

    Returns the drawn series as ``(study index, artist)`` pairs: the row key
    names the solid (study A) and dashed (study B) line styles, and
    ``build_figure`` checks that those swatches are what the panel draws.
    """
    d = _col_sel(allk, col)
    series: list[tuple[int, object]] = []
    for gi, (_, units) in enumerate(d.groupby("group", sort=True)):
        ls = "-" if gi == 0 else GLORI_LS
        tags = sorted(units.unit.unique())
        for tag, off in zip(tags, _offsets(len(tags), 0.13)):
            u = units[units.unit == tag].sort_values("k")
            series.append((gi, ax.plot(u.k + off, 100 * u.union_recall,
                                       color=UNION_C, ls=ls, lw=1.1, marker="o",
                                       ms=2.6, zorder=3)[0]))
            series.append((gi, ax.plot(u.k + off, 100 * u.isect_precision,
                                       color=ISECT_C, ls=ls, lw=1.0, marker="s",
                                       ms=2.4, mfc="white", mew=0.8, zorder=2)[0]))
    ax.set_xticks(range(1, 6))
    ax.set_ylim(0, 105)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlabel("tools combined ($k$)", fontsize=FS["label"])
    ax.set_ylabel("recall / PPV (%)", fontsize=FS["label"])
    return {"series": series}


def panel_c_own(ax, col: dict, single: pd.DataFrame, toolq: pd.DataFrame,
                tools: list[str], show_labels: bool = True) -> None:
    """C, left block: what each of the 13 configurations achieves on its own.

    ``show_labels`` is true for column 0 only: the three columns share the same
    thirteen rows, and one set of names is 0.69 in wide -- a second and third set
    would be printed on the effect matrix of the column to the left (measured
    2026-09-21, 20 px of overlap per label).
    """
    index = {t: i for i, t in enumerate(tools)}
    d = _col_sel(toolq, col)
    for t in tools:
        y = index[t]
        s = d[d.tool == t]
        ends: list[tuple[float, float, float]] = []      # (recall, ppv, y)
        for gi, (_, sub) in enumerate(s.groupby("group", sort=True)):
            mk, mc = ("o" if gi == 0 else "s"), _gcolour(col, gi)
            # two studies in one column would land on exactly the same row (one
            # unit each, so no jitter): offset them inside the row instead
            dy = 0.0 if len(col["groups"]) == 1 else (0.18 if gi == 0 else -0.18)
            for off, (_, r) in zip(_offsets(len(sub), 0.26), sub.iterrows()):
                ax.scatter([100 * r.frac_glori], [y + dy + off], marker=mk, s=8,
                           color=mc, linewidth=0, zorder=4)
                ax.scatter([100 * r.recall_in_universe], [y + dy + off],
                           marker=mk, s=8, facecolor="white", edgecolor=mc,
                           linewidth=0.6, zorder=3)
            m = _col_sel(single, col)
            m = m[m.tool.eq(t)]
            if len(m):
                ppv = 100 * _mean(m, "union_precision_mean")
                rec = 100 * _mean(m, "union_recall_mean")
                ax.plot([rec, ppv], [y + dy, y + dy], color="0.60", lw=0.6,
                        zorder=2)
                ax.plot([ppv - 4, ppv + 4], [y + dy, y + dy], color="0.20",
                        lw=1.3, zorder=5, solid_capstyle="butt")
                ends.append((rec, ppv, y + dy))
        # the two studies of a mouse column: join their recall and their PPV ends
        if len(ends) == 2:
            (r0, p0, _), (r1, p1, _) = ends
            ax.plot([r0, r1], [y + 0.18, y - 0.18], color="0.65", lw=0.6,
                    zorder=2)
            ax.plot([p0, p1], [y + 0.18, y - 0.18], color="0.65", lw=0.6,
                    zorder=2)
    ax.set_xlim(-4, 104)
    ax.set_ylim(len(tools) - 0.4, -0.6)          # row 0 = strongest tool on top
    ax.set_xticks([0, 50, 100])
    # the row positions stay on every column (the 13 rows must line up across the
    # figure); the names themselves are printed once, in column 0
    ax.set_yticks(range(len(tools)))
    ax.set_yticklabels(tools if show_labels else [], fontsize=FS["tool"])
    ax.set_xlabel("recall / PPV (%)", fontsize=FS["tool"], labelpad=1.5)


def panel_c_matrix(ax, col: dict, members: pd.DataFrame, effects: pd.DataFrame,
                   tools: list[str], lim: float) -> None:
    """C, right block: main effect of each configuration at fixed k (cell colour)
    and the configurations the criterion actually selects (markers)."""
    index = {t: i for i, t in enumerate(tools)}
    cmap = plt.get_cmap("RdBu_r")
    d = _col_sel(effects, col)
    for _, r in d.iterrows():
        v = 100 * float(r.effect_recall)
        shade = cmap(0.5 + 0.5 * float(np.clip(v / lim, -1, 1)))
        ax.add_patch(mpl.patches.Rectangle((int(r.k) - 0.5, index[r.tool] - 0.40),
                                           1.0, 0.80, facecolor=shade,
                                           edgecolor="white", lw=0.35, zorder=1))
    for gi, (_, sub) in enumerate(_col_sel(members, col).groupby("group",
                                                                 sort=True)):
        dx = 0.0 if len(col["groups"]) == 1 else (-0.20 if gi == 0 else 0.20)
        mk = "o" if gi == 0 else "s"
        hit = sub[sub.in_combination]
        ax.scatter(hit.k.astype(int) + dx, [index[t] for t in hit.tool],
                   marker=mk, s=13, color=_gcolour(col, gi), edgecolor="white",
                   linewidth=0.45, zorder=3)
    ax.set_xlim(0.5, 5.5)
    ax.set_xticks(range(1, 6))
    ax.set_ylim(len(tools) - 0.4, -0.6)
    ax.set_yticks([])
    ax.set_xlabel("tools ($k$)", fontsize=FS["tool"], labelpad=1.5)


def panel_d(ax, col: dict, space: pd.DataFrame) -> dict:
    """The criterion: every enumerated combination, p0, and the optima.

    The y range is guaranteed to contain p0 (2026-09-30): the dashed chance-level
    line was drawn in every column but sat below the axis floor in Arabidopsis
    (p0 1.50 % < 2.18 %) and mouse (0.70 % < 0.89 %), so the key named a line two
    of the three panels did not show.  Returns p0 and the limits it was drawn in
    so ``build_figure`` can assert the contract.
    """
    d = space[(space.species == col["sp"]) & space.group.isin(col["groups"])]
    for gi, (_, sub) in enumerate(d.groupby("group", sort=True)):
        alpha = 0.75 if gi == 0 else 0.45
        for k, colour in zip(range(1, 6), K_COLOURS):
            s = sub[sub.k == k]
            ax.scatter(100 * s.union_recall_mean, 100 * s.union_precision_mean,
                       s=2.0, color=colour, alpha=alpha, linewidths=0,
                       rasterized=True, zorder=2)
        sel = sub[sub.selected].sort_values("k")
        ls = "-" if gi == 0 else (0, (3, 2))
        ax.plot(100 * sel.union_recall_mean, 100 * sel.union_precision_mean,
                color=HILITE, lw=1.1, ls=ls, zorder=5)
        rings = ax.scatter(100 * sel.union_recall_mean,
                           100 * sel.union_precision_mean, s=48,
                           facecolors="none", edgecolors=HILITE, linewidths=1.2,
                           zorder=6)
        # white halo: lifts the ring off the 2 379 semi-transparent cloud dots
        rings.set_path_effects([pe.withStroke(linewidth=2.4,
                                              foreground="white")])
    p0 = 100 * float(d.p0_chance_precision.iloc[0])
    ax.axhline(p0, color="0.30", lw=0.9, ls=(0, (3, 2)), zorder=1)
    ax.set_yscale("log")
    #: the chance level is part of the selection criterion, so the axis floor is
    #: lowered until the dashed line sits inside the panel half a decade above it
    #: -- the same visible margin in all three columns, and low enough that the
    #: row shares one tick set (1-5-20) instead of a per-column one
    cloud_lo = max(0.6, 100 * float(d.union_precision_mean.min()) * 0.6)
    lo = min(cloud_lo, 0.5 * p0)
    hi = min(60.0, 100 * float(d.union_precision_mean.max()) * 2.0)
    ax.set_ylim(lo, hi)
    #: 2026-09-30: five decade ticks (1-2-5-10-20) no longer fit: with the floor
    #: below 1 the bottom pair of labels touched (caught by layout_report), and a
    #: 0.66 in panel cannot hold them; 1-5-20 keeps every label 13 pt apart
    ticks = [t for t in (1, 5, 20, 100) if lo <= t <= hi]
    ax.set_yticks(ticks)
    ax.set_yticklabels([str(t) for t in ticks])
    ax.minorticks_off()
    ax.set_xlim(0, 100 * float(d.union_recall_mean.max()) * 1.12)
    # fixed ticks: the automatic ones crowd the right edge of the row-D band
    ax.set_xticks([0, 20, 40, 60, 80])
    ax.set_xlabel("mean union recall (%)", fontsize=FS["label"])
    ax.set_ylabel("mean PPV (%)", fontsize=FS["label"])
    return {"p0": p0, "lo": lo, "hi": hi, "yticks": ticks}


def panel_e1(ax, col: dict, hist: pd.DataFrame) -> None:
    """Coverage CDF of every site set: thin = one unit, thick = mean CDF.

    Reading: a curve to the right is better covered; the union marginal sites are
    as well covered as the best single tool's calls, the GLORI reference is the
    worst covered set of all.
    """
    d = _col_sel(hist, col)
    grid = np.sort(hist.bin_hi.unique())               # common x for every unit
    for i, name in enumerate(QUAL_SETS):
        s = d[d.set == name]
        curves = []
        for _, sg in s.groupby(["group", "unit"], sort=True):
            u = sg.set_index("bin_hi").n_sites.reindex(grid, fill_value=0)
            n = float(u.sum())
            if n <= 0:
                continue
            y = np.cumsum(u.to_numpy(dtype=float)) / n
            ax.plot(grid, y, color=QUAL_COLOURS[i], lw=0.6, alpha=0.55, zorder=2)
            curves.append(y)
        if curves:
            ax.plot(grid, np.mean(np.vstack(curves), axis=0),
                    color=QUAL_COLOURS[i], lw=1.5, zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(4, 5000)
    ax.set_ylim(0, 1.0)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_yticklabels(["0", "0.5", "1"])
    ax.minorticks_off()
    ax.set_xticks([10, 100, 1000])
    ax.set_xticklabels(["10", "100", "1000"])
    ax.set_ylabel("cumulative\nfraction", fontsize=FS["label"], linespacing=1.1)
    ax.set_xlabel("coverage (reads)", fontsize=FS["label"], labelpad=1.5)


def panel_e2(ax, col: dict, ctl: pd.DataFrame) -> None:
    """Price of the selected sets over the whole path k = 1..5 (both controls)."""
    d = _col_sel(ctl, col, key="chosen_for")
    for i, cname in enumerate(CONTROLS):
        sub = d[d.control == cname]
        for gi, (_, sg) in enumerate(sub.groupby("chosen_for", sort=True)):
            mk, mc = ("o" if gi == 0 else "s"), _gcolour(col, gi)
            mean_x, mean_y = [], []
            for k in range(1, 6):
                s = sg[sg.k == k]
                off = -0.11 if i == 0 else 0.11
                for _, r in s.iterrows():
                    ax.plot([k + off], [r.fp_per_10kb], ls="none", marker=mk,
                            ms=2.6, mfc=mc if i == 0 else "white", mec=mc,
                            mew=0.6, zorder=3)
                if len(s):
                    mean_x.append(k + off)
                    mean_y.append(float(s.fp_per_10kb.mean()))
            ax.plot(mean_x, mean_y, color=mc if i == 0 else "0.60",
                    lw=1.2 if i == 0 else 0.8, ls="-" if i == 0 else (0, (3, 2)),
                    zorder=4)
    vals = pd.to_numeric(d.fp_per_10kb, errors="coerce").dropna()
    vals = vals[vals > 0]
    lo = max(1e-2, 10 ** np.floor(np.log10(vals.min()) - 0.1)) if len(vals) else 1e-2
    hi = 10 ** np.ceil(np.log10(vals.max() + 1e-9) + 0.1) if len(vals) else 1e3
    ax.set_yscale("log")
    ax.set_ylim(lo, hi)
    ticks = _decades(lo, hi, 4)
    ax.set_yticks(ticks)
    ax.set_yticklabels([_sci_label(t) for t in ticks])
    ax.minorticks_off()
    ax.set_xlim(0.5, 5.5)
    ax.set_xticks(range(1, 6))
    ax.set_xlabel("tools combined ($k$)", fontsize=FS["label"], labelpad=1.5)
    ax.set_ylabel("FP / 10 kb", fontsize=FS["label"])


# ----------------------------------------------------------------- legend ---- #
def handles() -> list:
    """Every legend entry of the page, in the order the rows introduce them.

    Kept as one flat list (the old bottom band used all ten); ``row_handles``
    picks the subset each row needs.
    """
    return [
        Line2D([], [], color=UNION_C, marker="o", ms=3.2, lw=1.3,
               label="union recall"),
        #: 2026-09-30: the label no longer claims a line style; the swatch *is*
        #: the style drawn in every column (``GLORI_LS``), mouse study B included
        Line2D([], [], color=(GREY_SINGLE[0],) * 3, lw=0.9, ls=GLORI_LS,
               label=GLORI_LABEL),
        Line2D([], [], color=UNION_C, marker="o", ms=3.2, mfc="white", mew=0.9,
               lw=1.2, label="union PPV"),
        Line2D([], [], color=ISECT_C, marker="s", ms=2.8, mfc="white", mew=0.9,
               lw=1.1, label="intersection PPV"),
        Line2D([], [], color=K_COLOURS[2], marker="o", ms=2.2, ls="none",
               label="all enumerated combinations"),
        Line2D([], [], color=HILITE, marker="o", ms=4.0, mfc="none", mew=1.2,
               lw=1.1, path_effects=[pe.withStroke(linewidth=2.4,
                                                   foreground="white")],
               label="selected optimum (k = 1-5)"),

        #: "all GLORI sites" series, and the two must be told apart at a glance;
        #: the supplementary companion (Figure S7A) dots this same quantity
        Line2D([], [], color=(GREY_SINGLE[1],) * 3, lw=0.9, ls=BEST_LS,
               label="best single tool"),
        Line2D([], [], color="0.30", lw=0.9, ls=(0, (3, 2)),
               label="chance level p0"),
        Line2D([], [], color="0.45", marker="o", ms=2.8, ls="none", mfc="0.45",
               mec="0.45", label="IVT controls"),
        #: the mouse column carries two independent studies in every row; where a
        #: line style can carry them, study B is the dashed one
        Line2D([], [], color=(STUDY_KEY_GREY,) * 3, lw=1.3, ls="-",
               label="mouse study A"),
        Line2D([], [], color=(STUDY_KEY_GREY,) * 3, lw=1.3, ls=GLORI_LS,
               label="mouse study B"),
    ]


#: which entry belongs to which row: "all the legends in one band at the bottom"
#: forced the reader to travel to the foot of the page for every panel

#: colours -- the union recall (blue, filled circles) and the intersection PPV
#: (orange, open squares) -- and its own key now names those two quantities, as
#: the other rows do, instead of describing what a dot is.  That a point is one
#: independent sequencing unit is stated in the caption (Fig. 6B).

#: the mouse column is the only one that carries two independent studies, and a
#: dashed line there means study B, in the rows where line style is the encoding.
#: Row A carries the shortest key that still names every mark it draws: the two
#: study entries ride in row B's key (and in the caption), because row A's line
#: would otherwise run past the page edge once the grey curve was renamed.
ROW_LEGEND = {"A": ["union recall", GLORI_LABEL, "union PPV",
                    "intersection PPV", "best single tool"],
              "B": ["union recall", "intersection PPV",
                    "mouse study A", "mouse study B"],
              "D": ["all enumerated combinations", "selected optimum (k = 1-5)",
                    "chance level p0"]}


def qual_handles() -> list:
    """Row E: the five site sets of the coverage CDF."""
    return [Line2D([], [], color=QUAL_COLOURS[i], lw=1.5,
                   label=lbl.replace("\n", " "))
            for i, lbl in enumerate(QUAL_LABELS)]


def control_handles() -> list:
    """Row F: filled = Curlcake IVT, open = HeLa IVT (mean line solid/dashed)."""
    return [Line2D([], [], color="0.35", marker="o", ms=2.8, ls="-", lw=1.2,
                   mfc="0.35", mec="0.35", label="Curlcake IVT"),
            Line2D([], [], color="0.45", marker="s", ms=2.6, ls=(0, (3, 2)),
                   lw=0.8, mfc="white", mec="0.45", label="Human IVT")]


def row_legend(fig, handles_: list, y_in: float, fig_h: float):
    """One legend line, left-aligned *below* the row of panels it belongs to.

    2026-09-30 : the key used to
    sit in the row's title band, above the panels; it is now anchored under the
    row's own panels, ``KEY_BELOW`` inches beneath their lower edge so it clears
    the x ticks and the x axis label.  ``check_keys_below_panels`` asserts that
    every key ends up below its panels.
    """
    return fig.legend(handles=handles_, loc="upper left",
                      bbox_to_anchor=(0.035, y_in / fig_h),
                      ncol=len(handles_), frameon=False, fontsize=FS["legend"],
                      handlelength=1.1, columnspacing=0.85, handletextpad=0.30,
                      labelspacing=0.25, borderpad=0.0)


# ------------------------------------------------------------------ checks --- #
def layout_report(fig, tools: list[str] | None = None) -> dict:
    """Printed-size audit: font sizes, text collisions, legend overflow.

    With ``tools`` given, the report also counts the row-C tool-name labels and
    checks that none of them is printed over a panel it does not belong to (the
    failure mode of an earlier draft: all three columns carried the names, and
    the 0.69 in labels of columns 1-2 landed on the effect matrix to their left).
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    items = []
    for t in fig.findobj(Text):
        if not t.get_visible() or not str(t.get_text()).strip():
            continue
        try:
            bb = t.get_window_extent(renderer=rend)
        except Exception:                                      # pragma: no cover
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        items.append({"s": str(t.get_text())[:32], "fs": float(t.get_fontsize()),
                      "bb": [float(bb.x0), float(bb.y0), float(bb.x1),
                             float(bb.y1)]})
    overlaps = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            inter = Bbox.intersection(Bbox.from_extents(*items[i]["bb"]),
                                      Bbox.from_extents(*items[j]["bb"]))
            if inter is not None and inter.width > 1.2 and inter.height > 1.2:
                overlaps.append((items[i]["s"], items[j]["s"]))
    fonts: dict[float, int] = {}
    for it in items:
        fonts[it["fs"]] = fonts.get(it["fs"], 0) + 1
    legs = []
    for lg in fig.legends:
        bb = lg.get_window_extent(renderer=rend)
        legs.append([float(bb.x0), float(bb.x1), float(bb.y0), float(bb.y1)])
    w_px = float(fig.get_size_inches()[0]) * float(fig.dpi)
    h_px = float(fig.get_size_inches()[1]) * float(fig.dpi)
    # text running off the canvas is silently truncated when the PDF is placed
    # in the manuscript -- the failure that clipped the tool labels of an earlier
    # draft of this figure
    clipped = [(it["s"], round(it["bb"][0], 1), round(it["bb"][2], 1))
               for it in items
               if it["bb"][0] < -0.5 or it["bb"][2] > w_px + 0.5]
    n_tool_labels, intrusion = 0, []
    if tools:
        names = set(tools)
        axes_bb = [(ax, ax.get_window_extent(renderer=rend)) for ax in fig.axes]
        n_tool_labels = sum(1 for it in items if it["s"] in names)
        for t in fig.findobj(Text):
            if not t.get_visible() or str(t.get_text()) not in names:
                continue
            bb = t.get_window_extent(renderer=rend)
            for ax, ab in axes_bb:
                if ax is t.axes:                 # its own axes does not count
                    continue
                inter = Bbox.intersection(bb, ab)
                if inter is not None and inter.width > 1.0 and inter.height > 1.0:
                    intrusion.append((str(t.get_text()), round(inter.width, 1),
                                      round(inter.height, 1)))
    return {"n_text": len(items),
            "min_fontsize": min((it["fs"] for it in items), default=float("nan")),
            "fonts": dict(sorted(fonts.items())), "overlaps": overlaps[:20],
            "legends": legs, "clipped": clipped, "n_tool_labels": n_tool_labels,
            "tool_label_intrusion": intrusion,
            "canvas_px": [w_px, h_px]}


def check_tables(t: dict[str, pd.DataFrame], tools: list[str], logger) -> None:
    """Anchors: the figure may only draw what the frozen tables state."""
    sel, allk, mem = t["sel"], t["allk"], t["members"]
    assert len(tools) == 13, f"expected 13 configurations, got {len(tools)}"
    for col in COLS:
        for g in col["groups"]:
            s = sel[sel.group == g]
            assert sorted(s.k) == [1, 2, 3, 4, 5], f"{g}: selected k incomplete"
            assert s.union_recall_all_mean.notna().all(), \
                f"{g}: union_recall_all missing (re-run 40 with the new column)"
            m = mem[mem.group == g]
            assert len(m) == 65, f"{g}: membership matrix has {len(m)} rows"
            a = allk[allk.group == g]
            for k in range(1, 6):
                u = a[a.k == k]
                want = int(s[s.k == k].n_units.iloc[0])
                assert len(u) == want, f"{g} k={k}: {len(u)} vs {want} unit rows"
                got = 100 * float(u.union_recall.mean())
                ref = 100 * float(s[s.k == k].union_recall_mean.iloc[0])
                assert abs(got - ref) < 1e-3, \
                    f"{g} k={k}: per-unit mean {got:.4f} != selected {ref:.4f}"
                assert u.union_recall_all.notna().all(), f"{g} k={k}: recall_all NaN"
            e = t["effects"][t["effects"].group == g]
            assert len(e) == 65, f"{g}: tool effects has {len(e)} rows"
            assert not e.effect_recall.isna().any(), f"{g}: NaN tool effect"
            c = t["ctl"][t["ctl"].chosen_for == g]
            assert sorted(set(c.k)) == [1, 2, 3, 4, 5], f"{g}: control k incomplete"
            q = t["toolq"][t["toolq"].group == g]
            assert set(q.tool) == set(tools), f"{g}: per-tool quality incomplete"
            s1 = t["single"][t["single"].group == g]
            assert set(s1.tool) == set(tools), f"{g}: single-tool plane incomplete"
        h = t["hist"][t["hist"].species == col["sp"]]
        assert set(h.set) >= set(QUAL_SETS), f"{col['sp']}: coverage sets missing"
        sp = t["space"][t["space"].species == col["sp"]]
        assert len(sp) in (2379, 4758), f"{col['sp']}: {len(sp)} search rows"
        assert set(t["qual"][t["qual"].species == col["sp"]].set) >= set(QUAL_SETS)
    logger.info("table anchors OK: 4 groups x k = 1..5 (both denominators), "
                "membership 65 rows/group, per-unit means match to 1e-3 pp, "
                "search space complete")


def _ls_key(art) -> tuple:
    """Comparable dash pattern of a Line2D.

    matplotlib reports every custom dash tuple back as the string ``"--"`` and
    keeps the pattern itself in ``_unscaled_dash_pattern``; the whole point of
    the checks below is that the two grey reference lines are *not*
    interchangeable, so the pattern has to be compared rather than its name.
    """
    pat = getattr(art, "_unscaled_dash_pattern", None)
    if pat is None:                                  # pragma: no cover (mpl move)
        return (str(art.get_linestyle()),)
    off, onoff = pat
    return (round(float(off or 0.0), 3),
            tuple(round(float(v), 3) for v in (onoff or ())))


def _style_key(ls) -> tuple:
    """``_ls_key`` of a line *style* (a dash tuple, ``"-"``, ``"none"`` ...)."""
    probe = Line2D([], [])
    probe.set_linestyle(ls)
    return _ls_key(probe)


def check_key_vs_artwork(panel_a_refs: list[dict], panel_b_series: dict,
                         d_contracts: list[dict], all_h: dict) -> None:
    """The row keys may only claim what the panels actually draw (2026-09-30).

    Every assertion here answers a complaint of the 2026-09-30 review: the key
    said "all GLORI (dashed)" while two of the three columns drew that line solid
    and the two greys differed only in dash density; the row-D key named a chance
    level that two of the three panels had clipped below the axis; and nothing on
    the page said what a dashed line in the mouse column means.
    """
    glori, best, series = [], [], []
    for r in panel_a_refs:
        glori += r["all_glori"]
        best += r["best_single"]
        series += r["series"]
    series += list(panel_b_series["series"])
    assert glori and best and series, "panel A/B drew no reference line"
    for art in glori:
        assert _ls_key(art) == _style_key(GLORI_LS), \
            (f'"all GLORI sites" drawn {art.get_linestyle()} / {_ls_key(art)}, '
             f"the key swatch is {GLORI_LS}")
    for art in best:
        assert _ls_key(art) == _style_key(BEST_LS), \
            (f'"best single tool" drawn {art.get_linestyle()} / {_ls_key(art)}, '
             f"the key swatch is {BEST_LS}")
    assert _ls_key(glori[0]) != _ls_key(best[0]), \
        "the two grey reference lines share one line style again"
    # the study key: study A solid, study B dashed, in rows A and B alike
    for gi, art in series:
        want = "-" if gi == 0 else GLORI_LS
        assert _ls_key(art) == _style_key(want), \
            (f"study {gi} series drawn {art.get_linestyle()}, the "
             f'"mouse study {"AB"[gi]}" key swatch is {want}')
    for label, want in ((GLORI_LABEL, GLORI_LS),
                        ("best single tool", BEST_LS),
                        ("mouse study A", "-"), ("mouse study B", GLORI_LS)):
        assert all_h[label].get_label() == label, \
            f'key entry "{label}" is labelled {all_h[label].get_label()!r}'
        assert _ls_key(all_h[label]) == _style_key(want), \
            f'key entry "{label}" is {all_h[label].get_linestyle()}, not {want}'
    for c in d_contracts:
        assert c["lo"] < c["p0"] < c["hi"], \
            (f"chance level p0 = {c['p0']:.2f} % outside the drawn range "
             f"{c['lo']:.2f}..{c['hi']:.2f} % of its panel")
        assert c["yticks"] == d_contracts[0]["yticks"], \
            f"row D panels carry different y ticks: {c['yticks']}"


def check_keys_below_panels(fig, keyed: list[tuple[object, float]]) -> None:
    """Every row key must be printed *below* the panels it belongs to.

    2026-09-30 : the keys moved
    from the title band into the space under their own row.  The requirement is
    geometric, so it is asserted against the rendered bounding boxes rather than
    trusted to the anchoring arithmetic.
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    dpi = fig.dpi
    for legend, panel_bottom in keyed:
        bb = legend.get_window_extent(renderer=rend)
        top, bottom = bb.y1 / dpi, bb.y0 / dpi
        assert top <= panel_bottom, \
            (f"row key reaches {top:.2f} in, above its own panel bottom "
             f"{panel_bottom:.2f} in")
        assert bottom >= panel_bottom - KEY_BELOW - 0.16, \
            (f"row key ends {panel_bottom - bottom:.2f} in below its panels, "
             f"more than the {KEY_BELOW + 0.16:.2f} in reserved for it")


# -------------------------------------------------------------------- main --- #
def build_figure(t: dict[str, pd.DataFrame],
                 tools: list[str]) -> tuple[plt.Figure, dict]:
    apply_style()
    xpos = {tool: float(i) for i, tool in enumerate(tools)}
    fig_h, rects, letters = build_stack()
    assert fig_h <= MAX_H, f"canvas {fig_h:.2f} in exceeds the {MAX_H} in budget"
    fig = plt.figure(figsize=(CANVAS_W, fig_h))

    ax_a = panel_axes(fig, rects, "A", 0, CANVAS_W, fig_h)
    ax_b = panel_axes(fig, rects, "B", 0, CANVAS_W, fig_h)
    ax_d = panel_axes(fig, rects, "D", 0, CANVAS_W, fig_h)
    ax_e1 = panel_axes(fig, rects, "E", 0, CANVAS_W, fig_h)
    ax_e2 = panel_axes(fig, rects, "E", 1, CANVAS_W, fig_h)

    # row C is split per column: own performance of every configuration (left)
    # and its main effect at fixed k plus the selected members (right)
    y0, y1 = rects["C"][0]
    wcol = (CANVAS_W - L_LEFT - L_RIGHT - 2 * L_GAP) / 3
    w1 = wcol * 0.56
    w2 = wcol - w1 - 0.10
    ax_c1, ax_c2 = [], []
    for i, _ in enumerate(COLS):
        x0 = L_LEFT + i * (wcol + L_GAP)
        ax_c1.append(fig.add_axes(_rect(x0, y0, w1, y1 - y0, CANVAS_W, fig_h)))
        ax_c2.append(fig.add_axes(_rect(x0 + w1 + 0.10, y0, w2, y1 - y0,
                                        CANVAS_W, fig_h)))
    eff = 100 * t["effects"].effect_recall.to_numpy(dtype=float)
    lim = float(np.nanmax(np.abs(eff))) if eff.size else 1.0

    refs_a: list[dict] = []
    refs_b: dict = {"series": []}
    contracts_d: list[dict] = []
    for j, col in enumerate(COLS):
        refs_a.append(panel_a(ax_a[j], col, t["sel"]))
        refs_b["series"] += panel_b(ax_b[j], col, t["allk"])["series"]
        panel_c_own(ax_c1[j], col, t["single"], t["toolq"], tools,
                    show_labels=(j == 0))
        panel_c_matrix(ax_c2[j], col, t["members"], t["effects"], tools, lim)
        contracts_d.append(panel_d(ax_d[j], col, t["space"]))
        panel_e1(ax_e1[j], col, t["hist"])
        panel_e2(ax_e2[j], col, t["ctl"])


    # every row, so a reader who lands on any one of them can key the three columns
    # without travelling back to row A.  Row C is the one row made of two blocks per
    # column, so its title is centred over the whole column by hand (x in axes
    # fractions of the left block: half a column divided by that block's width).
    for axs in (ax_a, ax_b, ax_d, ax_e1, ax_e2):
        for ax, col in zip(axs, COLS):
            ax.set_title(col["title"], fontsize=FS["title"], fontweight="bold",
                         pad=3)
    for i, col in enumerate(COLS):
        ax_c1[i].set_title(col["title"], x=(wcol / 2) / w1, fontsize=FS["title"],
                           fontweight="bold", pad=3)

    # one y axis label per row and column 0 only: the columns differ by species,
    # and a second label would eat into the neighbouring panel's data area
    for axs in (ax_a, ax_b, ax_d, ax_e1, ax_e2, ax_c1, ax_c2):
        for ax in axs[1:]:
            ax.set_ylabel("")

    for row, y_in in letters.items():
        fig.text(0.008, y_in / fig_h, row, fontsize=FS["letter"],
                 fontweight="bold", va="top", ha="left")

    # ---- one key per row, below that row's own panels --------------------- #

    # to sit in the title band, above the panels, next to the row letter.
    all_h = {h.get_label(): h for h in handles()}
    keyed: list[tuple[object, float]] = []          # (legend, panel bottom edge)
    for row, labels in ROW_LEGEND.items():
        bottom = rects[row][-1][0]
        keyed.append((row_legend(fig, [all_h[l] for l in labels],
                                 bottom - KEY_BELOW, fig_h), bottom))
    # rows E and F are the two panels of the "E" entry of ROWS
    for index, hh in ((0, qual_handles()), (1, control_handles())):
        bottom = rects["E"][index][0]
        keyed.append((row_legend(fig, hh, bottom - KEY_BELOW, fig_h), bottom))
    check_key_vs_artwork(refs_a, refs_b, contracts_d, all_h)
    check_keys_below_panels(fig, keyed)
    return fig, {"fig_h": fig_h, "tools": tools, "xpos": xpos,
                 "effect_lim": lim, "contracts_d": contracts_d, "keyed": keyed}


def main() -> None:
    logger = setup_logger("65_fig6_combination_page", log_dir=LOG)
    apply_style()
    t = load()
    tools = tool_order(t["members"], t["space"])
    check_tables(t, tools, logger)
    fig, meta = build_figure(t, tools)
    rep = layout_report(fig, tools)
    logger.info("canvas %.2f x %.2f in (printed 1:1 at 0.95\\textwidth = %.2f in)",
                CANVAS_W, meta["fig_h"], CANVAS_W)
    logger.info("tool order: %s", " > ".join(meta["tools"]))
    logger.info("row C effect colour scale: +/- %.2f pp (mean recall at fixed k)",
                meta["effect_lim"])
    logger.info("row D chance level inside every panel "
                "(p0, floor, ceiling, ticks): %s",
                [(*[round(c[k], 2) for k in ("p0", "lo", "hi")], c["yticks"])
                 for c in meta["contracts_d"]])
    logger.info('row A/B key <-> artwork: "%s" dashed, best single tool dotted, '
                "mouse study A solid / study B dashed (asserted)", GLORI_LABEL)
    logger.info("row keys below their own panels: %s (in below each panel edge)",
                [round(b - lg.get_window_extent().y0 / fig.dpi, 2)
                 for lg, b in meta["keyed"]])
    logger.info("text objects %d; font sizes %s; smallest %.2f pt; overlaps %d",
                rep["n_text"], rep["fonts"], rep["min_fontsize"],
                len(rep["overlaps"]))
    for a, b in rep["overlaps"]:
        logger.warning("text overlap: %r / %r", a, b)
    assert rep["min_fontsize"] >= MIN_PT, \
        f"smallest text {rep['min_fontsize']:.2f} pt < {MIN_PT} pt"
    assert not rep["overlaps"], f"{len(rep['overlaps'])} overlapping text pairs"
    assert not rep["clipped"], f"text outside the canvas: {rep['clipped'][:6]}"
    assert rep["n_tool_labels"] == len(tools), \
        (f"{rep['n_tool_labels']} tool-name labels on the canvas; expected "
         f"{len(tools)} (row C, left margin of column 0 only)")
    assert not rep["tool_label_intrusion"], \
        (f"tool names printed over another panel: "
         f"{rep['tool_label_intrusion'][:6]}")
    logger.info("row C tool names: %d labels, printed once; intrusion %d",
                rep["n_tool_labels"], len(rep["tool_label_intrusion"]))
    logger.info("no text runs off the %.0f x %.0f px canvas",
                rep["canvas_px"][0], rep["canvas_px"][1])
    for i, (x0, x1, _, _) in enumerate(rep["legends"]):
        assert x0 >= -1.0 and x1 <= rep["canvas_px"][0] + 1.0, \
            (f"legend {i} spans {x0:.0f}..{x1:.0f} px on a "
             f"{rep['canvas_px'][0]:.0f} px canvas (clipped at the page edge)")
    logger.info("legends inside the canvas: %s",
                [[round(v, 1) for v in b] for b in rep["legends"]])

    FIG.mkdir(parents=True, exist_ok=True)
    stem = FIG / "Figure6_rev"
    # exact canvas: no tight bounding box, so 6.66 in stays 479.5 pt on paper
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches=None, pad_inches=0)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches=None,
                pad_inches=0)
    plt.close(fig)
    logger.info("wrote %s.pdf + %s.png", stem, stem)
    logger.info("output: %s", FIG)


if __name__ == "__main__":
    main()
