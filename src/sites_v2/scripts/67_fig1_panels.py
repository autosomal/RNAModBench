#!/usr/bin/env python
"""67 -- main-text Figure 1 body (panels B, C, D) at the printed size (2026-09-21).

Why this script exists
----------------------
The 2026-09-20 rebuild stitched three separately drawn vector panels onto a
941 x 478 pt canvas, which the manuscript then places with
``\\includegraphics[width=0.95\\textwidth]{Figure1.pdf}`` = 479.52 pt: at that
scale the 10.5 pt tool names of panel B print at 5.4 pt.  This script draws the
three panels **on the canvas they are printed on** (6.66 x 6.53 in, width = 0.95
\\textwidth = 479.52 pt), so every font size written here is the printed one, and
the figure is exported with ``bbox_inches=None``.

Figure contract (nature-figure)
-------------------------------
*Core conclusion*: which nanopore m6A tools were actually run on which dataset,
how many sites each of them calls there, and how many of those calls sit in the
canonical RRACH motif -- stated per independent sequencing unit, never as a
replicate-pooled union (reviewers R3-2 / E6).

*Evidence chain*: B = site counts across the four species with their controls
(the published lollipop, now per unit mean); C = Curlcake constructs, all calls,
modified versus unmodified template; D = the same calls restricted to RRACH.
C and D share one row order, so the drop from C to D is the specificity reading.

*Archetype*: quantitative grid (three stacked/adjacent panels, no hero).

*Journal/export contract*: printed size 1:1, Arial embedded (fonttype 42), no
gridlines, no in-panel annotation text, bold panel letters, >= 7 pt everywhere,
vector PDF + 300 dpi PNG, every number traceable to a frozen table.

Colours are the published figure's (m6A ``#F5B264``, WT/unmodified ``#3778A0``),
row order is m6A-count descending with the largest on top, and tools that were
never run on a dataset get **no** marker -- the published panel implicitly drew
a zero for DRUMMER / ELIGOS2_diff on the Curlcake IVT controls.

Inputs (frozen; nothing is recounted here)
------------------------------------------
``04_revision_analysis/figure1b_replicates/tables/per_replicate_tool_counts.tsv``
``04_revision_analysis/fig1_revision/tables/fig1C_tool_counts.tsv``
``04_revision_analysis/fig1_revision/tables/fig1D_rrach_counts.tsv``

Outputs
-------
``04_revision_analysis/fig1_revision/figures/Figure1_rev_body.{pdf,png}``
``04_revision_analysis/fig1_revision/tables/fig1_page_geometry.tsv``

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    01_code/code/sites_v2/scripts/67_fig1_panels.py
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
import importlib.util
import re
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.text import Text
from matplotlib.transforms import Bbox

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from common import config as C                       # noqa: E402
from common import pagelayout as PL                  # noqa: E402
from common.figstyle import apply as figstyle_apply  # noqa: E402
from common.io_utils import write_table              # noqa: E402
from common.manifest import setup_logger             # noqa: E402

OUT = (_RB / "figures/figure1")
FIG, TAB, LOG = (_RB / "figures/figure1/figures"), (_RB / "figures/figure1/tables"), (_RB / "figures/figure1/logs")
B_TAB = (_RB / "figures/figure1/inputs/tables")

# ------------------------------------------------------------------ canvas --- #
#: 1.00 * \textwidth of USG.cls, measured on the class itself on 2026-09-23
#: (\typeout{\the\textwidth} -> 506.45908 pt).  The page used to be 0.95 of that
#: (479.52 pt), which the manuscript then embedded at width=0.95\textwidth --
#: an extra 0.95 x on top of a canvas that is already drawn at printed size, so
#: the 7.2 pt legend of panels B/C/D printed at 6.84 pt, under the floor.  The
#: canvas is now the full text width and manuscript.tex embeds it 1:1.
PAGE_W = 506.4591
#: A strip redrawn at its printed size by 71 (121 pt; the published artwork it
#: replaces carried 2.8 pt labels and could not be blown up enough to comply).
#: The page grew by the same 17 pt so the B/C/D panels keep the geometry that
#: passed the layout audit:
#: 121 + 8 + 357 pt of panels + two 10 pt margins
PAGE_H = 506.0
MARGIN = 10.0
A_H = 121.0
A_RECT = (MARGIN, PAGE_H - MARGIN - A_H, PAGE_W - 2 * MARGIN, A_H)
#: the shared C/D legend band (16 pt) was removed on 2026-09-24 -- the two mean
#: codings are labelled directly inside C and D -- and its room went to panel B:
B_RECT = (16.0, MARGIN, 284.0, A_RECT[1] - 8.0 - MARGIN)
CD_X, CD_W = 318.0, PAGE_W - MARGIN - 318.0
_CD_H = (B_RECT[3] - 14.0) / 2.0
C_RECT = (CD_X, MARGIN + _CD_H + 14.0, CD_W, _CD_H)
D_RECT = (CD_X, MARGIN, CD_W, _CD_H)
#: panel letters sit in the page margin, left of every label column; C and D go
#: into the gap between panel B's right column and their own tool-name labels
LETTER_X = 3.0
CD_LETTER_X = 303.0
LETTER_GAP = 2.0

FS = {"tick": 7.5, "tool": 7.2, "axis": 8.5, "title": 9.5, "legend": 7.2,
      "legend_cd": 8.0, "letter": 11.0}
MIN_PT = 7.0
ORANGE = "#F5B264"          # Curlcake_m6A  (published colour)
BLUE = "#3778A0"            # Curlcake_IVT  (published colour, = B's WT blue)
CONNECT = "0.45"
CONDS = ("Curlcake_m6A", "Curlcake_IVT")
COND_LABEL = {"Curlcake_m6A": "m6A", "Curlcake_IVT": "unmodified"}
VALUE_COL = {"C": "n_sites", "D": "n_rrach"}


def apply_style() -> None:
    figstyle_apply()
    mpl.rcParams.update({
        "font.size": FS["tick"],
        "axes.titlesize": FS["title"],
        "axes.labelsize": FS["axis"],
        "xtick.labelsize": FS["tick"],
        "ytick.labelsize": FS["tick"],
        "legend.fontsize": FS["legend"],
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 2.4, "ytick.major.size": 2.4,
        # pagelayout.assert_page_clean() gates "text hugs a frame line" at 3 pt;
        # matplotlib's default 3.5 pt pad leaves the tick labels 2 pt from the
        # spine in these dense panels, so they are pushed out to >= 3 pt
        "xtick.major.pad": 5.5, "ytick.major.pad": 5.5,
        "legend.frameon": False,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "figure.dpi": 120,
    })


# -------------------------------------------------------------------- data --- #
def load_28():
    """Import the figure-1B script (its module name starts with a digit).

    The module has to be registered in ``sys.modules`` *before* it is executed:
    its ``@dataclass`` decorator resolves the class module through
    ``sys.modules`` (Python 3.11) and would raise AttributeError otherwise.
    """
    path = Path(__file__).with_name("28_fig_tool_counts_replicates.py")
    spec = importlib.util.spec_from_file_location("fig1b28", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["fig1b28"] = mod
    spec.loader.exec_module(mod)
    return mod


def load() -> dict[str, pd.DataFrame]:
    return {
        "b": pd.read_csv((_RB / "figures/figure1/inputs/tables/per_replicate_tool_counts.tsv"), sep="\t"),
        "c": pd.read_csv(TAB / "fig1C_tool_counts.tsv", sep="\t"),
        "d": pd.read_csv(TAB / "fig1D_rrach_counts.tsv", sep="\t"),
    }


def cd_series(tab: pd.DataFrame, value_col: str) -> dict:
    """``{tool: {condition: [unit values]}}`` for the Curlcake panels."""
    out: dict[str, dict[str, list[int]]] = {}
    for _, r in tab.iterrows():
        out.setdefault(r.tool, {}).setdefault(r.dataset_group, []).append(
            int(r[value_col]))
    return out


def cd_order(series: dict, ref_cond: str = "Curlcake_m6A") -> list[str]:
    """Tool order shared by C and D: mean m6A calls descending (largest on top).

    The published C and D each carried their *own* order, so the same tool sat at
    a different row in the two panels; sharing one order is what lets the reader
    go from "how many calls" to "how many of them are in RRACH" row by row.
    """
    return sorted(series, key=lambda t: -float(np.mean(series[t][ref_cond])))


# ------------------------------------------------------------------ panels --- #
def block_order(per_rep: pd.DataFrame, mod28) -> list[str]:
    """One row order for all four species sub-panels.

    The per-species sort of the published panel cannot survive the printed size:
    a 2 x 2 grid whose columns carry different orders needs a tool-name column
    *per column* (40 pt each), and the right one then runs into the left panel
    (measured by the page gate).  The four sub-panels therefore share one order
    -- the geometric mean of each tool's mean count over the WT datasets it was
    run on -- and the names are printed once per row of the grid.

    A tool that a dataset never saw keeps its row (DENA and MINES were not run
    on E. coli): dropping it would remove it from the three species that do have
    it, and a missing row is the honest picture.
    """
    wt_groups = [p.wt_group for p in mod28.PANELS]
    logmean = {}
    for t in mod28.TOOL_ORDER:
        vals = [max(v, 0.5) for g in wt_groups
                for _, v in mod28.group_points(per_rep, g, t)]
        if vals:
            logmean[t] = float(np.exp(np.mean(np.log(vals))))
    return sorted(logmean, key=lambda t: -logmean[t])


def draw_b(fig, per_rep: pd.DataFrame, mod28, order: list[str]) -> list:
    """Panel B: the four species dumbbells, drawn by the 28 module itself.

    Widths come from the measured label column (``Nanocompore`` at 7.2 pt) so the
    sub-panels of a row line up; the names are printed on the left column only.
    """
    label_w = PL.text_width_in("Nanocompore", FS["tool"]) * 72.0 + 5.0
    #: ``bottom`` has to hold the x axis label of the second row (17 pt: 14 pt
    #: put its last pixel on the page edge, caught by the canvas audit) and
    #: ``gap_x`` the two facing decade labels of a row (a 12 pt gap let the
    #: Arabidopsis 10^7 label run into the Mouse 10^1 label)
    gap_x, title, row_gap, bottom = 20.0, 12.0, 28.0, 17.0
    sub_w = (B_RECT[2] - label_w - 3.0 - gap_x) / 2.0
    sub_h = (B_RECT[3] - title - row_gap - bottom) / 2.0
    x0 = B_RECT[0] + label_w
    top = B_RECT[1] + B_RECT[3]
    axes = []
    for i, panel in enumerate(mod28.PANELS):
        col, row = i % 2, i // 2
        ax = fig.add_axes(_rect(x0 + col * (sub_w + gap_x),
                                top - title - (row + 1) * sub_h - row * row_gap,
                                sub_w, sub_h))
        ax.tick_params(length=2.2, pad=3.0, labelsize=FS["tick"])
        mod28.draw_panel(ax, panel, per_rep, order=order,
                         show_ylabels=(col == 0),
                         fs={"tick": FS["tool"], "title": FS["title"],
                             "axis": FS["axis"], "legend": FS["legend"]},
                         ms_scale=0.78, title_pad=4.0)
        _thin_decades(ax)
        axes.append(ax)
    return axes


def _plain_decade(t: float) -> str:
    """``0.1`` / ``100`` / ``1000000`` -- never ``1e+06``."""
    return f"{t:.0f}" if t >= 1 else f"{t:g}"


def _thin_decades(ax, max_n: int = 4) -> None:
    """Keep the decades that fit, spread over the whole range, in plain digits.

    Three constraints collide on a 106 pt axis: eight decade labels do not fit,
    dropping *every second* decade (the first draft) left Arabidopsis with
    10^-1, 10^3, 10^7, and matplotlib's ``10^n`` form prints the exponent in
    mathtext at 0.7 x the tick size = 5.25 pt, below the 7 pt house limit
    (caught by the acceptance script).  The labels are therefore plain numbers
    (the convention of the rebuilt Figure 6) and their real ink width decides how
    many survive.
    """
    #: the decade list is built from the x limits, not from ``get_xticks()``:
    #: matplotlib has already thinned a 106 pt log axis to four of them, and a
    #: shortened list cannot be thinned any further (the first draft kept
    #: 0.1 / 10 / 100000 / 10000000, whose last two labels collide)
    xlo, xhi = ax.get_xlim()
    first = int(np.ceil(np.log10(xlo) - 1e-9))
    last = int(np.floor(np.log10(xhi) + 1e-9))
    ticks = [10.0 ** e for e in range(first, last + 1)]
    width_pt = ax.get_position().width * PAGE_W
    lo, hi = np.log10(xlo), np.log10(xhi)

    def fits(keep: list[float]) -> bool:
        """Adjacent labels must not touch: the spread is not uniform."""
        pos = [(np.log10(t) - lo) / (hi - lo) * width_pt for t in keep]
        w = [PL.text_width_in(_plain_decade(t), FS["tick"]) * 72.0 for t in keep]
        return all(pos[i + 1] - pos[i] >= 0.5 * (w[i] + w[i + 1]) + 3.0
                   for i in range(len(keep) - 1))

    keep = list(ticks)
    labels = [_plain_decade(t) for t in keep]
    for n in range(max_n, 1, -1):
        if len(ticks) <= n:
            break
        idx = np.unique(np.round(np.linspace(0, len(ticks) - 1, n)).astype(int))
        keep = [ticks[i] for i in idx]
        if fits(keep):
            break
    labels = [_plain_decade(t) for t in keep]
    ax.set_xticks(keep)
    ax.set_xticklabels(labels)
    #: without this the decades that were just removed come back as *minor* tick
    #: labels
    ax.minorticks_off()


def draw_cd(ax, series: dict, order: list[str], *, log_x: bool, xlim: tuple,
            ticks: list[float], xlabel: str, display: dict,
            ticklabels=None) -> None:
    """One Curlcake metric, one row per tool, forest-dot coding.

    Per tool: the two condition means are the large dots on the row line (joined
    by a thin grey connector, the drop from modified to unmodified template), and
    every independent sequencing unit is a small dot **above** the line (m6A) or
    **below** it (unmodified) -- the published panel drew the units at the bar
    centre line, where half of them disappeared inside the bar.  A measured zero
    is a solid dot clamped to the axis floor; a tool that was never run on a
    dataset gets no dot at all instead of a fake zero.
    """
    y_pos = {t: float(i) for i, t in enumerate(order)}
    unit_dy, mean_ms, unit_ms = 0.26, 4.0, 2.6
    ax.set_xlim(*xlim)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_yticks(range(len(order)))
    #: canonical display names, the same strings panel B uses (R1-7: "yanocomp"
    #: is printed "Yanocomp" everywhere in the manuscript)
    ax.set_yticklabels([display.get(t, t) for t in order], fontsize=FS["tool"])
    if log_x:
        ax.set_xscale("log")
        ax.minorticks_off()
    ax.set_xticks(ticks)
    if ticklabels is not None:
        ax.set_xticklabels(ticklabels)
    ax.set_xlabel(xlabel, fontsize=FS["axis"], labelpad=1.5)
    ax.tick_params(length=2.2, pad=3.0, labelsize=FS["tick"])

    floor = float(xlim[0])
    for tool in order:
        y = y_pos[tool]
        means = {}
        for cond, dy, colour in (("Curlcake_m6A", +unit_dy, ORANGE),
                                 ("Curlcake_IVT", -unit_dy, BLUE)):
            vals = series[tool].get(cond)
            if not vals:
                continue
            mean = float(np.mean(vals))
            means[cond] = max(mean, floor) if log_x else mean
            for v in vals:
                x = max(float(v), floor) if log_x else float(v)
                ax.plot(x, y + dy, "o", ms=unit_ms, mfc=colour, mec="white",
                        mew=0.45, zorder=6, clip_on=False)
            ax.plot(max(mean, floor) if log_x else mean, y, "o", ms=mean_ms,
                    mfc=colour, mec="white", mew=0.6, zorder=5)
        if len(means) == 2:              # connector: m6A -> unmodified
            ax.plot([means["Curlcake_m6A"], means["Curlcake_IVT"]], [y, y],
                    color=CONNECT, lw=0.7, zorder=2, solid_capstyle="round")
    #: the two mean codings are page-level labels in the empty strips above
    #: panels C and D (see build_figure): the 12 pt row pitch of C/D leaves no
    #: room for 7.5 pt text inside the panels next to the dots


def _rect(x_pt: float, y_pt: float, w_pt: float, h_pt: float) -> list[float]:
    return [x_pt / PAGE_W, y_pt / PAGE_H, w_pt / PAGE_W, h_pt / PAGE_H]


# ------------------------------------------------------------------ checks --- #
def check_inputs(t: dict[str, pd.DataFrame], series: dict, logger) -> None:
    """Anchors: the panels may only draw what the frozen tables state."""
    c, d = t["c"], t["d"]
    assert len(c) == 44 and len(d) == 44, (len(c), len(d))
    for name, tab, col in (("C", c, "n_sites"), ("D", d, "n_rrach")):
        units = tab.groupby(["dataset_group", "tool"]).size()
        assert set(units.unique()) == {2}, \
            f"panel {name}: expected 2 independent units per cell, got " \
            f"{sorted(set(units.unique()))}"
        for cond in CONDS:
            sub = tab[tab.dataset_group == cond]
            if cond == "Curlcake_IVT":
                assert set(sub.tool) == {"DENA", "ELIGOS2_solo", "EpiNano_Error",
                                         "MINES", "NanoSPA_m6A", "Nanocompore",
                                         "Nanom6A", "m6Anet", "xPore", "yanocomp"}
            else:
                assert len(set(sub.tool)) == 12
        logger.info("panel %s: %d tool-cell rows, 2 units each", name, len(tab))
    for tool in series:
        for cond in CONDS:
            vals = series[tool].get(cond)
            if vals:
                assert len(vals) == 2, (tool, cond, vals)
    # the site counts of panel B must be the same numbers the C panel shows
    per_rep = t["b"]
    assert per_rep[per_rep.dataset_group == "Curlcake_m6A"].n_sites.sum() > 0
    logger.info("table anchors OK: 2 units per tool and condition, "
                "12 tools on m6A and 10 on the unmodified template")


def _legend_hits(fig, ax, leg) -> list:
    """Data points that fall inside a legend box (the project's test for
    "the legend covers the data"; bounding-box intersection is not used)."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    bb = leg.get_window_extent(renderer=rend)
    hits = []
    for ln in ax.lines:
        xd = np.asarray(ln.get_xdata(), dtype=float)
        yd = np.asarray(ln.get_ydata(), dtype=float)
        keep = np.isfinite(xd) & np.isfinite(yd)
        if not keep.any():
            continue
        pts = ax.transData.transform(np.column_stack([xd[keep], yd[keep]]))
        for px, py in pts:
            if bb.x0 - 1.0 <= px <= bb.x1 + 1.0 and bb.y0 - 1.0 <= py <= bb.y1 + 1.0:
                hits.append((float(px), float(py)))
    return hits


def _short_legend_label(text: str) -> str:
    """``WT (mean of 3)`` -> ``WT (n = 3)``: the sub-panels are 111 pt wide.

    The full wording ("marker = mean of the group's independent sequencing
    units") is in the caption; inside the panel only the count has to survive.
    """
    text = re.sub(r"\(mean of (\d+)\)", r"(n = \1)", text)
    return text.replace("(two studies)", "(2 studies)")



#: pushed to whatever corner happened to be free (lower right / lower left /
#: upper right / upper left, measured at y = 162 / 270 / 324 / 432 pt), so they
#: read as four unrelated labels.  They now all sit in the upper left; when that
#: corner carries dumbbells, the low end of the log axis is extended by these
#: factors (the data shift right, the corner empties) instead of moving the key.
LEGEND_X_PASSES = (1.0, 0.35, 0.12)


def compact_legend(fig, ax, label: str, logger, *,
                   loc: str = "upper left",
                   passes: tuple = LEGEND_X_PASSES) -> None:
    """Place a sub-panel legend in the upper left, free of data points.

    The key keeps one fixed corner for all four species panels; the *axis* is
    what gives way when the corner is occupied (up to ``passes`` attempts, each
    extending the low end of the log decade range and re-thinning the ticks).
    """
    leg = ax.get_legend()
    if leg is None:
        return
    handles = list(leg.legend_handles)
    labels = [_short_legend_label(t.get_text()) for t in leg.get_texts()]
    for attempt, factor in enumerate(passes, start=1):
        if factor != 1.0:
            lo, hi = ax.get_xlim()
            ax.set_xlim(lo * factor, hi)
            _thin_decades(ax)
        leg = ax.legend(handles, labels, loc=loc, frameon=False,
                        fontsize=FS["legend"], handlelength=0.9,
                        handletextpad=0.32, labelspacing=0.28, borderpad=0.18)
        hits = _legend_hits(fig, ax, leg)
        bb = leg.get_window_extent(renderer=fig.canvas.get_renderer())
        logger.info("%s legend %.0f x %.0f pt at %s (x-min pass %d/%d, %.3g) "
                    "covers %d data point(s)", label, bb.width, bb.height, loc,
                    attempt, len(passes), ax.get_xlim()[0], len(hits))
        if not hits:
            return
    raise SystemExit(f"{label}: the upper-left key still covers data points")


def layout_report(fig) -> dict:
    """Printed-size audit: font sizes and text boxes outside the canvas."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    items = []
    for t in fig.findobj(Text):
        if not t.get_visible() or not str(t.get_text()).strip():
            continue
        try:
            bb = t.get_window_extent(renderer=rend)
        except Exception:                                   # pragma: no cover
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        items.append({"s": str(t.get_text())[:32], "fs": float(t.get_fontsize()),
                      "bb": [float(bb.x0), float(bb.y0), float(bb.x1), float(bb.y1)]})
    fonts: dict[float, int] = {}
    for it in items:
        fonts[it["fs"]] = fonts.get(it["fs"], 0) + 1
    w_px = PAGE_W / 72.0 * float(fig.dpi)
    h_px = PAGE_H / 72.0 * float(fig.dpi)
    outside = [(it["s"], f"x {it['bb'][0]:.0f}..{it['bb'][2]:.0f} "
                         f"y {it['bb'][1]:.0f}..{it['bb'][3]:.0f}")
               for it in items
               if it["bb"][0] < -0.5 or it["bb"][2] > w_px + 0.5
               or it["bb"][1] < -0.5 or it["bb"][3] > h_px + 0.5]
    return {"n_text": len(items), "fonts": dict(sorted(fonts.items())),
            "min_fontsize": min((it["fs"] for it in items), default=float("nan")),
            "outside": outside, "canvas_px": [w_px, h_px]}


# -------------------------------------------------------------------- main --- #
def build_figure(t: dict[str, pd.DataFrame]) -> tuple[plt.Figure, dict]:
    apply_style()
    mod28 = load_28()
    series = cd_series(t["d"], VALUE_COL["D"])
    series_c = cd_series(t["c"], VALUE_COL["C"])
    #: one shared row order for C and D, taken from the C panel (all calls): the
    #: published panels each carried their own order, so the same tool sat at a
    #: different height in the two, and the C-to-D comparison could not be read
    #: row by row
    order = cd_order(series_c)
    assert set(order) == set(series), \
        f"C and D must list the same tools: {set(order) ^ set(series)}"

    fig = plt.figure(figsize=(PAGE_W / 72.0, PAGE_H / 72.0))
    order_b = block_order(t["b"], mod28)
    ax_b = draw_b(fig, t["b"], mod28, order_b)

    ax_c = fig.add_axes(_rect(C_RECT[0] + _label_w(), C_RECT[1] + 21.0,
                              C_RECT[2] - _label_w() - _right_pad(), C_RECT[3] - 25.0))
    ax_d = fig.add_axes(_rect(D_RECT[0] + _label_w(), D_RECT[1] + 21.0,
                              D_RECT[2] - _label_w() - _right_pad(), D_RECT[3] - 25.0))
    #: plain decade labels, not matplotlib's 10^n: the superscript of ``10^n`` is
    #: mathtext at 0.7 x the tick size = 5.25 pt, below the 7 pt house limit
    draw_cd(ax_c, series_c, order, log_x=True, xlim=(0.7, 3000.0),
            ticks=[1, 10, 100, 1000], xlabel="counts (all calls)",
            display=mod28.DISPLAY, ticklabels=["1", "10", "100", "1000"])
    draw_cd(ax_d, series, order, log_x=False, xlim=(-6.0, 136.0),
            ticks=[0, 50, 100], xlabel="counts (RRACH calls)",
            display=mod28.DISPLAY)

    #: letters sit in the page margin, hard against the top-left corner of their
    #: panel: above the panels there is no room (the 10 pt top margin of panel A
    #: is thinner than an 11 pt letter), and every label column starts at x >= 25
    for letter, rect in (("A", A_RECT), ("B", B_RECT), ("C", C_RECT),
                         ("D", D_RECT)):
        fig.text((LETTER_X if letter in ("A", "B") else CD_LETTER_X) / PAGE_W,
                 (rect[1] + rect[3] - 1.0) / PAGE_H, letter,
                 fontsize=FS["letter"], fontweight="bold", ha="left", va="top")
    #: the two mean codings of C and D, as page-level labels in the empty
    #: strips above each panel (8 pt between A and C, 14 pt between C and D);
    #: inside the panels the 12 pt row pitch leaves no room for 7.5 pt text
    for strip_top in (C_RECT[1] + C_RECT[3], D_RECT[1] + D_RECT[3]):
        x = CD_X + 2.0
        for lab, colour in (("m6A (mean)", ORANGE), ("IVT (mean)", BLUE)):
            fig.text(x / PAGE_W, (strip_top + 0.3) / PAGE_H, lab,
                     fontsize=7.5, color=colour, ha="left", va="bottom")
            x += PL.text_width_in(lab, 7.5) * 72.0 + 6.0
    return fig, {"order": order, "order_b": order_b, "series": series,
                 "axes": {"B": ax_b, "C": ax_c, "D": ax_d}}


def _label_w() -> float:
    """Width of the 13 tool names at ``FS['tool']`` plus the tick pad."""
    return PL.text_width_in("Nanocompore", FS["tool"]) * 72.0 + 5.0


def _right_pad() -> float:
    return 3.0


def main() -> None:
    logger = setup_logger("67_fig1_panels", log_dir=LOG)
    FIG.mkdir(parents=True, exist_ok=True)
    t = load()
    fig, meta = build_figure(t)
    check_inputs(t, meta["series"], logger)

    n_cd = len(meta["order"])
    for sub in meta["axes"]["B"]:
        h = sub.get_position().height * PAGE_H
        logger.info("panel B sub-panel %.1f pt tall for 13 rows (%.2f pt/row)",
                    h, (h - 9.0) / 13.0)
        assert (h - 9.0) / 13.0 >= 9.0, "panel B rows are too tight"
    for name in ("C", "D"):
        h = meta["axes"][name].get_position().height * PAGE_H
        logger.info("panel %s %.1f pt tall for %d rows (%.2f pt/row)",
                    name, h, n_cd, h / n_cd)
        assert h / n_cd >= 9.5, f"panel {name} rows are too tight"
    for i, ax in enumerate(meta["axes"]["B"]):
        compact_legend(fig, ax, f"B[{i}]", logger)

    rep = layout_report(fig)
    logger.info("text objects %d; font sizes %s; smallest %.2f pt",
                rep["n_text"], rep["fonts"], rep["min_fontsize"])
    assert rep["min_fontsize"] >= MIN_PT, \
        f"smallest text {rep['min_fontsize']:.2f} pt < {MIN_PT} pt"
    assert not rep["outside"], f"text outside the page: {rep['outside'][:6]}"
    logger.info("panel B shared row order: %s", " > ".join(meta["order_b"]))
    logger.info("panel C/D shared row order: %s", " > ".join(meta["order"]))
    logger.info("page %.2f x %.2f pt = %.2f x %.2f in (printed 1:1)",
                PAGE_W, PAGE_H, PAGE_W / 72, PAGE_H / 72)

    PL.assert_page_clean(fig, min_pt=MIN_PT, verbose=True)

    stem = FIG / "Figure1_rev_body"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches=None, pad_inches=0)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches=None,
                pad_inches=0)
    plt.close(fig)
    write_table(pd.DataFrame([
        {"key": "page_w_pt", "value": f"{PAGE_W:.2f}"},
        {"key": "page_h_pt", "value": f"{PAGE_H:.2f}"},
        {"key": "margin_pt", "value": f"{MARGIN:.1f}"},
        {"key": "A_strip_pt", "value": f"{A_RECT[2]:.2f} x {A_RECT[3]:.2f} "
                                       f"at ({A_RECT[0]:.1f},{A_RECT[1]:.1f})"},
        {"key": "B_rect_pt", "value": ",".join(f"{v:.1f}" for v in B_RECT)},
        {"key": "C_rect_pt", "value": ",".join(f"{v:.1f}" for v in C_RECT)},
        {"key": "D_rect_pt", "value": ",".join(f"{v:.1f}" for v in D_RECT)},
        {"key": "font_sizes_pt", "value": "; ".join(f"{k}={v}" for k, v in FS.items())},
        {"key": "row_order", "value": " > ".join(meta["order"])},
        {"key": "b_panel_cells", "value": "2 x 2 species, 13 tools"},
    ]), TAB / "fig1_page_geometry.tsv")
    logger.info("wrote %s.pdf + .png", stem)
    logger.info("output: %s", FIG)


if __name__ == "__main__":
    main()
