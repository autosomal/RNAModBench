#!/usr/bin/env python
"""42 -- Figure S4 revision (panels A--D): one panel at a time, then assemble.

The submitted Figure S4 was one A4 portrait page whose rows each split into the
three species columns (Arabidopsis | Mouse | Human):

    A  calls per tool: WT-common counts (log), purified counts (log) and the
       purified/WT ratio, three stacked sub-rows per species
    B  the matching-window sweep (0--50 bp), upper facet: PPV vs. GLORI
    C  the same sweep, lower facet: the exact-nucleotide fraction
    D  PPV (2 bp) of the full WT-common call set vs the purified subset
       (lettered C until 2026-09-21)

Panel identity is the one the manuscript cites (S4A counts, S4B window-sweep
PPV, S4C exact-nucleotide fraction, S4D precision); the three panels the
reviewers asked for on top of the original figure (R3-7: GLORI overlap, DRACH
fraction, stoichiometry of the three site groups) are the separate
Supplementary Figure S10 (``48_figS5_validation.py``).

The PR-AUC row of the submitted page (PR-AUC vs. the matching window) was
0.017--0.106, i.e. it repeats the window dependence of PPV that the window sweep
already carries, and nothing in the manuscript cites it.  Its frozen tables
(``figS4_auprc_per_unit.tsv`` / ``figS4_auprc_summary.tsv``) and the anchors
that guard them stay in place (see :func:`write_anchors`), but no panel is drawn
and the row is gone from the page.  **Note** that this dropped row is *not* the
panel D of 2026-09-21: D is the purified-site comparison, which was lettered C
before the window sweep was split.

The window sweep **became two stacked facets on 2026-09-21** ("B panel is too dense 
 split into two facets "), and each facet **carries a panel letter of its own**
(" each cell keeps its own letter "): the upper facet is **B** (PPV vs. GLORI), the lower one
**C** (the exact-nucleotide fraction that used to be a dashed overlay on the
same 0--1 axis).  Both facets keep the per-unit thin lines and the group mean of
their own quantity, both share one x axis, and the purified-site panel moved
from C to **D**; the manuscript and the response letter were re-lettered
accordingly (S4A / S4B / S4C / S4D).  The height for the second facet comes from
panel A (4.00 -> 3.52 in; its counts sub-rows still hold three power-of-ten
labels), so the sub-panels stay close to square -- B/C 2.16 x 1.83 in, D
2.16 x 1.90 in -- which is what " BC" asked for:
0.05 + 3.52 + 0.08 + 0.36 + 0.08 + 4.80 (BC) + 0.08 + 2.62 (D) + 0.11 = 11.70 in.

Layout (2026-09-21, third rework -- "draw each panel, then assemble")
--------------------------------------------------------------------
Two earlier reworks failed because they fit four rows into one page *and* one
canvas: the panels ended up 3.2 x 0.85 in letterbox strips, and the sub-row
titles of panel A shared a band with the tick labels of the row above, so the
text printed over the axis line.  This version inverts the workflow:

* every panel (A--D) is drawn on its **own canvas
  at its final print size** (``common/panelpage.py``), with its own page margin
  -- panel letter, y axis label, tick labels -- so its internal layout is
  never squeezed by the page budget;
* each canvas passes ``pagelayout.assert_page_clean`` on its own (text-text,
  text-in-foreign-panel, **text-on-a-spine**, grid, minimum font size);
* the A4 portrait page is then assembled by translating the panel PDFs onto a
  blank page with pypdf (vector preserved); the placement rectangles are
  asserted to be inside the page and pairwise disjoint, so a cross-panel
  collision is impossible by construction;
* the 300 dpi PNG is rasterised from the assembled PDF (``pdftoppm``), so the
  bitmap cannot drift from the vector page.

Panel A gets dedicated title bands (0.30 in) between the sub-rows: a band has to
hold the title (0.03--0.19 in above its own axes) *and* the lowest y tick label
of the sub-row above (which reaches 0.08 in below that axes), which is exactly
the overlap that was visible before.

Legend placement (2026-09-27, sixth rework -- "BCBC"): the 13
tool colours are drawn once, **inside the BC canvas, in the band between the two
facets** (B above, C below), which is where they are read and where nothing is
drawn behind them.  Until then they were a page-level stripe between panels A
and B (``render_tool_stripe``, now unused).  The two keys that belong to a single
panel stay on their own line inside that panel's canvas: the species /
mouse-study / mean key at the top of A, the line-type key (which explains the
two C entries as well) at the top of B.  ``pagelayout.assert_legend_clear``
re-measures each in-panel legend against the drawn vertices and aborts if a key
hides data.

Contract
--------
* claim: the purified subset is only modestly better (descriptive -- R3-7
  circularity) and the apparent gain of a wide matching window is dominated by
  localisation tolerance rather than detection;
* evidence chain: A counts/ratio -> B window sweep (two facets: PPV above,
  exact-nucleotide fraction below) -> C WT vs purified;
* archetype: quantitative grid (three block-letter rows x three species columns;
  B is a two-facet row sharing one x axis, and every species sub-panel of A, B
  and C keeps its own y tick labels);
* mouse: the two cross-study datasets (SRP166020 = study A, SRP357195 = study
  B) are drawn separately in every panel and are NEVER averaged;
* export: A4 portrait, vector PDF with embedded Arial (``pdf.fonttype`` 42)
  plus a 300 dpi PNG, exact page size (no ``bbox_inches='tight'``).  House
  rules: Arial, no grid, no annotation text inside a panel, bold panel letters
  in the page margin, ticks/axis labels >= 8.5 pt.

Outputs
-------
``figures/figureS4/figures/panels/FigureS4_A.pdf`` (A),
``FigureS4_BC.pdf`` (B and C on one canvas), ``FigureS4_D.pdf`` and
``FigureS4_legend_tools.pdf``
``figures/figureS4/figures/FigureS4_rev.{pdf,png}``
``figures/figureS4/logs/42_figS4_anchors.tsv``

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/figures/figureS4/src/42_figS4_figure.py
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
import argparse
import sys
import time
from datetime import datetime
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.legend import Legend
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator

HERE = Path(__file__).resolve()
sys.path.insert(0, str(_RB / "src/harmonisation"))

from common import pagelayout, panelpage  # noqa: E402
from common.config import SITES_ROOT, TABLE_DIR, WINDOWS  # noqa: E402
from common.figstyle import apply  # noqa: E402

TAB = (_RB / "figures/figureS4/tables")
FIG = (_RB / "figures/figureS4/figures")
PANEL_DIR = FIG / "panels"
LOG = (_RB / "figures/figureS4/logs")

#: the page of the submitted sup4.pdf (and of every SI figure page)
#: 2026-09-23: A4 width, 0.36 in taller than A4 -- the upright tool names of
#: panel A and panel D's own key line added 0.42 in to the budget.  The SI page
#: scales this canvas by 0.80 (was 0.82), which keeps every label above 7 pt.
PAGE_IN = (pagelayout.A4_PORTRAIT_IN[0], 12.16)

# --------------------------------------------------------------------------- #
# geometry, in inches (the canvas of a panel *is* its final print size)
# --------------------------------------------------------------------------- #
#: every panel spans the full page width; the left strip holds the panel letter
#: (0.02--0.16 in), the y axis label (0.16--0.45) and the y tick labels
#: (0.46--0.70), so the three species columns start at LEFT_IN.  The column gap
#: has to hold the *neighbouring* column's y tick labels (0.26 in at 8 pt) plus
#: a visible clearance: at 0.38 in the labels of column 2 sat 1.4 pt from the

#: leaves 0.14 in and still buys the 45-degree tool names a 0.166 in slot.
PANEL_W_IN = PAGE_IN[0]
LEFT_IN = 0.74
RIGHT_IN = 0.06
COL_GAP_IN = 0.50
NCOL = 3
COL_W_IN = (PANEL_W_IN - LEFT_IN - RIGHT_IN - (NCOL - 1) * COL_GAP_IN) / NCOL

#: panel A: the three sub-rows and their dedicated title bands.  Compressed on
#: 2026-09-21 (4.00 -> 3.52 in) to pay for the second facet of the window sweep
#: (panels B/C): the counts sub-rows keep the three power-of-ten labels (0.50 in
#: leaves 12.4 pt between two 8 pt labels) and the ratio row keeps its 0.5 tick.
#: A_H_IN = 0.82 + 0.34 + (0.27 + 0.50) * 2 + 0.52 + 0.30 = 3.52 in
A_LABEL_IN = 1.06     # upright tool names: 5 pt pad + the 0.99 in of
                      # "NanoSPA_m6A" at 9.5 pt + ink
A_RATIO_IN = 0.34
A_BAND_IN = 0.27      # title band: title + the upper row's lowest tick label
A_SUB_IN = 0.50
A_HEAD_IN = 0.52      # species name (14 pt) on its own line + "WT counts"
#: each sub-row title starts at the top edge of that sub-row, i.e. next to the


A_TITLE_INDENT_IN = 0.07
#: the species / mouse-study / mean key lives on its own line inside the A
#: canvas (the ratio sub-rows are full to 1.0, so no sub-row has room for it)
A_KEY_IN = 0.30
A_H_IN = (A_LABEL_IN + A_RATIO_IN + A_BAND_IN + A_SUB_IN + A_BAND_IN
          + A_SUB_IN + A_HEAD_IN + A_KEY_IN)

#: the window sweep is split into **two stacked facets, one letter each**

#: upper facet is **B** = the PPV vs. GLORI curve, the lower one **C** = the
#: exact-nucleotide fraction; both carry the per-unit thin lines and the group
#: mean of their own quantity, and both share one x axis (only the lower facet
#: prints the window ticks and the axis name).  They stay on ONE canvas so each

#: in the page margin, each next to the top edge of its own facet.
#: 0.05 + 3.52 + 0.08 + 0.36 + 0.08 + 4.80 (BC) + 0.08 + 2.62 (D) + 0.11 = 11.70


#: the line-type key (Mean of units / Single unit / Study A / Study B), so the
#: band that used to sit above B is now empty and only keeps a hair of clearance
B_KEY_IN = 0.06       # clearance only: both keys live in the B/C band now
B_TITLE_IN = 0.28     # species column titles, printed once for the block
B_FACET_IN = 1.83     # one facet: upper = PPV, lower = exact fraction
#: the band between the facets carries the line-type key (one row) above the
#: tool key (two rows of seven entries); its height covers the clearance the "0"
#: tick label of the upper facet needs (0.073 in below its spine plus 3 pt), the
#: two keys and the top tick label of the lower facet
B_GAP_IN = 0.70
B_BOTTOM_IN = 0.44    # window ticks + "Matching window (bp)"
B_H_IN = (B_KEY_IN + B_TITLE_IN + 2 * B_FACET_IN + B_GAP_IN + B_BOTTOM_IN)

#: panel D (the purified-site comparison, lettered C until 2026-09-21): title
#: band, drawing area, bottom band (2.16 x 1.90 in sub-panels)
D_TITLE_IN = 0.46     # species titles + panel D's own key line
D_AXES_IN = 1.90
D_BOTTOM_IN = 0.44
D_H_IN = D_TITLE_IN + D_AXES_IN + D_BOTTOM_IN

#: the y axis name of the lower B facet; the full definition (fraction of the
#: window-matched calls that sit on the exact reference nucleotide, i.e. the
#: main-text quantity of Fig. 5D) lives in the legend text
B_EXACT_LABEL = "Exact-nucleotide fraction"

#: the 13 tool colours are explained once, in a centred stripe between panels A

#: full page width, no panel slot, no data behind it
TOOL_STRIPE_IN = 0.36

#: page assembly: top/bottom margin and the gaps between the stacked pieces
TOP_IN = 0.05
BOTTOM_IN = 0.11
GAP_IN = 0.08

TOOL_ORDER = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
              "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
              "NanoSPA_m6A", "xPore", "yanocomp"]
DISPLAY = {"yanocomp": "Yanocomp"}
TOOL_COLOR = {t: matplotlib.colormaps["tab20"].colors[i % 20]
              for i, t in enumerate(TOOL_ORDER)}

#: colours of the independent units (also used for the two mouse studies)
GROUP_COLOR = {"Arabidopsis": "#4b81b8", "Mouse": "#e8a76b", "HeLa": "#7aa25c"}
#: column titles follow the submitted figure / Fig. 5 ("Human", not "HeLa")
COLUMN_TITLE = {"Arabidopsis": "Arabidopsis", "Mouse": "Mouse", "HeLa": "Human"}

#: (column key, species label used in the tables, pairs in marker order).
#: The mouse column holds one pair per cross-study dataset; the two studies are
#: kept separate everywhere and are never averaged (R3-2 / E6).
COLUMNS: list[tuple[str, str, list[str]]] = [
    ("Arabidopsis", "Arabidopsis", ["Ath_rep1", "Ath_rep2", "Ath_rep3"]),
    ("Mouse", "Mouse", ["Mouse_studyA", "Mouse_studyB"]),
    ("HeLa", "Human", ["HeLa_rep1", "HeLa_rep2", "HeLa_rep3"]),
]
#: samples of the frozen localisation curve that belong to each column
CURVE_SAMPLES = {
    "Arabidopsis": ["Arabidopsis_WT_rep1", "Arabidopsis_WT_rep2",
                    "Arabidopsis_WT_rep3"],
    "Mouse": ["mESCs_Mettl3_WT", "mES_WT"],
    "HeLa": ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"],
}
#: mouse study -> (curve sample, line style, marker filled?)
MOUSE_STUDIES = [("Mouse_studyA", "mESCs_Mettl3_WT", "-", True),
                 ("Mouse_studyB", "mES_WT", "--", False)]

#: key line at the top of panel A: the species / mouse-study / mean semantics of
#: the dots (the species colours themselves are named by the column titles)
A_KEY: list[tuple[str, dict]] = [
    ("Arabidopsis (3 pairs)", dict(color=GROUP_COLOR["Arabidopsis"], marker="o",
                                   ls="none", ms=5)),
    ("Human (3 pairs)", dict(color=GROUP_COLOR["HeLa"], marker="o", ls="none",
                             ms=5)),
    ("Mouse study A", dict(color=GROUP_COLOR["Mouse"], marker="o", ls="none",
                           ms=5)),
    ("Mouse study B", dict(color="white", marker="o", ls="none", ms=5,
                           markeredgecolor=GROUP_COLOR["Mouse"])),
    ("Mean across pairs", dict(color="black", lw=2.0)),
]

#: line-type key of the B/C canvas, drawn on its own line at the top of it: the
#: two layers every facet draws (thin = single unit, thick = group mean), the
#: two mouse study styles, and the two panel-D marks.  The quantity of each
#: facet is named by its own y axis label, so the old "Exact nucleotide" entry
#: is gone: since 2026-09-21 the exact fraction is a facet of its own instead
#: of a dashed overlay.
LINE_KEY_B: list[tuple[str, dict]] = [
    ("Mean of units", dict(color="black", lw=2.0)),
    ("Single unit", dict(color="0.60", lw=0.8)),
    ("Study A", dict(color=GROUP_COLOR["Mouse"], lw=2.0)),
    ("Study B", dict(color=GROUP_COLOR["Mouse"], lw=2.0, ls="--")),
]

#: panel D's own key line: it used to be the last two entries of the B/C key,
#: i.e. printed on a canvas two panels away from the marks it explains
D_KEY: list[tuple[str, dict]] = [
    ("Per-tool pair", dict(color="0.55", lw=0.9)),
    ("Mean +- SD", dict(color="black", marker="o", ls="none", ms=4)),
]

#: the 13 tool colours, drawn in the centred page stripe (A -> B/C canvas)
TOOL_KEY: list[tuple[str, dict]] = [
    (DISPLAY.get(t, t), dict(color=TOOL_COLOR[t], lw=2.4)) for t in TOOL_ORDER
]

#: font sizes -- "one notch up" scale (2026-09-21).  Figure S10 uses a *larger*
#: scale of its own (its panels are taller, so the same point size reads
#: smaller); this dict is the S4 scale only.
#: 2026-09-23: tick_small / legend 8.0 -> 9.0 (at the SI scale of 0.82 the old
#: 8.0 printed at 6.6 pt).  The 45-degree tool names stay at 9.5 pt: 10 pt made
#: the neighbouring names touch in the 0.171 in slot (page gate aborts).
FS = {"tick": 9.5, "tick_small": 9.0, "axis": 11.5, "subrow": 9.5,
      "column": 14.0, "letter": 20.0, "legend": 9.0}

#: rotation of the 13 tool names on the panel-A x axis.  45 degrees is the
#: angle of the submitted figure; the gate measures the printed ink, and at
#: a 0.171 in slot the two neighbouring tilted names keep a ~1.5 pt gap
#: (see pagelayout self-test).  ``--rotation 90`` renders the upright variant
#: 2026-09-23: the tilted variant cannot carry a 7 pt floor -- two neighbouring
#: 45-degree names sit 11.3 pt apart, so the rotated text height (and with it the
#: own size) is capped at ~8 pt, i.e. 6.6 pt on the printed SI page.  The upright
#: variant is pitch-limited at 11.3 pt instead of 8 pt, so the names can be set
#: at 9.5 pt and stay clear of each other (page gate verifies).
TICK_ROTATION = 90


# --------------------------------------------------------------------------- #
# canvas helpers
# --------------------------------------------------------------------------- #
def col_x_in(ci: int) -> float:
    """Left edge of species column ``ci`` inside a panel canvas (inches)."""
    return LEFT_IN + ci * (COL_W_IN + COL_GAP_IN)


def col_center_in(ci: int) -> float:
    return col_x_in(ci) + COL_W_IN / 2.0


def _rect(height_in: float, *, x_in: float, y_in: float, w_in: float,
          h_in: float) -> list[float]:
    """``[left, bottom, width, height]`` figure fractions for ``add_axes``."""
    return [x_in / PANEL_W_IN, y_in / height_in, w_in / PANEL_W_IN,
            h_in / height_in]


def _title(fig: plt.Figure, height_in: float, x_in: float, y_in: float,
           text: str, *, ha: str = "left", fontsize: float | None = None,
           weight: str = "normal") -> None:
    """Panel text placed in its own band (never overlapping an axis line)."""
    fig.text(x_in / PANEL_W_IN, y_in / height_in, text, ha=ha, va="bottom",
             fontsize=FS["subrow"] if fontsize is None else fontsize,
             fontweight=weight)


def panel_letter(fig: plt.Figure, letter: str, *, y_in: float) -> None:
    """Bold panel letter in the panel's own left margin."""
    pagelayout.margin_letter(fig, letter, y=y_in / fig.get_figheight(),
                             x=0.02 / PANEL_W_IN, fontsize=FS["letter"])


def tool_ticks(ax: plt.Axes, rotation: int | None = None) -> None:
    """Full tool names on the bottom sub-row of panel A.

    The submitted figure tilted these names by 45 degrees (``ha="right"``, so
    the name runs down-left from its tick); 90 degrees is kept as the upright
    variant.  ``rotation_mode`` must stay at its default: with ``"anchor"`` the
    alignment is applied to the *unrotated* string, which centres the rotated
    name on the tick and pushes half of it up through the axis line (that was
    the visible "text over the frame" defect of panel A).

    The pad has to clear the same axis' lowest y tick label (``10^0``, which
    reaches ~5 pt below the axis) as well as the 3.5 pt tick marks, hence 8 pt
    for the tilted and 5 pt for the upright variant.
    """
    rotation = TICK_ROTATION if rotation is None else rotation
    ax.set_xticks(range(len(TOOL_ORDER)))
    ax.set_xticklabels([DISPLAY.get(t, t) for t in TOOL_ORDER],
                       rotation=rotation,
                       ha="center" if rotation == 90 else "right",
                       va="top", fontsize=FS["tick_small"])
    ax.tick_params(axis="x", pad=5.0 if rotation == 90 else 8.0)
    ax.set_xlim(-0.65, len(TOOL_ORDER) - 0.35)


def pair_rows(df: pd.DataFrame, tool: str,
              pairs: list[str]) -> list[pd.Series]:
    """Rows of ``df`` for one tool, one row per pair, in pair order."""
    sub = df[(df["tool"] == tool) & (df["pair_id"].isin(pairs))]
    out: list[pd.Series] = []
    for pair in pairs:
        hit = sub[sub["pair_id"] == pair]
        if not hit.empty:
            out.append(hit.iloc[0])
    return out


def pair_markers(ax: plt.Axes, x: float, values: list[float], filled: list[bool],
                 color: str, log_floor: float | None = None) -> None:
    """One marker per independent pair, slightly offset around ``x``."""
    if not values:
        return
    offsets = np.linspace(-0.20, 0.20, len(values))
    for value, offset, solid in zip(values, offsets, filled):
        y = value if log_floor is None else max(value, log_floor)
        ax.plot([x + offset], [y], marker="o", ms=4.0,
                mfc=color if solid else "white", mec=color, mew=0.8,
                ls="none", zorder=5)


def mean_tick(ax: plt.Axes, x: float, values: list[float], *,
              half_width: float = 0.28, log_floor: float | None = None) -> None:
    """Black tick marking the mean of the independent pairs of one column."""
    finite = [v for v in values if v is not None and np.isfinite(v)]
    if not finite:
        return
    y = float(np.mean(finite))
    if log_floor is not None:
        y = max(y, log_floor)
    ax.plot([x - half_width, x + half_width], [y, y], color="black", lw=2.0,
            zorder=4, solid_capstyle="butt")


def _key_handles(spec: list[tuple[str, dict]]) -> list[Line2D]:
    """Legend handles from a ``(label, kwargs)`` spec list."""
    return [Line2D([], [], label=label, **kw) for label, kw in spec]


def assert_key_fits(fig: plt.Figure, legend: Legend, what: str,
                    *, margin_in: float = 0.04) -> None:
    """Abort when a key line runs off its canvas.

    The key lines of this page are figure legends on their own band, so a key
    that grew (panel B gained two mouse-study entries when it was split into
    two facets) silently overflows the canvas instead of colliding with
    something the layout gate would notice.  Measure the real extent, never
    estimate it from the label lengths.
    """
    fig.canvas.draw()
    box = legend.get_window_extent(fig.canvas.get_renderer())
    right_in = box.x1 / fig.dpi
    left_in = box.x0 / fig.dpi
    width_in = box.width / fig.dpi
    print(f"[key] {what}: {width_in:.2f} in wide, spans "
          f"{left_in:.2f}--{right_in:.2f} in of {fig.get_figwidth():.2f} in")
    if left_in < -0.5 / fig.dpi or right_in > fig.get_figwidth() - margin_in:
        raise SystemExit(f"{what}: key line runs off the canvas "
                         f"({left_in:.2f}--{right_in:.2f} in of "
                         f"{fig.get_figwidth():.2f} in)")


def _style_axes(ax: plt.Axes, *, labelsize: float | None = None) -> None:
    ax.tick_params(direction="out", length=3.5, width=1.0,
                   labelsize=FS["tick"] if labelsize is None else labelsize)


# --------------------------------------------------------------------------- #
# data helpers
# --------------------------------------------------------------------------- #
def load_auprc_per_unit() -> pd.DataFrame:
    """Long form of ``figS4_auprc_per_unit.tsv`` (one row per unit/window)."""
    raw = pd.read_csv(TAB / "figS4_auprc_per_unit.tsv", sep="\t")
    parts = []
    for window in WINDOWS:
        column = f"pr_auc_w{window}"
        if column not in raw.columns:
            continue
        part = raw[["species_group", "unit", "tool", column]].copy()
        part.columns = ["species_group", "unit", "tool", "pr_auc"]
        part["window"] = window
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


def _assert_unit_means(df_unit: pd.DataFrame, df_summary: pd.DataFrame,
                       tol: float = 1e-5) -> None:
    """The per-unit mean must reproduce the summary table's ``pr_auc_mean``.

    Tolerance is 1e-5 relative because the frozen tables store floats as
    ``%.6g`` (same convention as the reconciliation in ``41_figS4_tables.py``).
    """
    left = (df_unit.groupby(["species_group", "tool", "window"])["pr_auc"]
            .mean().rename("unit_mean").reset_index())
    right = df_summary[["species_group", "tool", "window", "pr_auc_mean"]]
    merged = left.merge(right, on=["species_group", "tool", "window"],
                        how="inner")
    bad = merged[np.abs(merged["unit_mean"] - merged["pr_auc_mean"])
                 > tol * merged["pr_auc_mean"].abs().clip(lower=1e-12)]
    if not bad.empty:
        raise SystemExit("per-unit PR-AUC does not reproduce the summary table:\n"
                         + bad.head(10).to_string(index=False))


# --------------------------------------------------------------------------- #
# panel A -- counts, purified counts and the purified/WT ratio
# --------------------------------------------------------------------------- #
def render_panel_a(df_counts: pd.DataFrame, *, stem: str = "FigureS4_A",
                   gate: bool = True) -> panelpage.Panel:
    """Counts / purified counts / purified-WT ratio, three species columns.

    Dots are the individual independent pairs and the black tick is their mean
    (the bar version was rejected by the user).  The mouse column gets no
    combined mean: the two cross-study datasets are drawn as two dots joined by
    a thin line, so nothing is ever averaged.

    The sub-rows are separated by dedicated title bands: a band carries the
    title of the sub-row *below* it and has to stay clear of the lowest y tick
    label of the sub-row *above* it, which is the overlap the user saw.
    """
    fig = panelpage.new_panel(PANEL_W_IN, A_H_IN)
    y = 0.0
    y_ratio = y + A_LABEL_IN
    y_band1 = y_ratio + A_RATIO_IN
    y_pur = y_band1 + A_BAND_IN
    y_band2 = y_pur + A_SUB_IN
    y_wt = y_band2 + A_BAND_IN
    y_wt_title = y_wt + A_SUB_IN            # "WT counts (log)" line
    # species names on their own line, with a visible gap above the sub-row
    # title (0.22 in read as "text on text" at the bold 14 pt size)
    y_species = y_wt_title + 0.24
    for ci, (key, species, pairs) in enumerate(COLUMNS):
        color = GROUP_COLOR[key]
        mouse = key == "Mouse"
        x_in = col_x_in(ci)
        ax_wt = fig.add_axes(_rect(A_H_IN, x_in=x_in, y_in=y_wt,
                                   w_in=COL_W_IN, h_in=A_SUB_IN))
        ax_pu = fig.add_axes(_rect(A_H_IN, x_in=x_in, y_in=y_pur,
                                   w_in=COL_W_IN, h_in=A_SUB_IN),
                             sharex=ax_wt)
        ax_ra = fig.add_axes(_rect(A_H_IN, x_in=x_in, y_in=y_ratio,
                                   w_in=COL_W_IN, h_in=A_RATIO_IN),
                             sharex=ax_wt)
        sub = df_counts[df_counts["species"] == species]
        for ti, tool in enumerate(TOOL_ORDER):
            rows = pair_rows(sub, tool, pairs)
            wt = [float(r["n_wt_common"]) for r in rows]
            pu = [float(r["n_purified"]) for r in rows]
            ratio = [p / w if w > 0 else np.nan for p, w in zip(pu, wt)]
            filled = [r["pair_id"] != "Mouse_studyB" for r in rows]
            for ax, vals in ((ax_wt, wt), (ax_pu, pu), (ax_ra, ratio)):
                floor = 0.7 if ax is not ax_ra else None
                if mouse:  # two studies: dots + connector, never a mean
                    if len(vals) == 2:
                        ax.plot([ti - 0.20, ti + 0.20],
                                [max(vals[0], floor) if floor else vals[0],
                                 max(vals[1], floor) if floor else vals[1]],
                                color="0.70", lw=0.8, zorder=2)
                else:
                    mean_tick(ax, ti, vals, log_floor=floor)
            pair_markers(ax_wt, ti, wt, filled, color, log_floor=0.7)
            pair_markers(ax_pu, ti, pu, filled, color, log_floor=0.7)
            pair_markers(ax_ra, ti, ratio, filled, color)
        for ax in (ax_wt, ax_pu):
            ax.set_yscale("log")
            ax.set_ylim(0.6, 4e5)
            # three power-of-ten labels: the 0.50 in sub-rows keep the labels
            # ~12 pt apart (the 0.46 in sub-rows of the four-row version, i.e.
            # before the PR-AUC row was dropped, only fit two)
            ax.set_yticks([1, 100, 10000])
            ax.set_yticklabels([r"$10^{0}$", r"$10^{2}$", r"$10^{4}$"],
                               fontsize=FS["tick_small"])
            ax.yaxis.set_minor_locator(NullLocator())
            ax.tick_params(axis="x", which="both", labelbottom=False)
        ax_ra.set_ylim(-0.03, 1.03)
        ax_ra.set_yticks([0.0, 0.5, 1.0])
        ax_ra.set_yticklabels(["0", "0.5", "1.0"], fontsize=FS["tick_small"])
        tool_ticks(ax_ra)
        for ax in (ax_wt, ax_pu, ax_ra):
            _style_axes(ax, labelsize=FS["tick_small"])
        # titles live in their own bands, never on an axis line: the 0.03 in
        # offset keeps the ink (descenders included) off the spine below and
        # the indent keeps it clear of the top of the left spine
        t_x = x_in + A_TITLE_INDENT_IN
        _title(fig, A_H_IN, t_x, y_wt_title + 0.03, "WT counts (log)")
        _title(fig, A_H_IN, t_x, y_band2 + 0.03, "Purified counts (log)")
        _title(fig, A_H_IN, t_x, y_band1 + 0.03, "Purified/WT ratio")
        _title(fig, A_H_IN, col_center_in(ci), y_species,
               COLUMN_TITLE[key], ha="center", fontsize=FS["column"],
               weight="bold")
    # the species / mouse-study / mean key sits on its own line at the top of
    # this canvas: the three sub-rows run full to 1.0, so none has room
    # centred key: at "upper left" it read as belonging to the first column only
    legend = fig.legend(handles=_key_handles(A_KEY), loc="upper center",
                        bbox_to_anchor=(0.5, 1 - 0.05 / A_H_IN),
                        ncol=len(A_KEY), frameon=False, fontsize=FS["legend"],
                        handlelength=1.5, columnspacing=1.6, labelspacing=0.2,
                        borderaxespad=0.0)
    assert_key_fits(fig, legend, "panel A")
    panel_letter(fig, "A", y_in=A_H_IN - 0.04)
    return panelpage.save_panel(fig, PANEL_DIR / f"{stem}.pdf", gate=gate)


# --------------------------------------------------------------------------- #
# panels B/C -- PPV and single-nucleotide accuracy vs. the matching window
# --------------------------------------------------------------------------- #
def render_panels_bc(df_curve: pd.DataFrame, *, stem: str = "FigureS4_BC",
                     gate: bool = True) -> panelpage.Panel:
    """Two stacked facets against the matching window -- **panels B and C**.

    The submitted panel B overlaid two different quantities on one 0--1 axis
    (13 tools x three line families), which the user judged too crowded
 ("B panel is too dense split into two facets "). Every species column now holds two
    facets that share the x axis, and *each facet carries a panel letter of its
 own* (2026-09-21, " each cell keeps its own letter "):

    * upper facet = **B** = ``hit_rate``, the PPV vs. GLORI of the window;
    * lower facet = **C** = ``localization_accuracy``, the fraction of the
      window-matched calls that fall on the exact reference nucleotide -- the
      relationship between extension distance and single-nucleotide accuracy
      the reviewers asked for (R1-3 / E2 / R3-3).  Read top-to-bottom the
      column now shows directly that the apparent PPV rises with the window
      while the exact fraction falls.

    Both facets carry the same two layers per series (the per-replicate policy:
    the thin lines are the individual sequencing units and the thick line is
    their group mean), and every column keeps its own y tick labels.  The two
    letters live in the page margin, next to their own facet; the canvas stays
    a single one so the facets keep their near-square 2.16 x 1.83 in shape.

    The two mouse studies are never merged: each cross-study dataset
    contributes exactly one sequencing unit, so the mouse column draws one
    thick line per study and facet (study A solid, study B dashed) and no thin
    layer.
    """
    canvas_h = B_H_IN
    fig = panelpage.new_panel(PANEL_W_IN, canvas_h)
    y_lo = B_BOTTOM_IN
    y_up = B_BOTTOM_IN + B_FACET_IN + B_GAP_IN

    #: this canvas into the band between the two facets, where it explains the
    #: upper facet (B) and the lower one (C) at once and hides no data
    #: both keys of this canvas print in the band between the two facets, one
    #: above the other: the line-type key first (it explains the thin/thick and
    #: solid/dashed strokes of both facets), the 13 tool colours below it
    band = B_BOTTOM_IN + B_FACET_IN
    line_legend = fig.legend(handles=_key_handles(LINE_KEY_B), loc="center",
                             bbox_to_anchor=(0.5, (band + 0.51) / canvas_h),
                             ncol=len(LINE_KEY_B), frameon=False,
                             fontsize=FS["legend"], handlelength=1.6,
                             columnspacing=1.4, labelspacing=0.2,
                             borderaxespad=0.0)
    tool_legend = fig.legend(handles=_key_handles(TOOL_KEY), loc="center",
                             bbox_to_anchor=(0.5, (band + 0.235) / canvas_h),
                             ncol=7, frameon=False, fontsize=FS["legend"],
                             handlelength=1.2, columnspacing=0.9,
                             labelspacing=0.30, borderaxespad=0.0)
    #: (table column, facet index); facet 0 is the lower one, 1 the upper one
    facet_spec = (("localization_accuracy", 0), ("hit_rate", 1))
    for ci, (key, _species, _pairs) in enumerate(COLUMNS):
        facets = [fig.add_axes(_rect(canvas_h, x_in=col_x_in(ci), y_in=y_lo,
                                     w_in=COL_W_IN, h_in=B_FACET_IN)),
                  fig.add_axes(_rect(canvas_h, x_in=col_x_in(ci), y_in=y_up,
                                     w_in=COL_W_IN, h_in=B_FACET_IN))]
        facets[0].sharex(facets[1])          # one window axis for both facets
        sub = df_curve[df_curve["sample"].isin(CURVE_SAMPLES[key])]
        for tool in TOOL_ORDER:
            t = sub[sub["tool"] == tool]
            if t.empty:
                continue
            colour = TOOL_COLOR[tool]
            for column, fi in facet_spec:
                ax = facets[fi]
                if key == "Mouse":  # one line per study, never merged
                    for _study, sample, ls, _filled in MOUSE_STUDIES:
                        st = t[t["sample"] == sample].sort_values("window")
                        if st.empty:
                            continue
                        ax.plot(st["window"], st[column], color=colour, lw=2.0,
                                ls=ls, zorder=4)
                else:
                    for _unit, u in t.groupby("sample"):
                        u = u.sort_values("window")
                        ax.plot(u["window"], u[column], color=colour, lw=0.7,
                                alpha=0.25, zorder=2)
                    mean = t.groupby("window")[column].mean()
                    ax.plot(mean.index, mean.values, color=colour, lw=2.0,
                            zorder=4)
        for ax in facets:
            ax.set_xlim(-1, 51)
            # 3 % of headroom: a marker at fraction 1.0 was clipped by the frame
            ax.set_ylim(0, 1.03)
            ax.set_yticks([0.0, 0.5, 1.0])
            ax.set_yticklabels(["0", "0.5", "1.0"], fontsize=FS["tick"])
            _style_axes(ax)
        # the facets share the x axis: the ticks and the axis name belong to the
        # lower facet, the upper one keeps only its frame
        facets[0].set_xticks([0, 20, 50])
        facets[0].set_xlabel("Matching window (bp)", fontsize=FS["axis"],
                             labelpad=3)
        facets[1].tick_params(axis="x", which="both", bottom=False,
                              labelbottom=False)
        if ci == 0:
            facets[1].set_ylabel("PPV vs. GLORI", fontsize=FS["axis"],
                                 labelpad=3)
            facets[0].set_ylabel(B_EXACT_LABEL, fontsize=FS["axis"], labelpad=3)
        _title(fig, canvas_h, col_center_in(ci), y_up + B_FACET_IN + 0.02,
               COLUMN_TITLE[key], ha="center", fontsize=FS["column"],
               weight="bold")
    #: both keys are measured, so neither can grow past the band unnoticed
    assert_key_fits(fig, line_legend, "line-type key between panels B and C")
    assert_key_fits(fig, tool_legend, "tool key between panels B and C")
    # two letters for two panels: each one sits in the page margin next to the
    # top edge of its own facet (the block letter convention of rows A and D is
    # kept, but the window sweep is deliberately split into B and C)
    panel_letter(fig, "B", y_in=y_up + B_FACET_IN - 0.04)
    panel_letter(fig, "C", y_in=y_lo + B_FACET_IN - 0.04)
    return panelpage.save_panel(fig, PANEL_DIR / f"{stem}.pdf", gate=gate)


# --------------------------------------------------------------------------- #
# panel D -- PPV of the WT-common calls vs their purified subset
# --------------------------------------------------------------------------- #
def _dumbbells(df_prec: pd.DataFrame, species: str,
               pairs: list[str]) -> list[tuple[str, float, float]]:
    sub = df_prec[df_prec["species"] == species]
    out: list[tuple[str, float, float]] = []
    for tool in TOOL_ORDER:
        for row in pair_rows(sub, tool, pairs):
            out.append((tool, float(row["precision_wt_common_w2"]),
                        float(row["precision_purified_w2"])))
    out.sort(key=lambda item: item[2] - item[1])
    return out


def render_panel_d(df_prec: pd.DataFrame, *, stem: str = "FigureS4_D",
                   gate: bool = True, kind: str = "dumbbell") -> panelpage.Panel:
    """PPV (2 bp) of the WT-common calls vs their purified subset -- **panel D**.

    Dumbbells connect the same tool in the same pair before and after
    purification (colour = tool, ordered by the shift so the connectors do not
    cross); the black point with the error bar is the mean +- SD over all
    tool-pair combinations.  Descriptive only -- conditioning on absence from a
    deficient call set is circular for accuracy assessment (R3-7).

    This panel was lettered C until 2026-09-21; it became D when the two facets
    of the window sweep were given letters of their own (B and C).
    """
    fig = panelpage.new_panel(PANEL_W_IN, D_H_IN)
    for ci, (key, species, pairs) in enumerate(COLUMNS):
        ax = fig.add_axes(_rect(D_H_IN, x_in=col_x_in(ci), y_in=D_BOTTOM_IN,
                                w_in=COL_W_IN, h_in=D_AXES_IN))
        items = _dumbbells(df_prec, species, pairs)
        if kind == "box":
            for gi, tag in enumerate(("wt", "pu")):
                vals = [item[1 + gi] for item in items]
                ax.boxplot([vals], positions=[gi], widths=0.55,
                           showfliers=False, patch_artist=True, zorder=2,
                           medianprops=dict(color="black", lw=1.8),
                           boxprops=dict(facecolor="white", edgecolor="0.20",
                                         lw=1.1),
                           whiskerprops=dict(color="0.20", lw=1.1),
                           capprops=dict(color="0.20", lw=1.1))
                rng = np.random.default_rng(20260921)
                for tool, y0, y1 in items:
                    ax.plot([gi + rng.uniform(-0.22, 0.22)], [y0 if gi == 0 else y1],
                            "o", ms=2.0, mfc=TOOL_COLOR[tool], mec="none",
                            alpha=0.45, ls="none", zorder=3)
                sd = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
                ax.errorbar([gi + 0.22], [float(np.mean(vals))], yerr=[sd],
                            fmt="D", ms=5, color="black", ecolor="black",
                            elinewidth=1.4, capsize=0, zorder=5)
            ax.set_xlim(-0.6, 1.6)
        else:
            for tool, y0, y1 in items:
                ax.plot([0, 1], [y0, y1], color=TOOL_COLOR[tool], lw=0.6,
                        alpha=0.35, zorder=2)
                ax.plot([0, 1], [y0, y1], marker="o", ms=2.4, ls="none",
                        mfc=TOOL_COLOR[tool], mec=TOOL_COLOR[tool], mew=0.0,
                        alpha=0.85, zorder=3)
            for x, column in ((0, 1), (1, 2)):
                vals = [item[column] for item in items]
                if not vals:
                    continue
                sd = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
                ax.errorbar([x], [float(np.mean(vals))], yerr=[sd], fmt="o",
                            ms=6, color="black", ecolor="black", elinewidth=1.5,
                            capsize=3.5, capthick=1.5, zorder=5)
            ax.set_xlim(-0.35, 1.35)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["All WT-common", "Purified subset"],
                           fontsize=FS["tick"])
        ax.set_ylim(0, 1.03)          # headroom: PPV = 1.0 markers were clipped
        ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1.0"],
                           fontsize=FS["tick"])
        _style_axes(ax)
        if ci == 0:
            ax.set_ylabel("PPV vs. GLORI (2 bp)", fontsize=FS["axis"],
                          labelpad=3)
        _title(fig, D_H_IN, col_center_in(ci), D_BOTTOM_IN + D_AXES_IN + 0.02,
               COLUMN_TITLE[key], ha="center", fontsize=FS["column"],
               weight="bold")
    # panel D carries its own key line (it used to sit on the B/C canvas).

    # full-width canvas, so it floated in the gutter between the Arabidopsis and
    # Mouse panels.  It is now aligned with the first D panel's left edge, in the
    # same title band as the species names (which sit lower in the band).
    legend = fig.legend(handles=_key_handles(D_KEY), loc="upper left",
                        bbox_to_anchor=(LEFT_IN / PANEL_W_IN,
                                        1 - 0.04 / D_H_IN),
                        ncol=len(D_KEY),
                        frameon=False, fontsize=FS["legend"], handlelength=1.6,
                        columnspacing=1.4, labelspacing=0.2, borderaxespad=0.0)
    assert_key_fits(fig, legend, "panel D")
    panel_letter(fig, "D", y_in=D_H_IN - 0.04)
    return panelpage.save_panel(fig, PANEL_DIR / f"{stem}.pdf", gate=gate)


# --------------------------------------------------------------------------- #
# the tool colour key: its own centred stripe, so it cannot eat a panel slot
# --------------------------------------------------------------------------- #
def render_tool_stripe(*, stem: str = "FigureS4_legend_tools",
                       gate: bool = True) -> panelpage.Panel:
    """Superseded 2026-09-27: kept for reference, no longer on the page.

    Until 2026-09-27 the 13 tool colours were a centred stripe between panel A
    and the B/C canvas.  They are now drawn inside the BC canvas, between its
    two facets, so this page piece is unused.

    It is a canvas of its own (full page width, :data:`TOOL_STRIPE_IN` deep):
    the page budget reserves a band for it, so no panel loses drawing area and
    the key sits at the middle of the page, where the reader looks for it after
    panel A.  Two rows of seven entries at the key font size; nothing is drawn
    behind it, so it can never hide data.
    """
    fig = panelpage.new_panel(PANEL_W_IN, TOOL_STRIPE_IN)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_axis_off()
    ax.legend(handles=_key_handles(TOOL_KEY), loc="center", ncol=7,
              frameon=False, fontsize=FS["legend"], handlelength=1.6,
              columnspacing=1.6, labelspacing=0.30, borderaxespad=0.0)
    return panelpage.save_panel(fig, PANEL_DIR / f"{stem}.pdf",
                                ignore_axes=[ax], gate=gate)


# --------------------------------------------------------------------------- #
# page assembly
# --------------------------------------------------------------------------- #
def page_placements(panels: dict[str, panelpage.Panel],
                    ) -> list[tuple[panelpage.Panel, float, float]]:
    """Stack the pieces top-down; returns ``[(panel, x_in, y_in)]``.

    ``y_in`` is the lower-left corner on the page.  The pieces are laid out
    from a budget in inches, and the leftover is asserted to be non-negative.
    The tool key is no longer a page piece: since 2026-09-27 it is drawn inside
    the BC canvas, between its two facets (B above, C below).
    """
    order = [("A", GAP_IN), ("BC", GAP_IN), ("D", 0.0)]
    cursor = PAGE_IN[1] - TOP_IN
    out: list[tuple[panelpage.Panel, float, float]] = []
    for key, gap in order:
        panel = panels[key]
        bottom = cursor - panel.height_in
        if bottom < -1e-9:
            raise SystemExit(f"{key}: page budget overflows by {-bottom:.3f} in")
        out.append((panel, 0.0, bottom))
        cursor = bottom - gap
    if cursor < BOTTOM_IN - 1e-9:
        raise SystemExit(f"bottom margin {cursor:.3f} in < {BOTTOM_IN:.2f} in")
    return out


# --------------------------------------------------------------------------- #
# numeric anchors
# --------------------------------------------------------------------------- #
def write_anchors(df_curve: pd.DataFrame, df_prec: pd.DataFrame,
                  df_auprc: pd.DataFrame) -> pd.DataFrame:
    """Log the numbers the legend text quotes (checked against the tables).

    The ``pr_auc_*`` rows are the one block that no longer feeds a panel: they
    keep the archived PR-AUC tables of the dropped PR-AUC row under the same
    1e-6 check as everything else, so the numbers cannot drift unnoticed if that
    row is ever brought back.
    """
    rows: list[dict[str, object]] = []
    for key, _species, _pairs in COLUMNS:
        sub = df_curve[df_curve["sample"].isin(CURVE_SAMPLES[key])]
        for window in (0, 2, 50):
            rows.append({
                "anchor": f"ppv_w{window}_{key}",
                "value": sub.loc[sub["window"] == window, "hit_rate"].mean(),
                "definition": "mean PPV vs. GLORI over tools and independent "
                              "units (frozen localisation curve)"})
    for column, label in (("precision_wt_common_w2", "ppv_wt_common"),
                          ("precision_purified_w2", "ppv_purified")):
        rows.append({"anchor": label, "value": df_prec[column].mean(),
                     "definition": "mean over the 104 tool-pair combinations "
                                   "(2 bp window)"})
    delta = df_prec["precision_purified_w2"] - df_prec["precision_wt_common_w2"]
    rows.append({"anchor": "n_tool_pair_up",
                 "value": int((delta > 0).sum()),
                 "definition": "tool-pair combinations whose PPV rises after "
                               "purification"})
    for group in ("Arabidopsis", "Mouse_studyA", "Mouse_studyB", "HeLa"):
        sub = df_auprc[df_auprc["species_group"] == group]
        for window in (0, 50):
            rows.append({
                "anchor": f"pr_auc_w{window}_{group}",
                "value": sub.loc[sub["window"] == window, "pr_auc_mean"].mean(),
                "definition": "mean PR-AUC over tools (universe ranking)"})
    out = pd.DataFrame(rows)
    LOG.mkdir(parents=True, exist_ok=True)
    out.to_csv(LOG / "42_figS4_anchors.tsv", sep="\t", index=False)
    return out


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main() -> None:
    global TICK_ROTATION
    t0 = time.time()
    parser = argparse.ArgumentParser(description="Figure S4 (A-D) renderer.")
    parser.add_argument("--rotation", type=int, default=TICK_ROTATION,
                        choices=(45, 90),
                        help="rotation of the panel-A tool labels (default 45)")
    parser.add_argument("--panel-d", "--panel-c", dest="panel_d",
                        default="dumbbell", choices=("dumbbell", "box"),
                        help="rendering of panel D (the old --panel-c alias "
                             "is kept for the pre-2026-09-21 lettering)")
    parser.add_argument("--no-overlap-check", action="store_true",
                        help="skip the per-panel layout gate (tuning only)")
    parser.add_argument("--out", default="FigureS4_rev",
                        help="output stem inside figures/ (default FigureS4_rev)")
    args = parser.parse_args()
    TICK_ROTATION = args.rotation
    gate = not args.no_overlap_check
    apply()
    plt.rcParams.update({
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "axes.linewidth": 1.2,
    })
    FIG.mkdir(parents=True, exist_ok=True)

    df_counts = pd.read_csv(TAB / "figS4_purified_sites.tsv", sep="\t")
    df_prec = pd.read_csv(TAB / "figS4_precision_on_purified.tsv", sep="\t")
    df_auprc = pd.read_csv(TAB / "figS4_auprc_summary.tsv", sep="\t")
    df_curve = pd.read_csv(TABLE_DIR / "m6a_localization_curve.tsv", sep="\t")
    df_curve = df_curve[(df_curve["platform"] == "RNA002")
                        & (df_curve["tool"].isin(TOOL_ORDER))]
    # the two frozen PR-AUC tables of the dropped PR-AUC row are still read
    # once: the per-unit reshape has to reproduce the summary table (no panel is
    # drawn from them any more -- they fed the fourth row that was dropped on
    # 2026-09-21, which is *not* the panel D of the B/C re-lettering; the
    # archived numbers stay consistent either way)
    df_unit = load_auprc_per_unit()
    _assert_unit_means(df_unit, df_auprc)

    anchors = write_anchors(df_curve, df_prec, df_auprc)

    longest = pagelayout.rotated_label_in([DISPLAY.get(t, t) for t in TOOL_ORDER],
                                          FS["tick_small"])
    tilt = np.deg2rad(TICK_ROTATION)
    pad_in = (5.0 if TICK_ROTATION == 90 else 8.0) / 72.0
    layout_h = 1.117 * FS["tick_small"] / 72.0     # ascent + descent of Arial
    need = pad_in + longest * np.sin(tilt) + layout_h * np.cos(tilt)
    print(f"[layout] panel {PANEL_W_IN:.2f} in wide, column {COL_W_IN:.2f} in "
          f"(B/C facets {COL_W_IN:.2f} x {B_FACET_IN:.2f} in each, two per "
          f"column; D {COL_W_IN:.2f} x {D_AXES_IN:.2f} in), slot "
          f"{COL_W_IN / len(TOOL_ORDER) * 72:.2f} pt; panel-A names at "
          f"{TICK_ROTATION} deg need {need:.3f} in of {A_LABEL_IN:.2f} in")
    if need > A_LABEL_IN:
        raise SystemExit(f"panel-A label band too small: {need:.3f} > "
                         f"{A_LABEL_IN:.2f} in")

    panels = {
        "A": render_panel_a(df_counts, gate=gate),

        #: the page has three pieces, not four
        "BC": render_panels_bc(df_curve, gate=gate),
        "D": render_panel_d(df_prec, gate=gate, kind=args.panel_d),
    }
    if args.panel_d == "box":
        panels["D"] = render_panel_d(df_prec, stem="FigureS4_Dbox",
                                     gate=gate, kind="box")

    out_pdf = FIG / f"{args.out}.pdf"
    page = panelpage.compose_page(PAGE_IN, page_placements(panels), out_pdf)
    panelpage.pdf_to_png(page, FIG / f"{args.out}.png")

    print(f"saved -> {out_pdf} + {FIG / (args.out + '.png')} "
          f"({time.time() - t0:.1f} s)")
    print(f"anchors -> {LOG / '42_figS4_anchors.tsv'} "
          f"({datetime.now():%Y-%m-%d %H:%M})")
    with pd.option_context("display.width", 120, "display.max_columns", 8):
        print(anchors.to_string(index=False))


if __name__ == "__main__":
    main()
