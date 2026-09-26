#!/usr/bin/env python
"""S8 style tokens -- page geometry, type scale, palette.

Since 2026-09-25 (user) the page is an **A4 portrait sheet with three rows**:
row 1 carries the three bar panels A | B | C side by side, row 2 the full-width
window sweep D (its 13-entry key on the right), row 3 the merged effect-size
panel E (Pearson r and slope sharing one label column).

Layout arithmetic (inches, 1:1 print size)::

    page   595.44 x 842.4 pt  =  8.270 x 11.700
    margins 0.10 all round, gutters 0.14 (vertical and between the three columns)
    row heights 4.350 (A|B|C) + 3.900 (D) + 2.970 (E) + 2*0.14 = 11.500
    column width (8.27 - 0.20 - 2*0.14) / 3 = 2.5967

Every panel is drawn on its own print-size canvas and the page is composed 1:1
with pypdf (no scaling), so the type sizes below are the printed sizes.

House rules honoured here: Arial only, no grid, no in-panel annotation text
(only bold panel letters and >= 7 pt entries), vector PDF + 300 dpi PNG at the
exact page size (never ``bbox_inches="tight"``).
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

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# --------------------------------------------------------------------------- #
# paths
# --------------------------------------------------------------------------- #
PROJECT = Path(str(_RB))

#: the shared versioned helpers (layout gate, panel/page composition)
sys.path.insert(0, str((_RB / "src/sites_v2")))
OUT = (_RB / "figures/figureS9")
TABLES = (_RB / "figures/figureS9/tables")
FIGS = (_RB / "figures/figureS9/figures")
PANELS = (_RB / "figures/figureS9/figures/panels")
LOGS = (_RB / "figures/figureS9/logs")

#: frozen evidence tables owned by the sibling Figure-8 workflow (read-only)
SRC_TABLES = (_RB / "figures/figureS9/tables/source")

# --------------------------------------------------------------------------- #
# page geometry (inches)
# --------------------------------------------------------------------------- #
PAGE_IN = (842.4 / 72.0, 595.44 / 72.0)         # A4 landscape (2026-09-25 layout)

M_LEFT = M_RIGHT = 0.100
M_TOP = M_BOTTOM = 0.100
GUT_H = 0.140                                   # between the two columns
GUT_V = 0.140                                   # between the rows

N_COL = 3                                       # row 1 carries A | B | C
W_FULL = PAGE_IN[0] - M_LEFT - M_RIGHT          # 11.500 (landscape)
COL_W = (W_FULL - (N_COL - 1) * GUT_H) / N_COL  # 3.740
ROW1_W = (COL_W, COL_W, COL_W)


#: (D, E) on the second row and the effect sizes are two facets (F, G) on the
#: third, each facet with its own panel letter and the key on the right.  Row 1
#: shrank to the 23-entry count matrix it has to hold (23 rows at an 8.4 pt
#: pitch) and the freed height went into the two lower rows, which is what makes
#: the four facets read as panels instead of letterbox strips.
#: 3.600 + 2.100 + 2.090 + 2 * 0.140 = 8.070 = page height minus the margins.
ROW_H = {"r1": 3.600, "r2": 2.100, "r3": 2.090}

_X_L = M_LEFT
_X_M = M_LEFT + ROW1_W[0] + GUT_H
_X_R = _X_M + ROW1_W[1] + GUT_H
_Y_R3 = M_BOTTOM
_Y_R2 = _Y_R3 + ROW_H["r3"] + GUT_V
_Y_R1 = _Y_R2 + ROW_H["r2"] + GUT_V


#: facet gets its own canvas and its own letter -- row 2 carries D | E (the
#: window sweep: PPV, exact fraction) and row 3 carries F | G (the effect sizes:
#: Pearson r, slope).  The right-hand canvas of each row also carries the key.
ROW2_W = {"D": 3.95, "E": W_FULL - 3.95}
ROW3_W = {"F": 4.90, "G": W_FULL - 4.90}

#: panel -> (width, height) of its own canvas, in inches (print size)
PIECE = {
    "A": (ROW1_W[0], ROW_H["r1"]), "B": (ROW1_W[1], ROW_H["r1"]),
    "C": (ROW1_W[2], ROW_H["r1"]),
    "D": (ROW2_W["D"], ROW_H["r2"]), "E": (ROW2_W["E"], ROW_H["r2"]),
    "F": (ROW3_W["F"], ROW_H["r3"]), "G": (ROW3_W["G"], ROW_H["r3"]),
}

#: panel -> lower-left corner on the composed page, in inches
PLACE = {
    "A": (_X_L, _Y_R1), "B": (_X_M, _Y_R1), "C": (_X_R, _Y_R1),
    "D": (M_LEFT, _Y_R2), "E": (M_LEFT + ROW2_W["D"], _Y_R2),
    "F": (M_LEFT, _Y_R3), "G": (M_LEFT + ROW3_W["F"], _Y_R3),
}

# --------------------------------------------------------------------------- #
# type scale (points, printed 1:1); hard floor 7 pt (the layout gate enforces it)
# --------------------------------------------------------------------------- #
FS = {
    "tick": 10.5,       # ticks of the wide axes (window sweep, effect size)
    "tick_sm": 9.5,     # ticks of the narrow blocks (counts, false-positive rate)
    "tick_xs": 9.0,     # ticks of the specificity axis
    "axis": 12.5,       # axis titles
    "axis_sm": 11.0,    # axis titles of the narrow blocks
    "row": 8.0,         # row labels (landscape rows are shallow: the 23-row
    "row_a": 8.0,       # matrix of A would else collide; floor is 7 pt)
    "legend": 10.5,
    "letter": 24.0,
    "title": 10.5,
    "stripe": 9.0,
}
MIN_PT = 7.0

matplotlib.rcParams.update({
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
})

WT = "#4B81B8"
IVT = "#E8A76B"
CURLAKE = "#8C8C8C"
FAM = {
    "m6A_DRACH": "#E83C1E",
    "m6A": "#F5B264",
    "inosine_m6A": "#F18F01",
    "pseU": "#A23B72",
    "m5C": "#2E86AB",
    "tool": "#7F7F7F",
}
FAM_LABEL = {
    "m6A_DRACH": "m6A DRACH",
    "m6A": "m6A (non-DRACH)",
    "inosine_m6A": "inosine + m6A",
    "pseU": "pseU",
    "m5C": "m5C",
    "tool": "ORCA",
}
MEAN = "#222222"
#: the benchmark dot colour of the main figures (Figure 8 / S6 convention)
DOT_COLOR = "#1E888B"
GUIDE = "0.72"
AXIS_LW = 0.9

def apply() -> None:
    """Arial everywhere, no grid, vector-friendly fonts."""
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": FS["tick"],
        "axes.linewidth": AXIS_LW,
        "axes.grid": False,
        "xtick.labelsize": FS["tick"],
        "ytick.labelsize": FS["tick"],
        "legend.fontsize": FS["legend"],
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def bare(ax: plt.Axes) -> None:
    """House axis furniture: no top/right spine, ticks out, no minor ticks."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(AXIS_LW)
    ax.tick_params(direction="out", length=2.4, width=0.8, pad=3.6)
    ax.set_axisbelow(True)
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())


def letter(fig: plt.Figure, key: str, x: float = 0.006) -> None:
    """Bold panel letter in the top-left corner of a panel canvas.

    ``x`` (figure fraction) lets a canvas that carries two facets label the
    second one at its own left edge (2026-09-25: the window sweep and the
    effect-size row are two lettered facets each).
    """
    fig.text(x, 0.994, key, fontsize=FS["letter"], fontweight="bold",
             ha="left", va="top")


def dot_column_axes(ax: plt.Axes, n_rows: int) -> None:
    """A label-free dot column: no spines, no ticks, rows top-down."""
    bare(ax)
    ax.set_xlim(-0.55, 0.55)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xticks([])
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)


def log_ticks(ax: plt.Axes, ticks: list[float], labels: list[str] | None = None,
              axis: str = "x") -> None:
    """Explicit log tick set, no minor ticks, no 10^0-style auto labels."""
    loc = matplotlib.ticker.FixedLocator(ticks)
    target = ax.xaxis if axis == "x" else ax.yaxis
    target.set_major_locator(loc)
    if labels is not None:
        target.set_major_formatter(matplotlib.ticker.FixedFormatter(labels))
    target.set_minor_locator(matplotlib.ticker.NullLocator())


DECADES = {1e-4: "$10^{-4}$", 1e-3: "$10^{-3}$", 1e-2: "$10^{-2}$",
           1e-1: "0.1", 1: "1", 10: "10", 100: "100", 1000: "1000",
           10000: "10000", 100000: "100000"}


def tick_labels(ticks: list[float]) -> list[str]:
    return [DECADES.get(t, f"{t:g}") for t in ticks]
