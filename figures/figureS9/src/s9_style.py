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
sys.path.insert(0, str((_RB / "src/harmonisation")))
OUT = (_RB / "figures/figureS9")
TABLES = (_RB / "figures/figureS9/tables")
FIGS = (_RB / "figures/figureS9/figures")
PANELS = (_RB / "figures/figureS9/figures/panels")
LOGS = (_RB / "figures/figureS9/logs")

#: frozen evidence tables owned by the companion Figure-8 workflow (read-only)
SRC_TABLES = (_RB / "figures/figureS9/tables/source")

# --------------------------------------------------------------------------- #
# page geometry (inches)
# --------------------------------------------------------------------------- #
PAGE_IN = (842.4 / 72.0, 595.44 / 72.0)         # A4 landscape (restored layout)

M_LEFT = M_RIGHT = 0.100
M_TOP = M_BOTTOM = 0.100
GUT_H = 0.060                                   # between the three columns
GUT_V = 0.060                                   # between the three rows

N_COL = 3
W_FULL = PAGE_IN[0] - M_LEFT - M_RIGHT          # 11.500 (landscape)
COL_W = (W_FULL - (N_COL - 1) * GUT_H) / N_COL  # 3.740
ROW1_W = (COL_W, COL_W, COL_W)


#: three rows** (A | B | C, then D | E | H, then F | G | I).  The first row is
#: taller (3.600 in) because the 23-entry count matrix of A needs the height; the
#: two lower rows are 2.095 in.  Every panel draws its data in a **square box** of
#: the same side (BOX_SIDE_IN) on the same bottom offset, so the eight small
#: panels line up row by row on one horizontal plane; A keeps a rectangular box,
#: because a square of that side cannot hold its 23 rows.
#: 3.600 + 2.095 + 2.095 + 2 * 0.140 = 8.070 = page height minus the margins.
ROW_H = {"r1": 3.480, "r2": 2.235, "r3": 2.235}

_X_L = M_LEFT
_X_M = M_LEFT + COL_W + GUT_H
_X_R = _X_M + COL_W + GUT_H
_Y_R3 = M_BOTTOM
_Y_R2 = _Y_R3 + ROW_H["r3"] + GUT_V
_Y_R1 = _Y_R2 + ROW_H["r2"] + GUT_V

#: panel -> (width, height) of its own canvas, in inches (print size)
PIECE = {
    "A": (COL_W, ROW_H["r1"]), "B": (COL_W, ROW_H["r1"]),
    "C": (COL_W, ROW_H["r1"]),
    "D": (COL_W, ROW_H["r2"]), "E": (COL_W, ROW_H["r2"]),
    "H": (COL_W, ROW_H["r2"]),
    "F": (COL_W, ROW_H["r3"]), "G": (COL_W, ROW_H["r3"]),
    "I": (COL_W, ROW_H["r3"]),
}

#: panel -> lower-left corner on the composed page, in inches
PLACE = {
    "A": (_X_L, _Y_R1), "B": (_X_M, _Y_R1), "C": (_X_R, _Y_R1),
    "D": (_X_L, _Y_R2), "E": (_X_M, _Y_R2), "H": (_X_R, _Y_R2),
    "F": (_X_L, _Y_R3), "G": (_X_M, _Y_R3), "I": (_X_R, _Y_R3),
}

# --------------------------------------------------------------------------- #
# type scale (points, printed 1:1); hard floor 7 pt (the layout gate enforces it)
# --------------------------------------------------------------------------- #
FS = {
    "tick": 9.5,
    "tick_sm": 9.5,
    "tick_xs": 9.0,
    "axis": 10.5,
    "axis_sm": 10.5,

    #: longest model name sets -- so A prints its 23 names at the 7 pt floor
    "row": 8.0,
    "level": 8.0,      
    "row_a": 8.0,
    "legend": 9.0,
    "letter": 24.0,
    "title": 10.5,
    "stripe": 9.0,
    #: 2026-09-27: the thirteen-entry sweep key has its own cell (J) and the
    #: species key of the stability cells its own (K)
    "key": 8.0,
    "key_sm": 8.0,
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
    """House axis furniture: a closed box frame, ticks out, no minor ticks.

 2026-09-27 (user): every panel is a **boxed square** (""),
    so the top and right spines stay visible -- the earlier landscape layout kept
    them off.
    """
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(AXIS_LW)
    ax.tick_params(direction="out", length=2.4, width=0.8, pad=3.6)
    ax.set_axisbelow(True)
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())



#: all sit on this bottom offset, so a row of panels reads as one horizontal line.
#: The six cells of rows 2 and 3 draw the same square and leave the same strip
#: under it, because **a key prints under the panel it belongs to** (user rule):
#: H and I carry the species key, G the effect key, each one line under its own
#: box, inside its own canvas -- never above the box and never in a neighbour's
#: free space.  The strip has room for ticks + rotated-nothing axis title + one
#: key line: 0.580 (bottom offset) - 0.37 (furniture) - 0.17 (key) > 0.

#: is 3.740 in wide and 2.095 in (rows 2 and 3) or 3.600 in (row 1) tall; the box
#: runs from the panel's label column to the cell's right margin and from the
#: bottom offset up to a fixed pad under the cell top, so the nine boxes share
#: their right edge and their row-wise top and bottom lines.  A *square* box in a
#: 3.740 x 2.095 cell can never fill it -- it leaves ~2.4 in of white on the right,
#: which is what the squares of the first draft looked like.  The strip under a
#: lower box (0.520 in) carries the x tick labels, the x title and, where the
#: panel has one, **its own key** (H and I species, G effect, E the sweep key).
BOX_BOTTOM_IN = 0.700            
                               # under the box drops from 0.72 in to 0.28 in, so the box fills its cell and the
                               # white band between two rows shrinks from ~0.94 in to ~0.50 in
BOX_TOP_PAD_IN = 0.075
BOX_H_LOW_IN = 1.460              # the six cells of rows 2 and 3 (= 2.195 - 0.280 - 0.075)
ROW1_BOX_BOTTOM_IN = 0.414        # row 1: the 23-row matrix of A sets the band,
ROW1_BOX_H_IN = 2.620             # and C's second (FP-rate) axis caps it: its
                                  # ticks and title need 0.33 in above the box
BOX_SIDE_A_IN = 3.420

#: canvas y (inches) of the lower edge of a key printed under a panel
KEY_BOTTOM_IN = 0.020


def cell_rect(fig: plt.Figure, left_in: float, *, height_in: float,
              bottom_in: float = BOX_BOTTOM_IN,
              top_pad_in: float = BOX_TOP_PAD_IN,
              right_in: float = 0.060) -> list[float]:
    """The plotting area of a cell-filling panel, in figure fractions.

    It runs from ``left_in`` (the end of the label column) to the canvas' right
    margin and from ``bottom_in`` to ``top_pad_in`` under the canvas top, clamped
    to the canvas so a box can never leave it.
    """
    win, hin = fig.get_figwidth(), fig.get_figheight()
    width = max(win - left_in - right_in, 0.4)
    height = min(height_in, hin - bottom_in - top_pad_in)
    return [left_in / win, bottom_in / hin, width / win, height / hin]


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
