#!/usr/bin/env python
"""Figure S9 style tokens -- page geometry, type scale, palette.

Since 2026-09-30 (author decision) the page is an **enlarged 3 x 3 sheet**: A | B | C on
top, then D | E | H, then F | G | I.  The revision asked for the six lower panels to reach the height of the top three, for
a little more width per cell, and, once the sheet had grown, for the type to catch up.  The A4 landscape sheet could not hold any of that, so the sheet is
enlarged instead of shrinking row 1, and every size except ``row_a`` is one step
larger than the 2026-09-27 A4 page (see ``FS``).

2026-10-09 (author decision), of this page: " the nine panels can be a little wider, and square " -- the nine
**plotting areas are squares** now.  The box height cannot go down (A's 23-row
label pitch pins it at ``BOX_H_LOW_IN`` = 2.675 in), so the square is built by
widening: the shared label column is **frozen** as ``LABEL_COL_IN`` and the
column width -- and with it the page width -- is derived from it instead of the
other way round.  The sheet is composed 1:1, so the page itself keeps every type
size; in the SI album (7.008 in of printed width) the scale moves 0.565 -> 0.497,
the trade-off the author accepted on 2026-10-09.

Layout arithmetic (inches, 1:1 print size)::

    page   1015.776 x 792.0 pt  =  14.108 x 11.000
    margins 0.10 all round, gutters 0.060 (vertical and between the three columns)
    row heights 3.580 (A|B|C) + 3.550 (D|E|H) + 3.550 (F|G|I) + 2*0.060 = 11.000
    column width  frozen label column 1.861 + 0.060 + box side 2.675 = 4.596
    nine boxes  2.675 x 2.675 in -- **squares** (2026-10-09)

Every panel is drawn on its own print-size canvas and the page is composed 1:1
with pypdf (no scaling), so the type sizes below are the printed sizes.

House rules honoured here: Arial only, no grid, no in-panel annotation text
(only bold panel letters and entries at 8 pt or larger since the 2026-09-30 type
scale), vector PDF + 300 dpi PNG at the exact page size (never
``bbox_inches="tight"``).
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

#: (inside which rows 2 and 3 could never give the six lower panels more than
#: ~1.8 in of box height) to 12.400 x 10.700, and the same day deepened twice
#: with the type scale -- 10.900 for the two key strips (0.1 in each) and 11.000
#: for row 1, whose band carries C's FP-rate twin axis above its box and C's own
#: x title below it.

#: that the nine boxes come out square -- see ``LABEL_COL_IN`` / ``BOX_SIDE_IN``.
M_LEFT = M_RIGHT = 0.100
M_TOP = M_BOTTOM = 0.100
GUT_H = 0.060                                   # between the three columns
GUT_V = 0.060                                   # between the three rows
N_COL = 3


#: the value measured on 2026-09-30** -- A's 23 model names at ``FS["row_a"]``
#: plus 0.09 in of pad need 1.8562 in, F's longest tool name at ``FS["row"]``
#: plus the 0.32 in panel-letter provision needs 1.8607 in.  It is a constant
#: (not a measurement) because the page width is derived from it:
#: ``70_figs9_panels.column_left`` re-measures both needs at build time and fails
#: loudly when either outgrows this number -- that is the event which forces a
#: re-derivation of the page.
LABEL_COL_IN = 1.861
#: right margin kept inside every cell canvas: a box ends 0.060 in short of the
#: cell's right edge (one number, three users -- ``cell_rect``'s default
#: ``right_in``, the inline rect of the row-1 bar panels, and
#: ``draw_sweep_key``'s box centre)
BOX_RIGHT_IN = 0.060
#: side of the nine **square** boxes.  **The floor is A's 23-row pitch**: 23
#: labels at 8 pt need 2.675 in of box height, so a square is only reachable by
#: widening the sheet, never by lowering the box (``_check_geometry`` asserts
#: ``BOX_SIDE_IN == BOX_H_LOW_IN``).
BOX_SIDE_IN = 2.675


#: it the page width, follows from the frozen label column and the box side
COL_W = LABEL_COL_IN + BOX_RIGHT_IN + BOX_SIDE_IN           # 4.596
W_FULL = N_COL * COL_W + (N_COL - 1) * GUT_H                # 13.908
PAGE_IN = (M_LEFT + W_FULL + M_RIGHT, 792.0 / 72.0)         # 14.108 x 11.000


#: letters are anchored to their **own canvas' left edge** at this single
#: inset, so the page reads as a 3 x 3 grid (the canvases are equal per
#: column and the columns line up across the three rows).  This replaces the
#: per-panel placements of 2026-09-28/09-30, which pushed every letter just
#: left of its own label block (and D/E next to their boxes), leaving the
#: column x positions 0.02 .. 0.99 in apart.  The visual intent of those
#: rounds survives -- each letter still sits in its panel's top-left corner,
#: clear of the ink: the measured headroom in the letter band is large (the
#: closest label block starts 1.31 in in, so up to 0.97 in would fit), and
#: 0.030 in leaves C/F/G where they already were.
LETTER_IN = 0.030                                        # inches, per canvas


#: (A | B | C, then D | E | H, then F | G | I).  Row 1 is the tallest (3.580) and
#: draws its band in ``s9_panels._draw_bar_panel`` (``ROW1_BOTTOM_FRAC`` /
#: ``ROW1_BOX_FRAC`` below): the band holds C's FP-rate twin axis above its box
#: and C's own x title below it.  Rows 2 and 3 carry 3.550 in each (the key strip
#: deepened with the type scale), so their boxes come out 2.675 in tall
#: (BOX_H_LOW_IN, = r2 - 0.800 - 0.075).  Since 2026-10-09 row 1 draws a box of
#: exactly the same 2.675 in, so all nine boxes are equal squares -- their width
#: comes from ``COL_W`` -- while every key still prints under its own box.
#: 0.100 + 3.580 + 0.060 + 3.550 + 0.060 + 3.550 + 0.100 = 11.000.
ROW_H = {"r1": 3.580, "r2": 3.550, "r3": 3.550}

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
# type scale (points, printed 1:1)
# --------------------------------------------------------------------------- #

#: enlargement the sheet had grown by ~20 % while the type had not, so every size
#: **except row_a** was raised by ~20 %.  A's 23-row labels stay at 8 pt because
#: their row pitch (2.675 in / 23 = 8.37 pt, the square side of 2026-10-09)
#: leaves no room at all, and the letters of C/F/G end up at the canvas' left
#: edge as a consequence (see ``letter_near``).  ``"title"``/``"stripe"`` (panel
#: titles, page stripe) were never used and are gone; the layout gate carries its
#: own 7 pt default floor, while the smallest size the page actually prints is 8 pt
#: (the sweep key).
FS = {
    "tick": 11.5,
    "tick_sm": 11.5,
    "tick_xs": 10.5,
    "axis": 12.5,      #: unused on this page (every caller passes axis_sm), kept in step
    "axis_sm": 12.5,
    "row": 9.5,
    "level": 9.5,      
    #: locked: 2.675 in / 23 rows = 8.37 pt of pitch -- 9.5 pt would overlap
    "row_a": 8.0,
    "legend": 9.5,     #: only the grid_key default; every key of this page passes its size
    "letter": 28.0,
    "key": 9.5,
    "key_sm": 9.5,
}
#: the smallest size the page prints (the thirteen-entry sweep key); documents the
#: drawing code -- the layout gate's own default floor is still 7 pt
MIN_PT = 8.0

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

 2026-09-27 (author decision): every panel is a **boxed square** (""),
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



#: is 4.596 in wide (2026-10-09) and 3.550 in (rows 2 and 3) or 3.580 in (row 1)
#: tall; the box runs from the panel's label column to the cell's right margin and
#: from the bottom offset up to a fixed pad under the cell top, so the nine boxes
#: share their right edge and their row-wise top and bottom lines.  The strip
#: under a lower box carries its x tick labels, its x title and, where the panel
#: has one, **its own key** (H and I the species, G the effect, D/E the two halves
#: of the sweep key), because **a key prints under the panel it belongs to**
#: (house rule).
#: Row 1 draws its band in ``s9_panels._draw_bar_panel`` at the canvas fractions
#: ``ROW1_BOTTOM_FRAC`` / ``ROW1_BOX_FRAC`` (bottom 0.474 in, box 2.675 in, 0.431
#: in kept above the box for C's FP-rate twin axis; 2026-10-09: the box is the
#: same square as the six lower ones, so the strip above it keeps growing with
#: the band).

#: enlarged type -- it holds the x ticks (11.5 pt), the x title (12.5 pt) and the
#: sweep key, whose three rows need 0.024 + 3*0.117 in at 8 pt
BOX_BOTTOM_IN = 0.800
BOX_TOP_PAD_IN = 0.075
#: the six cells of rows 2 and 3 -- **derived, never a literal again**: a stale
#: literal is silently clamped by ``cell_rect`` (a small box under a fat white
#: band) and neither the layout gate nor the verifier would see it.
BOX_H_LOW_IN = ROW_H["r2"] - BOX_BOTTOM_IN - BOX_TOP_PAD_IN          # 2.675


#: the six lower ones -- the ceiling above it (C's FP-rate twin axis plus its
#: title) is what ``_check_geometry`` guards
ROW1_BOTTOM_IN = 0.474
ROW1_BOTTOM_FRAC = ROW1_BOTTOM_IN / ROW_H["r1"]
ROW1_BOX_FRAC = BOX_SIDE_IN / ROW_H["r1"]                            # 0.74721

#: canvas y (inches) of the lower edge of a key printed under a panel
KEY_BOTTOM_IN = 0.020


def _check_geometry() -> None:
    """Fail loudly when the page budget stops adding up (2026-09-30).

    The layout was re-derived by hand through many rounds (2026-09-25 …
    2026-09-30, and the square of 2026-10-09), and the one constant that was
    forgotten in a later edit is silently swallowed by ``cell_rect``'s clamp.
    These checks make a stale row height, column width or box side impossible to
    write down -- including the one that ties the nine boxes to **squares**.
    """
    assert ROW_H["r2"] == ROW_H["r3"], "rows 2 and 3 must share one height"
    total = M_TOP + M_BOTTOM + ROW_H["r1"] + ROW_H["r2"] + ROW_H["r3"] + 2 * GUT_V
    assert abs(total - PAGE_IN[1]) < 1e-9, (total, PAGE_IN[1])
    want = (PAGE_IN[0] - M_LEFT - M_RIGHT - (N_COL - 1) * GUT_H) / N_COL
    assert abs(COL_W - want) < 1e-12, (COL_W, want)
    assert abs(BOX_H_LOW_IN - (ROW_H["r2"] - BOX_BOTTOM_IN - BOX_TOP_PAD_IN)) < 1e-12
    assert BOX_BOTTOM_IN > 0 and BOX_TOP_PAD_IN > 0

    #: ``PAGE_IN`` that agrees with itself (the check above), so the square has
    #: its own ties: the column is exactly the frozen label column plus the right
    #: margin plus the box side, the box side equals the rows-2/3 box height, and
    #: row 1's band draws the same side inside its taller cell.
    assert abs(COL_W - (LABEL_COL_IN + BOX_RIGHT_IN + BOX_H_LOW_IN)) < 1e-12, \
        (COL_W, LABEL_COL_IN + BOX_RIGHT_IN + BOX_H_LOW_IN)
    assert abs(BOX_SIDE_IN - BOX_H_LOW_IN) < 1e-12, (BOX_SIDE_IN, BOX_H_LOW_IN)
    assert abs(ROW1_BOX_FRAC * ROW_H["r1"] - BOX_H_LOW_IN) < 1e-12
    assert abs(ROW1_BOTTOM_FRAC * ROW_H["r1"] - ROW1_BOTTOM_IN) < 1e-12
    #: C's FP-rate twin axis prints in the strip above row 1's box (its title,
    #: 12.5 pt, plus the 11.5 pt tick labels): the strip must not shrink below
    #: the 0.42 in it has needed since 2026-09-30
    strip = ROW_H["r1"] - ROW1_BOTTOM_IN - BOX_SIDE_IN
    assert strip >= 0.42, f"row 1 keeps only {strip:.3f} in above its box"


_check_geometry()


def cell_rect(fig: plt.Figure, left_in: float, *, height_in: float,
              bottom_in: float = BOX_BOTTOM_IN,
              top_pad_in: float = BOX_TOP_PAD_IN,
              right_in: float = BOX_RIGHT_IN) -> list[float]:
    """The plotting area of a cell-filling panel, in figure fractions.

    It runs from ``left_in`` (the end of the label column) to the canvas' right
    margin and from ``bottom_in`` to ``top_pad_in`` under the canvas top.

    2026-10-09 (author decision): a box that does not fit its canvas is now a **hard
    failure**.  Until this round the two clamps below silently shrank the box
    (``max(..., 0.4)`` / ``min(...)``), which no gate and no verifier could see:
    a stale column width would have printed nine 0.4 in boxes under fat white
    bands and looked like a styling accident.
    """
    win, hin = fig.get_figwidth(), fig.get_figheight()
    width = win - left_in - right_in
    avail = hin - bottom_in - top_pad_in
    if width < 0.4:
        raise ValueError(f"box would be only {width:.3f} in wide on a "
                         f"{win:.3f} in canvas (left_in {left_in:.3f})")
    if height_in > avail + 1e-9:
        raise ValueError(f"box wants {height_in:.3f} in of height but the "
                         f"canvas keeps {avail:.3f} in")
    return [left_in / win, bottom_in / hin, width / win, height_in / hin]


def letter(fig: plt.Figure, key: str, x: float | None = None) -> None:
    """Bold panel letter in the top-left corner of a panel canvas.

    ``x`` (figure fraction) is only for a canvas that carries two lettered
    facets (2026-09-25 design; no panel of this page uses it today).  All nine
    panels take the default -- one ``LETTER_IN`` inset from their own canvas,
    which is what makes the columns line up (2026-10-09 revision request).
    """
    fig.text((LETTER_IN / COL_W) if x is None else x, 0.994, key,
             fontsize=FS["letter"], fontweight="bold", ha="left", va="top")


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
