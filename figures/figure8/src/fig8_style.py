#!/usr/bin/env python3
"""Shared style helpers for the per-panel Figure 8 / Figure S9 construction.

Every panel is drawn on its own canvas at its final print size and exported as a
standalone piece (``figures/panels/<name>.pdf`` + ``.png``); ``60_assemble_fig8.py``
and ``62_assemble_figS9.py`` then place the pieces on the page with pypdf.  Sizes
here are therefore *printed* point sizes -- no scaling happens at assembly time.
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
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = Path(str(_RB))
OUT = (_RB / "figures/figure8")
TABLES = (_RB / "figures/figure8/tables")
FIGS = (_RB / "figures/figure8/figures")
PANELS = (_RB / "figures/figure8/figures/panels")

# ---- palette (REPLICATION_SPEC.md) ----------------------------------------- #
WT_COLOR, IVT_COLOR = "#4B81B8", "#E8A76B"
#: "inosine" used to be #F18F01, a second orange that read as the same series as
#: the non-DRACH m6A (#F5B264) in panels C and D; it is now the only green
FAM_COLOR = {"m6A_DRACH": "#E83C1E", "m6A": "#F5B264", "m5C": "#2E86AB",
             "pseU": "#A23B72", "inosine": "#4E9A51"}
FAM_LABEL = {"m6A_DRACH": "m6A DRACH", "m6A": "m6A (non-DRACH)",
             "m5C": "m5C", "pseU": "pseU", "inosine": "inosine + m6A"}
DOT_COLOR = "#1E888B"
GUIDE = "0.72"

STYLE = dict(font=8.5, tick=8.5, row=8.8, label=10.0, value=12.0, letter=15.0,
             legend=8.5, title=11.5, rod_alpha=0.40)

# ---- page geometry, inches (final print size, no scaling at assembly) -------- #
PAGE_IN = (6.95, 8.15)          # 500.4 x 586.8 pt, printed 1:1 (\textwidth)
PIECE = {"A": (3.30, 3.02), "B": (3.30, 3.02), "C": (3.30, 3.02),
         "D": (3.30, 3.02), "E": (6.85, 2.05)}

# ---- panel B: two sub-axes side by side in the same 3.30 x 3.02 in piece ----- #
# B1 (left) = unmodified Curlcake threshold scan, B2 (right) = unmodified HeLa
# IVT at the delivered >= 90 % cutoff; one shared legend sits under both.  They
# were stacked, which read as one panel with two unrelated y scales; each half
# is now 1.24 in wide and the shared y title sits at their common left edge.

#: 0.36 in to the left of its axes and were printed over B1's plot area -- the
#: two blocks were 0.18 in apart.  The gutter is now 0.69 in.



B_BOX = {"B1": [0.105, 0.280, 0.300, 0.595],
         "B2": [0.615, 0.280, 0.360, 0.595]}
#: the slanted family labels of B2 need ~0.10 in below the axes, which is the
#: band the shared legend would otherwise be printed on
B_YSHARE = [0.280, 0.875]        # common vertical extent of the two sub-axes
B_YLIM = {"B1": (0.35, 2.5e3), "B2": (3e-5, 12.0)}
B_YTICKS = {"B1": [1, 10, 100, 1000],
            "B2": [1e-4, 1e-3, 1e-2, 1e-1, 1.0]}
#: B2's labels in 10^-n form: two characters narrower each, which the 0.69 in
#: gutter absorbs with room to spare
B_YLABELS = {"B1": ["1", "10", "100", "1,000"],
             "B2": ["$10^{-4}$", "$10^{-3}$", "$10^{-2}$", "$10^{-1}$", "1"]}
B_FLOOR = {"B1": 0.5, "B2": 6e-5}   # open marker for a series with no call


def apply_style() -> None:
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": STYLE["font"],
        "axes.labelsize": STYLE["label"],
        "axes.titlesize": STYLE["title"],
        "xtick.labelsize": STYLE["tick"],
        "ytick.labelsize": STYLE["tick"],
        "legend.fontsize": STYLE["legend"],
        "axes.linewidth": 0.9,
        "xtick.major.width": 0.9,
        "ytick.major.width": 0.9,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "axes.grid": False,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "savefig.dpi": 300,
    })


def letter(fig, s: str, x: float = 0.004, y: float = 0.997) -> None:
    """Bold panel letter in the piece's top-left corner (figure coordinates)."""
    fig.text(x, y, s, fontsize=STYLE["letter"], fontweight="bold", ha="left",
             va="top")


def title(ax, s: str, pad: float = 6.0) -> None:
    ax.set_title(s, fontsize=STYLE["title"], fontweight="bold", pad=pad, loc="left")


def _fit(fig, ax, side: str, pad_in: float = 0.02,
         keep_right: float | None = None) -> None:
    """Shrink ``ax`` so tick/axis labels stay inside the canvas (no tight bbox).

    ``Axes.get_tightbbox`` clips to the figure, so the requirement is measured
    from the tick/axis label extents instead.  ``keep_right`` pins the right
    edge (figure fraction) after the left side has moved, so long row labels
    can never push the plotting area past the canvas edge.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    pos = ax.get_position()
    w_in, h_in = fig.get_size_inches()
    x0, y0 = pos.x0, pos.y0
    x1, y1 = pos.x1, pos.y1

    def wid(t) -> float:
        return t.get_window_extent(r).width / fig.dpi

    def hei(t) -> float:
        return t.get_window_extent(r).height / fig.dpi

    if side == "left":
        need = max([wid(t) for t in ax.get_yticklabels() if t.get_text()] or [0])
        need += 0.03                                        # tick pad
        lab = ax.yaxis.get_label()
        if lab.get_text():                                  # rotated y title
            need += wid(lab) + 0.02
        # get_window_extent() is measured on the open figure, whose font metrics
        # are not the PDF's, and it ignores labelpad: the DRACH row labels of
        # panel A and the y title of panel D ran 2-5 pt off the canvas.  The
        # estimate below is therefore corrected against the real tight bbox.
        need *= 1.06
        x0 = min(max((need + pad_in) / w_in, 0.02), 0.75)
        ax.set_position([x0, y0, x1 - x0, y1 - y0])
        fig.canvas.draw()
        left = ax.get_tightbbox(r).x0 / (fig.get_figwidth() * fig.dpi)
        if left < pad_in:
            x0 = min(x0 + (pad_in - left), 0.75)
        if keep_right is not None:
            x1 = max(x0 + 0.20, keep_right)
    elif side == "bottom":
        need = max([hei(t) for t in ax.get_xticklabels() if t.get_text()] or [0])
        need += 0.05                                        # tick length + pad
        lab = ax.xaxis.get_label()
        if lab.get_text():
            need += hei(lab) + 0.02
        y0 = min(max((need + pad_in) / h_in, 0.02), 0.6)
    elif side == "top":
        y1 = max(min(1 - pad_in, 0.995), 0.4)
    ax.set_position([x0, y0, x1 - x0, y1 - y0])
    fig.canvas.draw()


def fit_labels(fig, ax, left=True, bottom=True, top=False, pad_in: float = 0.02,
               keep_right: float | None = None) -> None:
    """Apply :func:`_fit` to the requested sides (labels must not be clipped)."""
    for side, on in (("left", left), ("bottom", bottom), ("top", top)):
        if on:
            _fit(fig, ax, side, pad_in, keep_right=keep_right)


def save_piece(fig, name: str, width: float, height: float, dpi: int = 300,
               fit: tuple = (), keep_right: float | None = None) -> None:
    """Export one panel piece at its exact print size (vector PDF + 300 dpi PNG).

    ``fit`` lists the sides whose tick/axis labels must be kept inside the canvas;
    the figure is resized first so the fit is computed for the final canvas.
    """
    PANELS.mkdir(parents=True, exist_ok=True)
    fig.set_size_inches(width, height)
    if fit:
        ax = fig.axes[0]
        fit_labels(fig, ax, left="left" in fit, bottom="bottom" in fit,
                   top="top" in fit, keep_right=keep_right)
    for ext in ("pdf", "png"):
        fig.savefig(PANELS / f"{name}.{ext}", dpi=dpi if ext == "png" else None)
    plt.close(fig)
    print(f"[piece] {PANELS / name}.pdf (+png)  {width} x {height} in", flush=True)
