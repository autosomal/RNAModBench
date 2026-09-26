#!/usr/bin/env python
"""71 -- redraw Figure 1 panel A (the dataset schematic) at the printed size.

Why a redraw and not a crop
---------------------------
68 cuts panel A out of the *published* Figure 1 as vectors.  That strip is
515 pt wide on a 595 pt page and has to be placed into the 459.52 pt wide A
strip of the rebuilt page, i.e. it is scaled by 0.89 -- and its own text was
already drawn small (3.3 - 8 pt on the full page).  Printed, the smallest
glyphs land at **3.3 pt**, which is both unreadable (it reads as mojibake ) and
below the 7 pt floor.  No placement can fix that: the artwork would have to be
scaled by ~2x to comply, which does not fit the page.

So panel A is redrawn here, in vectors, **at the size it is printed**:
every string is >= 7.5 pt, nothing is scaled afterwards (69 places it 1:1).
This is also what keeps the strip fully vector -- converting the published
artwork's Lab/CMYK colours with Ghostscript instead was measured on 2026-09-23
to keep the vector paths but (a) rasterise 12 spots (transparency/shading
fallbacks), (b) leave its WinAnsi fonts without ToUnicode (" mojibake " persists) and
(c) still print at 2.5 pt; see tables/fig1a_type_palette.tsv.
The redraw keeps the published schematic's information structure, measured off
the published page with ``pdftotext -bbox``:

* section headers "Datasets" / "Epigenetic modification" /
  "Pre-processing" / "Tools comparison";
* three datasets, each wild type (blue) versus the modification-deficient
  treatment (orange) -- A. thaliana, E. coli, HeLa;
* the four modification classes (m6A, m5C, Psi/m1Psi, Nm);
* pre-processing: the RRACH sequence ACTGGACTCTCGAGGA (the DRACH adenosine
  highlighted in the m6A orange) and the raw-current trace (Current / Time).

The published colours are kept (m6A orange ``#F5B264``, unmodified/wild type
blue ``#3778A0``); unlike the crop the strip is written in RGB, which also
removes the CMYK shift the PDF version showed.

Outputs
-------
``figures/figure1/figures/panels/Figure1_A_redrawn.{pdf,png}``
(pagesize exactly ``STRIP_W x STRIP_H`` pt, so 69 can place it with scale 1)

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    src/harmonisation/scripts/71_fig1_A_redraw.py
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
from matplotlib.patches import FancyBboxPatch, Circle

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from common import config as C                        # noqa: E402
from common.figstyle import apply as figstyle_apply   # noqa: E402
from common.manifest import setup_logger              # noqa: E402

OUT = (_RB / "figures/figure1")
FIG, LOG = (_RB / "figures/figure1/figures"), (_RB / "figures/figure1/logs")
PANEL_DIR = FIG / "panels"
A_PDF = PANEL_DIR / "Figure1_A_redrawn.pdf"
A_PNG = PANEL_DIR / "Figure1_A_redrawn.png"

#: must match 67's A_RECT width/height -- the strip is placed 1:1 by 69.
#: 2026-09-23: the page grew to the full text width (506.4591 pt) and the strip
#: with it, 486.4591 x 121.0 pt; the strings keep their printed size.
STRIP_W, STRIP_H = 486.4591, 121.0
MIN_PT = 7.5

ORANGE = "#F5B264"      # m6A / modification-deficient treatment
BLUE = "#3778A0"        # wild type / unmodified
GREY = "#4D4D4D"
FAM = {"m6A": ORANGE, "m5C": "#2E86AB", "\u03a8": "#A23B72",
       "m1\u03a8": "#6A994E", "Nm": "#7F7F7F"}
SEQ = "ACTGGACTCTCGAGGA"
SEQ_DRACH_A = 5         # the adenosine of GGACT

FS_HEAD = 8.5           # section headers
FS_TEXT = 7.5           # everything else (printed size, no scaling)


def style() -> None:
    figstyle_apply()
    mpl.rcParams.update({
        "font.size": FS_TEXT,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "figure.dpi": 150,
    })


def pill(ax, x, y, w, h, text, *, fc="white", ec=GREY, fs=FS_TEXT,
         italic=False, tc="black"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=2.2",
                                linewidth=0.7, edgecolor=ec, facecolor=fc,
                                mutation_scale=1.0, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=tc, style="italic" if italic else "normal", zorder=3)


def dot(ax, x, y, color, r=2.6):
    ax.add_patch(Circle((x, y), r, facecolor=color, edgecolor="none", zorder=3))


def datasets(ax, x0, y0, y1):
    """three datasets, wild type (blue) versus treatment (orange)"""
    ax.text(x0, y1 + 4.0, "Datasets", ha="left", va="bottom", fontsize=FS_HEAD,
            fontweight="bold", color="black")
    names = ["A. thaliana", "E. coli", "HeLa"]
    h, gap = 15.0, 7.0
    ys = [y1 - 2.0 - i * (h + gap) - h for i in range(len(names))]
    w = 70.0
    for y, name in zip(ys, names):
        pill(ax, x0, y, w, h, name, italic=True)
        dot(ax, x0 + w + 7.0, y + h / 2, BLUE)
        dot(ax, x0 + w + 17.0, y + h / 2, ORANGE)
    # shared key, directly under the column
    ky = y0 + 1.0
    dot(ax, x0 + 3.0, ky + 3.0, BLUE)
    ax.text(x0 + 8.0, ky, "wild type", ha="left", va="bottom", fontsize=FS_TEXT)
    dot(ax, x0 + 56.0, ky + 3.0, ORANGE)
    ax.text(x0 + 61.0, ky, "deficient", ha="left", va="bottom", fontsize=FS_TEXT)


def modifications(ax, x0, y0, y1):
    ax.text(x0, y1 + 4.0, "Epigenetic modification", ha="left", va="bottom",
            fontsize=FS_HEAD, fontweight="bold", color="black")
    labels = list(FAM.items())                      # (name, colour)
    w, h, gap = 50.0, 15.0, 7.0
    for i, (name, colour) in enumerate(labels):
        col, row = i % 2, i // 2
        x = x0 + col * (w + gap)
        y = y1 - 2.0 - row * (h + gap) - h
        pill(ax, x, y, w, h, name, fc=colour, ec="none", tc="white",
             fs=FS_TEXT)


def preprocessing(ax, x0, y0, y1, w):
    ax.text(x0, y1 + 4.0, "Pre-processing", ha="left", va="bottom",
            fontsize=FS_HEAD, fontweight="bold", color="black")
    # --- RRACH sequence: one box per base, the DRACH A in the m6A orange --- #
    bw, bh = 7.0, 13.0
    seq_w = bw * len(SEQ)
    sy = y1 - 2.0 - bh
    for i, base in enumerate(SEQ):
        hit = i == SEQ_DRACH_A
        ax.add_patch(FancyBboxPatch((x0 + i * bw, sy), bw, bh,
                                    boxstyle="round,pad=0.3,rounding_size=1.0",
                                    linewidth=0.5,
                                    edgecolor=ORANGE if hit else "#BFBFBF",
                                    facecolor=ORANGE if hit else "white",
                                    mutation_scale=1.0, zorder=2))
        ax.text(x0 + i * bw + bw / 2, sy + bh / 2, base, ha="center",
                va="center", fontsize=FS_TEXT,
                color="black" if hit else GREY, zorder=3)
    # --- raw current trace: Current (y) against Time (x) ------------------- #
    tx, ty, tw, th = x0 + 4.0, y0 + 10.0, min(w - 14.0, 66.0), 30.0
    t = np.linspace(0, 1, 220)
    sig = (0.30 * np.sin(2 * np.pi * 3.2 * t)
           + 0.12 * np.sin(2 * np.pi * 11.0 * t + 0.7)
           + 0.05 * np.random.default_rng(20260923).normal(size=t.size))
    ax.plot(tx + t * tw, ty + th / 2 + sig * th * 1.6, color=BLUE, lw=0.8,
            zorder=3)
    ax.plot([tx, tx], [ty, ty + th], color=GREY, lw=0.7)
    ax.plot([tx, tx + tw], [ty, ty], color=GREY, lw=0.7)
    ax.text(tx - 2.5, ty + th / 2, "Current", ha="right", va="center",
            fontsize=FS_TEXT, rotation=90)
    ax.text(tx + tw / 2, ty - 3.0, "Time", ha="center", va="top",
            fontsize=FS_TEXT)


def draw() -> plt.Figure:
    style()
    fig = plt.figure(figsize=(STRIP_W / 72.0, STRIP_H / 72.0))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, STRIP_W)
    ax.set_ylim(0, STRIP_H)
    ax.set_axis_off()

    head_y = STRIP_H - 11.0
    body_y0, body_y1 = 1.0, head_y - 1.0

    datasets(ax, 0.0, body_y0, body_y1)
    modifications(ax, 116.0, body_y0, body_y1)
    preprocessing(ax, 240.0, body_y0, body_y1, 124.0)

    # "Tools comparison" heads the tool panels drawn below the strip (67)
    ax.text(392.0, head_y + 4.0, "Tools comparison", ha="left", va="bottom",
            fontsize=FS_HEAD, fontweight="bold", color="black")
    ax.annotate("", xy=(392.0 + 44.0, 2.0), xytext=(392.0 + 44.0, body_y1 - 6.0),
                arrowprops=dict(arrowstyle="-|>", color=GREY, lw=0.8))
    return fig


def main() -> None:
    logger = setup_logger("71_fig1_A_redraw", log_dir=LOG)
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    fig = draw()
    fig.savefig(A_PDF, bbox_inches=None, pad_inches=0.0)
    fig.savefig(A_PNG, dpi=300, bbox_inches=None, pad_inches=0.0)
    plt.close(fig)

    from pypdf import PdfReader
    page = PdfReader(str(A_PDF)).pages[0]
    w, h = float(page.mediabox.width), float(page.mediabox.height)
    if abs(w - STRIP_W) > 0.2 or abs(h - STRIP_H) > 0.2:
        raise SystemExit(f"A strip is {w:.2f} x {h:.2f} pt, expected "
                         f"{STRIP_W} x {STRIP_H}")
    smallest = min(FS_HEAD, FS_TEXT)
    if smallest < MIN_PT:
        raise SystemExit(f"font {smallest} pt is below the {MIN_PT} pt floor")
    logger.info("wrote %s  %.2f x %.2f pt  (min font %.1f pt, printed 1:1)",
                A_PDF.name, w, h, smallest)
    logger.info("wrote %s", A_PNG.name)


if __name__ == "__main__":
    main()
