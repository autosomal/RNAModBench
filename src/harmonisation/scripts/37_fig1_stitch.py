#!/usr/bin/env python
"""Stitch the rebuilt Figure 1 panels B + C + D into one vector page (pypdf).

B = `figure1b_replicates/figures/Fig1B_tool_counts_replicates.pdf` (9.0 x 8.2 in)
C = `figures/figure1/figures/panels/Figure1_rev_C_curlcake_counts.pdf` (6.8 x 4.0 in)
D = `figures/figure1/figures/panels/Figure1_rev_D_rrach_counts.pdf`   (6.8 x 4.0 in)

Layout (pt; 1 in = 72 pt), one landscape page, B left and C/D stacked right:

    page 1070 x 670
    B  at ( 36, 114)  scale .75  ->  486 x 443
    C  at (544, 344)
    D  at (544,  36)
    letters B/C/D (Helvetica bold 18 pt) just above each panel's top-left corner

Vector composition via pypdf (conda-forge, installed into `benchmark-revision`
on 2026-09-20 after both TeX envs turned out to have broken format files); the
letter overlay is a transparent matplotlib page.  A 300 dpi PNG is rendered
from the composed page with pdftoppm.

Usage
-----
conda run -n benchmark-revision --no-capture-output python scripts/37_fig1_stitch.py
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
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

import matplotlib.pyplot as plt  # noqa: E402
from pypdf import PdfReader, PdfWriter, Transformation  # noqa: E402

from common.manifest import setup_logger  # noqa: E402

PROJECT = Path(str(_RB))
OUT = (_RB / "figures/figure1")
WORK = (_RB / "figures/figure1/logs")
FIGDIR = (_RB / "figures/figure1/figures")

B_PDF = ((_RB / "figures/figure1/inputs/figures/Fig1B_tool_counts_replicates.pdf"))
C_PDF = (_RB / "figures/figure1/figures/panels/Figure1_rev_C_curlcake_counts.pdf")
D_PDF = (_RB / "figures/figure1/figures/panels/Figure1_rev_D_rrach_counts.pdf")
PDF = (_RB / "figures/figure1/figures/Figure1_rev_BCD.pdf")
PNG = (_RB / "figures/figure1/figures/Figure1_rev_BCD.png")

PAGE_W, PAGE_H = 941, 478
LETTERS = (("B", 38, 446), ("C", 503.5, 446), ("D", 503.5, 237))


def letter_overlay(path: Path) -> None:
    """Transparent page carrying the B/C/D panel letters."""
    fig = plt.figure(figsize=(PAGE_W / 72, PAGE_H / 72))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    for txt, x, y in LETTERS:
        ax.text(x / PAGE_W, y / PAGE_H, txt, fontsize=18, weight="bold",
                family="Arial", ha="left", va="bottom")
    fig.savefig(path, transparent=True)
    plt.close(fig)


def main() -> None:
    logger = setup_logger("37_fig1_stitch")
    WORK.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    for f in (B_PDF, C_PDF, D_PDF):
        if not f.exists():
            raise SystemExit(f"missing panel: {f}")

    writer = PdfWriter()
    page = writer.add_blank_page(width=PAGE_W, height=PAGE_H)

    def place(src: Path, x: float, y: float, scale: float = 1.0) -> None:
        src_page = PdfReader(src).pages[0]
        op = Transformation().scale(scale, scale).translate(x, y)
        page.merge_transformed_page(src_page, op)

    place(B_PDF, 36, 36, 0.8362)
    place(C_PDF, 501.5, 258.9)
    place(D_PDF, 501.5, 50)

    overlay = (_RB / "figures/figure1/logs/fig1bcd_letters.pdf")
    letter_overlay(overlay)
    page.merge_page(PdfReader(overlay).pages[0])

    with open(PDF, "wb") as fh:
        writer.write(fh)
    subprocess.run(["pdftoppm", "-png", "-r", "300", "-singlefile",
                    str(PDF), str(PNG.with_suffix(""))], check=True)
    logger.info("wrote %s (+ .png)", PDF)


if __name__ == "__main__":
    main()
