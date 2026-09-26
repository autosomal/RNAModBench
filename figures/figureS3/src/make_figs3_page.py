#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Assemble the rebuilt Figure S3 as a single page.

Page geometry is identical to the submitted ``sup3.pdf`` (595.276 x 633.598
pt).  Rows:
  A  GLORI replicate-overlap venns (unified criterion, both replicates > 0.1)
  B  PPV vs. GLORI (2 bp), per independent unit (revision evaluation layer)
  C  per-replicate modification-ratio agreement (the companion analysis recipe,
     re-rendered at page geometry from their cached tables)
Every row is drawn by the same code that produces its standalone panel (the
``page=True`` drawing mode places the row at its absolute page position), so
the page can never drift from the panels.

Output: figures/figureS3/FigureS3_rev.pdf / .png
(the submitted sup3.pdf is never touched).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import matplotlib.pyplot as plt
import pandas as pd

import make_s3a_venn
import make_s3b_ppv
import make_s3c_reuse_other
import s3_common as sc


def main() -> int:
    #: the sibling module sets its own rcParams at import -> our style last
    other = make_s3c_reuse_other.load_other_recipe()
    sc.apply_page_style()

    print("[compute] panel A overlaps ...")
    overlaps = sc.compute_overlap()
    sc.check_overlap(overlaps)
    print("[compute] panel B per-unit PPV ...")
    table_b = make_s3b_ppv.load_confusion()
    print("[compute] panel C (the companion analysis cached tables) ...")
    sites = pd.read_csv(make_s3c_reuse_other.MATCHED, sep="\t")
    summary = pd.read_csv(make_s3c_reuse_other.SUMMARY, sep="\t")
    tools, colors = list(other.TOOLS), dict(other.TOOL_COLOR)

    fig = plt.figure(figsize=(sc.PAGE_W / 72.0, sc.PAGE_H / 72.0))

    # ---- row A -----------------------------------------------------------
    ax_a = fig.add_axes((0.0, 1.0 - 185.0 / sc.PAGE_H, 1.0, 185.0 / sc.PAGE_H))
    make_s3a_venn.draw_panel_a(ax_a)

    # ---- rows B and C at their absolute page positions -------------------
    make_s3b_ppv.draw_row_b(fig, table_b, y_top=200.0, row_h=245.0,
                            page=True, titles=True)
    make_s3c_reuse_other.draw_row_c(fig, sites, summary, tools, colors,
                                    titles=True, y_top=420.0,
                                    row_h=sc.PAGE_H - 420.0, page=True)

    # ---- panel letters ----------------------------------------------------
    for letter, (x, y) in sc.LETTER_POS.items():
        sc.panel_letter(fig, letter, x, y)

    sc.OUT_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(sc.OUT_DIR / f"FigureS3_rev.{ext}", facecolor="white")
    plt.close(fig)
    print(f"[write] {sc.OUT_DIR / 'FigureS3_rev.pdf'} (+png)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
