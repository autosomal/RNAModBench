#!/usr/bin/env python3
"""68 -- assemble the four Figure 8 pieces on an A4-landscape page (1:1, vector)."""

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

from pypdf import PdfReader, PdfWriter, Transformation

FIGS = Path(str(_RB / "figures/figure8/figures"))
PANELS = (_RB / "figures/figure8/figures/panels")

#: 3.02 in at y = 5.24 and 2.16) and the page keeps its 6.95 in width (the height
#: gained 0.09 in on the third pass, below), but each of
#: the two upper rows is three equal columns of 2.18333 in with a 0.15 in gutter
#: (fig8_style.COL_IN / GUTTER_IN): row 1 = A | C's PPV dot plot | C's calls bar,
#: row 2 = B's Curlcake sub-plot | B's Human sub-plot | D.  The pieces span one
#: column (A, D) or two (B, C); their internal boxes put the sub-axes on the
#: column edges, so the row reads as three equal thirds.  (2026-10-08 second pass:
#: inside C the PPV|calls split was rebalanced to ~1:1.5 and left the column
#: boundary; C's own edges and letter did not move, and row 2 is untouched.)

#: equal thirds.  A's piece is ROW1_A_IN = 2.78 in and C's ROW1_C_IN = 3.92 in
#: (fig8_style), so C starts at MARGIN + ROW1_A_IN + GUT = 2.98 in and its right
#: edge still lands on the page's right margin (2.98 + 3.92 = 6.90 in).  Row 2 and
#: E keep the equal columns, which is why their PLACE entries still use COL.
from fig8_style import (COL_IN as COL, GUTTER_IN as GUT,      # noqa: E402
                        PAGE_MARGIN_IN as MARGIN, ROW1_A_IN)

#: clears the panel letters (0.05 -> 0.14 in); every PLACE below stays untouched.
PAGE_PT = (6.95 * 72, 8.40 * 72)
PLACE = {"fig8A_counts": (MARGIN, 5.24),
         "fig8C_ppv": (MARGIN + ROW1_A_IN + GUT, 5.24),
         "fig8B_fpr": (MARGIN, 2.16),
         "fig8D_tradeoff": (MARGIN + 2 * (COL + GUT), 2.16),
         "fig8E_guitar": (MARGIN, 0.05)}


def main() -> None:
    writer = PdfWriter()
    page = writer.add_blank_page(*PAGE_PT)
    for name, (x_in, y_in) in PLACE.items():
        src = PdfReader(str(PANELS / f"{name}.pdf"))
        box = src.pages[0].mediabox
        page.merge_transformed_page(
            src.pages[0], Transformation().translate(tx=x_in * 72, ty=y_in * 72))
        print(f"[assemble] {name}: {box.width / 72:.3f} x {box.height / 72:.3f} in "
              f"at ({x_in:.2f}, {y_in:.2f})")
    out = (_RB / "figures/figure8/figures/Figure8_rev.pdf")
    with open(out, "wb") as fh:
        writer.write(fh)
    print(f"[assemble] wrote {out}  page {PAGE_PT[0]} x {PAGE_PT[1]} pt")


if __name__ == "__main__":
    main()
