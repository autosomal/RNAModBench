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

#: unchanged (6.95 x 8.31 in; rows of 3.02 in at y = 5.24 and 2.16), but each of
#: the two upper rows is three equal columns of 2.18333 in with a 0.15 in gutter
#: (fig8_style.COL_IN / GUTTER_IN): row 1 = A | C's PPV dot plot | C's calls bar,
#: row 2 = B's Curlcake sub-plot | B's Human sub-plot | D.  The pieces span one
#: column (A, D) or two (B, C); their internal boxes put the sub-axes on the
#: column edges, so the row reads as three equal thirds.
MARGIN, COL, GUT = 0.05, 2.18333, 0.15
PAGE_PT = (6.95 * 72, 8.31 * 72)
PLACE = {"fig8A_counts": (MARGIN, 5.24),
         "fig8C_ppv": (MARGIN + COL + GUT, 5.24),
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
