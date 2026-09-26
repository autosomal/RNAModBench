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
PAGE_PT = (6.95 * 72, 8.15 * 72)
PLACE = {"fig8A_counts": (0.05, 5.08), "fig8B_fpr": (3.60, 5.08),
         "fig8C_ppv": (0.05, 2.01), "fig8D_tradeoff": (3.60, 2.01),
         "fig8E_guitar": (0.05, 0.05)}


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
