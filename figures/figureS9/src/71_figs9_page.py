#!/usr/bin/env python
"""71 -- compose the nine panels onto the A4 landscape (842.4 x 595.44 pt) S8 page.

The panels were drawn at their final print size on a 3 x 3 grid, so pypdf only
translates them (no scaling): the point sizes in the panel PDFs are the printed
sizes.  ``compose_page`` asserts that every placement stays inside the page and
that no two panels overlap.

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/analysis/figS8_rebuild/scripts/71_figs9_page.py
"""
from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import s9_style as S  # noqa: E402  (also puts src/harmonisation on sys.path)
from common.panelpage import Panel, compose_page, pdf_to_png  # noqa: E402

PAGE = S.FIGS / "FigureS8_rev.pdf"
PNG = S.FIGS / "FigureS8_rev.png"


def main() -> None:
    placements = []
    for key in S.PIECE:
        w, h = S.PIECE[key]
        path = S.PANELS / f"figS8{key}.pdf"
        if not path.exists():
            raise SystemExit(f"missing panel {path} -- run 70_figs9_panels.py")
        x, y = S.PLACE[key]
        placements.append((Panel(name=key, path=path, width_in=w, height_in=h),
                           x, y))
    compose_page(S.PAGE_IN, placements, PAGE)
    pdf_to_png(PAGE, PNG, dpi=300)
    stamp = datetime.now().strftime("%H:%M:%S")
    lines = [f"[{stamp}] page {S.PAGE_IN[0]:.4f} x {S.PAGE_IN[1]:.4f} in "
             f"({S.PAGE_IN[0] * 72:.3f} x {S.PAGE_IN[1] * 72:.3f} pt)"]
    for key in S.PIECE:
        w, h = S.PIECE[key]
        x, y = S.PLACE[key]
        lines.append(f"  {key}: piece {w:.4f} x {h:.4f} in at "
                     f"({x:.4f}, {y:.4f}) in")
    lines.append(f"  wrote {PAGE.relative_to(S.PROJECT)} and "
                 f"{PNG.relative_to(S.PROJECT)}")
    (S.LOGS / "71_page.log").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
