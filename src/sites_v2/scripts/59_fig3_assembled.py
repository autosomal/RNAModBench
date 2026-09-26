#!/usr/bin/env python
"""Assemble Figure 3 (A-D) into one page, at the final printed size.

The four rows come from two backends and are merged as **vector** pages:

  A  Fig3A_tool_similarity_mds.pdf     matplotlib  (56_fig3a_tool_similarity_mds.py)
  B  Fig3B_modratio_wt_treatment.pdf   matplotlib  (57_fig3b_modratio_wt_treatment.py)
  C  Fig3C_metagene_species.pdf        R / Guitar  (58_fig3cd_guitar.R)
  D  Fig3D_metagene_tools.pdf          R / Guitar  (58_fig3cd_guitar.R)

Every row is drawn at exactly 6.66 in = 0.95 x \\textwidth of the Wiley USG
layout, so the merged page is 6.66 in wide and ~7.8 in tall -- the printed
height of the published Figure 3 (8.34 x 11.59 in artwork shown at
0.8 x \\textwidth = 5.61 x 7.79 in), i.e. replacing the file does not reflow the
manuscript beyond the \\includegraphics width.

pypdf places each row with a pure translation (no scaling, no rasterisation), so
text stays real text with embedded Arial in the merged PDF.

Checks reported here (and written to tables/fig3_assembly_geometry.tsv):
  * page size and the position of every row;
  * every row's width equals the page width (no hidden scaling);
  * the smallest ``Tf`` font size found in the merged page (must be >= 7.2 pt);
  * pdffonts embedding column (must not be "no" for any font).

Outputs -> 04_revision_analysis/fig3_revision/figures/Figure3_rev.{pdf,png}

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/src/sites_v2/scripts/59_fig3_assembled.py
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
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
from pypdf import PdfReader, PdfWriter, Transformation

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from common import config as C                                        # noqa: E402
from common.io_utils import write_table                               # noqa: E402
from common.manifest import setup_logger                              # noqa: E402

OUT = (_RB / "figures/figure3")
FIG, TAB, LOG = OUT / "figures", OUT / "tables", OUT / "logs"

ROWS = [
    ("A", "Fig3A_tool_similarity_mds.pdf", "matplotlib (56_fig3a)"),
    ("B", "Fig3B_modratio_boxplot.pdf", "matplotlib (57_fig3b)"),
    ("C", "Fig3C_metagene_species.pdf", "R/Guitar (58_fig3cd)"),
    ("D", "Fig3D_metagene_tools.pdf", "R/Guitar (58_fig3cd)"),
    ("L", "Fig3_legend_band.pdf", "matplotlib (61_legend_band)"),
]
W_IN = 6.66                # 0.95 x \textwidth
GAP_IN = 0.05              # default vertical gap between rows
# per-transition gaps: B-C and C-D sit closer, as requested (upper
# letters keep a little more air)
ROW_GAPS = {"A": 0.05, "B": 0.02, "C": 0.02, "D": 0.05}
MIN_FONT_PT = 7.2
PT = 72.0


def min_font_size(page) -> float:
    """Smallest effective text size in a content stream, in points.

    Two conventions are handled, because the rows come from two backends:
      * matplotlib  -- ``/F1 7.2 Tf`` with an unscaled text matrix;
      * cairo (R)   -- ``/f-1-0 1 Tf`` with the size carried by the text matrix
                       (``8.598633 0 0 -8.598633 x y Tm`` = 8.6 pt).
    The effective size of every ``Tj`` is therefore ``Tf_size * |Tm_scale|``,
    and the minimum over all text objects is returned.
    """
    data = page.get_contents()
    if data is None:
        return float("nan")
    txt = data.get_data().decode("latin-1", errors="ignore")
    tf, tm_scale = 1.0, 1.0
    sizes: list[float] = []
    for m in re.finditer(r"/([A-Za-z0-9#\-.]+)\s+([0-9]*\.?[0-9]+)\s+Tf"
                         r"|([0-9.\-]+)\s+([0-9.\-]+)\s+([0-9.\-]+)\s+([0-9.\-]+)"
                         r"\s+([0-9.\-]+)\s+([0-9.\-]+)\s+Tm"
                         r"|(Tj|TJ|'|\")", txt):
        if m.group(1) is not None:                       # Tf
            tf = float(m.group(2))
        elif m.group(3) is not None:                     # Tm
            a, b = float(m.group(3)), float(m.group(4))
            tm_scale = (a * a + b * b) ** 0.5
        else:                                            # a text-showing operator
            sizes.append(tf * tm_scale)
    sizes = [s for s in sizes if s > 0]
    return min(sizes) if sizes else float("nan")


def main() -> None:
    logger = setup_logger("59_fig3_assembled", log_dir=LOG)
    w_pt = W_IN * PT
    rows = []
    for letter, name, backend in ROWS:
        p = FIG / name
        if not p.exists():
            raise SystemExit(f"missing panel: {p} -- run its script first")
        pg = PdfReader(str(p)).pages[0]
        rows.append({"letter": letter, "file": name, "backend": backend,
                     "page": pg,
                     "w": float(pg.mediabox.width), "h": float(pg.mediabox.height),
                     "font_pt": min_font_size(pg)})
    # R (cairo) writes 479.0 pt for 6.66 in where matplotlib writes 479.52 pt:
    # a 0.07 % rounding difference, absorbed by centring each row, never by scaling.
    for r in rows:
        if abs(r["w"] - w_pt) > 2.0:
            raise SystemExit(f"row {r['letter']} is {r['w']:.1f} pt wide, expected "
                             f"{w_pt:.1f} pt -- re-run its script (no scaling here)")
        r["x_pt"] = (w_pt - r["w"]) / 2.0
        if r["font_pt"] < MIN_FONT_PT - 0.05:
            raise SystemExit(f"row {r['letter']} contains {r['font_pt']:.2f} pt text "
                             f"(< {MIN_FONT_PT} pt)")

    h_pt = sum(r["h"] for r in rows) + PT * sum(
        ROW_GAPS.get(r["letter"], GAP_IN) for r in rows[:-1])
    writer = PdfWriter()
    page = writer.add_blank_page(width=w_pt, height=h_pt)
    y = h_pt                                    # cursor from the top edge
    for r in rows:
        y -= r["h"]
        r["y_bottom_pt"] = y
        page.merge_transformed_page(r["page"],
                                    Transformation().translate(tx=r["x_pt"], ty=y))
        logger.info("placed %s at y=%.1f pt (h=%.1f pt, min font %.2f pt)",
                    r["letter"], y, r["h"], r["font_pt"])
        y -= PT * ROW_GAPS.get(r["letter"], GAP_IN)

    FIG.mkdir(parents=True, exist_ok=True)
    out_pdf = FIG / "Figure3_rev.pdf"
    with open(out_pdf, "wb") as fh:
        writer.write(fh)
    logger.info("wrote %s (%.2f x %.2f in)", out_pdf.name, w_pt / PT, h_pt / PT)

    out_png = FIG / "Figure3_rev.png"
    subprocess.run(["pdftoppm", "-png", "-r", "300", "-singlefile",
                    str(out_pdf), str(out_png.with_suffix(""))], check=True)
    logger.info("wrote %s (300 dpi)", out_png.name)

    # font embedding audit on the merged file
    fonts = subprocess.run(["pdffonts", str(out_pdf)], capture_output=True,
                           text=True).stdout.strip().splitlines()
    emb = []
    for line in fonts[2:]:
        parts = line.split()
        if len(parts) >= 5:
            emb.append((parts[0], parts[-5]))          # name, emb column
    bad = [n for n, e in emb if e.lower() == "no"]
    if bad:
        logger.error("fonts NOT embedded: %s", ", ".join(bad))
    else:
        logger.info("pdffonts: %d fonts, all embedded", len(emb))

    merged = PdfReader(str(out_pdf)).pages[0]
    merged_font = min_font_size(merged)
    logger.info("merged page: %d text objects, minimum font %.2f pt",
                len(re.findall(r"Tj", merged.get_contents().get_data()
                               .decode("latin-1", errors="ignore"))), merged_font)
    if merged_font < MIN_FONT_PT - 0.05:
        logger.error("merged page contains %.2f pt text (< %.1f pt)",
                     merged_font, MIN_FONT_PT)
    geo = [{"item": "page_in", "value": f"{w_pt / PT:.2f} x {h_pt / PT:.2f}"},
           {"item": "page_pt", "value": f"{w_pt:.1f} x {h_pt:.1f}"},
           {"item": "text_width_fraction", "value": f"{W_IN / 7.01:.2f}"}]
    for r in rows:
        geo.append({"item": f"row_{r['letter']}",
                    "value": f"{r['file']} ({r['backend']}), "
                             f"{r['w'] / PT:.2f} x {r['h'] / PT:.2f} in, "
                             f"y_bottom {r['y_bottom_pt'] / PT:.2f} in, "
                             f"min font {r['font_pt']:.2f} pt"})
    geo.append({"item": "min_font_pt",
                "value": f"{min(r['font_pt'] for r in rows):.2f} (per row); "
                         f"{merged_font:.2f} (merged page)"})
    geo.append({"item": "fonts_embedded", "value": f"{len(emb) - len(bad)}/{len(emb)}"})
    geo.append({"item": "vector", "value": "pypdf translation merge (no rasterisation)"})
    write_table(pd.DataFrame(geo), TAB / "fig3_assembly_geometry.tsv")
    logger.info("output: %s", OUT)


if __name__ == "__main__":
    main()
