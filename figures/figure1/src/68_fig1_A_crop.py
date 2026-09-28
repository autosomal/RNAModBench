#!/usr/bin/env python
"""68 -- cut panel A of Figure 1 out of the submitted PDF, as vectors (2026-09-21).

Why a real crop
---------------
The revised Figure 1 keeps everything of the published panel A (the dataset /
pre-processing / tools-comparison schematic was drawn by hand and has no
editable source in the repository) but redraws B, C and D at the printed size.
Painting a white box over the published B/C/D would leave their text *inside*
the file -- ``pdftotext`` would still return "RRACH Tool Counts Comparison" and
"Tool Frequency Comparison", and a vector editor would show the hidden page.
So the A strip is physically cut out of the page with ghostscript's ``CropBox``
pdfmark, which keeps the artwork vector (no rasterisation) and deletes every
object outside the box.

What this script measures
-------------------------
The published page is 595.276 x 532.615 pt (origin bottom-left).  Panel A is the
wide strip at the top, its three section headers sit at y ~ 473 and the lowest
piece of its artwork ("Time", the current trace) at y ~ 383, while the panel
letter "B" of the row below is at y = 345.  Everything is measured from a 300 dpi
render instead of being hard-coded:

* ``letter``  -- ink of the published "A" (x < 45 pt): the rebuild draws its own
  bold letters, so the old one must not survive inside the crop;
* ``artwork`` -- ink of the strip at x > 45 pt, i.e. the schematic itself;
* ``strip``   -- the artwork box grown by ``PAD_PT``, clipped to stay below the
  letter (when there is a clean white row between the two) and above the
  B/C/D row.

The crop is verified afterwards by re-rendering it and comparing it with the
same region of the original page pixel by pixel: a ghostscript rewrite that
moves the artwork by a point or loses a font would otherwise go unnoticed.

Outputs
-------
``figures/figure1/figures/panels/Figure1_A_from_published.pdf``
``figures/figure1/tables/fig1A_geometry.tsv``   (crop box, size, aspect)

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    figures/figure1/src/68_fig1_A_crop.py
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

import numpy as np
from PIL import Image, ImageFilter

HERE = _RB / "src/harmonisation"
sys.path.insert(0, str(HERE))

from common import config as C                        # noqa: E402
from common.io_utils import write_table               # noqa: E402
from common.manifest import setup_logger              # noqa: E402

OUT = (_RB / "figures/figure1")
FIG, TAB, LOG = (_RB / "figures/figure1/figures"), (_RB / "figures/figure1/tables"), (_RB / "figures/figure1/logs")
PANEL_DIR = FIG / "panels"
A_PDF = PANEL_DIR / "Figure1_A_from_published.pdf"
PUB = (_XB / "submission/manuscript/Figure1.pdf")

DPI = 300.0
#: page of the submitted figure, measured with pdfinfo (pt, origin bottom-left)
PAGE_W, PAGE_H = 595.276, 532.615
#: everything above this y is panel A; the B/C/D row starts at y ~ 345
BAND_TOP, BAND_BOTTOM = 525.0, 355.0
#: window that holds the published panel letter ("A" is at x 13.2..22.6,
#: y 499.7..508.9 per pdftotext, and it is the topmost ink of the page)
LETTER_WINDOW = (0.0, 40.0, 490.0, 520.0)
#: white space kept around the artwork
PAD_PT = 2.0
#: a pixel is "ink" below this grey level
INK = 250


def render(pdf: Path, png: Path, dpi: float = DPI) -> None:
    subprocess.run(["pdftoppm", "-png", "-singlefile", "-r", str(int(dpi)),
                    str(pdf), str(png.with_suffix(""))], check=True)
    assert png.exists(), f"pdftoppm wrote no file for {pdf}"


def grey(png: Path) -> np.ndarray:
    return np.asarray(Image.open(png).convert("L"), dtype=np.uint8)


def ink_box(img: np.ndarray, x0: float, x1: float, y0: float, y1: float,
            dpi: float = DPI) -> tuple[float, float, float, float] | None:
    """Ink bounding box of a page window, returned in PDF points.

    ``x0..y1`` are PDF coordinates (origin bottom-left) and the render is
    top-left, so the row window is flipped before it is sliced.
    """
    scale = dpi / 72.0
    px0, px1 = int(round(x0 * scale)), int(round(x1 * scale))
    r0 = int(round((PAGE_H - y1) * scale))
    r1 = int(round((PAGE_H - y0) * scale))
    win = img[max(r0, 0):r1, max(px0, 0):px1]
    if win.size == 0:
        return None
    mask = win < INK
    rows, cols = np.where(mask)
    if len(rows) == 0:
        return None
    return (px0 / scale + cols.min() / scale, PAGE_H - (r0 + rows.max()) / scale,
            px0 / scale + (cols.max() + 1) / scale,
            PAGE_H - (r0 + rows.min()) / scale)


def main() -> None:
    logger = setup_logger("68_fig1_A_crop", log_dir=LOG)
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    if not PUB.exists():
        raise SystemExit(f"published figure not found: {PUB}")

    raw = LOG / "fig1_published_page_300dpi.png"
    render(PUB, raw)
    img = grey(raw)
    logger.info("rendered %s at %d dpi -> %s", PUB.name, DPI, img.shape)

    lx0, lx1, ly0, ly1 = LETTER_WINDOW
    letter = ink_box(img, lx0, lx1, ly0, ly1)
    if letter is None:
        raise SystemExit("the published panel letter was not found -- the page "
                         "layout changed, re-measure before cropping")
    logger.info("published letter A: %s", _fmt(letter))

    #: the schematic is everything below the letter (a rectangle cannot exclude
    #: the letter vertically, because the artwork reaches further left than it)
    art = ink_box(img, 0.0, PAGE_W, BAND_BOTTOM, letter[1] - 1.0)
    if art is None:
        raise SystemExit("no artwork ink found in the panel-A band")
    logger.info("panel-A artwork   : %s", _fmt(art))

    x0 = max(0.0, art[0] - PAD_PT)
    x1 = min(PAGE_W, art[2] + PAD_PT)
    y0 = max(BAND_BOTTOM, art[1] - PAD_PT)
    y1 = min(letter[1] - 1.0, art[3] + PAD_PT)
    box = (x0, y0, x1, y1)
    gap = letter[1] - y1
    logger.info("crop box: x %.2f..%.2f, y %.2f..%.2f  (%.2f x %.2f pt); "
                "%.2f pt of white above it keep the published letter out",
                x0, x1, y0, y1, x1 - x0, y1 - y0, gap)
    if gap < 1.0:
        raise SystemExit(f"only {gap:.2f} pt between the artwork and the "
                         f"published letter: the rebuild would have to paint "
                         f"the letter over -- re-measure the page")
    stray = ink_box(img, 0.0, PAGE_W, y1, letter[1] - 0.5)
    if stray is not None and not (stray[0] >= letter[0] - 1.0
                                  and stray[2] <= letter[2] + 1.0):
        raise SystemExit(f"ink between the crop top and the letter: {_fmt(stray)}")

    _gs_crop(PUB, A_PDF, box)
    size = subprocess.run(["pdfinfo", str(A_PDF)], capture_output=True,
                          text=True, check=True).stdout
    page = [ln for ln in size.splitlines() if ln.startswith("Page size")][0]
    logger.info("gs output: %s (page origin moved to the crop corner)", page.strip())
    _check_text(A_PDF, logger)
    _reconcile(img, A_PDF, box, logger)

    write_table(_geometry_table(box, letter, art), TAB / "fig1A_geometry.tsv")
    logger.info("aspect (w/h) = %.4f; at 459.52 pt page width -> %.2f pt tall",
                (x1 - x0) / (y1 - y0), 459.52 / ((x1 - x0) / (y1 - y0)))
    logger.info("output: %s", A_PDF)


def _fmt(box: tuple[float, float, float, float] | None) -> str:
    if box is None:
        return "-"
    return (f"x {box[0]:.1f}..{box[2]:.1f}, y {box[1]:.1f}..{box[3]:.1f} "
            f"({box[2] - box[0]:.1f} x {box[3] - box[1]:.1f} pt)")


def _gs_crop(src: Path, dst: Path, box: tuple[float, float, float, float]) -> None:
    """Vector crop that really drops the rest of the page.

    A ``[/CropBox ...]`` pdfmark alone only *hides* the other panels: the media
    box stays 595 x 533 pt and ``pdftotext`` still returns the published B/C/D
    labels (measured 2026-09-21).  A fixed-size device plus a ``PageOffset``
    writes a page whose origin *is* the crop corner, and ghostscript clips every
    object outside it.
    """
    x0, y0, x1, y1 = box
    cmd = ["gs", "-q", "-o", str(dst), "-sDEVICE=pdfwrite",
           "-dPDFSETTINGS=/prepress", "-dEmbedAllFonts=true", "-dFIXEDMEDIA",
           f"-dDEVICEWIDTHPOINTS={x1 - x0:.3f}",
           f"-dDEVICEHEIGHTPOINTS={y1 - y0:.3f}",
           "-c", f"<</PageOffset [{-x0:.3f} {-y0:.3f}]>> setpagedevice",
           "-f", str(src)]
    subprocess.run(cmd, check=True)


def _check_text(pdf: Path, logger) -> None:
    """The dropped panels must be gone from the file, panel A must be intact."""
    text = subprocess.run(["pdftotext", "-q", str(pdf), "-"], capture_output=True,
                          text=True, check=True).stdout
    gone = ["Tool Frequency Comparison", "RRACH Tool Counts Comparison",
            "Curlcake_IVT", "Nanocompore"]
    kept = ["Datasets", "Pre-processing", "Tools comparison"]
    still = [s for s in gone if s in text]
    if still:
        raise SystemExit(f"the crop still contains the published B/C/D text: {still}")
    missing = [s for s in kept if s not in text]
    if missing:
        raise SystemExit(f"the crop lost panel-A text: {missing}")
    logger.info("text check: %d published B/C/D labels gone, all %d panel-A "
                "headers kept", len(gone), len(kept))


def _reconcile(orig: np.ndarray, crop_pdf: Path,
               box: tuple[float, float, float, float], logger) -> None:
    """Re-render the cut and compare it with the same window of the original.

    Two nuisance effects have to be separated from real loss:

    * ghostscript rounds the crop corner to the pixel grid, so the two renders
      are compared over the best of the +/- 2 px alignments;
    * a sub-pixel offset repaints every edge of the artwork, which shows up as a
      band of 30-60/255 differences along every 1 px line -- harmless, but it
      hides a dropped glyph.  Both renders are therefore blurred (sigma 1.2 px)
      before the comparison, and the gate asks for *content*: a lost letter, a
      re-routed arrow or a dropped icon moves whole regions by far more than
      64/255 and cannot be blurred away.
    """
    scale = DPI / 72.0
    x0, y0, x1, y1 = box
    c0, c1 = int(round(x0 * scale)), int(round(x1 * scale))
    r0, r1 = int(round((PAGE_H - y1) * scale)), int(round((PAGE_H - y0) * scale))
    ref = orig[max(r0 - 3, 0):r1 + 3, max(c0 - 3, 0):c1 + 3].astype(np.uint8)
    #: where the crop corner sits inside the padded reference window
    off_r, off_c = r0 - max(r0 - 3, 0), c0 - max(c0 - 3, 0)
    out_png = LOG / "fig1A_crop_300dpi.png"
    render(crop_pdf, out_png)
    got = grey(out_png)

    best = None
    for dy in range(-2, 3):
        for dx in range(-2, 3):
            sub = ref[off_r + dy:off_r + dy + got.shape[0],
                      off_c + dx:off_c + dx + got.shape[1]]
            g = got[:sub.shape[0], :sub.shape[1]]
            if sub.size == 0:
                continue
            d = np.abs(g.astype(np.int16) - sub.astype(np.int16))
            score = float((d > 64).mean())
            if best is None or score < best[0]:
                best = (score, sub, g, dx, dy)
    if best is None:
        raise SystemExit("the cropped page could not be compared with the original")
    _, sub, g, dx, dy = best
    blur = ImageFilter.GaussianBlur(1.2)
    rb = np.asarray(Image.fromarray(sub).filter(blur)).astype(np.int16)
    gb = np.asarray(Image.fromarray(g).filter(blur)).astype(np.int16)
    diff = np.abs(gb - rb)
    hard, soft = float((diff > 64).mean()), float((diff > 32).mean())
    logger.info("pixel reconciliation at the best alignment (%+d, %+d) px: "
                "blurred mean |diff| %.2f/255, %.4f%% of pixels move by more "
                "than 64/255, %.2f%% by more than 32/255",
                dx, dy, diff.mean(), 100 * hard, 100 * soft)
    if hard > 0.0005:
        raise SystemExit(f"cropped panel A lost content ({100 * hard:.3f}% of "
                         f"pixels move by more than 64/255) -- check the gs rewrite")
    logger.info("no panel-A content lost in the ghostscript rewrite")


def _geometry_table(box, letter, art):
    import pandas as pd

    x0, y0, x1, y1 = box
    rows = [
        {"key": "source", "value": str(PUB)},
        {"key": "source_page_pt", "value": f"{PAGE_W:.3f} x {PAGE_H:.3f}"},
        {"key": "crop_x0_pt", "value": f"{x0:.3f}"},
        {"key": "crop_y0_pt", "value": f"{y0:.3f}"},
        {"key": "crop_x1_pt", "value": f"{x1:.3f}"},
        {"key": "crop_y1_pt", "value": f"{y1:.3f}"},
        {"key": "crop_w_pt", "value": f"{x1 - x0:.3f}"},
        {"key": "crop_h_pt", "value": f"{y1 - y0:.3f}"},
        {"key": "aspect_w_over_h", "value": f"{(x1 - x0) / (y1 - y0):.6f}"},
        {"key": "artwork_box_pt", "value": f"{art[0]:.2f},{art[1]:.2f},"
                                           f"{art[2]:.2f},{art[3]:.2f}"},
        {"key": "letter_box_pt", "value": f"{letter[0]:.2f},{letter[1]:.2f},"
                                          f"{letter[2]:.2f},{letter[3]:.2f}"},
        {"key": "letter_kept_out_by_pt", "value": f"{letter[1] - y1:.2f}"},
        {"key": "gs_cmd", "value": "gs -sDEVICE=pdfwrite -dFIXEDMEDIA "
                                   "-dDEVICEWIDTHPOINTS=w -dDEVICEHEIGHTPOINTS=h "
                                   "-c '<</PageOffset [-x0 -y0]>> setpagedevice'"},
    ]
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
