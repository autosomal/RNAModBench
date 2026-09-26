#!/usr/bin/env python
"""72 -- panel A from the Illustrator RGB artwork (user supplied 2026-09-23).

Why this script exists
----------------------
67/69 reserve a 486.4591 x 121.0 pt strip for panel A and place it 1:1 (the
only scaling in the chain happens at drawing time).  The redraw (71) satisfied
the 7 pt floor but not the published look, so the user re-exported the original
panel from Illustrator with RGB process colours (``03_figures/Fig1a.pdf``,
read-only).  This script turns that artwork into the strip, in vectors only:

* the panel letter "A" is a standalone text run (x 13.3-22.7, y 6.8-24.9 from
  the page top) and is **deleted from the content stream** -- cropping it away
  is not an option, the column headers start further right at the same height;
* the strip is then trimmed to the remaining ink (measured on a 300 dpi render
  in memory, +1.5 pt of margin) so the page's empty margins do not travel;
* the trimmed box is scaled by ``min(W/ink_w, H/ink_h)`` (~0.888) and centred,
  and that transform is **baked into the strip's page box** (486.4591 x 121.0
  pt) so 69 keeps its 1:1 invariant and 70's page-box assertion still holds.

Accepted consequences (user, 2026-09-23)
----------------------------------------
* panel A is exempt from the 7 pt floor: the artwork's own sizes print as-is
  (body labels ~6.6 pt, small print ~2.5-4.5 pt) -- the printed table below is
  the record of that exception;
* two of the three fonts carry no ToUnicode, so a few words extract fragmented
  ("la", "He", "li co"): that is the artwork's own encoding, kept on purpose;
* the six DIC spot inks stay ``/Separation`` with Lab alternates, so panel A's
  palette stays its own (measured distance to the B/C/D constants: 1-44).

Outputs
-------
``04_revision_analysis/fig1_revision/figures/panels/Figure1_A_from_ai_srgb.pdf``
``04_revision_analysis/fig1_revision/figures/panels/Figure1_A_from_ai_srgb.png``
``04_revision_analysis/fig1_revision/tables/fig1a_ai_printed_type.tsv``
``04_revision_analysis/_figfix_20260923/_audit_tmp/fig1a_noletter.pdf`` (scratch)

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    01_code/code/sites_v2/scripts/72_fig1_A_from_ai.py
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
import io
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import DecodedStreamObject, NameObject

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from common import config as C               # noqa: E402
from common.io_utils import write_table      # noqa: E402
from common.manifest import setup_logger     # noqa: E402

OUT = (_RB / "figures/figure1")
FIG, TAB, LOG = (_RB / "figures/figure1/figures"), (_RB / "figures/figure1/tables"), (_RB / "figures/figure1/logs")
PANEL = FIG / "panels"
STRIP = PANEL / "Figure1_A_from_ai_srgb.pdf"
PREVIEW = PANEL / "Figure1_A_from_ai_srgb.png"
SRC = (_XB / "figures_original/Fig1a.pdf")          # user artwork, read-only
SCRATCH = ((_RB / "analysis/_figfix_20260923/_audit_tmp/fig1a_noletter.pdf"))

#: the strip 67 reserves -- the single source of truth is its A_RECT, so read it
def strip_size() -> tuple[float, float]:
    import importlib.util
    path = Path(__file__).with_name("67_fig1_panels.py")
    spec = importlib.util.spec_from_file_location("fig1_panels_72", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["fig1_panels_72"] = mod
    spec.loader.exec_module(mod)
    return float(mod.A_RECT[2]), float(mod.A_RECT[3])


MARGIN_PT = 1.5          # ink margin kept around the trimmed artwork
DPI = 300                # render used for the ink measurement
#: pdftotext word boxes are font metric boxes (ascent+descent), not em sizes;
#: 1.117 is the Arial ratio measured on the redrawn strip (70 derives it too)
BOX_PER_PT = 1.117
LETTER_A_MAX_X = 30.0    # PDF pt from the left edge: the letter lives below it
LETTER_A_MIN_Y = 120.0   # PDF pt from the bottom: and above it

BT_ET = re.compile(rb"BT(.*?)ET", re.S)
TM = re.compile(rb"([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+"
                rb"([-\d.]+)\s+([-\d.]+)\s+Tm")
SHOW = re.compile(rb"\((?:[^()\\]|\\.)*\)\s*Tj|\[(?:[^\[\]\\]|\\.)*\]\s*TJ")
LITERAL = re.compile(rb"\((?:[^()\\]|\\.)*\)")
WORD = re.compile(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" '
                  r'yMax="([\d.]+)">(.*?)</word>')


def literal_text(raw: bytes) -> str:
    """Characters of the first literal string of a show operator."""
    m = LITERAL.search(raw)
    if not m:
        return ""
    body = m.group(0)[1:-1]
    body = re.sub(rb"\\([()\\])", rb"\1", body)
    return body.decode("latin-1")


def drop_letter_a(page, logger) -> None:
    """Delete the standalone panel letter "A" from a copied page's content."""
    data = page.get_contents().get_data()
    keep, dropped = [], []
    pos = 0
    for m in BT_ET.finditer(data):
        block = m.group(0)
        show = SHOW.search(block)
        tm = TM.search(block)
        if show and tm:
            text = literal_text(show.group(0))
            x, y = float(tm.group(5)), float(tm.group(6))
            if (text.strip() == "A" and x < LETTER_A_MAX_X
                    and y > LETTER_A_MIN_Y):
                dropped.append((text, x, y))
                keep.append(data[pos:m.start()])
                keep.append(b" ")           # a space keeps the operators apart
                pos = m.end()
                continue
        logger.info("text run kept: %r at (%.1f, %.1f)",
                    literal_text(show.group(0)) if show else "", *(
                        (float(tm.group(5)), float(tm.group(6))) if tm
                        else (float("nan"), float("nan"))))
    keep.append(data[pos:])
    if len(dropped) != 1:
        raise SystemExit(f"expected exactly one letter-A run, found {dropped} "
                         f"-- adjust LETTER_A_MAX_X / LETTER_A_MIN_Y")
    logger.info("dropped the panel letter: %r at (%.1f, %.1f)", *dropped[0])
    stream = DecodedStreamObject()
    stream.set_data(b"".join(keep))
    add = getattr(page, "replace_contents", None)
    if callable(add):                       # pypdf >= 3.5
        page.replace_contents(stream)
    else:                                   # pragma: no cover
        raise SystemExit("pypdf too old: no PageObject.replace_contents")


def ink_box(pdf: Path, logger) -> tuple[float, float, float, float]:
    """Ink bounding box in PDF points (origin bottom-left), margin included."""
    png = subprocess.run(["/usr/bin/pdftocairo", "-png", "-r", str(DPI),
                          "-singlefile", str(pdf), "-"],
                         capture_output=True, check=True).stdout
    import matplotlib.image as mpimg
    img = mpimg.imread(io.BytesIO(png), format="png")
    ink = (img[:, :, :3].min(axis=2) < 0.98)
    rows, cols = np.where(ink)
    if not len(rows):
        raise SystemExit(f"{pdf} renders blank")
    box = PdfReader(str(pdf)).pages[0].mediabox
    w_pt, h_pt = float(box.width), float(box.height)
    scale = 72.0 / DPI
    x0, x1 = cols.min() * scale, (cols.max() + 1) * scale
    top, bot = rows.min() * scale, (rows.max() + 1) * scale
    logger.info("rendered ink %d x %d px = %.2f x %.2f pt at %d dpi",
                cols.max() - cols.min() + 1, rows.max() - rows.min() + 1,
                x1 - x0, bot - top, DPI)
    return (x0, h_pt - bot, x1, h_pt - top)


def main() -> None:
    logger = setup_logger("72_fig1_A_from_ai", log_dir=LOG)
    strip_w, strip_h = strip_size()
    logger.info("target strip %.4f x %.2f pt (from 67.A_RECT)", strip_w, strip_h)
    if not SRC.exists():
        raise SystemExit(f"missing user artwork: {SRC}")

    # ---- 1) letter A out, in vectors -------------------------------------- #
    writer = PdfWriter(clone_from=str(SRC))
    drop_letter_a(writer.pages[0], logger)
    PANEL.mkdir(parents=True, exist_ok=True)
    SCRATCH.parent.mkdir(parents=True, exist_ok=True)
    with open(SCRATCH, "wb") as fh:
        writer.write(fh)

    # ---- 2) trim to the ink, scale, centre -------------------------------- #
    x0, y0, x1, y1 = ink_box(SCRATCH, logger)
    x0, y0 = x0 - MARGIN_PT, y0 - MARGIN_PT
    x1, y1 = x1 + MARGIN_PT, y1 + MARGIN_PT
    ink_w, ink_h = x1 - x0, y1 - y0
    s = min(strip_w / ink_w, strip_h / ink_h)
    tx = (strip_w - ink_w * s) / 2.0 - x0 * s
    ty = (strip_h - ink_h * s) / 2.0 - y0 * s
    logger.info("crop (%.2f, %.2f)-(%.2f, %.2f) = %.2f x %.2f pt -> scale %.4f "
                "-> %.2f x %.2f pt, offset (%.2f, %.2f)",
                x0, y0, x1, y1, ink_w, ink_h, s, ink_w * s, ink_h * s, tx, ty)

    strip = PdfWriter()
    page = strip.add_blank_page(width=strip_w, height=strip_h)
    page.merge_transformed_page(writer.pages[0],
                               Transformation().scale(s, s).translate(tx, ty))
    with open(STRIP, "wb") as fh:
        strip.write(fh)

    # ---- 3) gates --------------------------------------------------------- #
    box = PdfReader(str(STRIP)).pages[0].mediabox
    got_w, got_h = float(box.width), float(box.height)
    if abs(got_w - strip_w) > 0.05 or abs(got_h - strip_h) > 0.05:
        raise SystemExit(f"strip page is {got_w:.2f} x {got_h:.2f} pt, "
                         f"expected {strip_w:.2f} x {strip_h:.2f}")
    listed = subprocess.run(["pdfimages", "-list", str(STRIP)],
                            capture_output=True, text=True).stdout
    images = [ln for ln in listed.splitlines()[2:] if ln.strip()]
    if images:
        raise SystemExit(f"the strip is not fully vector: {len(images)} images")
    fonts = subprocess.run(["pdffonts", str(STRIP)], capture_output=True,
                           text=True).stdout
    logger.info("strip fonts (emb/uni are the artwork's own):\n%s", fonts.strip())

    xml = subprocess.run(["/usr/bin/pdftotext", "-bbox", str(STRIP), "-"],
                         capture_output=True, text=True, check=True).stdout
    words = [(float(a), float(b), float(c), float(d), t)
             for a, b, c, d, t in WORD.findall(xml)]
    lone = [t for *_, t in words if t.strip() == "A"]
    if lone:
        raise SystemExit(f"the panel letter is still on the strip: {lone}")
    for must in ("Datasets", "Pre-processing", "Tools"):
        if not any(must in t for *_, t in words):
            raise SystemExit(f"the strip lost its header {must!r}")
    sizes = sorted(round((d - b) / BOX_PER_PT, 2) for _, b, _, d, _ in words)
    logger.info("printed type of panel A: min %.2f pt | median %.2f pt | "
                "max %.2f pt (%d words) -- the artwork's own sizes, panel A is "
                "exempt from the 7 pt floor by the user's decision",
                sizes[0], float(np.median(sizes)), sizes[-1], len(sizes))

    subprocess.run(["/usr/bin/pdftoppm", "-png", "-r", "300", "-singlefile",
                    str(STRIP), str(PREVIEW.with_suffix(""))], check=True)
    rows = [{"key": "source_pdf", "value": str(SRC)},
            {"key": "source_page_pt", "value": "595.276 x 158.261"},
            {"key": "letter_a_dropped", "value": "yes (vector, content stream)"},
            {"key": "crop_pt", "value": f"{x0:.2f},{y0:.2f},{x1:.2f},{y1:.2f}"},
            {"key": "scale", "value": f"{s:.6f}"},
            {"key": "printed_pt_min", "value": f"{sizes[0]:.2f}"},
            {"key": "printed_pt_median", "value": f"{float(np.median(sizes)):.2f}"},
            {"key": "printed_pt_max", "value": f"{sizes[-1]:.2f}"},
            {"key": "strip_pt", "value": f"{got_w:.4f} x {got_h:.2f}"},
            {"key": "images_in_strip", "value": "0 (fully vector)"},
            {"key": "note", "value": "panel A: artwork look and own type sizes; "
                                     "DIC spot inks kept as Separation/Lab"}]
    write_table(pd.DataFrame(rows), TAB / "fig1a_ai_printed_type.tsv")
    write_table(pd.DataFrame([{"word": t, "box_h_pt": round(d - b, 2),
                              "printed_pt": round((d - b) / BOX_PER_PT, 2)}
                             for _, b, _, d, t in
                             sorted(words, key=lambda v: v[3] - v[1])]),
                TAB / "fig1a_ai_words.tsv")
    logger.info("wrote %s (+ .png preview) and %s", STRIP,
                TAB / "fig1a_ai_printed_type.tsv")


if __name__ == "__main__":
    main()
