#!/usr/bin/env python
"""69 -- compose the replaceable Figure 1 (published panel A + rebuilt B/C/D).

The submitted page is assembled from exactly two vector files:

* ``figures/panels/Figure1_A_redrawn.pdf`` -- the dataset schematic, redrawn at
  its printed size by 71 (the published crop is 3.3 pt at this width, so it can
  neither be read nor meet the 7 pt floor);
* ``figures/Figure1_rev_body.pdf`` -- panels B, C and D, drawn at the printed
  size by 67, with the A strip left blank.

pypdf places the A crop in that blank strip scaled to the measured aspect (it
must be the *only* scaling in the chain: the panels are already at their printed
size, so 1 + 1 means the fonts land on paper as written).  The composed page is
then re-rendered at 300 dpi and audited: page size, panel count, that the
published B/C/D text is gone, and that panel A still carries its own content.

Outputs
-------
``04_revision_analysis/fig1_revision/figures/Figure1_rev.{pdf,png}``
``04_revision_analysis/fig1_revision/tables/fig1_page_files.tsv`` (md5 + geometry)

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    01_code/code/sites_v2/scripts/69_fig1_page.py
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
import hashlib
import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd
from pypdf import PdfReader, PdfWriter, Transformation

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from common import config as C            # noqa: E402
from common.io_utils import write_table   # noqa: E402
from common.manifest import setup_logger  # noqa: E402

OUT = (_RB / "figures/figure1")
FIG, TAB, LOG = (_RB / "figures/figure1/figures"), (_RB / "figures/figure1/tables"), (_RB / "figures/figure1/logs")
BODY = FIG / "Figure1_rev_body.pdf"
#: 2026-09-23: the delivered panel A is the user's Illustrator RGB export of the
#: original artwork, turned into this strip by 72 (letter "A" deleted in the
#: content stream, trimmed to the ink, scaled and centred with the transform
#: baked into the page box).  Point this back at Figure1_A_redrawn.pdf to fall
#: back to the redraw, which is the variant that meets the 7 pt floor.
A_PDF = FIG / "panels" / "Figure1_A_from_ai_srgb.pdf"
PAGE = FIG / "Figure1_rev.pdf"
PNG = FIG / "Figure1_rev.png"


def load_panels_module():
    """The page geometry lives in 67 (its module name starts with a digit)."""
    path = Path(__file__).with_name("67_fig1_panels.py")
    spec = importlib.util.spec_from_file_location("fig1_panels", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["fig1_panels"] = mod
    spec.loader.exec_module(mod)
    return mod


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pdfinfo(path: Path, key: str = "Page size") -> str:
    out = subprocess.run(["pdfinfo", str(path)], capture_output=True, text=True,
                         check=True).stdout
    return next(ln.split(":", 1)[1].strip() for ln in out.splitlines()
                if ln.startswith(key))


def page_text(path: Path) -> str:
    return subprocess.run(["pdftotext", "-q", str(path), "-"], capture_output=True,
                          text=True, check=True).stdout


def main() -> None:
    logger = setup_logger("69_fig1_page", log_dir=LOG)
    mod67 = load_panels_module()
    for f in (BODY, A_PDF):
        if not f.exists():
            raise SystemExit(f"missing input: {f} (run 67 / 68 first)")

    a_page = PdfReader(str(A_PDF)).pages[0]
    a_w, a_h = float(a_page.mediabox.width), float(a_page.mediabox.height)
    # 71 draws the strip at exactly the size it is printed, so the only legal
    # placement is 1:1 -- a scale here would silently shrink its fonts again.
    if abs(a_w - mod67.A_RECT[2]) > 0.2 or abs(a_h - mod67.A_RECT[3]) > 0.2:
        raise SystemExit(f"A strip is {a_w:.2f} x {a_h:.2f} pt but the page "
                         f"reserves {mod67.A_RECT[2]:.2f} x {mod67.A_RECT[3]:.2f} pt")
    body_page = PdfReader(str(BODY)).pages[0]
    b_w, b_h = float(body_page.mediabox.width), float(body_page.mediabox.height)
    assert abs(b_w - mod67.PAGE_W) < 0.05 and abs(b_h - mod67.PAGE_H) < 0.05, \
        f"body page is {b_w:.2f} x {b_h:.2f} pt, expected {mod67.PAGE_W} x {mod67.PAGE_H}"

    scale = 1.0
    writer = PdfWriter()
    page = writer.add_blank_page(width=mod67.PAGE_W, height=mod67.PAGE_H)
    page.merge_transformed_page(
        body_page, Transformation().scale(1.0, 1.0).translate(0.0, 0.0))
    page.merge_transformed_page(
        a_page, Transformation().scale(scale, scale)
        .translate(mod67.A_RECT[0], mod67.A_RECT[1]))
    logger.info("A crop %.2f x %.2f pt scaled by %.4f -> %.2f x %.2f pt at "
                "(%.1f, %.1f)", a_w, a_h, scale, a_w * scale, a_h * scale,
                mod67.A_RECT[0], mod67.A_RECT[1])

    FIG.mkdir(parents=True, exist_ok=True)
    with open(PAGE, "wb") as fh:
        writer.write(fh)
    subprocess.run(["pdftoppm", "-png", "-r", "300", "-singlefile",
                    str(PAGE), str(PNG.with_suffix(""))], check=True)

    size = pdfinfo(PAGE)
    logger.info("composed page: %s (%d page(s))", size, len(PdfReader(str(PAGE)).pages))
    # compare numerically: pdfinfo prints as many decimals as the value needs
    # (506.4591 -> "506.459"), so a substring match on "%.2f" is not a check
    got_w, got_h = (float(part.split()[0]) for part in size.split(" x ")[:2])
    if abs(got_w - mod67.PAGE_W) > 0.02 or abs(got_h - mod67.PAGE_H) > 0.02:
        raise SystemExit(f"page size drifted: {size}, expected "
                         f"{mod67.PAGE_W:.4f} x {mod67.PAGE_H:.4f} pt")

    text = page_text(PAGE)
    gone = ["Tool Frequency Comparison", "RRACH Tool Counts Comparison"]
    kept = ["Datasets", "Pre-processing", "Tools comparison", "Nanocompore",
            "m6A (mean)", "IVT (mean)"]
    still = [s for s in gone if s in text]
    missing = [s for s in kept if s not in text]
    if still or missing:
        raise SystemExit(f"composed page text is wrong: still {still}, missing {missing}")
    logger.info("text check: published B/C/D headers gone, panel A and the "
                "rebuilt panels present")

    write_table(pd.DataFrame([
        {"key": "figure", "value": str(PAGE)},
        {"key": "page_pt", "value": size},
        {"key": "A_source", "value": str(A_PDF)},
        {"key": "A_scale", "value": f"{scale:.6f}"},
        {"key": "body_source", "value": str(BODY)},
        {"key": "md5_final_pdf", "value": md5(PAGE)},
        {"key": "md5_final_png", "value": md5(PNG)},
        {"key": "md5_body_pdf", "value": md5(BODY)},
        {"key": "md5_A_pdf", "value": md5(A_PDF)},
    ]), TAB / "fig1_page_files.tsv")
    logger.info("md5 final = %s", md5(PAGE))
    logger.info("wrote %s (+ .png)", PAGE)


if __name__ == "__main__":
    main()
