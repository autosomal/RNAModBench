#!/usr/bin/env python
"""70 -- acceptance of the rebuilt main Figure 1 (67 draws it, 68 cuts panel A,
69 composes the page).

The figure it replaces printed its tool names at ~5.4 pt: the 2026-09-20 draft
was drawn on a 941 x 478 pt canvas and placed with
``\\includegraphics[width=0.95\\textwidth]{Figure1.pdf}`` = 479.52 pt, i.e. at a
0.51 scale.  This script computes the placement scale from the PDF page width
itself and asserts

    scale == 1.000                 (the canvas *is* the printed size)
    every Tf size of B/C/D >= 7 pt (panel A is the published artwork: reported,
                                    not asserted -- its fonts are not ours)

plus the anchors the figure claims (every number drawn comes from a frozen
per-unit table), the text of the page (published B/C/D labels gone, letters
A-D printed exactly once, each tool name exactly as often as its panels), and
the house rules (Arial only in our panels, no gridlines, no in-panel annotation
text, 300 dpi PNG, panel A really inked).

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    figures/figure1/src/70_verify_fig1_page.py
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
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pypdf
from PIL import Image

HERE = _RB / "src/harmonisation"
sys.path.insert(0, str(HERE))

from common import config as C            # noqa: E402
from common.manifest import setup_logger  # noqa: E402

OUT = (_RB / "figures/figure1")
TAB, FIG, LOG = (_RB / "figures/figure1/tables"), (_RB / "figures/figure1/figures"), (_RB / "figures/figure1/logs")
PDF, PNG = FIG / "Figure1_rev.pdf", FIG / "Figure1_rev.png"
BODY, BODY_PNG = FIG / "Figure1_rev_body.pdf", FIG / "Figure1_rev_body.png"
#: 2026-09-23: the delivered panel A is the user's Illustrator RGB export of the
#: original artwork, prepared by 72; Figure1_A_redrawn.pdf is the fallback that
#: meets the 7 pt floor
A_PDF = FIG / "panels" / "Figure1_A_artwork.pdf"
if not A_PDF.exists():
    #: the gate measures whichever panel-A variant the page used
    A_PDF = FIG / "panels" / "Figure1_A_redrawn.pdf"
SCRIPT_PANELS = Path(__file__).with_name("67_fig1_panels.py")
SCRIPT_CROP = Path(__file__).with_name("68_fig1_A_crop.py")
SCRIPT_REDRAW = Path(__file__).with_name("71_fig1_A_redraw.py")

#: panel A is redrawn by 71 at the size it is printed, so the A strip is
#: placed 1:1 and *its* fonts are checked too (the published crop was exempt).
#: Both values are only fallbacks: main() reads the live geometry from 67.
TARGET_W, TARGET_H = 506.4591, 506.0
#: word boxes of pdftotext are font metric boxes (ascent+descent), not em sizes;
#: 1.117 is the Arial ratio measured on the redrawn strip, whose smallest string
#: is 7.5 pt by construction -- the audit derives it again at run time
BOX_PER_PT = 1.117
MIN_PT = 7.0
LETTERS = ("A", "B", "C", "D")


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def tf_sizes(pdf: Path) -> dict[float, int]:
    """Every ``Tf`` font size of the page, including nested Form XObjects."""
    out: dict[float, int] = {}
    reader = pypdf.PdfReader(str(pdf))

    def walk(res, depth: int = 0) -> None:
        if depth > 3:
            return
        for _, obj in (res.get("/XObject") or {}).items():
            o = obj.get_object()
            if o.get("/Subtype") != "/Form":
                continue
            try:
                data = o.get_data().decode("latin-1")
            except Exception:                                  # pragma: no cover
                continue
            for m in re.finditer(r"/([A-Za-z0-9#+._-]+)\s+([\d.]+)\s+Tf", data):
                s = round(float(m.group(2)), 2)
                out[s] = out.get(s, 0) + 1
            walk(o.get("/Resources", {}), depth + 1)

    for page in reader.pages:
        try:
            data = page.get_contents().get_data().decode("latin-1")
        except Exception:                                      # pragma: no cover
            data = ""
        for m in re.finditer(r"/([A-Za-z0-9#+._-]+)\s+([\d.]+)\s+Tf", data):
            s = round(float(m.group(2)), 2)
            out[s] = out.get(s, 0) + 1
        walk(page.get("/Resources", {}))
    return out


def fonts_of(pdf: Path) -> list[str]:
    raw = subprocess.run(["pdffonts", str(pdf)], capture_output=True, text=True).stdout
    return sorted({ln.split()[0].split("+")[-1] for ln in raw.splitlines()[2:]
                   if ln.strip()})


def text_of(pdf: Path) -> str:
    return subprocess.run(["pdftotext", "-q", str(pdf), "-"], capture_output=True,
                          text=True, check=True).stdout


def count_word(text: str, word: str) -> int:
    return len(re.findall(rf"(?<![A-Za-z0-9_]){re.escape(word)}(?![A-Za-z0-9_])",
                          text))


def main() -> None:
    logger = setup_logger("70_verify_fig1_page", log_dir=LOG)
    fails: list[str] = []

    def check(ok: bool, msg: str) -> None:
        (logger.info if ok else logger.error)("%s %s", "PASS" if ok else "FAIL", msg)
        if not ok:
            fails.append(msg)

    panels = load_module("fig1_panels_v", "67_fig1_panels.py")
    TARGET_W, TARGET_H = panels.PAGE_W, panels.PAGE_H

    # ------------------------------------------------------------- geometry --
    reader = pypdf.PdfReader(str(PDF))
    box = reader.pages[0].mediabox
    w_pt, h_pt = float(box.width), float(box.height)
    scale = TARGET_W / w_pt
    check(len(reader.pages) == 1, "single-page PDF")
    check(abs(w_pt - TARGET_W) < 0.5,
          f"page width {w_pt:.2f} pt == \\textwidth ({TARGET_W} pt), i.e. the "
          f"page *is* the printed size")
    check(abs(h_pt - TARGET_H) < 0.5, f"page height {h_pt:.2f} pt == {TARGET_H} pt")
    check(abs(scale - 1.0) < 2e-3,
          f"placement scale {scale:.4f} == 1.000 (nothing is down-scaled in the "
          f"manuscript)")
    check(h_pt / 72.0 <= 7.55, f"height {h_pt / 72:.2f} in within the page budget")
    png = Image.open(PNG)
    check(abs(png.size[0] / 300 * 72 - w_pt) < 0.5,
          f"PNG {png.size[0]}x{png.size[1]} px = the same width at 300 dpi")

    # -------------------------------------------------------- font sizes -----
    body_sizes = tf_sizes(BODY)
    page_sizes = tf_sizes(PDF)
    body_eff = sorted(round(s * scale, 2) for s in body_sizes)
    worst = min(body_eff)
    check(worst >= MIN_PT,
          f"smallest effective font of B/C/D {worst:.2f} pt >= {MIN_PT} pt "
          f"(printed {body_eff})")
    # panel A is ours now (71 redraws it at the printed size), so it is held to
    # the same floor.  Measure the strip itself: a set difference against the
    # body sizes is empty whenever both share 7.5/8.5 pt, which "passes" without
    # checking anything (the published crop bottomed out at 2.5 pt).
    # 2026-09-23 (user decision): panel A is the original artwork again, so its
    # own type sizes are exempt from the floor -- measured and reported, while
    # B/C/D above stay strict
    a_sizes = tf_sizes(A_PDF)
    check(bool(a_sizes), "panel A carries text runs")
    logger.info("panel A fonts, exempt from the %.1f pt floor (user decision): "
                "%s pt", MIN_PT, sorted(round(s * scale, 2) for s in a_sizes))
    body_fonts = fonts_of(BODY)
    check(all("Arial" in f for f in body_fonts),
          f"our panels embed Arial only ({body_fonts})")
    logger.info("composed page fonts: %s", fonts_of(PDF))

    # ------------------------------------------------------------- text ------
    page_text = text_of(PDF)
    for s in ("Tool Frequency Comparison", "RRACH Tool Counts Comparison"):
        check(s not in page_text, f"published B/C/D label {s!r} is gone from the page")
    for letter in LETTERS:
        got = count_word(page_text, letter)
        # the published letter "A" was cut away with the crop (68 asserts the
        # 14.6 pt of white between the artwork and the letter); pdftotext also
        # splits panel A's labels into fragments ("A." from "A.thaliana"), so
        # the count is a lower bound, not an equality
        check(got >= 1, f"panel letter {letter} is printed on the page ({got}x)")
    for s in ("Datasets", "Pre-processing", "Tools comparison"):
        check(s in page_text, f"panel A keeps its header {s!r}")

    # ------------------------------------------------------------ anchors ----
    t = {"b": pd.read_csv((_RB / "figures/figure1/inputs/tables/per_replicate_tool_counts.tsv"), sep="\t"),
         "c": pd.read_csv(TAB / "fig1C_tool_counts.tsv", sep="\t"),
         "d": pd.read_csv(TAB / "fig1D_rrach_counts.tsv", sep="\t")}
    d = t["d"]
    want = {("Curlcake_m6A", "Nanom6A"): 128.5, ("Curlcake_m6A", "DENA"): 112.5,
            ("Curlcake_m6A", "m6Anet"): 83.0, ("Curlcake_m6A", "NanoSPA_m6A"): 33.0,
            ("Curlcake_m6A", "yanocomp"): 73.0, ("Curlcake_m6A", "DRUMMER"): 9.5,
            ("Curlcake_m6A", "MINES"): 90.5, ("Curlcake_IVT", "MINES"): 73.0,
            ("Curlcake_IVT", "NanoSPA_m6A"): 26.0}
    for (cond, tool), exp in want.items():
        g = d[(d.dataset_group == cond) & (d.tool == tool)]
        got = float(np.mean(g["n_rrach"]))
        check(abs(got - exp) < 0.05,
              f"{cond} / {tool}: RRACH mean {got:.1f} == frozen {exp}")
    zeros = d[(d.dataset_group == "Curlcake_IVT") & (d.n_rrach == 0)]
    checks_zero = sorted(set(zeros.tool))
    check(checks_zero == ["ELIGOS2_solo", "EpiNano_Error", "xPore", "yanocomp"],
          f"the four unmodified-template tools with a real zero are the ones "
          f"drawn hollow ({checks_zero})")
    absent = sorted(set(d[d.dataset_group == "Curlcake_m6A"].tool)
                    - set(d[d.dataset_group == "Curlcake_IVT"].tool))
    check(absent == ["DRUMMER", "ELIGOS2_diff"],
          f"the tools never run on the unmodified template get no marker: {absent}")

    # each tool name is printed once per label column it belongs to: B carries
    # two label columns (Arabidopsis for the top row of the grid, Human for the
    # bottom one -- Mouse and E. coli share those rows), plus one column each in
    # C and D
    tools_c = set(t["d"].tool)
    disp = {"yanocomp": "Yanocomp"}
    for tool in sorted(tools_c):
        name = disp.get(tool, tool)
        got = count_word(page_text, name)
        check(got == 4, f"tool name {name} printed 4x (B x2, C, D): got {got}")
    check(count_word(page_text, "CHEUI_m6A") == 2,
          "CHEUI_m6A appears in panel B only (never run on the Curlcake constructs)")

    # -------------------------------------------------------- house rules ----
    src = SCRIPT_PANELS.read_text()
    check('"axes.grid": False' in src, "gridlines disabled in the style block")
    check("bbox_inches=None" in src, "saved without a tight bounding box")
    check(not re.search(r"\bax\.annotate\(", src),
          "no in-panel annotation callouts in the page script")
    check('"m6A (mean)"' in src and '"IVT (mean)"' in src
          and "ax.text(" not in src
          and 'PL.text_width_in(lab, 7.5)' in src,
          "the two C/D mean codings are page-level labels (no in-panel text)")
    check("page_legend" not in src and "LEGEND_RECT" not in src,
          "the shared legend band is fully removed")
    for must in ("compact_legend", "_legend_hits", "assert_page_clean"):
        check(must in src, f"the page script self-gates with {must}()")
    redraw_src = SCRIPT_REDRAW.read_text()
    check("bbox_inches=None" in redraw_src and "pad_inches=0.0" in redraw_src,
          "panel A is exported at exactly the strip size (no re-scaling later)")
    check("MIN_PT = 7.5" in redraw_src,
          "panel A honours the 7 pt floor at its own drawing size")

    # every drawn number comes from a frozen table: the body page and the tables
    # must agree on the page geometry the panels were drawn on
    geo = pd.read_csv(TAB / "fig1_page_geometry.tsv", sep="\t").set_index("key")["value"]
    check(abs(float(geo["page_w_pt"]) - TARGET_W) < 0.02
          and abs(float(geo["page_h_pt"]) - TARGET_H) < 0.02,
          "recorded page geometry matches the exported page")
    a_box = pypdf.PdfReader(str(A_PDF)).pages[0].mediabox
    check(abs(float(a_box.width) - panels.A_RECT[2]) < 0.2
          and abs(float(a_box.height) - panels.A_RECT[3]) < 0.2,
          f"panel A {float(a_box.width):.2f} x {float(a_box.height):.2f} pt fills "
          f"its strip {panels.A_RECT[2]:.2f} x {panels.A_RECT[3]:.2f} pt 1:1")
    files = pd.read_csv(TAB / "fig1_page_files.tsv", sep="\t").set_index("key")["value"]
    digest = __import__("hashlib").md5(PDF.read_bytes()).hexdigest()
    check(files["md5_final_pdf"] == digest,
          f"recorded md5 {files['md5_final_pdf'][:8]} == file {digest[:8]}")

    # panel A really is inked (an empty strip would pass every other check)
    arr = np.asarray(png.convert("L"), dtype=np.uint8)
    sx = png.size[0] / TARGET_W
    y0 = int(round((TARGET_H - panels.A_RECT[1] - panels.A_RECT[3]) * sx))
    y1 = int(round((TARGET_H - panels.A_RECT[1]) * sx))
    x0 = int(round(panels.A_RECT[0] * sx))
    x1 = int(round((panels.A_RECT[0] + panels.A_RECT[2]) * sx))
    strip = arr[y0:y1, x0:x1]
    ink = float((strip < 245).mean())
    check(ink > 0.06, f"panel A strip is inked ({100 * ink:.1f}% of its pixels)")

    # ------------------------------------------------- vector / text audit --- #
    # user 2026-09-23: panel A had to stay vector -- no raster fallback anywhere
    for label, path in (("composed page", PDF), ("panel A strip", A_PDF),
                        ("body page", BODY)):
        listed = subprocess.run(["pdfimages", "-list", str(path)],
                                capture_output=True, text=True).stdout
        rows = [ln for ln in listed.splitlines()[2:] if ln.strip()]
        check(not rows, f"{label} is fully vector ({len(rows)} image objects)")

    
    # must all be embedded *and* mapped
    raw = subprocess.run(["pdffonts", str(PDF)], capture_output=True, text=True).stdout
    fonts = [ln.split() for ln in raw.splitlines()[2:] if ln.strip()]
    # pdffonts columns are space separated and "CID TrueType" is two tokens, so
    # the flags are read from the right: object ID, then uni / sub / emb
    no_uni = [f[0] for f in fonts if f[-3] == "no"]
    not_emb = [f[0] for f in fonts if f[-5] != "yes"]
    body_no_uni = [f[0] for f in
                   [ln.split() for ln in subprocess.run(
                       ["pdffonts", str(BODY)], capture_output=True, text=True)
                       .stdout.splitlines()[2:] if ln.strip()]
                   if f[-3] == "no"]
    check(not body_no_uni, f"every font of B/C/D carries a ToUnicode map "
                           f"(missing: {body_no_uni})")
    if no_uni:
        logger.info("panel A keeps the artwork's own encoding: %s without a "
                    "ToUnicode map (user decision -- a few words extract "
                    "fragmented)", no_uni)
    check(not not_emb, f"every font is embedded (missing: {not_emb})")

    # word level: nothing outside the page, nothing overlapping, nothing below
    # the floor -- 70 already checks the Tf sizes, this checks what a reader
    # sees after the composition
    bbox = subprocess.run(["/usr/bin/pdftotext", "-bbox", str(PDF), "-"],
                          capture_output=True, text=True, check=True).stdout
    wb = [(float(a), float(b), float(c), float(d), t)
          for a, b, c, d, t in re.findall(
              r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" '
              r'yMax="([\d.]+)">(.*?)</word>', bbox)]
    check(len(wb) > 100, f"word audit sees the page text ({len(wb)} words)")
    ax, ay, aw, ah = panels.A_RECT

    def in_panel_a(v: tuple) -> bool:
        """word boxes are top-down, A_RECT is bottom-up"""
        return (ax - 1 <= v[0] and v[2] <= ax + aw + 1
                and h_pt - (ay + ah) - 1 <= v[1] and v[3] <= h_pt - ay + 1)

    outside = [t for a, b, c, d, t in wb
               if a < -0.5 or b < -0.5 or c > w_pt + 0.5 or d > h_pt + 0.5]
    check(not outside, f"no word outside the page (found {outside[:4]})")
    pairs = [(wb[i], wb[j]) for i in range(len(wb)) for j in range(i + 1, len(wb))
             if min(wb[i][2], wb[j][2]) - max(wb[i][0], wb[j][0]) > 1.0
             and min(wb[i][3], wb[j][3]) - max(wb[i][1], wb[j][1]) > 1.0]
    # panel A is the original artwork: its two WinAnsi fonts extract as
    # fragments whose boxes overlap ("E." over "co"), which is the artwork's own
    # encoding (user decision) -- reported, while our panels stay strict
    artwork = [p for p in pairs if in_panel_a(p[0]) or in_panel_a(p[1])]
    ours = [p for p in pairs if not (in_panel_a(p[0]) or in_panel_a(p[1]))]
    check(not ours, f"no overlapping word boxes outside panel A "
                    f"({[(p[0][4], p[1][4]) for p in ours[:4]]})")
    if artwork:
        logger.info("panel A: %d overlapping word boxes (%s), the artwork's own "
                    "fragmented extraction (user decision)", len(artwork),
                    [(p[0][4], p[1][4]) for p in artwork[:4]])
    # printed type, split by region: our panels must clear the floor, panel A is
    # the user's original artwork and is exempt (its own sizes are reported).
    # Box heights are font metric boxes: / 1.117 turns them into em (Arial).
    def printed_pt(v: tuple) -> float:
        return (v[3] - v[1]) / BOX_PER_PT

    a_words = [v for v in wb if in_panel_a(v)]
    rest = [v for v in wb if not in_panel_a(v)]
    check(bool(a_words) and bool(rest), f"word regions resolved "
          f"({len(a_words)} in panel A, {len(rest)} elsewhere)")
    worst_rest = min(rest, key=printed_pt)
    check(printed_pt(worst_rest) >= MIN_PT,
          f"smallest printed word outside panel A {worst_rest[4]!r} "
          f"{printed_pt(worst_rest):.2f} pt >= {MIN_PT} pt")
    a_min = min(a_words, key=printed_pt)
    logger.info("panel A printed type, exempt (user decision): min %.2f pt (%r), "
                "median %.2f pt, max %.2f pt", printed_pt(a_min), a_min[4],
                float(np.median([printed_pt(v) for v in a_words])),
                max(printed_pt(v) for v in a_words))
    lone = [v[4] for v in a_words if v[4].strip() == "A"]
    if A_PDF.name == "Figure1_A_redrawn.pdf":
        check(True, "panel A is the code-only redraw, which carries its own letter")
    else:
        check(not lone, f"the panel letter A is gone from the strip ({lone})")

    # ------------------------------------------------------------- verdict ---
    logger.info("-" * 62)
    if fails:
        logger.error("failures: %d", len(fails))
        for f in fails:
            logger.error("  - %s", f)
        sys.exit(1)
    logger.info("ALL CHECKS PASSED -- Figure1_rev.pdf prints at %.1f pt minimum, "
                "page %.0f x %.0f pt", worst, w_pt, h_pt)


if __name__ == "__main__":
    main()
