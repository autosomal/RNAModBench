#!/usr/bin/env python
"""66 -- acceptance of the rebuilt main Figure 6 (65_fig6_combination_page.py).

The 2026-09-20 main figure failed exactly one test that nobody was running: its
**effective printed font size**.  It was drawn on a 380 x 367 mm canvas and placed
with ``\\includegraphics[width=0.95\\textwidth]{Figure6.pdf}`` (479.5 pt = 169 mm),
i.e. at a 0.445 scale, so its 8-14 pt labels printed at 3.7-6.2 pt.  This script
therefore computes the placement scale from the PDF page width itself
(``scale = 0.95 * \\textwidth / page width``) and asserts

    scale == 1.0            (the canvas *is* the printed size)
    every Tf size * scale >= 7.0 pt

plus the anchors the figure claims (recall/PPV/p0/membership/effects/controls) and
the geometry rules of the house style (Arial-only embedding, no gridlines, no
in-panel annotation text, 300 dpi PNG).

2026-09-30 
---------------------------------------------------
Adds the key-vs-artwork checks: the row-A reference lines are read back from the
page as dash patterns (thin dashed = recall against every GLORI site, thin dotted
= best single tool, in *every* column), the chance level of row D must be drawn
in all three columns *and* visible inside each panel (counted as 0.30-grey pixels
of the rendered PNG), the retired "all GLORI (dashed)" wording must be gone, and
the grey curve is named "share of all GLORI sites".

2026-09-30 
------------------------------------------------------------------------------------
Adds the placement checks: the six panel bands are read back from the page's clip
rectangles and one distinctive word of each row key must sit below its own band
and above the next one; the row-E key must not print the k of its sets (the values
live in the caption).  Every one of the six rows must print the three species
titles, and the page-budget constants below are re-measured on the delivered
manuscript's Figure 6 page instead of the older page-20 estimate.

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    figures/figure6/src/66_verify_fig6_page.py
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

import numpy as np
import pandas as pd
import pypdf
from PIL import Image

HERE = _RB / "src/harmonisation"
sys.path.insert(0, str(HERE))

from common import config as C                    # noqa: E402
from common.manifest import setup_logger          # noqa: E402

OUT = (_RB / "figures/figure6")
TAB, FIG, LOG = OUT / "tables", OUT / "figures", OUT / "logs"
PDF = FIG / "Figure6_rev.pdf"
PNG = FIG / "Figure6_rev.png"
SCRIPT = Path(__file__).with_name("65_fig6_combination_page.py")

#: 0.95 * \textwidth of USG.cls (178 mm text width); the manuscript places the
#: figure with exactly this width, so this is the printed width in points
TARGET_W = 479.52
MIN_PT = 7.0
#: panel letters of the rebuilt figure, in order
LETTERS = ("A", "B", "C", "D", "E", "F")


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


def dashed_strokes(pdf: Path) -> list[tuple[float, float, str]]:
    """``(clip x0, clip y0, dash pattern)`` of every dashed stroke on the page.

    matplotlib writes a dash tuple scaled by the line width, right beside the
    clip rectangle of the axes it belongs to, so the two grey reference lines of
    row A -- the styles the 2026-09-30 review was about -- can be read back from
    the PDF instead of being taken on trust from the source.
    """
    data = pypdf.PdfReader(str(pdf)).pages[0].get_contents().get_data()
    pat = re.compile(r"([\d.]+) ([\d.]+) [\d.]+ [\d.]+ re W n[^\n]*?"
                     r"\[ ([\d.]+) ([\d.]+) \] 0 d")
    return [(float(m.group(1)), float(m.group(2)), f"{m.group(3)} {m.group(4)}")
            for m in pat.finditer(data.decode("latin-1"))]


def clip_rects(pdf: Path) -> list[tuple[float, float, float, float]]:
    """``(x0, y0, x1, y1)`` of every axes clip rectangle, in PDF points.

    The panel boxes are the only rectangles matplotlib clips a panel to, so the
    row bands of the page can be read back from them (bottom-left origin, as the
    PDF stores them).
    """
    data = pypdf.PdfReader(str(pdf)).pages[0].get_contents().get_data()
    out = []
    for m in re.finditer(r"([\d.]+) ([\d.]+) ([\d.]+) ([\d.]+) re W n",
                         data.decode("latin-1")):
        x0, y0, w, h = (float(m.group(i)) for i in range(1, 5))
        out.append((x0, y0, x0 + w, y0 + h))
    return out


def page_words(pdf: Path) -> list[tuple[str, float, float, float, float]]:
    """``(word, x0, y0, x1, y1)`` of every word, in points from the page's top left."""
    bbox = subprocess.run(["pdftotext", "-bbox", str(pdf), "-"],
                          capture_output=True, text=True).stdout
    return [(m.group(5), float(m.group(1)), float(m.group(2)),
             float(m.group(3)), float(m.group(4)))
            for m in re.finditer(
                r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" '
                r'yMax="([\d.]+)">([^<]*)</word>', bbox)]


def main() -> None:
    logger = setup_logger("66_verify_fig6_page", log_dir=LOG)
    fails: list[str] = []

    checks_run: list[int] = []          #: so the summary counts what actually ran

    def check(ok: bool, msg: str) -> None:
        checks_run.append(1)
        (logger.info if ok else logger.error)("%s %s", "PASS" if ok else "FAIL",
                                              msg)
        if not ok:
            fails.append(msg)

    # ---------------------------------------------------------------- geometry
    reader = pypdf.PdfReader(str(PDF))
    box = reader.pages[0].mediabox
    w_pt, h_pt = float(box.width), float(box.height)
    scale = TARGET_W / w_pt
    check(len(reader.pages) == 1, "single-page PDF")
    check(abs(w_pt - TARGET_W) < 0.5,
          f"page width {w_pt:.2f} pt == 0.95\\textwidth ({TARGET_W} pt)")
    check(abs(scale - 1.0) < 2e-3,
          f"placement scale {scale:.4f} == 1.000 (no down-scaling in the "
          f"manuscript: {TARGET_W / w_pt * 100:.1f} %)")

    #: binding constraint is the text block minus the caption the figure carries.
    #: 2026-09-30: both numbers were re-measured on the delivered manuscript's
    #: Figure 6 page (page 21 of 23, US Letter 792 pt, margins 16 mm -> text block
    #: 9.74 in; `pdftotext -bbox` gives the caption block y 617.7 -> 734.8 pt =
    #: 1.63 in, with 11.8 pt left to the bottom margin).  The old pair (234 mm /
    #: 1.05 in, page 20) described an earlier caption and was 0.6 in optimistic.
    #: The canvas may therefore be 9.74 - 1.63 = 8.11 in; it is 8.07 in.
    TEXTBLOCK_IN = 9.74
    CAPTION_IN = 1.63
    budget_in = TEXTBLOCK_IN - CAPTION_IN
    check(h_pt / 72.0 <= budget_in + 1e-6,
          f"height {h_pt / 72:.2f} in within the page budget (<= {budget_in:.2f} in: "
          f"US Letter text block {TEXTBLOCK_IN:.2f} in minus the caption "
          f"{CAPTION_IN:.2f} in measured on page 21)")

    # ------------------------------------------------------ effective font size

    #: scale mathtext gives sub- and superscripts (``$p_0$``, ``$10^6$`` ...), so it
    #: is typography by construction rather than a small label.  Every size must be
    #: at least MIN_PT, or be that sub/superscript scale of a size that is.
    sizes = tf_sizes(PDF)
    eff = {s: round(s * scale, 2) for s in sizes}
    worst = min(eff.values()) if eff else float("nan")
    bases = [s for s in sizes if s >= MIN_PT - 1e-6]
    unexplained = [s for s in sorted(sizes)
                   if s < MIN_PT - 1e-6
                   and not any(abs(s - 0.7 * b) < 0.02 * b for b in bases)]
    check(not unexplained,
          f"smallest effective font {worst:.2f} pt: every drawn size is >= {MIN_PT} pt "
          f"or a mathtext sub/superscript of one (sizes {sorted(sizes)}, "
          f"unexplained {unexplained})")

    png = Image.open(PNG)
    check(abs(png.size[0] / 300 * 25.4 - TARGET_W / 72 * 25.4) < 0.5,
          f"PNG {png.size[0]}x{png.size[1]} px = same width at 300 dpi")

    # --------------------------------------------------- style / house rules
    text = SCRIPT.read_text()
    check('"axes.grid": False' in text, "gridlines disabled in the style block")
    check("bbox_inches=None" in text, "saved without a tight bounding box")
    check(not re.search(r"\bax\.annotate\(|\bax\.text\(", text),
          "no in-panel annotation text / value callouts in the script")
    raw = __import__("subprocess").run(["pdffonts", str(PDF)],
                                       capture_output=True, text=True).stdout
    fonts = {ln.split()[0].split("+")[-1] for ln in raw.splitlines()[2:]
             if ln.strip()}
    check(all("Arial" in f for f in fonts),
          f"embedded fonts are Arial only ({sorted(fonts)})")

    # ------------------------------------------------------------- anchors
    sel = pd.read_csv(TAB / "fig6_combination_selected.tsv", sep="\t")
    mem = pd.read_csv(TAB / "fig6_selected_members.tsv", sep="\t")
    allk = pd.read_csv(TAB / "fig6_per_unit_allk.tsv", sep="\t")
    eff_t = pd.read_csv(TAB / "fig6_tool_effects.tsv", sep="\t")
    hist = pd.read_csv(TAB / "figS6_coverage_hist.tsv", sep="\t")
    ctl = pd.read_csv(TAB / "fig6_negative_control_fp.tsv", sep="\t")
    want_k5 = {"Arabidopsis_WT": (46.7, 14.1), "studyA": (80.1, 6.8),
               "studyB": (74.0, 6.3), "HeLa_WT": (56.8, 7.5)}
    for g, (rec, ppv) in want_k5.items():
        r = sel[(sel.group == g) & (sel.k == 5)].iloc[0]
        check(abs(100 * r.union_recall_mean - rec) < 0.1
              and abs(100 * r.union_precision_mean - ppv) < 0.1,
              f"{g} k=5 anchor: recall {100 * r.union_recall_mean:.1f} % / "
              f"PPV {100 * r.union_precision_mean:.1f} %")
    for g in want_k5:
        s = sel[sel.group == g]
        p0 = 100 * float(s.p0_chance_precision.iloc[0])
        exp = 1.5 if g == "Arabidopsis_WT" else (0.85 if g == "HeLa_WT" else 0.7)
        check(abs(p0 - exp) < 0.06, f"{g} chance level p0 = {p0:.2f} %")
        check(s.union_recall_all_mean.notna().all(),
              f"{g}: recall against every GLORI site present (both denominators)")
        check(len(mem[mem.group == g]) == 65, f"{g}: membership 65 rows")
        check(len(eff_t[eff_t.group == g]) == 65, f"{g}: tool effects 65 rows")
        check(sorted(set(ctl[ctl.chosen_for == g].k)) == [1, 2, 3, 4, 5],
              f"{g}: control burden covers k = 1..5")
        a = allk[allk.group == g]
        got = 100 * a[a.k == 5].union_recall.mean()
        ref = 100 * float(s[s.k == 5].union_recall_mean.iloc[0])
        check(abs(got - ref) < 1e-3,
              f"{g}: per-unit mean of k=5 matches the selected table "
              f"({got:.3f} vs {ref:.3f})")
    check(set(hist.set) >= {"single (k=1)", "union marginal", "union (k=5)",
                            "intersection (k=2)", "GLORI reference"},
          "coverage CDF covers all five site sets")
    check(len(hist[hist.set == "union marginal"]) > 0
          and len(hist[hist.set == "intersection (k=5)"]) > 0,
          "coverage histograms exist for the marginal and intersection sets")
    lim = 100 * float(np.nanmax(np.abs(eff_t.effect_recall)))
    check(lim > 1.0, f"tool-effect colour scale is finite (+/- {lim:.1f} pp)")

    # ------------------------------- 2026-09-21 readability fixes ------------
    # (a) row C prints the 13 tool names once, in the left margin of column 0:
    #     the three columns share those rows, and a second and third set of 0.69
    #     in labels landed on the effect matrix of the column to their left
    body = __import__("subprocess").run(["pdftotext", str(PDF), "-"],
                                        capture_output=True, text=True).stdout
    names = sorted(set(mem.tool))
    seen = {n: len(re.findall(rf"(?<![A-Za-z0-9_]){re.escape(n)}"
                              rf"(?![A-Za-z0-9_])", body)) for n in names}
    check(all(v == 1 for v in seen.values()),
          f"each of the {len(names)} tool names is printed exactly once (row C, "
          f"column 0): "
          f"{ {k: v for k, v in seen.items() if v != 1} or 'all 1'}")
    # (b) row D marks the selected optimum of every k with an accent ring and a
    #     white halo, and that accent has to survive into the printed artwork
    src = SCRIPT.read_text()
    check('HILITE = "#b02418"' in src and "withStroke" in src,
          "row D draws the selected optima as accent rings with a white halo")
    arr = np.asarray(png.convert("RGB"), dtype=np.int16)
    n_accent = int((np.abs(arr - np.array([0xB0, 0x24, 0x18])).sum(axis=2)
                    <= 24).sum())
    check(n_accent > 400,
          f"accent-coloured optima actually printed ({n_accent} px at 300 dpi)")

    # -------------------------------- 2026-09-30: key <-> artwork, one to one --
    # (user review) The row-A key said "all GLORI (dashed)" while that line was
    # drawn solid in two of the three columns, its swatch was told apart from the
    # best-single-tool swatch by dash density alone, and the row-D key named a
    # chance level that two of the three panels had clipped below the axis.
    #: the patterns are unique to their line on this page, so they are counted
    #: page-wide: the rows moved down when the keys moved below their panels
    strokes = [p for _, _, p in dashed_strokes(PDF)]
    n_glori_dash = strokes.count("3.6 1.8")
    check(n_glori_dash == 4,
          f'"share of all GLORI sites" drawn as a 0.9 x (4 2) dashed line in '
          f"every column ({n_glori_dash} strokes, expected 4: Arabidopsis, mouse "
          f"study A, mouse study B, human)")
    n_best_dot = strokes.count("1.44 1.44")
    check(n_best_dot == 4,
          f'"best single tool" drawn as a 0.9 x (1.6 1.6) dotted line in every '
          f"column ({n_best_dot} strokes, expected 4)")
    check("2 1.2" not in strokes,
          "the retired grey dash of the best-single-tool line is gone (it was "
          "the style indistinguishable from all GLORI in the key)")
    n_p0 = strokes.count("2.7 1.8")
    check(n_p0 == 3, f"row D draws the chance level p0 in all three columns "
                     f"({n_p0} strokes)")
    #: drawn is not enough -- p0 must be *visible*: the dashed line was clipped
    #: away in Arabidopsis and mouse until 2026-09-30
    grey = np.abs(arr - 76).max(axis=2) <= 8        # the 0.30 grey of that line
    for name, x0 in (("Arabidopsis", 57.6), ("Mouse", 204.0), ("Human", 350.4)):
        xa = int(round(x0 / TARGET_W * png.size[0])) + 3
        xb = int(round((x0 + 124.8) / TARGET_W * png.size[0])) - 3
        ya = int(round((h_pt - (192.96 + 47.52)) / h_pt * png.size[1])) + 3
        yb = int(round((h_pt - 192.96) / h_pt * png.size[1])) - 3
        n_px = int(grey[ya:yb, xa:xb].sum())
        check(n_px > 150,
              f"row D {name}: the dashed chance level is inside the panel "
              f"({n_px} px of 0.30 grey)")
    check("share of all GLORI sites" in body and "all GLORI (dashed)" not in body,
          "the row-A key names the grey curve by what it measures ('share of all "
          "GLORI sites'), not by a line style or a bare denominator")
    check("union marginal" in body and "k = 5" not in body and "k = 2" not in body,
          "the site-set key names the sets, not their k (the k values live in the "
          "caption: 5 tools for the union, 2 for the intersection)")
    n_a, n_b = body.count("mouse study A"), body.count("mouse study B")
    check(n_a == 1 and n_b == 1,
          f"the mouse study convention is stated once, in row B's key "
          f"(A x{n_a}, B x{n_b})")


    #: printed below the panels it belongs to.  Read the six panel bands back from
    #: the page's clip rectangles and locate one distinctive word of each key.
    bands = sorted({(round(h_pt - y1, 3), round(h_pt - y0, 3))
                    for x0, y0, x1, y1 in clip_rects(PDF)
                    if x1 - x0 > 40 and y1 - y0 > 21})        # a panel, not a stray
    check(len(bands) == 6,
          f"the page carries six rows of panels ({len(bands)} bands found)")
    words = page_words(PDF)
    for token, index in (("share", 0), ("study", 1), ("chance", 3),
                         ("marginal", 4), ("Curlcake", 5)):
        hit = [w for w in words if w[0].strip(".,") == token]
        check(bool(hit), f"key word {token!r} is printed ({len(hit)} times)")
        top, bottom = bands[index] if index < len(bands) else (0.0, 0.0)
        below = bands[index + 1][0] if index + 1 < len(bands) else h_pt
        for _, x0, y0, x1, y1 in hit:
            check(top < y0 and y1 < below,
                  f"{token!r} at {y0:.1f}..{y1:.1f} pt sits below its own panels "
                  f"({top:.1f}..{bottom:.1f}) and above the next row ({below:.1f})")
    titles = {n: len(re.findall(rf"(?<![A-Za-z]){n}(?![A-Za-z])", body))
              for n in ("Arabidopsis", "Mouse", "Human")}
    check(titles == {"Arabidopsis": 6, "Mouse": 6, "Human": 7},
          f"species titles printed on all six rows: {titles} (six rows x three "
          f"columns, plus the row-F key entry 'Human IVT')")

    # --------------------------------------------------------------- verdict
    tokens = set(__import__("subprocess").run(
        ["pdftotext", str(PDF), "-"], capture_output=True, text=True
    ).stdout.split())
    for letter in LETTERS:
        check(letter in tokens, f"panel letter {letter} is printed in the figure")
    logger.info("-" * 62)
    if fails:
        logger.error("failures: %d", len(fails))
        for f in fails:
            logger.error("  - %s", f)
        sys.exit(1)

    #: mathtext sub/superscript scale (5.2 pt = 7.5 x 0.7) and read like a small
    #: label.  It now quotes the smallest *base* size and says what the smaller one
    #: is.  (The check count is derived from the checks themselves, not hard-wired.)
    logger.info("ALL CHECKS PASSED (%d) -- Figure6_rev.pdf base type >= %.1f pt, "
                "mathtext sub/superscript %.1f pt",
                len(checks_run), min(bases) if bases else float("nan"), worst)


if __name__ == "__main__":
    main()
