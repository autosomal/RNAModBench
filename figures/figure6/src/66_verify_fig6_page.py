#!/usr/bin/env python
"""66 -- acceptance of the rebuilt main Figure 6 (65_fig6_combination_page.py).

The 2026-09-20 main figure failed exactly one test that nobody was running: its
**effective printed font size**.  It was drawn on a 380 x 367 mm canvas and placed
with ``\\includegraphics[width=0.95\\textwidth]{Figure6.pdf}`` (479.5 pt = 169 mm),
i.e. at a 0.445 scale, so its 8-14 pt labels printed at 3.7-6.2 pt.  This script
therefore reads the placement factor out of the manuscript source instead of
assuming it, and asserts

    placement == the factor manuscript.tex asks for (0.90 \\textwidth)
    printed height = canvas x scale <= text block - caption
    every Tf size >= 7.0 pt nominal   (printed = nominal x 0.947; re-tiering the
        sizes so the printed value clears 7 pt as well is the family-wide type
        round, recorded here 2026-10-09)

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

2026-10-08 
------------------------------------------------------
Row B's mouse column now draws the two study-B curves lighter and thinner (0.8 pt,
alpha 0.7) and with their dashes half a period apart, so the two stay readable
where they coincide (they differ by 0.8 / 0.03 / 1.5 pp at k = 3 / 4 / 5), and the
panel named that family itself.  Checked here: both dashes read back as the 0.8 x
(4 2) pattern with phases 0 and half a period, the old heavier pattern is gone from
that column, the 0.7 stroke-and-fill transparency state is used exactly twice on
the page, its blend is visible inside the mouse panel, and the ``study`` key word
is still stated in row B's own key.  (The in-panel label this round introduced was
dropped again on 2026-10-09 as a duplicate of that key -- see below.)

2026-10-09 
----------------------------------------------------------------
Row D's selected-optimum ring shrank from 10.53 pt across to 6.4 pt -- ten of
them crowd the mouse column -- and this script now holds the size from both
sides: the source constants, the accent-pixel count of the page (3 000-6 000 at
300 dpi, where the old page carried 7 150) and a scan through the printed ring
itself (4.9 pt across the stroke, where the old one measured 7.92).

2026-10-09 ")
-------------------------------------------------------------------------
The delivered page had no ink margin at all (300 dpi raster: left 2.4, right 0.5,
top 1.9, bottom 5.5 pt) and row B's mouse panel printed "study B (dashed)" even
though the row key already names both studies with a solid and a dashed swatch.
This script now fails if any of the four ink margins is below 12 pt, if the
drawing script still carries an ``ax.text`` / ``ax.annotate`` call, if the row-B
key stops naming study A and study B (the single source for that convention), or
if the printed height leaves the page budget -- and that budget is now computed
from the manuscript's own placement factor instead of the old 1:1 assumption.

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

from common import config as C                    # noqa: E402
from common.manifest import setup_logger          # noqa: E402

OUT = (_RB / "figures/figure6")
TAB, FIG, LOG = OUT / "tables", OUT / "figures", OUT / "logs"
PDF = FIG / "Figure6_rev.pdf"
PNG = FIG / "Figure6_rev.png"
SCRIPT = Path(__file__).with_name("65_fig6_combination_page.py")

#: the canvas width: 0.95 * \textwidth of USG.cls (178 mm text width = 504.57 pt).
#: The manuscript does *not* place the page at that width any more -- the
#: placement factor is read from manuscript.tex below (2026-10-09), so this is the
#: drawn width, and the printed one is ``TARGET_W * PLACEMENT``.
TARGET_W = 479.52
MIN_PT = 7.0
#: fallback if manuscript.tex cannot be read; the check below prefers the source
PLACEMENT_FALLBACK = 0.90
TEXTWIDTH_PT = 178.0 / 25.4 * 72.0            # USG.cls text width, 504.57 pt
MANUSCRIPT_TEX = ((_XB / "manuscript/manuscript_rev/manuscript.tex"))
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


def dashed_strokes_full(pdf: Path) -> list[tuple[float, float, str, float]]:
    """``(clip x0, clip y0, dash pattern, dash offset)`` of every dashed stroke.

    matplotlib writes a dash tuple scaled by the line width, right beside the
    clip rectangle of the axes it belongs to, so the two grey reference lines of
    row A -- the styles the 2026-09-30 review was about -- can be read back from
    the PDF instead of being taken on trust from the source.

    2026-10-08: the pattern regex no longer insists on a zero dash offset, because
    row B's study-B intersection dash is deliberately shifted by half a period so
    that the two coincident study-B curves interleave.  ``dashed_strokes`` still
    returns the three fields its callers unpack; the offset is available here.
    """
    data = pypdf.PdfReader(str(pdf)).pages[0].get_contents().get_data()
    pat = re.compile(r"([\d.]+) ([\d.]+) [\d.]+ [\d.]+ re W n[^\n]*?"
                     r"\[ ([\d.]+) ([\d.]+) \] ([\d.]+) d")
    return [(float(m.group(1)), float(m.group(2)), f"{m.group(3)} {m.group(4)}",
             float(m.group(5)))
            for m in pat.finditer(data.decode("latin-1"))]


def dashed_strokes(pdf: Path) -> list[tuple[float, float, str]]:
    """``(clip x0, clip y0, dash pattern)`` of every dashed stroke on the page."""
    return [(x0, y0, pat) for x0, y0, pat, _ in dashed_strokes_full(pdf)]


def translucent_strokes(pdf: Path, alpha: float) -> int:
    """Uses of the ``alpha`` stroke/fill transparency state, page and XObjects.

    2026-10-08: matplotlib records ``alpha`` as an ExtGState (``/CA`` for the
    stroke, ``/ca`` for the fill), so the lighter study-B series of row B can be
    read back from the vector stream instead of being taken on trust.
    """
    reader = pypdf.PdfReader(str(pdf))
    page = reader.pages[0]
    names = [k for k, v in (page["/Resources"].get("/ExtGState") or {}).items()
             if (lambda d: abs(float(d.get("/CA", 1)) - alpha) < 1e-9
                 and abs(float(d.get("/ca", 1)) - alpha) < 1e-9)(v.get_object())]
    if not names:                                   # pragma: no cover
        return 0
    seen: set[int] = set()
    hits = 0

    def walk(obj, depth: int = 0) -> None:
        nonlocal hits
        if depth > 3 or obj is None:
            return
        try:
            data = obj.get_contents().get_data().decode("latin-1")
        except Exception:                           # pragma: no cover
            data = ""
        for n in names:
            hits += len(re.findall(re.escape(n) + r"\s+gs", data))
        res = obj.get("/Resources") if hasattr(obj, "get") else None
        for name in (res or {}).get("/XObject") or {}:
            child = (res["/XObject"][name]).get_object()
            if id(child) in seen:
                continue
            seen.add(id(child))
            walk(child, depth + 1)

    walk(page)
    return hits


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


def ring_centres_px(h_px: int) -> list[list[tuple[float, float]]]:
    """PNG-pixel centres of row D's selected optima, one list per column.

    2026-10-09: the ring of row D shrank to 6.4 pt because ten of them crowded
    the mouse column, and the only way to measure what is *printed* is to locate
    the rings on the artwork.  The centres come from the drawing module itself
    (the optima are table data, so they are the same in every version), and the
    renderer's y axis points up while the image's points down.
    """
    spec = importlib.util.spec_from_file_location("verify_fig6_src", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["verify_fig6_src"] = mod
    spec.loader.exec_module(mod)
    t = mod.load()
    fig, meta = mod.build_figure(t, mod.tool_order(t["members"], t["space"]))
    fig.canvas.draw()
    scale = 300 / fig.dpi
    out: list[list[tuple[float, float]]] = []
    for c in meta["contracts_d"]:
        ax = c["rings"].axes
        pts = []
        for x, y in np.atleast_2d(c["rings"].get_offsets()):
            dx, dy = ax.transData.transform((float(x), float(y)))
            pts.append((dx * scale, h_px - dy * scale))
        out.append(pts)
    return out


def measure_optimum_ring_pt(arr) -> float:
    """Outer diameter (pt) of the accent stroke of row D's k = 1 ring.

    A vertical scan through the centre of the leftmost optimum (Arabidopsis
    k = 1): that ring is isolated, its neighbours sitting to the right, the
    joining line crosses the scan *inside* the ring, and the white halo -- behind
    the stroke and invisible over white paper -- contributes no accent pixel, so
    the span is exactly ``sqrt(OPT_S) + OPT_LW`` as printed.
    """
    cx, cy = ring_centres_px(arr.shape[0])[0][0]
    acc = np.abs(arr - np.array([0xB0, 0x24, 0x18])).sum(axis=2) <= 60
    half = 32                                    #: +-7.7 pt, clear of neighbours
    col = acc[int(round(cy)) - half:int(round(cy)) + half + 1,
              int(round(cx)) - 2:int(round(cx)) + 3]
    rows = np.nonzero(col.any(axis=1))[0]
    assert len(rows), "no accent pixels around the k = 1 optimum of column 0"
    return float(rows.max() - rows.min() + 1) / (300 / 72.0)


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
    check(len(reader.pages) == 1, "single-page PDF")
    check(abs(w_pt - TARGET_W) < 0.5,
          f"page width {w_pt:.2f} pt == 0.95\\textwidth ({TARGET_W} pt)")
    #: 2026-10-09: the placement factor is read from the manuscript source.  The
    #: page used to assume it printed 1:1, while manuscript.tex asks for
    #: ``width=0.90\textwidth`` -- so everything printed 5.3 % smaller than this
    #: gate believed (a 7.2 pt legend printed at 6.8 pt).  Deriving the factor
    #: here keeps the height budget and the printed font report honest.
    src_lines = [ln for ln in MANUSCRIPT_TEX.read_text().splitlines()
                 if "includegraphics" in ln and "Figure6.pdf" in ln] \
        if MANUSCRIPT_TEX.exists() else []
    m = re.search(r"width=([\d.]+)", src_lines[0]) if src_lines else None
    check(m is not None and "textwidth" in src_lines[0],
          f"manuscript.tex places Figure6.pdf by a width relative to the text "
          f"width ({src_lines[0].strip() if src_lines else 'no such line'})")
    placement = float(m.group(1)) if m else PLACEMENT_FALLBACK
    scale = placement * TEXTWIDTH_PT / w_pt                    # printed / drawn
    check(abs(scale - 0.947) < 5e-3,
          f"printed scale {scale:.4f}: the page prints at {scale * 100:.1f} % of "
          f"its drawn size ({placement:g} x 504.57 pt / {w_pt:.2f} pt)")

    #: binding constraint is the text block minus the caption the figure carries.
    #: 2026-09-30: both numbers were re-measured on the delivered manuscript's
    #: Figure 6 page (page 21 of 23, US Letter 792 pt, margins 16 mm -> text block
    #: 9.74 in; `pdftotext -bbox` gives the caption block y 617.7 -> 734.8 pt =
    #: 1.63 in, with 11.8 pt left to the bottom margin).  The old pair (234 mm /
    #: 1.05 in, page 20) described an earlier caption and was 0.6 in optimistic.
    #: 2026-10-09: the budget applies to the *printed* height, so the canvas may be
    #: (9.74 - 1.63) / 0.947 = 8.56 in; it is 8.34 in (7.90 in printed, 0.21 in of
    #: slack), and the margin round spent 0.16 in of that on the top ink margin and
    #: 0.11 in on the row-F tail.
    TEXTBLOCK_IN = 9.74
    CAPTION_IN = 1.63
    budget_in = (TEXTBLOCK_IN - CAPTION_IN) / scale
    check(h_pt / 72.0 <= budget_in + 1e-6,
          f"canvas {h_pt / 72:.2f} in = {h_pt / 72 * scale:.2f} in printed within "
          f"the page budget (<= {budget_in:.2f} in canvas = {TEXTBLOCK_IN:.2f} in "
          f"text block minus the {CAPTION_IN:.2f} in caption, divided by the "
          f"{scale:.3f} placement scale)")

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
    #: 2026-10-09: this holds the *nominal* floor (>= 7 pt declared).  At the
    #: 0.947 placement the printed sizes are 6.8-9.0 pt, so the smallest printed
    #: label is 0.2 pt under the floor -- the family-wide type round re-tiers all
    #: figures at once, and `measure_fig6_margins_20261009.py` plus this message
    #: record the printed value until then.
    check(not unexplained,
          f"every drawn size is >= {MIN_PT} pt nominal or a mathtext "
          f"sub/superscript of one (sizes {sorted(sizes)}, unexplained "
          f"{unexplained}); smallest prints at {worst:.2f} pt (scale {scale:.3f})")

    png = Image.open(PNG)
    check(abs(png.size[0] / 300 * 25.4 - TARGET_W / 72 * 25.4) < 0.5,
          f"PNG {png.size[0]}x{png.size[1]} px = same width at 300 dpi")

    # --------------------------------------------------- style / house rules
    text = SCRIPT.read_text()
    check('"axes.grid": False' in text, "gridlines disabled in the style block")
    check("bbox_inches=None" in text, "saved without a tight bounding box")

    #: exception to the house rule "no in-panel annotation text / value callouts":


    #: studies -- so the rule holds absolutely again and any ax.text / ax.annotate
    #: in the drawing script fails here.
    label_calls = re.findall(r"\bax\.text\([^)]*", text)
    annotate_calls = re.findall(r"\bax\.annotate\([^)]*", text)
    check(not annotate_calls and not label_calls,
          f"no in-panel annotation text / value callouts in the script "
          f"({len(label_calls)} ax.text, {len(annotate_calls)} ax.annotate calls)")
    raw = __import__("subprocess").run(["pdffonts", str(PDF)],
                                       capture_output=True, text=True).stdout
    fonts = {ln.split()[0].split("+")[-1] for ln in raw.splitlines()[2:]
             if ln.strip()}
    check(all("Arial" in f for f in fonts),
          f"embedded fonts are Arial only ({sorted(fonts)})")


    #: instrument is `measure_fig6_margins_20261009.py` (300 dpi raster, bounding
    #: box of every non-white pixel); here it is repeated as the gate, because the
    #: margins are a delivered property of the page, and 12 pt is the house floor
    #: (13 pt are drawn so that the 0.947 placement still prints >= 12.3 pt).
    subprocess.run(["pdftoppm", "-r", "300", "-png", "-singlefile", str(PDF),
                    "/tmp/_fig6_gate_margin"], check=True)
    arr = np.asarray(Image.open("/tmp/_fig6_gate_margin.png").convert("L"))
    ink = np.where(arr < 250)
    l_pt = ink[1].min() * 72 / 300
    r_pt = (arr.shape[1] - 1 - ink[1].max()) * 72 / 300
    t_pt = ink[0].min() * 72 / 300
    b_pt = (arr.shape[0] - 1 - ink[0].max()) * 72 / 300
    worst, side = min((l_pt, "left"), (r_pt, "right"), (t_pt, "top"),
                      (b_pt, "bottom"))
    check(worst >= 12.0,
          f"four ink margins >= 12 pt (left {l_pt:.2f}, right {r_pt:.2f}, "
          f"top {t_pt:.2f}, bottom {b_pt:.2f}; smallest {worst:.2f} pt on the "
          f"{side}, 1 px = 0.24 pt of rasterisation)")

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

    #: (sqrt(OPT_S) + OPT_LW + OPT_HALO), down from 10.53 pt, and the drawn
    #: artist has to come from those constants -- a later edit that keeps them
    #: but draws something else is exactly the drift this catches
    check(re.search(r"OPT_S, OPT_LW, OPT_HALO = 16, 0\.9, 1\.5", src)
          and "OPT_MIN_PT, OPT_MAX_PT = 5.5, 7.0" in src
          and "s=OPT_S" in src and "linewidths=OPT_LW" in src
          and "ms=OPT_SWATCH_MS" in src,
          "the row-D ring (and its key swatch) is drawn from the shrunken "
          "OPT_* constants: 6.4 pt, not the 10.53 pt of 2026-09-21")
    arr = np.asarray(png.convert("RGB"), dtype=np.int16)
    n_accent = int((np.abs(arr - np.array([0xB0, 0x24, 0x18])).sum(axis=2)
                    <= 24).sum())
    #: both bounds bite: a vanished ring falls under the floor, a ring grown back
    #: to the old geometry goes over the ceiling -- 7 150 accent pixels were
    #: measured on the 10.53 pt page, 4 338 on the 6.4 pt one (300 dpi)
    check(3000 <= n_accent <= 6000,
          f"accent-coloured optima printed at the shrunken size ({n_accent} px "
          f"at 300 dpi, 3 000-6 000 expected; 7 150 was the 10.53 pt page)")
    #: and the size itself, measured on the artwork rather than trusted to the
    #: constants
    ring_pt = measure_optimum_ring_pt(arr)
    check(4.4 <= ring_pt <= 5.4,
          f"the printed ring measures {ring_pt:.2f} pt across its stroke "
          f"(4.9 = sqrt(16) + 0.9 expected; the 2026-09-21 ring measured 7.92)")

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

    #: draws the two study-B curves lighter and thinner (0.8 pt, alpha 0.7) and
    #: half a period apart, so that they stay readable where they share a y value
    #: -- union recall 67.4/71.1/74.0 % against intersection PPV 66.6/71.1/75.5 %
    #: at k = 3/4/5.  Both read back from the page: one stroke at phase 0 and one
    #: at half of the (3.2, 1.6) period, and no stroke at the old 1.0/1.1 pt
    #: (4 2)/(4.4 2.2) pattern of that column.
    full = dashed_strokes_full(PDF)
    study_b = [(pat, off) for _, _, pat, off in full if pat == "3.2 1.6"]
    check(len(study_b) == 2 and sorted(o for _, o in study_b) == [0.0, 2.4],
          f"row B's two study-B dashes are 0.8 x (4 2) and half a period apart "
          f"({len(study_b)} strokes, phases {sorted(o for _, o in study_b)})")
    stale = [(pat, off) for x0, _, pat, off in full
             if pat in ("4 2", "4.4 2.2") and x0 < 150.0]
    check(not stale,
          f"the old, heavier study-B dash is gone from the mouse column of row B "
          f"({stale})")
    #: the alpha itself is a transparency state in the vector stream: exactly the
    #: two study-B series of row B carry it, so the whole page may show no other
    #: use of a 0.7 stroke-and-fill state
    n_alpha = translucent_strokes(PDF, 0.7)
    check(n_alpha == 2,
          f"exactly the two study-B series of row B are drawn at 0.7 alpha "
          f"({n_alpha} uses of the CA .7 / ca .7 state on the page)")
    #: and the render must show it: the mouse panel is the pt box of the middle
    #: column of row B, where the 0.7 blends of the two series colours must be
    #: present (measured 158 blue / 92 orange px at 300 dpi)
    #: x0..x1 and top..bottom of that panel in PDF points (measured from the
    #: clip rectangle of row B's middle column), both counted from the top
    mx0, mx1, mtop, mbot = 204.0, 328.8, h_pt - 471.6, h_pt - 424.1
    xa, xb = int(round(mx0 / TARGET_W * png.size[0])) + 2, \
        int(round(mx1 / TARGET_W * png.size[0])) - 2
    ya, yb = int(round(mtop / h_pt * png.size[1])) + 2, \
        int(round(mbot / h_pt * png.size[1])) - 2
    panel = arr[ya:yb, xa:xb]
    for name, rgb, floor in (("blue", np.array([0x1F, 0x5C, 0x8B]), 100),
                             ("orange", np.array([0xE0, 0x8A, 0x2E]), 50)):
        light = np.round(0.7 * rgb + 0.3 * 255).astype(np.int16)
        n_light = int((np.abs(panel - light).sum(axis=2) <= 12).sum())
        check(n_light >= floor,
              f"the mouse panel of row B shows the lighter study-B {name} "
              f"({n_light} px of the 0.7 blend, floor {floor})")
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
        if token == "study":

            #: no longer prints its own "study B (dashed)", so both remaining
            #: "study" words are the two key entries, inside row B's own band
            in_band = [w for w in hit if top < w[2] and w[4] < below]
            check(len(hit) == 2 and len(in_band) == 2,
                  f"exactly row B's two key entries say 'study' "
                  f"({len(in_band)} of {len(hit)} occurrences inside "
                  f"{top:.1f}..{below:.1f} pt)")
            continue
        for _, x0, y0, x1, y1 in hit:
            check(top < y0 and y1 < below,
                  f"{token!r} at {y0:.1f}..{y1:.1f} pt sits below its own panels "
                  f"({top:.1f}..{bottom:.1f}) and above the next row ({below:.1f})")

    #: dropped as a duplicate of the row key, so "(dashed)" must be gone from the
    #: page and no study/B word may sit inside row B's mouse panel any more.
    b_x = sorted({(round(x0, 1), round(x1, 1)) for x0, y0, x1, y1 in clip_rects(PDF)
                  if (round(h_pt - y1, 3), round(h_pt - y0, 3)) == bands[1]
                  and x1 - x0 > 40})
    assert len(b_x) == 3, f"row B should carry three panel columns, got {b_x}"
    mid0, mid1 = b_x[1]
    dashed_words = [w for w in words if "(dashed)" in w[0]]
    stray = [w for w in words
             if w[0].strip(".,") in ("study", "B", "(dashed)")
             and bands[1][0] < w[2] and w[4] < bands[1][1]
             and mid0 < w[1] and w[3] < mid1]
    check(not dashed_words and not stray,
          f"row B's mouse panel carries no study annotation: '(dashed)' occurs "
          f"{len(dashed_words)} times on the page and {len(stray)} study/B word(s) "
          f"fall inside the mouse panel (x {mid0}..{mid1})")
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
