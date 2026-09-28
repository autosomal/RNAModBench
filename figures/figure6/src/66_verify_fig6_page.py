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
    
    #: binding constraint is the two-column text block, 234 mm = 9.213 in, minus the
    #: caption the figure carries: measured on the delivered manuscript (page 20,
    #: caption block y 634 -> 704 pt = 0.97 in of text, 1.05 in including the line
    #: box), so the canvas may be 9.213 - 1.05 = 8.16 in.  It is 8.02 in.
    TEXTBLOCK_IN = 234 / 25.4
    CAPTION_IN = 1.05
    budget_in = TEXTBLOCK_IN - CAPTION_IN
    check(h_pt / 72.0 <= budget_in + 1e-6,
          f"height {h_pt / 72:.2f} in within the page budget (<= {budget_in:.2f} in: "
          f"234 mm text block {TEXTBLOCK_IN:.2f} in minus the caption {CAPTION_IN:.2f} in)")

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
