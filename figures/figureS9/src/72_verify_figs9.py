#!/usr/bin/env python
"""72 -- acceptance checks for the rebuilt Figure S9 (run after 70/71).

Contract of the 2026-10-09 page: **1015.776 x 792.0 pt** (14.108 x 11.000 in)
with the 3 x 3 grid A | B | C / D | E | H / F | G | I and nine **square** boxes
(2.675 x 2.675 in, all derived from ``s9_style`` -- A's 23-row pitch pins the
side, the sheet was widened to reach it) -- **RNA004 data only** (no RNA002 run
is drawn anywhere), every panel's numbers reproduce the frozen tables, and the
drawing code keeps the house rules (no grid, no scatter clouds, Arial only,
nothing below 8 pt -- the sweep key; the layout gate's own default floor is
7 pt).
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import s9_data as D  # noqa: E402
import s9_panels as P  # noqa: E402
import s9_style as S  # noqa: E402

PDF = S.FIGS / "FigureS9_rev.pdf"

#: 792.0 pt since 2026-10-09, so this check can never fall out of step with the
#: builder
PAGE_PT = (S.PAGE_IN[0] * 72.0, S.PAGE_IN[1] * 72.0)
KEYS = list("ABCDEFGHI")                        # 2026-09-27: three-by-three grid
                                                # (sweep D|E, effect F|G, the two
                                                # stability cells H|I)
fails: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"[{'ok' if ok else 'FAIL'}] {name}{' ' + detail if detail else ''}")
    if not ok:
        fails.append(name)


def text_of(path: Path) -> str:
    return subprocess.run(["/usr/bin/pdftotext", str(path), "-"],
                          capture_output=True, text=True).stdout


def word_boxes(pdf: Path) -> list[tuple[float, float, float, float, str]]:
    """(x0, y0, x1, y1, text) of every word box, in pt (``pdftotext -bbox``)."""
    raw = subprocess.run(["/usr/bin/pdftotext", "-bbox", str(pdf), "-"],
                         capture_output=True, text=True).stdout
    return [(float(a), float(b), float(c), float(d), e) for a, b, c, d, e in re.findall(
        r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">(.*?)</word>',
        raw)]


def main() -> None:
    info = subprocess.run(["/usr/bin/pdfinfo", str(PDF)], capture_output=True,
                          text=True).stdout
    size = [ln for ln in info.splitlines() if ln.startswith("Page size")][0]
    got = tuple(float(x) for x in re.findall(r"[\d.]+", size)[:2])
    check("page size", abs(got[0] - PAGE_PT[0]) < 0.01
          and abs(got[1] - PAGE_PT[1]) < 0.01, f"{got}")
    check("single page", "Pages:           1" in info)

    
    #: squares.  Nothing else can see a box that came out a hair narrow: the
    #: layout gate measures text, pdftotext sees no geometry, and ``cell_rect``
    #: used to clamp silently.  This check re-derives the square from the style
    #: module (column width == frozen label column + right margin + box side;
    #: row 1's band draws the same side) and, second, that every rendered panel
    #: canvas still is the width the placement table promises.
    box_w = S.COL_W - S.LABEL_COL_IN - S.BOX_RIGHT_IN
    check("nine plotting areas are 2.675 in squares",
          abs(box_w - S.BOX_H_LOW_IN) < 1e-9
          and abs(S.BOX_SIDE_IN - S.BOX_H_LOW_IN) < 1e-9
          and abs(S.ROW1_BOX_FRAC * S.ROW_H["r1"] - S.BOX_H_LOW_IN) < 1e-9
          and abs(S.COL_W - (S.LABEL_COL_IN + S.BOX_RIGHT_IN + S.BOX_SIDE_IN)) < 1e-9,
          f"box {box_w:.4f} x {S.BOX_H_LOW_IN:.4f} in on a "
          f"{S.COL_W:.4f} in cell, page {S.PAGE_IN[0]:.4f} in")
    panel_bad = []
    for key, (w, h) in S.PIECE.items():
        p_info = subprocess.run(["/usr/bin/pdfinfo",
                                 str(S.PANELS / f"figS9{key}.pdf")],
                                capture_output=True, text=True).stdout
        p_size = [ln for ln in p_info.splitlines() if ln.startswith("Page size")]
        p_got = tuple(float(x) for x in re.findall(r"[\d.]+", p_size[0])[:2]) \
            if p_size else (0.0, 0.0)
        if abs(p_got[0] - w * 72.0) > 0.01 or abs(p_got[1] - h * 72.0) > 0.01:
            panel_bad.append(f"{key} {p_got} != ({w * 72:.2f}, {h * 72:.2f})")
    check("panel canvases match their print size", not panel_bad,
          str(panel_bad))

    fonts = subprocess.run(["/usr/bin/pdffonts", str(PDF)],
                           capture_output=True, text=True).stdout
    names = {ln.split()[0].split("+")[-1] for ln in fonts.splitlines()[2:] if ln}
    check("Arial only", all(n.startswith("Arial") for n in names),
          str(sorted(names)))

    # each panel PDF carries its own letter (the composed page merges the glyph
    # into a neighbouring text run during pypdf composition)
    page = text_of(PDF)
    got_letters = []
    for key in KEYS:
        piece = text_of(S.PANELS / f"figS9{key}.pdf")
        found = [t for t in piece.split() if len(t) == 1 and t.isupper()]
        got_letters.append(found[0] if found else "?")
    check("panel letters A-I once each", got_letters == KEYS, str(got_letters))

    
    #: anchored to its own canvas at ``S.LETTER_IN``, so the three letters of a
    #: column share one x -- and that x is the canvas' left edge plus that
    #: constant.  The composed page carries the glyphs as text, but pdftotext
    #: splits the rotated tool labels into 8-9 pt fragments, so a panel letter is
    #: a lone A-I whose box is taller than 25 pt (the 28 pt glyph measures 31.3).
    _pos = {t: (x0, y0) for x0, y0, x1, y1, t in word_boxes(PDF)
            if len(t) == 1 and t in KEYS and (y1 - y0) > 25.0}
    check("nine panel letters visible in the page text layer",
          sorted(_pos) == sorted(KEYS), f"found {sorted(_pos)}")
    if sorted(_pos) == sorted(KEYS):
        _off = {k: _pos[k][0] - (S.PLACE[k][0] + S.LETTER_IN) * 72.0 for k in KEYS}
        _worst = max(_off.items(), key=lambda kv: abs(kv[1]))
        check("panel letters sit at their canvas' left edge + LETTER_IN",
              abs(_worst[1]) <= 0.5, f"worst {_worst[0]} {_worst[1]:+.2f} pt (tol 0.5)")
        _col = max(max(_pos[k][0] for k in g) - min(_pos[k][0] for k in g)
                   for g in ("ADF", "BEG", "CHI"))
        check("letters of one column share one x", _col <= 0.5,
              f"spread {_col:.2f} pt (tol 0.5)")
        _row = max(max(_pos[k][1] for k in g) - min(_pos[k][1] for k in g)
                   for g in ("ABC", "DEH", "FGI"))
        check("letters of one row share one y", _row <= 0.5,
              f"spread {_row:.2f} pt (tol 0.5)")

    banned = ("hit rate", "near-perfect", "jaccard", "percent_modified",
              "rna002", "chemistry", "nested subset", "unit mean",
              "3'utr", "5'utr", "cross-replicate", "dumbbell")
    bad = [w for w in banned if w in page.lower()]
    check("no retired wording, no RNA002 on the page", not bad, str(bad))
    check("no main-figure duplicate", "PPV vs. GLORI (2 bp)" not in page)
    check("false-positive-rate axis on the page", "FP rate" in page,
          "the C top axis prints the false-positive rate in scientific notation"
          if "FP rate" in page else "axis label missing")
    check("coverage-normalised axis on C", "candidates" in page)
    check("count facets titled by their key only",
          "ORCA" in page and page.count("WT") >= 2 and page.count("IVT") >= 2)
    check("no unused 002 tool name on the page",
          all(t not in page for t in ("DENA", "EpiNano_Error", "MINES",
                                      "Nanocompore", "Nanom6A", "xPore",
                                      "Yanocomp", "CHEUI")))
    
    #: by name (a marker is a drawn entry, so the house rule wants it named), and
    #: each key prints **under its own panel**, so every species name appears twice
    #: on the page -- once under H and once under I.  What has to stay out is a
    #: *sample* name printed inside a panel.
    #: 2026-09-30: the marker key prints the species names of the 2026-09-28
    
    #: Arabidopsis / mouse / Human -- and a stray "HeLa" on the sheet is itself a
    #: failure now (the caption keeps its own wording; this check reads the figure).
    check("no panel title, species named only by the two marker keys",
          "Curlcake" not in page and "HeLa" not in page
          and all(page.count(s) == 2 for s in ("Arabidopsis", "mouse", "Human")),
          f"A {page.count('Arabidopsis')} M {page.count('mouse')} "
          f"H {page.count('Human')}")

    sh = D.panel_stability("H")
    si = D.panel_stability("I")
    check("anchor coverage-stability cells",
          len(sh) == 12 and abs(sh["rho"].min() - 0.801) < 1e-3
          and abs(sh["rho"].max() - 1.000) < 1e-3,
          f"{len(sh)} cells, rho {sh['rho'].min():.3f}-{sh['rho'].max():.3f}")
    check("anchor ratio-stability cells",
          len(si) == 9 and abs(si["rho"].min() - 0.478) < 1e-3
          and abs(si["rho"].max() - 0.904) < 1e-3
          and int((~si["top1_same"]).sum()) == 6,
          f"{len(si)} cells, rho {si['rho'].min():.3f}-{si['rho'].max():.3f}, "
          f"top-ranked tool differs in {int((~si['top1_same']).sum())} of "
          f"{len(si)} cells")
    check("no anonymous entry group on the page",
          "other entries" not in page)
    sweep = D.panel_window()
    # the printed names come from the panel's own palette (one source of truth)
    names = [P.D_ENTRY[str(t)][0] for t in sweep["tool"].unique()]
    check("every sweep entry named", len(names) == 13
          and all(n in page for n in names), f"{len(names)} entries")
    cols = [c for _t, _n, c in P.window_key(sweep)]
    check("thirteen distinct sweep colours", len(set(cols)) == 13,
          f"{len(set(cols))} colours")
    
    #: note were dropped; the axis now reads plain "Slope"
    check("effect-size panel axis", "Slope" in page and "1.0 = proportional" not in page)
    check("window panel axes", "Matching window" in page
          and "Exact fraction" in page and "PPV" in page)
    check("effect-size panel axes", "Pearson r" in page and "Slope" in page)
    check("fitted-lines panel retired (F is the effect size now)",
          "GLORI ratio (%)" not in page and "Predicted ratio" not in page)
    #: 2026-09-30: the "m6A (non-DRACH)" family key of the first drafts went away
    
    
    #: any built sheet -- old or new.  What the page must still name is the four
    #: modified families its rows belong to, one word each.
    check("family key on the page",
          all(w in page for w in ("m6A DRACH", "inosine", "pseU", "m5C")))

    # ---- drawing-code house rules -----------------------------------------
    src = "\n".join(p.read_text() for p in sorted(HERE.glob("s9_*.py")))
    check("no grid in the drawing code",
          ".grid(" not in src and "axes.grid': True" not in src)
    check("no RNA002 branch left in the data layer",
          '"RNA002"' not in src and "'RNA002'" not in src)
    check("no scatter clouds", "scatter(" not in src
          and "MAX_DOTS_PER_CALL" in src)
    log = (S.LOGS / "70_panels.log").read_text()
    m = re.search(r"dots drawn in total: (\d+)", log)
    dots = int(m.group(1)) if m else -1
    check("dot count within the guard", 0 < dots < 400, f"dots in 70 log: {dots}")

    # ---- frozen-table anchors ---------------------------------------------
    a = D.panel_a()
    a_idx = a.set_index("label")
    check("anchor m6Anet HeLa WT sites",
          float(a_idx.loc["m6Anet", "WT"]) == 30327.0)
    check("anchor count matrix rows", len(a) == 23, str(len(a)))
    dr = a[a["label"].str.contains("DRACH")]
    check("anchor DRACH IVT range",
          (int(dr["IVT"].min()), int(dr["IVT"].max())) == (2, 10),
          f"{int(dr['IVT'].min())}-{int(dr['IVT'].max())}")

    b = D.panel_orca()
    check("anchor ORCA channels", len(b) == 8, str(len(b)))
    check("anchor ORCA WT range",
          (int(b["WT"].min()), int(b["WT"].max())) == (39, 14784))

    c = D.panel_curlcake()
    c_idx = c.set_index("label")
    dr_c = c[c["block"] == "DRACH"]
    nd_c = c[c["block"] == "non-DRACH"]
    check("anchor Curlcake entries", len(c) == 10, str(len(c)))
    check("anchor Curlcake DRACH calls",
          (int(dr_c["calls"].min()), int(dr_c["calls"].max())) == (0, 4))
    check("anchor Curlcake non-DRACH calls",
          (int(nd_c["calls"].min()), int(nd_c["calls"].max())) == (10, 13))
    check("anchor Curlcake m6Anet rate",
          abs(float(c_idx.loc["m6Anet", "per_10kb"]) - 4.9334) < 1e-3)
    check("anchor Curlcake NanoSPA_m6A rate",
          abs(float(c_idx.loc["NanoSPA_m6A", "per_10kb"]) - 16.7736) < 1e-3)
    check("anchor specificity range",
          c["specificity"].min() > 0.99 and c["specificity"].max() == 1.0,
          f"{c['specificity'].min() * 100:.2f}-{c['specificity'].max() * 100:.2f} %")

    e = D.panel_effect()
    e_idx = e.set_index("label")
    check("anchor nine ratio tools", len(e) == 9, str(len(e)))
    check("anchor m6Anet r", abs(float(e_idx.loc["m6Anet", "r"]) - 0.6526) < 1e-3)
    check("anchor m6Anet slope",
          abs(float(e_idx.loc["m6Anet", "slope"]) - 0.6115) < 1e-3)

    check("anchor nine effect sizes", len(e) == 9, str(len(e)))

    print(f"\nfailures: {len(fails)}" + (f" {fails}" if fails else " (none)"))
    (S.LOGS / "72_verify.log").write_text("\n".join(fails) or "ALL CHECKS PASSED\n")
    if fails:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
