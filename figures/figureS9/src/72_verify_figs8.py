#!/usr/bin/env python
"""72 -- acceptance checks for the rebuilt Figure S8 (run after 70/71).

Contract of the 2026-09-25 page: A4 **landscape** 842.4 x 595.44 pt, three rows
-- A | B | C, then the full-width window sweep D (its key on the right), then
the merged effect-size panel E -- **RNA004 data only** (no RNA002 run is drawn
anywhere), every panel's numbers reproduce the frozen tables, and the drawing
code keeps the house rules (no grid, no scatter clouds, Arial only, nothing
below 7 pt).
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import s8_data as D  # noqa: E402
import s8_panels as P  # noqa: E402
import s8_style as S  # noqa: E402

PDF = S.FIGS / "FigureS8_rev.pdf"
PAGE_PT = (842.4, 595.44)                       # A4 landscape since 2026-09-25
KEYS = list("ABCDEFG")                          # 2026-09-25: the sweep and the
                                                # effect sizes are two facets each
fails: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"[{'ok' if ok else 'FAIL'}] {name}{' ' + detail if detail else ''}")
    if not ok:
        fails.append(name)


def text_of(path: Path) -> str:
    return subprocess.run(["/usr/bin/pdftotext", str(path), "-"],
                          capture_output=True, text=True).stdout


def main() -> None:
    info = subprocess.run(["/usr/bin/pdfinfo", str(PDF)], capture_output=True,
                          text=True).stdout
    size = [ln for ln in info.splitlines() if ln.startswith("Page size")][0]
    got = tuple(float(x) for x in re.findall(r"[\d.]+", size)[:2])
    check("page size", abs(got[0] - PAGE_PT[0]) < 0.01
          and abs(got[1] - PAGE_PT[1]) < 0.01, f"{got}")
    check("single page", "Pages:           1" in info)

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
        piece = text_of(S.PANELS / f"figS8{key}.pdf")
        found = [t for t in piece.split() if len(t) == 1 and t.isupper()]
        got_letters.append(found[0] if found else "?")
    check("panel letters A-G once each", got_letters == KEYS, str(got_letters))

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
    check("no panel title (names live in the legend)",
          "HeLa" not in page and "Curlcake" not in page)
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
    check("family key on the page",
          "m6A DRACH" in page and "inosine" in page
          and "m6A (non-DRACH)" in page)

    # ---- drawing-code house rules -----------------------------------------
    src = "\n".join(p.read_text() for p in sorted(HERE.glob("s8_*.py")))
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
    dr_c = c[c["block"] == "m6A DRACH"]
    nd_c = c[c["block"] == "m6A (non-DRACH)"]
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
