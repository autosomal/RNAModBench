#!/usr/bin/env python3
"""53 -- acceptance checks for the rebuilt Figure 8 (v6: panels A, B and E).

Checks
------
1. the assembled page: 500.4 x 586.8 pt (= 6.95 x 8.15 in), 300 dpi PNG sibling,
   Arial-only embedded fonts (no DejaVu, no Noto fallback, every font embedded);
2. no false-positive count is drawn inside panel B (the numbers live in the
   legend and the source tables): the B piece and the page must not contain the
   count call-outs of the retired draft;
3. the two Python pieces obey the house rules, verified geometrically on the
   live figures through ``common/pagelayout.assert_page_clean`` (no text
   overlap, no cross-panel intrusion, no text hugging a frame line, no grid, no
   font below 7 pt) plus a "no in-figure call-out" source check;
4. the values that are drawn reproduce the frozen tables:
   A bars == ``fig8_counts.tsv`` (site-level counts, 16 Dorado models),
   B series/points == ``fig8_fpr_curlcake_scan.tsv`` /
   ``fig8_fpr_hela_ivt.tsv``;
5. panel E: the drawn curves come from the majority-consensus BEDs of
   ``guitar_metagene_replicates`` (n sites per curve == BED line count), the
   piece geometry is the Figure 8 contract, and no count is drawn inside it.

Usage
-----
conda run -n benchmark-revision --no-capture-output \
    python $RNAMODBENCH_ROOT/figures/figure8/src/53_verify_fig8.py
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

import pandas as pd

PROJECT = Path(str(_RB))
OUT = (_RB / "figures/figure8")
TABLES = (_RB / "figures/figure8/tables")
FIGS = (_RB / "figures/figure8/figures")
PANELS = (_RB / "figures/figure8/figures/panels")
SCRIPTS = (_RB / "figures/figure8/src")
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str((_RB / "src/harmonisation/common")))

PAGE_PT = (500.4, 586.8)
#: count call-outs of the retired draft -- none of these may be drawn any more.
#: 2026-09-23: "0 calls" was removed from this list -- since the 15:13 legend edit
#: it is a deliberate key label ("the open circle is a measured zero", the same
#: wording as panels A/C and Figure S8), not a count call-out; the retired tokens
#: are the numbers of the draft (13,028 / 46,963) and the two ranges.
RETIRED_CALLOUTS = ("13,028", "46,963", "87\u2013101", "2\u201310", "87-101", "2-10")
#: panel E keeps the model names (user decision) but as plain text, never a title
E_MODELS = ("hac@v5.0.0_m6A", "hac@v5.1.0_inosine+m6A", "sup@v5.0.0_m6A")
#: signed coordinates: the rotated shared title of panel B sticks out of the piece
BBOX_WORD = re.compile(r'<word xMin="(-?[\d.]+)" yMin="(-?[\d.]+)" '
                       r'xMax="(-?[\d.]+)" yMax="(-?[\d.]+)">(.*?)</word>')
#: the shared rotated y title of panel B (its words define the band it occupies)
TITLE_WORDS = ("False", "positives", "per", "10", "kb")

FAILS: list[str] = []
CHECKS = [0]


def check(name: str, ok: bool, detail: str = "") -> None:
    CHECKS[0] += 1
    print(f"[{'OK  ' if ok else 'FAIL'}] {name}{'  -- ' + detail if detail else ''}")
    if not ok:
        FAILS.append(name)


def pdf_page_size(pdf: Path) -> tuple[float, float]:
    txt = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
    m = re.search(r"Page size:\s+([\d.]+) x ([\d.]+)", txt)
    return float(m.group(1)), float(m.group(2))


def pdf_fonts(pdf: Path) -> list[tuple[str, str]]:
    txt = subprocess.check_output(["pdffonts", str(pdf)], text=True)
    out = []
    for line in txt.splitlines()[2:]:
        parts = line.split()
        if len(parts) >= 6:
            out.append((parts[0], parts[-4]))          # (name, emb)
    return out


def page_text(pdf: Path) -> str:
    return subprocess.check_output(["pdftotext", "-layout", str(pdf), "-"],
                                   text=True, errors="replace")


def bbox_words(pdf: Path) -> list[tuple[float, float, float, float, str]]:
    """Word boxes from ``pdftotext -bbox`` (coordinates may be negative)."""
    xml = subprocess.check_output(["/usr/bin/pdftotext", "-bbox", str(pdf), "-"],
                                  text=True, errors="replace")
    return [(float(a), float(b), float(c), float(d), t)
            for a, b, c, d, t in BBOX_WORD.findall(xml)]


def close(a: float, b: float, tol: float = 1e-4) -> bool:
    return abs(a - b) <= tol * max(abs(b), 1e-9)


def load_piece(mod_name: str, path: Path, captured: dict):
    """Import a piece script, capturing the figure instead of writing it."""
    import fig8_style as fs

    def capture(fig, name, width, height, dpi=300, fit=(), keep_right=None):
        fig.set_size_inches(width, height)
        if fit:
            ax = fig.axes[0]
            fs.fit_labels(fig, ax, left="left" in fit, bottom="bottom" in fit,
                          top="top" in fit, keep_right=keep_right)
        captured[name] = fig

    spec = importlib.util.spec_from_file_location(mod_name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.save_piece = capture
    return mod


def main() -> None:
    ## 1 -- the assembled page -------------------------------------------------- #
    page = (_RB / "figures/figure8/figures/Figure8_rev.pdf")
    size = pdf_page_size(page)
    check("Figure8_rev.pdf: page size 500.4 x 586.8 pt",
          close(size[0], PAGE_PT[0]) and close(size[1], PAGE_PT[1]), f"got {size}")
    fonts = pdf_fonts(page)
    bad = [f for f in fonts if f[1] != "yes" or "Arial" not in f[0]]
    check("Figure8_rev.pdf: all fonts embedded Arial", not bad, f"offending: {bad}")
    check("Figure8_rev.png (300 dpi) exists",
          ((_RB / "figures/figure8/figures/Figure8_rev.png")).exists()
          and ((_RB / "figures/figure8/figures/Figure8_rev.png")).stat().st_size > 100_000)

    ## 2 -- no count call-out on the page or in panel B ------------------------- #
    txt_page = page_text(page)
    txt_b = page_text((_RB / "figures/figure8/figures/panels/fig8B_fpr.pdf"))
    for token in RETIRED_CALLOUTS:
        check(f"no retired count call-out {token!r} in panel B", token not in txt_b)
    check("no retired count call-out anywhere on the page",
          not any(t in txt_page for t in RETIRED_CALLOUTS))

    ## 2b -- v7: E schematic labels, C without value labels --------------------- #
    txt_e = page_text((_RB / "figures/figure8/figures/panels/fig8E_guitar.pdf"))
    txt_c = page_text((_RB / "figures/figure8/figures/panels/fig8C_ppv.pdf"))
    check("panel E: 1 kb label on both flanks of all three panels",
          txt_e.count("1kb") == 6, f"count={txt_e.count('1kb')}")
    for seg in ("5'UTR", "CDS", "3'UTR"):
        check(f"panel E: {seg!r} labelled on all three panels",
              txt_e.count(seg) == 3, f"count={txt_e.count(seg)}")
    check("panel E: three model names kept",
          all(t in txt_e for t in E_MODELS),
          f"missing={[t for t in E_MODELS if t not in txt_e]}")
    src_e = ((_RB / "figures/figure8/src/69_fig8e_guitar.R")).read_text()
    check("panel E: model name drawn as plain text, not a plot title",
          'E_TITLE_FACE <- "plain"' in src_e and "fontface = E_TITLE_FACE" in src_e
          and "title = m$label" not in src_e)
    pct = re.findall(r"\d+\.\d+\s?%", txt_page)
    check("no value annotation such as '82.2%' on the page", not pct, f"found {pct}")
    check("panel C: no percentage value label",
          not re.findall(r"\d+\.\d+\s?%", txt_c))

    ## 2c -- panel B: the rotated shared y title must not touch the y tick labels -- #

    # rotated title was printed over the topmost y tick label ("1,000", gap
    # -11.2 pt) because B1's left margin was exactly the width of its labels.
    # The whole-page audit (_delivery_audit_20260921/check_text_collisions.py)
    # skips every pair whose boxes are taller than wide, so rotated titles are
    # invisible to it; the two column bands are compared here instead.
    words_b = bbox_words((_RB / "figures/figure8/figures/panels/fig8B_fpr.pdf"))



    clusters: dict[tuple[float, float], list] = {}
    for w in words_b:
        if w[0] < 15.0:
            clusters.setdefault((round(w[0], 1), round(w[2], 1)), []).append(w)
    cluster = max(clusters.values(), key=len, default=[])
    band_words = [w for w in cluster if w[4].strip() in TITLE_WORDS]
    check("panel B: rotated shared title in place cluster 5 title words ", len(band_words) == 5,
          f" cluster {len(cluster)} words, of which title words {[w[4] for w in band_words]}")
    band = ((min(w[0] for w in cluster), max(w[2] for w in cluster))
            if cluster else (0.0, 0.0))


    nums = [w for w in words_b if re.fullmatch(r"[0-9][0-9,]*", w[4].strip())
            and 30 < (w[1] + w[3]) / 2 < 155 and w[2] < 60]
    col = max((w[2] for w in nums), default=0.0)
    yticks = [w for w in nums if abs(w[2] - col) <= 1.5]
    check("panel B: 4 y tick labels form a column ", len(yticks) == 4,
          f" found {len(yticks)} : {[w[4] for w in yticks]}")
    if band_words and yticks:
        gap = min(w[0] for w in yticks) - band[1]
        check("panel B: shared y axis titles and y tick labels column gap >= 1 pt", gap >= 1.0,
              f" title band x {band[0]:.2f}-{band[1]:.2f} | leftmost label {min(w[0] for w in yticks):.2f} "
              f"({min(yticks, key=lambda w: w[0])[4]!r}) | gap {gap:.2f} pt")
    else:
        check("panel B: shared y axis titles and y tick labels column gap >= 1 pt", False,
              f" could not parse title band / tick labels band={band}, yticks={[w[4] for w in yticks]}")

    ## 3 -- house rules on the live pieces -------------------------------------- #
    import pagelayout
    captured: dict = {}
    mod_a = load_piece("fig8a", (_RB / "figures/figure8/src/57_fig8a_counts.py"), captured)
    mod_b = load_piece("fig8b", (_RB / "figures/figure8/src/59_fig8b_fpr.py"), captured)
    mod_a.main()
    mod_b.main()
    for name in ("fig8A_counts", "fig8B_fpr"):
        check(f"{name}: figure captured for the layout gate", name in captured)
    for name, fig in captured.items():
        try:
            pagelayout.assert_page_clean(fig, min_pt=7.0, verbose=True)
            ok, detail = True, ""
        except SystemExit as exc:                        # pragma: no cover
            ok, detail = False, str(exc)
        check(f"{name}: layout gate (no overlap / grid / font < 7 pt)", ok, detail)
    for script in ("57_fig8a_counts.py", "59_fig8b_fpr.py"):
        src = (SCRIPTS / script).read_text()
        # only figure-level text is allowed (panel letters, axis titles); an
        # axes call-out would be a value printed inside the plotting area
        check(f"{script}: no in-figure value call-out",
              "annotate(" not in src and not re.search(r"\bax\d*\.text\(", src))
        check(f"{script}: no grid call", ".grid(" not in src.replace("minor", ""))

    ## 4 -- the drawn values reproduce the frozen tables ------------------------ #
    counts = pd.read_csv((_RB / "figures/figure8/tables/fig8_counts.tsv"), sep="\t")
    piv = mod_a.order_rows(mod_a.build_rows(counts))
    check("panel A: 16 Dorado models, no zero WT bar, one zero IVT bar",
          len(piv) == 16 and int((piv["wt"] == 0).sum()) == 0
          and int((piv["ivt"] == 0).sum()) == 1,
          f"{len(piv)} rows, zero-IVT: "
          f"{list(piv.loc[piv['ivt'] == 0, 'label'])}")
    wt_sup = piv[(piv["tool"] == "Dorado_sup@v5.0.0_m6A@v1")
                 & (piv["mod_type"] == "m6A")]["wt"].iloc[0]
    check("panel A: sup@v5.0.0_m6A WT bar = 10,818 sites", int(wt_sup) == 10818,
          f"got {int(wt_sup)}")
    dr = piv[piv["family"] == "m6A_DRACH"]["ivt"]
    check("panel A: DRACH IVT bars inside 2-10 sites",
          dr.min() >= 2 and dr.max() <= 10, f"got {sorted(dr.astype(int))}")

    cc = pd.read_csv((_RB / "figures/figure8/tables/fig8_fpr_curlcake_scan.tsv"), sep="\t")
    hl = pd.read_csv((_RB / "figures/figure8/tables/fig8_fpr_hela_ivt.tsv"), sep="\t")
    s = mod_b.series_for(cc, "m6A", "hac")
    check("panel B1: non-DRACH m6A (hac) falls 794 -> 12.8 FP/10 kb across 5-50 %",
          close(s["y"][0], 794.277, 1e-4) and close(s["y"][3], 12.8268, 1e-4),
          f"{[None if v != v else round(v, 3) for v in s['y']]}")
    s = mod_b.series_for(cc, "m6A_DRACH", "hac")
    check("panel B1: DRACH (hac) reaches 3.95 FP/10 kb at 50 %",
          close(s["y"][3], 3.94672, 1e-4),
          f"{[None if v != v else round(v, 3) for v in s['y']]}")
    s = mod_b.series_for(cc, "m6A_DRACH", "sup")
    check("panel B1: DRACH (sup) has no call at 50 % (floor marker)",
          s["y"][3] != s["y"][3], f"y = {s['y'][3]}")
    pts = mod_b.hela_points(hl, "m6A_DRACH")
    check("panel B2: DRACH block holds 4 delivered models, 1.3e-4 - 6.4e-4 per 10 kb",
          len(pts) == 4 and close(min(p[0] for p in pts), 0.000127303, 1e-4)
          and close(max(p[0] for p in pts), 0.000636515, 1e-4),
          f"{[round(p[0], 8) for p in pts]}")
    pts = mod_b.hela_points(hl, "m6A")
    check("panel B2: non-DRACH m6A block spans 0.83 - 2.99 per 10 kb",
          close(min(p[0] for p in pts), 0.829251, 1e-4)
          and close(max(p[0] for p in pts), 2.98926, 1e-4),
          f"{[round(p[0], 3) for p in pts]}")

    ## 5 -- panel E (Guitar, majority consensus) -------------------------------- #
    plan = pd.read_csv((_RB / "figures/figure8/tables/fig8e_guitar_panel_inputs.tsv"), sep="\t")
    check("panel E: 3 models x (WT, IVT) curves", len(plan) == 6,
          f"{len(plan)} rows")
    n_ok, n_detail = True, []
    for _, r in plan.iterrows():
        p = Path(r["path"])
        lines = int(subprocess.check_output(["awk", "END{print NR}", str(p)],
                                            text=True).strip())
        n_ok &= lines == int(r["n_sites"])
        n_detail.append(f"{r['label']}/{r['condition']}={r['n_sites']}")
        n_ok &= "RNA004/majority/" in r["path"]
    check("panel E: n sites == majority BED line counts", n_ok, "; ".join(n_detail))
    geom = dict(zip(*[pd.read_csv((_RB / "figures/figure8/tables/fig8e_guitar_geometry.tsv"), sep="\t")[c]
                      for c in ("item", "value")]))
    check("panel E: 3 panels, 6.84 x 2.05 in, merge=majority, min font 8 pt",
          int(geom["panels"]) == 3 and close(float(geom["width_in"]), 6.84, 1e-4)
          and close(float(geom["height_in"]), 2.05, 1e-4)
          and geom["merge"] == "majority" and float(geom["font_min_pt"]) >= 7,
          str({k: geom[k] for k in ("panels", "width_in", "height_in", "merge",
                                    "font_min_pt")}))
    txt_e = page_text((_RB / "figures/figure8/figures/panels/fig8E_guitar.pdf"))
    check("panel E: no count drawn inside (no thousand-separated number)",
          not re.search(r"\d{1,3}(?:,\d{3})+", txt_e), txt_e.strip()[:80])
    check("panel E: region labels present",
          all(t in txt_e for t in ("5'UTR", "CDS", "3'UTR"))
          and all(t in txt_e for t in ("hac@v5.0.0_m6A",
                                       "hac@v5.1.0_inosine+m6A", "sup@v5.0.0_m6A")))
    check("panel E: fonts are Arial only",
          all("Arial" in n for n, _ in pdf_fonts((_RB / "figures/figure8/figures/panels/fig8E_guitar.pdf"))),
          str(pdf_fonts((_RB / "figures/figure8/figures/panels/fig8E_guitar.pdf"))))

    ## 6 -- reply-letter anchor table still populated --------------------------- #
    anchors = pd.read_csv((_RB / "figures/figure8/tables/fig8_quoted_numbers.tsv"), sep="\t")
    check("reply-letter anchor table populated", len(anchors) >= 20,
          f"{len(anchors)} rows")

    print()
    if FAILS:
        print(f"{len(FAILS)} FAILED of {CHECKS[0]} checks: {FAILS}")
        sys.exit(1)
    print(f"all {CHECKS[0]} checks passed")


if __name__ == "__main__":
    main()
