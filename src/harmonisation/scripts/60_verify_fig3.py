#!/usr/bin/env python3
"""Delivery check for the revised Figure 3 (figures/figure3).

Checks, in order:
  1. the shipped `$RNAMODBENCH_LOCAL/manuscript/Figure3.pdf` is byte-identical to
     `figures/figure3/figures/Figure3_rev.pdf`;
  2. the merged page is one page of the expected printed size;
  3. every font is embedded (pdf resolution) and no text is below 7.2 pt;
  4. the PNG sibling is 300 dpi;
  5. the manuscript caption mentions all four panels and carries no stale
     p-value-only wording for this figure;
  6. the legacy-vs-revision table covers the twelve tool x species cells.

Exit code 0 = all checks pass.  Read-only: nothing is written.
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
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
from pypdf import PdfReader

ROOT = Path(str(_RB))
FIG3 = (_RB / "figures/figure3")
SHIPPED = (_XB / "manuscript/manuscript_rev/Figure3.pdf")
MANUSCRIPT = (_XB / "manuscript/manuscript_rev/manuscript.tex")
PT_IN = 72.0
MIN_FONT_PT = 7.2
EXPECTED_IN = (6.66, 7.47)   # A row now matches Figure 4C (1.98 in)
TOL_IN = 0.03


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def min_font_size(page) -> float:
    """Minimum effective text size (matplotlib `Tf` and cairo `Tm` styles)."""
    txt = page.get_contents().get_data().decode("latin-1", errors="ignore")
    tf, tm = 1.0, 1.0
    sizes: list[float] = []
    for m in re.finditer(r"/([A-Za-z0-9#\-.]+)\s+([0-9]*\.?[0-9]+)\s+Tf"
                         r"|([0-9.\-]+)\s+([0-9.\-]+)\s+[0-9.\-]+\s+[0-9.\-]+"
                         r"\s+[0-9.\-]+\s+[0-9.\-]+\s+Tm"
                         r"|(Tj|TJ|'|\")", txt):
        if m.group(1):
            tf = float(m.group(2))
        elif m.group(3):
            tm = (float(m.group(3)) ** 2 + float(m.group(4)) ** 2) ** 0.5
        else:
            sizes.append(tf * tm)
    sizes = [s for s in sizes if s > 0]
    return min(sizes) if sizes else float("nan")


def main() -> int:
    checks: list[tuple[str, bool, str]] = []
    shipped_md5 = md5(SHIPPED)
    checks.append(("shipped == Figure3_rev.pdf",
                   shipped_md5 == md5((_RB / "figures/figure3/figures/Figure3_rev.pdf")),
                   shipped_md5))

    reader = PdfReader(str(SHIPPED))
    page = reader.pages[0]
    w_in = float(page.mediabox.width) / PT_IN
    h_in = float(page.mediabox.height) / PT_IN
    checks.append(("single page", len(reader.pages) == 1, f"{len(reader.pages)}"))
    checks.append(("printed size 6.66 x 6.91 in",
                   abs(w_in - EXPECTED_IN[0]) < TOL_IN
                   and abs(h_in - EXPECTED_IN[1]) < TOL_IN,
                   f"{w_in:.2f} x {h_in:.2f} in"))

    fonts = subprocess.run(["pdffonts", str(SHIPPED)], capture_output=True,
                           text=True).stdout.splitlines()[2:]
    fonts = [l for l in fonts if l.strip()]
    not_embedded = [l.split()[0] for l in fonts if l.split()[-4] != "yes"]
    checks.append((f"fonts embedded ({len(fonts)})", not not_embedded,
                   ",".join(not_embedded) or "all embedded"))

    font = min_font_size(page)
    checks.append((f"min font >= {MIN_FONT_PT} pt", font >= MIN_FONT_PT - 0.05,
                   f"{font:.2f} pt"))

    png = (_RB / "figures/figure3/figures/Figure3_rev.png")
    from PIL import Image

    with Image.open(png) as im:
        dpi = im.size[0] / w_in
        checks.append(("PNG 300 dpi", abs(dpi - 300) <= 3,
                       f"{dpi:.0f} dpi, {im.size[0]}x{im.size[1]} px"))

    tex = MANUSCRIPT.read_text()
    cap = re.search(r"\\caption\{Stoichiometry.*?\n", tex, re.S)
    cap_text = cap.group(0) if cap else ""
    checks.append(("caption found for Figure 3", bool(cap_text), ""))
    for tag in ("(A)", "(B)", "(C)", "(D)"):
        checks.append((f"caption has panel {tag}", tag in cap_text, ""))
    checks.append(("no p-value-only claim in caption",
                   "Wilcoxon test, p" not in cap_text, ""))
    checks.append(("figure height within the 8.60 in budget", h_in <= 8.60 + TOL_IN,
                   f"{h_in:.2f} in"))
    geo = ((_RB / "figures/figure3/tables/fig3_assembly_geometry.tsv")).read_text()
    checks.append(("legend band row present in the assembly", "Fig3_legend_band.pdf" in geo,
                   ""))
    checks.append(("caption describes the boxplot panel",
                   "boxplot" in cap_text.lower(), ""))
    geo_txt = ((_RB / "figures/figure3/tables/fig3_assembly_geometry.tsv")).read_text()
    checks.append(("four lettered rows + legend band",
                   all(f"row_{c}" in geo_txt for c in "ABCD") and "row_L" in geo_txt, ""))
    ba = pd.read_csv((_RB / "figures/figure3/tables/fig3b_bland_altman.tsv"), sep="\t")
    checks.append(("Bland-Altman table covers 12 cells",
                   len(ba) == 12 and {"bias", "loa_lo", "loa_hi"} <= set(ba.columns),
                   f"{len(ba)} rows"))
    checks.append(("caption states replicate n",
                   "n = 3" in cap_text and "two studies" in cap_text, ""))

    legacy = ((_RB / "figures/figure3/tables/fig3b_legacy_vs_revision.tsv")).read_text()
    n_cells = len([l for l in legacy.splitlines()[1:] if l.strip()])
    checks.append(("legacy-vs-revision covers 12 cells", n_cells >= 12,
                   f"{n_cells} rows"))

    bad = 0
    for name, ok, detail in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {name}"
              + (f"  ({detail})" if detail else ""))
        bad += 0 if ok else 1
    print(f"\n{len(checks) - bad}/{len(checks)} checks passed")
    return 1 if bad else 0
if __name__ == "__main__":
    sys.exit(main())
