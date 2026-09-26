#!/usr/bin/env python
"""Clipping audit for the redrawn Figure S1 (read-only inspection).

Checks every standalone panel plus the assembled page: no text span and no
vector-drawing rectangle may extend beyond the page box.  Writes
figures/figureS1/tables/figS1_clipping_check.tsv; exits non-zero if any file
fails, so the render loop can bump margins and retry.

Usage:
  conda run -n viz python 23d_check_clipping.py
"""

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
import glob
import os
import sys

import pymupdf

ROOT = str(_RB / "figures/figureS1")
OUT = os.path.join(ROOT, "tables", "figS1_clipping_check.tsv")
TOL = 0.6  # pt of slack for hairline rounding

targets = sorted(glob.glob(os.path.join(ROOT, "figures", "panels", "*.pdf")))
targets.append(os.path.join(ROOT, "figures", "FigureS1_rev.pdf"))

rows, n_fail = [], 0
for f in targets:
    page = pymupdf.open(f)[0]
    R = page.rect
    text_out, draw_out, min_sz, examples = [], [], 1e9, []
    for b in page.get_text("dict")["blocks"]:
        if b["type"]:
            continue
        for line in b.get("lines", []):
            for s in line["spans"]:
                t = s["text"].strip()
                if not t:
                    continue
                x0, y0, x1, y1 = s["bbox"]
                min_sz = min(min_sz, round(s["size"], 2))
                if x0 < R.x0 - TOL or x1 > R.x1 + TOL or \
                   y0 < R.y0 - TOL or y1 > R.y1 + TOL:
                    text_out.append(t)
                    if len(examples) < 3:
                        examples.append(f"{t}@({x0:.0f},{y0:.0f},{x1:.0f},{y1:.0f})")
    for d in page.get_drawings():
        r = d["rect"]
        if r.x1 > R.x1 + TOL or r.y1 > R.y1 + TOL or \
           r.x0 < R.x0 - TOL or r.y0 < R.y0 - TOL:
            draw_out.append(f"({r.x0:.0f},{r.y0:.0f},{r.x1:.0f},{r.y1:.0f})")
    ok = not text_out and not draw_out
    n_fail += 0 if ok else 1
    rows.append((os.path.relpath(f, ROOT), f"{R.width:.1f}", f"{R.height:.1f}",
                 len(text_out), len(draw_out), min_sz,
                 "PASS" if ok else "FAIL", " | ".join(examples or draw_out[:3])))

with open(OUT, "w") as fh:
    fh.write("file\tpage_w_pt\tpage_h_pt\tn_text_out\tn_draw_out\tmin_font_pt\t"
             "status\texamples\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")

print(f"checked {len(rows)} files -> {OUT}")
for r in rows:
    if r[6] == "FAIL":
        print("FAIL:", r)
print(f"status: {'ALL PASS' if n_fail == 0 else f'{n_fail} FAILED'}")
sys.exit(0 if n_fail == 0 else 1)
