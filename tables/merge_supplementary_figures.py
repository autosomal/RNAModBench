#!/usr/bin/env python3
"""Merge FigureS1-S10 into one Supplementary_Figures.pdf with the caption under each figure.

The previously submitted `Supplementary_Figures.pdf` carried a caption beneath every
figure (e.g. "Figure S1. Metagene plots revealing ..."), so the rebuilt file does too.
Captions are read from the delivered legends (`figures/figureS*/
figures/*_legends.md`); S3 has no legend file yet, so its caption is written here.

Usage: conda run -n benchmark-revision --no-capture-output python merge_supplementary_figures.py
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
import re
import subprocess
from pathlib import Path

RA = Path(str(_RB / "analysis"))
NEW = Path(str(_XB / "submission/new_submission"))
LEG = {
    1: (_RB / "figures/figureS1/figures/FigS1_legends.md"),
    2: (_RB / "figures/figureS2/figures/FigS2_legends.md"),
    4: (_RB / "figures/figureS4/figures/FigS4_legends.md"),
    5: (_RB / "figures/figureS5/figures/FigS5_legends.md"),
    6: (_RB / "figures/figureS6/figures/FigS6_legends.md"),
    7: (_RB / "figures/figureS7/figures/FigS7_legends.md"),
    8: (_RB / "figures/figureS8/figures/FigS8_legends.md"),
    9: (_RB / "figures/figureS9/figures/FigS9_legends.md"),
    10: (_RB / "figures/figureS10/figures/FigS10_legends.md"),
}
S3 = (
    "**Figure S3.** GLORI reference construction, tool performance against it and "
    "modification-ratio agreement. (A) Overlap between the two GLORI biological "
    "replicates per species, with the high-confidence reference defined as the sites "
    "modified in both replicates at a modification ratio > 0.1 (Arabidopsis 80,624; "
    "mouse 41,961; human 112,451 high-confidence m6A sites). (B) Positive predictive "
    "value (PPV) against GLORI within the explicit candidate-site set at the "
    "primary 2-bp matching window, computed per independent sequencing unit (Arabidopsis "
    "and HeLa, three replicates; mouse, its two independent studies, never pooled); bars "
    "are group means and error bars the SD across units. (C) Agreement between the "
    "modification ratio reported by each tool and the GLORI ratio at shared sites, one "
    "panel per species; per-replicate Pearson r and Lin's concordance correlation "
    "coefficient (CCC) are printed for every tool."
)


CAP = Path(__file__).resolve().parent / "sup_figure_captions.md"


def caption(i: int) -> str:
    blocks = [b.strip() for b in re.split(r"\n\s*\n", CAP.read_text())]
    hit = [b for b in blocks if re.match(r"\*{0,2}Figure S%d\." % i, b.lstrip("> ").strip())]
    if hit:
        raw = hit[0]
    elif i == 3:
        raw = S3
    else:
        text = LEG[i].read_text()
        blocks = [b.strip() for b in re.split(r"\n\s*\n", text)]
        cand = [b for b in blocks
                if not b.lstrip().startswith("#") and re.match(r"\*{0,2}Figure S%d\b" % i, b.lstrip())]
        if not cand:
            cand = [b for b in blocks
                    if not b.lstrip().startswith("#") and "Figure S" in b]
        if not cand:
            cand = [b for b in blocks if re.search(r"\*{0,2}Figure S%d\b" % i, b)]
        raw = cand[0]
    raw = raw.replace("$", " ")
    for a, b in (("%", r"\%"), ("&", r"\&"), ("_", r"\_"), ("#", r"\#"),
                 ("{", "("), ("}", ")"), ("~", r"\textasciitilde{}"),
                 ("^", r"\textasciicircum{}")):
        raw = raw.replace(a, b)
    raw = re.sub(r"\*\*(.+?)\*\*", lambda m: "\\textbf{" + m.group(1) + "}", raw, flags=re.S)
    raw = raw.replace("\u2013", "--").replace("\u2014", "---").replace("\u2212", "-")
    for a, b in UNI.items():
        raw = raw.replace(a, b)
    raw = raw.replace("*", "")
    raw = re.sub(r"(^|\s)>\s*", " ", raw)
    raw = re.sub(r"\s+", " ", raw).strip()
    return raw


UNI = {
    "\u207b": r"\textsuperscript{$-$}", "\u2076": r"\textsuperscript{6}",
    "\u2075": r"\textsuperscript{5}", "\u2074": r"\textsuperscript{4}",
    "\u00b2": r"\textsuperscript{2}", "\u00b9": r"\textsuperscript{1}",
    "\u00b3": r"\textsuperscript{3}", "\u03a8": r"$\Psi$", "\u00b1": r"$\pm$",
    "\u00d7": r"$\times$", "\u2248": r"$\approx$", "\u2264": r"$\leq$",
    "\u2265": r"$\geq$", "\u0394": r"$\Delta$", "\u03b4": r"$\delta$",
    "\u03b5": r"$\epsilon$", "\u03c1": r"$\rho$", "\u03b1": r"$\alpha$",
    "\u03c4": r"$\tau$", "\u03bc": r"$\mu$", "\u201c": "``", "\u201d": "''",
    "\u2018": "`", "\u2019": "'",
    "\u2078": r"\textsuperscript{8}", "\u2080": r"\textsubscript{0}",
    "\u2081": r"\textsubscript{1}", "\u2082": r"\textsubscript{2}",
    "\u2086": r"\textsubscript{6}", "\u2229": r"$\cap$", "\u222a": r"$\cup$",
    "\u2208": r"$\in$", "\u2260": r"$\neq$", "\u00b7": r"$\cdot$",
}
parts = [r"\documentclass[a4paper]{article}", r"\usepackage[margin=1.6cm]{geometry}",
         r"\usepackage[T1]{fontenc}", r"\usepackage[utf8]{inputenc}",
         r"\usepackage{graphicx}", r"\pagestyle{plain}", r"\begin{document}"]
for i in range(1, 11):
    parts += [
        r"\begin{figure}[p]\centering",
        r"\includegraphics[width=\textwidth,height=0.80\textheight,keepaspectratio]"
        r"{sup/FigureS%d_rev.pdf}" % i,
        r"\par\vspace{6pt}\begin{minipage}{\textwidth}\footnotesize",
        caption(i),
        r"\end{minipage}\end{figure}\clearpage",
    ]
parts.append(r"\end{document}")
tex = (_XB / "submission/new_submission/merge_supplementary_figures.tex")
tex.write_text("\n".join(parts) + "\n")
r = subprocess.run([pdflatex,
                    "-interaction=nonstopmode", tex.name],
                   cwd=NEW, capture_output=True, text=True)
print("pdflatex exit", r.returncode)
for f in ("merge_supplementary_figures.pdf", "merge_supplementary_figures.aux",
          "merge_supplementary_figures.log"):
    if (NEW / f).exists() and f.endswith(".pdf"):
        ((_XB / "submission/new_submission/Supplementary_Figures.pdf")).write_bytes((NEW / f).read_bytes())
        print("wrote Supplementary_Figures.pdf")
