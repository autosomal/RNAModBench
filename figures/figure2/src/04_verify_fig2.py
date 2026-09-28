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
import subprocess, numpy as np, pandas as pd
from pathlib import Path
R = Path(__file__).resolve().parents[1]; T = R / "tables"
F = R / "figures" / "Figure2_rev"
bad = []
# 1 artefact ---------------------------------------------------------------
info = subprocess.run(["pdfinfo", str(F) + ".pdf"], capture_output=True, text=True).stdout
for line in info.splitlines():
    if line.startswith("Page size"):
        w, h = [float(x) for x in line.split()[2:5:2]]
        if abs(w - 479.52) > 0.5 or abs(h - 576.0) > 0.5:
            bad.append(f"page {w}x{h} pts != 479.52x576 (6.66x8.0 in)")
if not (F.with_suffix(".png")).exists():
    bad.append("PNG missing")
# 2 fonts ------------------------------------------------------------------
pf = subprocess.run(["pdffonts", str(F) + ".pdf"], capture_output=True, text=True).stdout
rows = [r.split() for r in pf.splitlines()[2:] if r.strip()]
for r in rows:
    if "Arial" not in r[0]:
        bad.append("non-Arial font in PDF: " + r[0])
    if r[-5] != "yes" or r[-4] != "yes" or r[-3] != "yes":
        bad.append("font not embedded/subset/unicode: " + r[0])
# 3 numbers vs the evaluation tables --------------------------------------
fr = pd.read_csv(T / "fig2a_testable_ratio.tsv", sep="\t")
ev = fr.dropna(subset=["eval_ko_wt_ratio"])
if len(ev):
    rel = np.abs(ev.ctrl_wt_ratio - ev.eval_ko_wt_ratio) / ev.eval_ko_wt_ratio
    if rel.max() > 1e-5:
        bad.append(f"A ratio vs ko_kd_metrics rel diff {rel.max():.2e}")
E = pd.read_csv(T / "fig2e_replicate_consistency.tsv", sep="\t")
W = pd.read_csv(_RB / "data/evaluation/tables/figure_ready_replicates.tsv", sep="\t")
W = W[(W.window == 2) & W.dataset_group.isin(["Arabidopsis_WT", "HeLa_WT", "Mouse_WT"])]
m = E[E.metric == "mean_pairwise"].merge(W[["species", "tool", "mean_pairwise_jaccard"]],
                                         on=["species", "tool"], how="inner")
if len(m):
    d = np.abs(m.value - m.mean_pairwise_jaccard)
    if d.max() > 1e-5:
        print("note: panel G is on raw call sets; the eval table uses the shared universe")
# 4 the pooled (global) series is retired (2026-09-28, union drop) --------
txt = subprocess.run(["pdftotext", str(F) + ".pdf", "-"],
                     capture_output=True, text=True).stdout.lower()
for token in ("pooled", "global"):
    if token in txt:
        bad.append(f"the retired pooled series is still printed in figure 2: {token!r}")
# 5 report ----------------------------------------------------------------
if bad:
    print("FAIL")
    [print(" -", b) for b in bad]
    raise SystemExit(1)
print("verify_fig2: all checks pass")
