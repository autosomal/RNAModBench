#!/usr/bin/env python
"""R3-6 six-panel figure, house-style redo (2026-09-19).

Visual evidence that motif profiles are dominated by the tool rather than by the
species.  The plotted content is unchanged w.r.t. the first version; only the
presentation follows the house rules confirmed by the author on 2026-09-19:

* no figure title and no sentence panel titles -- the only headings inside the
  figure are the bold panel letters A-F; every statistic that used to sit in a
  panel title is written to ``../figures/FigR3_6_legends.md``, which this script
  regenerates from the very tables it plots (legend cannot drift from the figure);
* every text element is >= 12 pt (in-cell KL values, eta^2 bar labels, 5-mer
  labels, tick labels, legends);
* no gridlines anywhere (heatmap cell borders removed, no ax.grid);
* the 13 tool + 3 species swatches live in a dedicated column on the right of
  the canvas instead of being crammed into the clustering panel;
* the m6Anet panel is a species x rank grid (colour = relative frequency) so the
  5-mer labels can no longer overlap each other.

The same code path renders both the legacy Top-5 inputs and the full-count
(callsets) inputs, selected purely by the CLI arguments.

Outputs: ../figures/<out>.pdf, ../figures/<out>.png
         ../figures/FigR3_6_legends.md  (section <out> replaced in place)

Usage (conda env benchmark-revision, from this directory):
  # legacy Top-5 inputs
  python r3_6_figures.py \
      --kl $RNAMODBENCH_LOCAL/legacy/next_postprocessing/Total_new/combined_kl_data.csv \
      --top5 $RNAMODBENCH_LOCAL/legacy/next_postprocessing/Total_new/All_Species_Top5_5mer_Data.csv \
      --eta kl_variance_decomposition.tsv --ggac ggac_agac_within_tool.tsv \
      --out FigR3_6_motif_tool_bias --dlab "top-5"
  # full-count inputs
  python r3_6_figures.py --kl fullcount_kl_pooled.tsv --top5 fullcount_top5.tsv \
      --eta fullcount_variance.tsv --ggac fullcount_ggac.tsv \
      --out FigR3_6_motif_tool_bias_fullcount --dlab "all called sites"
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
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.text as mtext
from matplotlib.patches import Patch
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import cophenet, dendrogram, linkage
from scipy.spatial.distance import squareform

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.normpath(os.path.join(HERE, "..", "figures"))
os.makedirs(FIGDIR, exist_ok=True)

# shared house style (Arial, no gridlines, pdf.fonttype 42, save() -> pdf+png)
sys.path.insert(0, str(_RB / "src/harmonisation"))
from common.figstyle import apply as apply_style, save  # noqa: E402

# style floor of this figure: nothing below 12 pt, panel letters much larger
TEXT_FLOOR = 12.0
RC = {
    "font.size": 13.5,
    "axes.labelsize": 15,
    "axes.titlesize": 20,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 12.5,
}

SP3 = ["Arabidopsis", "Mouse", "Human (HeLa)"]
pal = {"Arabidopsis": "#2ca02c", "Mouse": "#4f81bd", "Human (HeLa)": "#c1272d"}
LEGEND_MD = os.path.join(FIGDIR, "FigR3_6_legends.md")


ap = argparse.ArgumentParser()
ap.add_argument("--kl", required=True, help="KL table (legacy csv or fullcount tsv)")
ap.add_argument("--top5", required=True, help="Top-5 5-mer table")
ap.add_argument("--eta", required=True, help="variance-decomposition table")
ap.add_argument("--ggac", required=True, help="AGAC/GGAC fractions per tool x species")
ap.add_argument("--out", default="FigR3_6_motif_tool_bias")
ap.add_argument("--dlab", default="top-5", help="scope of the AGAC-GGAC counts (x label)")
A = ap.parse_args()


def read_auto(path, index_col=False):
    """Read a csv/tsv without being told which separator it uses."""
    with open(path) as fh:
        head = fh.readline()
    sep = "\t" if head.count("\t") > head.count(",") else ","
    return pd.read_csv(path, sep=sep, index_col=index_col)


def panel_letter(ax, ch):
    """The only heading a panel is allowed to carry."""
    ax.set_title(ch, loc="left", fontsize=20, fontweight="bold", pad=6)


def audit_figure(fig, floor=TEXT_FLOOR):
    """Report every text element below the font floor and any visible gridline."""
    small = []
    for t in fig.findobj(mtext.Text):
        s = (t.get_text() or "").strip()
        if not s:
            continue
        fs = t.get_fontsize()
        if fs is not None and float(fs) < floor:
            small.append((round(float(fs), 1), s[:46]))
    grids = []
    for ax in fig.axes:
        vis = [ln for ln in ax.get_xgridlines() + ax.get_ygridlines()
               if ln.get_visible() and ln.get_linestyle() not in ("None", " ")]
        if vis:
            grids.append(len(vis))
    print(f"[style-check] text elements below {floor:g} pt: {len(small)}")
    for fs, s in small:
        print(f"[style-check]     {fs:>5} pt  {s!r}")
    print(f"[style-check] axes with visible gridlines: {len(grids)} {grids}")
    return small, grids


def update_legend_md(key, body):
    """Replace (or append) the ``## <key>`` section of the legend file."""
    header = (
        "# R3-6 figure legends\n\n"
        "Auto-generated by `analysis/r3_6_figures.py`; every number below is computed from the\n"
        "same tables that are plotted, so the legend cannot drift from the figure.  The panels\n"
        "themselves carry only the letters A-F (house rule: no titles inside a figure).\n"
    )
    order, blocks, cur = [], {}, None
    if os.path.exists(LEGEND_MD):
        with open(LEGEND_MD) as fh:
            for line in fh:
                if line.startswith("## "):
                    cur = line[3:].strip()
                    order.append(cur)
                    blocks[cur] = []
                elif cur is not None:
                    blocks[cur].append(line)
    blocks[key] = [ln + "\n" for ln in body.rstrip().split("\n")]
    if key not in order:
        order.append(key)
    with open(LEGEND_MD, "w") as fh:
        fh.write(header)
        for k in order:
            fh.write(f"\n## {k}\n\n")
            fh.writelines(blocks[k])


# --------------------------------------------------------------------------- #
# data (read-only)
# --------------------------------------------------------------------------- #
kl = read_auto(A.kl)
top = read_auto(A.top5)
top = top[top.Species.isin(SP3)]
eta = read_auto(A.eta, index_col=0)
gg = read_auto(A.ggac, index_col=0)

apply_style()
plt.rcParams.update(RC)

fig = plt.figure(figsize=(12.6, 7.8), constrained_layout=True)
gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 0.55])
axes = [[fig.add_subplot(gs[r, c]) for c in range(3)] for r in range(2)]
lax = fig.add_subplot(gs[:, 3])
lax.axis("off")

# ---------------- A: KL heatmap ----------------
ax = axes[0][0]
km = kl.pivot(index="Tool", columns="Species_Name", values="KL_Divergence")[SP3]
sns.heatmap(km, annot=True, fmt=".1f", annot_kws={"fontsize": 12.5},
            cmap="rocket_r", linewidths=0.0, ax=ax,
            cbar_kws={"label": "KL vs RRACH", "shrink": 0.85, "pad": 0.02})
cb = ax.collections[0].colorbar
cb.ax.tick_params(labelsize=12)
cb.ax.yaxis.label.set_size(13)
ax.set_xlabel("")
ax.set_ylabel("")
plt.setp(ax.get_xticklabels(), rotation=25, ha="right", rotation_mode="anchor")
panel_letter(ax, "A")

# ---------------- B: variance decomposition ----------------
ax = axes[0][1]
vals = [float(eta.loc["Tool", "eta_sq"]), float(eta.loc["Species", "eta_sq"]),
        float(eta.loc["Interaction+Residual", "eta_sq"])]
terms = ["Tool", "Species", "Tool\u00d7Species\n+ residual"]
bars = ax.barh(terms[::-1], vals[::-1],
               color=["#c1272d", "#4f81bd", "#bbbbbb"][::-1])
for b, v in zip(bars, vals[::-1]):
    ax.text(v + 0.02, b.get_y() + b.get_height() / 2, f"{v:.2f}",
            va="center", ha="left", fontsize=13)
ax.set_xlim(0, 1.12)
ax.set_xlabel("variance explained ($\\eta^2$)")
panel_letter(ax, "B")

# ---------------- C: top-5 overlap by pair type ----------------
ax = axes[0][2]
sets = {k: set(g["5mer"]) for k, g in top.groupby(["Species", "Tool"], sort=False)}
recs = []
for (s1, t1), a in sets.items():
    for (s2, t2), b in sets.items():
        if (s1, t1) >= (s2, t2):
            continue
        if s1 == s2 and t1 != t2:
            recs.append(("cross tool\n(same species)", len(a & b)))
        elif t1 == t2 and s1 != s2:
            recs.append(("same tool\n(cross species)", len(a & b)))
ov = pd.DataFrame(recs, columns=["kind", "n"])
g1 = ov[ov.kind == "same tool\n(cross species)"].n
g2 = ov[ov.kind == "cross tool\n(same species)"].n
u = stats.mannwhitneyu(g1, g2, alternative="greater")
sns.boxplot(data=ov, x="kind", y="n", hue="kind",
            palette={"same tool\n(cross species)": "#c1272d",
                     "cross tool\n(same species)": "#4f81bd"},
            legend=False, ax=ax, width=0.55)
sns.stripplot(data=ov, x="kind", y="n", color="k", size=3.0, alpha=0.35, ax=ax)
plt.setp(ax.get_xticklabels(), rotation=20, ha="right", rotation_mode="anchor")
ax.set_ylabel("Top-5 5-mer overlap")
ax.set_xlabel("")
panel_letter(ax, "C")

# ---------------- D: within-tool AGAC-GGAC flip ----------------
ax = axes[1][0]
cols = {s: (f"{s}_AGAC", f"{s}_GGAC") for s in SP3}
rows = []
for tool, r in gg.iterrows():
    for sp in SP3:
        a, b = r[cols[sp][0]], r[cols[sp][1]]
        if a + b > 0:
            rows.append({"Tool": tool, "Species": sp, "diff": a - b})
dd = pd.DataFrame(rows)
sns.barplot(data=dd, x="diff", y="Tool", hue="Species", palette=pal, ax=ax,
            dodge=True, errorbar=None)
if ax.get_legend() is not None:          # species colours are in the right column
    ax.get_legend().remove()
ax.axvline(0, color="k", lw=1.0)
ax.set_xlabel(f"freq(AGAC{{A,U}}) \u2212 freq(GGAC{{A,U}})\n({A.dlab})")
ax.set_ylabel("")
ax.set_xlim(-0.9, 0.9)
panel_letter(ax, "D")

# ---------------- E: clustering of tool x species motif profiles ----------------
ax = axes[1][1]
keys = list(sets.keys())
tool_list = sorted({t for _, t in keys})
tool_pal = {t: matplotlib.colors.to_hex(c)
            for t, c in zip(tool_list, sns.color_palette("tab20", len(tool_list)))}
mer = sorted(top["5mer"].unique())
M = pd.DataFrame(0.0, index=range(len(keys)), columns=mer)
for i, k in enumerate(keys):
    for m in sets[k]:
        M.loc[i, m] = 1.0
I = M.dot(M.T).values
r = M.sum(axis=1).values
J = 1 - I / (r[:, None] + r[None, :] - I)
np.fill_diagonal(J, 0)
Z = linkage(squareform(J, checks=False), method="average")
dn = dendrogram(Z, ax=ax, no_labels=True, color_threshold=0)
leaforder = [keys[i] for i in dn["leaves"]]
trans = ax.get_xaxis_transform()
for i, (sp, t) in enumerate(leaforder):
    ax.add_patch(plt.Rectangle((i - 0.45, -0.045), 0.9, 0.035,
                               transform=trans, color=tool_pal[t], clip_on=False))
    ax.add_patch(plt.Rectangle((i - 0.45, -0.095), 0.9, 0.035,
                               transform=trans, color=pal[sp], clip_on=False))
cond = squareform(J)
cc = stats.pearsonr(cophenet(Z, cond)[1], cond)[0]
pairs_wt = [J[i, j] for i in range(len(keys)) for j in range(i + 1, len(keys))
            if keys[i][1] == keys[j][1]]
pairs_bt = [J[i, j] for i in range(len(keys)) for j in range(i + 1, len(keys))
            if keys[i][1] != keys[j][1] and keys[i][0] == keys[j][0]]
uE = stats.mannwhitneyu(pairs_wt, pairs_bt, alternative="less")
ax.set_ylabel("Jaccard distance (Top-5 sets)")
ax.set_xticks([])
ax.tick_params(axis="x", length=0)
ax.margins(x=0.02)
panel_letter(ax, "E")

# ---------------- F: m6Anet exemplar, species x rank grid ----------------
ax = axes[1][2]
m6 = top[top.Tool == "m6Anet"]
val = np.full((5, len(SP3)), np.nan)
txt = np.full((5, len(SP3)), "", dtype=object)
for j, sp in enumerate(SP3):
    s = m6[m6.Species == sp].sort_values("Rank").head(5)
    for i, (_, row) in enumerate(s.iterrows()):
        val[i, j] = row.Relative_Frequency
        txt[i, j] = row["5mer"]
cmap = matplotlib.colormaps["pink"]  # white (0) -> dark (1): darker cell = more frequent
norm = matplotlib.colors.Normalize(vmin=0.0, vmax=float(np.nanmax(val)))
ax.imshow(val, cmap=cmap, norm=norm, aspect="auto")
for i in range(5):
    for j in range(len(SP3)):
        if txt[i, j]:
            lum = float(np.dot(cmap(norm(val[i, j]))[:3], (0.299, 0.587, 0.114)))
            ax.text(j, i, txt[i, j], ha="center", va="center", fontsize=12.5,
                    color="white" if lum < 0.55 else "black")
ax.set_xticks(range(len(SP3)))
ax.set_xticklabels(SP3)
ax.set_yticks(range(5))
ax.set_yticklabels([f"#{k}" for k in range(1, 6)])
ax.set_xlabel("")
ax.set_ylabel("motif rank")
plt.setp(ax.get_xticklabels(), rotation=20, ha="right", rotation_mode="anchor")
cax = ax.inset_axes([0.30, -0.32, 0.40, 0.045])
smF = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
cbF = matplotlib.colorbar.Colorbar(cax, smF, orientation="horizontal")
cbF.set_label("relative frequency", fontsize=13, labelpad=2)
cbF.ax.tick_params(labelsize=12)
panel_letter(ax, "F")

# ---------------- right-hand legend column ----------------
leg_tools = lax.legend(
    handles=[Patch(color=tool_pal[t], label=t) for t in tool_list],
    loc="upper left", frameon=False, fontsize=12.5,
    title="tool  (upper bar of E)", title_fontsize=12.5,
    handlelength=1.1, handleheight=1.1, handletextpad=0.6, labelspacing=0.62)
lax.add_artist(leg_tools)
lax.legend(handles=[Patch(color=pal[s], label=s) for s in SP3],
           loc="lower left", frameon=False, fontsize=12.5,
           title="species  (lower bar of E;\ncolours of D and F)", title_fontsize=12.5,
           handlelength=1.1, handleheight=1.1, handletextpad=0.6, labelspacing=0.62)

# ---------------- legend markdown (replaces the deleted panel titles) ---------
n_tool_kl = km.shape[0]
kl_range = km.max(axis=1) - km.min(axis=1)
gg_pos = {sp: int((gg[f"{sp}_AGAC_minus_GGAC"] > 0).sum()) for sp in SP3}
gg_neg = {sp: int((gg[f"{sp}_AGAC_minus_GGAC"] < 0).sum()) for sp in SP3}
m6_lists = {sp: " / ".join(txt[i, j] for i in range(5) if txt[i, j])
            for j, sp in enumerate(SP3)}
body = f"""Inputs: `{os.path.basename(A.kl)}`, `{os.path.basename(A.top5)}`, \
`{os.path.basename(A.eta)}`, `{os.path.basename(A.ggac)}` ({A.dlab} counts; \
{n_tool_kl} tools with KL, {len(keys)} tool x species profiles, {len(SP3)} species).

- **A** KL divergence of each tool's site 5-mer profile from RRACH, per tool (rows)
  x species (columns); values are printed in the cells.  Within a row the three
  species differ by at most {kl_range.max():.1f} KL units (median within-tool range
  {kl_range.median():.1f}), whereas between tools the range spans
  {km.values.min():.0f}-{km.values.max():.0f}; i.e. the distance from RRACH is a
  property of the caller, not of the species.
- **B** Two-way variance decomposition of the same KL values: tool
  $\\eta^2$ = {vals[0]:.3f} (F = {float(eta.loc['Tool', 'F']):.1f},
  p = {float(eta.loc['Tool', 'p']):.1e}), species $\\eta^2$ = {vals[1]:.3f}
  (F = {float(eta.loc['Species', 'F']):.1f}, p = {float(eta.loc['Species', 'p']):.2f},
  not significant), tool x species + residual $\\eta^2$ = {vals[2]:.3f}.
- **C** Top-5 5-mer overlap per tool pair: same tool across two species (red,
  n = {len(g1)} pairs) vs two tools within one species (blue, n = {len(g2)} pairs);
  one-sided Mann-Whitney p = {u.pvalue:.1e}, median
  {g1.median():.0f} vs {g2.median():.0f} of 5.
- **D** Within-tool frequency difference freq(AGAC{{A,U}}) - freq(GGAC{{A,U}})
  ({A.dlab}): AGAC exceeds GGAC in {gg_pos['Arabidopsis']}/{gg.shape[0]} tools in
  Arabidopsis, while the reverse holds
  for human ({(gg[f'Human (HeLa)_AGAC_minus_GGAC'] < 0).sum()}/{gg.shape[0]} negative)
  and mouse ({(gg['Mouse_AGAC_minus_GGAC'] < 0).sum()}/{gg.shape[0]} negative) --
  the plant/animal direction survives tool control.
- **E** Hierarchical clustering (average linkage) of the {len(keys)} tool x species
  profiles by Jaccard distance of their Top-5 5-mer sets; the two bars under the
  leaves code the tool (upper, colours per the legend) and the species (lower).
  Mean distance {np.mean(pairs_wt):.2f} within the same tool vs {np.mean(pairs_bt):.2f}
  between tools of the same species (one-sided Mann-Whitney p = {uE.pvalue:.1e};
  cophenetic correlation r = {cc:.2f}).
- **F** m6Anet only: the Top-5 5-mers per species (cell label) with their relative
  frequency (cell shading, colour bar below).  Arabidopsis: {m6_lists['Arabidopsis']};
  mouse: {m6_lists['Mouse']}; human HeLa: {m6_lists['Human (HeLa)']} -- the same
  GAC-type motifs top the ranking in all three species, i.e. the human training
  prior travels with the tool.

(No p-value in this figure is corrected for multiple testing; the variance
decomposition treats the 13 tools as a methodological family, not as independent
replicates.)"""

update_legend_md(A.out, body)

audit_figure(fig)
save(fig, os.path.join(FIGDIR, A.out))
print(f"[stats] eta2 tool={vals[0]:.4f} F={float(eta.loc['Tool', 'F']):.3f} "
      f"p={float(eta.loc['Tool', 'p']):.3e} | species={vals[1]:.4f} "
      f"F={float(eta.loc['Species', 'F']):.3f} p={float(eta.loc['Species', 'p']):.4f}")
print(f"[stats] MWU overlap (same tool > cross tool) p={u.pvalue:.3e}; "
      f"MWU Jaccard (within tool < cross tool) p={uE.pvalue:.3e}; cophenetic r={cc:.3f}")
print(f"[stats] mean Jaccard within tool={np.mean(pairs_wt):.3f} "
      f"cross tool={np.mean(pairs_bt):.3f}; KL overall min={km.values.min():.1f} "
      f"max={km.values.max():.1f}")
print(f"[out] {os.path.join(FIGDIR, A.out)}.pdf/.png; legend -> {LEGEND_MD}")
