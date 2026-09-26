#!/usr/bin/env python3
"""47 -- Figure 7 (revised): non-m6A tools, paper narrative A-D."""

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
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve()
PROJECT = _RB
sys.path.insert(0, str((_RB / "src/harmonisation/common")))

import figstyle  # noqa: E402

EV = (_RB / "analysis/nonm6a_false_positives/evidence")
OUT = (_RB / "figures/figure7")
FIG = (_RB / "figures/figure7/figures")
LOG = (_RB / "figures/figure7/logs")
DENS = (_RB / "figures/figure7/tables/fig7_metagene_density.tsv")

#: restrained palette (paper-like): one blue WT, one orange unmodified control,
#: one neutral guide, one accent for the synthetic-control row.
C_WT = "#4B81B8"
C_IVT = "#E8A76B"
C_GREY = "0.45"
C_ACC = "#B3261E"

TOOL_LABEL = {
    "NanoNm": "NanoNm",
    "NanoMUD_psi": "NanoMUD-\u03a8",
    "NanoPsu": "NanoPsu",
    "NanoSPA_psU": "NanoSPA-\u03a8",
    "NanoMUD_m1psi": "NanoMUD-m1\u03a8",
    "CHEUI_m5C": "CHEUI-m5C",
}
#: three-class answer to R3-9, in the order used by the manuscript.
CLASSES = [("FP-dominated", ["CHEUI_m5C", "NanoMUD_psi", "NanoMUD_m1psi"]),
           ("Intermediate", ["NanoNm"]),
           ("Sparse", ["NanoPsu", "NanoSPA_psU"])]
A_TOOLS = [t for _, ts in CLASSES for t in ts]
A_POS = {t: i for i, t in enumerate(A_TOOLS)}
#: tools carrying an external reference (NanoMUD_m1Psi has none)
REF_TOOLS = ["CHEUI_m5C", "NanoMUD_psi", "NanoNm", "NanoPsu", "NanoSPA_psU"]
#: consecutive slots inside the enrichment panel (no empty slot for m1Psi)
REF_POS = {t: i for i, t in enumerate(REF_TOOLS)}

SEGS = ["five_prime_flank", "five_prime_UTR", "CDS", "three_prime_UTR",
        "three_prime_flank"]
GRID = 200
#: panel-D tool order -- same three-class order as panel A
F_ORDER = list(A_TOOLS)
COND_OF = {"WT": "HeLa_WT", "IVT": "HeLa_IVT"}
COND_COLOR = {"WT": C_WT, "IVT": C_IVT}

#: sizes are chosen for a 180 mm canvas that LaTeX scales to 0.95 \textwidth
#: (169 mm, x0.85), so the printed minimum stays above 7 pt.
##: The six y axis variables of rows A-C, in panel order.  They are named in
##: the caption because an in-figure y title is 90-degree type, which the
##: user forbade for these figures (2026-09-23).
Y_TITLES = False
Y_TITLES_TEXT = {"calls": "Calls per replicate",
                 "enrich": "Enrichment over chance",
                 "fp": "FP per 10$^{6}$\ncandidate sites",
                 "density": "Density",
                 "jaccard": "Global Jaccard",
                 "auc": "AUC, WT vs.\nunmodified IVT"}

STYLE = dict(font=9.6, label=10.6, tick=9.4, legend=9.0, letter=15.0)


def apply_style() -> None:
    figstyle.apply()
    mpl.rcParams.update({
        "font.size": STYLE["font"],
        "axes.labelsize": STYLE["label"],
        "axes.titlesize": STYLE["label"],
        "xtick.labelsize": STYLE["tick"],
        "ytick.labelsize": STYLE["tick"],
        "legend.fontsize": STYLE["legend"],
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.9,
        "xtick.major.width": 0.9,
        "ytick.major.width": 0.9,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "figure.dpi": 120,
    })


def letter(ax, s: str, dx: float = -0.005, dy: float = 1.055) -> None:
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=STYLE["letter"],
            fontweight="bold", va="bottom", ha="left")


def jitter(n: int, width: float = 0.055) -> np.ndarray:
    return np.linspace(-width, width, n) if n > 1 else np.zeros(1)


def pow_ticks(ax, lo: float, hi: float, fmt: str = "10$^{%d}$") -> None:
    """Log axis with decade ticks only (never '1e+06'-style labels)."""
    ax.set_yscale("log")
    ax.set_ylim(lo, hi)
    exps = [e for e in range(-12, 13) if lo <= 10.0 ** e <= hi]
    ax.set_yticks([10.0 ** e for e in exps])
    ax.set_yticklabels([fmt % e for e in exps])


def load() -> dict[str, pd.DataFrame]:
    """Five frozen evidence tables (read-only; the single source of numbers)."""
    return {
        "rep": pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/hela_wt_ivt_per_replicate.tsv"), sep="\t"),
        "sum": pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/hela_wt_ivt_summary.tsv"), sep="\t"),
        "cc": pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/curlcake_ivt_fp_per_construct.tsv"), sep="\t"),
        "truth": pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/truth_precision_per_replicate.tsv"), sep="\t"),
        "scor": pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/score_distributions.tsv"), sep="\t"),
        "eco": pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/ecoli_gse271571_prob.tsv"), sep="\t"),
    }


def load_density() -> pd.DataFrame:
    d = pd.read_csv(DENS, sep="\t")
    d = d[d["merge"].eq("majority") | d["merge"].str.startswith("unit:")].copy()
    d["curve"] = [np.fromiter(map(float, s.split("|")), float, GRID)
                  for s in d["density"]]
    return d


def spine(d: pd.DataFrame, condition: str, tool: str,
          merge: str) -> np.ndarray | None:
    """Concatenate the five segments of one curve (GRID points each)."""
    out = np.zeros(len(SEGS) * GRID)
    seen = False
    for i, name in enumerate(SEGS):
        r = d[(d.condition == condition) & (d.tool == tool)
              & (d["merge"] == merge) & (d.kind == name)]
        if len(r):
            out[i * GRID:(i + 1) * GRID] = r.curve.iloc[0]
            seen = True
    return out if seen else None


def cat_axis(ax, tools: list[str], classes: bool = False,
             local: bool = False) -> None:
    """Tool axis with rotated tick labels and an optional class-name row.

    ``local=True`` packs the given tools into consecutive slots (used by the
    enrichment panel, which has five tools with an external reference and must
    not inherit the empty sixth slot of the call-count panel).
    """
    if local:
        pos = {t: i for i, t in enumerate(tools)}
        ax.set_xlim(-0.62, len(tools) - 0.38)
    else:
        pos = A_POS
        ax.set_xlim(-0.62, len(A_TOOLS) - 0.38)
    if classes:
        for cname, group in CLASSES:
            idx = [A_POS[t] for t in group]
            ax.text(float(np.mean(idx)), -1.12, cname,
                    transform=ax.get_xaxis_transform(), ha="center", va="top",
                    fontsize=STYLE["tick"], color="0.15", fontweight="bold")
    ax.set_xticks([pos[t] for t in tools])
    ax.set_xticklabels([TOOL_LABEL[t] for t in tools], rotation=30, ha="right",
                       rotation_mode="anchor")


def panel_a1(ax, rep: pd.DataFrame) -> None:
    """Per-replicate calls: unmodified IVT libraries call as much as WT."""
    for tl in A_TOOLS:
        sub = rep[rep["tool"] == tl]
        x = A_POS[tl]
        for grp, off, filled, col in (("WT", -0.17, True, C_WT),
                                      ("IVT", 0.17, False, C_IVT)):
            rows = sub[sub["group"] == grp].sort_values("sample")
            vals = rows["n_calls_in_universe"].to_numpy(float)
            if not vals.size:
                continue
            ax.plot(x + off + jitter(len(vals), 0.05), vals, linestyle="none",
                    marker="o", markersize=3.9,
                    markerfacecolor=col if filled else "white",
                    markeredgecolor=col, markeredgewidth=0.9, zorder=3)
            ax.plot([x + off - 0.15, x + off + 0.15], [vals.mean()] * 2,
                    color=col, linewidth=1.4, zorder=4, solid_capstyle="butt")
    pow_ticks(ax, 1e2, 2e5)
    cat_axis(ax, A_TOOLS, classes=True)
    
    ## is a rotated label by definition, so every axis variable is named in
    ## the caption instead (panels A-F all carry their metric there, see
    ## manuscript.tex, Figure 7).  Set Y_TITLES=True to bring the six back.
    if Y_TITLES:
        ax.set_ylabel(Y_TITLES_TEXT['calls'])
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=4.6,
                      markerfacecolor=C_WT, markeredgecolor=C_WT, label="HeLa WT"),
               Line2D([], [], marker="o", linestyle="none", markersize=4.6,
                      markerfacecolor="white", markeredgecolor=C_IVT,
                      markeredgewidth=0.9, label="HeLa IVT")]
    place_legend(ax, handles, slot="upper right", handletextpad=0.35,
              borderpad=0.15, labelspacing=0.25)


def panel_a2(ax, truth: pd.DataFrame) -> None:
    """Truth-anchored enrichment over the chromosome-stratified permutation null."""
    t = truth[(truth["window_bp"] == 1)
              & (truth["reference"].isin(["RMBase+DirectRMDB", "NGS"]))]
    markers = {"RMBase+DirectRMDB": "o", "NGS": "s"}
    floor = 1e-2
    for tl in REF_TOOLS:
        x = REF_POS[tl]
        for ref, mk in markers.items():
            sub = t[(t["tool"] == tl) & (t["reference"] == ref)]
            if not len(sub):
                continue
            for grp, off, filled, col in (("WT", -0.17, True, C_WT),
                                          ("IVT", 0.17, False, C_IVT)):
                rows = sub[sub["group"] == grp]
                if not len(rows):
                    continue
                ys = rows["enrichment"].to_numpy(float)
                ys = np.where(np.isfinite(ys) & (ys > 0), ys, floor)
                ax.plot(x + off + jitter(len(rows), 0.05), ys, linestyle="none",
                        marker=mk, markersize=3.9,
                        markerfacecolor=col if filled else "white",
                        markeredgecolor=col, markeredgewidth=0.9, zorder=3)
    ax.axhline(1.0, color=C_GREY, linewidth=1.0, linestyle=(0, (2, 2)), zorder=1)
    pow_ticks(ax, 3e-2, 4e2)
    cat_axis(ax, REF_TOOLS, local=True)
    if Y_TITLES:          # axis variable named in the caption
        ax.set_ylabel(Y_TITLES_TEXT['enrich'])
    # No legend for this panel: every corner carries data points, so the corner
    # search kept escaping above the axes, and the two keys are one clause of
    # the caption ("circles, RMBase + DirectRMDB compilation; squares,
    # orthogonal NGS reference sets; dashed line, chance").  Filled WT / open
    # unmodified IVT is panel A's key and the caption's lead-in.


def panel_b1(ax, cc: pd.DataFrame) -> None:
    """Normalised false-positive density on the unmodified Curlcake constructs."""
    tools = [t for t in A_TOOLS if t in set(cc["tool"])]
    ind = cc.loc[cc["construct_role"] == "independent", "fp_per_1e6_candidates"]
    floor = float(ind[ind > 0].min()) * 0.02
    pos = {t: i for i, t in enumerate(tools)}   # tight slots, no empty column
    for tl in tools:
        sub = cc[cc["tool"] == tl]
        x = pos[tl]
        for role, filled in (("independent", True),
                             ("depth_matched_subset_of_rep3", False)):
            rows = sub[sub["construct_role"] == role]
            if not len(rows):
                continue
            vals = rows["fp_per_1e6_candidates"].to_numpy(float)
            ax.plot(x + jitter(len(vals), 0.11), np.maximum(vals, floor),
                    linestyle="none", marker="o", markersize=4.0,
                    markerfacecolor=C_ACC if filled else "white",
                    markeredgecolor=C_ACC, markeredgewidth=0.9, zorder=3)
        vals = sub.loc[sub["construct_role"] == "independent",
                       "fp_per_1e6_candidates"].to_numpy(float)
        if vals.size:
            ax.plot([x - 0.22, x + 0.22], [max(vals.mean(), floor)] * 2,
                    color="0.2", linewidth=1.5, zorder=4, solid_capstyle="butt")
    pow_ticks(ax, floor * 0.5, 2e5)
    ax.set_yticks([floor] + [10.0 ** e for e in (2, 3, 4)])
    ax.set_yticklabels(["0", "10$^{2}$", "10$^{3}$", "10$^{4}$"])
    cat_axis(ax, tools, local=True)
    if Y_TITLES:          # axis variable named in the caption
        ax.set_ylabel(Y_TITLES_TEXT['fp'])
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=4.6,
                      markerfacecolor=C_ACC, markeredgecolor=C_ACC,
                      label="construct"),
               Line2D([], [], marker="o", linestyle="none", markersize=4.6,
                      markerfacecolor="white", markeredgecolor=C_ACC,
                      markeredgewidth=0.9,
                      label="depth-matched subset"),
               Line2D([], [], color="0.2", linewidth=1.5,
                      label="mean of constructs")]
    
    #: panel carries a single legend line instead of a stacked three-row box;
    #: the corner search keeps it collision-free
    place_legend(ax, handles, slot="lower right", ncol=3, handletextpad=0.35,
              borderpad=0.15, labelspacing=0.25, columnspacing=0.9)


def panel_b2(ax, eco: pd.DataFrame) -> None:
    """Third-party GSE271571: CHEUI probabilities on unmodified E. coli RNA."""
    bins = eco[eco["row_type"] == "prob_bin"]
    styles = {"m5C": "-", "m6A": (0, (3, 2))}
    for mod, ls in styles.items():
        for cond, col in (("WT", C_WT), ("IVT", C_IVT)):
            d = bins[(bins["mod_type"] == mod) & (bins["condition"] == cond)]
            d = d.sort_values("bin_lo")
            centres = d["bin_lo"].to_numpy(float) + 0.01
            ax.plot(centres, d["density"].to_numpy(float), color=col,
                    linestyle=ls, linewidth=1.3, zorder=3)
    ax.set_xlim(0, 1)
    ax.set_ylim(bottom=0)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["0", "0.5", "1"])
    ax.set_xlabel("CHEUI probability")
    if Y_TITLES:          # axis variable named in the caption
        ax.set_ylabel(Y_TITLES_TEXT['density'])
    # two entries only: the modification is the line style, the condition is
    # the colour (blue WT, orange unmodified IVT -- stated in the caption)
    handles = [Line2D([], [], color="0.25", linewidth=1.4,
                      label="m$^{5}$C"),
               Line2D([], [], color="0.25", linewidth=1.4, linestyle=(0, (3, 2)),
                      label="m$^{6}$A")]
    place_legend(ax, handles, slot="upper right", handletextpad=0.45,
              borderpad=0.15, labelspacing=0.25)


def panel_c1(ax, summ: pd.DataFrame) -> None:
    """Replicate consistency: global vs mean pairwise Jaccard (log-log)."""
    tools = [t for t in A_TOOLS if t in set(summ["tool"])]
    ax.plot([3e-4, 0.6], [3e-4, 0.6], color="0.8", linewidth=0.8,
            linestyle=(0, (2, 2)), zorder=1)
    for tl in tools:
        row = summ[summ["tool"] == tl].iloc[0]
        for grp, filled, col in (("wt", True, C_WT), ("ivt", False, C_IVT)):
            gx = float(row[f"{grp}_mean_pairwise_jaccard_uni"])
            gy = float(row[f"{grp}_global_jaccard_uni"])
            ax.plot([gx], [gy], linestyle="none", marker="o", markersize=4.6,
                    markerfacecolor=col if filled else "white",
                    markeredgecolor=col, markeredgewidth=0.9, zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(3e-4, 1.2)
    ax.set_ylim(3e-4, 1.2)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_ticks([1e-3, 1e-2, 1e-1, 1.0])
        axis.set_ticklabels(["10$^{-3}$", "10$^{-2}$", "10$^{-1}$", "1"],
                            fontsize=STYLE["tick"])
    ax.set_xlabel("Mean pairwise Jaccard")
    if Y_TITLES:          # axis variable named in the caption
        ax.set_ylabel(Y_TITLES_TEXT['jaccard'])
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=4.6,
                      markerfacecolor=C_WT, markeredgecolor=C_WT, label="HeLa WT"),
               Line2D([], [], marker="o", linestyle="none", markersize=4.6,
                      markerfacecolor="white", markeredgecolor=C_IVT,
                      markeredgewidth=0.9, label="HeLa IVT")]
    place_legend(ax, handles, slot="upper right", handletextpad=0.35,
              borderpad=0.15, labelspacing=0.25)


def panel_c2(ax, scor: pd.DataFrame) -> None:
    """Discrimination of each tool's own score between WT and unmodified IVT."""
    d = scor[scor["sample"] == "__discrimination__"]
    tools = [t for t in A_TOOLS if t in set(d["tool"])]
    ax.axhline(0.5, color=C_GREY, linewidth=1.0, linestyle=(0, (2, 2)), zorder=1)
    for tl in tools:
        auc = float(d[d["tool"] == tl]["auc"].iloc[0])
        ax.plot([A_POS[tl]], [auc], linestyle="none", marker="o", markersize=4.6,
                markerfacecolor="0.25", markeredgecolor="0.25", zorder=3)
    ax.set_ylim(0.15, 0.72)
    ax.set_yticks([0.2, 0.4, 0.6])
    cat_axis(ax, tools)
    if Y_TITLES:          # axis variable named in the caption
        ax.set_ylabel(Y_TITLES_TEXT['auc'])
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=4.6,
                      markerfacecolor="0.25", markeredgecolor="0.25",
                      label="Score discrimination"),
               Line2D([], [], color=C_GREY, linewidth=1.0, linestyle=(0, (2, 2)),
                      label="Chance")]
    
    #: empty band directly above the panel
    place_legend(ax, handles, slot="lower right", ncol=2, handletextpad=0.35,
              borderpad=0.15, labelspacing=0.25,
              bbox_to_anchor=(1.0, 1.03))


def metagene_frame(ax, ylabel: bool) -> None:
    """Paper-style frame: segment labels, dotted separators, gene-model bar."""
    n = len(SEGS)
    ax.set_xticks(np.arange(n) + 0.5)
    ax.set_xticklabels([figstyle.SEGMENT_LABELS[s] for s in SEGS],
                       fontsize=STYLE["tick"])
    ax.set_xlim(0, n)
    ax.set_ylim(bottom=0)
    ax.minorticks_off()
    for i in range(1, n):
        ax.axvline(i, color="k", lw=1.0, ls=(0, (1, 1.6)), zorder=1)
    if ylabel:
        ax.set_ylabel("Density")
    tr = ax.get_xaxis_transform()
    ax.plot([0, n], [-0.030, -0.030], color="0.35", lw=1.0, transform=tr,
            clip_on=False, zorder=6)
    ax.plot([1, n - 1], [-0.008, -0.008], color="0.80", lw=6.0, transform=tr,
            clip_on=False, zorder=5, solid_capstyle="butt")
    for a, b in ((0, 1), (n - 1, n)):
        ax.plot([a, b], [-0.008, -0.008], color="k", lw=3.2, transform=tr,
                clip_on=False, zorder=6, solid_capstyle="butt")


def panel_d(f_axes, d: pd.DataFrame) -> None:
    """Replicate-aware metagene: majority consensus + per-replicate curves."""
    for idx, (ax, tool) in enumerate(zip(f_axes.ravel(), F_ORDER)):
        metagene_frame(ax, ylabel=idx % 3 == 0)
        tops = []
        for cond in ("WT", "IVT"):
            col = COND_COLOR[cond]
            main = spine(d, COND_OF[cond], tool, "majority")
            if main is None:
                continue
            tops.append(float(main.max()))
            x = (np.arange(len(main)) + 0.5) / GRID
            reps = sorted(m.split(":", 1)[1]
                          for m in d.loc[d.condition == COND_OF[cond], "merge"]
                          if m.startswith("unit:"))
            for u in reps:
                c = spine(d, COND_OF[cond], tool, f"unit:{u}")
                if c is not None:
                    ax.plot(x, c, color=col, lw=0.6, alpha=0.55,
                            ls=(0, (2.5, 2)), zorder=2)
            ax.fill_between(x, 0, main, color="0.86" if cond == "WT" else "#F6DCC3",
                            lw=0, zorder=3)
            ax.plot(x, main, color=col, lw=1.4, zorder=4)
        ax.set_ylim(0, (max(tops) if tops else 1.0) * 1.1)
        handles = [Line2D([], [], color=COND_COLOR[c], linewidth=4.0, alpha=0.9,
                          label=f"{TOOL_LABEL[tool]}-{c}") for c in ("WT", "IVT")]
        place_legend(ax, handles, slot="upper left", bbox_to_anchor=(0.0, -0.26),
                  ncol=1, handlelength=1.0, handleheight=0.9, handletextpad=0.3,
                  borderpad=0.05, columnspacing=0.55, labelspacing=0.15)


def _legend_hits(ax, leg, fig) -> list[str]:
    """Labels of the data/text artists overlapped by a legend box."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    lb = leg.get_window_extent(rend)
    hits: list[str] = []
    for coll in list(ax.collections) + list(ax.lines) + list(ax.patches):
        if not coll.get_visible():
            continue
        bb = coll.get_window_extent(rend)
        if bb.width * bb.height > 0 and lb.overlaps(bb):
            hits.append(type(coll).__name__)
    others = (list(ax.texts) + list(ax.get_xticklabels())
              + list(ax.get_yticklabels()))
    for name in ("xaxis", "yaxis"):
        lab = getattr(ax, name).get_label()
        if lab.get_visible() and str(lab.get_text()).strip():
            others.append(lab)
    for t in others:
        if t.get_visible() and lb.overlaps(t.get_window_extent(rend)):
            hits.append("text")
    for o in ax.figure.axes:               # other panels' keys must not collide
        if o is ax:
            continue
        ol = o.get_legend()
        if ol is not None and lb.overlaps(ol.get_window_extent(rend)):
            hits.append("companion-legend")
    return hits


def relayout_legends(fig) -> None:
    """Re-run the legend search once every panel letter / label is on the canvas."""
    for _ in range(2):                     # second pass sees all companion keys
        for ax in fig.axes:
            spec = getattr(ax, "_fig7_legend", None)
            if spec is None:
                continue
            if ax.get_legend() is not None:
                ax.get_legend().remove()
            place_legend(ax, spec[0], **spec[1])


def place_legend(ax, handles, **kw):
    """Legend in the first corner that overlaps no data artist or text.

    ``frameon=False`` is kept (house style), so a legend sitting on top of a
    marker is a visible collision -- hence the automated corner search.
    """
    kw.setdefault("frameon", False)
    ax._fig7_legend = (handles, dict(kw))      # kept for the final relayout
    kw.pop("loc", None)          # the slot is chosen by the search below
    ncols = int(kw.pop("ncol", 1))
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    fbox = fig.bbox
    sib = [o.get_window_extent(rend) for o in fig.axes
           if o is not ax and o.get_visible()]
    sib_text = []
    for o in fig.axes:
        if o is ax or not o.get_visible():
            continue
        for t in (list(o.texts) + list(o.get_xticklabels())
                  + list(o.get_yticklabels())):
            if t.get_visible() and str(t.get_text()).strip():
                sib_text.append(t)
    sib_boxes = [t.get_window_extent(rend) for t in sib_text]
    nc_variants = sorted({ncols, max(1, ncols // 2),
                          min(2 * ncols, len(handles))})
    cands: list[tuple[object, int, object]] = [
        ("corner", loc, nc)
        for loc in ("upper right", "upper left", "lower right", "lower left")
        for nc in nc_variants]
    for y, nc in ((1.005, ncols), (-0.14, ncols), (-0.30, ncols),
                  (1.005, min(2 * ncols, len(handles))),
                  (1.12, ncols), (1.12, min(2 * ncols, len(handles))),
                  (1.24, min(2 * ncols, len(handles)))):
        cands.append(("outside", (0.5, y), nc))
    slot = kw.pop("slot", "inside")
    if slot == "below":                 # under the rotated tick labels
        cands = [("below", (0.5, y), nc) for y in (-0.80, -0.95, -0.66)
                 for nc in nc_variants] + cands
    elif slot == "above":
        cands = [("outside", (0.5, y), nc) for y in (1.30, 1.42, 1.18)
                 for nc in nc_variants] + cands
    elif slot in ("upper right", "lower right"):
        cands = [("corner", slot, nc) for nc in nc_variants] + cands
    leg = None
    for kind, anchor, nc in cands:
        if leg is not None:
            leg.remove()
        if kind == "corner":
            leg = ax.legend(handles=handles, loc=anchor, ncol=nc, **kw)
        else:
            loc = "upper center" if kind == "below" else "lower center"
            leg = ax.legend(handles=handles, loc=loc, ncol=nc,
                            bbox_to_anchor=anchor, **kw)
        if _legend_hits(ax, leg, fig):
            continue
        bb = leg.get_window_extent(rend)
        if (bb.x0 < fbox.x0 or bb.x1 > fbox.x1
                or bb.y0 < fbox.y0 or bb.y1 > fbox.y1):
            continue
        if any(bb.overlaps(s) for s in sib):
            continue
        if any(bb.overlaps(s) for s in sib_boxes):
            continue
        return leg
    # nothing clean: keep the candidate with the smallest overlap area
    best = (float("inf"), leg)
    for kind, anchor, nc in cands:
        if leg is not None:
            leg.remove()
        loc = ("upper center" if kind == "below" else "lower center")
        leg = (ax.legend(handles=handles, loc=anchor, ncol=nc, **kw)
               if kind == "corner" else
               ax.legend(handles=handles, loc=loc, ncol=nc,
                         bbox_to_anchor=anchor, **kw))
        bb = leg.get_window_extent(rend)
        area = 0.0
        for s in sib + sib_boxes:
            dx = min(bb.x1, s.x1) - max(bb.x0, s.x0)
            dy = min(bb.y1, s.y1) - max(bb.y0, s.y0)
            if dx > 1.0 and dy > 1.0:
                area += dx * dy
        for hit_ax in (ax,):
            for coll in list(hit_ax.collections) + list(hit_ax.lines):
                cb = coll.get_window_extent(rend)
                dx = min(bb.x1, cb.x1) - max(bb.x0, cb.x0)
                dy = min(bb.y1, cb.y1) - max(bb.y0, cb.y0)
                if dx > 1.0 and dy > 1.0:
                    area += dx * dy
        if area < best[0]:
            best = (area, leg)
    return best[1]


def print_preview(src: Path, dst: Path, width_mm: float = 169.0) -> None:
    """Down-scale to the printed width (0.95 * textwidth) for a human check."""
    from PIL import Image

    im = Image.open(src)
    target = int(round(width_mm / 25.4 * 300))
    im.resize((target, int(round(im.height * target / im.width))),
              Image.LANCZOS).save(dst, dpi=(300, 300))
    print("print preview ->", dst, f"({width_mm:.0f} mm wide @ 300 dpi)")


def layout_report(fig) -> dict:
    """Measure the printed text objects: min font size and overlap pairs."""
    import json
    import matplotlib.text as mtext

    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    items = []
    for t in fig.findobj(mtext.Text):
        if not t.get_visible() or not str(t.get_text()).strip():
            continue
        try:
            bb = t.get_window_extent(renderer=rend)
        except Exception:                                    # pragma: no cover
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        items.append({"s": str(t.get_text())[:40], "fs": float(t.get_fontsize()),
                      "bb": [float(bb.x0), float(bb.y0), float(bb.x1), float(bb.y1)]})
    # axes-level check: the tight bbox of every panel includes its ticks, axis
    # labels, titles and legends, so disjoint tight bboxes rule out the
    # E-label / F-title style collisions between neighbouring panels.
    axes = [ax for ax in fig.axes if ax.get_visible()]
    # plotting rectangles only: the tight bbox would include an outside legend
    # and a legend that merely sits above its own panel is not a panel clash
    tbs = [ax.get_window_extent(rend) for ax in axes]
    overlaps = []
    for i in range(len(axes)):
        for j in range(i + 1, len(axes)):
            a, b = tbs[i], tbs[j]
            dx = min(a.x1, b.x1) - max(a.x0, b.x0)
            dy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if dx > 1.0 and dy > 1.0:
                overlaps.append({"kind": "axes", "a": f"axes{i}", "b": f"axes{j}",
                                 "dx": round(dx, 1), "dy": round(dy, 1)})
    # free-standing text (panel letters, class labels) vs the other panels
    for i, ax in enumerate(axes):
        for t in ax.texts:
            if not t.get_visible() or not str(t.get_text()).strip():
                continue
            tb = t.get_window_extent(renderer=rend)
            own = (list(ax.get_xticklabels()) + list(ax.get_yticklabels())
                   + [ax.xaxis.get_label(), ax.yaxis.get_label()])
            for o in own:                      # own tick labels (rotated rows)
                if not o.get_visible() or not str(o.get_text()).strip():
                    continue
                ob = o.get_window_extent(renderer=rend)
                if (min(tb.x1, ob.x1) - max(tb.x0, ob.x0) > 1.0
                        and min(tb.y1, ob.y1) - max(tb.y0, ob.y0) > 1.0):
                    overlaps.append({"kind": "text-tick",
                                     "a": str(t.get_text())[:24],
                                     "b": str(o.get_text())[:24], "dx": 0, "dy": 0})
            for j, other in enumerate(axes):
                if i == j:
                    continue
                o = tbs[j]
                dx = min(tb.x1, o.x1) - max(tb.x0, o.x0)
                dy = min(tb.y1, o.y1) - max(tb.y0, o.y0)
                if dx > 1.0 and dy > 1.0:
                    overlaps.append({"kind": "text", "a": str(t.get_text())[:24],
                                     "b": f"axes{j}", "dx": round(dx, 1),
                                     "dy": round(dy, 1)})
    # legends: box + every legend text against data artists, texts and neighbours
    legend_rows: list[dict] = []
    for i, ax in enumerate(axes):
        leg = ax.get_legend()
        if leg is None:
            continue
        lb = leg.get_window_extent(renderer=rend)
        for coll in list(ax.collections) + list(ax.lines) + list(ax.patches):
            if not coll.get_visible():
                continue
            bb = coll.get_window_extent(renderer=rend)
            if bb.width * bb.height <= 0:
                continue
            dx = min(lb.x1, bb.x1) - max(lb.x0, bb.x0)
            dy = min(lb.y1, bb.y1) - max(lb.y0, bb.y0)
            if dx > 1.0 and dy > 1.0:
                overlaps.append({"kind": "legend", "a": f"legend{i}",
                                 "b": type(coll).__name__, "dx": round(dx, 1),
                                 "dy": round(dy, 1)})
        texts = list(ax.texts) + list(ax.get_xticklabels()) + list(ax.get_yticklabels())
        for o in axes:                     # companion labels must stay clear too
            if o is ax:
                continue
            texts += (list(o.texts) + list(o.get_xticklabels())
                      + list(o.get_yticklabels()))
        for t in texts:
            if not t.get_visible():
                continue
            tb = t.get_window_extent(renderer=rend)
            dx = min(lb.x1, tb.x1) - max(lb.x0, tb.x0)
            dy = min(lb.y1, tb.y1) - max(lb.y0, tb.y0)
            if dx > 1.0 and dy > 1.0:
                overlaps.append({"kind": "legend", "a": f"legend{i}",
                                 "b": str(t.get_text())[:24], "dx": round(dx, 1),
                                 "dy": round(dy, 1)})
        legend_rows.append({"panel": ax.get_label(), "loc": leg._loc,
                            "x0": round(lb.x0, 1), "y0": round(lb.y0, 1),
                            "w": round(lb.width, 1), "h": round(lb.height, 1)})
    if legend_rows:
        LOG.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(legend_rows).to_csv((_RB / "figures/figure7/logs/47_layout_report.tsv"), sep="\t",
                                         index=False)
    ticks = []
    for i, ax in enumerate(axes):
        labs = [t for t in ax.get_xticklabels() + ax.get_yticklabels()
                if t.get_visible() and str(t.get_text()).strip()]
        ticks += [t.get_text() for t in labs]
        for m in range(len(labs)):                  # tick labels must not touch
            for n in range(m + 1, len(labs)):
                a, b = labs[m], labs[n]
                if abs(a.get_rotation()) < 1e-6:     # horizontal: exact bbox test
                    ba = a.get_window_extent(renderer=rend)
                    bb_ = b.get_window_extent(renderer=rend)
                    dx = min(ba.x1, bb_.x1) - max(ba.x0, bb_.x0)
                    dy = min(ba.y1, bb_.y1) - max(ba.y0, bb_.y0)
                    hit = dx > 1.0 and dy > 1.0
                else:                                # rotated: parallel-baseline gap
                    pa = a.get_transform().transform(a.get_position())
                    pb = b.get_transform().transform(b.get_position())
                    th = np.deg2rad(a.get_rotation())
                    perp = abs((pa[0] - pb[0]) * np.sin(th)) + \
                        abs((pa[1] - pb[1]) * np.cos(th))
                    hit = perp < 0.8 * float(a.get_fontsize())
                if hit:
                    overlaps.append({"kind": "tick", "a": a.get_text()[:16],
                                     "b": b.get_text()[:16], "ax": f"axes{i}"})
    
    # instead of shipping one, whatever produced it (ylabel, tick, annotation)
    rotated = sorted({str(x.get_text())[:24] for x in fig.findobj(mtext.Text)
                      if x.get_visible() and str(x.get_text()).strip()
                      and abs(float(x.get_rotation()) % 180.0 - 90.0) < 1e-6})
    if rotated:
        raise SystemExit("90-degree type inside the figure: "
                         + ", ".join(rotated[:8]))
    report = {"min_fontsize": min((it["fs"] for it in items), default=0.0),
              "n_text": len(items), "n_axes": len(axes),
              "n_overlap": len(overlaps), "overlaps": overlaps[:40],
              "tick_labels": ticks}
    LOG.mkdir(parents=True, exist_ok=True)
    ((_RB / "figures/figure7/logs/47_layout.json")).write_text(json.dumps(report, indent=1))
    print(f"layout: {report['n_text']} text objects, min font "
          f"{report['min_fontsize']:.1f} pt, {report['n_overlap']} overlaps")
    return report


def build_figure() -> tuple[plt.Figure, dict[str, pd.DataFrame]]:
    """Assemble the four-row figure (no file output; reused by the verifier)."""
    apply_style()
    data = load()
    dens = load_density()

    # only rows A-C live here since v3: row D (metagene) is drawn by the Guitar
    # R pipeline (48_fig7d_guitar.R) and the two rows are stitched afterwards.
    # drawn 193 mm wide and scaled to 180 mm by the assembly (scale 0.933), so
    # the rotated tick labels and the wide reference legend stay inside the page
    fig = plt.figure(figsize=(7.60, 6.10))           # rows A-C
    outer = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.06, 1.06],
                             hspace=2.05, left=0.115, right=0.975,
                             top=0.900, bottom=0.115)
    gs_a = outer[0].subgridspec(1, 2, wspace=0.38)
    gs_b = outer[1].subgridspec(1, 2, wspace=0.38)
    gs_c = outer[2].subgridspec(1, 2, wspace=0.38)

    ax_a1 = fig.add_subplot(gs_a[0, 0])
    ax_a2 = fig.add_subplot(gs_a[0, 1])
    ax_b1 = fig.add_subplot(gs_b[0, 0])
    ax_b2 = fig.add_subplot(gs_b[0, 1])
    ax_c1 = fig.add_subplot(gs_c[0, 0])
    ax_c2 = fig.add_subplot(gs_c[0, 1])
    panel_a1(ax_a1, data["rep"])
    panel_a2(ax_a2, data["truth"])
    panel_b1(ax_b1, data["cc"])
    panel_b2(ax_b2, data["eco"])
    panel_c1(ax_c1, data["sum"])
    panel_c2(ax_c2, data["scor"])

    letter(ax_a1, "A")             # per sub-panel letters A-F (G is on the R side)
    letter(ax_a2, "B")
    letter(ax_b1, "C")
    letter(ax_b2, "D")
    letter(ax_c1, "E")
    letter(ax_c2, "F")
    relayout_legends(fig)          # panel letters now exist: re-check legends
    rows = []
    for name, ax in (("A-left", ax_a1), ("A-right", ax_a2),
                     ("B-left", ax_b1), ("B-right", ax_b2)):
        for x, lab in zip(ax.get_xticks(), ax.get_xticklabels()):
            rows.append({"panel": name, "x": float(x),
                         "label": lab.get_text()})
    LOG.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv((_RB / "figures/figure7/logs/47_axis_report.tsv"), sep="\t", index=False)

    return fig, data


def main() -> None:
    fig, data = build_figure()
    FIG.mkdir(parents=True, exist_ok=True)
    stem = (_RB / "figures/figure7/figures/Figure7_rev_AC")
    # the row is stitched with the Guitar row later, so keep the nominal canvas:
    # bbox_inches="tight" would change the page width and the two rows would be
    # scaled by different factors when the page is assembled.
    for suffix, dpi in (("", 300), ("_600dpi", 600)):
        fig.savefig(f"{stem}{suffix}.pdf")
        fig.savefig(f"{stem}{suffix}.png", dpi=dpi)
    print("wrote", stem, "(pdf + png 300/600 dpi)")
    print_preview(f"{stem}.png", (_RB / "figures/figure7/figures/Figure7_print_preview.png"))
    layout_report(fig)
    plt.close(fig)


if __name__ == "__main__":
    main()
