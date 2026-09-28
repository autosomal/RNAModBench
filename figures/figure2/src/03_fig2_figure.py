#!/usr/bin/env python3
"""03 -- Figure 2 (revision): the assembled page.

Panels (all replicate-aware; mouse = two cross-study datasets, never pooled):

A  control x wild-type detection on the SHARED TESTABLE UNIVERSE
   (A1: log-log scatter of the two counts per unit pair, parity line;
    A2: control/WT ratio per unit pair + mean, parity line at 1)
B  inter-tool support: fraction of a unit's sites supported by >= k of the 13
   m6A tools (thin = one independent unit, thick = species mean)
C  species-native GO:BP enrichment of the majority high-confidence sites
   (top terms per species, x = fold enrichment)
D  wild-type x control agreement and the explicit control-side false-positive
   burden (D1: Jaccard per unit pair + mean; D2: control calls not within 2 bp
   of GLORI, per 10^6 candidate positions, log scale)
E  replicate consistency per tool: mean pairwise Jaccard vs the pooled
   ("global") Jaccard |intersection| / |union| over the same independent units

House rules enforced here (see src/harmonisation/common/figstyle.py):
real embedded Arial, no gridlines, NO annotation text or numbers inside a
panel, English only, vector PDF + 300 dpi PNG drawn at the printed size
(6.66 x 7.4 in = 0.95 x \textwidth of the Wiley USG layout).

Usage
-----
conda run -n benchmark-revision --no-capture-output python 03_fig2_figure.py
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
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import pandas as pd
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

HERE = Path(__file__).resolve()
REV = HERE.parents[1]
sys.path.insert(0, str(_RB / "src/harmonisation"))
from common import figstyle  # noqa: E402

#: published species palette (kept so the revision reads as the same figure)
COLOR = {"Arabidopsis": "#1E888B", "Mouse": "#F5B264", "Human": "#3778A0"}
MARKER = {"Arabidopsis": "o", "Mouse": "s", "Human": "^"}
BLOCKS = ["Arabidopsis", "Mouse", "Human"]
TOOLS = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
         "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
         "NanoSPA_m6A", "xPore", "yanocomp"]
#: R1-7 standardised labels
LABEL = {"yanocomp": "Yanocomp", "xPore": "xPore", "Nanom6A": "Nanom6A"}

FS_TICK, FS_LABEL, FS_PANEL, FS_LEGEND = 7.5, 9.0, 12.0, 8.5
PAGE_W, PAGE_H = 6.66, 8.00
#: tool-name tick rotation (2026-09-24: 88 deg, the largest tilt with zero
#: measured collisions -- qa_tick_overlaps: 88 deg "none", 85 deg 3+3 pairs,
#: 80 deg 7-12 pairs; 45 deg needs 2.3 pt type, far under the 7 pt floor)
TICK_ROT = 88


def lab(tool: str) -> str:
    return LABEL.get(tool, tool)


def tool_ticks(ax, order, rot=None) -> None:
    """Tool names, upright Arial, rotated by ``rot`` (default TICK_ROT = 90).

    A slot is only 13.0 pt (D/E/F) to 15.8 pt (A2) wide, while the widest label
    ("NanoSPA_m6A", 51.4 pt at 7.5 pt type) needs 42.5 pt of horizontal room at
    45 deg: 13/13 labels overflow their slot and the text collides (measured
    98 colliding pairs), and a clean 45 deg slant would need 2.3 pt type.  The
    angle stays switchable for experiments: ``--tick-rot 45``.
    """
    rot = TICK_ROT if rot is None else rot
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([lab(t) for t in order], rotation=rot,
                       ha="center" if rot >= 60 else "right",
                       rotation_mode="default" if rot >= 60 else "anchor",
                       fontsize=FS_TICK)


def shorten(text: str, width: int = 25) -> str:
    """One-line GO term label; the full string stays in the tables."""
    text = str(text)
    if len(text) <= width:
        return text
    cut = text[:width].rsplit(" ", 1)[0]
    return (cut or text[:width]) + "\u2026"


def style(ax) -> None:
    ax.tick_params(direction="out", length=3, labelsize=FS_TICK, width=0.8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_linewidth(0.9)
    ax.spines["bottom"].set_linewidth(0.9)
    ax.grid(False)
    ax.minorticks_off()


#: Deferred panel letters: (axes, letter, dx, dy) queued while the panels are
#: built and drawn by ``place_letters()`` once every tick label exists.
_LETTERS: list = []


def letter(ax, ch: str, dx=0.0, dy=0.007) -> None:
    """Queue a panel letter; it is drawn at the panel's own top-left corner.

    2026-09-27 (user): the letter used to sit at the top-left of the *axes*
    bounding box, i.e. directly over the y axis, so it read as part of the axis
    rather than of the panel (and E's letter collided with C's x tick labels
    one row above).  It must sit at the top-left of the whole panel, left of the
    y tick labels.  The tick labels are unknown while a panel is being built, so
    the letters are queued here and placed after the figure has been drawn.
    """
    _LETTERS.append((ax, ch, dx, dy))


def place_letters(fig) -> None:
    """Draw the queued letters at each panel's own top-left corner."""
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    for ax, ch, dx, dy in _LETTERS:
        #: tight bbox = axes + tick labels + axis labels (legends below the axes
        #: only push the box downwards, so the top stays the panel's top)
        bb = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(inv)
        fig.text(bb.x0 + dx, bb.y1 + dy, ch, fontsize=FS_PANEL,
                 fontweight="bold", va="bottom", ha="left")


# --------------------------------------------------------------------------- #
def panel_a(fig, rect, tab: Path) -> None:
    """A: the two counts on the shared testable universe (log-log, parity)."""
    ratio = pd.read_csv(tab / "fig2a_testable_ratio.tsv", sep="\t")
    ax1 = fig.add_axes(rect)

    # ---- A1: the two counts on the shared universe (log-log, parity line) ----
    for sp in BLOCKS:
        d = ratio[ratio.species == sp]
        ax1.scatter(d["n_wt_common"], d["n_ctrl_common"], s=7, marker=MARKER[sp],
                    facecolor=COLOR[sp], edgecolor="none", alpha=0.75,
                    label=sp, zorder=3)
    lo, hi = 10.0, 2e5
    ax1.plot([lo, hi], [lo, hi], color="0.35", lw=0.9, ls=(0, (4, 2)), zorder=1)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(lo, hi)
    ax1.set_ylim(lo, hi)
    ax1.set_xlabel("Control-condition calls", fontsize=FS_LABEL)
    ax1.set_ylabel("WT calls in the\ncandidate-site set", fontsize=FS_LABEL)
    ax1.legend(frameon=False, fontsize=FS_LEGEND, loc="upper left",
               handletextpad=0.3, borderpad=0.1, labelspacing=0.25,
               markerscale=1.8)
    style(ax1)
    letter(ax1, "A")


def panel_b_ratio(fig, rect, tab: Path) -> None:
    """B: control/WT ratio per tool, horizontal 13-row layout (2026-09-24: the
    old A2 facet is promoted to its own letter)."""
    ratio = pd.read_csv(tab / "fig2a_testable_ratio.tsv", sep="\t")
    ax = fig.add_axes(rect)
    rng = np.random.default_rng(20260914)
    ords = order_by(ratio.groupby("tool")["ctrl_wt_ratio"].mean(), asc=False)
    n = len(ords)
    for i, sp in enumerate(BLOCKS):
        d = ratio[ratio.species == sp]
        for j, tl in enumerate(ords):
            v = d[d.tool == tl]["ctrl_wt_ratio"].to_numpy()
            if v.size == 0:
                continue
            jitter = (i - 1) * 0.10 + rng.uniform(-0.022, 0.022, size=v.size)
            ax.scatter(v, np.full(v.size, j) + jitter, s=7, marker=MARKER[sp],
                       facecolor=COLOR[sp], edgecolor="none", alpha=0.65, zorder=3)
            ax.scatter([np.nanmean(v)], [j + (i - 1) * 0.10], s=26,
                       marker=MARKER[sp], facecolor=COLOR[sp], edgecolor="white",
                       linewidth=0.5, zorder=4)
    ax.axvline(1.0, color="0.35", lw=0.9, ls=(0, (4, 2)), zorder=1)
    ax.set_xlim(-0.05, 2.9)
    ax.set_xticks([0, 1, 2])
    ax.set_xlabel("Control / WT ratio", fontsize=FS_LABEL, labelpad=2)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_yticks(range(n))
    ax.set_yticklabels([lab(t) for t in ords], fontsize=FS_TICK)
    ax.tick_params(axis="y", length=2.0, pad=2.5)
    handles = [Patch(facecolor=COLOR[sp], edgecolor="white", linewidth=0.3,
                     label=sp) for sp in BLOCKS]
    ax.legend(handles=handles, frameon=False, fontsize=7.0, ncol=3,
              handlelength=0.7, handletextpad=0.25, columnspacing=0.6,
              borderpad=0.0, borderaxespad=0.0, labelspacing=0.15,
              loc="upper left", bbox_to_anchor=(0.45, -0.215))
    style(ax)
    letter(ax, "B")


# --------------------------------------------------------------------------- #
def panel_c_curves(fig, rect, tab: Path) -> None:
    ax = fig.add_axes(rect)
    d = pd.read_csv(tab / "fig2b_support_curves.tsv", sep="\t")
    s = pd.read_csv(tab / "fig2b_support_summary.tsv", sep="\t")
    for sp in BLOCKS:
        for _, u in d[d.species == sp].groupby("sample"):
            ax.plot(u["k"], u["frac_support_ge_k"], color=COLOR[sp], lw=0.7,
                    alpha=0.40, zorder=2)
    for sp in BLOCKS:
        u = s[s.species == sp]
        ax.plot(u["k"], u["frac_mean"], color=COLOR[sp], lw=2.0, zorder=4,
                label=sp, solid_capstyle="round")
    ax.set_yscale("log")
    ax.set_ylim(1e-7, 2.0)
    ax.set_yticks([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
    ax.set_yticklabels(["10$^{-6}$", "10$^{-5}$", "10$^{-4}$", "10$^{-3}$",
                        "10$^{-2}$", "10$^{-1}$", "1"])
    ax.set_xticks(range(2, 14))
    ax.set_xlim(0.7, 13.3)
    #: 2026-09-24: shortened so the centred label cannot reach the bold "D"
    #: letter of the panel below (the letter sits at the panel's left edge)
    ax.set_xlabel("Minimum supporting tools $k$", fontsize=FS_LABEL)
    ax.set_ylabel("Fraction of sites", fontsize=FS_LABEL)
    ax.legend(frameon=False, fontsize=FS_LEGEND, loc="lower left",
              handlelength=1.4, handletextpad=0.4, labelspacing=0.25,
              borderpad=0.1)
    style(ax)
    letter(ax, "C")


# --------------------------------------------------------------------------- #
def panel_d_gobp(fig, rects, tab: Path) -> None:
    for i, sp in enumerate(BLOCKS):
        ax = fig.add_axes(rects[i])
        f = tab / f"fig2c_gobp_{sp}.tsv"
        if not f.exists():
            ax.axis("off")
            continue
        d = pd.read_csv(f, sep="\t").head(5)
        d = d.sort_values("FoldEnrichment", ascending=False).iloc[::-1]
        # the term labels sit *left* of the axis, so they have to fit in the
        # Cx margin; 26 chars is what the 1.46 in column holds at 7.5 pt
        names = [textwrap.fill(t, 26) for t in d["Description"]]
        ax.barh(np.arange(len(d)), d["FoldEnrichment"], height=0.55,
                color=COLOR[sp], edgecolor="none", zorder=3)
        ax.set_yticks(np.arange(len(d)))
        #: 2026-09-24 (user, "text overlap"): with the 34-char wrap the labels
        #: are at most two lines, and the line spacing went 1.05 -> 1.15 so the
        #: two lines of one label no longer read as overprinted type
        ax.set_yticklabels(names, fontsize=FS_TICK, linespacing=1.15)
        ax.set_xlim(0, max(2.0, float(d["FoldEnrichment"].max()) * 1.25))
        #: 2026-09-24 (user, "is D's x axis pressing into the text?"): the three
        #: blocks share one x scale, so only the bottom one prints numbers.  The
        #: upper ones used to put "0" 4.4 pt to the right of the first row label
        #: of the block below, with 2.1 pt of vertical overlap, so the tick read
        #: as part of that label; ticks stay, numbers go.
        if i != len(BLOCKS) - 1:
            ax.tick_params(axis="x", labelbottom=False)
        if i == len(BLOCKS) - 1:          # one shared axis label for the block
            ax.set_xlabel("Fold enrichment", fontsize=FS_LABEL)
        style(ax)
        for t in ax.get_yticklabels():
            t.set_fontweight("normal")
        if i == 0:
            letter(ax, "D")
    
    #: same species patches as B/E, in the service strip under the shared x
    #: label of the bottom block (rects[-1] = the third GO block).
    handles = [Patch(facecolor=COLOR[sp], edgecolor="white", linewidth=0.3,
                     label=sp) for sp in BLOCKS]
    fig.legend(handles=handles, frameon=False, fontsize=7.0, ncol=3,
               handlelength=0.7, handletextpad=0.25, columnspacing=0.6,
               borderpad=0.0, loc="upper left",
               bbox_to_anchor=(rects[-1][0] + 0.037, rects[-1][1] - 0.055))


def order_by(series, asc: bool) -> list:
    """One panels own tool order, best first."""
    vals = {t: float(series.get(t, np.nan)) for t in TOOLS}
    good = [t for t in TOOLS if not np.isnan(vals[t])]
    good.sort(key=lambda t: vals[t], reverse=not asc)
    return good + [t for t in TOOLS if t not in good]


def _mean_panel(fig, rect, df, col, order, ylabel, *, logy=False, ylim=None,
                yticks=None, legend_ncol=0, legend_loc="upper right",
                letter_ch="", horizontal=False, legend_dy=-0.215):
    """Per-unit scatter + per-species mean for one metric.

    ``horizontal=False`` keeps the original vertical layout (tools on x, names
    rotated).  ``horizontal=True`` (2026-09-24, panels D/E) puts the tools on
    y -- a 13-row stack whose names are plain horizontal text, so no tick
    rotation is needed at all; ``ylim``/``yticks`` then describe the value
    axis, which is now x.
    """
    ax = fig.add_axes(rect)
    rng = np.random.default_rng(7)
    n = len(order)
    for i, sp in enumerate(BLOCKS):
        means = []
        for j, tl in enumerate(order):
            v = df[(df.species == sp) & (df.tool == tl)][col].to_numpy()
            v = v[np.isfinite(v)]
            if v.size == 0:
                means.append(np.nan)
                continue
            if horizontal:
                ax.scatter(v, np.full(v.size, j) + rng.uniform(-0.16, 0.16, v.size),
                           s=6, marker=MARKER[sp], facecolor=COLOR[sp],
                           edgecolor="none", alpha=0.55, zorder=3)
            else:
                ax.scatter(np.full(v.size, j) + rng.uniform(-0.16, 0.16, v.size),
                           v, s=6, marker=MARKER[sp], facecolor=COLOR[sp],
                           edgecolor="none", alpha=0.55, zorder=3)
            means.append(float(np.mean(v)))
        if horizontal:
            ax.plot(means, range(n), color=COLOR[sp], lw=1.1, zorder=4, label=sp)
            ax.scatter(means, range(n), s=18, marker=MARKER[sp],
                       facecolor=COLOR[sp], edgecolor="white", linewidth=0.4,
                       zorder=5)
        else:
            ax.plot(range(n), means, color=COLOR[sp], lw=1.1, zorder=4, label=sp)
            ax.scatter(range(n), means, s=18, marker=MARKER[sp],
                       facecolor=COLOR[sp], edgecolor="white", linewidth=0.4,
                       zorder=5)
    if horizontal:
        if logy:
            ax.set_xscale("log")
        if ylim:
            ax.set_xlim(*ylim)
        if yticks is not None:
            ax.set_xticks(yticks[0])
            ax.set_xticklabels(yticks[1])
        ax.set_xlabel(ylabel, fontsize=FS_LABEL, labelpad=2)
        ax.set_ylim(n - 0.5, -0.5)
        ax.set_yticks(range(n))
        ax.set_yticklabels([lab(t) for t in order], fontsize=FS_TICK)
        ax.tick_params(axis="y", length=2.0, pad=2.5)
    else:
        if logy:
            ax.set_yscale("log")
        if ylim:
            ax.set_ylim(*ylim)
        if yticks is not None:
            ax.set_yticks(yticks[0])
            ax.set_yticklabels(yticks[1])
        ax.set_ylabel(ylabel, fontsize=FS_LABEL)
        tool_ticks(ax, order)
    if legend_ncol:
        #: species colour-block legend (2026-09-24): patches, not mean lines;
        #: horizontal panels put it in the empty strip above the axes (the
        #: 13 rows of unit scatter leave no free corner inside either panel)
        handles = [Patch(facecolor=COLOR[sp], edgecolor="white", linewidth=0.3,
                         label=sp) for sp in BLOCKS]
        if horizontal:
            #: 2026-09-24: the legend lives in the 12 pt service strip below
            #: the x label -- inside the panel it sat on the row connectors
            #: (markers are not the only ink), above it sat on the
            #: neighbouring panel's tick numbers
            ax.legend(handles=handles, frameon=False, fontsize=7.0,
                      ncol=3, handlelength=0.7, handletextpad=0.25,
                      columnspacing=0.6, borderpad=0.0, borderaxespad=0.0,
                      labelspacing=0.15, loc="upper left",
                      bbox_to_anchor=(0.12, legend_dy))
        else:
            ax.legend(handles=handles, frameon=False, fontsize=FS_LEGEND,
                      loc=legend_loc, ncol=legend_ncol, handlelength=1.0,
                      handletextpad=0.35, columnspacing=0.9, borderpad=0.1,
                      labelspacing=0.2)
    style(ax)
    if letter_ch:
        letter(ax, letter_ch)
    return ax


def panel_e_jaccard(fig, rect, tab: Path) -> None:
    jac = pd.read_csv(tab / "fig2d_wt_ctrl_jaccard.tsv", sep="\t")
    order = order_by(jac.groupby("tool")["jaccard"].mean(), asc=False)
    _mean_panel(fig, rect, jac, "jaccard", order, "WT vs. control Jaccard",
                ylim=(-0.03, 0.85),
                yticks=([0.2, 0.4, 0.6, 0.8], ["0.2", "0.4", "0.6", "0.8"]),
                legend_ncol=1, legend_loc="lower right",
                letter_ch="E", horizontal=True)


def panel_f_fp(fig, rect, tab: Path) -> None:
    fp = pd.read_csv(tab / "fig2d_ctrl_fp_by_unit.tsv", sep="\t")
    order = order_by(fp.groupby("tool")["fp_per_1e6_candidates"].mean(),
                     asc=True)
    #: 2026-09-24: its long xlabel shares the default legend line, so this key
    #: drops one line lower (the F/G gap has the room)
    _mean_panel(fig, rect, fp, "fp_per_1e6_candidates", order,
                "Control-side FP per 10$^{6}$ candidates", logy=True,
                legend_dy=-0.268,
                ylim=(1.0, 3e4),
                yticks=([10, 100, 1000, 10000],
                        ["10", "100", "1000", "10$^{4}$"]),
                legend_ncol=1, legend_loc="upper right", letter_ch="F",
                horizontal=True)


def panel_g_replicate(fig, rect, tab: Path) -> None:
    """G: replicate consistency per tool, horizontal 13-row layout."""
    ax = fig.add_axes(rect)
    d = pd.read_csv(tab / "fig2e_replicate_consistency.tsv", sep="\t")
    mean = d[d.metric == "mean_pairwise"][["species", "tool", "value"]]
    #: 2026-09-28 (union drop): the pooled (global) series was withdrawn from
    #: the text and the caption; only the per-unit-pair mean is drawn here.
    order = order_by(mean.groupby("tool")["value"].mean(), asc=False)
    n = len(order)
    for i, sp in enumerate(BLOCKS):
        m = mean[mean.species == sp].set_index("tool")["value"]
        ys, hi = [], []
        for j, tl in enumerate(order):
            if tl not in m.index:
                continue
            ys.append(j)
            hi.append(m[tl])
        ax.scatter(hi, ys, s=16, marker="o", facecolor=COLOR[sp],
                   edgecolor="white", linewidth=0.4, zorder=4)
    ax.set_xlim(0, 0.92)
    ax.set_xticks([0, 0.25, 0.5, 0.75])
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_yticks(range(n))
    ax.set_yticklabels([lab(t) for t in order], fontsize=FS_TICK)
    ax.tick_params(axis="y", length=2.0, pad=2.5)
    ax.set_xlabel("Jaccard between independent units", fontsize=FS_LABEL,
                  labelpad=2)
    handles = [plt.Line2D([], [], marker="o", ls="none", markersize=4.5,
                          markerfacecolor="0.45", markeredgecolor="white",
                          label="mean pairwise")]
    
    #: colours every point by species (COLOR[sp]) but its key names the value
    #: and the species
    handles += [plt.Line2D([], [], marker=MARKER[sp], ls="none", markersize=4.5,
                           markerfacecolor=COLOR[sp], markeredgecolor="white",
                           label=sp) for sp in BLOCKS]
    #: 2026-09-28 (union drop): with the pooled series gone the key has four
    #: entries; it keeps the lower-right corner, whose rows stop below x = 0.2
    #: (``--qa`` reports 0 data points covered by the legend box).
    ax.legend(handles=handles, frameon=False, fontsize=7.0, ncol=2,
              loc="lower right", bbox_to_anchor=(0.995, 0.020),
              handlelength=0.9, handletextpad=0.35, columnspacing=0.9,
              borderpad=0.0)
    style(ax)
    letter(ax, "G")


# --------------------------------------------------------------------------- #
def main() -> None:
    global TICK_ROT
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", default=str(REV / "tables"))
    ap.add_argument("--out", default=str(REV / "figures" / "Figure2_rev"))
    ap.add_argument("--qa", action="store_true")
    ap.add_argument("--tick-rot", type=int, default=TICK_ROT,
                    help="rotation of the tool-name tick labels (default 90)")
    args = ap.parse_args()
    TICK_ROT = args.tick_rot
    tab, out = Path(args.tables), Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    figstyle.apply()
    plt.rcParams.update({
        "font.size": FS_TICK,
        # mathtext must be pinned to Arial: the default (DejaVu) leaks a second
        # font family into the PDF wherever a tick label carries $...$
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "mathtext.default": "regular",
    })
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))
    # Explicit rectangles: with bbox_inches=None the canvas IS the page, so every
    # label must fit inside it (a GridSpec crop would silently cut panel C's
    # term labels).  [x0, y0, width, height] in figure fractions.
    L = 0.150      # left column axes
    Wl = 0.395
    R = 0.665      # right column axes
    Wr = 0.325
    
    #: start 0.045 closer to the edge while the right edge stays on 0.545 like
    #: A and B.  The label column is 1.46 in wide: measured at 300 dpi a 7.5 pt
    #: Arial line runs ~3.7 pt per char, so the 26-char wrap needs ~96 pt and
    #: still clears the 10 pt page margin (the 34-char attempt was clipped).
    Cx = 0.240
    Wc = 0.305

    
    #: every tool-named panel is horizontal (names on y, no tick rotation);
    #: letters A-G run down the left column (A, B, D) and the right (C, E, F, G)
    ax_a = [L, 0.845, Wl, 0.130]
    ax_b = [L, 0.601, Wl, 0.190]
    #: 2026-09-27 (user, fourth pass): the four right-column panels are spaced
    #: evenly now.  The first attempt (0.845 / 0.622 / 0.363) left C's x label and
    #: E's bold letter in the same band (~0 pt clear), the second (0.850 / 0.530)
    #: opened 43 pt, which read as too airy; this spacing leaves ~20 pt between
    #: C's x label and E's letter and ~12 pt between the other pairs.
    ax_c = [R, 0.860, Wr, 0.115]
    ax_e = [R, 0.584, Wr, 0.185]
    ax_f = [R, 0.330, Wr, 0.185]
    
    #: full service band above the bottom margin
    ax_g = [R, 0.082, Wr, 0.190]
    
    #: (the old strip left ~85 pt of dead space below D)
    ax_d1 = [Cx, 0.395, Wc, 0.137]
    ax_d2 = [Cx, 0.246, Wc, 0.137]
    ax_d3 = [Cx, 0.097, Wc, 0.137]

    panel_a(fig, ax_a, tab)
    panel_b_ratio(fig, ax_b, tab)
    panel_c_curves(fig, ax_c, tab)
    panel_d_gobp(fig, [ax_d1, ax_d2, ax_d3], tab)
    panel_e_jaccard(fig, ax_e, tab)
    panel_f_fp(fig, ax_f, tab)
    panel_g_replicate(fig, ax_g, tab)

    
    #: at each panel's own top-left corner (see place_letters()).
    place_letters(fig)

    # NOT figstyle.save(): its bbox_inches="tight" crop changes the page size,
    # and the page size IS the contract here -- the manuscript includes the PDF
    # at 0.95\textwidth, so any crop would scale the type below 7 pt on paper.
    if args.qa:
        import _qa
        _qa.qa_tick_overlaps(fig)
        _qa.qa_overlaps(fig, out)
    fig.savefig(str(out) + ".pdf", bbox_inches=None)
    fig.savefig(str(out) + ".png", bbox_inches=None, dpi=300)
    plt.close(fig)
    print(f"wrote {out}.pdf / {out}.png  (page {PAGE_W} x {PAGE_H} in)")


if __name__ == "__main__":
    main()
