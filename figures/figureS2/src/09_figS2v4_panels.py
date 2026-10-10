#!/usr/bin/env python
"""Figure S2 v3 - per-panel drawing.

Each panel is drawn with the SAME functions twice: once alone (origin at
(0, 0), so the panel PDF/PNG can be inspected or replaced on its own) and once
inside the assembled page (origin at the panel's slot) - see
``07_figS2_assemble.py``.  Panel geometry is therefore expressed in page
inches with a top-left origin, and every panel receives its own offset.

Panels
  A  ten tools (rows, two groups of five) x three species (columns) 5-mer
     sequence logos; one shared 0-2 bits axis per group, labelled
     "Information content (bits)" (R3-m4).  Rows 0.42 in since 2026-10-09 (was
     0.36), so that the tool name and its KL value clear each other.
  B  top-5 preferred 5-mers, four species blocks (2 x 2), single-hue blue
     scale, RRACH 5-mers in bold orange, shared horizontal colour bar.
  C  cross-species AGAC vs GGAC preference (R3-6): per tool, three species
     markers at freq(AGAC) - freq(GGAC); positive = plant-type AGAC, negative
     = animal-type GGAC; note on the right states that the tool factor
     dominates the variance (95.3% vs 1.0%).

Outputs figures/panelA_logos.{pdf,png}, panelB_top5.{pdf,png},
        panelC_agac_ggac.{pdf,png}
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
import os
import sys
from pathlib import Path

import logomaker
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TAB, FIG = ROOT / "analysis", ROOT / "figures"
sys.path.insert(0, str(_RB / "src/harmonisation"))
from common.figstyle import apply as apply_style  # noqa: E402

PAGE_W, PAGE_H = 595.276 / 72.0, 841.89 / 72.0
MIN_PT = 7.0
FS = {"panel": 14.0, "tool": 9.5, "kl": 9.0, "species": 9.0, "block": 10.0,
      "rank": 9.0, "kmer": 9.5, "cbar": 9.0, "ctick": 9.0, "axis": 9.0,
      "tick": 9.0, "note": 9.0}
INK = {"dark": "#272727", "text": "#1A1A1A", "muted": "#6B6B6B",
       "rule": "#C9C9C9", "rule_soft": "#EFEFEF", "blue": "#4b81b8",
       "blue_dark": "#2c5f8a", "orange": "#e8a76b"}
FREQ_CMAP = LinearSegmentedColormap.from_list(
    "s2_blue", ["#FFFFFF", "#DCE9F4", "#9DBEDC", "#4b81b8", "#2c5f8a"])
VMAX = 0.45
WHITE_TEXT_NORM = 0.78
LOGO_COLOR = {"A": "#33A02C", "C": "#1F78B4", "G": "#E58606", "U": "#7F8C99"}
BASES = ["A", "U", "C", "G"]
SP_A = ["Arabidopsis", "Mouse", "Human"]
SP_COLOR = {"Arabidopsis": INK["blue"], "Mouse": INK["orange"],
            "Human": "#4D4D4D"}
#: shape as well as colour, so the three species stay apart in panel D (and in
#: a grey print) -- colour alone was the one thing the panel's key was missing
SP_MARKER = {"Arabidopsis": "o", "Mouse": "s", "Human": "^"}
ORDER10 = ["CHEUI_m6A", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo", "EpiNano_Error",
           "NanoSPA_m6A", "Nanocompore", "Yanocomp", "m6Anet", "xPore"]
ORDER13 = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
           "EpiNano_Error", "MINES", "NanoSPA_m6A", "Nanocompore",
           "Nanom6A", "Yanocomp", "m6Anet", "xPore"]
SP_DISP = {"Arabidopsis": "Arabidopsis", "Mouse": "Mouse", "Human": "Human",
           "E.coli": "E. coli"}
RANKS = ["1st", "2nd", "3rd", "4th", "5th"]
# standardized nomenclature (Table S11): ELIGOS2_diff / ELIGOS2_solo everywhere.
# The frozen tables still carry the legacy display names, so normalize on load.
RENAME = {"ELIGOS_diff": "ELIGOS2_diff", "ELIGOS_solo": "ELIGOS2_solo"}
DISPLAY = {"yanocomp": "Yanocomp"}



#: budget is re-cut instead -- A and B grow, the C/D/E row is compressed, and the
#: two defects the re-layout exposed are fixed with it:
#:   * A's right block printed its last logo axes at 8.355 in on an 8.268 in page
#:     (the Human column was clipped at the page edge; 433 ink pixels sat in the
#:     last 12 px column of the 300 dpi render).  ``A_CELL`` 0.96 -> 0.91 and
#:     ``A_GAP`` 0.20 -> 0.19 bring the right edge to 8.055 in.
#:   * every row of A stacked its tool name and its "KL = x" line 7.2 pt apart
#:     while the two text boxes are 10.6 / 10 pt tall, so all 30 word pairs
#:     overlapped by 3.1-3.2 pt (measured with pdftotext -bbox).  The taller rows
#:     (``A_ROW`` 0.360 -> 0.410) let the two lines sit 10.9 pt apart.
#: Panel D's 13-tool row pitch was hardcoded as 0.180 in two places; it is one
#: constant now (``D_ROW``), because ``panel_heights`` has to follow it.
#: The freed space is split A +0.25 / B +0.29 in of band height, paid for by the
#: C/D/E band (-0.55 in): A 2.36 -> 2.61, B 5.25 -> 5.54, C/D/E 3.33 -> 2.78.
#: The species key under the C/D/E band stays where it was (the sheet still ends
#: 0.42 in of clear space below the last band).
A_LAB, A_CELL, A_GAP, A_HEAD, A_ROW = 0.90, 0.91, 0.19, 0.28, 0.410
A_PAD = 0.28
#: the two label lines of an A row: name below, KL value above, 10.9 pt apart
#: (their boxes are 10.6 / 10 pt tall, so this is what clears them; at the old
#: spacing of 7.2 pt all 30 pairs overlapped)
A_NAME_DY, A_KL_DY = -0.076, 0.076
B_LAB, B_CELL, B_GAP = 0.78, 0.55, 0.26
B_TITLE, B_HEAD, B_ROW = 0.16, 0.13, 0.1545
B_PAD, B_BLOCK_GAP, B_CBAR_PAD = 0.30, 0.10, 0.18
C_LAB, C_ROW = 0.95, 0.130
C_PAD, C_BOT = 0.34, 0.40
D_ROW = 0.160         #: row pitch of panel D's 13 tools (was a hardcoded 0.180)
PANEL_BOX = 13 * D_ROW    #: common content box height for panels C, D and E
X0, X1 = 0.45, PAGE_W - 0.30
#: the right edge every band has to stay inside (2026-10-09): the C/D/E row ends
#: at 8.10 in and A's cells do too now, so the page keeps a real right margin
RIGHT_LIMIT = 8.10


def page_axes(fig, w=PAGE_W, h=PAGE_H):
    """A helper axes covering a page-sized box, with a top-left origin."""
    pax = fig.add_axes([0, 0, 1, 1], zorder=0)
    pax.set_xlim(0, w)
    pax.set_ylim(h, 0)
    pax.axis("off")
    return pax


def fig_top(ax_top: float, h: float, page_h: float = PAGE_H) -> float:
    return 1.0 - (ax_top + h) / page_h


def bits_matrix(pwm: np.ndarray) -> pd.DataFrame:
    h = -(pwm * np.log2(np.clip(pwm, 1e-12, None))).sum(axis=1)
    ic = np.clip(2.0 - h, 0.0, None)
    return pd.DataFrame(pwm * ic[:, None], index=range(-2, 3), columns=BASES)


def is_rrach(k: str) -> bool:
    return (len(k) == 5 and k[2] == "A" and k[3] == "C" and k[0] in "AG"
            and k[1] in "AG" and k[4] in "ACU")


def _standard_names(df):
    """Map the legacy Tool_display values onto the standardized nomenclature."""
    if "Tool_display" in df.columns:
        df = df.copy()
        df["Tool_display"] = df["Tool_display"].replace(RENAME)
    return df


def load_inputs():
    pooled = _standard_names(pd.read_csv(TAB / "figS2_kl_pooled.tsv", sep="\t"))
    top5 = _standard_names(pd.read_csv(TAB / "figS2_top5_full.tsv", sep="\t"))
    pwm = _standard_names(pd.read_csv(TAB / "figS2_pwm_mean.tsv", sep="\t"))
    agac = _standard_names(pd.read_csv(TAB / "figS2_panelC_agac.tsv", sep="\t"))
    return pooled, top5, pwm, agac


def draw_a(fig, pooled, pwm, org=(0.0, 0.0), page_w=PAGE_W, page_h=PAGE_H,
           panel_w=PAGE_W):
    """Panel A: logos, ten tools as rows in two groups of five."""
    x_org, y_org = org
    pax = page_axes(fig, page_w, page_h)
    klx = {(t.Tool_display, t.Species): t.KL_mean for t in pooled.itertuples()}
    top = y_org + A_PAD
    y_head = top + A_HEAD
    for g in (0, 1):
        x_lab = x_org + X0 + g * (A_LAB + 3 * A_CELL + A_GAP)
        x_cells = x_lab + A_LAB + 0.17
        tools = ORDER10[g * 5:(g + 1) * 5]
        for j, sp in enumerate(SP_A):
            pax.text(x_cells + (j + 0.5) * A_CELL, top + 0.10, sp,
                     ha="center", va="center", fontsize=FS["species"],
                     color=INK["dark"])
        for j in range(1, 3):
            xj = x_cells + j * A_CELL
            pax.plot([xj, xj], [y_head, y_head + 5 * A_ROW],
                     color=INK["rule_soft"], lw=0.4, zorder=1)
        ax_x = x_cells - 0.085
        pax.plot([ax_x, ax_x], [y_head, y_head + 5 * A_ROW], color=INK["muted"],
                 lw=0.6)
        for v in (1, 2):
            yv = y_head + 5 * A_ROW - (v / 2.0) * (5 * A_ROW)
            pax.plot([ax_x - 0.03, ax_x], [yv, yv], color=INK["muted"], lw=0.6)
        if g == 0:
            pax.text(x_org + 0.18, y_head + 2.5 * A_ROW,
                     "Information content (bits)", rotation=90, ha="center",
                     va="center", fontsize=FS["axis"], color=INK["muted"])
        for i, tool in enumerate(tools):
            yc = y_head + (i + 0.5) * A_ROW
            if i % 2 == 1:                      # light row banding for guidance
                x_band = x_lab - 0.04
                pax.add_patch(Rectangle((x_band, y_head + i * A_ROW),
                                        x_cells + 3 * A_CELL - x_band, A_ROW,
                                        facecolor="#FAFAFA", edgecolor="none",
                                        zorder=0))
            kbar = float(np.nanmean([klx.get((tool, s), np.nan) for s in SP_A]))
            pax.text(x_lab, yc + A_NAME_DY, tool, ha="left", va="center",
                     fontsize=FS["tool"], fontweight="bold", color=INK["dark"])
            pax.text(x_lab, yc + A_KL_DY, f"KL = {kbar:.1f}", ha="left",
                     va="center", fontsize=FS["kl"], color=INK["blue_dark"])
            for j, sp in enumerate(SP_A):
                sub = pwm[(pwm.Tool_display == tool) & (pwm.Species == sp)]
                sub = sub.sort_values("Pos")
                axx = fig.add_axes([(x_cells + j * A_CELL + 0.025) / page_w,
                                    fig_top(yc - 0.5 * A_ROW + 0.03,
                                            A_ROW - 0.06, page_h),
                                    (A_CELL - 0.05) / page_w,
                                    (A_ROW - 0.06) / page_h])
                logo = logomaker.Logo(
                    bits_matrix(sub[BASES].to_numpy(dtype=float)), ax=axx,
                    color_scheme=LOGO_COLOR)
                logo.ax.set_ylim(0, 2.45)
                logo.ax.set_xlim(-2.45, 2.45)
                logo.ax.set_xticks([])
                logo.ax.set_yticks([])
                for sp_ in logo.ax.spines.values():
                    sp_.set_visible(False)
    pax.text(x_org + X0, y_org + 0.10, "A", fontsize=FS["panel"],
             fontweight="bold", color=INK["dark"], ha="left", va="center")


def draw_b(fig, top5, org=(0.0, 0.0), page_w=PAGE_W, page_h=PAGE_H,
           panel_w=PAGE_W):
    """Panel B: top-5 5-mers in four species blocks (2 x 2) + colour bar."""
    x_org, y_org = org
    pax = page_axes(fig, page_w, page_h)
    y0 = y_org + B_PAD
    row0_h = B_TITLE + B_HEAD + 13 * B_ROW
    blocks = [("Arabidopsis", 0, 0), ("Mouse", 0, 1),
              ("Human", 1, 0), ("E.coli", 1, 1)]
    for sp, r, c in blocks:
        bx = x_org + X0 + c * (B_LAB + 5 * B_CELL + B_GAP)
        by = y0 + r * (row0_h + 0.14)
        sub = top5[top5.Species == sp]
        tools = [t for t in ORDER13 if t in set(sub.Tool_display)]
        pax.text(bx + (B_LAB + 5 * B_CELL) / 2, by + B_TITLE / 2, SP_DISP[sp],
                 ha="center", va="center", fontsize=FS["block"],
                 fontweight="bold", color=INK["dark"])
        cells_x = bx + B_LAB
        for j, rk in enumerate(RANKS):
            pax.text(cells_x + (j + 0.5) * B_CELL, by + B_TITLE + B_HEAD / 2,
                     rk, ha="center", va="center", fontsize=FS["rank"],
                     color=INK["muted"])
        pax.plot([cells_x, cells_x + 5 * B_CELL],
                 [by + B_TITLE + B_HEAD] * 2, color=INK["rule"], lw=0.4)
        for i, tool in enumerate(tools):
            yc = by + B_TITLE + B_HEAD + (i + 0.5) * B_ROW
            pax.text(cells_x - 0.06, yc, DISPLAY.get(tool, tool), ha="right",
                     va="center", fontsize=FS["kmer"], color=INK["dark"])
            row = sub[sub.Tool_display == tool].set_index("Rank")
            for j in range(1, 6):
                km = row.loc[j, "Kmer"]
                f = float(row.loc[j, "RelFreq_repmix"])
                norm_v = min(f / VMAX, 1.0)
                axc = fig.add_axes([(cells_x + (j - 1) * B_CELL + 0.008) / page_w,
                                    fig_top(yc - B_ROW / 2 + 0.006,
                                            B_ROW - 0.012, page_h),
                                    (B_CELL - 0.016) / page_w,
                                    (B_ROW - 0.012) / page_h])
                axc.axis("off")
                axc.add_patch(Rectangle((0, 0), 1, 1, transform=axc.transAxes,
                                        facecolor=FREQ_CMAP(norm_v),
                                        edgecolor="none"))
                rr = is_rrach(km)
                txtcol = ("white" if norm_v > WHITE_TEXT_NORM
                          else (INK["orange"] if rr else INK["text"]))
                axc.text(0.5, 0.5, km, transform=axc.transAxes, ha="center",
                         va="center", fontsize=FS["kmer"], color=txtcol,
                         fontweight="bold" if rr else "normal")
    # colour bar + direct legend, below the bottom block row
    yb = y0 + 2 * row0_h + B_BLOCK_GAP + 0.30
    cax = fig.add_axes([(x_org + 1.75) / page_w, fig_top(yb, 0.10, page_h),
                        1.7 / page_w, 0.10 / page_h])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, VMAX),
                                            cmap=FREQ_CMAP), cax=cax,
                      orientation="horizontal", ticks=[0.0, 0.2, 0.4])
    cb.ax.tick_params(labelsize=FS["ctick"], length=2, width=0.4, pad=1.5,
                      color=INK["muted"], labelcolor=INK["muted"])
    cb.outline.set_visible(False)
    pax.text(x_org + X0, yb + 0.05, "Relative frequency", ha="left",
             va="center", fontsize=FS["cbar"], color=INK["dark"])
    # legend: two sample cells with plain labels, generous spacing
    lx, ly, lw, lh = x_org + X0 + 3.60, yb - 0.02, 0.42, 0.16
    pax.add_patch(Rectangle((lx, ly), lw, lh, facecolor=FREQ_CMAP(0.30),
                            edgecolor="none"))
    pax.text(lx + lw / 2, ly + lh / 2, "GGACU", ha="center", va="center",
             fontsize=FS["kmer"], color=INK["orange"], fontweight="bold")
    pax.text(lx + lw + 0.10, ly + lh / 2, "RRACH", ha="left", va="center",
             fontsize=FS["note"], color=INK["dark"])
    lx2 = lx + lw + 0.72
    pax.add_patch(Rectangle((lx2, ly), lw, lh, facecolor="white",
                            edgecolor=INK["rule"], lw=0.5))
    pax.text(lx2 + lw / 2, ly + lh / 2, "UGUCA", ha="center", va="center",
             fontsize=FS["kmer"], color=INK["text"])
    pax.text(lx2 + lw + 0.10, ly + lh / 2, "other", ha="left", va="center",
             fontsize=FS["note"], color=INK["dark"])
    pax.text(x_org + X0, y_org + 0.10, "B", fontsize=FS["panel"],
             fontweight="bold", color=INK["dark"], ha="left", va="center")


def draw_d(fig, agac, org=(0.0, 0.0), page_w=PAGE_W, page_h=PAGE_H,
           panel_w=PAGE_W):
    """Panel C: cross-species AGAC vs GGAC preference (R3-6)."""
    x_org, y_org = org
    pax = page_axes(fig, page_w, page_h)
    top = y_org + C_PAD
    row0 = top
    order = (agac.pivot(index="Tool_display", columns="Species", values="delta")
             .sort_values("Arabidopsis", ascending=False))
    x_left, x_right = x_org + 1.35, x_org + panel_w - 0.15
    lo, hi = -0.42, 0.35

    def xmap(v):
        return x_left + (v - lo) / (hi - lo) * (x_right - x_left)

    # axis + zero line
    y_bot = row0 + len(order) * D_ROW
    pax.plot([x_left, x_right], [y_bot + 0.02] * 2, color=INK["muted"], lw=0.6)
    pax.plot([xmap(0)] * 2, [top - 0.04, y_bot + 0.02], color=INK["muted"],
             lw=0.6)
    for v in (-0.4, -0.2, 0.0, 0.2):
        pax.plot([xmap(v)] * 2, [y_bot + 0.02, y_bot + 0.07], color=INK["muted"],
                 lw=0.6)
        #: 2026-10-09: the tick row and the axis title sit 0.08 / 0.25 in below the
        #: axis so that both 9 pt lines fit C_BOT = 0.40 in without touching
        #: (0.08-0.227 and 0.25-0.397 in; the old 0.12 / 0.28 pair fitted 0.42 in)
        pax.text(xmap(v), y_bot + 0.08, f"{v:+.1f}", ha="center", va="top",
                 fontsize=FS["tick"], color=INK["muted"])
    pax.text((x_left + x_right) / 2, y_bot + 0.25,
             "freq(AGAC) − freq(GGAC)", ha="center", va="top",
             fontsize=FS["axis"], color=INK["dark"])
    for i, (tool, row) in enumerate(order.iterrows()):
        yc = row0 + (i + 0.5) * D_ROW
        pax.text(x_left - 0.20, yc, tool, ha="right", va="center",
                 fontsize=FS["tool"], color=INK["dark"])
        for k, sp in enumerate(("Arabidopsis", "Mouse", "Human")):
            pax.scatter([xmap(float(row[sp]))], [yc + (k - 1) * 0.045], s=19,
                        color=SP_COLOR[sp], marker=SP_MARKER[sp], zorder=3,
                        edgecolor="white", linewidth=0.4)
    pax.text(x_org + 0.10, y_org + 0.10, "D", fontsize=FS["panel"],
             fontweight="bold", color=INK["dark"], ha="left", va="center")


def draw_c(fig, over, summary, org=(0.0, 0.0), page_w=PAGE_W, page_h=PAGE_H,
           panel_w=PAGE_W):
    """Panel C: top-5 identity - within-species tools vs same tool across species."""
    x_org, y_org = org
    pax = page_axes(fig, page_w, page_h)
    top = y_org + 0.34
    d = over[over.scope == summary.scope.iloc[0]]
    groups = [("different\ntools", "within species, different tools",
               INK["blue"]),
              ("same tool\ndifferent species", "same tool, different species",
               INK["orange"])]
    ax = fig.add_axes([(x_org + 1.05) / page_w, fig_top(top, PANEL_BOX, page_h),
                       (panel_w - 1.35) / page_w, PANEL_BOX / page_h])
    rng = np.random.default_rng(20260921)
    for gi, (glabel, gkey, gcol) in enumerate(groups):
        sub = d[d.group == gkey]
        y = sub.overlap.to_numpy(float) + rng.uniform(-0.16, 0.16, len(sub))
        ax.scatter(np.full(len(y), gi) + rng.uniform(-0.16, 0.16, len(y)), y,
                   s=7, color=gcol, alpha=0.55, linewidths=0)
        med = float(np.median(sub.overlap))
        ax.plot([gi - 0.28, gi + 0.28], [med, med], color=INK["dark"], lw=1.4,
                zorder=5)
        ax.text(gi, 5.55, f"median {med:.0f}", ha="center", va="bottom",
                fontsize=FS["tick"], color=INK["dark"])
        ax.text(gi, -0.85, f"n = {len(sub)}", ha="center", va="top",
                fontsize=FS["tick"], color=INK["muted"])
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(-1.4, 6.0)
    ax.set_xticks([0, 1])
    #: 2026-10-09: the two lines of these labels met (0.6 pt of overlap in the
    #: printed page, measured with pdftotext -bbox); 1.45 linespacing separates
    #: them without moving the axis
    ax.set_xticklabels(["different\ntools", "same\ntool"], fontsize=8.0,
                       linespacing=1.45, color=INK["dark"])
    ax.set_yticks([0, 1, 2, 3, 4, 5])
    ax.set_ylabel("Top-5 5-mer overlap", fontsize=FS["axis"], labelpad=1.5)
    ax.tick_params(labelsize=FS["tick"], length=2, pad=1.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    p = float(summary.p_one_sided_less.iloc[0])
    pax.text(x_org + 0.10, y_org + 0.10, "C", fontsize=FS["panel"],
             fontweight="bold", color=INK["dark"], ha="left", va="center")


def draw_e(fig, org=(0.0, 0.0), page_w=PAGE_W, page_h=PAGE_H, panel_w=PAGE_W):
    """Panel E: two-way variance decomposition of the KL divergence."""
    x_org, y_org = org
    pax = page_axes(fig, page_w, page_h)
    top = y_org + 0.34
    var = pd.read_csv(TAB / "figS2v4_variance.tsv", sep="\t")
    var = var[var.scope == "cell means, panel-A 10 tools, three species"]
    eta = {r.Source: r.eta_sq * 100 for r in var.itertuples()}
    names = ["Tool", "Species", "Residual"]
    vals = [eta["Tool"], eta["Species"], eta["Interaction+Residual"]]
    cols = [INK["blue"], INK["orange"], INK["rule"]]
    ax = fig.add_axes([(x_org + 1.05) / page_w, fig_top(top, PANEL_BOX, page_h),
                       (panel_w - 1.30) / page_w, PANEL_BOX / page_h])
    yy = np.arange(3)[::-1]
    ax.barh(yy, vals, height=0.55, color=cols, edgecolor="none")
    for v, yv in zip(vals, yy):
        ax.text(v + 2.0, yv, f"{v:.0f}%", va="center", ha="left",
                fontsize=FS["axis"], color=INK["dark"])
    ax.set_yticks(yy)
    ax.set_yticklabels(names, fontsize=FS["axis"], color=INK["dark"])
    ax.set_xlim(0, 112)
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel("KL variance (%)", fontsize=FS["axis"], labelpad=1.5)
    ax.tick_params(labelsize=FS["tick"], length=2, pad=1.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    pax.text(x_org + 0.10, y_org + 0.10, "E", fontsize=FS["panel"],
             fontweight="bold", color=INK["dark"], ha="left", va="center")


def text_objects(fig):
    """Every printed text object of the page (axes texts, tick labels, labels).

    ``ax.texts`` always counts: an axes switched off with ``axis("off")`` still
    prints the annotations drawn into it (the panel letters, the species key).
    Tick labels only count where the axes keeps its frame on -- the sheet-wide
    helper axes of ``page_axes()`` has default ticks at 0,1,2,... that are never
    drawn, and counting them made the audit report text outside the page.
    """
    out = []
    for ax in fig.axes:
        out += [t for t in ax.texts if t.get_text() and t.get_visible()]
        if ax.axison:
            out += [t for t in ax.get_xticklabels() + ax.get_yticklabels()
                    if t.get_text() and t.get_visible()]
            for lbl in (ax.xaxis.label, ax.yaxis.label, ax.title):
                if lbl.get_text() and lbl.get_visible():
                    out.append(lbl)
    out += [t for t in fig.texts if t.get_text() and t.get_visible()]
    return out


def assert_fonts(fig, min_pt=MIN_PT):
    sizes = [float(t.get_fontsize()) for t in text_objects(fig)]
    smallest = float(min(sizes))
    assert smallest >= min_pt, f"text {smallest:.2f} pt < {min_pt} pt"
    return smallest


def assert_page_geometry(fig, page_w=PAGE_W, page_h=PAGE_H, min_pt=MIN_PT,
                         overlap_tol=1.0, edge_tol=0.5, right_limit=RIGHT_LIMIT):
    """Printed-size audit: fonts, page overflow, and text-on-text collisions.

    2026-10-09: the re-layout of S2 (A and B taller, C/D/E shorter) was the first
    pass that could tighten two pitches that had never been checked, and the two
    defects it exposed were only visible in the *printed* page:
      * A's right block ran 0.087 in past the page edge, so the Human column of
        the second group was clipped (433 ink pixels in the last 12 px column of
        the 300 dpi render);
      * the tool name and the "KL = x" line of every A row overlapped by 3.1-3.2
        pt (all 30 word pairs), and panel C's two two-line tick labels by 0.6 pt.
    Both are now asserted against the drawn objects: every text box inside the
    page, no pair overlapping more than ``overlap_tol`` pt in either axis, and
    every axes/rectangle inside the page and clear of ``right_limit``.  The check
    runs on the same figure the PDF is written from, so a later edit cannot put
    either defect back without failing the build.
    """
    smallest = assert_fonts(fig, min_pt)
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    dpi = float(fig.dpi)
    texts = text_objects(fig)
    boxes = []
    for t in texts:
        bb = t.get_window_extent(renderer=rend)
        if bb.width > 0 and bb.height > 0:
            boxes.append((str(t.get_text())[:24], bb))
    outside = [(s, round(bb.x0 / dpi, 3), round(bb.x1 / dpi, 3),
                round(bb.y0 / dpi, 3), round((page_h * dpi - bb.y1) / dpi, 3))
               for s, bb in boxes
               if bb.x0 < edge_tol or bb.x1 > page_w * dpi - edge_tol
               or bb.y0 < edge_tol or bb.y1 > page_h * dpi - edge_tol]
    assert not outside, f"text outside the page (x0, x1, y_bottom, y_top): {outside[:6]}"
    hits = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            dx = min(a.x1, b.x1) - max(a.x0, b.x0)
            dy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if dx > overlap_tol and dy > overlap_tol:
                hits.append((boxes[i][0], boxes[j][0], round(dx / dpi * 72, 2),
                             round(dy / dpi * 72, 2)))
    assert not hits, (f"{len(hits)} text pairs overlap by more than "
                      f"{overlap_tol} pt (name, name, dx pt, dy pt): {hits[:6]}")
    over = []
    for ax in fig.axes:
        bb = ax.get_window_extent(renderer=rend)
        if bb.x1 > page_w * dpi + edge_tol or bb.x0 < -edge_tol \
                or bb.y0 < -edge_tol or bb.y1 > page_h * dpi + edge_tol:
            over.append((round(bb.x0 / dpi, 3), round(bb.x1 / dpi, 3),
                         round(bb.y0 / dpi, 3), round(bb.y1 / dpi, 3)))
    assert not over, f"axes outside the page: {over}"
    #: the widest inked *element* is what ran off the page (a logo axes ending at
    #: 8.355 in on an 8.268 in page, clipping the Human column of A's second
    #: group).  Every drawing axes narrower than the page is measured, plus every
    #: inked patch/collection (A's row-banding strips, B's sample cells) except the
    #: sheet-wide helper axes and its background, which span the page by design.
    full = 0.99 * page_w * dpi
    right = 0.0
    for ax in fig.axes:
        bb = ax.get_window_extent(renderer=rend)
        if bb.width < full:
            right = max(right, bb.x1 / dpi)
        for art in list(ax.patches) + list(ax.collections):
            if not art.get_visible():
                continue
            ab = art.get_window_extent(renderer=rend)
            if 0 < ab.width < full:
                right = max(right, ab.x1 / dpi)
    assert right <= right_limit + 0.02, \
        (f"the right-most drawing box ends at {right:.3f} in, past the "
         f"{right_limit} in right limit of this page")
    return {"n_text": len(boxes), "min_pt": smallest, "overlaps": len(hits),
            "max_element_x1": round(right, 3)}


def panel_heights():
    row0_h = B_TITLE + B_HEAD + 13 * B_ROW
    return {"A": A_PAD + A_HEAD + 5 * A_ROW,
            "B": B_PAD + 2 * row0_h + B_BLOCK_GAP + 0.30 + 0.20,
            "C": C_PAD + PANEL_BOX + C_BOT, "D": C_PAD + 13 * D_ROW + C_BOT,
            "E": C_PAD + PANEL_BOX + C_BOT}


def render_panel(name: str, draw, data, height: float, page_w=PAGE_W,
                 width: float | None = None, skirt: float = 0.0):
    """Draw one panel alone on its own canvas.

    ``height`` is the panel's band height on the assembled page; ``skirt`` adds
    canvas below it for the furniture that is allowed to hang into the gap
    between bands there.  The B colour bar's tick labels live 0.04 in below the
    band edge, so without a skirt the delivered ``panelB_top5.pdf`` used to have
    them cut off by the canvas edge (287 ink pixels on the last 300 dpi row).
    """
    page_w = width or page_w
    canvas_h = height + skirt
    fig = plt.figure(figsize=(page_w, canvas_h))
    draw(fig, *data, org=(0.0, 0.0), page_w=page_w, page_h=canvas_h,
         panel_w=page_w)
    rep = assert_page_geometry(fig, page_w=page_w, page_h=canvas_h,
                               right_limit=min(RIGHT_LIMIT, page_w - 0.02))
    smallest = rep["min_pt"]
    out = FIG / f"panel{name}"
    fig.savefig(str(out) + ".pdf", bbox_inches=None)
    fig.savefig(str(out) + ".png", dpi=300, bbox_inches=None)
    plt.close(fig)
    print(f"panel{name}: {page_w:.3f} x {height:.3f} in, min text "
          f"{smallest:.2f} pt, texts {rep['n_text']} (overlaps "
          f"{rep['overlaps']}, right-most box {rep['max_element_x1']:.3f} in)"
          f" -> {out}.pdf/.png")


def load_v4():
    over = pd.read_csv(TAB / "figS2v4_overlap.tsv", sep="\t")
    summ = pd.read_csv(TAB / "figS2v4_overlap_summary.tsv", sep="\t")
    return over, summ


def main():
    apply_style()
    FIG.mkdir(parents=True, exist_ok=True)
    pooled, top5, pwm, agac = load_inputs()
    over, summ = load_v4()
    h = panel_heights()
    print("panel heights (in):", {k: round(v, 3) for k, v in h.items()})
    render_panel("A_logos", draw_a, (pooled, pwm), h["A"])
    render_panel("B_top5", draw_b, (top5,), h["B"], skirt=0.12)
    render_panel("C_identity", draw_c, (over, summ), h["C"], width=PAGE_W / 2 - 0.15)
    render_panel("D_agac", draw_d, (agac,), h["D"], width=3.90)
    render_panel("E_variance", draw_e, (), h["E"], width=2.55)


if __name__ == "__main__":
    main()
