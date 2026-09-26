#!/usr/bin/env python3
"""45 -- SUPERSEDED (2026-09-20): old six-panel A-F assembly of Figure 7.

Kept for history only; the delivered Figure 7 is now the four-row A-D
rebuild in ``47_fig7_rebuild.py`` (paper narrative order).

Original docstring: Figure 7 (revised): non-m6A tool evaluation, six panels A-F.

Panels A-E reuse the panel functions of the finalised R3-9 analysis figure
(``R3-9_nonm6a_fp_analysis/analysis/r39_figure.py``) so the main figure and the
response-letter figure stay consistent and read the same evidence tables
(``evidence/*.tsv`` -- the single source of numbers).  Panel F is new: the
replicate-aware non-m6A metagene rebuilt from ``harmonisation/callsets`` with the
Ensembl region model (script 43), drawn in the FigS1 manner -- majority
consensus thick line + fill, individual replicates as thin dashed lines.

House rules: no grid, no in-figure value callouts; bold panel letters only;
vector PDF + 300/600 dpi PNG.

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    src/harmonisation/scripts/45_fig7_assembled.py
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
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve()
PROJECT = _RB
sys.path.insert(0, str((_RB / "src/harmonisation/common")))
sys.path.insert(0, str((_RB / "analysis/nonm6a_false_positives/analysis")))

import figstyle                        # noqa: E402
import r39_figure as r39               # noqa: E402  (panels A-E, palette, ORDER)

OUT = (_RB / "figures/figure7")
TAB = (_RB / "figures/figure7/tables")
FIG = (_RB / "figures/figure7/figures")

GRID = 200
#: panel-F tool order -- grouped by modification (Psi, m1Psi, Nm, m5C)
F_ORDER = ["NanoPsu", "NanoSPA_psU", "NanoMUD_psi", "NanoMUD_m1psi", "NanoNm", "CHEUI_m5C"]
WT_BLUE = "#4b81b8"
IVT_ORANGE = "#e8a76b"
COND_OF = {"WT": "HeLa_WT", "IVT": "HeLa_IVT"}
CONDITION_COLOR = {"WT": WT_BLUE, "IVT": IVT_ORANGE}
SEGS = ["five_prime_flank", "five_prime_UTR", "CDS", "three_prime_UTR",
        "three_prime_flank"]


def load_density() -> pd.DataFrame:
    d = pd.read_csv((_RB / "figures/figure7/tables/fig7_metagene_density.tsv"), sep="\t")
    d = d[d["merge"].eq("majority") | d["merge"].str.startswith("unit:")].copy()
    d["curve"] = [np.fromiter(map(float, s.split("|")), float, GRID) for s in d.density]
    return d


def spine(d: pd.DataFrame, condition: str, tool: str, merge: str) -> np.ndarray | None:
    out = np.zeros(len(SEGS) * GRID)
    seen = False
    for i, name in enumerate(SEGS):
        r = d[(d.condition == condition) & (d.tool == tool)
              & (d["merge"] == merge) & (d.kind == name)]
        if len(r):
            out[i * GRID:(i + 1) * GRID] = r.curve.iloc[0]
            seen = True
    return out if seen else None


def guitar_small(ax, title: str, *, tick_fs: float = 8.5, title_fs: float = 11.0) -> None:
    """Compact GUITAR frame (segment labels, dotted separators, schematic bar)."""
    figstyle.guitar_panel(ax, SEGS, title)
    ax.set_title(title, fontweight="bold", fontsize=title_fs, pad=3)
    ax.set_xticklabels([figstyle.SEGMENT_LABELS[s] for s in SEGS], fontsize=tick_fs)
    ax.tick_params(length=2.5, labelsize=tick_fs)
    ax.set_ylabel("")


def panel_f(f_axes, d: pd.DataFrame) -> None:
    """Replicate-aware non-m6A metagene: majority consensus + per-replicate lines."""
    for ax, tool in zip(f_axes.ravel(), F_ORDER):
        guitar_small(ax, r39.TOOL_LABEL.get(tool, tool))
        have = False
        tops = []
        for cond in ("WT", "IVT"):
            colour = CONDITION_COLOR[cond]
            main = spine(d, COND_OF[cond], tool, "majority")
            if main is None:
                continue
            have = True
            tops.append(float(main.max()))
            x = (np.arange(len(main)) + 0.5) / GRID
            reps = sorted(m.split(":", 1)[1]
                          for m in d.loc[d.condition == COND_OF[cond], "merge"]
                          if m.startswith("unit:"))
            for u in reps:
                c = spine(d, COND_OF[cond], tool, f"unit:{u}")
                if c is not None:
                    ax.plot(x, c, color=colour, lw=0.8, alpha=0.55,
                            ls=(0, (2, 1.5)), zorder=2)
            ax.fill_between(x, 0, main, color=colour, alpha=0.16, lw=0, zorder=3)
            ax.plot(x, main, color=colour, lw=1.9, zorder=4)
        if not have:
            ax.text(0.5, 0.5, "no calls", transform=ax.transAxes, ha="center",
                    va="center", fontsize=9, color="0.4")
        ax.set_ylim(0, (max(tops) if tops else 1.0) * 1.12)
    for ax in f_axes[0]:                     # top row: hide segment labels
        ax.set_xticklabels([])


def main() -> None:
    figstyle.apply()
    r39.apply_style()
    # mathtext (e.g. the 10^6 in panel B) must render in Arial too, otherwise
    # matplotlib silently falls back to DejaVu Sans inside the vector PDF
    plt.rcParams.update({
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
    })
    data = r39.load()
    dens = load_density()

    fig = plt.figure(figsize=(7.087, 9.9))            # 180 mm wide
    gs = fig.add_gridspec(
        7, 6,
        height_ratios=[0.92, 1.02, 0.58, 0.58, 0.78, 0.90, 0.90],
        hspace=1.12, wspace=1.05,
        left=0.085, right=0.985, top=0.972, bottom=0.048)

    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0:2])
    ax_c1 = fig.add_subplot(gs[1, 2:4])
    ax_c2 = fig.add_subplot(gs[1, 4:6])
    d_axes = np.array([[fig.add_subplot(gs[2, 0:2]), fig.add_subplot(gs[2, 2:4]),
                        fig.add_subplot(gs[2, 4:6])],
                       [fig.add_subplot(gs[3, 0:2]), fig.add_subplot(gs[3, 2:4]),
                        fig.add_subplot(gs[3, 4:6])]])
    ax_e1 = fig.add_subplot(gs[4, 0:3])
    ax_e2 = fig.add_subplot(gs[4, 3:6])
    f_axes = np.array([[fig.add_subplot(gs[5, 0:2]), fig.add_subplot(gs[5, 2:4]),
                        fig.add_subplot(gs[5, 4:6])],
                       [fig.add_subplot(gs[6, 0:2]), fig.add_subplot(gs[6, 2:4]),
                        fig.add_subplot(gs[6, 4:6])]])

    r39.panel_a(ax_a, data["truth"])
    r39.panel_b(ax_b, data["cc"])
    # headroom for the B legend: extend the log axis so all bars sit in the
    # lower 60% and the key can occupy the empty top-left band
    _ind = data["cc"].loc[data["cc"].construct_role == "independent",
                          "fp_per_1e6_candidates"].to_numpy(float)
    _floor = float(_ind[_ind > 0].min()) * 0.12
    r39.fake_zero_axis(ax_b, _floor, [1e2, 1e4, 1e6])
    ax_b.set_ylim(_floor * 0.55, 1e8)
    r39.panel_c1(ax_c1, data["rep"])
    r39.panel_c2(ax_c2, data["sum"])
    r39.panel_d(d_axes, data["hist"], data)
    r39.panel_e([ax_e1, ax_e2], data["eco"])
    panel_f(f_axes, dens)

    # share the tool axis across B/C1/C2: labels only under the rightmost panel
    for ax in (ax_b, ax_c1):
        ax.set_xticklabels([])
    # legend placement fixes (avoid covering data)
    old = ax_b.get_legend()
    if old is not None:
        old.remove()
    ax_b.legend(handles=[Line2D([], [], marker="o", linestyle="none", markersize=8,
                                markerfacecolor=r39.C_HL, markeredgecolor=r39.C_HL,
                                label="Curlcake IVT"),
                         Line2D([], [], marker="o", linestyle="none", markersize=8,
                                markerfacecolor="white", markeredgecolor=r39.C_HL,
                                markeredgewidth=1.4, label="depth-matched subset"),
                         Line2D([], [], color=r39.C_WT, linewidth=2.4,
                                label="mean (independent)")],
                loc="upper left", bbox_to_anchor=(0.02, 0.99), frameon=False,
                ncol=1, handletextpad=0.3,
                borderpad=0.05, labelspacing=0.10, columnspacing=0.9, fontsize=8.5)
    old = ax_c2.get_legend()
    if old is not None:
        old.remove()
    ax_c2.legend(handles=[Line2D([], [], marker="o", linestyle="none", markersize=8,
                                 markerfacecolor=r39.C_WT, markeredgecolor=r39.C_WT,
                                 label="WT"),
                          Line2D([], [], marker="o", linestyle="none", markersize=8,
                                 markerfacecolor="white", markeredgecolor=r39.C_IVT,
                                 markeredgewidth=1.4, label="IVT (unmodified)"),
                          Line2D([], [], marker="s", linestyle="none", markersize=7,
                                 markerfacecolor=r39.C_CHANCE, markeredgecolor=r39.C_CHANCE,
                                 label="mean pairwise")],
                 loc="upper right", bbox_to_anchor=(0.99, 1.0), frameon=False,
                 ncol=1, handletextpad=0.5,
                 borderpad=0.1, labelspacing=0.22, fontsize=8.5)
    ax_c2.set_ylim(-0.02, 0.78)

    # shared keys (legends only, no value callouts)
    d_axes[0, 0].legend(
        handles=[Line2D([], [], color=r39.C_WT, lw=6, alpha=0.5, label="WT"),
                 Line2D([], [], color=r39.C_IVT, lw=2.2, label="IVT (unmodified)")],
        loc="upper right", frameon=False, handletextpad=0.5, borderpad=0.2,
        labelspacing=0.3, fontsize=10)
    ax_e1.set_ylabel("density", fontsize=13)
    f_axes[0, 0].set_ylabel("Density", fontsize=12)
    f_axes[1, 0].set_ylabel("Density", fontsize=12)

    r39.letter(ax_a, "A", dx=0.0, dy=1.14)
    r39.letter(ax_b, "B", dy=1.30)
    r39.letter(ax_c1, "C")
    r39.letter(d_axes[0, 0], "D", dx=-0.02, dy=1.30)
    r39.letter(ax_e1, "E", dx=-0.02, dy=1.22)
    r39.letter(f_axes[0, 0], "F", dx=-0.02, dy=1.34)

    FIG.mkdir(parents=True, exist_ok=True)
    stem = (_RB / "figures/figure7/figures/Figure7_rev")
    for suffix, dpi in (("", 300), ("_600dpi", 600)):
        fig.savefig(f"{stem}{suffix}.pdf", bbox_inches="tight", pad_inches=0.02)
        fig.savefig(f"{stem}{suffix}.png", bbox_inches="tight",
                    pad_inches=0.02, dpi=dpi)
    plt.close(fig)
    print("wrote", stem, "(pdf + png 300/600 dpi)")


if __name__ == "__main__":
    main()
