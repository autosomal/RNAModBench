#!/usr/bin/env python3
"""R3-9b -- figure for the depth-vs-replicate diagnosis of the RNA002 non-m6A
WT vs unmodified-IVT comparison.

House rules (src/harmonisation/common/figstyle.py): no figure or panel
titles, no grid lines, no in-figure numbers or callouts -- values live in the
evidence tables; only axis labels, ticks, legends and bold panel letters.
Vector PDF (editable text, fonttype 42) + PNG, 180 mm double-column canvas.

Panels
------
A  per-site call rate against coverage, per tool, WT solid/filled vs IVT
   open/dashed -> the coverage-driven detection floor shared by both conditions
B  IVT/WT ratio from raw counts, from per-site density, and after coverage-band
   matching (bootstrap CI over the three units per condition)
C  median coverage of the candidate sites of every sequencing unit on a log axis
   -> the synthetic Curlcake null sits two to three orders of magnitude deeper
   than any HeLa library, so its density per candidate is not on the same footing

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    analysis/nonm6a_replicate_depth/analysis/rd3_figure.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve()
PKG = HERE.parents[1]
PROJECT = HERE.parents[3]
sys.path.insert(0, str(PROJECT / "src/harmonisation/common"))
import figstyle  # noqa: E402

EV = PKG / "evidence"
C_WT = "#1F3B64"
C_IVT = "#E08A2E"
C_CHANCE = "#8A8A8A"
STYLE = dict(font_size=12.5, label_size=15.0, tick_size=12.0, legend_size=10.5,
             letter_size=18.0)
TOOL_LABEL = {
    "NanoNm": "NanoNm", "NanoMUD_psi": "NanoMUD-\u03a8", "NanoPsu": "NanoPsu",
    "NanoSPA_psU": "NanoSPA-\u03a8", "NanoMUD_m1psi": "NanoMUD-m1\u03a8",
    "CHEUI_m5C": "CHEUI-m5C",
}
ORDER = ["CHEUI_m5C", "NanoMUD_m1psi", "NanoMUD_psi", "NanoNm",
         "NanoPsu", "NanoSPA_psU"]
TOOL_COLOR = {
    "CHEUI_m5C": "#B3261E", "NanoMUD_m1psi": "#E08A2E", "NanoMUD_psi": "#6C8EBF",
    "NanoNm": "#2F6B4F", "NanoPsu": "#7A4FA3", "NanoSPA_psU": "#4A4A4A",
}
BANDS = ["10-19", "20-29", "30-49", "50-99", "100-199", "200-399", "400-99999"]
BAND_MID = np.sqrt(np.array([10 * 20, 20 * 30, 30 * 50, 50 * 100, 100 * 200,
                             200 * 400, 400 * 1000], dtype=float))


def apply_style() -> None:
    figstyle.apply()
    mpl.rcParams.update({
        "font.size": STYLE["font_size"], "axes.labelsize": STYLE["label_size"],
        "xtick.labelsize": STYLE["tick_size"], "ytick.labelsize": STYLE["tick_size"],
        "legend.fontsize": STYLE["legend_size"], "axes.grid": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 1.2, "xtick.major.width": 1.2, "ytick.major.width": 1.2,
        "xtick.major.size": 4.5, "ytick.major.size": 4.5,
        "pdf.fonttype": 42, "ps.fonttype": 42, "figure.dpi": 120})


def letter(ax, s, dx=-0.14, dy=1.05):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=STYLE["letter_size"],
            fontweight="bold", va="bottom", ha="left")


def panel_a(ax, bands: pd.DataFrame) -> None:
    # rate per 1e6 sites at risk, pooled over the three units of each condition
    strat = pd.read_csv(EV / "cov_strata.tsv", sep="\t")
    strat = strat[strat.dataset == "HeLa"]
    for tool in ORDER:
        g = strat[strat.tool == tool]
        piv_c = g.pivot_table(index="unit", columns="stratum", values="calls",
                              aggfunc="sum", fill_value=0).reindex(columns=BANDS, fill_value=0)
        piv_s = g.pivot_table(index="unit", columns="stratum", values="sites_at_risk",
                              aggfunc="sum", fill_value=0).reindex(columns=BANDS, fill_value=0)
        for cond, units, ls, mk in (("WT", ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"], "-", "o"),
                                    ("IVT", [f"HeLa_IVT_rep{i}" for i in (1, 2, 3)], "--", "s")):
            c = piv_c.loc[units].to_numpy().sum(axis=0).astype(float)
            s = piv_s.loc[units].to_numpy().sum(axis=0).astype(float)
            rate = np.where(s > 0, c / np.maximum(s, 1) * 1e6, np.nan)
            ok = np.isfinite(rate) & (rate > 0)
            ax.plot(BAND_MID[ok], rate[ok], ls, marker=mk, ms=5.5, lw=1.5,
                    mfc=TOOL_COLOR[tool] if cond == "WT" else "white",
                    mec=TOOL_COLOR[tool], mew=1.4, label=None)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(8, 1200)
    ax.set_ylim(0.3, 1.6e6)
    ax.set_xlabel("coverage of candidate sites")
    ax.set_ylabel("calls per $10^{6}$ sites at risk")
    handles = [Line2D([], [], color=TOOL_COLOR[t], lw=1.6, label=TOOL_LABEL[t])
               for t in ORDER]
    handles += [Line2D([], [], color=C_CHANCE, lw=1.5, marker="o", ms=5.5,
                       label="wild type"),
                Line2D([], [], color=C_CHANCE, lw=1.5, ls="--", marker="s", ms=5.0,
                       mfc="white", label="unmodified IVT")]
    ax.legend(handles=handles, loc="upper left", frameon=False, ncol=1,
              handletextpad=0.4, labelspacing=0.28, columnspacing=0.9,
              fontsize=8.5, borderaxespad=0.15)
    letter(ax, "A")


def panel_b(ax, band: pd.DataFrame) -> None:
    y = np.arange(len(ORDER))[::-1]
    for i, tool in enumerate(ORDER):
        r = band[band.tool == tool].iloc[0]
        yy = y[i]
        raw = float(r.raw_density_ratio)
        mat = float(r.band_matched_rate_ratio)
        lo, hi = float(r.boot_matched_ci_lo), float(r.boot_matched_ci_hi)
        ax.plot([lo, hi], [yy, yy], color=C_CHANCE, lw=1.4, zorder=1)
        ax.plot([mat, mat], [yy - 0.20, yy + 0.20], color=C_IVT, lw=3.0, zorder=3)
        ax.plot([raw], [yy], marker="o", ms=7.5, mfc=C_WT, mec=C_WT, zorder=4)
        ax.plot([mat], [yy], marker="D", ms=6.5, mfc=C_IVT, mec=C_IVT, zorder=4)
    ax.axvline(1.0, color="black", lw=1.1, ls=":", zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels([TOOL_LABEL[t] for t in ORDER])
    ax.set_ylim(-0.7, len(ORDER) - 0.3)
    ax.tick_params(axis="y", labelsize=11.0)
    ax.set_xlabel("unmodified-IVT / wild-type ratio")
    ax.set_ylabel("")
    ax.set_xscale("log")
    handles = [Line2D([], [], marker="o", ls="", color=C_WT, ms=7.5,
                      label="per-site density (all depths)"),
               Line2D([], [], marker="D", ls="", color=C_IVT, ms=6.5,
                      label="coverage-band matched"),
               Line2D([], [], color=C_CHANCE, lw=1.4, label="95% CI (unit bootstrap)")]
    ax.figure._rd1_b_keys = handles
    letter(ax, "B")


def panel_c(ax, bands: pd.DataFrame) -> None:
    """Share of candidate sites per coverage band, WT vs IVT (composition shift)."""
    strat = pd.read_csv(EV / "cov_strata.tsv", sep="\t")
    strat = strat[strat.dataset == "HeLa"]
    tot = strat.pivot_table(index="unit", columns="stratum", values="sites_at_risk",
                            aggfunc="sum", fill_value=0).reindex(columns=BANDS, fill_value=0)
    wt = ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"]
    iv = [f"HeLa_IVT_rep{i}" for i in (1, 2, 3)]
    share_wt = (tot.loc[wt].sum() / tot.loc[wt].sum().sum() * 100).reindex(BANDS)
    share_iv = (tot.loc[iv].sum() / tot.loc[iv].sum().sum() * 100).reindex(BANDS)
    x = np.arange(len(BANDS))
    ax.bar(x - 0.2, share_wt.to_numpy(), width=0.38, color=C_WT, label="wild type")
    ax.bar(x + 0.2, share_iv.to_numpy(), width=0.38, color=C_IVT,
           label="unmodified IVT")
    ax.set_xticks(x)
    ax.set_xticklabels(BANDS, rotation=45, ha="right", fontsize=9.5)
    ax.set_xlabel("coverage of candidate sites")
    ax.set_ylabel("share of candidate sites (%)")
    ax.legend(loc="upper right", frameon=False, fontsize=9.5)
    letter(ax, "C")


def panel_d(ax, units: pd.DataFrame) -> None:
    h = units[units.dataset == "HeLa"].copy()
    h["order"] = np.where(h.condition_class == "WT", 0, 1)
    c = units[units.dataset == "Curlcake"].copy()
    groups = [("HeLa WT", h[h.order == 0], C_WT),
              ("HeLa IVT", h[h.order == 1], C_IVT),
              ("Curlcake IVT", c, "#2F6B4F")]
    pos, labels, keys = [], [], []
    for name, sub, col in groups:
        for _, r in sub.sort_values("median_cov_all_ge10").iterrows():
            pos.append(len(pos))
            labels.append(str(r.unit))
            keys.append(col)
            ax.plot(len(pos) - 1, float(r.median_cov_all_ge10), marker="o",
                    ms=8, color=col, zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(0.7, 1e6)
    ax.set_xticks(pos)
    ax.set_xticklabels(labels, rotation=90, fontsize=8.0)
    ax.set_ylabel("median coverage of\ncandidate sites")
    for x in (2.5, 5.5):
        ax.axvline(x, color=C_CHANCE, lw=0.9)
    handles = [Line2D([], [], marker="o", ls="", color=col, ms=8, label=name)
               for name, _, col in groups]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=9.0)
    letter(ax, "D")


def main() -> None:
    apply_style()
    band = pd.read_csv(EV / "band_matched_rate_ratios.tsv", sep="\t")
    units = pd.read_csv(EV / "units_depth_all.tsv", sep="\t")
    fig, axes = plt.subplots(2, 2, figsize=(18.0 / 2.54, 11.6 / 2.54),
                             gridspec_kw={"width_ratios": [1.15, 1.0]})
    panel_a(axes[0, 0], band)
    panel_b(axes[0, 1], band)
    panel_c(axes[1, 0], band)
    panel_d(axes[1, 1], units)
    keys = getattr(axes[0, 1], "_rd1_b_keys", [])
    axes[0, 1].legend(handles=keys, loc="lower right", frameon=False, fontsize=8.5,
                      borderaxespad=0.25)
    fig.subplots_adjust(left=0.115, right=0.985, top=0.935, bottom=0.115,
                        wspace=0.62, hspace=0.50)
    figstyle.save(fig, str(EV / "FigR3_9b_depth_vs_replicate"))
    print("written", EV / "FigR3_9b_depth_vs_replicate.pdf")


if __name__ == "__main__":
    main()
