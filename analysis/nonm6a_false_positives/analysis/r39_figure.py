#!/usr/bin/env python3
"""R3-9 -- main figure for the non-m6A false-positive analysis (2026-09-19).

House rules enforced here (see ``01_code/code/sites_v2/common/figstyle.py``):
no figure/panel titles, no grid lines, no in-figure annotations or value
callouts -- numbers live in the evidence tables; only axis labels, ticks,
legends and bold panel letters are allowed.  Fonts are as large as the 180 mm
(double column) canvas allows.  Vector PDF (editable text) + PNG.

Panels
------
A  truth-anchored specificity: RMBase + DirectRMDB / orthogonal NGS precision, WT vs IVT
B  unmodified Curlcake controls: FP density per 10^6 candidate sites
C  HeLa WT vs IVT per replicate: call counts (left) and Jaccard (right)
D  reported score/stoichiometry density, WT vs IVT (six tool-mod pairs)
E  third-party GSE271571 (E. coli): CHEUI probability density, WT vs IVT

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    04_revision_analysis/R3-9_nonm6a_fp_analysis/analysis/r39_figure.py
"""

from __future__ import annotations

import argparse
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
sys.path.insert(0, str(PROJECT / "src/sites_v2/common"))

import figstyle  # noqa: E402

EV = PKG / "evidence"

#: house palette (restrained): one neutral, one signal, one accent
C_WT = "#1F3B64"
C_IVT = "#E08A2E"
C_CHANCE = "#8A8A8A"
C_HL = "#B3261E"
C_III = "#6C8EBF"

TOOL_LABEL = {
    "NanoNm": "NanoNm",
    "NanoMUD_psi": "NanoMUD-Ψ",
    "NanoPsu": "NanoPsu",
    "NanoSPA_psU": "NanoSPA-Ψ",
    "NanoMUD_m1psi": "NanoMUD-m1Ψ",
    "CHEUI_m5C": "CHEUI-m5C",
}
ORDER = ["NanoNm", "NanoMUD_psi", "NanoPsu", "NanoSPA_psU", "NanoMUD_m1psi", "CHEUI_m5C"]
MOD_OF = {"NanoNm": "Nm", "NanoMUD_psi": "Psi", "NanoPsu": "Psi",
          "NanoSPA_psU": "Psi", "NanoMUD_m1psi": "m1Psi", "CHEUI_m5C": "m5C"}

STYLE = dict(font_size=13.0, label_size=15.5, tick_size=12.5, legend_size=12.0,
             letter_size=19.0)


def apply_style() -> None:
    figstyle.apply()
    mpl.rcParams.update({
        "font.size": STYLE["font_size"],
        "axes.labelsize": STYLE["label_size"],
        "axes.titlesize": STYLE["label_size"],
        "xtick.labelsize": STYLE["tick_size"],
        "ytick.labelsize": STYLE["tick_size"],
        "legend.fontsize": STYLE["legend_size"],
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "xtick.major.size": 4.5,
        "ytick.major.size": 4.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 120,
    })


def letter(ax, s: str, dx: float = -0.16, dy: float = 1.06,
           ha: str = "left") -> None:
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=STYLE["letter_size"],
            fontweight="bold", va="bottom", ha=ha)


def fake_zero_axis(ax, floor: float, ticks: list[float]) -> None:
    """Log axis with an explicit zero row drawn at ``floor``."""
    ax.set_yscale("log")
    ax.set_ylim(floor * 0.55, 1.6)
    ax.set_yticks(ticks + [floor])
    ax.set_yticklabels([f"{t:g}" for t in ticks] + ["0"])


def jitter(n: int, width: float = 0.055) -> np.ndarray:
    return np.linspace(-width, width, n) if n > 1 else np.zeros(1)


def load() -> dict[str, pd.DataFrame]:
    return {
        "cc": pd.read_csv(EV / "curlcake_ivt_fp_per_construct.tsv", sep="\t"),
        "rep": pd.read_csv(EV / "hela_wt_ivt_per_replicate.tsv", sep="\t"),
        "sum": pd.read_csv(EV / "hela_wt_ivt_summary.tsv", sep="\t"),
        "truth": pd.read_csv(EV / "truth_precision_per_replicate.tsv", sep="\t"),
        "scor": pd.read_csv(EV / "score_distributions.tsv", sep="\t"),
        "hist": pd.read_csv(EV / "score_histograms.tsv", sep="\t"),
        "eco": pd.read_csv(EV / "ecoli_gse271571_prob.tsv", sep="\t"),
    }


# --------------------------------------------------------------------------- #
def panel_a(ax, truth: pd.DataFrame) -> None:
    """Truth-anchored specificity: overlap with ORCA / NGS vs the permutation null.

    The y-axis is the enrichment over the chromosome-stratified permutation
    expectation, i.e. a currency that is comparable across references (each
    reference has its own, much smaller, chance overlap).  A value of 1 means
    "indistinguishable from random candidate sites"; the unmodified IVT library
    must sit at 1 for a tool whose calls carry modification-specific information.
    """
    t = truth[(truth["window_bp"] == 1)
              & (truth["reference"].isin(["RMBase+DirectRMDB", "NGS"]))]
    tools = [tl for tl in ORDER if tl in set(t["tool"])]
    pos = {tl: i for i, tl in enumerate(tools)}
    vals = t["enrichment"].to_numpy(float)
    pos_vals = vals[(vals > 0) & np.isfinite(vals)]
    floor = float(pos_vals.min()) * 0.25 if len(pos_vals) else 0.01
    markers = {"RMBase+DirectRMDB": "o", "NGS": "s"}

    for tl in tools:
        x = pos[tl]
        for ref, marker in markers.items():
            sub = t[(t["tool"] == tl) & (t["reference"] == ref)]
            if not len(sub):
                continue
            for grp, off, filled, colour in (("WT", -0.15, True, C_WT),
                                             ("IVT", 0.15, False, C_IVT)):
                rows = sub[sub["group"] == grp]
                if not len(rows):
                    continue
                xs = x + off + jitter(len(rows), 0.045)
                ys = rows["enrichment"].to_numpy(float)
                ys = np.where(np.isfinite(ys) & (ys > 0), ys, floor)
                ax.plot(xs, ys, linestyle="none", marker=marker, markersize=7.2,
                        markerfacecolor=colour if filled else "white",
                        markeredgecolor=colour, markeredgewidth=1.4, zorder=3)
    ax.axhline(1.0, color=C_CHANCE, linewidth=1.6, linestyle=(0, (2, 2)), zorder=1)
    ax.set_yscale("log")
    ax.set_ylim(floor * 0.7, 2000)
    ax.set_yticks([floor, 1, 10, 100, 1000])
    ax.set_yticklabels(["0", "1", "10", "100", "1000"])
    ax.set_xticks(range(len(tools)))
    ax.set_xticklabels([TOOL_LABEL[tl] for tl in tools], fontsize=12.5)
    ax.set_ylabel("enrichment\nover chance")
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=8,
                      markerfacecolor=C_WT, markeredgecolor=C_WT, label="WT"),
               Line2D([], [], marker="o", linestyle="none", markersize=8,
                      markerfacecolor="white", markeredgecolor=C_IVT, markeredgewidth=1.4,
                      label="IVT"),
               Line2D([], [], marker="o", linestyle="none", markersize=7,
                      markerfacecolor="none", markeredgecolor="0.25",
                      label="RMBase + DirectRMDB"),
               Line2D([], [], marker="s", linestyle="none", markersize=7,
                      markerfacecolor="none", markeredgecolor="0.25", label="NGS gold standard")]
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, 1.005),
              frameon=False, ncol=4, handletextpad=0.4, borderpad=0.1,
              labelspacing=0.3, columnspacing=1.5)


def panel_b(ax, cc: pd.DataFrame) -> None:
    tools = [tl for tl in ORDER if tl in set(cc["tool"])]
    pos = {tl: i for i, tl in enumerate(tools)}
    vals = cc.loc[cc["construct_role"] == "independent", "fp_per_1e6_candidates"].to_numpy(float)
    pos_vals = vals[vals > 0]
    floor = float(pos_vals.min()) * 0.12 if len(pos_vals) else 1e-6
    for tl in tools:
        x = pos[tl]
        sub = cc[cc["tool"] == tl]
        for role, filled in (("independent", True), ("depth_matched_subset_of_rep3", False)):
            rows = sub[sub["construct_role"] == role]
            xs = x + jitter(len(rows), 0.13)
            ys = np.where(rows["fp_per_1e6_candidates"].to_numpy(float) > 0,
                          rows["fp_per_1e6_candidates"].to_numpy(float), floor)
            ax.plot(xs, ys, linestyle="none", marker="o", markersize=7.4,
                    markerfacecolor=C_HL if filled else "white",
                    markeredgecolor=C_HL, markeredgewidth=1.4, zorder=3)
        ind = sub[sub["construct_role"] == "independent"]["fp_per_1e6_candidates"].to_numpy(float)
        if ind.size:
            mean = float(np.mean(ind))
            ax.plot([x - 0.26, x + 0.26], [max(mean, floor)] * 2, color=C_WT,
                    linewidth=2.4, solid_capstyle="butt", zorder=4)
    ax.set_xticks(range(len(tools)))
    ax.set_xticklabels([TOOL_LABEL[tl] for tl in tools], rotation=45,
                       ha="right", rotation_mode="anchor", fontsize=11)
    fake_zero_axis(ax, floor, [1e2, 1e3, 1e4, 1e5])
    ax.set_ylim(floor * 0.55, 2e6)
    ax.set_ylabel("FP per 10$^6$\ncandidates")
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=8,
                      markerfacecolor=C_HL, markeredgecolor=C_HL, label="Curlcake IVT"),
               Line2D([], [], marker="o", linestyle="none", markersize=8,
                      markerfacecolor="white", markeredgecolor=C_HL, markeredgewidth=1.4,
                      label="depth-matched subset"),
               Line2D([], [], color=C_WT, linewidth=2.4,                       label="mean (independent)")]
    ax.legend(handles=handles, loc="upper left", frameon=False, ncol=1,
              handletextpad=0.3, borderpad=0.1,
              labelspacing=0.2, columnspacing=0.9, fontsize=9.5)


def panel_c1(ax, rep: pd.DataFrame) -> None:
    tools = [tl for tl in ORDER if tl in set(rep["tool"])]
    pos = {tl: i for i, tl in enumerate(tools)}
    for tl in tools:
        sub = rep[rep["tool"] == tl]
        x = pos[tl]
        for grp, off, filled, colour in (("WT", -0.16, True, C_WT), ("IVT", 0.16, False, C_IVT)):
            rows = sub[sub["group"] == grp].sort_values("sample")
            xs = x + off + jitter(len(rows), 0.05)
            ax.plot(xs, rows["n_calls_in_universe"].to_numpy(float), linestyle="none",
                    marker="o", markersize=7.4, markerfacecolor=colour if filled else "white",
                    markeredgecolor=colour, markeredgewidth=1.4, zorder=3)
    ax.set_yscale("log")
    ax.set_xticks(range(len(tools)))
    ax.set_xticklabels([TOOL_LABEL[tl] for tl in tools], rotation=45,
                       ha="right", rotation_mode="anchor", fontsize=11)
    ax.set_ylabel("calls")


def panel_c2(ax, summ: pd.DataFrame) -> None:
    tools = [tl for tl in ORDER if tl in set(summ["tool"])]
    pos = {tl: i for i, tl in enumerate(tools)}
    for tl in tools:
        row = summ[summ["tool"] == tl].iloc[0]
        x = pos[tl]
        for grp, off, filled, colour in (("wt", -0.16, True, C_WT), ("ivt", 0.16, False, C_IVT)):
            ax.plot([x + off], [row[f"{grp}_global_jaccard_uni"]], linestyle="none",
                    marker="o", markersize=7.4,
                    markerfacecolor=colour if filled else "white",
                    markeredgecolor=colour, markeredgewidth=1.4, zorder=3)
            ax.plot([x + off], [row[f"{grp}_mean_pairwise_jaccard_uni"]], linestyle="none",
                    marker="s", markersize=6.6,
                    markerfacecolor=colour if filled else "white",
                    markeredgecolor=colour, markeredgewidth=1.4, zorder=3)
        ax.plot([x - 0.16, x + 0.16],
                [row["wt_global_jaccard_uni"], row["ivt_global_jaccard_uni"]],
                color=C_CHANCE, linewidth=1.2, zorder=1)
    ax.set_ylim(-0.02, 0.52)
    ax.set_xticks(range(len(tools)))
    ax.set_xticklabels([TOOL_LABEL[tl] for tl in tools], rotation=45,
                       ha="right", rotation_mode="anchor", fontsize=11)
    ax.set_ylabel("Jaccard")
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=8,
                      markerfacecolor=C_WT, markeredgecolor=C_WT, label="WT"),
               Line2D([], [], marker="o", linestyle="none", markersize=8,
                      markerfacecolor="white", markeredgecolor=C_IVT, markeredgewidth=1.4,
                      label="IVT"),
               Line2D([], [], marker="s", linestyle="none", markersize=7,
                      markerfacecolor=C_CHANCE, markeredgecolor=C_CHANCE, label="mean pairwise")]
    ax.legend(handles=handles, loc="upper left", frameon=False, ncol=1,
              handletextpad=0.5, borderpad=0.2, labelspacing=0.35)


def panel_d(axes, hist: pd.DataFrame, res: dict[str, pd.DataFrame]) -> None:
    tools = [tl for tl in ORDER if tl in set(hist["tool"])]
    for ax, tl in zip(axes.ravel(), tools):
        sub = hist[(hist["tool"] == tl)].sort_values("bin_lo")
        mid = (sub["bin_lo"].to_numpy(float) + sub["bin_hi"].to_numpy(float)) / 2
        binmax = sub.groupby("bin_lo", as_index=False)["density"].max()
        bins = (sub.drop_duplicates("bin_lo")[["bin_lo", "bin_hi"]]
                .merge(binmax, on="bin_lo").sort_values("bin_lo"))
        centres_all = ((bins["bin_lo"].to_numpy(float) + bins["bin_hi"].to_numpy(float)) / 2)
        dens_any = bins["density"].to_numpy(float)
        for grp, colour, filled in (("WT", C_WT, True), ("IVT", C_IVT, False)):
            parts = []
            for _, ssub in sub[sub["group"] == grp].groupby("sample"):
                d = ssub.sort_values("bin_lo")
                parts.append(d["density"].to_numpy(float))
            if not parts:
                continue
            dens = np.mean(np.vstack(parts), axis=0)
            centres = centres_all[: dens.size]
            ax.fill_between(centres, 0, dens, step="mid", color=colour,
                            alpha=0.38 if filled else 0.14, linewidth=0)
            ax.step(centres, dens, where="mid", color=colour, linewidth=1.8)
        # zoom x to bins that actually carry density (quantised scores pile at 1.0)
        nz = dens_any > 0
        lo_all = float(bins["bin_lo"].to_numpy(float)[nz].min()) if nz.any() else float(bins["bin_lo"].min())
        hi_all = float(bins["bin_hi"].to_numpy(float)[nz].max()) if nz.any() else float(bins["bin_hi"].max())
        span = hi_all - lo_all
        pad = max(0.02, span * 0.06)
        ax.set_xlim(lo_all - pad, hi_all + pad)
        ax.set_ylim(bottom=0)
        kind = "mod. ratio" if tl == "CHEUI_m5C" else "prob."
        ax.set_xlabel(f"{TOOL_LABEL[tl]}\n{kind}", fontsize=12.5, labelpad=2)
        if span < 0.15:                        # quantised scores: zoomed tick row
            ticks = [round(lo_all - pad, 3), round((lo_all + hi_all) / 2, 3), 1.0]
            ax.set_xticks(ticks)
            ax.set_xticklabels([f"{t:g}" for t in ticks], fontsize=10.5)
        else:
            ax.set_xticks([0, 0.5, 1.0])
            ax.set_xticklabels(["0", "0.5", "1"], fontsize=10.5)
        ax.set_yticks([])
    for ax in axes.ravel()[len(tools):]:
        ax.axis("off")


def panel_e(axes, eco: pd.DataFrame) -> None:
    bins = eco[eco["row_type"] == "prob_bin"]
    for ax, mod in zip(axes, ("m5C", "m6A")):
        for cond, colour, filled in (("WT", C_WT, True), ("IVT", C_IVT, False)):
            d = bins[(bins["mod_type"] == mod) & (bins["condition"] == cond)].sort_values("bin_lo")
            centres = d["bin_lo"].to_numpy(float) + 0.01
            dens = d["density"].to_numpy(float)
            ax.fill_between(centres, 0, dens, step="mid", color=colour,
                            alpha=0.38 if filled else 0.14, linewidth=0)
            ax.step(centres, dens, where="mid", color=colour, linewidth=1.8)
        ax.set_xlim(0, 1)
        ax.set_ylim(bottom=0)
        ax.set_xlabel(f"CHEUI-{mod}\nprobability", fontsize=13, labelpad=2)
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_xticklabels(["0", "0.5", "1"], fontsize=10.5)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stem", type=Path, default=EV / "FigR3_9_nonm6a_fp")
    args = ap.parse_args()
    apply_style()
    data = load()

    fig = plt.figure(figsize=(7.087, 8.9))          # 180 mm wide
    gs = fig.add_gridspec(5, 6, height_ratios=[0.95, 1.0, 0.72, 0.72, 0.85],
                          hspace=0.88, wspace=1.15,
                          left=0.088, right=0.985, top=0.985, bottom=0.055)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0:2])
    ax_c1 = fig.add_subplot(gs[1, 2:4])
    ax_c2 = fig.add_subplot(gs[1, 4:6])
    d_axes = np.array([[fig.add_subplot(gs[2, 0:2]), fig.add_subplot(gs[2, 2:4]),
                        fig.add_subplot(gs[2, 4:6])],
                       [fig.add_subplot(gs[3, 0:2]), fig.add_subplot(gs[3, 2:4]),
                        fig.add_subplot(gs[3, 4:6])]])
    e_axes = [fig.add_subplot(gs[4, 0:3]), fig.add_subplot(gs[4, 3:6])]

    panel_a(ax_a, data["truth"])
    panel_b(ax_b, data["cc"])
    panel_c1(ax_c1, data["rep"])
    panel_c2(ax_c2, data["sum"])
    panel_d(d_axes, data["hist"], data)
    panel_e(e_axes, data["eco"])

    letter(ax_a, "A", dx=0.0, dy=1.16)
    letter(ax_b, "B", dy=1.18)
    letter(ax_c1, "C")
    letter(d_axes[0, 0], "D", dx=-0.02, dy=1.26)
    letter(e_axes[0], "E", dx=-0.02, dy=1.18)
    e_axes[0].set_ylabel("density")
    d_axes[0, 0].set_ylabel("density")
    d_axes[1, 0].set_ylabel("density")
    handles = [Line2D([], [], color=C_WT, linewidth=6, alpha=0.5, label="WT"),
               Line2D([], [], color=C_IVT, linewidth=2.2, label="IVT")]
    d_axes[0, 0].legend(handles=handles, loc="upper right", frameon=False,
                        handletextpad=0.5, borderpad=0.2, labelspacing=0.3)

    for stem_suffix, dpi in (("", 300), ("_600dpi", 600)):
        fig.savefig(f"{args.stem}{stem_suffix}.pdf", bbox_inches="tight", pad_inches=0.02)
        fig.savefig(f"{args.stem}{stem_suffix}.png", bbox_inches="tight",
                    pad_inches=0.02, dpi=dpi)
    plt.close(fig)
    print("wrote", args.stem, "(pdf + png 300/600 dpi)")


if __name__ == "__main__":
    main()
