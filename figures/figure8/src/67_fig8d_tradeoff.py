#!/usr/bin/env python3
"""Figure 8 panel D: PPV against GLORI (2 bp) vs. false positives per 10 kb."""
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from fig8_style import (FAM_COLOR, FAM_LABEL, PIECE, STYLE, TABLES, apply_style,
                        fit_labels, letter, save_piece, title)
import matplotlib.pyplot as plt

FAM_ORDER = ["m6A_DRACH", "m6A", "m5C", "pseU", "inosine"]


def fam_of(name: str) -> str:
    for key in ("m6A_DRACH", "m5C", "pseU", "inosine"):
        if str(name).startswith(key):
            return key
    return "m6A"


def main() -> None:
    apply_style()
    df = pd.read_csv(TABLES / "fig8_tradeoff.tsv", sep="\t")
    df = df[df["platform"] == "RNA004"].copy()
    df["fam"] = [fam_of(f) for f in df["family"]]
    fig = plt.figure()
    ax = fig.add_axes([0.135, 0.330, 0.845, 0.580])
    ax.set_xlim(1e-5, 6.0)
    ax.fill_between([1e-6, 0.12], 0.66, 1.02, color="#E6F2E6", zorder=0)
    ax.plot([1e-6, 0.12], [0.66, 0.66], color="#7FAF7F", lw=1.0, zorder=0)
    ax.plot([0.12, 0.12], [0.66, 1.02], color="#7FAF7F", lw=1.0, zorder=0)
    ax.set_xscale("log")
    for fam in FAM_ORDER:
        sub = df[df["fam"] == fam]
        if sub.empty:
            continue
        # No error bar: every tool was run once on this library (one sequencing
        # unit), so no interval is drawn anywhere in Fig. 8 (user decision).
        ax.scatter(sub["fp_per_10kb"], sub["ppv_glori_w2"], s=62,
                   facecolor=FAM_COLOR[fam], edgecolor="black", linewidth=0.8,
                   zorder=3)
    # two lines: the single-line label is wider than the panel and was clipped
    ax.set_xlabel("False positives per 10 kb\n(unmodified HeLa IVT)",
                  fontsize=STYLE["label"])
    ax.set_ylabel("PPV vs. GLORI (2 bp)", fontsize=STYLE["label"])
    ax.set_ylim(-0.02, 1.0)
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "0.25", "0.50", "0.75", "1.00"])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # only the families this panel actually plots: the RNA004 trade-off table
    # carries no m5C or pseU row (the GLORI reference is m6A-specific), so those
    # two keys were dangling -- a legend entry whose points are not on the axes
    plotted = [f for f in FAM_ORDER if not df[df["fam"] == f].empty]
    handles = [Line2D([], [], marker="o", ls="none", ms=7, mec="k",
                      mfc=FAM_COLOR[f], label=FAM_LABEL[f]) for f in plotted]
    # two columns: three keys in one row ran 2 pt past the piece's right edge
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 0.012), handletextpad=0.4,
               columnspacing=0.9, labelspacing=0.28,
               fontsize=STYLE["legend"])

    # no panel title (letters only)
    letter(fig, "D")
    save_piece(fig, "fig8D_tradeoff", *PIECE["D"], fit=("left",))


if __name__ == "__main__":
    main()
