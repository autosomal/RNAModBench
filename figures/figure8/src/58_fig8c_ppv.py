#!/usr/bin/env python3
"""Figure 8 panel C as a standalone piece: PPV vs. GLORI (2 bp)."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fig8_style import (DOT_COLOR, GUIDE, PIECE, STYLE, TABLES, apply_style,
                        letter, save_piece, title)


def label_of(tool: str) -> str:
    return str(tool).replace("Dorado_", "").replace("_otherMod", "").replace("@v1", "")


def main() -> None:
    apply_style()
    df = pd.read_csv(TABLES / "fig8_ppv_glori.tsv", sep="\t")
    df = df.sort_values("ppv_glori_w2").reset_index(drop=True)

    fig = plt.figure()
    ax = fig.add_axes([0.385, 0.150, 0.430, 0.755])
    ax2 = fig.add_axes([0.855, 0.150, 0.115, 0.755])
    y = np.arange(len(df))
    fam_col = {"m6A_DRACH": "#E83C1E", "m6A": "#F5B264", "m5C": "#2E86AB",
               "pseU": "#A23B72", "inosine": "#F18F01",
               "inosine_m6A": "#F18F01", "m5C@v1": "#2E86AB",
               "pseU@v1": "#A23B72"}
    cols = [fam_col.get(str(f), DOT_COLOR) for f in df["family"]]
    # No error bar: every tool was run once on the RNA004 HeLa WT library (single
    # sequencing unit), so a replicate interval cannot be drawn.  The site-level
    # bootstrap interval is reported in tables/fig8_ppv_glori.tsv, never in the
    # figure (author decision, 2026-09-21).
    sizes = 18 + 62 * (df["n_calls_total"] / df["n_calls_total"].max())
    ax.scatter(df["ppv_glori_w2"], y, s=sizes, c=cols, edgecolor="black",
               linewidth=0.8, zorder=3)
    # No value labels: the bar length carries the PPV, the numbers live in the
    # caption and in tables/fig8_ppv_glori.tsv (house rule: no in-figure values).
    ax.set_yticks(y)
    ax.set_yticklabels([label_of(t) for t in df["tool"]], fontsize=STYLE["row"])
    ax.tick_params(axis="y", length=0, pad=2.5)
    ax2.barh(y, df["n_calls_total"], height=0.62, color=cols,
             edgecolor="black", linewidth=0.6, zorder=2)
    ax2.set_xscale("log")
    ax2.set_xlim(1, 1.4e5)
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.set_ylim(-0.6, len(df) - 0.4)
    ax2.set_xticks([1, 1000])
    ax2.set_xticklabels(["1", "1k"])
    ax2.set_xlabel("calls", fontsize=STYLE["font"])
    ax2.tick_params(axis="y", length=0)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)

    ax.set_xlim(0, 1.10)          # no value labels any more: bars fill the box
    ax.set_ylim(-0.6, len(df) - 0.4)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["0%", "50%", "100%"])
    ax.set_xlabel("PPV vs. GLORI (2 bp)", fontsize=STYLE["label"])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    # no panel title (letters only); the caption names the sample
    letter(fig, "C")
    save_piece(fig, "fig8C_ppv", *PIECE["C"], fit=("left",))


if __name__ == "__main__":
    main()
