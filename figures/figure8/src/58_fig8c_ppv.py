#!/usr/bin/env python3
"""Figure 8 panel C as a standalone piece: PPV vs. GLORI (2 bp)."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fig8_style import (DOT_COLOR, GUIDE, PIECE, STYLE, TABLES, apply_style,
                        letter, save_piece, title)

#: 2026-10-08 (author decision, three equal columns per row): the row labels stay verbatim and
#: drop to 7.5 pt, so the PPV dot plot keeps ~0.9 in inside its 2.18333 in column.
ROW_PT = 7.5


def label_of(tool: str) -> str:
    return str(tool).replace("Dorado_", "").replace("_otherMod", "").replace("@v1", "")


def main() -> None:
    apply_style()
    df = pd.read_csv(TABLES / "fig8_ppv_glori.tsv", sep="\t")
    df = df.sort_values("ppv_glori_w2").reset_index(drop=True)

    fig = plt.figure()
    #: 2026-10-08 (author decision, three equal columns per row): C spans the row's second and
    #: third columns.  The PPV dot plot fills the second (its labels inside it, the
    #: left fit moves x0), the calls bar fills the third; both land on the column
    #: edges, so the row reads as three equal thirds.
    #: y0 is 0.205, not 0.150: the folded x title needs two lines below the axes.
    ax = fig.add_axes([0.2435, 0.205, 0.2354, 0.700])
    ax2 = fig.add_axes([0.5210, 0.205, 0.4746, 0.700])
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
    ax.set_yticklabels([label_of(t) for t in df["tool"]], fontsize=ROW_PT)
    #: 2026-10-08: 3.6 pt, not 2.5 -- C is inside the layout gate now and its
    #: clearance requirement is 3.0 pt; the y floor moves to -0.75 for the same
    #: reason on the bottom spine.
    ax.tick_params(axis="y", length=0, pad=3.6)
    ax2.barh(y, df["n_calls_total"], height=0.62, color=cols,
             edgecolor="black", linewidth=0.6, zorder=2)
    ax2.set_xscale("log")
    ax2.set_xlim(1, 1.4e5)
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.set_ylim(-0.75, len(df) - 0.4)
    ax2.set_xticks([1, 1000])
    ax2.set_xticklabels(["1", "1k"])
    ax2.set_xlabel("calls", fontsize=STYLE["font"])
    ax2.tick_params(axis="y", length=0)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)

    ax.set_xlim(0, 1.10)          # no value labels any more: bars fill the box
    ax.set_ylim(-0.75, len(df) - 0.4)
    ax.set_xticks([0, 0.5, 1.0])
    #: 2026-10-08: the three labels share ~0.9 in now, so they print at 7.5 pt
    ax.set_xticklabels(["0%", "50%", "100%"], fontsize=7.5)
    ax.set_xlabel("PPV vs. GLORI\n(2 bp)", fontsize=STYLE["label"])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    # no panel title (letters only); the caption names the sample
    letter(fig, "C")
    save_piece(fig, "fig8C_ppv", *PIECE["C"], fit=("left",))


if __name__ == "__main__":
    main()
