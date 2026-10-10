#!/usr/bin/env python3
"""Figure 8 panel C as a standalone piece: PPV dot plot | calls bar.

2026-10-08 (author decision): the piece spans the row's second and third columns; the two plot
areas inside it were rebalanced twice the same day.  The dot plot's left edge is
pinned by the left fit at the measured row-label width (1.379 in), so only its
right edge can move: PPV right edge 2.164 -> 2.380 -> 3.079 in (pass 2: "C
 ", pass 3: " ").

2026-10-10 (author note: "abcd a
"): the piece narrows from 4.51667 to **3.92 in** (fig8_style.ROW1_C_IN) and the
pair inside it is re-cut to **dot 1.30 in : calls 0.99 in** (gap 0.19 in), so the
dot plot no longer runs past the old second-column boundary; the dot's right edge
goes 3.079 -> 2.679 in and the calls bar 3.269-4.497 -> 2.869-3.859 in.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fig8_style import (DOT_COLOR, GUIDE, PIECE, STYLE, TABLES, apply_style,
                        letter, save_piece, title)

#: 2026-10-08 (author decision, three equal columns per row): the row labels stay verbatim and
#: drop to 7.5 pt; at that size the widest one measures 1.379 in, which the left fit
#: turns into the pinned left edge of the dot plot.  Rebalanced twice the same day:
#: the dot plot's right edge 2.164 -> 2.380 -> 3.079 in (1.70 in wide) while the
#: calls bar keeps its own right edge and shrinks to 1.23 in -- the piece itself,
#: the panel letter and the label column do not move.

#: the pair inside is re-cut to dot 1.30 in : calls 0.99 in (gap 0.19 in), i.e.
#: fractions below are x/3.92 of the narrower canvas.
ROW_PT = 7.5


def label_of(tool: str) -> str:
    return str(tool).replace("Dorado_", "").replace("_otherMod", "").replace("@v1", "")


def main() -> None:
    apply_style()
    df = pd.read_csv(TABLES / "fig8_ppv_glori.tsv", sep="\t")
    df = df.sort_values("ppv_glori_w2").reset_index(drop=True)

    fig = plt.figure()
    #: 2026-10-08 (author decision, three equal columns per row; rebalanced twice more the same
    #: day): C spans the row's second and third columns, and the two plot areas
    #: inside it move by their **right/left edges** -- the PPV left edge is pinned
    #: by the left fit at the measured row-label width (1.379 in from the piece
    #: left, whatever x0 says: the fit overrides it), so widening the dot plot means
    #: moving its right edge only.  The x1 values therefore read: PPV right edge
    #: 2.164 (thirds) -> 2.380 (pass 2, real width only 1.00 in) -> **3.079 in**
    #: (1.70 in wide); the calls bar keeps its right edge (4.497 in) and starts



    #: y0 is 0.205, not 0.150: the folded x title needs two lines below the axes.
    #: 2026-10-10: fractions of the 3.92 in canvas -- dot 1.379-2.679 in
    #: (0.35179 + 0.33163), calls 2.869-3.859 in (0.73189 + 0.25255), 0.061 in of
    #: white to the piece's right edge.  The dot's x0 is only a seed: the left fit
    #: overrides it with the measured label width, so the right edge is the knob.
    ax = fig.add_axes([0.35179, 0.205, 0.33163, 0.700])
    ax2 = fig.add_axes([0.73189, 0.205, 0.25255, 0.700])
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
    #: 2026-10-08: the three labels shared ~0.9 in in the thirds pass and print at
    #: 7.5 pt; the rebalances handed the axis 1.70 in, the size is kept as it is
    ax.set_xticklabels(["0%", "50%", "100%"], fontsize=7.5)
    ax.set_xlabel("PPV vs. GLORI\n(2 bp)", fontsize=STYLE["label"])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    # no panel title (letters only); the caption names the sample
    letter(fig, "C")
    save_piece(fig, "fig8C_ppv", *PIECE["C"], fit=("left",))


if __name__ == "__main__":
    main()
