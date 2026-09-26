"""Shared publication figure style.

House rules:
* English only, Arial for every text element;
* no gridlines inside the plotting area;
* no small annotation text inside a figure (2026-09-19): explanatory strings,
  value callouts and cell numbers below the tick-label size do not belong in a
  panel -- they go into the panel title, the legend, or the tables
  (``figure1b_replicates`` set the precedent: the n-counts moved from a grey
  line above the axes into the panel title).  The audit list of the remaining
  offenders is ``$RNAMODBENCH_LOCAL/figure_text_audit/``;
* the data area should dominate the canvas -- tight layouts, no wasted margin,
  legends placed inside the axes when they would otherwise buy empty space.
"""

from __future__ import annotations

from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np

STYLE = {
    "font.family": "Arial",
    "font.size": 13,
    "axes.titlesize": 15,
    "axes.labelsize": 14,
    "axes.linewidth": 1.1,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
    "axes.grid": False,
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,   # editable text in the vector PDF
    "ps.fonttype": 42,
}


def apply() -> None:
    rcParams.update(STYLE)


def save(fig, stem: str) -> None:
    """Write ``<stem>.pdf`` and ``<stem>.png`` side by side."""
    fig.savefig(stem + ".pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(stem + ".png", bbox_inches="tight", pad_inches=0.02, dpi=300)
    plt.close(fig)


#: short x-axis labels of the concatenated metagene axis
SEGMENT_LABELS = {"five_prime_flank": "1kb", "five_prime_UTR": "5'UTR",
                  "CDS": "CDS", "three_prime_UTR": "3'UTR",
                  "three_prime_flank": "1kb", "ncRNA_body": "ncRNA"}


def guitar_panel(ax, segs: list[str], title: str, *, flank_color: str = "#4b81b8",
                 ylabel: str = "Density") -> None:
    """Frame a metagene axis like the published Guitar panels."""
    n = len(segs)
    ax.set_title(title, fontweight="bold", fontsize=16)
    ax.set_xticks(np.arange(n) + 0.5)
    ax.set_xticklabels([SEGMENT_LABELS[s] for s in segs], fontsize=13)
    ax.set_xlim(0, n)
    ax.set_ylim(bottom=0)
    ax.tick_params(direction="out", length=4, labelsize=13)
    ax.minorticks_off()
    for i in range(1, n):
        ax.axvline(i, color="k", lw=1.8, ls=":", zorder=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.set_ylabel(ylabel, fontsize=15)
    # transcript schematic below the axis: baseline, dark flank bars, grey body
    tr = ax.get_xaxis_transform()
    ax.plot([0, n], [-0.030, -0.030], color=flank_color, lw=1.6,
            transform=tr, clip_on=False, zorder=6)
    ax.plot([1, n - 1], [-0.008, -0.008], color="0.80", lw=9,
            transform=tr, clip_on=False, zorder=5, solid_capstyle="butt")
    for a, b in ((0, 1), (n - 1, n)):
        ax.plot([a, b], [-0.008, -0.008], color="k", lw=5,
                transform=tr, clip_on=False, zorder=6, solid_capstyle="butt")
    for y in (0.014, -0.030):
        ax.plot([1, n - 1], [y, y], color="k", lw=1.0, transform=tr,
                clip_on=False, zorder=7)
