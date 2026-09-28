#!/usr/bin/env python3
"""Figure 8 panel B: strict FPR on the two unmodified controls (two sub-axes).

B1 (top)  unmodified Curlcake IVT - false positives per 10 kb against the
          Dorado percent-modified threshold (5/10/20/50 %), one series per
          model family x call mode (sup: solid line, filled circle; hac:
          dashed line, open triangle).
B2 (below) unmodified HeLa IVT at the delivered >= 90 % cutoff - one open
          diamond per delivered model, family mean as a horizontal bar.

No number is written inside the figure: the false-positive counts live in the
legend and in ``tables/fig8_fpr_curlcake_scan.tsv`` /
``tables/fig8_fpr_hela_ivt.tsv`` (a series with zero calls at a threshold
is an open marker on the axis floor, never a fallen line).  Piece size =
3.30 x 3.02 in.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from fig8_style import (B_BOX, B_FLOOR, B_YLIM, B_YSHARE, B_YTICKS, B_YLABELS,
                        FAM_COLOR, STYLE, TABLES, PIECE, apply_style, letter,
                        save_piece, title)

SLOT_OF = {"m6A_DRACH": "m6A_DRACH", "pseU_m6A": "m6A", "inosine_m6A": "inosine"}
SLOTS = ["m6A_DRACH", "m6A", "inosine"]
FAM_SHORT = {"m6A_DRACH": "m6A DRACH", "m6A": "non-DRACH m6A",
             "inosine": "inosine+m6A"}
B2_XLAB = {"m6A_DRACH": "DRACH", "m6A": "non-DRACH", "inosine": "inosine"}
HELA_FAM = {"m6A_DRACH": ["m6A_DRACH"], "m6A": ["m6A"],
            "inosine": ["inosine_m6A"]}
MODE_STYLE = {"sup": dict(ls="-", marker="o", mfc="full"),
              "hac": dict(ls=(0, (3.2, 1.6)), marker="^", mfc="white")}
B1_XLAB = ["5", "10", "20", "50"]
THR = [5, 10, 20, 50]
TICK_PAD_PT = 3.6               # > layout-gate clearance (3.0 pt)
TITLE_TEXT = "False positives per 10 kb"    # shared rotated y title (both sub-axes)
TITLE_X = 0.020                 # its anchor, figure fraction (panel's left edge)
TITLE_GAP_PT = 1.5              # required clearance between title band and tick labels


def fit_b1_left(fig, ax1, title_artist, pad_pt: float = TICK_PAD_PT,
                tick_len_pt: float = 2.4, gap_pt: float = TITLE_GAP_PT,
                k: float = 1.15):
    """Move B1 right so its y tick labels clear the shared rotated title.

    2026-09-23 : B1's
    left margin (``B_BOX["B1"]``, 24.9 pt) was exactly the width of its tick
    labels (24.8 pt including the pad), so the shared title -- a ``fig.text``
    that ``fig8_style._fit`` cannot see, because that helper only measures an
    axis' own ``ylabel`` -- was printed on top of the topmost label ("1,000",
    gap -11.2 pt).  The margin is therefore measured here instead of hardcoded:

        half the rotated title band  +  clearance  +  tick length + tick pad
                                     +  widest y tick label

    Two corrections are applied to the renderer measurements, both checked
    against the exported piece with ``pdftotext -bbox``: the rotated band is
    widened by ``k = 1.15`` (the renderer reports 9.8 pt where the PDF prints
    11.2 pt) and the margin is expressed in the **final** piece width
    (``PIECE["B"]``) because ``save_piece`` sets the canvas size only at the end.

    Returns ``(new_x0, title_band_pt, widest_tick_pt)``.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    px2pt = 72.0 / fig.dpi
    band_pt = title_artist.get_window_extent(renderer).width * px2pt * k
    tick_pt = max((t.get_window_extent(renderer).width * px2pt
                   for t in ax1.get_yticklabels() if t.get_text()), default=0.0)
    panel_pt = PIECE["B"][0] * 72.0
    need_pt = TITLE_X * panel_pt + band_pt / 2 + gap_pt + tick_len_pt + pad_pt + tick_pt
    x0 = min(need_pt / panel_pt, 0.75)
    pos = ax1.get_position()
    ax1.set_position([x0, pos.y0, pos.width, pos.height])
    fig.canvas.draw()
    return x0, band_pt, tick_pt


def series_for(cc: pd.DataFrame, slot: str, mode: str) -> dict:
    """FP per 10 kb at THR for one family slot and call mode."""
    sub = cc[(cc["mod_type"] == "m6A") & (cc["mode"] == mode)
             & cc["threshold_pct"].notna()]
    fams = [k for k, v in SLOT_OF.items() if v == slot]
    sub = sub[sub["family"].isin(fams)]
    vals, nfp = [], []
    for t in THR:
        hit = sub[sub["threshold_pct"] == t].sort_values("fp_per_10kb")
        v = hit["fp_per_10kb"].astype(float)
        vals.append(float(v.mean()) if len(v) and float(v.mean()) > 0 else np.nan)
        nfp.append([int(x) for x in hit["n_fp"]])
    return {"y": vals, "n_fp": nfp}


def hela_points(hl: pd.DataFrame, slot: str) -> list:
    """FP per 10 kb of every HeLa IVT model belonging to one family slot."""
    ivt = hl[(hl["sample"] == "HeLa_RNA004_IVT") & (hl["mod_type"] == "m6A")
             & hl["tool"].str.startswith("Dorado")]
    ivt = ivt[ivt["family"].isin(HELA_FAM[slot])]
    return [(float(r["fp_per_10kb"]), int(r["n_fp"])) for _, r in ivt.iterrows()]


def styles(ax, yticks, ylabels) -> None:
    ax.set_yscale("log")
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=STYLE["tick"])
    ax.tick_params(axis="y", length=2.4, width=0.9, pad=TICK_PAD_PT)
    ax.tick_params(axis="x", length=2.4, width=0.9, pad=TICK_PAD_PT)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def main() -> None:
    apply_style()
    cc = pd.read_csv(TABLES / "fig8_fpr_curlcake_scan.tsv", sep="\t")
    hl = pd.read_csv(TABLES / "fig8_fpr_hela_ivt.tsv", sep="\t")

    fig = plt.figure()

    ## ---- B1: Curlcake threshold scan ----------------------------------------- #
    ax1 = fig.add_axes(B_BOX["B1"])
    x = np.arange(len(THR))
    for slot in SLOTS:
        for mode in ("sup", "hac"):
            s = series_for(cc, slot, mode)
            st = MODE_STYLE[mode]
            mfc = FAM_COLOR[slot] if st["mfc"] == "full" else "white"
            ax1.plot(x, s["y"], ls=st["ls"], lw=0.9, color=FAM_COLOR[slot],
                     marker=st["marker"], ms=4.6, mfc=mfc, mec="black",
                     mew=0.85, zorder=3)
            for i, v in enumerate(s["y"]):
                if np.isnan(v):
                    ax1.plot([x[i]], [B_FLOOR["B1"]], marker="o", ms=4.0,
                             mfc="white", mec=FAM_COLOR[slot], mew=1.05,
                             ls="none", zorder=4)
    ax1.set_ylim(*B_YLIM["B1"])
    styles(ax1, B_YTICKS["B1"], B_YLABELS["B1"])
    ax1.set_xlim(-0.55, len(THR) - 0.45)
    ax1.set_xticks(x)
    ax1.set_xticklabels(B1_XLAB, fontsize=STYLE["tick"])
    ax1.set_xlabel("threshold (%)", fontsize=STYLE["label"], labelpad=2.5)
    ax1.set_title("Curlcake IVT", fontsize=STYLE["label"],
                  fontweight="normal", loc="left", pad=4)

    # the shared rotated y title is created *before* B1 is positioned: its band
    # width is what fixes B1's left margin (see fit_b1_left)
    ylab = fig.text(TITLE_X, (B_YSHARE[0] + B_YSHARE[1]) / 2, TITLE_TEXT,
                    rotation=90, ha="center", va="center",
                    fontsize=STYLE["label"])
    x0, band_pt, tick_pt = fit_b1_left(fig, ax1, ylab)
    print(f"[fig8B] shared title band {band_pt:.2f} pt | widest y tick label {tick_pt:.2f} pt | "
          f"B1 left edge {B_BOX['B1'][0]:.3f} -> {x0:.3f}", flush=True)

    ## ---- B2: HeLa IVT at the delivered cutoff -------------------------------- #
    ax2 = fig.add_axes(B_BOX["B2"])
    slots_x = np.arange(len(SLOTS), dtype=float)
    for i, slot in enumerate(SLOTS):
        pts = hela_points(hl, slot)
        vals = [p[0] for p in pts]
        jit = np.linspace(-0.17, 0.17, len(vals))
        ax2.scatter(slots_x[i] + jit, vals, s=30, facecolor="white",
                    edgecolor=FAM_COLOR[slot], linewidth=1.15, marker="D",
                    zorder=4)
        ax2.plot([slots_x[i] - 0.24, slots_x[i] + 0.24], [float(np.mean(vals))] * 2,
                 color=FAM_COLOR[slot], lw=1.3, zorder=3)
    ax2.set_ylim(*B_YLIM["B2"])
    styles(ax2, B_YTICKS["B2"], B_YLABELS["B2"])
    ax2.set_xlim(-0.55, len(SLOTS) - 0.45)
    ax2.set_xticks(slots_x)
    # three family names in 1.47 in: upright they collided, so they are slanted
    ax2.set_xticklabels([B2_XLAB[s] for s in SLOTS], fontsize=STYLE["tick"],
                        rotation=30, ha="right", rotation_mode="anchor")
    ax2.set_title("Human IVT (\u2265 90 %)", fontsize=STYLE["label"],
                  fontweight="normal", loc="left", pad=4)

    # (the shared y title itself is created above, before B1 is positioned; its
    #  anchor spans the shared extent of the two sub-axes)

    handles = [Patch(facecolor=FAM_COLOR[SLOTS[0]], edgecolor="black",
                     linewidth=0.55, label=FAM_SHORT[SLOTS[0]]),
               Line2D([], [], color="0.1", ls=MODE_STYLE["sup"]["ls"],
                      marker="o", ms=4.6, label="sup"),
               Patch(facecolor=FAM_COLOR[SLOTS[1]], edgecolor="black",
                     linewidth=0.55, label=FAM_SHORT[SLOTS[1]]),
               Line2D([], [], color="0.1", ls=MODE_STYLE["hac"]["ls"],
                      marker="^", ms=4.6, mfc="white", label="hac"),
               Patch(facecolor=FAM_COLOR[SLOTS[2]], edgecolor="black",
                     linewidth=0.55, label=FAM_SHORT[SLOTS[2]]),
               Line2D([], [], color="0.1", ls="none", marker="o", ms=4.0,
                      mfc="white", label="0 calls")]   # the same "measured
                      # zero" wording as panels A/C and Figure S8: the open
                      # circle is a zero, never a missing measurement
    # matplotlib fills a legend column by column, so the interleaved order above
    # puts the three model families on the first row and the keys on the second.
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.49, 0.002), handletextpad=0.35, columnspacing=0.55,
               labelspacing=0.22, borderaxespad=0.0, fontsize=8.3)

    # no panel title (letters only): the two block labels inside the piece already
    # name the controls
    letter(fig, "B")
    save_piece(fig, "fig8B_fpr", *PIECE["B"], fit=())


if __name__ == "__main__":
    main()
