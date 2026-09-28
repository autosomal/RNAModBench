#!/usr/bin/env python3
"""Figure 8 panel A as a standalone piece: Dorado model counts, WT vs IVT.

Grouped horizontal bars on a log10 axis: for every Dorado model the WT count
(full-strength family colour) and the IVT count (same colour at 45 % alpha)
are drawn as one pair; rows are clustered by model family (m6A DRACH,
non-DRACH m6A, m5C, pseU, inosine + m6A) and sorted by WT count inside the
cluster.  A model with no call in a library gets no bar: an open left-pointing
triangle sits on the axis floor instead (never a fake zero-length bar, never a
number inside the figure).  Site-level counts from ``tables/fig8_counts.tsv``.

Piece size = 3.30 x 3.02 in (final print size; assembled 1:1 by
``68_assemble_a4.py``).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

from fig8_style import (FAM_COLOR, FAM_LABEL, STYLE, TABLES, PIECE, apply_style,
                        letter, save_piece, title)
import matplotlib.pyplot as plt

WT_SAMPLE, IVT_SAMPLE = "HeLa_RNA004_WT", "HeLa_RNA004_IVT"
FAM_ORDER = ["m6A_DRACH", "m6A", "m5C", "pseU", "inosine"]
FLOOR = 0.75                     # x position of "no call" (log axis floor)
ROW_PT = 8.2                     # 16 rows: the "@" descender sets the row pitch
TICK_PAD_PT = 3.6                # > layout-gate clearance (3.0 pt)


def family_of(row: pd.Series) -> str:
    fam = str(row["family"])
    if fam.startswith("m6A_DRACH"):
        return "m6A_DRACH"
    if fam.startswith("m5C"):
        return "m5C"
    if fam.startswith("pseU"):
        return "pseU"
    if fam.startswith("inosine"):
        return "inosine"
    return "m6A"


#: the v5.1.0 run that loads the whole model set is called "..._all" in the
#: tables; its m6A channel is the DRACH-equivalent model of that release (the
#: same reading as CURLAKE_FAMILY / s9_data.ALL_RUN).  No row of this figure
#: carries that name today; the mapping keeps the code name off the page if one
#: ever does.
ALL_RUN = {"all": "m6A DRACH", "all_Psi": "\u03a8", "all_m5C": "m5C"}


def label_of(row: pd.Series) -> str:
    tool = (str(row["tool"]).replace("Dorado_", "").replace("_otherMod", "")
            .replace("@v1", ""))
    for raw, printed in ALL_RUN.items():
        if tool.endswith("_" + raw):
            tool = tool[: -len(raw)] + printed
    if tool.endswith("_inosine_m6A"):
        base, chan = tool[: -len("_inosine_m6A")], "inosine"
    else:
        base, chan = tool, None
    if row["mod_type"] == "m6A" and chan == "inosine":
        return f"{base}#m6A"
    if row["mod_type"] == "inosine":
        return f"{base}#inosine"
    return base


def build_rows(counts: pd.DataFrame) -> pd.DataFrame:
    df = counts[(counts["sample"].isin([WT_SAMPLE, IVT_SAMPLE]))
                & counts["tool"].str.startswith("Dorado")].copy()
    df["family"] = [family_of(r) for _, r in df.iterrows()]
    df["label"] = [label_of(r) for _, r in df.iterrows()]
    wide = df.pivot_table(index=["mod_type", "tool", "family", "label"],
                          columns="sample", values="n_sites", aggfunc="sum").reset_index()
    wide["wt"] = wide.get(WT_SAMPLE, 0).fillna(0).astype(float)
    wide["ivt"] = wide.get(IVT_SAMPLE, 0).fillna(0).astype(float)
    return wide.sort_values("wt", ascending=False).reset_index(drop=True)


def order_rows(piv: pd.DataFrame) -> pd.DataFrame:
    """Cluster by family (fixed family order), then by WT count descending."""
    rank = {f: i for i, f in enumerate(FAM_ORDER)}
    piv = piv.copy()
    piv["frank"] = piv["family"].map(rank).fillna(len(FAM_ORDER))
    piv = piv.sort_values(["frank", "wt"], ascending=[True, False])
    return piv.reset_index(drop=True)


def main() -> None:
    apply_style()
    counts = pd.read_csv(TABLES / "fig8_counts.tsv", sep="\t")
    piv = order_rows(build_rows(counts))

    fig = plt.figure()
    ax = fig.add_axes([0.420, 0.290, 0.565, 0.645])

    y = np.arange(len(piv), dtype=float)
    fam_cols = list(piv["family"].map(FAM_COLOR))
    wt = piv["wt"].to_numpy(float)
    ivt = piv["ivt"].to_numpy(float)
    called_wt, called_ivt = wt > 0, ivt > 0
    zero_wt, zero_ivt = wt == 0, ivt == 0          # measured zero: no call
    miss_wt, miss_ivt = ~np.isfinite(wt), ~np.isfinite(ivt)   # never measured

    ax.barh(y[ called_wt] - 0.19, wt[called_wt], height=0.30, color=np.array(fam_cols)[called_wt],
            edgecolor="black", linewidth=0.55, zorder=3)
    ax.barh(y[called_ivt] + 0.19, ivt[called_ivt], height=0.30, color=np.array(fam_cols)[called_ivt],
            edgecolor="black", linewidth=0.55, alpha=0.45, zorder=3)

    
    # A measured zero (the model made no call in that library) is an open
    # circle; an entry that was never measured stays a hollow triangle.  Neither
    # is a zero-length bar -- a log axis cannot show 0.
    for slot, zero, miss in ((-0.19, zero_wt, miss_wt), (0.19, zero_ivt, miss_ivt)):
        for i in np.where(zero)[0]:
            ax.plot([FLOOR * 1.04], [y[i] + slot], marker="o", ms=4.0,
                    mfc="white", mec=fam_cols[i], mew=0.95, ls="none", zorder=4,
                    clip_on=False)
        for i in np.where(miss)[0]:
            ax.plot([FLOOR * 1.04], [y[i] + slot], marker="<", ms=4.4,
                    mfc="white", mec=fam_cols[i], mew=0.95, ls="none", zorder=4,
                    clip_on=False)

    # light separators between the family clusters
    fam = list(piv["family"])
    for i in range(len(fam) - 1):
        if fam[i] != fam[i + 1]:
            
            #: the three model families is no longer drawn.
            pass

    ax.set_xscale("log")
    ax.set_xlim(FLOOR * 0.90, max(wt.max(), ivt.max()) * 1.35)
    ax.xaxis.set_major_locator(FixedLocator([1, 100, 10000]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_ylim(len(piv) - 0.40, -0.60)
    ax.set_yticks(y)
    ax.set_yticklabels(piv["label"], fontsize=ROW_PT)
    ax.tick_params(axis="y", length=0, pad=TICK_PAD_PT)
    ax.set_xlabel("Detected sites (log scale)", fontsize=STYLE["label"])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    handles = [Patch(facecolor="0.45", edgecolor="k", linewidth=0.55, label="WT"),
               Patch(facecolor="0.45", edgecolor="k", linewidth=0.55, alpha=0.45,
                     label="IVT"),
               Line2D([], [], color="0.45", ls="none", marker="o", ms=4.0,
                      mfc="white", label="0 calls")]
    if bool(np.isnan(piv[["wt", "ivt"]].to_numpy(float)).any()):
        # only when such a row really exists (never today)
        handles.append(Line2D([], [], color="0.45", ls="none", marker="<",
                              ms=4.4, mfc="white", label="not measured"))
    handles += [Patch(facecolor=FAM_COLOR[f], edgecolor="black", linewidth=0.55,
                      label=FAM_LABEL[f]) for f in FAM_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 0.000), labelspacing=0.18, columnspacing=0.9,
               handletextpad=0.35, borderaxespad=0.0, fontsize=8.3)

    # no panel title (author decision: letters only; the caption carries the text)
    letter(fig, "A")

    save_piece(fig, "fig8A_counts", *PIECE["A"], fit=("left",), keep_right=0.985)


if __name__ == "__main__":
    main()
