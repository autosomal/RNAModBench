#!/usr/bin/env python
"""S8 panel drawing -- every panel on its own print-size canvas.

Six panels, no panel titles (the figure legend carries the sample names, the
Figure 8 convention), every key next to the panel it belongs to, and every entry
that is drawn is named somewhere on the page:

A  detection counts as grouped horizontal bars on a log axis: one pair per
   entry, WT in the project's WT blue, IVT in the project's IVT orange; rows
   clustered in three blocks (Dorado m6A models, other m6A tools, other
   modification models).  An entry with no call gets an open left triangle on
   the axis floor instead of a zero-length bar.  Its WT / IVT key sits in the
   free strip above the axes.
B  the eight ORCA channels, same layout as A with its own WT / IVT key.
C  the unmodified RNA004 Curlcake control: false positives per 10^6 candidate
   adenosines as one bar per entry (IVT orange -- it is an unmodified IVT
   library) on a log axis, with the matching specificity on the top axis.
D  matching window (0-50 bp) against PPV (upper tile) and the exact-nucleotide
   fraction (lower tile).  All thirteen entries are drawn and all thirteen are
   named in the key under the tiles; the best DRACH model carries the project's
   "selected" red, the 2-bp working point is a dashed guide.
E  Pearson r against GLORI with its Fisher-z 95 % CI, one lollipop per ratio
   tool (the row with the highest r highlighted).
F  the effect size: the OLS slope of the tool ratio on the GLORI ratio with its
   bootstrap 95 % CI, one lollipop per ratio tool, reference line at 1.0.

Panels of one column share the width of their label column (``left_in``), so
their plotting areas start at the same x and the page reads as a grid.

House rules: Arial only, no grid, no rules inside a panel, no small in-panel
annotation text, no number printed inside a panel, keys in the free strip of
their own panel (never over the data), panel letter in the panel's left margin.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

import s8_style as S
from s8_style import FS, FAM, IVT, WT, GUIDE

matplotlib.use("Agg")

#: the key of panel D -- one entry, one colour, no colour used twice, so any
#: curve of the sweep can be traced back to its name.  The Dorado models take
#: well separated tab20 entries, the best DRACH model the project's "selected"
#: red, the two benchmark tools that the main figures single out their own hue.
D_ENTRY: dict[str, tuple[str, str]] = {
    "Dorado_hac@v5.0.0_m6A@v1": ("hac v5.0.0 m6A", "#1f77b4"),
    "Dorado_hac@v5.0.0_m6A_DRACH@v1": ("hac v5.0.0 m6A DRACH", "#9467bd"),
    "Dorado_hac@v5.1.0_inosine_m6A": ("hac v5.1.0 inosine+m6A", "#bcbd22"),
    "Dorado_hac@v5.1.0_m6A_DRACH@v1": ("hac v5.1.0 m6A DRACH", "#e377c2"),
    "Dorado_sup@v5.0.0_m6A@v1": ("sup v5.0.0 m6A", "#17becf"),
    "Dorado_sup@v5.0.0_m6A_DRACH@v1": ("sup v5.0.0 m6A DRACH", "#b02418"),
    "Dorado_sup@v5.1.0_inosine_m6A": ("sup v5.1.0 inosine+m6A", "#c49c94"),
    "Dorado_sup@v5.1.0_m6A_DRACH@v1": ("sup v5.1.0 m6A DRACH", "#aec7e8"),
    "DRUMMER": ("DRUMMER", "#ff7f0e"),
    "ELIGOS2_diff": ("ELIGOS2_diff", "#ffbb78"),
    "ELIGOS2_solo": ("ELIGOS2_solo", "#2ca02c"),
    "NanoSPA_m6A": ("NanoSPA_m6A", "#8c564b"),
    "m6Anet": ("m6Anet", "#1E888B"),
}
#: key order: the eight Dorado models first (hac, then sup), then the tools
D_ORDER: list[str] = list(D_ENTRY)

#: the best DRACH model carries the project's "selected" colour
#: the working window of the primary analysis (dashed guide in panel D)
WORK_POINT_BP = 2.0

#: bar geometry: a pair per row, +/- this offset in row units, this height
BAR_OFFSET = 0.19
BAR_HEIGHT = 0.30
BAR_EDGE_LW = 0.5

#: hard cap on the number of points any single call may draw (no scatter clouds)
MAX_DOTS_PER_CALL = 120
DOTS_DRAWN = 0


def _dots(ax, xs, ys, **kw):
    """Draw isolated points, refusing clouds (the house rule bans raw scatter)."""
    global DOTS_DRAWN
    n = len(np.atleast_1d(xs))
    if n > MAX_DOTS_PER_CALL:
        raise SystemExit(f"{n} points in one _dots call exceeds "
                         f"{MAX_DOTS_PER_CALL} -- no scatter clouds allowed")
    DOTS_DRAWN += n
    return ax.plot(xs, ys, ls="none", **kw)


def _text_width_in(text: str, fontsize: float) -> float:
    from common.pagelayout import text_width_in
    return text_width_in(text, fontsize)


#: allowance on a measured label width: ``text_width_in`` returns the ink of the
#: string, while matplotlib lays the label out with its side bearings as well
_LABEL_ALLOWANCE = 1.08
_LABEL_PAD_IN = 0.02


def label_width_in(labels: list[str], fontsize: float) -> float:
    """Widest label in inches, incl. the layout allowance; blank rows ignored.

    Public because the renderer measures the three panels of a column and hands
    the common width back as ``left_in``.
    """
    vals = [_text_width_in(t, fontsize) for t in labels if t and t.strip()]
    return max(vals, default=0.4) * _LABEL_ALLOWANCE + _LABEL_PAD_IN



#: readings.  A *measured* zero (the entry made no call on that library -- in the
#: Curlcake panel a zero false-positive result) is now an open circle; an entry
#: that was never measured on that library stays a hollow left triangle, so the
#: two can never be confused.
def _top_key(fig: plt.Figure, handles: list, ncol: int, *, x: float = 0.075,
             y: float = 0.995) -> None:
    """Key in the free strip above the axes, right of the panel letter."""
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(x, y),
               ncol=ncol, frameon=False, fontsize=FS["legend"],
               handlelength=1.3, columnspacing=1.2, handletextpad=0.45,
               borderaxespad=0.0, borderpad=0.0)


# --------------------------------------------------------------------------- #
# A / B / C -- horizontal bar panels
# --------------------------------------------------------------------------- #
def _rows(df: pd.DataFrame, columns: list[str],
          ) -> tuple[list[str], list[int], np.ndarray]:
    """(labels, header indices, values) of a bar panel.

    Every block contributes one bold header row that carries no bar, so the
    blocks are separated by a blank line instead of a rule.
    """
    labels: list[str] = []
    heads: list[int] = []
    values: list[list[float]] = []
    for block, g in df.groupby("block", sort=False, observed=True):
        heads.append(len(labels))
        labels.append(str(block))
        values.append([np.nan] * len(columns))
        for _, r in g.iterrows():
            labels.append(str(r["label"]))
            values.append([float(r[c]) for c in columns])
    return labels, heads, np.array(values, dtype=float)


def bar_labels(df: pd.DataFrame, paired: bool = True) -> list[str]:
    """Row labels a bar panel would print (for the shared label column)."""
    return _rows(df, ["WT", "IVT"] if paired else ["per1e6"])[0]


def _draw_bar_panel(fig: plt.Figure, df: pd.DataFrame, key: str, row_fs: float,
                    left_in: float, *, paired: bool, xlabel: str, floor: float,
                    ticks: list[float], y0: float, height: float,
                    sec_ticks: dict | None = None) -> None:
    """One horizontal-bar panel; ``paired`` draws the WT / IVT bar pair.

    ``floor`` is the x position of the "no call" triangle; the axis starts a
    little left of it and ends 30 % above the largest value, so the log scale
    covers exactly the range the data need.
    """
    columns = ["WT", "IVT"] if paired else ["per1e6"]
    labels, heads, mat = _rows(df, columns)
    n = mat.shape[0]
    y = np.arange(n, dtype=float)
    w_in = fig.get_figwidth()
    left = left_in / w_in
    ax = fig.add_axes([left, y0, 0.975 - left, height])
    S.bare(ax)
    ax.set_xscale("log")
    ax.set_xlim(floor * 0.90, np.nanmax(mat) * 1.35)
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.tick_params(axis="x", labelsize=FS["tick_sm"], pad=3.6)
    ax.set_xlabel(xlabel, fontsize=FS["axis_sm"], labelpad=2.4)
    ax.set_ylim(n - 0.45, -0.55)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=row_fs)
    for k, t in enumerate(ax.get_yticklabels()):
        if k in heads:
            t.set_fontweight("bold")
    ax.tick_params(axis="y", length=0, pad=3.6)

    slots = ((-BAR_OFFSET, 0, WT), (BAR_OFFSET, 1, IVT)) if paired \
        else ((0.0, 0, IVT),)
    for slot, j, colour in slots:
        value = mat[:, j]
        called = np.isfinite(value) & (value > 0)
        if called.any():
            ax.barh(y[called] + slot, value[called], height=BAR_HEIGHT,
                    color=colour, edgecolor="black", linewidth=BAR_EDGE_LW,
                    zorder=3)
        for i in np.where(np.isfinite(value) & (value == 0))[0]:
            # measured zero: open circle on the axis floor, never a zero-length
            # bar (a log axis cannot show 0)
            ax.plot([floor * 1.06], [y[i] + slot], marker="o", ms=4.0,
                    mfc="white", mec=colour, mew=0.9, ls="none", zorder=4)
        
        # table, so its hollow-triangle convention was removed here and from
        # the Figure 8 caption; measured zeros stay open circles above.

    
    # caption, not in the key; the key keeps only the condition patches.
    if paired:
        entries = [Patch(facecolor=WT, edgecolor="black",
                         linewidth=BAR_EDGE_LW, label="WT"),
                   Patch(facecolor=IVT, edgecolor="black",
                         linewidth=BAR_EDGE_LW, label="IVT")]
        _top_key(fig, entries, ncol=len(entries))
    
    #: single-entry key ("unmodified IVT") was redundant and is gone -- the axis
    #: label and the caption name the library
    if sec_ticks:
        sec = ax.twiny()
        sec.set_xscale("log")
        sec.set_xlim(ax.get_xlim())
        S.log_ticks(sec, list(sec_ticks), list(sec_ticks.values()))
        sec.tick_params(axis="x", labelsize=FS["tick_xs"], pad=3.4)
        #: the twin only carries the top FP-rate axis; its inherited y ticks sit
        #: on top of the row labels and the layout gate reads them as frame
        #: marks, so every twin y tick line and label is switched off
        #: (2026-09-25; ``yaxis.set_visible`` does not reach these painters)
        for _tick in sec.yaxis.get_major_ticks():
            _tick.tick1line.set_visible(False)
            _tick.tick2line.set_visible(False)
            _tick.label1.set_visible(False)
            _tick.label2.set_visible(False)
        
        #: tick label (the strip above is free), so it clears it
        sec.set_xlabel("FP rate", fontsize=FS["axis_sm"], labelpad=12.0)
        sec.spines["top"].set_linewidth(S.AXIS_LW)
        for side in ("left", "right", "bottom"):
            sec.spines[side].set_visible(False)
    S.letter(fig, key)


#: x value of the false-positive-rate ticks of panel C: rate -> printed label.

#: duplicated number, so the axis now prints the complementary error rate in
#: scientific notation, where the order of magnitude is immediate (mathtext is
#: already set to Arial in the style module, so the glyphs match the figure).
#: The axis caps at max(rate) * 1.35 ~ 4657; the right tick sits at 3000 so the
#: wider scientific-notation label keeps ~40 pt of clearance to the axis end
#: (at 4000 the label pushed 0.049 in past the canvas).
SPEC_TICKS = {100.0: "$10^{-4}$", 3000.0: "$3\\times10^{-3}$"}


def draw_a(fig: plt.Figure, df: pd.DataFrame, left_in: float) -> None:
    """A -- detection counts per entry, WT | IVT (23 entries)."""
    #: row 1 shrank to 3.60 in on 2026-09-25, so the x label needs a slightly
    #: higher floor than the 0.105 of the 4.08 in canvas
    _draw_bar_panel(fig, df, "A", FS["row_a"], left_in, paired=True,
                    xlabel="Detected sites (log scale)", floor=0.75,
                    ticks=[1, 100, 10000], y0=0.120, height=0.800)


def draw_orca(fig: plt.Figure, df: pd.DataFrame, left_in: float) -> None:
    """B -- the eight ORCA channels, same layout as A."""
    _draw_bar_panel(fig, df, "B", FS["row"], left_in, paired=True,
                    xlabel="Detected sites (log scale)", floor=10.0,
                    ticks=[10, 100, 1000, 10000], y0=0.120, height=0.800)


def draw_curlcake(fig: plt.Figure, df: pd.DataFrame, left_in: float) -> None:
    """C -- false positives per 10^6 candidates on the unmodified control.

    One bar per entry (the control is an unmodified IVT library, so the bars
    carry the IVT colour) with the matching false-positive rate on the top axis.
    """
    
    _draw_bar_panel(fig, df, "C", FS["row"], left_in, paired=False,
                    xlabel="FP per 10$^6$ candidates", floor=100.0,
                    ticks=[100, 1000], y0=0.120, height=0.720,
                    sec_ticks=SPEC_TICKS)


# --------------------------------------------------------------------------- #
# D -- matching window against detection quality
# --------------------------------------------------------------------------- #
def grid_key(fig: plt.Figure, items: list[tuple[str, str, object]], *,
             x_in: float, y_top_in: float, ncol: int = 2,
             fontsize: float | None = None, sw_in: float = 0.185,
             row_h_in: float = 0.152, col_gap_in: float = 0.18,
             centre: bool = False) -> None:
    """A hand-laid key on a fixed grid, so its columns line up exactly.

    ``matplotlib`` legend columns race each other (each column is only as wide
    as its own longest label and the entries are centred in it).  Here every
    swatch of a column starts at the same x and every row sits on the same
    baseline, which is what makes the key read as a table.
    """
    fs = FS["legend"] if fontsize is None else fontsize
    win, hin = fig.get_figwidth(), fig.get_figheight()
    nrow = -(-len(items) // ncol)
    col_w = [sw_in + 0.05 + max(_text_width_in(name, fs)
                                for _t, name, _c in items[c * nrow:(c + 1) * nrow])
             for c in range(ncol)]
    if centre:                     # centre the whole block on the canvas
        total = sum(col_w) + (ncol - 1) * col_gap_in
        x_in = max(x_in, (win - total) / 2.0)
    xs, x = [], x_in
    for w in col_w:
        xs.append(x)
        x += w + col_gap_in
    for i, (_tool, name, colour) in enumerate(items):
        col, row = divmod(i, nrow)
        yf = (y_top_in - (row + 0.5) * row_h_in) / hin
        x0 = xs[col]
        fig.add_artist(Line2D([x0 / win, (x0 + sw_in) / win], [yf, yf],
                              color=colour, lw=1.5, transform=fig.transFigure))
        fig.add_artist(Line2D([(x0 + sw_in / 2.0) / win], [yf], color=colour,
                              marker="o", ms=3.2, mfc=colour, mec="white",
                              mew=0.3, ls="none", transform=fig.transFigure))
        fig.text((x0 + sw_in + 0.05) / win, yf, name, fontsize=fs,
                 va="center", ha="left")


def window_key(loc: pd.DataFrame) -> list[tuple[str, str, object]]:
    """(tool, printed name, colour) of every entry of the window sweep.

    Built from the palette above; a tool that the table draws but the palette
    does not know is a hard error, never a silently grey curve.
    """
    present = {str(t) for t in loc["tool"].unique()}
    missing = sorted(present - set(D_ENTRY))
    if missing:
        raise SystemExit(f"panel D: no key entry for {missing}")
    return [(tool, *D_ENTRY[tool]) for tool in D_ORDER if tool in present]


def _window_axes(fig: plt.Figure, rect: list[float], loc: pd.DataFrame,
                 value: str, ylabel: str) -> plt.Axes:
    """One window-sweep facet: every entry drawn against the matching window."""
    ax = fig.add_axes(rect)
    S.bare(ax)
    ax.set_xlim(-1.5, 52.5)
    ax.set_ylim(0, 100)
    ax.set_xticks([0, 10, 20, 30, 40, 50])
    ax.set_yticks([0, 50, 100])
    ax.tick_params(labelsize=FS["tick_sm"], pad=3.6)
    ax.axvline(WORK_POINT_BP, color=GUIDE, lw=0.9, ls=(0, (3.0, 2.5)),
               zorder=1)
    ax.set_xlabel("Matching window (bp)", fontsize=FS["axis_sm"], labelpad=2.4)
    ax.set_ylabel(ylabel, fontsize=FS["axis_sm"], labelpad=3.0)
    colour_of = {tool: colour for tool, _name, colour in window_key(loc)}
    for tool, g in loc.groupby("tool", sort=False):
        g = g.sort_values("window")
        colour = colour_of[str(tool)]
        ax.plot(g["window"], g[value] * 100.0, color=colour, lw=1.5,
                marker="o", ms=3.2, mfc=colour, mec="white", mew=0.3, zorder=3)
    return ax


def draw_window_ppv(fig: plt.Figure, loc: pd.DataFrame, left_in: float) -> None:
    """D -- PPV against the matching window (left facet of the sweep row)."""
    win = fig.get_figwidth()
    w = (win - left_in - 0.12) / win
    _window_axes(fig, [left_in / win, 0.205, w, 0.680], loc, "hit_rate",
                 "PPV vs. GLORI (%)")
    S.letter(fig, "D")


def draw_window_exact(fig: plt.Figure, loc: pd.DataFrame,
                      left_in: float) -> None:
    """E -- exact-nucleotide fraction (right facet) with the 13-entry key.

    2026-09-25 (user): the sweep is two numbered facets (D, E) instead of one
    full-width row, and the key of the sweep sits on the right of E.
    """
    win = fig.get_figwidth()
    key_w = 3.35
    w = (win - left_in - key_w - 0.55) / win
    _window_axes(fig, [left_in / win, 0.205, w, 0.680], loc,
                 "localization_accuracy", "Exact fraction (%)")
    grid_key(fig, window_key(loc), x_in=left_in + w * win + 0.45,
             y_top_in=1.66, ncol=2, fontsize=FS["tick_sm"], sw_in=0.140,
             col_gap_in=0.100, row_h_in=0.176)
    S.letter(fig, "E")


# --------------------------------------------------------------------------- #
# E / F -- correlation and effect size as lollipops
# --------------------------------------------------------------------------- #
def effect_labels(rows: pd.DataFrame) -> list[str]:
    """Row labels of the two lollipop panels (for the shared label column)."""
    return [str(v) for v in rows["label"]]


def _lollipop_axes(fig: plt.Figure, rect: list[float], rows: pd.DataFrame,
                   value: str, lo: str, hi: str, xlabel: str,
                   show_labels: bool, ref: float | None = None,
                   series2: dict | None = None) -> plt.Axes:
    """One row per tool: a capped CI whisker and a dot.

    ``series2`` (2026-09-26) adds a second estimate on the same row -- used by
    the RNA004 panel G, which reports the calibration slope and Lin's concordance
    correlation coefficient side by side.  The two estimates are dodged within
    the row (filled dot above, open dot below) so neither interval hides the
    other; no in-panel text is added.
    """
    labels = effect_labels(rows)
    n = len(labels)
    y = np.arange(n, dtype=float)
    ax = fig.add_axes(rect)
    S.bare(ax)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(n - 0.45, -0.55)
    ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
    ax.set_yticks(y)
    ax.set_yticklabels(labels if show_labels else [""] * n, fontsize=FS["row"])
    ax.tick_params(axis="x", labelsize=FS["tick_sm"], pad=3.6)
    #: 2026-09-25: pad 5.0 -- at 8 pt the longest row label cleared the invisible
    #: tick mark by only 2.8 pt, under the 3 pt the layout gate asks for
    ax.tick_params(axis="y", length=0, pad=5.0)
    ax.set_xlabel(xlabel, fontsize=FS["axis_sm"], labelpad=2.4)
    if ref is not None:
        ax.axvline(ref, color=GUIDE, lw=1.0, ls=(0, (4.0, 3.0)), zorder=1)
    cap = 0.16
    dodge = 0.17 if series2 else 0.0
    for i_, r in rows.iterrows():
        
        #: (m6Anet in both panels) reads from its own estimate and interval
        colour = S.DOT_COLOR
        yy = y[i_] - dodge
        lo_v, hi_v = max(r[lo], 0.0), r[hi]
        ax.plot([lo_v, hi_v], [yy] * 2, color=colour, lw=1.5, zorder=3)
        
        #: to vanish completely behind the 5 pt dot; brackets make it readable,
        #: and a CI narrower than the dot gets a smaller dot so the bracket
        #: still shows on both sides (the dot size carries no data)
        ax.plot([lo_v, lo_v], [yy - cap, yy + cap], color=colour, lw=1.6,
                zorder=3)
        ax.plot([hi_v, hi_v], [yy - cap, yy + cap], color=colour, lw=1.6,
                zorder=3)
        ms = 5.0 if (hi_v - lo_v) >= 0.045 else 3.0
        _dots(ax, [r[value]], [yy], marker="o", ms=ms, mfc=colour,
              mec="white", mew=0.8, zorder=4)
        if series2:
            yy2 = y[i_] + dodge
            lo2, hi2 = max(r[series2["lo"]], 0.0), r[series2["hi"]]
            ax.plot([lo2, hi2], [yy2] * 2, color=colour, lw=1.2, zorder=3)
            ax.plot([lo2, lo2], [yy2 - cap, yy2 + cap], color=colour, lw=1.3,
                    zorder=3)
            ax.plot([hi2, hi2], [yy2 - cap, yy2 + cap], color=colour, lw=1.3,
                    zorder=3)
            ms2 = 4.6 if (hi2 - lo2) >= 0.045 else 2.8
            _dots(ax, [r[series2["value"]]], [yy2], marker="o", ms=ms2,
                  mfc="white", mec=colour, mew=1.2, zorder=4)
    return ax


def _effect_key(fig: plt.Figure, ax: plt.Axes,
                ref_label: str | None = None,
                series2_label: str | None = None) -> None:
    """Corner key of the two lollipop panels (E/F had none).

    An axes legend, not a figure legend: inside the axes it is audited by
    ``assert_legend_clear`` (no data point may sit under it), which is exactly
    the check the empty lower-right corner of these two panels has to pass.
    """
    
    #: estimate handle (an optional reference-line handle stays available);
    #: 2026-09-26: panel G adds the Lin's CCC handle (open dot)
    handles = [
        Line2D([], [], color=S.DOT_COLOR, marker="o", ms=5.0, mfc=S.DOT_COLOR,
               mec="white", mew=0.8, ls="none", label="estimate, 95 % CI"),
    ]
    if series2_label:
        handles.append(Line2D([], [], color=S.DOT_COLOR, marker="o", ms=4.6,
                              mfc="white", mec=S.DOT_COLOR, mew=1.2, ls="none",
                              label=series2_label))
    if ref_label:
        handles.append(Line2D([], [], color=GUIDE, lw=1.0, ls=(0, (4.0, 3.0)),
                              label=ref_label))
    
    #: **right** of the two facets (a figure legend, so it is not boxed by the
    #: Slope axes and cannot collide with a CI end)
    fig.legend(handles=handles, loc="center right", bbox_to_anchor=(0.968, 0.55),
               ncol=1, frameon=False, fontsize=FS["legend"],
               handlelength=1.3, handletextpad=0.45, columnspacing=1.4,
               labelspacing=0.30, borderaxespad=0.0, borderpad=0.0)


def draw_effect_r(fig: plt.Figure, rows: pd.DataFrame, left_in: float) -> None:
    """F -- Pearson r against GLORI (left facet of the effect-size row).

    The two nine-row lollipop facets share one row grid; the labels are printed
    once, on F.  The dashed "slope = 1" reference stays dropped (2026-09-24):
    every tool sits well below proportional tracking.
    """
    win = fig.get_figwidth()
    width = (win - left_in - 0.12) / win
    _lollipop_axes(fig, [left_in / win, 0.205, width, 0.680], rows,
                   "r", "r_lo", "r_hi", "Pearson r", show_labels=True)
    S.letter(fig, "F")


def draw_effect_slope(fig: plt.Figure, rows: pd.DataFrame,
                      left_in: float) -> None:
    """G -- calibration slope and Lin's CCC of the tool ratio, key on the right.

    2026-09-26 (user): the facet carries two estimates per tool -- the
    least-squares calibration slope (filled dot) and the Lin's concordance
    correlation coefficient (open dot) -- so that association, agreement and
    calibration are all reported for the RNA004 library.
    """
    win = fig.get_figwidth()
    width = 3.65 / win                   # the same facet width as F
    ax_s = _lollipop_axes(fig, [left_in / win, 0.205, width, 0.680], rows,
                          "slope", "slope_lo", "slope_hi", "Slope or Lin's CCC",
                          show_labels=False,
                          series2={"value": "ccc", "lo": "ccc_lo", "hi": "ccc_hi"})
    _effect_key(fig, ax_s, series2_label="Lin's CCC, 95 % CI")
    S.letter(fig, "G")
