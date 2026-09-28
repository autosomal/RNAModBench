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
import re
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

import s8_style as S
from s8_style import FS, FAM, IVT, WT, DOT_COLOR, GUIDE

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
def _top_key(fig: plt.Figure, handles: list, ncol: int, *, x: float = 0.5,
             y: float = 0.995) -> None:
    """Key in the free strip above the axes, right of the panel letter.

    2026-09-27 (user): the canvases are 2.5967 in squares now, so the top key
    prints at the small-legend size -- at 10.5 pt the key of C ran past the right
    edge of its cell.
    """

    #: is centred on its own cell instead of hugging the panel letter.
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(x, y),
               ncol=ncol, frameon=False, fontsize=FS["key_sm"],
               handlelength=1.3, columnspacing=1.2, handletextpad=0.45,
               borderaxespad=0.0, borderpad=0.0)


# --------------------------------------------------------------------------- #
# A / B / C -- horizontal bar panels
# --------------------------------------------------------------------------- #
def letter_near(fig: plt.Figure, key: str, left_in: float,
                labels: list[str] | None, fs: float,
                box_h_in: float, pad_in: float = 0.34) -> None:
    """Panel letter immediately left of this panel's own topmost label.

 2026-09-28 (user): " " -- the letter used to sit
    at the left edge of the cell, which for a panel with short row labels (or with
    none at all) left it a whole label column away from the ink.  It is now placed
    0.30 in left of the top row's own label, so every letter reads as part of its
    panel; the layout gate checks that it clears that label.
    """
    win = fig.get_figwidth()
    if not labels:
        #: a panel with no row labels (the sweep cells): the tick labels and the
        #: rotated y title are its leftmost ink, measured at 0.86 in
        x_in = max(left_in - pad_in - 0.52, 0.020)
    else:

        #: letter is 0.37 in tall, so it shares its band with the first two to
        #: three label rows; the width that matters is the widest of *those*, not
        #: the top one.  F/G put a 1.0 in name in their second row (the letter was
        #: sitting on it) and H/I have short two-row blocks (the letter floated a
        #: row too far left).  The band is now measured from the row pitch.
        pitch = max(box_h_in / max(len(labels), 1), 0.02)
        rows = max(1, min(len(labels), int(-(-0.30 // pitch))))
        width = label_width_in(list(labels[:rows]), fs)
        x_in = max(left_in - width - pad_in, 0.020)
    S.letter(fig, key, x=x_in / win)


def _rows(df: pd.DataFrame, columns: list[str],
          ) -> tuple[list[str], list[int], np.ndarray]:
    """(labels, header indices, values) of a bar panel.

 2026-09-28 (user): "Dorado m6A modelABC ORCA (one tool)" -- with
    the family bands gone the header rows only repeated a grouping the reader
    already gets from the caption and from the model names, and an empty row
    ("ORCA (one tool)") looked like a stray label.  The headers are no longer
    inserted: every row is a model or channel, and the families stay together
    because the blocks keep their order.
    """
    labels: list[str] = []
    heads: list[int] = []
    values: list[list[float]] = []
    for _block, g in df.groupby("block", sort=False, observed=True):
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
                    sec_ticks: dict | None = None,
                    square: bool = True) -> None:
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


    #: span the same width as the eight other panels (canvas minus the shared label
    #: column and the 0.06 in right margin) instead of stopping at 0.900 of the
    #: cell and shrinking to the square where the box was taller than the square.
    ax = fig.add_axes([left, y0,
                       (w_in - 0.060 - left_in) / w_in,
                       height])
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

    #: group name stays a row label but in grey italics, and a light grey band
    #: spans exactly the rows of that group, so the grouping is carried by the
    #: band instead of by bold type.
    for k, t in enumerate(ax.get_yticklabels()):
        if k in heads:
            t.set_fontweight("normal")
            t.set_fontstyle("italic")
            t.set_color("#6f6f6f")

    #: every printed row is a model or a channel.
    _ = heads
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

        #: tick label (the strip above is free), so it clears it; 2026-09-27: the
        #: box now fills the row-1 band, so the strip above it is only 0.33 in and
        #: the title keeps a tight 2 pt pad (its ticks stay clear at 3.4 pt)
        sec.set_xlabel("FP rate", fontsize=FS["axis_sm"], labelpad=2.0)
        sec.spines["top"].set_linewidth(S.AXIS_LW)
        for side in ("left", "right", "bottom"):
            sec.spines[side].set_visible(False)
    letter_near(fig, key, left_in, labels, row_fs,
                fig.get_figheight() - S.BOX_BOTTOM_IN - S.BOX_TOP_PAD_IN)


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

                    #: (3.000 in of the 3.600 in cell), so the row reads as a row
                    ticks=[1, 100, 10000], y0=0.1218, height=0.7706, square=False)


def draw_orca(fig: plt.Figure, df: pd.DataFrame, left_in: float) -> None:
    """B -- the eight ORCA channels, same layout as A."""
    _draw_bar_panel(fig, df, "B", FS["row"], left_in, paired=True,
                    xlabel="Detected sites (log scale)", floor=10.0,

                    ticks=[10, 100, 1000, 10000], y0=0.1218, height=0.7706)


def draw_curlcake(fig: plt.Figure, df: pd.DataFrame, left_in: float) -> None:
    """C -- false positives per 10^6 candidates on the unmodified control.

    One bar per entry (the control is an unmodified IVT library, so the bars
    carry the IVT colour) with the matching false-positive rate on the top axis.
    """


    #: the row share their top and bottom lines
    _draw_bar_panel(fig, df, "C", FS["row"], left_in, paired=False,
                    xlabel="FP per 10$^6$ candidates", floor=100.0,
                    ticks=[100, 1000], y0=0.1218, height=0.7706,
                    sec_ticks=SPEC_TICKS)


# --------------------------------------------------------------------------- #
# D -- matching window against detection quality
# --------------------------------------------------------------------------- #
def grid_key(fig: plt.Figure, items: list[tuple[str, str, object]], *,
             x_in: float, y_top_in: float, ncol: int = 2,
             fontsize: float | None = None, sw_in: float = 0.185,
             row_h_in: float = 0.152, col_gap_in: float = 0.18,
             centre: bool = False, centre_in: float | None = None) -> None:
    """A hand-laid key on a fixed grid, so its columns line up exactly.

    ``matplotlib`` legend columns race each other (each column is only as wide
    as its own longest label and the entries are centred in it).  Here every
    swatch of a column starts at the same x and every row sits on the same
    baseline, which is what makes the key read as a table.
    """
    fs = FS["legend"] if fontsize is None else fontsize
    win, hin = fig.get_figwidth(), fig.get_figheight()

    #: and a half that holds fewer than four entries used to leave an empty column
    #: behind -- ``max()`` then had nothing to measure and the build died.  The
    #: column count is clamped to the entry count and an empty column is worth
    #: zero width, so a half with one or two entries still lays out.
    ncol = max(1, min(int(ncol), len(items)))
    nrow = -(-len(items) // ncol)
    col_w = [sw_in + 0.05 + max((_text_width_in(name, fs) * _LABEL_ALLOWANCE
                                 for _t, name, _c in items[c * nrow:(c + 1) * nrow]),
                                default=0.0)
             for c in range(ncol)]
    if centre_in is not None:      # centre the block on a box, not on the cell
        total = sum(col_w) + (ncol - 1) * col_gap_in
        x_in = centre_in - total / 2.0
        #: a key wider than its box would run off the canvas, so it is kept inside

        #: the key is wider than the canvas (hi < lo), which is how the sweep key left
        #: the cell; the block is now pinned to the left margin instead.
        lo, hi = 0.012, win - 0.012 - total
        x_in = max(lo, min(x_in, hi)) if hi >= lo else lo
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


    #: x title keeps a tight labelpad and leaves the key strip its clearance
    ax.set_xlabel("Matching window (bp)", fontsize=FS["axis_sm"], labelpad=0.8)
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

    _window_axes(fig, S.cell_rect(fig, left_in, height_in=S.BOX_H_LOW_IN), loc,
                 "hit_rate", "PPV vs. GLORI (%)")

    #: no row labels, so the letter sat a whole label column (1.55 in) away from
    #: its box; it now hugs the box, clear of the rotated y title and the tick
    #: labels by 0.36 in (the gate checks both).

    #: no row labels, so the letter sat a whole label column (1.86 in) away from
    #: its box; it now sits in the free band just left of the box.
    S.letter(fig, "D", x=min(left_in - 0.86, 2.40) / fig.get_figwidth())


def draw_window_exact(fig: plt.Figure, loc: pd.DataFrame,
                      left_in: float) -> None:
    """E -- exact-nucleotide fraction (right facet of the sweep).

    2026-09-27 (user): the box fills its cell like D's; the thirteen-entry key of
    the sweep prints under the pair (E draws its second half).
    """
    _window_axes(fig, S.cell_rect(fig, left_in, height_in=S.BOX_H_LOW_IN), loc,
                 "localization_accuracy", "Exact fraction (%)")

    #: no row labels, so the letter sat a whole label column (1.55 in) away from
    #: its box; it now hugs the box, clear of the rotated y title and the tick
    #: labels by 0.36 in (the gate checks both).

    #: no row labels, so the letter sat a whole label column (1.86 in) away from
    #: its box; it now sits in the free band just left of the box.
    S.letter(fig, "E", x=min(left_in - 0.86, 2.40) / fig.get_figwidth())


# --------------------------------------------------------------------------- #
# H / I -- rank stability of the tool ordering (2026-09-27)
# --------------------------------------------------------------------------- #
#: marker per species: the shape carries the species, so no extra hue enters the
#: page (the sweep already spends thirteen colours)
STAB_MARKER = {"Arabidopsis": "o", "mouse": "s", "HeLa": "^"}

#: the library key (HeLa), the printed legend names the species (Human).
SPECIES_DISPLAY = {"Arabidopsis": "Arabidopsis", "mouse": "mouse", "HeLa": "Human"}


def species_key(fig: plt.Figure, printed: list[str], *, centre_x: float) -> None:
    """One row naming the three species markers, **under the panel that draws them**.

    2026-09-27 (user): a key belongs under its own panel -- the three markers are
    drawn in both stability cells, so each of them prints this line under its own
    box (``centre_x`` is the centre of that box, in figure fractions), never in a
    neighbour's free space and never above the axes.
    """
    handles = [Line2D([], [], color=DOT_COLOR, marker=STAB_MARKER[p], ms=4.6,
                      mfc=DOT_COLOR, mec="white", mew=0.6, ls="none",
                      label=SPECIES_DISPLAY.get(p, p))
               for p in printed]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(centre_x, S.KEY_BOTTOM_IN / fig.get_figheight()),
               ncol=len(printed), frameon=False, fontsize=FS["key_sm"],
               handlelength=1.0, columnspacing=1.1, handletextpad=0.35,
               borderaxespad=0.0, borderpad=0.0)


def _stability_axes(fig: plt.Figure, rect: list[float], rows: pd.DataFrame,
                    xlabel: str) -> plt.Axes:
    """One row per level, one marker per species, against the primary ranking.

    The dashed guide marks rho = 1.0, i.e. the level reproduces the primary
    ordering exactly; the three species of a level share the row.
    """
    levels = list(dict.fromkeys(str(v) for v in rows["level"]))
    idx = {lv: i for i, lv in enumerate(levels)}
    n = len(levels)
    ax = fig.add_axes(rect)
    S.bare(ax)
    ax.set_xlim(0.40, 1.02)
    ax.set_ylim(n - 0.45, -0.55)
    ax.set_xticks([0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels(["0.4", "0.6", "0.8", "1"])
    ax.set_yticks(np.arange(n))
    #: 2026-09-27: the level names print at the 7 pt floor -- at 7.5 pt the
    #: longest one ran 0.02 in past the left edge of the 2.5967 in square

    #: unified row size instead of the 7 pt floor.
    ax.set_yticklabels(levels, fontsize=FS["level"])
    ax.tick_params(axis="x", labelsize=FS["tick_sm"], pad=3.6)
    ax.tick_params(axis="y", length=0, pad=5.0)

    #: sits tight against its ticks
    ax.set_xlabel(xlabel, fontsize=FS["axis_sm"], labelpad=0.8)

    for _, r in rows.iterrows():
        _dots(ax, [float(r["rho"])], [idx[str(r["level"])]],
              marker=STAB_MARKER[str(r["printed"])], ms=4.6, mfc=DOT_COLOR,
              mec="white", mew=0.6, zorder=4)
    return ax


def draw_stability_coverage(fig: plt.Figure, rows: pd.DataFrame,
                            left_in: float) -> None:
    """H -- ranking stability against transcript coverage (a square cell).

    2026-09-27 (user): the cells are squares, so the thirteen-entry key of the
    window sweep left this cell; the species key it needs prints **under its own
    box** (the same line I also carries under I), not in any other cell.
    """
    rect = S.cell_rect(fig, left_in, height_in=S.BOX_H_LOW_IN)

    #: edge, so their box is inset a little (F keeps the full width).

    #: inset is gone, so every box spans the full cell width like F and A.
    _stability_axes(fig, rect, rows, "Spearman \u03c1 vs. ranking")
    species_key(fig, [str(v) for v in dict.fromkeys(rows["printed"])],
                centre_x=rect[0] + rect[2] / 2.0)
    _levels = list(dict.fromkeys(str(v) for v in rows["level"]))
    letter_near(fig, "H", left_in, _levels, FS["level"], S.BOX_H_LOW_IN)


def draw_stability_ratio(fig: plt.Figure, rows: pd.DataFrame,
                         left_in: float) -> None:
    """I -- ranking stability against the reference's own modification ratio.

    The reference-composition axis is where the ordering actually moves (the
    manuscript quotes rho = 0.48-0.90 with the top-ranked tool changing), so the
    cell is drawn on the same axes as H; its species key prints under its own box.
    """
    rect = S.cell_rect(fig, left_in, height_in=S.BOX_H_LOW_IN)

    #: edge, so their box is inset a little (F keeps the full width).

    #: inset is gone, so every box spans the full cell width like F and A.
    _stability_axes(fig, rect, rows, "Spearman \u03c1 vs. ranking")
    species_key(fig, [str(v) for v in dict.fromkeys(rows["printed"])],
                centre_x=rect[0] + rect[2] / 2.0)
    _levels = list(dict.fromkeys(str(v) for v in rows["level"]))
    letter_near(fig, "I", left_in, _levels, FS["level"], S.BOX_H_LOW_IN)


def draw_sweep_key(fig: plt.Figure, loc: pd.DataFrame, *, half: str,
                   left_in: float, row_h_in: float = 0.105,
                   rows: int = 3) -> None:
    """The thirteen-entry key of the window sweep, **under the pair that draws it**.

    2026-09-27 (user): a key prints under its own panel, and the sweep is the pair
    D | E -- but thirteen entries do not fit in one cell's strip (their labels
    need ~5 in of width and three 7 pt rows).  Each facet's canvas therefore
    prints one half of the key in the strip under its own box, centred; the two
    halves read as one band under the row.  ``half`` is ``"left"`` (D: the first
    six entries, two columns) or ``"right"`` (E: the remaining seven, three).
    """
    items = window_key(loc)

    #: shrink from 0.72 in to 0.42 in (the S8-like packing the user asked for).

    #: shave a row off the block, but seven entries at 7 pt need 5.06 in and a cell
    #: is 3.79 in wide, so the block left the canvas (the layout gate caught it).

    #: seven-entry half -- both come out three rows tall, and 3 x 1.24 in is the
    #: only arrangement that fits a 3.79 in cell (four columns need 5.06 in and
    #: three columns for the six-entry half need 3.9 in; the gate caught both).
    part, ncol = (items[:6], 2) if half == "left" else (items[6:], 3)

    #: the canvas edge, so the block starts 0.10 in up instead of at KEY_BOTTOM_IN.
    y_top_in = 0.024 + rows * 0.098
    box_centre = left_in + (fig.get_figwidth() - left_in - 0.060) / 2.0

    #: "sup v5.0.0 m6A DRACH") and the unified 8.5 pt needs more air between the
    #: columns, so the gap grows; the strip is 0.72 in, the three rows fit.

    #: the unified 8 pt): at 8 pt its three columns do not fit a 3.79 in cell -- the
    #: layout gate measured 4.59 in for the right half.
    grid_key(fig, part, x_in=0.100, y_top_in=y_top_in, ncol=ncol,
             fontsize=7.0, sw_in=0.075, row_h_in=0.098,
             col_gap_in=0.020, centre_in=box_centre)


# --------------------------------------------------------------------------- #
# E / F -- correlation and effect size as lollipops
# --------------------------------------------------------------------------- #
def effect_labels(rows: pd.DataFrame) -> list[str]:
    """Row labels of the two lollipop panels -- the **full** model names.

 2026-09-28 (user): "FG" -- F and G print "hac v5.1.0 m6A DRACH"
    again, not "hac m6A DRACH".  The shared label column holds them at the unified
    8 pt scale (the layout gate verifies the fit), and G borrows the column from F.
    """
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

    ax.set_xlabel(xlabel, fontsize=FS["axis_sm"], labelpad=0.8)

    #: used to be drawn here (`if ref is not None: ax.axvline(ref, ... dashed)`) is gone;
    #: the parameter stays so the callers that pass a reference keep working.
    _ = ref
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


def _effect_key(fig: plt.Figure, *, centre_x: float,
                series2_label: str | None = None) -> None:
    """Key of the two lollipop facets, in one row **under the box of its panel**.

    2026-09-27 (user): a key prints under the panel it belongs to.  The two
    handles describe marks that F and G draw alike and the pair shares one row,
    so the key sits under the right-hand facet of that row -- the same place the
    thirteen-entry key of the sweep takes in its own row.
    """

    #: estimate handle; 2026-09-26: the facet adds the Lin's CCC handle (open dot)
    handles = [
        Line2D([], [], color=DOT_COLOR, marker="o", ms=5.0, mfc=DOT_COLOR,
               mec="white", mew=0.8, ls="none", label="estimate, 95 % CI"),
    ]
    if series2_label:
        handles.append(Line2D([], [], color=DOT_COLOR, marker="o", ms=4.6,
                              mfc="white", mec=DOT_COLOR, mew=1.2, ls="none",
                              label=series2_label))

    #: the key of this one panel prints at the 7 pt floor
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(centre_x, S.KEY_BOTTOM_IN / fig.get_figheight()),
               ncol=len(handles), frameon=False, fontsize=FS["key"],
               handlelength=1.2, handletextpad=0.40, columnspacing=1.0,
               borderaxespad=0.0, borderpad=0.0)


def draw_effect_r(fig: plt.Figure, rows: pd.DataFrame, left_in: float) -> None:
    """F -- Pearson r against GLORI (left facet of the effect-size row).

    The two nine-row lollipop facets share one row grid; the labels are printed
    once, on F.  The dashed "slope = 1" reference stays dropped (2026-09-24):
    every tool sits well below proportional tracking.
    """

    _lollipop_axes(fig, S.cell_rect(fig, left_in, height_in=S.BOX_H_LOW_IN),
                   rows, "r", "r_lo", "r_hi", "Pearson r", show_labels=True)
    letter_near(fig, "F", left_in, effect_labels(rows), FS["row"],
                S.BOX_H_LOW_IN)


def draw_effect_slope(fig: plt.Figure, rows: pd.DataFrame,
                      left_in: float) -> None:
    """G -- calibration slope and Lin's CCC of the tool ratio.

    2026-09-26 (user): the facet carries two estimates per tool -- the
    least-squares calibration slope (filled dot) and the Lin's concordance
    correlation coefficient (open dot) -- so that association, agreement and
    calibration are all reported for the RNA004 library.  2026-09-27 (user): the
    facet is exactly as wide as F's (one grid column) and the key prints under
    this facet's own box.
    """


    #: Lin's CCC" title and the key together ran past the cell edge (0.28 in measured), so G's
    #: box is inset by that much (F keeps the full width).
    rect = S.cell_rect(fig, left_in, height_in=S.BOX_H_LOW_IN)

    #: inset is gone, so every box spans the full cell width like F and A.

    #: The 11:01 state of this file had them (G was self-contained), and dropping
    #: them to make the full model names fit was the wrong trade: the names are the
    #: only way to read the facet.  G keeps F's box left edge, so the two facets
    #: still line up, and the gate below fails if the labels leave the cell.
    _lollipop_axes(fig, rect, rows, "slope", "slope_lo", "slope_hi",
                   "Slope or Lin's CCC", show_labels=True,
                   series2={"value": "ccc", "lo": "ccc_lo", "hi": "ccc_hi"})
    #: the key is ~2.2 in wide against a 1.87 in box, so it is centred on the
    #: canvas instead of on the box -- that overflow is what forced the inset
    _effect_key(fig, centre_x=0.5,
                series2_label="Lin's CCC, 95 % CI")
    letter_near(fig, "G", left_in, effect_labels(rows), FS["row"],
                S.BOX_H_LOW_IN)
