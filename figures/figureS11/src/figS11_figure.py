#!/usr/bin/env python
"""NA-4b figure -- PPV and F1 of six representative tool variants as the calling
threshold is swept, on the manuscript (w = 2 bp, cov >= 10 universe, drawn
units of Fig. 5A).

Layout (since the 2026-10-10 transpose): **two metric rows x three species columns** -- row 1 is
PPV, row 2 is F1, and the columns are Arabidopsis, mouse and HeLa from left to right -- in the
house style of Figure 3 / Figure S9 / Figure S5: **every panel carries its own axis titles** (the
y title repeats its row's metric, the x title repeats below every panel) **and its own x tick
numbers**, so both rows read on their own even though the x limits are shared; one bold species
header per **column** is centred over that column's two panels; the bold row letters (A = the PPV
row, B = the F1 row) sit in the left margin; and **each row carries its own framed six-entry tool
key** below it.  A framed six-entry tool key sits below everything (author decisions 2026-10-10: letters
name a block, the key must not sit on the data and carries the house frame, one species name
serves both panels of its column, two-decimal tick values that fill the panel's own limit; the
"table" variant with shared column headers was rejected in favour of the other figures'
style).  ``x`` is the
mean fraction of the units' own calls that survive the threshold
(100 % = the full distributed callset = the published operating point); the documented
operating point of each tool is drawn as an open circle, so a reader can see how much is
gained or lost by moving the threshold away from the published one.  Absolute thresholds
are in ``tables/figS11_sweep_summary.tsv`` (they are tool-specific: probabilities for
m6Anet, mod-ratios for DENA / Nanom6A, GMM p-values for Nanocompore, adjusted p-values
for ELIGOS2).

2026-10-10 (R1-5 shipping round): Nanocompore joins the sweep as the second
differential-testing representative, the panels are lettered, and each panel is
also written as a single-panel file under ``figures/panels/`` so the delivery can
carry the same per-panel set as every other figure.

2026-10-10 (aggregation + layout round, on the author's request): every curve is now
**averaged over the independent sequencing units at each threshold** (x = mean retained
fraction, y = mean metric), so a plotted point really is the mean over the units the
caption names.  The first version grouped on the raw per-unit ``retained_fraction`` and
therefore drew an interleaving of single units' points -- measured in
``$RNAMODBENCH_LOCAL/submission/_build/s11_letters_legend_20261010``: 93.8 % of the drawn symbols
came from one unit while the caption promised means over units.  The per-tile letters A-F
became the two block letters A/B in the same round, the shared legend moved below the
panels, and ``_assert_style`` now guards the unit count behind every point.

2026-10-10 (margins + free-scale round, on the author's request: " 
F1 vs GLORI (2 bp) 0.20 "): two defects of the shipped page are
fixed.  *Margins*: the page used to leave 0 pt, measured from the word boxes of the saved PDF
-- left -0.2 pt, top -0.5 pt (the top title and the left y titles were shaved), right 1.4 pt,
bottom 7.4 pt -- because the house ``figstyle.save`` hard-codes ``pad_inches=0.02`` (1.4 pt);
this page now saves itself through ``_save_padded`` with ``MARGIN_PT = 10`` clear points on
every side, and ``FIGSIZE`` shrank by the same 20 pt so the product width stays 178 mm.  *Y
scales*: the limit used to be one value per metric column (``_vmax``: f1 -> 0.3686) while the
ticks were the fixed 0.00-0.20 set, so 46 % of the F1 axis carried no tick and the
Arabidopsis / Human F1 panels filled only 55 % / 64 % of their height; each panel now scales to
its own curves (``panel_vmax`` x ``Y_LIM_HEADROOM``) with the locator's nice two-decimal ticks,
and ``_draw_panel`` asserts the panel's curves reach at least 92 % of the drawn height.

2026-10-10 (x-labels round, on the author's request: "Tool calls retained (% 
"): the page used to leave the upper two rows without any x information -- ``sharex`` hid
their 0/25/50/75/100 numbers and the x title was printed only under the bottom row.  Now every
one of the six panels prints ``X_LABEL`` and re-enables its own bottom tick labels
(``tick_params(labelbottom=True)``; the shared limits stay linked), the two middle gaps carry the
previous row's numbers and title above that row's species header, and the canvas grew by 55 pt
so the panels themselves keep their size.

House rules: English only, Arial, no gridlines, no in-panel annotation text,
178 mm page width, minimum font 9 pt (asserted at draw time), pdf + png.

Outputs
-------
figures/FigureS11_rev.{pdf,png}         the six-panel page
figures/panels/FigureS11_panel{A1..A3,B1..B3}.{pdf,png}   single panels

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/figures/figureS11/src/figS11_figure.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FormatStrFormatter, MaxNLocator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "analysis/coverage_rank_stability/analysis"))
sys.path.insert(0, str(HERE))

from r21_common import PANELS, setup_logger                              # noqa: E402
from common.figstyle import apply as apply_style                         # noqa: E402

OUT = HERE.parent
TAB = OUT / "tables"
FIG = OUT / "figures"
PER_UNIT = TAB / "figS11_sweep_per_unit.tsv"
METRIC_COLS = [("precision", "PPV vs GLORI (2 bp)"), ("f1", "F1 vs GLORI (2 bp)")]
#: 2026-10-10 (R1-5 shipping round): Nanocompore joins as the second differential-testing
#: representative, so the page carries six curves and the panels are lettered A-F
TOOL_ORDER = ["m6Anet", "DENA", "Nanom6A", "Nanocompore",
              "ELIGOS2_solo", "ELIGOS2_diff"]
TOOL_COLOR = {"m6Anet": "#C44E52", "DENA": "#4C72B0", "Nanom6A": "#55A868",
              "Nanocompore": "#937860", "ELIGOS2_solo": "#DD8452",
              "ELIGOS2_diff": "#8172B3"}
TOOL_MARKER = {"m6Anet": "o", "DENA": "s", "Nanom6A": "^",
               "Nanocompore": "P", "ELIGOS2_solo": "D", "ELIGOS2_diff": "v"}

#: read top to bottom (Arabidopsis, mouse, HeLa): A = the PPV column, B = the F1 column.
#: The single-panel files keep a letter per tile, named after their block plus the row.
BLOCK_LETTERS = ("A", "B")
MIN_FONT = 9.0   # the page claims 9 pt and now enforces it

#: clear printed points on **every** side.  The house ``figstyle.save`` hard-codes
#: ``pad_inches=0.02`` (~1.4 pt) and is shared by 60+ scripts of the frozen ``$RNAMODBENCH_LOCAL/code`` tree,
#: so its default stays untouched; this page saves itself through ``_save_padded`` instead.
MARGIN_PT = 10.0
#: 2026-10-10 (labels round): 5.6 -> 6.0 in, because the species row headers need vertical
#: room between two rows of panels -- at 5.6 in they shared a baseline with the upper row's
#: bottom tick labels (measured: 3 px clearance, i.e. touching).
#: 2026-10-10 (margins round): the tight box grows by 2 x MARGIN_PT = 20 pt once the padding
#: above is applied, so the canvas shrinks by the same 20 pt to keep the **product** size --
#: 178 mm wide, 1:1 in the album, which places this page at ``\textwidth`` unscaled.  The
#: height went 5.76 -> 5.96 in after the key moved into its own reserved band (LEGEND_BAND):
#: that band is taken out of the panels, so the canvas has to grow to keep the page's aspect
#: ratio within 1 % of the shipped 1.0980 (measured 500.72 x 441.05 pt = 1.1353 without it).

#: upper rows now carry their own x numbers and x title as well, which costs about 2 x 27 pt of
#: panel space -- the canvas grows 5.96 -> 6.72 in so the **panels keep their size** and only
#: the two middle gaps get taller (5.96 + 55 pt/72 would be 6.72; the measured page height
#: lands near 510 pt, vs 455.306 before, i.e. +55 pt as the author accepted).  Measured at 6.72:
#: the panels came out 2.9 pt shorter than the shipped page (the distance between their 0.00 and
#: 1.00 y labels, 82.93 vs 85.81 pt), so the canvas adds 8.7 pt = 0.12 in to restore them.

#: metric rows instead of two metric columns -- so it gets **shorter** (two panel rows, one
#: middle gap, two keys) while the product width stays 178 mm.  The canvas is re-calibrated by
#: measurement; ``BETWEEN_ROWS_HSPACE`` carries the middle gap (x numbers + x title + row 0's
#: key), ``LETTER_BAND`` the left band of the bold row letters, ``LEGEND_BAND`` the last key.
FIGSIZE = (8.75, 5.07)
BETWEEN_ROWS_HSPACE = 0.15    # middle gap: x numbers + x title + the first row's key
#: 2026-10-10 (transpose round): the left margin of this page is fully taken by the y titles and
#: their five numbers (measured: they span x 8.3-46.3 pt), so the row letters need their own
#: reserved column -- hence a band wide enough for letters + y title + numbers, with the letters
#: pinned at ``LETTER_X`` (figure fraction) so they cannot creep into the labels.
LETTER_BAND = 0.045           # left band: bold row letters A/B, then the y titles and numbers
LETTER_X = 0.012              # where the letters are pinned, in figure fraction


def panel_labels() -> list[str]:
    """Per-tile labels in reading order, aligned with the page's block letters."""
    return [f"{BLOCK_LETTERS[j]}{i + 1}" for i in range(len(PANELS))
            for j in range(len(METRIC_COLS))]


#: 2026-10-10 (labels round, panel-labels round, x-labels round, transpose round): the page
#: follows the house style of Figure 3 / Figure S9 / Figure S5 -- **every panel carries its own
#: axis titles and its own tick numbers**.  Since the transpose the grid is row = metric, column
#: = species: one **column header** per species is centred over that column's two panels, the
#: bold row letters A/B sit in the left margin, and each row carries its own key.  The metric
#: names are *not* printed as extra column headers: the author rejected that "table" look and asked
#: for the style of the other figures (their own y titles already state the metric).
#:

#: each panel scales its **own** y axis (``panel_vmax`` x ``Y_LIM_HEADROOM``) and the ticks are
#: the locator's nice two-decimal values *inside* that limit.  The shipped page shared one limit
#: per metric column -- f1 took the global 0.3511 -> 0.3686 while the fixed ticks stopped at
#: 0.20, so 46 % of the axis height carried no tick, and the Arabidopsis / Human F1 panels
#: filled only 55 % / 64 % of their height.
Y_LIM_HEADROOM = 1.05         # a panel's limit = its own curve maximum x this
Y_TICK_NBINS = 5              # locator target for the number of tick intervals per panel
#: 2026-10-10 (margins round): the shared key sits at **fixed figure coordinates** inside a band
#: the layout engine reserves under the panels, instead of ``loc="outside lower center"``.
#: Measured on the first save of this round: the outside locator let the framed key end up
#: ~6 pt below the tight box that ``savefig(bbox_inches="tight")`` had computed (bottom margin
#: 4.08 pt instead of 10), because it is re-positioned while the save itself draws.  A fixed
#: anchor has nothing to re-position.  ``LEGEND_BAND_PANEL`` is taller because a single-panel
#: file wraps the same six entries into two rows of three.
LEGEND_BAND = 0.060           # figure fraction under the panels, six entries in one row
LEGEND_BAND_PANEL = 0.230     # ... and for a single-panel file, six entries in two rows
X_LABEL = "Tool calls retained (%)"
#: 2026-10-10 (transpose round): three species columns make each panel ~95 pt wide, and the
#: automatic x locator then thins the numbers out -- measured on the first 2x3 render, only
#: ``0`` and ``100`` survived.  The page promises these five, so they are pinned; the gate then
#: asserts their labels keep a real gap.
X_TICKS = (0, 25, 50, 75, 100)

#: numbers inherited the shared figstyle's ``font.size = 13``, i.e. ~13.5 pt as printed -- about
#: 1.5x what the house figures print (Figure 3 ~8.7 pt, S8 ~8.9, S10 ~7.3).  That oversized
#: number set is what made five y numbers collide inside a 62 pt panel: measured 15 word-vs-word
#: collisions on the delivered page.  They now use the house size.
TICK_FONTSIZE = 9.0
#: matplotlib's tick-label boxes are the tight glyph boxes while the PDF's nominal boxes are the
#: font's em boxes -- about 1.5 pt wider on each side.  Geometry gates therefore demand this much
#: clearance instead of mere non-overlap (the delivered page passed "no overlap" with 15 real
#: collisions), and the verifier re-checks the saved PDF word by word.
BOX_CLEARANCE_PT = 3.0
AXIS_LABEL_FONTSIZE = 9.5
HEADER_FONTSIZE = {"row": 13.0, "letter": 11.0}   # "row" = the species **column** headers here
HEADER_GAP = 0.012            # figure fraction between a column header and the panels it names


def species_label(i: int) -> str:
    """``Arabidopsis  (n = 3)`` -- the wording of manuscript Figure 5, per drawn condition."""
    species, _group, units = PANELS[i]
    return f"{species}  ({'n = 1 study' if len(units) == 1 else f'n = {len(units)}'})"


def _texts(fig) -> list:
    out = []
    for ax in fig.axes:
        out += [t for t in ax.texts if t.get_text()]
        out += [ax.title, ax.xaxis.label, ax.yaxis.label]
        out += list(ax.get_xticklabels()) + list(ax.get_yticklabels())
        if ax.get_legend() is not None:
            out += list(ax.get_legend().get_texts())
    #: 2026-10-10: the page carries a figure-level key, and it must obey the 9 pt floor too
    for lg in fig.legends:
        out += list(lg.get_texts())
    out += [t for t in fig.texts if t.get_text()]
    return out


def _overlaps(a, b) -> bool:
    """True when two display-space boxes intersect (a shared pixel is already a defect)."""
    return bool(a.x0 < b.x1 and b.x0 < a.x1 and a.y0 < b.y1 and b.y0 < a.y1)


def _clear(a, b, gap: float = 0.0) -> bool:
    """True when two boxes keep at least ``gap`` points of daylight between them."""
    return bool(a.x0 >= b.x1 + gap or b.x0 >= a.x1 + gap
                or a.y0 >= b.y1 + gap or b.y0 >= a.y1 + gap)


def _assert_style(fig, *, letters: list | None = None, below: list | None = None,
                  cols: list | None = None) -> float:
    """House-style gate; returns the smallest font size.

    ``letters`` are the bold row letters drawn as figure text -- exactly one per **metric row**
    (A = PPV, B = F1), each **outside every axes** and clear of every y title (2026-10-10: the
    the author's complaint was a letter drawn inside the plot area, on top of the data).  ``below``
    are figure-level artists -- the per-row keys -- each of which must sit strictly under **its
    own row's** axes, never on a panel, and must carry the house frame.  Passing ``cols`` turns
 on the page checks of the transposed layout (2026-10-10, author note "2*3 legend"): one
    species header per **column**, centred over that column's two panels; every panel naming its
    own y title (its row's metric) and x title, with its x numbers visible; no panel title.
    """
    for ax in fig.axes:
        for gl in list(ax.get_xgridlines()) + list(ax.get_ygridlines()):
            if gl.get_visible():
                raise AssertionError(f"gridline visible on axis {ax.get_ylabel()!r}")
        if len([t for t in ax.texts if t.get_text()]):
            raise AssertionError(f"in-panel annotation text on {ax.get_ylabel()!r}")
        #: 2026-10-10 (transpose round, tightened in the legibility round): five pinned x numbers
        #: must fit a panel's width and five y numbers its height -- with ``BOX_CLEARANCE_PT`` of
        #: daylight, not merely without overlapping.  The delivered page passed the old
        #: "no overlap" version with 15 real collisions, because matplotlib's boxes are the tight
        #: glyph boxes while the PDF's are the em boxes.
        for kind, labels, key in (("x", ax.get_xticklabels(), lambda b: b.x0),
                                  ("y", ax.get_yticklabels(), lambda b: b.y0)):
            boxes = sorted((t.get_window_extent() for t in labels if t.get_text()), key=key)
            for a, b in zip(boxes, boxes[1:]):
                if not _clear(a, b, BOX_CLEARANCE_PT):
                    raise AssertionError(
                        f"{kind} tick numbers have less than {BOX_CLEARANCE_PT:g} pt of daylight "
                        f"on a panel of {ax.get_window_extent().width:.1f} x "
                        f"{ax.get_window_extent().height:.1f} pt -- shrink the font or grow the "
                        f"panel")

        #: panel sticks out into the neighbouring label bands -- measured +15.7 pt on the
        #: delivered page.  Its span must stay inside the span of that panel's y numbers.
        if ax.get_ylabel():
            nb = [t.get_window_extent() for t in ax.get_yticklabels() if t.get_text()]
            if nb:
                tspan = ax.yaxis.label.get_window_extent()
                lo = min(b.y0 for b in nb)
                hi = max(b.y1 for b in nb)
                if tspan.y0 < lo - BOX_CLEARANCE_PT or tspan.y1 > hi + BOX_CLEARANCE_PT:
                    raise AssertionError(
                        f"the y title sticks out of its panel: title y {tspan.y0:.1f}.."
                        f"{tspan.y1:.1f} vs that panel's y numbers {lo:.1f}..{hi:.1f} pt")
    if letters is not None:
        labels = [t.get_text() for t in letters]
        if labels != list(BLOCK_LETTERS):
            raise AssertionError(f"expected block letters {list(BLOCK_LETTERS)}, "
                                 f"found {labels}")
        for j, t in enumerate(letters):
            lbox = t.get_window_extent()
            #: 2026-10-10 (letters round): the letter belongs to the **top-left corner** of its
            #: own row's block -- above that row's panels and left of every axes on the page
            row = fig.axes[j * len(PANELS):(j + 1) * len(PANELS)]
            for ax in row:
                if lbox.y0 < ax.get_window_extent().y1:
                    raise AssertionError(f"row letter {t.get_text()!r} is not above its own row "
                                         f"of panels")
            for ax in fig.axes:
                if lbox.x1 > ax.get_window_extent().x0:
                    raise AssertionError(f"row letter {t.get_text()!r} is not in the left margin "
                                         f"(left of every panel)")
                if _overlaps(lbox, ax.get_window_extent()):
                    raise AssertionError(f"block letter {t.get_text()!r} sits on an axes")
                if ax.get_ylabel() and _overlaps(lbox, ax.yaxis.label.get_window_extent()):
                    raise AssertionError(f"row letter {t.get_text()!r} touches a y title")
                #: 2026-10-10 (transpose round): on this narrow page the letters had to be moved
                #: into a reserved left column -- measured, they used to sit on top of the
                #: y numbers.  This gate makes that collision impossible to reintroduce.
                for tkl in ax.get_yticklabels():
                    if tkl.get_text() and _overlaps(lbox, tkl.get_window_extent()):
                        raise AssertionError(f"row letter {t.get_text()!r} touches the y numbers "
                                             f"{tkl.get_text()!r}")
                if ax.title.get_text() and _overlaps(lbox, ax.title.get_window_extent()):
                    raise AssertionError(f"block letter {t.get_text()!r} touches a title")
    if cols is not None:
        nrow, ncol = len(METRIC_COLS), len(PANELS)
        if len(cols) != ncol:
            raise AssertionError(f"expected {ncol} column headers, found {len(cols)}")
        tol = 2.0 / fig.bbox.width          # the 2 pt the plan allows, in figure fractions
        for i, t in enumerate(cols):
            left = fig.axes[i].get_position().x0                 # row 0 of column i
            right = fig.axes[(nrow - 1) * ncol + i].get_position().x1
            if abs(float(t.get_position()[0]) - (left + right) / 2.0) > tol:
                raise AssertionError(f"column header {t.get_text()!r} is not centred over its "
                                     f"two panels")

        #: wide species headers ran into the neighbouring column's rotated y title and into the
        #: top y numbers of their own column.  A header must keep daylight from every y title,
        #: every y number, every panel and its neighbour headers.
        for i, t in enumerate(cols):
            cbox = t.get_window_extent()
            for ax in fig.axes:
                if _overlaps(cbox, ax.get_window_extent()):
                    raise AssertionError(f"column header {t.get_text()!r} sits on a panel")
                if ax.get_ylabel() and _overlaps(cbox, ax.yaxis.label.get_window_extent()):
                    raise AssertionError(f"column header {t.get_text()!r} touches the y title "
                                         f"of another panel")
                for tkl in ax.get_yticklabels():
                    if tkl.get_text() and _overlaps(cbox, tkl.get_window_extent()):
                        raise AssertionError(f"column header {t.get_text()!r} touches the y "
                                             f"numbers {tkl.get_text()!r}")
            for k, other in enumerate(cols):
                if k != i and not _clear(cbox, other.get_window_extent(), BOX_CLEARANCE_PT):
                    raise AssertionError(f"column headers {t.get_text()!r} and "
                                         f"{other.get_text()!r} are closer than "
                                         f"{BOX_CLEARANCE_PT:g} pt")
        #: 2026-10-10 (panel-labels round, x-labels round, transpose round): the house style of
        #: Figure 3 / S8 / S10 -- every panel names **both** of its axes with its own numbers --
        #: and since the transpose the page is row = metric, column = species
        for j in range(nrow):
            for i in range(ncol):
                ax = fig.axes[j * ncol + i]
                if ax.get_ylabel() != METRIC_COLS[j][1]:
                    raise AssertionError(f"panel row {j} col {i}: y title "
                                         f"{ax.get_ylabel()!r} != {METRIC_COLS[j][1]!r}")
                if ax.get_xlabel() != X_LABEL:
                    raise AssertionError(f"panel row {j} col {i}: x title "
                                         f"{ax.get_xlabel()!r} != {X_LABEL!r} (every panel "
                                         f"carries it since 2026-10-10)")
                shown = [t for t in ax.xaxis.get_ticklabels() if t.get_text()]
                if not shown:
                    raise AssertionError(f"panel row {j} col {i}: no x tick numbers")
                if not all(t.get_visible() for t in shown):
                    raise AssertionError(f"panel row {j} col {i}: x tick numbers hidden -- "
                                         f"sharex hides every row but the last unless "
                                         f"tick_params(labelbottom=True)")
                if ax.get_title():
                    raise AssertionError(f"panel row {j} col {i}: unexpected panel title "
                                         f"{ax.get_title()!r}")
    for k, artist in enumerate(below or []):
        bbox = artist.get_window_extent()
        row = fig.axes[k * len(PANELS):(k + 1) * len(PANELS)]   # the axes this key belongs to
        for ax in row:
            if _overlaps(bbox, ax.get_window_extent()) or bbox.y1 > ax.get_window_extent().y0:
                raise AssertionError(f"key {k + 1} is not below the panels of its own row")
        for ax in fig.axes:                 # ... and clear of every other axes as well
            if _overlaps(bbox, ax.get_window_extent()):
                raise AssertionError(f"key {k + 1} overlaps a panel")
            #: 2026-10-10 (legibility round): measured on the single-panel file A1, the key was
            #: sitting on the x title and the x numbers -- "below the axes" is not enough, it has
            #: to clear the labels those axes draw underneath themselves
            if ax.get_xlabel() and not _clear(bbox, ax.xaxis.label.get_window_extent(),
                                              BOX_CLEARANCE_PT):
                raise AssertionError(f"key {k + 1} is closer than {BOX_CLEARANCE_PT:g} pt to the "
                                     f"x title of a panel")
            for tkl in ax.get_xticklabels():
                if tkl.get_text() and not _clear(bbox, tkl.get_window_extent(),
                                                 BOX_CLEARANCE_PT):
                    raise AssertionError(f"key {k + 1} is closer than {BOX_CLEARANCE_PT:g} pt to "
                                         f"the x numbers {tkl.get_text()!r}")
        if hasattr(artist, "get_frame_on") and not artist.get_frame_on():
            raise AssertionError("a shared key must carry the house frame (2026-09-30)")
        for other in (below or [])[:k]:     # the keys must not collide with each other
            if _overlaps(bbox, other.get_window_extent()):
                raise AssertionError(f"keys {k} and {k + 1} overlap")
    if letters is not None and cols is not None:
        for t in letters:                    # a row letter must clear the column headers too
            lbox = t.get_window_extent()
            for c in cols:
                if _overlaps(lbox, c.get_window_extent()):
                    raise AssertionError(f"row letter {t.get_text()!r} touches a column header")
    sizes = [float(t.get_fontsize()) for t in _texts(fig)
             if t.get_text() and t.get_fontsize()]
    if not sizes:
        raise AssertionError("no text found -- figstyle not applied?")
    if min(sizes) < MIN_FONT:
        raise AssertionError(f"smallest font {min(sizes):.2f} pt < {MIN_FONT} pt")
    #: 2026-10-10 (legibility round): the author explicitly accepted a "≥1.3x width" page so the
    #: three columns can breathe; the album prints it at \textwidth, which scales the drawing
    #: back down by ~0.77.  Bounded on both sides so a typo cannot silently move the print scale.
    width_pt = 72.0 * float(fig.get_size_inches()[0])
    if cols is not None and not 600.0 <= width_pt <= 690.0:
        raise AssertionError(f"canvas is {width_pt:.1f} pt wide, outside the agreed 1.30x band "
                             f"(600-690 pt before the 2 x {MARGIN_PT:g} pt page margin)")
    if cols is None and width_pt > 5.0 * 72.0:
        raise AssertionError(f"a single-panel file is {width_pt:.1f} pt wide -- too wide for a "
                             f"tile")
    return min(sizes)


def curve(per_unit: pd.DataFrame, species: str, metric: str) -> dict:
    """{tool: frame} -- one point per calling threshold, **averaged over the units**.

    ``x`` is the mean retained fraction of the units of that condition, ``y`` the mean
    metric over those units, so every plotted point is the mean over the independent
    sequencing units the caption names (three for Arabidopsis and Human, one for the single
    mouse study).  ``figS11_sweep.py`` writes ``retained_fraction`` **per unit**, so grouping
    on it -- what this function did until 2026-10-10 -- turned a curve into an interleaving
    of single units' points: measured in ``_build/s11_letters_legend_20261010``, 93.8 % of
    the drawn symbols came from one unit.  ``n_units``/``n_rows`` stay in the frame so
    ``_draw_panel`` can insist each point averages one row per unit of the condition.
    """
    sub = per_unit[per_unit["species"] == species]
    out: dict[str, pd.DataFrame] = {}
    for tool in TOOL_ORDER:
        t = sub[sub["tool"] == tool]
        if t.empty:
            continue
        g = (t.groupby("threshold", as_index=False)
             .agg(retained_fraction=("retained_fraction", "mean"),
                  value=(metric, "mean"),
                  is_default=("is_default", "max"),
                  n_units=("sample", "nunique"),
                  n_rows=("sample", "size"))
             .sort_values("retained_fraction"))
        out[tool] = g
    return out


def panel_vmax(per_unit: pd.DataFrame, species: str, metric: str) -> float:
    """Highest value drawn by **this panel's** six curves (2026-10-10: free y scales).

    The shipped page shared one limit per metric column, which is why the F1 panels topping out
    at 0.20 / 0.24 sat in a 0.3686-tall box whose top tick was 0.20 -- the blank band the author
    pointed at.  Every panel now fits its own curves; the 0.05 floor keeps a degenerate panel
    from collapsing to a zero-height axis.
    """
    vmax = 0.0
    for _tool, g in curve(per_unit, species, metric).items():
        vmax = max(vmax, float(pd.to_numeric(g["value"]).max()))
    return max(vmax, 0.05)


def _draw_panel(ax, per_unit: pd.DataFrame, i: int, j: int, *, legend: bool,
                standalone: bool) -> None:
    """One species x metric panel; ``standalone`` repeats the axis labels for a
    single-panel file, ``legend`` adds the tool key (panel A only, on both)."""
    metric, ylab = METRIC_COLS[j]
    species, _group, units = PANELS[i]
    want_units = len(units)
    for tool in TOOL_ORDER:
        g = curve(per_unit, species, metric).get(tool)
        if g is None:
            continue
        #: 2026-10-10: every point must average one row per unit of the condition -- this is
        #: the gate that makes "means over the independent sequencing units" a checked claim
        if not (g["n_units"] == want_units).all() or not (g["n_rows"] == want_units).all():
            bad = g[(g["n_units"] != want_units) | (g["n_rows"] != want_units)]
            raise AssertionError(f"{species}/{tool}/{metric}: a point averages "
                                 f"{sorted(set(bad['n_units']))} unit(s) instead of "
                                 f"{want_units} (thresholds {list(bad['threshold'])[:4]})")
        x = 100.0 * pd.to_numeric(g["retained_fraction"]).to_numpy(dtype=float)
        y = pd.to_numeric(g["value"]).to_numpy(dtype=float)
        ax.plot(x, y, "-", marker=TOOL_MARKER[tool], ms=3.6, lw=1.5,
                color=TOOL_COLOR[tool], label=tool)
        d = g["is_default"].to_numpy(dtype=bool)
        if d.any():
            ax.plot(x[d], y[d], marker=TOOL_MARKER[tool], ms=7.5, ls="none",
                    mfc="white", mec=TOOL_COLOR[tool], mew=1.4)
    top = panel_vmax(per_unit, species, metric)
    ylim = Y_LIM_HEADROOM * top
    ax.set_ylim(0, ylim)
    ax.set_xlim(0, 103)
    #: 2026-10-10 (transpose round): pin the five x numbers -- on the narrow 3-column panels the
    #: automatic locator dropped 25/50/75 (measured: only 0 and 100 were printed).
    #: 2026-10-10 (legibility round): and give every panel the house tick size; the shared
    #: figstyle leaves them at ``font.size = 13``, which is what collided on the delivered page.
    ax.set_xticks(X_TICKS)
    ax.tick_params(axis="both", labelsize=TICK_FONTSIZE)
    #: 2026-10-10 (free-scale round): the locator's nice values inside *this* panel's own limit,
    #: written to two decimals.  The old fixed set (0.20 on a 0.3686-tall F1 axis) left 46 % of
    #: the panel tick-free -- the blank band the author objected to -- so the assert below is that
    #: complaint's direct inverse: the panel's own curves must fill >= 92 % of the drawn height.
    ax.yaxis.set_major_locator(MaxNLocator(nbins=Y_TICK_NBINS, min_n_ticks=3))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ticks = [float(t) for t in ax.get_yticks() if t <= ylim + 1e-9]
    if not 3 <= len(ticks) <= 8:
        raise AssertionError(f"{species}/{metric}: {len(ticks)} y ticks inside "
                             f"[0, {ylim:.4f}]")
    if top / ylim < 0.92:
        raise AssertionError(f"{species}/{metric}: the curves reach only {100 * top / ylim:.0f} "
                             f"% of the axis height (limit {ylim:.4f})")
    #: 2026-10-10 (panel-labels round, author decision): every panel carries its own y title -- the style
    #: of Figure S9/S10, which repeat the same axis title panel by panel.

    #: the same holds for the **x** axis now -- all six panels print `X_LABEL` and their own
    #: 0/25/50/75/100 numbers, so no reader has to look two rows down for either.
    ax.set_ylabel(ylab, fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_xlabel(X_LABEL, fontsize=AXIS_LABEL_FONTSIZE)
    #: ``sharex`` hides every row's numbers but the last; the page wants all three, so each
    #: axes re-enables its own bottom labels (the shared limits stay linked as before)
    ax.tick_params(labelbottom=True)
    if standalone:
        #: a single-panel file is read on its own, so it repeats the species name as its title
        ax.set_title(species_label(i), fontweight="bold")
    if legend:
        #: the key never sits on the data: measured 2026-10-10, every panel of this page is
        #: crossed by curves (83-100 % of the axis height where a key would have to go), so it
        #: goes below the axes -- here for a single-panel file, once for the whole page in
        #: ``figure()``.  Fixed figure coordinates + a reserved band (see ``LEGEND_BAND_PANEL``)
        #: instead of an "outside" locator, which drifts during the save draw.
        handles, labels = ax.get_legend_handles_labels()
        #: centred in the reserved band (``loc="center"``, not "lower center"): anchored at its
        #: bottom edge the two-row key grew up into the x title -- measured on panel A1 in the
        #: legibility round, three collisions with `Tool calls retained (%)`
        ax.figure.legend(handles, labels, loc="center",
                         bbox_to_anchor=(0.5, 0.5 * LEGEND_BAND_PANEL),
                         bbox_transform=ax.figure.transFigure, ncol=3,
                         frameon=True, facecolor="white", edgecolor="0.45",
                         framealpha=1.0, borderpad=0.32, fontsize=9, handlelength=1.4,
                         columnspacing=0.9)


def _add_headers(fig, axes) -> tuple[list, list]:
    """The page's two row letters and its three species column headers.

    *Letters*: two bold capitals, one per **metric row** -- A = the PPV row, B = the F1 row --
    drawn in the left band the layout reserves (``LETTER_BAND``), outside every axes, so the
    caption's "(A)"/"(B)" keep naming what they named before the 2026-10-10 transpose.
    *Column headers*: one species label (with its ``n``) centred over that column's **two**
    panels, the transposed twin of the earlier row headers, so one name still serves both
    metrics of that species (the author's standing request).  Everything is figure text, which
    keeps the in-panel annotation gate meaningful (``ax.texts`` must stay empty), and it is laid
    out **after** ``fig.canvas.draw()`` because with constrained_layout the axes only reach
    their final position at draw time -- and because the shared keys must already be attached,
    since they move the panels.

    2026-10-10 history: one letter per tile inside every panel (rejected); two block letters in
    the top margin naming the metric columns (kept for several rounds); a "table" variant with
    metric column headers (rejected); per-panel axis titles with three species row headers; and
 on the author's "2*3 legend" the page turned 90 degrees -- so the letters name the
    metric **rows** and the species labels became the **column** headers.
    """
    fig.canvas.draw()                       # constrained_layout resolves here
    nrows, ncols = len(METRIC_COLS), len(PANELS)
    letters = []
    for j in range(nrows):

        #: at that block's **top-left corner** (the house wording since 2026-09-21), i.e. just
        #: above the first panel of its own row, still inside the leftmost margin band.  It used
        #: to sit at the vertical middle of the left margin, which read as a row label.
        bb = axes[j][0].get_position()
        letters.append(fig.text(LETTER_X, bb.y1 + HEADER_GAP, BLOCK_LETTERS[j],
                                fontsize=HEADER_FONTSIZE["letter"], fontweight="bold",
                                ha="left", va="bottom"))
    #: the column headers sit in a band above the first row, each centred over its column
    y = axes[0][0].get_position().y1 + HEADER_GAP
    cols = []
    for i in range(ncols):
        left = axes[0][i].get_position().x0
        right = axes[nrows - 1][i].get_position().x1
        cols.append(fig.text((left + right) / 2.0, y, species_label(i),
                             fontsize=HEADER_FONTSIZE["row"], fontweight="bold",
                             ha="center", va="bottom"))
    return letters, cols


def _save_padded(fig, stem: str) -> None:
    """``figstyle.save`` with this page's ``MARGIN_PT`` clear points on every side.

    2026-10-10 (margins round): the house helper hard-codes ``pad_inches=0.02`` (~1.4 pt) and
    is imported read-only from the frozen ``$RNAMODBENCH_LOCAL/code`` tree by 60+ scripts, so its default stays
    as it is; this page needs a real margin because the saved PDF measured left -0.2 pt / top
    -0.5 pt (ink past the page edge) with it.  PDF + 300 dpi PNG, exactly like the helper.
    """
    pad = MARGIN_PT / 72.0
    fig.savefig(stem + ".pdf", bbox_inches="tight", pad_inches=pad)
    fig.savefig(stem + ".png", bbox_inches="tight", pad_inches=pad, dpi=300)
    plt.close(fig)


def _add_legends(fig, axes, handles, labels) -> list:
    """One framed tool key **per metric row**, at fixed figure coordinates.

 2026-10-10 (transpose round, author note: "2*3 legend"): every row carries its own key, so
    a reader of the F1 strip never has to travel back up to the PPV strip's key.  The last row's
    key sits in the band the layout engine reserves (``LEGEND_BAND``); the one above it sits in
    the middle gap, between that row's x titles and the next row's panels, measured after
    ``fig.canvas.draw()``.  Fixed coordinates matter: in the margins round the "outside" locator
    was measured to drift ~6 pt during the save draw and eat the bottom page margin (4.08 pt
    instead of 10).
    """
    fig.canvas.draw()
    nrows = len(METRIC_COLS)
    #: centred on the panels, not on the figure -- the left band shifts the axes right
    x = 0.5 * (axes[0][0].get_position().x0
               + axes[0][len(PANELS) - 1].get_position().x1)
    inv = fig.transFigure.inverted()
    legends = []
    for j in range(nrows):
        if j == nrows - 1:
            y = 0.5 * LEGEND_BAND
        else:
            top = min(ax.xaxis.label.get_window_extent().y0 for ax in axes[j])
            bottom = max(ax.get_window_extent().y1 for ax in axes[j + 1])
            y = inv.transform((0.0, 0.5 * (top + bottom)))[1]
        legend = fig.legend(handles, labels, loc="center", bbox_to_anchor=(x, y),
                            bbox_transform=fig.transFigure, ncol=len(TOOL_ORDER),
                            frameon=True, facecolor="white", edgecolor="0.45",
                            framealpha=1.0, borderpad=0.32, fontsize=9, handlelength=1.4,
                            columnspacing=0.9)
        legend.get_frame().set_linewidth(0.7)
        legends.append(legend)
    return legends


def figure(per_unit: pd.DataFrame) -> float:
    apply_style()

    #: (PPV above F1) and columns the **species** (Arabidopsis, mouse, HeLa) -- the 3x2 page
    #: turned 90 degrees -- and every row carries its own key.  ``sharex`` still links each
    #: column, and every panel re-enables its own bottom tick labels (x-labels round).
    fig, axes = plt.subplots(len(METRIC_COLS), len(PANELS), figsize=FIGSIZE,
                             constrained_layout=True, sharex=True)
    #: the middle gap must hold, top to bottom: the upper row's x tick numbers and x title (both
    #: inside that axes' tight box, so constrained_layout reserves them automatically) and then
    #: that row's key (placed after the draw by ``_add_legends``, which is why hspace has to be
    #: generous).  ``rect`` additionally keeps a left band for the bold row letters and a bottom
    #: band for the last key, so both have fixed figure coordinates and cannot drift on save.
    fig.get_layout_engine().set(hspace=BETWEEN_ROWS_HSPACE,
                                rect=(LETTER_BAND, LEGEND_BAND, 1.0 - LETTER_BAND,
                                      1.0 - LEGEND_BAND))
    for j in range(len(METRIC_COLS)):
        for i in range(len(PANELS)):
            _draw_panel(axes[j][i], per_unit, i, j, legend=False, standalone=False)
    #: the per-panel axis titles are set in ``_draw_panel`` (house style); here the page adds
    #: the per-row keys -- framed, per the 2026-09-30 rule that a frame keeps a key from
    #: colliding with the artwork
    handles, labels = axes[0][0].get_legend_handles_labels()
    if list(labels) != TOOL_ORDER:
        raise AssertionError(f"legend entries {list(labels)} != TOOL_ORDER")
    legends = _add_legends(fig, axes, handles, labels)
    letters, cols = _add_headers(fig, axes)          # only now is the layout final
    min_font = _assert_style(fig, letters=letters, below=legends, cols=cols)
    FIG.mkdir(parents=True, exist_ok=True)
    stem = FIG / "FigureS11_rev"
    _save_padded(fig, str(stem))
    scales = ", ".join(f"{species}: " + "/".join(f"{m} <{Y_LIM_HEADROOM * panel_vmax(per_unit, species, m):.3f}"
                                                 for m, _lab in METRIC_COLS)
                       for species, _g, _u in PANELS)
    print(f"wrote {stem.with_suffix('.pdf').name} + .png (min font {min_font:.1f} pt; "
          f"{len(cols)} species column headers, {len(BLOCK_LETTERS)} row letters "
          f"{'/'.join(BLOCK_LETTERS)} outside the axes, per-panel axis titles with the x title "
          f"+ its own x tick numbers on every panel, {len(legends)} framed {len(labels)}-entry "
          f"keys (one per metric row), free per-panel y scales [{scales}], "
          f"{MARGIN_PT:g} pt page margin on every side)")
    return min_font


def panels(per_unit: pd.DataFrame) -> list[Path]:
    """The six panels as single-panel files, so the delivery can carry the same
    per-panel set as every other figure (``figures/panels/FigureS11_panelA1.pdf`` ...).

    A tile keeps a letter of its own, but named after its **block** plus the row, so the
    single-panel files map onto the page's A/B blocks (A1-A3 PPV, B1-B3 F1; author decision
    2026-10-10).  The first tile carries the tool key, below the axes like on the page.
    """
    out_dir = FIG / "panels"
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    labels = panel_labels()
    for i in range(len(PANELS)):
        for j in range(len(METRIC_COLS)):
            k = i * len(METRIC_COLS) + j
            apply_style()
            fig, ax = plt.subplots(1, 1, figsize=(3.4, 2.8), constrained_layout=True)
            #: reserve the band the two-row key sits in, exactly as on the page (see LEGEND_BAND)
            fig.get_layout_engine().set(rect=(0.0, LEGEND_BAND_PANEL, 1.0,
                                              1.0 - LEGEND_BAND_PANEL))
            _draw_panel(ax, per_unit, i, j, legend=(k == 0), standalone=True)
            #: a single-panel file must still name both axes (house style: Figure 3 / S8 / S10)
            if ax.get_ylabel() != METRIC_COLS[j][1] or ax.get_xlabel() != X_LABEL:
                raise AssertionError(f"panel {labels[k]}: titles {ax.get_ylabel()!r} / "
                                     f"{ax.get_xlabel()!r} do not match the house style")
            fig.text(0.02, 0.99, labels[k], fontsize=11, fontweight="bold",
                     ha="left", va="top")
            fig.canvas.draw()           # the gate reads display boxes, so draw first
            _assert_style(fig, below=list(fig.legends))
            stem = out_dir / f"FigureS11_panel{labels[k]}"
            _save_padded(fig, str(stem))    # same 10 pt margin as the six-panel page
            written.append(stem.with_suffix(".pdf"))
    print("wrote " + ", ".join(p.name for p in written))
    return written


def main() -> None:
    logger = setup_logger("figS11_figure", log_dir=OUT / "logs")
    per_unit = pd.read_csv(PER_UNIT, sep="\t")
    for species, _g, _u in PANELS:
        for metric, _lab in METRIC_COLS:
            for tool, g in curve(per_unit, species, metric).items():
                if g["value"].isna().all():
                    raise AssertionError(f"{species}/{tool}/{metric}: empty curve")
        if not (per_unit["retained_fraction"].between(0, 1).all()):
            raise AssertionError("retained_fraction outside [0, 1]")
    logger.info("data integrity ok: 3 species x %d tools x %d metrics, each point the mean "
                "over that condition's units at one threshold",
                len(TOOL_ORDER), len(METRIC_COLS))
    figure(per_unit)
    panels(per_unit)
    logger.info("wrote FigureS11_rev.pdf (no gridlines, no in-panel text; "
                "%d panels in %d metric rows x %d species columns, %d per-panel y titles with "
                "free per-panel y scales filling >= 92 %% of each axis, the x title and its own "
                "x tick numbers on every panel, %d species column headers centred over their two "
                "panels, %s row letters outside the axes, one framed key below each row, "
                "%g pt page margin on every side) and panels %s",
                len(PANELS) * len(METRIC_COLS), len(METRIC_COLS), len(PANELS),
                len(PANELS) * len(METRIC_COLS), len(PANELS),
                "/".join(BLOCK_LETTERS), MARGIN_PT, ", ".join(panel_labels()))
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
