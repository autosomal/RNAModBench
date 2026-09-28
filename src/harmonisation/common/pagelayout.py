"""Shared page layout for the Supporting-Information figure pages (Figure S4 and Figure S5).

Why this module exists
----------------------
The first portrait draft of Figure S4 was rejected twice for layout reasons,
and both defects were pure geometry that a hand-tuned layout cannot see:

* the panel letter "A" was drawn with ``ax.text(x=-0.18, y=1.30)`` in *axes*
  coordinates -- on a log axis the top-most y tick label sits in exactly that
  spot, so the letter was printed on top of ``10^4``;
* the y tick labels of the second and third species column ran into the
  drawing area of the column to their left, because the column gap was chosen
  by eye instead of from the measured label width.

This module removes the class of error instead of the instance:

* :data:`A4_LANDSCAPE_IN` -- the page every SI figure page uses now;
* :func:`text_width_in` / :func:`rotated_label_in` -- measure the real ink
  width of a label with the current rcParams font, so the space a row needs is
  computed, not guessed;
* :func:`wspace_for_gap` -- turn a required gap in inches into the GridSpec
  ``wspace`` that leaves exactly that gap;
* :func:`margin_letter` -- panel letters live in the page margin, where no
  tick label, title or axis label can ever be, instead of at an axes-relative
  offset;
* :func:`assert_page_clean` -- after drawing, every text artist is compared
  with every other one *across axes* and against the drawing area of every
  other axes, plus a grid and minimum-font-size audit.  A collision aborts the
  render.

The gate is deliberately blunt: it fails loudly rather than shipping a page
that only looks right on the author's screen.
"""
from __future__ import annotations

import itertools
import sys
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PathCollection, PolyCollection
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.textpath import TextPath
from matplotlib.transforms import Bbox

__all__ = [
    "A4_LANDSCAPE_IN",
    "A4_PORTRAIT_IN",
    "LETTER_X",
    "text_width_in",
    "rotated_label_in",
    "wspace_for_gap",
    "margin_letter",
    "text_quad",
    "rect_quad",
    "quads_intersect",
    "assert_legend_clear",
    "assert_page_clean",
]

#: page sizes in inches; every Supporting-Information page is A4 landscape
#: since 2026-09-21 (four content rows x three species columns do not fit a
#: portrait sheet once the tick labels are large enough to be printed)
A4_LANDSCAPE_IN: tuple[float, float] = (11.69, 8.27)
#: kept for reference -- the submitted (pre-revision) ``sup4.pdf`` page
A4_PORTRAIT_IN: tuple[float, float] = (8.27, 11.70)

#: x position of every panel letter, in page fractions.  The strip from the
#: page edge to the axes' left spine is reserved for the letter, the y axis
#: label and the y tick labels, in that order -- the letter has to stay clear
#: of the rotated y axis label (0.29 in wide on the Figure S4 and Figure S5 pages), so it sits
#: 0.105 in from the page edge.
LETTER_X: float = 0.009
#: how far above a row's top edge the letter starts (inches)
LETTER_GAP_IN: float = 0.02


# --------------------------------------------------------------------------- #
# measuring text
# --------------------------------------------------------------------------- #
def _font_prop(weight: str = "normal") -> FontProperties:
    family = plt.rcParams.get("font.family", ["Arial"])
    if isinstance(family, str):
        family = [family]
    return FontProperties(family=list(family), weight=weight)


def text_width_in(text: str, fontsize: float, weight: str = "normal") -> float:
    """Ink width of ``text`` in inches when printed at ``fontsize`` points.

    Measured with the same (Arial) font matplotlib will embed, so the numbers
    are directly comparable with the inches a page budget is written in.
    """
    path = TextPath((0.0, 0.0), text, size=fontsize, prop=_font_prop(weight))
    return float(path.get_extents().width) / 72.0


def rotated_label_in(labels: Iterable[str], fontsize: float) -> float:
    """Space a 90-degree rotated x tick label needs below its axis (inches)."""
    return max((text_width_in(s, fontsize) for s in labels), default=0.0)


def wspace_for_gap(fig_width_in: float, *, left: float, right: float, ncols: int,
                   gap_in: float) -> float:
    """``wspace`` value that leaves exactly ``gap_in`` inches between columns.

    ``wspace`` is expressed as a fraction of the average *axes* width, so it
    has to be computed from the page geometry instead of being tuned by eye.
    The gap has to hold the neighbouring column's y tick labels.
    """
    usable = (right - left) * fig_width_in
    column = (usable - (ncols - 1) * gap_in) / ncols
    if column <= 0:
        raise ValueError(f"gap {gap_in:.2f} in does not fit {ncols} columns "
                         f"into {usable:.2f} in")
    return gap_in / column


# --------------------------------------------------------------------------- #
# panel letters in the page margin
# --------------------------------------------------------------------------- #
def margin_letter(fig: plt.Figure, letter: str, *, y: float,
                  x: float = LETTER_X, fontsize: float = 20.0,
                  gap_in: float = LETTER_GAP_IN) -> Text:
    """Draw a bold panel letter in the left page margin.

    ``y`` is the page-fraction height of the row's top edge (typically
    ``ax.get_position().y1``).  The letter is placed left of the y axis label,
    so it can never land on a tick label (see the module docstring).
    """
    return fig.text(x, y + gap_in / fig.get_figheight(), letter, ha="left",
                    va="top", fontsize=fontsize, fontweight="bold")


# --------------------------------------------------------------------------- #
# rotated ink geometry
#
# A 45-degree x tick label is a *rotated* rectangle.  Its axis-aligned bounding
# box is a much larger square, so comparing bounding boxes would flag overlaps
# that are never printed (and, worse, would tempt anyone to weaken the gate).
# The gate therefore works on quads.
# --------------------------------------------------------------------------- #
def _layout_size(text: Text, renderer) -> tuple[float, float]:
    """Unrotated layout size of ``text`` in display pixels (ascent+descent)."""
    angle = text.get_rotation()
    if angle == 0:
        box = text.get_window_extent(renderer)
        return box.width, box.height
    text.set_rotation(0)                    # measuring only; restored below
    try:
        box = text.get_window_extent(renderer)
    finally:
        text.set_rotation(angle)
    return box.width, box.height


def _ink_rect(text: Text, renderer) -> tuple[float, float, float, float] | None:
    """Ink rectangle inside the layout box, in display-pixel coordinates.

    A 45-degree tick label is placed by its *layout* box (ascent + descent), but
    only the glyphs are printed.  For the dense tool names of Figure S4 panel A
    the difference decides whether two neighbouring labels merely touch their
    boxes or really collide, so the gate measures the ink when it can: the
    string is laid out with the same font through ``TextPath`` and mapped into
    the layout box (whose left edge is the text origin and whose baseline sits
    ``descent`` above its bottom edge).

    Returns ``None`` when the string cannot be measured this way (mathtext, an
    empty font); the caller then keeps the conservative layout box.
    """
    label = text.get_text().strip()
    if not label or "$" in label:            # mathtext / no ink path
        return None
    try:
        prop = text.get_fontproperties()
        extent = TextPath((0.0, 0.0), label, size=prop.get_size_in_points(),
                          prop=prop).get_extents()
    except Exception:                        # font without an outline: keep box
        return None
    scale = renderer.points_to_pixels(1.0)   # points -> display pixels
    # Arial's descent is 434/2048 em; the layout box is ascent + descent, so the
    # baseline sits ``descent`` above the box bottom
    baseline = 0.212 * prop.get_size_in_points() * scale
    return (float(extent.x0 * scale), float(extent.y0 * scale + baseline),
            float(extent.x1 * scale), float(extent.y1 * scale + baseline))


def text_quad(text: Text, renderer, *, ink: bool = True,
              ) -> tuple[tuple[float, float], ...]:
    """Corners of the rotated box of ``text``, in display pixels.

    matplotlib has two placement rules and both are used on these pages (tick
    labels and titles are ``"default"``, y axis labels are ``"anchor"``):

    * ``"default"`` rotates the string first and then aligns the *rotated*
      bounding box to the anchor, so the quad is the text rectangle rotated
      about the centre of that aligned rotated box;
    * ``"anchor"`` aligns the *unrotated* box to the anchor and then rotates it
      about the anchor -- that is what keeps a y axis label centred on its axis.

    Inside that rectangle the ink rectangle is used when it can be measured
    (see :func:`_ink_rect`); the layout box stays the fallback.
    """
    theta = np.deg2rad(text.get_rotation() % 360.0)
    cos, sin = float(np.cos(theta)), float(np.sin(theta))
    w, h = _layout_size(text, renderer)
    anchor = text.get_transform().transform(text.get_position())
    ha, va = text.get_ha(), text.get_va()
    centred = ("center", "center_baseline")
    if text.get_rotation_mode() == "anchor":
        x0 = anchor[0] - (w if ha == "right" else w / 2.0 if ha in centred else 0.0)
        y0 = anchor[1] - (h if va == "top" else h / 2.0 if va in centred else 0.0)
        cx, cy = anchor
    else:
        rot_w = abs(w * cos) + abs(h * sin)      # size of the rotated bbox
        rot_h = abs(w * sin) + abs(h * cos)
        bx0 = anchor[0] - (rot_w if ha == "right"
                           else rot_w / 2.0 if ha in centred else 0.0)
        by0 = anchor[1] - (rot_h if va == "top"
                           else rot_h / 2.0 if va in centred else 0.0)
        cx, cy = bx0 + rot_w / 2.0, by0 + rot_h / 2.0
        x0, y0 = cx - w / 2.0, cy - h / 2.0
    ink_box = _ink_rect(text, renderer) if ink else None
    if ink_box is None:
        rx0, ry0, rx1, ry1 = x0, y0, x0 + w, y0 + h
    else:
        rx0 = x0 + min(max(ink_box[0], 0.0), w)
        ry0 = y0 + min(max(ink_box[1], 0.0), h)
        rx1 = x0 + min(max(ink_box[2], 0.0), w)
        ry1 = y0 + min(max(ink_box[3], 0.0), h)
    corners = ((rx0, ry0), (rx1, ry0), (rx1, ry1), (rx0, ry1))
    return tuple((cx + (px - cx) * cos - (py - cy) * sin,
                  cy + (px - cx) * sin + (py - cy) * cos)
                 for px, py in corners)


def _data_vertices(ax: plt.Axes) -> list[tuple[float, float]]:
    """Display coordinates of every plotted vertex/point drawn in ``ax``.

    Used by :func:`assert_legend_clear`: the project's rule is that a legend
    hiding data is judged by *whether data vertices fall inside the legend box*,
    never by comparing bounding boxes (a long diagonal or a curve has a bounding
    box that overlaps almost everything).
    """
    points: list[tuple[float, float]] = []
    legend = ax.get_legend()
    for artist in ax.get_children():
        if artist is legend:
            continue
        children = getattr(artist, "lines", None) or ()
        for member in (artist, *children):
            try:
                if isinstance(member, Line2D):
                    data = member.get_xydata()
                    if data is None or len(data) == 0:
                        continue
                    display = member.get_transform().transform(
                        np.asarray(data, dtype=float)[:, :2])
                elif isinstance(member, PathCollection):
                    offsets = member.get_offsets()
                    if offsets is None or len(offsets) == 0:
                        continue
                    display = member.get_transform().transform(
                        np.asarray(offsets, dtype=float)[:, :2])
                elif isinstance(member, (LineCollection, PolyCollection)):
                    chunks = []
                    if isinstance(member, LineCollection):
                        chunks = [seg for seg in member.get_segments() if len(seg)]
                    if not chunks:
                        chunks = [path.vertices for path in member.get_paths()
                                  if len(path.vertices)]
                    if not chunks:
                        continue
                    display = np.vstack([
                        member.get_transform().transform(np.asarray(c, dtype=float)[:, :2])
                        for c in chunks])
                else:
                    continue
            except (AttributeError, TypeError, ValueError):
                continue
            points.extend((float(x), float(y)) for x, y in np.asarray(display))
    return points


def assert_legend_clear(ax: plt.Axes, legend=None, *, pad_pt: float = 2.0,
                        verbose: bool = True) -> None:
    """Abort the render when a legend box covers plotted data.

    ``pad_pt`` shrinks the legend box before the test, so a marker that merely
    touches the border is not a failure.  Only the axes that carry an in-panel
    legend need to be checked.
    """
    legend = legend if legend is not None else ax.get_legend()
    if legend is None:
        return
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    pad = pad_pt * renderer.points_to_pixels(1.0)
    box = legend.get_window_extent(renderer)
    inner = Bbox.from_extents(box.x0 + pad, box.y0 + pad, box.x1 - pad,
                              box.y1 - pad)
    points = _data_vertices(ax)
    inside = [(x, y) for x, y in points
              if inner.x0 <= x <= inner.x1 and inner.y0 <= y <= inner.y1]
    panel = ax.get_title() or "panel"
    if inside:
        raise SystemExit(f"legend over data: {len(inside)} of {len(points)} drawn "
                         f"vertices of {panel} fall inside the legend box "
                         f"({box.x0:.0f},{box.y0:.0f})-({box.x1:.0f},{box.y1:.0f})")
    if verbose:
        print(f"[legend] {panel}: legend box clear of all {len(points)} drawn "
              f"vertices")


def rect_quad(box: Bbox) -> tuple[tuple[float, float], ...]:
    """The four corners of an axis-aligned box (same type as ``text_quad``)."""
    return ((box.x0, box.y0), (box.x1, box.y0), (box.x1, box.y1),
            (box.x0, box.y1))


def _edge_normals(poly: Sequence[tuple[float, float]]):
    n = len(poly)
    for i in range(n):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % n]
        ex, ey = x1 - x0, y1 - y0
        length = float(np.hypot(ex, ey))
        if length > 0:
            yield (-ey / length, ex / length)


def quad_gap(a: Sequence[tuple[float, float]],
             b: Sequence[tuple[float, float]]) -> float:
    """Signed clearance between two convex quads in pixels.

    Positive = the widest separating direction (a real gap), negative = the
    shallowest overlap depth.  Reported by the geometry self-test so a layout
    can be judged by its margin instead of a yes/no answer.
    """
    worst = None
    for axis in itertools.chain(_edge_normals(a), _edge_normals(b)):
        ax, ay = axis
        pa = [x * ax + y * ay for x, y in a]
        pb = [x * ax + y * ay for x, y in b]
        sep = max(min(pa), min(pb)) - min(max(pa), max(pb))
        worst = sep if worst is None else max(worst, sep)
    return float(worst) if worst is not None else float("inf")


def quads_intersect(a: Sequence[tuple[float, float]],
                    b: Sequence[tuple[float, float]], pad: float = 0.0) -> bool:
    """True when two convex quads overlap, or stay closer than ``pad`` pixels.

    ``pad`` gives the same "touching is not overlap" tolerance the previous
    bounding-box comparison had, so a printed gap below ``pad`` is still a
    failure rather than a silent pass.
    """
    return quad_gap(a, b) <= pad


def _text_artists(fig: plt.Figure) -> list[tuple[str, int | None, Text]]:
    """``(owner, id(axes), artist)`` for every visible text on the page."""
    out: list[tuple[str, int | None, Text]] = []
    seen: set[int] = set()

    def add(owner: str, ax_id: int | None, texts: Iterable[Text]) -> None:
        for t in texts:
            if id(t) in seen:
                continue
            seen.add(id(t))
            if not t.get_visible() or not t.get_text().strip():
                continue
            out.append((owner, ax_id, t))

    for i, ax in enumerate(fig.axes):
        title = ax.get_title()
        owner = f"axes{i}" + (f" [{title}]" if title else "")
        if not ax.axison:
            # an axis-off stripe keeps its tick *objects* but never draws them;
            # only the legend it carries is part of the printed page
            leg = ax.get_legend()
            if leg is not None:
                add(f"{owner} legend", id(ax), list(leg.get_texts()))
            continue
        add(owner, id(ax), (ax.title, ax.xaxis.label, ax.yaxis.label))
        if ax.xaxis.get_visible():
            add(owner, id(ax), list(ax.get_xticklabels()))
        if ax.yaxis.get_visible():
            add(owner, id(ax), list(ax.get_yticklabels()))
        leg = ax.get_legend()
        if leg is not None:
            add(f"{owner} legend", id(ax), list(leg.get_texts()))
    add("page", None, list(fig.texts))
    for leg in fig.legends:
        add("page legend", None, list(leg.get_texts()))
    return out


def assert_page_clean(fig: plt.Figure, *, ignore_axes: Sequence[plt.Axes] = (),
                      check_spines: bool = True, pad_px: float = 0.5,
                      near_pt: float = 3.0, min_pt: float = 7.0,
                      verbose: bool = True) -> None:
    """Abort the render if the page is not printable.

    Four checks, all geometric (no eyeballing):

    1. no two text labels overlap anywhere on the page -- this covers the
       neighbouring-column tick labels *and* the cross-row cases (an x axis
       title running into the next row's panel title);
    2. no text enters the drawing area of a panel it does not belong to;
    3. no text sits on -- or *hugs* -- a spine or tick mark of another panel
 (``check_spines``): the defect the author called " font/box/line overlap ". A text
       that merely touches a frame line is invisible to an overlap test but is
       exactly what a reader sees, so the clearance has to be at least
       ``near_pt`` points (0.04 in at the default); the "next column's tick
       labels hugging the previous column's frame" case is caught here too;
    4. no visible grid line and no text below ``min_pt`` points.

    ``ignore_axes`` lists background stripes (axis-off legend bands) that must
    not act as keep-out zones; their own texts are still checked by (1).
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    items = _text_artists(fig)
    quads = [text_quad(t, renderer) for _o, _a, t in items]
    owners = [owner for owner, _a, _t in items]
    labels = [t.get_text() for _o, _a, t in items]
    owners_of_axes = [ax_id for _o, ax_id, _t in items]
    problems: list[str] = []

    for i in range(len(quads)):
        for j in range(i + 1, len(quads)):
            if quads_intersect(quads[i], quads[j], pad=pad_px):
                problems.append(f"text overlap: {owners[i]} {labels[i]!r}"
                                f" <-> {owners[j]} {labels[j]!r}")

    skip = {id(ax) for ax in ignore_axes}
    zones = [(ax, rect_quad(ax.get_window_extent(renderer))) for ax in fig.axes
             if id(ax) not in skip]
    for idx, quad in enumerate(quads):
        for ax, zone in zones:
            if id(ax) == owners_of_axes[idx]:
                continue
            if quads_intersect(quad, zone, pad=pad_px):
                problems.append(f"text enters a foreign panel: {owners[idx]} "
                                f"{labels[idx]!r} -> "
                                f"{ax.get_title() or 'panel'} area")

    for i, ax in enumerate(fig.axes):
        lines = list(ax.get_xgridlines()) + list(ax.get_ygridlines())
        if any(line.get_visible() for line in lines):
            problems.append(f"axes{i}: a grid line is visible")

    if check_spines:
        near = near_pt * renderer.points_to_pixels(1.0)
        for ax in fig.axes:
            if not ax.axison:
                continue
            marks: list[tuple[str, tuple[tuple[float, float], ...]]] = []
            for name, spine in ax.spines.items():
                if spine.get_visible():
                    marks.append((f"{name} spine",
                                  rect_quad(spine.get_window_extent(renderer))))
            ticks = list(ax.xaxis.get_major_ticks()) + list(ax.yaxis.get_major_ticks())
            for tick in ticks:
                for line in (tick.tick1line, tick.tick2line):
                    if line.get_visible():
                        marks.append(("tick mark",
                                      rect_quad(line.get_window_extent(renderer))))
            if not marks:
                continue
            panel = ax.get_title() or "panel"
            for idx, quad in enumerate(quads):
                for label, mark in marks:
                    if quads_intersect(quad, mark, pad=pad_px):
                        problems.append(f"text on a frame line: {owners[idx]} "
                                        f"{labels[idx]!r} -> {label} of {panel}")
                        continue
                    gap = quad_gap(quad, mark) / renderer.points_to_pixels(1.0)
                    if gap < near_pt:
                        problems.append(f"text hugs a frame line ({gap:.1f} pt): "
                                        f"{owners[idx]} {labels[idx]!r} -> "
                                        f"{label} of {panel}")

    sizes = [(owner, text, artist.get_fontsize())
             for (owner, _a, text), (_o, _ai, artist) in zip(
                 ((o, a, t.get_text()) for (o, a, t) in items), items)]
    for owner, text, size in sizes:
        if size < min_pt:
            problems.append(f"font below {min_pt:g} pt: {owner} {text!r} "
                            f"({size:.2f} pt)")

    if problems:
        for line in problems[:40]:
            print(f"  ! {line}", file=sys.stderr)
        raise SystemExit(f"page layout check failed: {len(problems)} problem(s)")
    if verbose:
        smallest = min(size for _o, _t, size in sizes)
        rotated = sum(1 for _o, _a, t in items if t.get_rotation() % 360.0)
        print(f"[layout] {len(quads)} text labels ({rotated} rotated): no "
              f"overlap, no cross-panel intrusion, smallest font "
              f"{smallest:.1f} pt")


# --------------------------------------------------------------------------- #
# self-test: the quad must reproduce matplotlib's own rotated extent, and the
# separation of two adjacent 45-degree tick labels must be measured on the ink
# (bounding boxes of rotated labels overlap even when the printed labels do not)
# --------------------------------------------------------------------------- #
def _selftest() -> int:
    """Run the geometry self-test; returns a process exit code."""
    failures: list[str] = []

    def expect(label: str, condition: bool, detail: str = "") -> None:
        print(f"[{'ok  ' if condition else 'FAIL'}] {label}"
              f"{(' -- ' + detail) if detail else ''}")
        if not condition:
            failures.append(label)

    rc = plt.rcParams
    rc["font.family"] = "Arial"
    for angle, size, slot_in in ((45.0, 8.0, 0.1710), (60.0, 8.5, 0.1646),
                                 (90.0, 8.5, 0.1646)):
        axes_w = slot_in * 13                       # the real panel-A width
        fig_w = axes_w + 0.9                        # room for the ylabel
        fig = plt.figure(figsize=(fig_w, 1.0), dpi=120)
        ax = fig.add_axes([0.9 / fig_w, 0.1, axes_w / fig_w, 0.8])
        ax.set_ylabel("PPV vs. GLORI (2 bp)", fontsize=11.0)   # anchor mode
        ax.set_xlim(0, 13)
        ax.set_xticks(range(13))
        labels = ["NanoSPA_m6A", "EpiNano_Error", "ELIGOS2_solo"] * 4 + ["xPore"]
        ax.set_xticklabels(labels, rotation=angle, ha="right", va="top",
                           fontsize=size)
        ax.tick_params(axis="x", pad=5.0)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()

        # 1. the layout quad must reproduce matplotlib's own rotated extent,
        #    and the ink quad must sit inside it
        worst = 0.0
        outside = 0
        for t in ax.get_xticklabels():
            outer = text_quad(t, renderer, ink=False)
            inner = text_quad(t, renderer)
            if min(p[0] for p in inner) < min(p[0] for p in outer) - 0.5 \
                    or max(p[0] for p in inner) > max(p[0] for p in outer) + 0.5 \
                    or min(p[1] for p in inner) < min(p[1] for p in outer) - 0.5 \
                    or max(p[1] for p in inner) > max(p[1] for p in outer) + 0.5:
                outside += 1
            quad = outer
            xs = [p[0] for p in quad]
            ys = [p[1] for p in quad]
            box = t.get_window_extent(renderer)
            worst = max(worst, abs(min(xs) - box.x0), abs(max(xs) - box.x1),
                        abs(min(ys) - box.y0), abs(max(ys) - box.y1))
        expect(f"{angle:g} deg: quad matches matplotlib's rotated extent",
               worst <= 1.0, f"max deviation {worst:.2f} px")
        expect(f"{angle:g} deg: the ink quad stays inside the layout quad",
               outside == 0, f"{outside} label(s) outside")

        # the y axis label is the "anchor" rotation mode -- the other rule
        ylab = ax.yaxis.label
        qy = text_quad(ylab, renderer, ink=False)
        by = ylab.get_window_extent(renderer)
        dev = max(abs(min(p[0] for p in qy) - by.x0),
                  abs(max(p[0] for p in qy) - by.x1),
                  abs(min(p[1] for p in qy) - by.y0),
                  abs(max(p[1] for p in qy) - by.y1))
        expect(f"{angle:g} deg: anchor-mode quad matches the y axis label",
               dev <= 1.0, f"max deviation {dev:.2f} px")

        # 2. adjacent labels must not be reported as overlapping, and the real
        #    clearance must be reported so the margin is visible
        quads = [text_quad(t, renderer) for t in ax.get_xticklabels()]
        gaps = [quad_gap(quads[i], quads[i + 1]) for i in range(len(quads) - 1)]
        hits = [g for g in gaps if g <= 0.5]
        expect(f"{angle:g} deg: {len(labels)} labels at {slot_in:.4f} in slots "
               f"do not overlap", not hits,
               f"tightest gap {min(gaps):.2f} px "
               f"({min(gaps) / rc['figure.dpi'] * 72:.2f} pt)")

        # 3. the SAT decision itself: a neighbour pushed 6 px closer must be
        #    flagged, a label 60 px away must not be (the tick machinery owns
        #    the real label positions, so the quads are shifted numerically)
        q0 = text_quad(ax.get_xticklabels()[0], renderer)
        q1 = text_quad(ax.get_xticklabels()[1], renderer)
        closer = tuple((x - 10.0, y) for x, y in q1)
        farther = tuple((x + 60.0, y) for x, y in q1)
        expect(f"{angle:g} deg: a deliberate collision is detected",
               quads_intersect(q0, closer, pad=0.5))
        expect(f"{angle:g} deg: a distant label is not flagged",
               not quads_intersect(q0, farther, pad=0.5))
        plt.close(fig)

    if failures:
        print(f"\n{len(failures)} failing self-test(s): {failures}")
        return 2
    print("\nlayout geometry self-test passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(_selftest())
