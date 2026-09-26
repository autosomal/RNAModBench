"""Per-panel rendering and page composition ("draw one panel at a time").

Why this module exists
----------------------
The first two reworks of Figure S4 tried to fit four rows x three species
columns into a single matplotlib page.  Every attempt squeezed the panels into
whatever space was left after the page budget, which is how the page ended up
with letterbox-shaped panels (3.22 x 0.85 in) and with sub-row titles sharing a
band with the tick labels of the row above.

This module inverts the workflow, as the user asked (" draw panels individually then compose --
draw each panel first, then assemble them"): every panel is drawn on its **own
canvas at its final print size**, keeps its own page margin (panel letter,
y axis label, tick labels) and passes the layout gate on its own.  The page is
then assembled by translating the panel PDFs onto one blank page with pypdf,
which keeps the artwork vectorial.

Because a panel never shares a canvas with its neighbours, cross-panel
collisions are impossible by construction: the only thing that can go wrong at
assembly time is the arithmetic of the placement rectangles, and
:func:`compose_page` asserts that they do not overlap.
"""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
from pypdf import PageObject, PdfReader, PdfWriter, Transformation

from common import pagelayout

__all__ = ["Panel", "new_panel", "save_panel", "compose_page", "pdf_to_png"]

#: PDF user-space units per inch
PT_PER_IN = 72.0


@dataclass(frozen=True)
class Panel:
    """One finished panel: its vector file and its exact print size."""

    name: str
    path: Path
    width_in: float
    height_in: float

    @property
    def size_in(self) -> tuple[float, float]:
        return (self.width_in, self.height_in)


def new_panel(width_in: float, height_in: float) -> plt.Figure:
    """A blank canvas at the exact final print size.

    The canvas is *not* cropped on export: a panel's page margin (panel letter,
    y axis label, tick labels) is part of its own width, so stacking panels is
    a pure translation and the left margins of all panels line up.
    """
    return plt.figure(figsize=(width_in, height_in))


def save_panel(fig: plt.Figure, path: Path, *, ignore_axes: Sequence[plt.Axes] = (),
               gate: bool = True, verbose: bool = True) -> Panel:
    """Run the layout gate on one panel and write it as a single-page PDF."""
    if gate:
        if verbose:
            print(f"[panel] {path.stem}")
        pagelayout.assert_page_clean(fig, ignore_axes=ignore_axes, verbose=verbose)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(path))
    size = (fig.get_figwidth(), fig.get_figheight())
    plt.close(fig)
    return Panel(name=path.stem, path=path, width_in=size[0], height_in=size[1])


def _overlaps(a: tuple[float, float, float, float],
              b: tuple[float, float, float, float]) -> bool:
    ax0, ay0, ax1, ay1 = a
    bx0, by0, bx1, by1 = b
    return ax0 < bx1 and ax1 > bx0 and ay0 < by1 and ay1 > by0


def compose_page(page_size_in: tuple[float, float],
                 placements: Iterable[tuple[Panel, float, float]],
                 out_pdf: Path, *, verbose: bool = True) -> Path:
    """Translate every panel onto one blank page (vector, no rasterisation).

    ``placements`` is ``[(panel, x_in, y_in)]`` where ``(x_in, y_in)`` is the
    lower-left corner of the panel on the page.  The placement rectangles are
    asserted to be inside the page and pairwise disjoint, so an assembly bug
    fails loudly instead of silently printing one panel on top of another.
    """
    page_w = page_size_in[0] * PT_PER_IN
    page_h = page_size_in[1] * PT_PER_IN
    page = PageObject.create_blank_page(width=page_w, height=page_h)
    boxes: list[tuple[str, tuple[float, float, float, float]]] = []

    for panel, x_in, y_in in placements:
        reader = PdfReader(str(panel.path))
        if len(reader.pages) != 1:
            raise SystemExit(f"{panel.name}: expected a single-page panel, "
                             f"got {len(reader.pages)}")
        src = reader.pages[0]
        src_w = float(src.mediabox.width)
        src_h = float(src.mediabox.height)
        for got, want, axis in ((src_w, panel.width_in * PT_PER_IN, "width"),
                                (src_h, panel.height_in * PT_PER_IN, "height")):
            if abs(got - want) > 1.0:
                raise SystemExit(f"{panel.name}: PDF {axis} {got:.1f} pt does not "
                                 f"match the declared {want:.1f} pt")
        x = x_in * PT_PER_IN
        y = y_in * PT_PER_IN
        box = (x, y, x + src_w, y + src_h)
        if box[0] < -0.5 or box[1] < -0.5 or box[2] > page_w + 0.5 \
                or box[3] > page_h + 0.5:
            raise SystemExit(f"{panel.name}: placement {box} leaves the "
                             f"{page_w:.0f} x {page_h:.0f} pt page")
        boxes.append((panel.name, box))
        page.merge_transformed_page(src, Transformation().translate(tx=x, ty=y))

    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            if _overlaps(boxes[i][1], boxes[j][1]):
                raise SystemExit(f"panels overlap on the page: {boxes[i][0]} "
                                 f"<-> {boxes[j][0]}")
    if not boxes:
        raise SystemExit("no panel placements given")

    writer = PdfWriter()
    writer.add_page(page)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    with open(out_pdf, "wb") as handle:
        writer.write(handle)
    if verbose:
        print(f"[page] {out_pdf.name}: {len(boxes)} panels on "
              f"{page_size_in[0]:.2f} x {page_size_in[1]:.2f} in")
    return out_pdf


def pdf_to_png(pdf: Path, png: Path, *, dpi: int = 300) -> None:
    """Rasterise the *assembled* page, so the PNG cannot drift from the PDF."""
    subprocess.run(["pdftoppm", "-png", "-r", str(dpi), "-singlefile",
                    str(pdf), str(png.with_suffix(""))], check=True)
