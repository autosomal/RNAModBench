#!/usr/bin/env python3
"""Render the Supplementary Tables S1-S11 into one A4-landscape PDF.

Layout mirrors the previously submitted `Supplementary_Table.pdf` (ReportLab,
A4 landscape 841.89 x 595.276 pt, centred bold "Table Sn: ..." heading,
horizontal rules only, very light row shading, Helvetica).

Inputs : tables/TableS*.tsv written by `build_supp_tables.py` (the first line is
         a `# title`, a `# note:` line may close the file).
Output : Supplementary_Tables.pdf (+ per-table PDFs in figures/)

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/tables/render_supp_tables.py
"""
from __future__ import annotations

from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.platypus import (BaseDocTemplate, Frame, PageTemplate, Paragraph,
                                Spacer, Table, TableStyle)

HERE = Path(__file__).resolve().parent
TAB = HERE / "tables"
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)

PAGE = landscape(A4)
MARGIN = 40
AVAIL = PAGE[0] - 2 * MARGIN

TITLE = ParagraphStyle("title", fontName="Helvetica-Bold", fontSize=11, leading=14,
                       alignment=1, spaceAfter=10)
HEAD = ParagraphStyle("head", fontName="Helvetica-Bold", fontSize=8, leading=9.5,
                      alignment=1)
CELL = ParagraphStyle("cell", fontName="Helvetica", fontSize=7.5, leading=9)
CELLC = ParagraphStyle("cellc", fontName="Helvetica", fontSize=7.5, leading=9,
                       alignment=1)
NOTE = ParagraphStyle("note", fontName="Helvetica-Oblique", fontSize=7, leading=9,
                      spaceBefore=8, alignment=0)
SHADE = colors.Color(0.96, 0.96, 0.96)


def read_table(path: Path) -> tuple[str, list[str], list[list[str]], str]:
    lines = path.read_text().split("\n")
    title = lines[0].lstrip("# ").strip()
    note = ""
    body: list[str] = []
    for ln in lines[1:]:
        if ln.startswith("# note:"):
            note = ln[len("# note:"):].strip()
        elif ln.strip():
            body.append(ln)
    header = body[0].split("\t")
    rows = [r.split("\t") for r in body[1:]]
    return title, header, rows, note


def _word_pt(text: str, font: str, size: float) -> float:
    """Width of the widest unbreakable word of a cell, in points."""
    return max((stringWidth(w, font, size) for w in text.split()), default=0.0)


def _cell_pt(text: str, font: str, size: float) -> float:
    """Width of a cell rendered on one line, in points."""
    return stringWidth(text, font, size)


def col_widths(header: list[str], rows: list[list[str]]) -> list[float]:
    """Widths that never split a word, then share the rest by content.

    2026-09-27 (user): with every cell printed in full the old
    proportional-to-length rule squeezed the narrow columns -- a tool name like
    ``CHEUI_m6A`` broke in two, and the header of the last column came out as
    "Coordina te harmo nisation".  Each column now first gets the width of its
    widest *word* (plus padding), and what is left is shared in proportion to the
    content, so a name is always readable and long prose still gets the room it
    needs.
    """
    n = len(header)
    if n == 1:
        return [AVAIL]
    pad = 9.0                                    #: LEFTPADDING + RIGHTPADDING + slack
    #: A token wider than this is a path, a URL or a long identifier, not a name:
    #: it is allowed to break, so its column does not reserve room for it.  A
    #: *name* (a tool, a mode, a level) never breaks.
    name_cap = 62.0
    floor = 40.0
    wants, musts = [], []
    for i in range(n):
        cells = [r[i] for r in rows if i < len(r)]
        head_word = _word_pt(header[i], "Helvetica-Bold", 8.0)
        cell_word = max((_word_pt(c, "Helvetica", 7.5) for c in cells), default=0.0)
        word = max(head_word, cell_word if cell_word + pad <= name_cap else 0.0) + pad
        musts.append(max(word, floor))
        widths = sorted(_cell_pt(c, "Helvetica", 7.5) for c in cells)
        wants.append(max(widths[int(0.8 * (len(widths) - 1))] + pad,
                         musts[-1]))
    rest = max(AVAIL - sum(musts), 0.0)
    weights = [max(w - m, 1.0) for w, m in zip(wants, musts)]
    total = sum(weights)
    return [m + rest * w / total for m, w in zip(musts, weights)]


def table_flowable(header: list[str], rows: list[list[str]]) -> Table:
    data = [[Paragraph(h, HEAD) for h in header]]
    for r in rows:
        data.append([Paragraph(c if c else "", CELL) for c in r])
    widths = col_widths(header, rows)
    t = Table(data, colWidths=widths, repeatRows=1, hAlign="CENTER")
    style = [
        ("GRID", (0, 0), (-1, -1), 0.25, colors.Color(0.75, 0.75, 0.75)),
        ("LINEABOVE", (0, 0), (-1, 0), 0.75, colors.black),
        ("LINEBELOW", (0, 0), (-1, 0), 0.5, colors.black),
        ("LINEBELOW", (0, -1), (-1, -1), 0.75, colors.black),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 2),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
    ]
    for i in range(1, len(data)):
        if i % 2 == 0:
            style.append(("BACKGROUND", (0, i), (-1, i), SHADE))
    t.setStyle(TableStyle(style))
    return t


def table_no(path: Path) -> int:
    """TableS10_... -> 10; plain alphabetical order puts S10/S11 before S1."""
    return int(path.name.split("_")[0][len("TableS"):])


def build() -> None:
    files = sorted(TAB.glob("TableS*_*.tsv"), key=table_no)
    assert files, "no TableS*.tsv found - run build_supp_tables.py first"
    out = HERE / "Supplementary_Tables.pdf"

    doc = BaseDocTemplate(str(out), pagesize=PAGE,
                          leftMargin=MARGIN, rightMargin=MARGIN,
                          topMargin=MARGIN, bottomMargin=MARGIN,
                          title="Supplementary Tables - RNA modification benchmark",
                          author="RNAModBench")
    frame = Frame(MARGIN, MARGIN, AVAIL, PAGE[1] - 2 * MARGIN, id="main")

    def footer(canv, d):
        canv.saveState()
        canv.setFont("Helvetica", 7.5)
        canv.drawCentredString(PAGE[0] / 2, 20, f"{canv.getPageNumber()}")
        canv.restoreState()

    doc.addPageTemplates([PageTemplate(id="all", frames=[frame], onPage=footer)])

    story = []
    for i, f in enumerate(files):
        title, header, rows, note = read_table(f)
        if i:
            story.append(Spacer(1, 6))
        story.append(Paragraph(title, TITLE))
        story.append(table_flowable(header, rows))
        if note:
            story.append(Paragraph("Note. " + note, NOTE))
        if i + 1 < len(files):
            from reportlab.platypus import PageBreak
            story.append(PageBreak())
    doc.build(story)
    print(f"wrote {out} ({out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    build()
