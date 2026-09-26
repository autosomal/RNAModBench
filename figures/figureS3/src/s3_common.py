#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared constants and helpers for the rebuilt Figure S3.

Everything on the page is drawn in *page points* (1/72 inch).  The page is
sized exactly like the submitted ``sup3.pdf`` (595.276 x 633.598 pt) so the
rebuilt page can replace it one-to-one.

Style constants were sampled from a 150-dpi raster of the submitted page
(``$RNAMODBENCH_LOCAL/submission/$RNAMODBENCH_LOCAL/figures_original_and_build_inputs/sup/sup3.pdf``):

* venn fills   -- orange #F1BA8A, steel blue #7D9EBB, lens sage #D7D7C2
* text ink     -- #231F20 (no pure black)
* hit-rate row -- species colours #1E888B / #F5B264 / #3778A0, markers o/s/^
* ratio row    -- m6Anet #58BFC0, Nanom6A #F5B264, DENA #A0D9F3, MINES #9FB7C5
  (the legacy notebook assigned these randomly per run; they are now frozen).

Venn geometry follows the submitted design: both circles have the same
radius and the centre separation is chosen so that the lens area equals the
replicates' Jaccard index (verified against the submitted panels: measured
d/r 0.553 -> J 0.484 vs. actual 0.498 for Arabidopsis; 0.29 -> 0.689 vs.
0.689 for Human).
"""

from __future__ import annotations

# --- RNAModBench path bootstrap (added when this file was deposited) ----------
import os as _rb_os, pathlib as _rb_pl


def _rb_find(start):
    for p in (start, *start.parents):
        if (p / "RNAMOD_BENCH_ROOT").exists():
            return p
    return start


_RB = _rb_pl.Path(_rb_os.environ.get("RNAMODBENCH_ROOT") or _rb_find(_rb_pl.Path(__file__).resolve().parent))
_XB = _rb_pl.Path(_rb_os.environ.get("RNAMODBENCH_LOCAL") or (_RB / "_local"))
# --------------------------------------------------------------------------- #
import math
from pathlib import Path

import matplotlib.pyplot as plt

# --------------------------------------------------------------------------- #
# paths
# --------------------------------------------------------------------------- #
PROJECT = Path(str(_RB))
GLORI_DIR = (_XB / "third_party/NGS/GLORI")
M6ASEQ_DIR = (_XB / "third_party/NGS/m6A-Seq")
LEGACY = (_XB / "archive/output_legacy_20260916")

#: moved by the user on 2026-09-19 (was $RNAMODBENCH_LOCAL/figures_original/NGS_visualization/S3_revision);
#: renamed S3_revision -> figures/figureS3 on 2026-09-21 (fig*_revision house style)
OUT_DIR = (_RB / "figures/figureS3")
PANEL_DIR = (_RB / "figures/figureS3/panels")
TABLE_DIR = (_RB / "figures/figureS3/tables")

#: panel-A replicate labels: unified "species + repN" style (user 2026-09-19)
REPLICATE_LABELS = {
    "Arabidopsis": ("Arabidopsis-rep1", "Arabidopsis-rep2"),
    "Mouse": ("Mouse-rep1", "Mouse-rep2"),
    "Human": ("HeLa-rep1", "HeLa-rep2"),
}

#: retired "GLORI hit rate" / "Hit Rate" -> house label (figures/figure5 README)
YLABEL_PPV = "PPV vs. GLORI (2 bp)"

SUP3_PDF = (_XB / "submission/manuscript/sup/sup3.pdf")

# --------------------------------------------------------------------------- #
# revision evaluation layer (harmonisation) -- the source of panels B and C
# --------------------------------------------------------------------------- #
SITES_V2 = (_RB / "data")
EVAL_TABLES = (_RB / "data/evaluation/tables")
SITES_CLEAN = (_RB / "data/callsets")
UNIVERSE_ROOT = (_XB / "harmonisation/universe")
CONFUSION_TSV = (_RB / "data/evaluation/tables/m6a_glori_confusion.tsv")

#: manuscript tool scope (13 m6A tools actually used in the paper)
M6A_TOOLS = (
    "CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
    "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
    "NanoSPA_m6A", "xPore", "yanocomp",
)


PPV_WINDOW = 2

#: independent units per species group (mouse = two independent studies)
UNITS_BY_SPECIES = {
    "Arabidopsis": ("Arabidopsis_WT", ["Arabidopsis_WT_rep1",
                                       "Arabidopsis_WT_rep2",
                                       "Arabidopsis_WT_rep3"]),
    "Mouse": ("Mouse_WT", ["mES_WT", "mESCs_Mettl3_WT"]),
    "Human": ("HeLa_WT", ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"]),
}

#: mouse study provenance (never averaged together -- user rule)
MOUSE_STUDY = {"mES_WT": "SRP357195", "mESCs_Mettl3_WT": "SRP166020"}

#: in-figure names of the two independent mouse WT mESC samples (user 2026-09-21:
#: the raw sample ids invited "different cell line / mutant" doubts -- both are WT).
#: study A = mESCs_Mettl3_WT (SRP166020); study B = mES_WT (SRP357195) -- never swap.
MOUSE_STUDY_LABEL = {"mESCs_Mettl3_WT": "mouse study A", "mES_WT": "mouse study B"}

GLORI_BEDS = {
    "Arabidopsis": (_XB / "third_party/NGS/GLORI/Arabidopsis_GLORI.bed"),
    "Mouse": (_XB / "third_party/NGS/GLORI/Mouse_GLORI_liftover.bed"),
    "Human": (_XB / "third_party/NGS/GLORI/Hela_GLORI.bed"),  # rebuilt 2026-09-19 (both > 0.1)
}

#: (chrom, start, ratio) allAs BED per replicate -- criterion: ratio > 0.1.
ALLAS_BEDS = {
    "Arabidopsis": (
        (_XB / "third_party/NGS/GLORI/GSM7873784_Arabidopsis_allAs_cov15_rep1.bed"),
        (_XB / "third_party/NGS/GLORI/GSM7873785_Arabidopsis_allAs_cov15_rep2.bed"),
    ),
    "Mouse": (
        (_XB / "third_party/NGS/GLORI/GSM7873782_mESC_allAs_cov15_rep1.bed"),
        (_XB / "third_party/NGS/GLORI/GSM7873783_mESC_allAs_cov15_rep2.bed"),
    ),
}

#: HeLa GLORI FDR txt per replicate (criterion unified to ratio > 0.1 each).
HELA_FDR = (
    (_XB / "third_party/NGS/m6A-Seq/GSM6432595_Hela-1_35bp_m2.totalm6A.FDR.txt"),
    (_XB / "third_party/NGS/m6A-Seq/GSM6432596_Hela-2_35bp_m2.totalm6A.FDR.txt"),
)

LEGACY_TOOLS_TXT = {
    "Arabidopsis": (_XB / "archive/output_legacy_20260916/Arabidopsis_WT/Tools.txt"),
    "Mouse": (_XB / "archive/output_legacy_20260916/Mouse_WT/Tools.txt"),
    "Human": (_XB / "archive/output_legacy_20260916/HeLa_WT/aggregated_data/m6A/Tools.txt"),
}
LEGACY_MODRATIO_DIR = {
    "Arabidopsis": (_XB / "archive/output_legacy_20260916/Arabidopsis_WT/mod_ratio"),
    "Mouse": (_XB / "archive/output_legacy_20260916/Mouse_WT/mod_ratio"),
    "Human": (_XB / "archive/output_legacy_20260916/HeLa_WT/tools/mod_ratio"),
}

# --------------------------------------------------------------------------- #
# page geometry (pt, origin = top-left of the page, y grows downward)
# --------------------------------------------------------------------------- #
PAGE_W = 595.276
PAGE_H = 633.598

#: panel letters (x = left edge, y = centre; measured on the submitted page)
LETTER_POS = {"A": (11.0, 16.5), "B": (11.0, 196.5), "C": (11.0, 411.0)}
LETTER_SIZE = 12.0

#: row A -- venns
A_TITLE_Y = 25.9          # species title centre
A_CIRCLE_Y = 100.1        # circle centre
A_R = 59.0                # circle radius (uniform)
A_LABEL_Y = 169.4         # replicate-label centre
A_PANEL_CX = (100.25, 303.85, 499.65)   # panel centre x (from the submitted page)

#: row B -- hit rate
B_AXES_X0 = 48.5          # left edge of the first axes
B_AXES_W = 143.5          # axes width
B_AXES_GAP = 45.0         # gap between axes
B_AXES_TOP = 212.2
B_AXES_BOTTOM = 368.6

#: row C -- modification-ratio agreement
C_AXES_TOP = 444.9
C_AXES_BOTTOM = 596.1

# --------------------------------------------------------------------------- #
# style
# --------------------------------------------------------------------------- #
ORANGE = "#F1BA8A"
BLUE = "#7D9EBB"
SAGE = "#D7D7C2"
INK = "#231F20"

SPECIES_ORDER = ("Arabidopsis", "Mouse", "Human")
SPECIES_COLORS = {"Arabidopsis": "#1E888B", "Mouse": "#F5B264", "Human": "#3778A0"}
SPECIES_MARKERS = {"Arabidopsis": "o", "Mouse": "s", "Human": "^"}

#: legacy tool -> colour mapping frozen from the submitted panel C
TOOL_COLORS = {
    "m6Anet": "#58BFC0",
    "Nanom6A": "#F5B264",
    "DENA": "#A0D9F3",
    "MINES": "#9FB7C5",
}

#: font sizes (pt) measured on the submitted page
F_PANEL_LETTER = 12.0
F_TITLE = 10.5          # species titles (rows A/B/C)
F_COUNT = 9.5           # venn counts
F_REPLABEL = 9.0        # replicate labels under the venns
F_AXIS_LABEL = 9.5      # "Hit Rate" / "GLORI Modification Ratio" (bold)
F_TICK = 7.0
F_LEGEND_B = 6.8
F_LEGEND_C = 6.3

LW_MAIN = 1.2           # hit-rate line
LW_REG = 1.4            # regression lines (panel C)
LW_SPINE = 0.9
LW_AVG = 0.9            # dashed average / identity lines


def apply_page_style() -> None:
    """rcParams shared by every S3 script (Arial, no grid, editable PDF text)."""
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.0,
            "axes.linewidth": LW_SPINE,
            "axes.grid": False,
            "axes.edgecolor": "black",
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "figure.dpi": 300,
            "savefig.dpi": 300,
        }
    )


# --------------------------------------------------------------------------- #
# replicate overlap under the unified criterion
# --------------------------------------------------------------------------- #
def _parse_allAs_bed(path: Path) -> dict[tuple[str, int], float]:
    """(chrom, start) -> modification ratio (column 5, index 4) of an allAs BED."""
    data: dict[tuple[str, int], float] = {}
    with path.open() as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            try:
                data[(p[0], int(p[1]))] = float(p[4])
            except ValueError:
                continue
    return data


def _parse_hela_fdr(path: Path) -> dict[tuple[str, int], float]:
    """HeLa FDR txt -> {(chrom, start=Sites-1): NormeRatio} (column 11, index 10)."""
    data: dict[tuple[str, int], float] = {}
    with path.open(encoding="utf-8") as fh:
        header = fh.readline()
        if not header.startswith("Chr"):
            raise ValueError(f"{path}: missing 'Chr' header")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 11:
                continue
            try:
                chrom = p[0][3:] if p[0].startswith("chr") else p[0]
                data[(chrom, int(p[1]) - 1)] = float(p[10])
            except ValueError:
                continue
    return data


def compute_overlap() -> dict[str, dict[str, float]]:
    """Replicate overlap per species under the unified criterion.

    Criterion (user decision 2026-09-19, applied to all three species):
    a replicate's site set = positions with modification ratio > 0.1 in that
    replicate; the reference intersection = positions present in both.
    """
    out: dict[str, dict[str, float]] = {}
    for species in SPECIES_ORDER:
        if species == "Human":
            d1 = _parse_hela_fdr(HELA_FDR[0])
            d2 = _parse_hela_fdr(HELA_FDR[1])
        else:
            d1 = _parse_allAs_bed(ALLAS_BEDS[species][0])
            d2 = _parse_allAs_bed(ALLAS_BEDS[species][1])
        s1 = {k for k, v in d1.items() if v > 0.1}
        s2 = {k for k, v in d2.items() if v > 0.1}
        inter = s1 & s2
        raw_inter = set(d1) & set(d2)
        mean01 = {k for k in raw_inter if (d1[k] + d2[k]) / 2.0 > 0.1}
        out[species] = {
            "rep1_sites": len(s1),
            "rep2_sites": len(s2),
            "only1": len(s1 - s2),
            "only2": len(s2 - s1),
            "intersection": len(inter),
            "union": len(s1 | s2),
            "jaccard": len(inter) / len(s1 | s2),
            "raw_intersection": len(raw_inter),   # no ratio filter (legacy figure)
            "mean_gt01": len(mean01),             # legacy HeLa file criterion
        }
    return out


#: counts independently verified on 2026-09-19 (awk / python recomputation)
EXPECTED_OVERLAP = {
    "Arabidopsis": {"only1": 40_461, "only2": 41_894, "intersection": 80_624},
    "Mouse": {"only1": 48_221, "only2": 59_580, "intersection": 41_961},
    "Human": {"only1": 23_686, "only2": 26_606, "intersection": 112_451},
}
EXPECTED_HUMAN_TRACE = {"raw_intersection": 113_606, "mean_gt01": 113_485}


def check_overlap(overlaps: dict[str, dict[str, float]]) -> None:
    """Hard assert: the unified-criterion counts must reproduce the verified numbers."""
    for species, exp in EXPECTED_OVERLAP.items():
        got = overlaps[species]
        for field, want in exp.items():
            if got[field] != want:
                raise AssertionError(
                    f"{species}.{field} = {got[field]}, expected {want}"
                )
    human = overlaps["Human"]
    for field, want in EXPECTED_HUMAN_TRACE.items():
        if human[field] != want:
            raise AssertionError(f"Human.{field} = {human[field]}, expected {want}")
    print("[overlap] all counts match the verified values "
          "(Arab 80,624 / Mouse 41,961 / Human 112,451; traces 113,606 / 113,485)")


# --------------------------------------------------------------------------- #
# venn geometry / drawing
# --------------------------------------------------------------------------- #
def jaccard_to_separation(jaccard: float) -> float:
    """Centre separation d (in units of the radius r) for a lens area = Jaccard.

    Two equal circles, separation d: lens/circle area = (2/pi) * (acos(t) -
    t*sqrt(1-t^2)) with t = d/(2r); the drawn Jaccard is lens/(2-lens).
    """
    from scipy.optimize import brentq

    s = 2.0 * jaccard / (1.0 + jaccard)          # lens area / circle area
    f = lambda t: (2.0 / math.pi) * (math.acos(t) - t * math.sqrt(1.0 - t * t)) - s
    return 2.0 * brentq(f, 1e-9, 1.0 - 1e-9)


def draw_venn_pair(
    ax: plt.Axes,
    cx: float,
    n1: int,
    n_inter: int,
    n2: int,
    jaccard: float,
    label1: str,
    label2: str,
    title: str,
) -> None:
    """Draw one species' two-circle venn, centred horizontally at ``cx``.

    Flat fills (no alpha): orange circle, blue circle, then a sage copy of the
    blue circle clipped to the orange one (exact lens).  The overlap encodes
    the replicates' Jaccard index.  Counts sit in the three regions; the
    species title is centred above, the replicate labels below.
    """
    from matplotlib.patches import Circle

    r = A_R
    sep = jaccard_to_separation(jaccard) * r
    c1 = cx - sep / 2.0
    c2 = cx + sep / 2.0

    circ1 = Circle((c1, A_CIRCLE_Y), r, facecolor=ORANGE, edgecolor="none", zorder=1)
    circ2 = Circle((c2, A_CIRCLE_Y), r, facecolor=BLUE, edgecolor="none", zorder=2)
    ax.add_patch(circ1)
    ax.add_patch(circ2)
    lens = Circle((c2, A_CIRCLE_Y), r, facecolor=SAGE, edgecolor="none", zorder=3)
    ax.add_patch(lens)
    lens.set_clip_path(circ1)

    ax.text(cx, A_TITLE_Y, title, ha="center", va="center",
            fontsize=F_TITLE, fontweight="bold", color=INK)

    for x, n in ((c1 - r + sep / 2.0, n1),
                 (c1 + sep / 2.0, n_inter),
                 (c1 + r + sep / 2.0, n2)):
        ax.text(x, A_CIRCLE_Y, f"{n:,}".replace(",", ""),
                ha="center", va="center", fontsize=F_COUNT,
                fontweight="bold", color=INK)

    # replicate labels: centred on their circle, pushed apart if they collide
    t1 = ax.text(c1, A_LABEL_Y, label1, ha="center", va="center",
                 fontsize=F_REPLABEL, fontweight="bold", color=INK)
    t2 = ax.text(c2, A_LABEL_Y, label2, ha="center", va="center",
                 fontsize=F_REPLABEL, fontweight="bold", color=INK)
    fig = ax.figure
    fig.canvas.draw()  # expire the renderer so window extents are valid
    b1, b2 = t1.get_window_extent(), t2.get_window_extent()
    inv = ax.transData.inverted()
    x1a, _ = inv.transform((b1.x0, b1.y0))
    x1b, _ = inv.transform((b1.x1, b1.y1))
    x2a, _ = inv.transform((b2.x0, b2.y0))
    x2b, _ = inv.transform((b2.x1, b2.y1))
    gap_pt = 4.0
    if x1b > x2a - gap_pt:  # collide -> push outward, keep the pair centred
        w1, w2 = x1b - x1a, x2b - x2a
        total = w1 + w2 + gap_pt
        left = cx - total / 2.0
        t1.set_position((left + w1 / 2.0, A_LABEL_Y))
        t2.set_position((left + w1 + gap_pt + w2 / 2.0, A_LABEL_Y))


def label_axes_page(ax: plt.Axes, y_top: float = 0.0, height: float = PAGE_H) -> None:
    """Configure a full-page / full-width axes in top-down point coordinates."""
    ax.set_xlim(0.0, PAGE_W)
    ax.set_ylim(y_top + height, y_top)
    ax.set_aspect("equal")
    ax.axis("off")


def panel_letter(fig: plt.Figure, letter: str, x: float, y: float) -> None:
    fig.text(x / PAGE_W, 1.0 - y / PAGE_H, letter, fontsize=F_PANEL_LETTER,
             fontweight="bold", color=INK, ha="left", va="center")
