#!/usr/bin/env python
"""48 -- Supplementary Figure S10: independent validation of purified sites.

Reviewer 3 (R3-7) asked how the ``purified`` sites -- called in WT but absent
from the matched modification-deficient (KD/KO/IVT) call set -- behave when
they are judged *outside* the circular definition.  The answer is this new,
stand-alone Supplementary figure (it is NOT part of Figure S4, which keeps the
layout of the submitted ``sup4.pdf``):

    A  overlap of each site group with the GLORI reference (2 bp window)
    B  DRACH-motif fraction of each site group against the per-sample universe
       background (dashed)
    C  stoichiometry-semantic score (mod-ratio / probability tools only) of
       purified vs. shared sites, in two halves of three tools each

The three site groups are ``def-only`` (called in the deficient sample only),
``shared`` (called in both) and ``purified`` (called in WT only), always inside
the pair's common testable universe.  The panel letters restart at A because
this is a new figure number (S10, appended after the existing S1--S9).

Layout (2026-09-21, third rework -- one panel at a time, then assemble)
----------------------------------------------------------------------
Same architecture as Figure S4 (``common/panelpage.py``): every panel is drawn
on its own canvas at its final print size, passes
``pagelayout.assert_page_clean`` (text-text, text-in-foreign-panel,
**text-on-a-spine**, grid, minimum font size) on its own, and the A4 portrait
page is assembled from the panel PDFs with pypdf.  Sub-panels are 2.14 x 2.88
in, i.e. close to square instead of the letterbox strips of the landscape
draft; the 300 dpi PNG is rasterised from the assembled PDF.

Contract
--------
* claim: purified sites are systematically enriched for GLORI support, DRACH
  context and stoichiometry-semantic signal relative to shared / def-only
  sites, i.e. the group is not an artefact of the WT-minus-deficient
  definition (R3-7 (iii)); no accuracy claim is made from it;
* archetype: quantitative grid (three species columns; C split in two halves);
* per independent unit: every statistic is first computed per pair (Arabidopsis
  and HeLa: three replicates; mouse: the two cross-study datasets SRP166020 and
  SRP357195, never averaged), dots are pair values and the black bar is the mean
  across pairs;
* export: A4 portrait, vector PDF with embedded Arial (``pdf.fonttype`` 42) and
  a 300 dpi PNG, exact page size (no ``bbox_inches='tight'``); house rules:
  Arial, no grid, no annotation text inside a panel, bold panel letters in the
  page margin, ticks/axis labels >= 8.5 pt.

Reads (read-only) ``figures/figureS4/tables/``:
``figS4_validation_groups.tsv``, ``figS4_universe_drach.tsv`` -- both written
by ``41_figS4_tables.py``; nothing is recomputed here.

Outputs
-------
``figures/figureS5/figures/panels/FigureS10_{A,B,C}.pdf``
``figures/figureS5/figures/FigureS10_rev.{pdf,png}``
``figures/figureS5/logs/48_figS10_anchors.tsv``

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/src/harmonisation/scripts/48_figS10_validation.py
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
import sys
import time
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common import pagelayout, panelpage  # noqa: E402
from common.config import SITES_ROOT  # noqa: E402
from common.figstyle import apply  # noqa: E402

S4_DIR = (_RB / "figures/figureS4")
OUT_DIR = (_RB / "figures/figureS5")
FIG = OUT_DIR / "figures"
PANEL_DIR = FIG / "panels"
LOG = OUT_DIR / "logs"

#: the page of the submitted sup4.pdf (and of every SI figure page)
PAGE_IN = pagelayout.A4_PORTRAIT_IN

# --------------------------------------------------------------------------- #
# geometry, in inches (a panel canvas *is* its final print size)
# --------------------------------------------------------------------------- #
PANEL_W_IN = PAGE_IN[0]
LEFT_IN = 0.80
RIGHT_IN = 0.12
#: the column gap holds the neighbouring column's y tick labels (0.23 in at
#: 12 pt) plus a visible clearance from the previous column's frame line
COL_GAP_IN = 0.48
NCOL = 3
COL_W_IN = (PANEL_W_IN - LEFT_IN - RIGHT_IN - (NCOL - 1) * COL_GAP_IN) / NCOL

#: bands of a panel: 17 pt column titles, the drawing area (sub-panel
#: 2.13 x 2.84 in) and the x tick labels below it (no x axis title on this page)
PANEL_TITLE_IN = 0.36
PANEL_AXES_IN = 3.01
PANEL_BOTTOM_IN = 0.38
PANEL_H_IN = PANEL_TITLE_IN + PANEL_AXES_IN + PANEL_BOTTOM_IN

#: panel A carries the shared key on its own line inside the canvas (the
#: per-group panels keep no free band: their y range is max x 1.15); 2 rows of
#: 4 entries at the key font size need ~0.45 in

#: above A is gone; its height went into the three plot areas instead
#: (PANEL_AXES_IN 2.84 -> 3.01) and only the panel letter needs head room here
A_KEY_IN = 0.10
TOP_IN = 0.05
BOTTOM_IN = 0.11
GAP_IN = 0.08

#: site groups in plotting order
GROUPS = ["def_only", "shared", "purified"]
GROUP_LABEL = {"def_only": "def-only", "shared": "shared",
               "purified": "purified"}
GROUP_COLOR = {"def_only": "#9a9a9a", "shared": "#4b81b8",
               "purified": "#e8a76b"}

#: (column key, table species label, pairs) -- mouse stays split by study
COLUMNS: list[tuple[str, str, list[str]]] = [
    ("Arabidopsis", "Arabidopsis", ["Ath_rep1", "Ath_rep2", "Ath_rep3"]),
    ("Mouse", "Mouse", ["Mouse_studyA", "Mouse_studyB"]),
    ("HeLa", "Human", ["HeLa_rep1", "HeLa_rep2", "HeLa_rep3"]),
]
COLUMN_TITLE = {"Arabidopsis": "Arabidopsis", "Mouse": "Mouse", "HeLa": "Human"}
#: WT samples whose universe supplies the DRACH background of each column
BACKGROUND_SAMPLES = {
    "Arabidopsis": ["Arabidopsis_WT_rep1", "Arabidopsis_WT_rep2",
                    "Arabidopsis_WT_rep3"],
    "Mouse": ["mESCs_Mettl3_WT", "mES_WT"],
    "HeLa": ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"],
}

#: block above panel A.  Three of them only repeated the x tick labels
#: (def-only / shared / purified) and the error-bar rule is a caption sentence,
#: so panel A now carries the three marks of its own drawing and panel B names
#: its one reference line.
KEY_ENTRIES: list[tuple[str, dict]] = [
    ("def-only sites", dict(color=GROUP_COLOR["def_only"], lw=3.0)),
    ("shared sites", dict(color=GROUP_COLOR["shared"], lw=3.0)),
    ("purified sites", dict(color=GROUP_COLOR["purified"], lw=3.0)),
    ("Mouse A", dict(color="black", marker="o", ls="none", ms=5, mec="white")),
    ("Mouse B", dict(color="black", marker="o", ls="none", ms=5, mfc="white")),
    ("Mean of pairs", dict(color="black", lw=2.0)),
    ("DRACH background", dict(color="black", lw=1.4, ls=(0, (4, 2)))),
    ("Error bar: SD / IQR", dict(color="0.55", lw=3.0)),
]

#: tools whose native score is a modification ratio / probability
SEMANTIC_TOOLS = ["CHEUI_m6A", "m6Anet", "MINES", "Nanom6A", "DENA",
                  "NanoSPA_m6A"]
#: panel C is split in two halves so that neither becomes a letterbox strip
SCORE_HALVES = [SEMANTIC_TOOLS[:3], SEMANTIC_TOOLS[3:]]
DISPLAY = {"yanocomp": "Yanocomp"}

#: font sizes -- a scale of its own, one full notch *above* Figure S4

#: ones, so the same point size read smaller here; +25 % on every element makes
#: the type proportional to the drawing again.  Figure S4 keeps its own scale.
FS = {"tick": 12.0, "axis": 14.0, "subrow": 12.0, "column": 17.0,
      "letter": 24.0, "legend": 10.5}


# --------------------------------------------------------------------------- #
# canvas helpers
# --------------------------------------------------------------------------- #
def col_x_in(ci: int) -> float:
    """Left edge of species column ``ci`` inside a panel canvas (inches)."""
    return LEFT_IN + ci * (COL_W_IN + COL_GAP_IN)


def col_center_in(ci: int) -> float:
    return col_x_in(ci) + COL_W_IN / 2.0


def _rect(height_in: float, *, x_in: float, y_in: float, w_in: float,
          h_in: float) -> list[float]:
    return [x_in / PANEL_W_IN, y_in / height_in, w_in / PANEL_W_IN,
            h_in / height_in]


def _title(fig: plt.Figure, height_in: float, x_in: float, y_in: float,
           text: str, *, ha: str = "left", fontsize: float | None = None,
           weight: str = "normal") -> None:
    """Panel text placed in its own band (never overlapping an axis line)."""
    fig.text(x_in / PANEL_W_IN, y_in / height_in, text, ha=ha, va="bottom",
             fontsize=FS["subrow"] if fontsize is None else fontsize,
             fontweight=weight)


def panel_letter(fig: plt.Figure, letter: str, *, y_in: float) -> None:
    """Bold panel letter in the panel's own left margin."""
    pagelayout.margin_letter(fig, letter, y=y_in / fig.get_figheight(),
                             x=0.02 / PANEL_W_IN, fontsize=FS["letter"])


def pair_markers(ax: plt.Axes, x: float, values: list[float],
                 filled: list[bool], color: str) -> None:
    """One marker per independent pair, slightly offset around ``x``."""
    if not values:
        return
    offsets = np.linspace(-0.18, 0.18, len(values))
    for value, offset, solid in zip(values, offsets, filled):
        ax.plot([x + offset], [value], marker="o", ms=4.5,
                mfc=color if solid else "white",
                mec="white" if solid else color, mew=0.9, ls="none", zorder=5)


def group_axis(ax: plt.Axes, df: pd.DataFrame, value_col: str, key: str,
               pairs: list[str], ylabel: str, show_ylabel: bool,
               y_max: float) -> None:
    """Per-pair dots + black mean tick for the three site groups of a species.

    Same visual language as Figure S4 (dots + mean tick, no bars); the error bar
    is the SD across the independent pairs and the y ticks are labelled on every
    column.
    """
    for gi, group in enumerate(GROUPS):
        values, filled = [], []
        for pair in pairs:
            sub = df[(df["pair_id"] == pair) & (df["group"] == group)]
            if sub.empty:
                continue
            values.append(float(sub[value_col].mean()))
            filled.append(pair != "Mouse_studyB")
        if values:
            mean = float(np.mean(values))
            ax.plot([gi - 0.28, gi + 0.28], [mean, mean], color="black",
                    lw=2.0, zorder=4, solid_capstyle="butt")
            if len(values) > 1:
                sd = float(np.std(values, ddof=1))
                ax.errorbar([gi], [mean], yerr=[sd], fmt="none", ecolor="black",
                            elinewidth=1.4, capsize=3.5, capthick=1.4,
                            zorder=5)
        pair_markers(ax, gi, values, filled, GROUP_COLOR[group])
    ax.set_xticks(range(len(GROUPS)))
    ax.set_xticklabels([GROUP_LABEL[g] for g in GROUPS], fontsize=FS["tick"])
    ax.set_xlim(-0.55, len(GROUPS) - 0.45)
    ax.set_ylim(0, y_max)
    ax.tick_params(direction="out", length=3.5, width=1.0, labelsize=FS["tick"])
    if show_ylabel:
        ax.set_ylabel(ylabel, fontsize=FS["axis"], labelpad=3)


# --------------------------------------------------------------------------- #
# panels
# --------------------------------------------------------------------------- #
def render_panel_glori(dfv: pd.DataFrame, *, stem: str = "FigureS10_A",
                       gate: bool = True) -> panelpage.Panel:
    """(A) GLORI (2 bp) overlap of the three site groups per species."""
    y_max = float(dfv["glori_hit_rate_w2"].max()) * 1.15
    canvas_h = PANEL_H_IN + A_KEY_IN
    fig = panelpage.new_panel(PANEL_W_IN, canvas_h)
    for ci, (key, _species, pairs) in enumerate(COLUMNS):
        ax = fig.add_axes(_rect(canvas_h, x_in=col_x_in(ci),
                                y_in=PANEL_BOTTOM_IN, w_in=COL_W_IN,
                                h_in=PANEL_AXES_IN))
        group_axis(ax, dfv, "glori_hit_rate_w2", key, pairs,
                   "GLORI overlap (2 bp)", show_ylabel=(ci == 0), y_max=y_max)
        _title(fig, canvas_h, col_center_in(ci),
               PANEL_BOTTOM_IN + PANEL_AXES_IN + 0.02, COLUMN_TITLE[key],
               ha="center", fontsize=FS["column"], weight="bold")
        
        #: the group mean and the error bar in every column, plus the Mouse A /
        #: Mouse B coding in the column that uses it; the top key band is gone
        wanted = ["Mean of pairs", "Error bar: SD"]
        if key == "Mouse":
            wanted = ["Mouse A", "Mouse B"] + wanted
        handles = [Line2D([], [], label=label, **kw)
                   for label, kw in KEY_ENTRIES if label in wanted]
        handles.append(Line2D([], [], color="0.55", lw=3.0,
                              label="Error bar: SD"))
        ax.legend(handles=handles, loc="upper left", frameon=False,
                  fontsize=FS["legend"], handlelength=1.5, handletextpad=0.45,
                  labelspacing=0.28, borderpad=0.0, borderaxespad=0.35)
    panel_letter(fig, "A", y_in=canvas_h - 0.04)
    return panelpage.save_panel(fig, PANEL_DIR / f"{stem}.pdf", gate=gate)


def render_panel_drach(dfv: pd.DataFrame, df_bg: pd.DataFrame, *,
                       stem: str = "FigureS10_B",
                       gate: bool = True) -> panelpage.Panel:
    """(B) DRACH fraction of the three site groups + universe background."""
    y_max = float(dfv["drach_rate"].max()) * 1.15
    fig = panelpage.new_panel(PANEL_W_IN, PANEL_H_IN)
    for ci, (key, _species, pairs) in enumerate(COLUMNS):
        ax = fig.add_axes(_rect(PANEL_H_IN, x_in=col_x_in(ci),
                                y_in=PANEL_BOTTOM_IN, w_in=COL_W_IN,
                                h_in=PANEL_AXES_IN))
        group_axis(ax, dfv, "drach_rate", key, pairs, "DRACH fraction",
                   show_ylabel=(ci == 0), y_max=y_max)
        for sample in BACKGROUND_SAMPLES[key]:
            hit = df_bg[df_bg["sample"] == sample]
            if hit.empty:
                continue
            rate = float(hit["drach_rate"].iloc[0])
            ax.plot([-0.45, len(GROUPS) - 0.55], [rate, rate], color="black",
                    lw=1.4, ls=(0, (4, 2)), zorder=4)
            
            #: line (and its error bar) inside the axes -- no column is left
            #: without a key, and nothing sits above the plot area
            ax.legend(handles=[Line2D([], [], color="black", lw=1.4,
                                      ls=(0, (4, 2)),
                                      label="DRACH background"),
                               Line2D([], [], color="0.55", lw=3.0,
                                      label="Error bar: SD")],
                      loc="upper right", frameon=False,
                      fontsize=FS["legend"], handlelength=1.5,
                      handletextpad=0.45, labelspacing=0.28, borderpad=0.0,
                      borderaxespad=0.4)
        _title(fig, PANEL_H_IN, col_center_in(ci),
               PANEL_BOTTOM_IN + PANEL_AXES_IN + 0.02, COLUMN_TITLE[key],
               ha="center", fontsize=FS["column"], weight="bold")
    panel_letter(fig, "B", y_in=PANEL_H_IN - 0.04)
    return panelpage.save_panel(fig, PANEL_DIR / f"{stem}.pdf", gate=gate)


def render_panel_score(dfv: pd.DataFrame, *, stem: str = "FigureS10_C",
                       gate: bool = True) -> panelpage.Panel:
    """(C) Stoichiometry-semantic score of purified vs. shared sites.

    Two halves of three tools each, so neither half becomes a letterbox strip
    (the draft put all six tools in one wide panel).  Dots are the mean of the
    per-pair score medians and the error bar the mean IQR of those per-pair
    distributions.
    """
    fig = panelpage.new_panel(PANEL_W_IN, PANEL_H_IN)
    gap_in = 0.60
    half_w = (PANEL_W_IN - LEFT_IN - RIGHT_IN - gap_in) / 2.0
    y_max = 1.05
    for hi, tools in enumerate(SCORE_HALVES):
        x_in = LEFT_IN + hi * (half_w + gap_in)
        ax = fig.add_axes(_rect(PANEL_H_IN, x_in=x_in, y_in=PANEL_BOTTOM_IN,
                                w_in=half_w, h_in=PANEL_AXES_IN))
        for ti, tool in enumerate(tools):
            for group, offset in (("shared", -0.16), ("purified", 0.16)):
                sub = dfv[(dfv["tool"] == tool) & (dfv["group"] == group)]
                median = sub["score_median"].dropna()
                iqr = (sub["score_q75"] - sub["score_q25"]).dropna()
                if median.empty:
                    continue
                half = float(iqr.mean()) / 2.0 if not iqr.empty else 0.0
                ax.errorbar([ti + offset], [float(median.mean())],
                            yerr=[[half], [half]], fmt="o", ms=7,
                            color=GROUP_COLOR[group],
                            ecolor=GROUP_COLOR[group], elinewidth=2.6,
                            capsize=4.0, capthick=2.4, zorder=3)
        ax.set_xlim(-0.6, len(tools) - 0.4)
        ax.set_ylim(0, y_max)
        ax.set_xticks(range(len(tools)))
        ax.set_xticklabels([DISPLAY.get(t, t) for t in tools],
                           fontsize=FS["tick"])
        ax.tick_params(direction="out", length=3.5, width=1.0,
                       labelsize=FS["tick"])
        if hi == 0:
            ax.set_ylabel("Score (mod. ratio / probability)", fontsize=FS["axis"],
                          labelpad=3)
        
        #: colour pair and the IQR error bar are now named inside each half
        ax.legend(handles=[
            Line2D([], [], color=GROUP_COLOR["shared"], marker="o", ms=6,
                   ls="none", label="shared"),
            Line2D([], [], color=GROUP_COLOR["purified"], marker="o", ms=6,
                   ls="none", label="purified"),
            Line2D([], [], color="0.55", lw=3.0, label="Error bar: IQR"),
        ], loc="upper left", frameon=False, fontsize=FS["legend"],
            handlelength=1.5, handletextpad=0.45, labelspacing=0.28,
            borderpad=0.0, borderaxespad=0.35)
    panel_letter(fig, "C", y_in=PANEL_H_IN - 0.04)
    return panelpage.save_panel(fig, PANEL_DIR / f"{stem}.pdf", gate=gate)


# --------------------------------------------------------------------------- #
# legend stripe
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# page assembly
# --------------------------------------------------------------------------- #
def page_placements(panels: dict[str, panelpage.Panel],
                    ) -> list[tuple[panelpage.Panel, float, float]]:
    """Stack the pieces top-down; returns ``[(panel, x_in, y_in)]``.

    Every key of this page sits inside the A canvas, so there is no legend
    stripe to place and the freed height went into the panels.  The last panel
    carries no trailing gap (``("C", 0.0)``), otherwise the budget would be
    charged 0.08 in for a gap that is never drawn.
    """
    order = [("A", GAP_IN), ("B", GAP_IN), ("C", 0.0)]
    cursor = PAGE_IN[1] - TOP_IN
    out: list[tuple[panelpage.Panel, float, float]] = []
    for key, gap in order:
        panel = panels[key]
        bottom = cursor - panel.height_in
        if bottom < -1e-9:
            raise SystemExit(f"{key}: page budget overflows by {-bottom:.3f} in")
        out.append((panel, 0.0, bottom))
        cursor = bottom - gap
    if cursor < BOTTOM_IN - 1e-9:
        raise SystemExit(f"bottom margin {cursor:.3f} in < {BOTTOM_IN:.2f} in")
    return out


# --------------------------------------------------------------------------- #
# numeric anchors
# --------------------------------------------------------------------------- #
def write_anchors(dfv: pd.DataFrame, df_bg: pd.DataFrame) -> pd.DataFrame:
    """Log the numbers the figure legend quotes (per group and per species)."""
    rows: list[dict[str, object]] = []
    for group in GROUPS:
        sub = dfv[dfv["group"] == group]
        rows.append({"anchor": f"glori_overlap_{group}",
                     "value": sub["glori_hit_rate_w2"].mean(),
                     "definition": "mean over the 104 tool-pair combinations "
                                   "(GLORI, 2 bp)"})
        rows.append({"anchor": f"drach_{group}",
                     "value": sub["drach_rate"].mean(),
                     "definition": "mean over the 104 tool-pair combinations"})
    for key, _species, pairs in COLUMNS:
        sub = dfv[dfv["pair_id"].isin(pairs)]
        for group in GROUPS:
            g = sub[sub["group"] == group]
            rows.append({
                "anchor": f"glori_overlap_{group}_{key}",
                "value": g["glori_hit_rate_w2"].mean(),
                "definition": "mean over the column's tools and pairs"})
        for sample in BACKGROUND_SAMPLES[key]:
            hit = df_bg[df_bg["sample"] == sample]
            if not hit.empty:
                rows.append({
                    "anchor": f"drach_background_{sample}",
                    "value": float(hit["drach_rate"].iloc[0]),
                    "definition": "DRACH fraction of the sample's universe"})
    for tool in SEMANTIC_TOOLS:
        sub = dfv[dfv["tool"] == tool]
        for group in ("shared", "purified"):
            rows.append({
                "anchor": f"score_{group}_{tool}",
                "value": sub[sub["group"] == group]["score_median"].mean(),
                "definition": "mean over pairs of the per-pair score median"})
    out = pd.DataFrame(rows)
    LOG.mkdir(parents=True, exist_ok=True)
    out.to_csv(LOG / "48_figS10_anchors.tsv", sep="\t", index=False)
    return out


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main() -> None:
    t0 = time.time()
    apply()
    plt.rcParams.update({"axes.linewidth": 1.2})
    FIG.mkdir(parents=True, exist_ok=True)

    dfv = pd.read_csv(S4_DIR / "tables" / "figS4_validation_groups.tsv",
                      sep="\t")
    df_bg = pd.read_csv(S4_DIR / "tables" / "figS4_universe_drach.tsv",
                        sep="\t")
    anchors = write_anchors(dfv, df_bg)

    panels = {
        "A": render_panel_glori(dfv),
        "B": render_panel_drach(dfv, df_bg),
        "C": render_panel_score(dfv),
    }
    out_pdf = FIG / "FigureS10_rev.pdf"
    page = panelpage.compose_page(PAGE_IN, page_placements(panels), out_pdf)
    panelpage.pdf_to_png(page, FIG / "FigureS10_rev.png")

    print(f"saved -> {page} + {FIG / 'FigureS10_rev.png'} "
          f"({time.time() - t0:.1f} s)")
    print(f"anchors -> {LOG / '48_figS10_anchors.tsv'} "
          f"({datetime.now():%Y-%m-%d %H:%M})")
    with pd.option_context("display.width", 120, "display.max_rows", 80):
        print(anchors.to_string(index=False))


if __name__ == "__main__":
    main()
