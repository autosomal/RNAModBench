#!/usr/bin/env python
"""Figure 5, assembled: the four revision rows on one page-width canvas.

Panels (rows), all read from the ``figures/figure5/tables`` summary tables -- no
metric is recomputed here, so the assembled figure cannot drift from the
per-panel figures:

  A  mean rank (per independent sequencing unit) of MCC / F1 / Precision /
     Recall (column order 2026-09-20, user rev.2: MCC then F1 first), 13 m6A tools
     x 3 species;      each species panel is ordered by the mean **MCC** rank,
     best first (user decision 2026-09-20, rev.2; rev.1 was mean F1 rank;
     originally mean Precision (PPV) rank per 2026-09-19).  Mouse = the
     single study ``mES_WT``;
  B  tool-reported modification ratio (10 bins) vs ``PPV vs. GLORI (2 bp)``,
     MINES / m6Anet / DENA / Nanom6A;
  C  matching window (0/1/2/5/10/20/50 bp) vs ``PPV vs. GLORI (2 bp)``;
  D  matching window vs the exact single-nucleotide localisation among the
     window matches (R1-3 / R3-3).

House rules applied here (2026-09-19): drawn at the **final printed size**
(0.95 * \\textwidth = 6.66 in of the Wiley USG layout, 178 mm text width), no
rescaling afterwards; every text element is asserted to be >= 7.0 pt at that
size; no gridlines; no small in-panel annotation text; Arial; vector PDF + 300
dpi PNG.

Outputs -> figures/figure5/figures/Figure5_rev.{pdf,png}

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/figures/figure5/src/39_fig5_assembled.py
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
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
import numpy as np
import pandas as pd
import seaborn as sns

HERE = _RB / "src/harmonisation"
sys.path.insert(0, str(HERE))

from common import config as C                                        # noqa: E402
from common.figstyle import apply as apply_style                     # noqa: E402
from common.manifest import setup_logger                             # noqa: E402

TAB = (_RB / "figures/figure5/tables")
FIG = (_RB / "figures/figure5/figures")

# ---------------------------------------------------------------- canvas ----
#: 0.95 * \textwidth of USG.cls (paper 210 mm, left/right 16 mm -> 178 mm)
CANVAS_W = 6.66
#: \textheight = 276 - 12 - 30 mm = 234 mm = 9.21 in; a full-page figure* must
#: still fit the caption, so the canvas stays below this
MAX_H = 8.95            # 2026-09-27: raised from 8.50 so the C row can keep its
# own 13-tool key in the service strip under the row (\textheight = 9.21 in)
MIN_PT = 7.0

TOP = 0.02
ROW_TITLE = 0.13           # space reserved above each row of axes for its title
#: the gap under a row must host that row's x tick labels and axis label:
#: A = four metric names rotated 45 deg (0.37 in + pad), B = ten bin labels
#: rotated 45 deg (0.30 in) plus the axis label, C/D = short window ticks plus
#: the axis label.  Nothing is ever shared with the next row's title band.
#: 2026-09-20: widened by ~0.10 in per gap (user request) so the four panel
#: rows read as clearly separated facets; canvas total 8.14 -> 8.45 in,
#: still below the 8.50 in full-page cap.

#: 13-tool key fits under the C row without touching its x axis label.
ROW_GAP = {"A": 0.55, "B": 0.58, "C": 0.72}
ROWS_H = {"A": 1.88, "B": 1.28, "C": 1.28, "D": 1.28}
LEGEND_H = 0.66            # hosts the D row's axis label + the 13-tool legend

#: geometry of the A row.  Its three panels each carry their own tool-name
#: labels, and the widest label ("NanoSPA_m6A" = 0.713 in at 7.5 pt Arial) must
#: fit inside the reserved gap.  Shrinking the font does NOT solve that (0.714 in
#: at 7.0 pt, and 6.5 pt would break the >= 7 pt floor), so the gap is what
#: grows -- asserted at run time by ``assert_label_fit``.
A_LEFT, A_W, A_LABEL_GAP, A_GAP = 0.82, 1.24, 0.80, 0.06
#: geometry of the B/C/D rows (tool labels only in the shared legends).

#: read as separated panels; each panel is 1.86 in wide, and the B-row 45 deg
#: bin-label clearance assertion still passes (0.132 in > 0.109 in needed).
L_LEFT, L_GAP, L_RIGHT = 0.42, 0.30, 0.06

SPECIES = ["Arabidopsis", "Mouse", "Human"]
#: species names only: the replicate counts (Arabidopsis n = 3, mouse n = 1
#: study, human n = 3) and the HeLa provenance belong to the caption
TITLE = {"Arabidopsis": "Arabidopsis", "Mouse": "Mouse", "Human": "Human"}

#: column order of the 5A heatmap; 2026-09-20 (user, rev.2): MCC / F1 first
METRICS = ["mcc", "f1", "precision", "recall"]
METRIC_LABEL = {"precision": "Precision", "recall": "Recall",
                "f1": "F1", "mcc": "MCC"}
#: row-ordering key of the 5A heatmap; 2026-09-20 (user, rev.2): MCC rank, best
#: first (rev.1 was the mean F1 rank; originally the mean Precision (PPV)
#: rank per user decision 2026-09-19) -- MCC row order matches the
#: "m6Anet highest MCC across species" narrative
SORT_METRIC = "mcc"

TOOL_ORDER = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
              "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
              "NanoSPA_m6A", "xPore", "yanocomp"]
DISPLAY = {"yanocomp": "Yanocomp"}
TOOL4 = ["MINES", "m6Anet", "DENA", "Nanom6A"]
TOOL4_COLOR = {"MINES": "#4C72B0", "m6Anet": "#DD8452",
               "DENA": "#55A868", "Nanom6A": "#C44E52"}
_pal = plt.get_cmap("tab20").colors
TOOL_COLOR = {t: _pal[i % 20] for i, t in enumerate(TOOL_ORDER)}

BINS = [f"{a:.1f}-{b:.1f}" for a, b in
        zip([0, .1, .2, .3, .4, .5, .6, .7, .8, .9],
            [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0])]
WINDOWS = C.WINDOWS

FS = {"tick": 7.5, "label": 8.0, "title": 8.5, "annot": 7.5,
      "legend": 7.2, "letter": 9.5, "cbar": 7.5}


# ------------------------------------------------------------------ data ----
def load_rank_matrix() -> dict[str, pd.DataFrame]:
    """species -> tools x metrics matrix of mean ranks (Precision-ordered)."""
    d = pd.read_csv(TAB / "fig5a_mean_rank_by_group_tool.tsv", sep="\t")
    out = {}
    for sp in SPECIES:
        s = d[d["species"] == sp]
        mat = s.pivot_table(index="tool", columns="metric", values="rank_mean")
        mat = mat.reindex(columns=METRICS)
        order = mat[SORT_METRIC].sort_values(kind="stable").index
        mat = mat.reindex(list(order) + [t for t in sorted(set(d["tool"]))
                                         if t not in order])
        assert len(mat) == len(TOOL_ORDER), f"{sp}: {len(mat)} tools"
        out[sp] = mat
    return out


def load_bins() -> pd.DataFrame:
    return pd.read_csv(TAB / "fig5b_bin_summary.tsv", sep="\t")


def load_window() -> pd.DataFrame:
    return pd.read_csv(TAB / "fig5cd_window_sweep_summary.tsv", sep="\t")


# ----------------------------------------------------------------- layout ---
def _rect(x_in: float, y_in: float, w_in: float, h_in: float,
          fig_w: float, fig_h: float) -> list[float]:
    return [x_in / fig_w, y_in / fig_h, w_in / fig_w, h_in / fig_h]


def build_layout(fig_h: float) -> dict[str, list[float]]:
    """Axes rectangles ``[y0, y1]`` in inches, top-down, for the 4 rows."""
    rects: dict[str, list[float]] = {}
    y = fig_h - TOP
    for i, row in enumerate(("A", "B", "C", "D")):
        y -= ROW_TITLE
        rects[row] = [y - ROWS_H[row], y]
        y -= ROWS_H[row]
        if i < 3:
            y -= ROW_GAP[row]
    return rects


# ------------------------------------------------------------- A row --------
def draw_row_a(fig, rects, mats, fig_w, fig_h, logger):
    y0, y1 = rects["A"]
    h = y1 - y0
    cmap = sns.diverging_palette(10, 240, as_cmap=True)     # rank 1 = red
    cmap.set_bad("0.86")
    last = None
    axes_a = []
    for i, sp in enumerate(SPECIES):
        mat = mats[sp].astype(float)
        x = A_LEFT + i * (A_W + A_LABEL_GAP + A_GAP)
        ax = fig.add_axes(_rect(x, y0, A_W, h, fig_w, fig_h))
        axes_a.append(ax)
        arr = np.ma.masked_invalid(mat.to_numpy(dtype=float))
        last = ax.imshow(arr, cmap=cmap, vmin=1, vmax=len(mat), aspect="auto")
        for r in range(mat.shape[0]):
            for c in range(mat.shape[1]):
                v = mat.iloc[r, c]
                if pd.isna(v):
                    continue
                ax.text(c, r, f"{v:.0f}" if abs(v - round(v)) < 1e-9 else f"{v:.1f}",
                        ha="center", va="center", fontsize=FS["annot"], color="black")
        ax.set_xticks(range(len(METRICS)))
        ax.set_xticklabels([METRIC_LABEL[m] for m in METRICS], rotation=45,
                           ha="right", fontsize=FS["tick"])
        ax.set_yticks(range(mat.shape[0]))
        ax.set_yticklabels([DISPLAY.get(t, t) for t in mat.index],
                           fontsize=FS["tick"])
        ax.tick_params(length=2, pad=1.5)
        for s in ax.spines.values():
            s.set_visible(True)
            s.set_linewidth(0.6)
        ax.set_title(TITLE[sp], fontsize=FS["title"], fontweight="bold", pad=3)
        # the ordered tool list is the panel's key: log it
        logger.info("A/%s order (%s rank): %s", sp, SORT_METRIC,
                    " > ".join(mat.index.astype(str)))
    cax = fig.add_axes(_rect(A_LEFT + 3 * A_W + 2 * (A_LABEL_GAP + A_GAP) + 0.02,
                             y0 + 0.10, 0.075, h - 0.20, fig_w, fig_h))
    cb = fig.colorbar(last, cax=cax)
    cb.set_label("Mean rank  (1 = best)", fontsize=FS["cbar"], labelpad=2)
    cb.ax.tick_params(labelsize=FS["tick"], length=2)
    cb.outline.set_linewidth(0.6)
    return axes_a


# --------------------------------------------------------- B / C / D --------
def _sp_axes(fig, rects, row, fig_w, fig_h):
    y0, y1 = rects[row]
    w = (CANVAS_W - L_LEFT - L_RIGHT - 2 * L_GAP) / 3
    return [fig.add_axes(_rect(L_LEFT + i * (w + L_GAP), y0, w, y1 - y0,
                               fig_w, fig_h)) for i in range(3)], w


def draw_row_b(fig, rects, bins, fig_w, fig_h):
    axes, _ = _sp_axes(fig, rects, "B", fig_w, fig_h)
    x = np.arange(len(BINS))
    for ax, sp in zip(axes, SPECIES):
        for tool in TOOL4:
            s = (bins[(bins["species"] == sp) & (bins["tool"] == tool)]
                 .set_index("mod_ratio_bin").reindex(BINS))
            y = s["hit_rate_mean"].to_numpy(dtype=float)
            lo = s["hit_rate_min"].to_numpy(dtype=float)
            hi = s["hit_rate_max"].to_numpy(dtype=float)
            ok = np.isfinite(y)
            if not ok.any():
                continue
            if np.isfinite(lo).sum() > 1:
                ax.fill_between(x[ok], lo[ok], hi[ok], color=TOOL4_COLOR[tool],
                                alpha=.15, lw=0)
            ax.plot(x[ok], y[ok], "-o", color=TOOL4_COLOR[tool], ms=2.6, lw=1.2,
                    label=tool)
        ax.set_xticks(x)
        ax.set_xticklabels(BINS, rotation=45, ha="right", fontsize=FS["tick"])
        ax.set_ylim(0, 1)
        ax.set_xlim(-0.4, len(BINS) - 0.6)
        ax.set_title(TITLE[sp], fontsize=FS["title"], fontweight="bold", pad=3)
        ax.tick_params(labelsize=FS["tick"], length=2, pad=1.5)
        # every species panel carries the x axis title (user 2026-09-19), like
        # the C and D rows do
        ax.set_xlabel("Tool-reported modification ratio", fontsize=FS["label"])
    axes[0].set_ylabel(f"PPV vs. GLORI ({C.PRIMARY_WINDOW} bp)",
                       fontsize=FS["label"])
    
    #: original figure carried one on each of the three; every panel gets its own.
    for ax in axes:
        ax.legend(frameon=False, fontsize=FS["legend"], loc="upper left",
                  handlelength=1.3, labelspacing=.25, borderpad=.1,
                  handletextpad=.4)
    return axes


def draw_window_row(fig, rects, curve, metric, ylabel, row, fig_w, fig_h):
    axes, _ = _sp_axes(fig, rects, row, fig_w, fig_h)
    handles = []
    for ax, sp in zip(axes, SPECIES):
        for tool in TOOL_ORDER:
            s = (curve[(curve["species"] == sp) & (curve["tool"] == tool)
                       & (curve["metric"] == metric)]
                 .set_index("window").reindex(WINDOWS))
            y = s["value_mean"].to_numpy(dtype=float)
            lo = s["value_min"].to_numpy(dtype=float)
            hi = s["value_max"].to_numpy(dtype=float)
            if not np.isfinite(y).any():
                continue
            if np.isfinite(lo).sum() > 1:
                ax.fill_between(WINDOWS, lo, hi, color=TOOL_COLOR[tool],
                                alpha=.12, lw=0)
            ln, = ax.plot(WINDOWS, y, "-o", color=TOOL_COLOR[tool], ms=2.4,
                          lw=1.1, label=DISPLAY.get(tool, tool))
            if sp == SPECIES[0]:
                handles.append(ln)
        ax.axvline(C.PRIMARY_WINDOW, color="0.45", lw=0.7, ls=(0, (3, 2)), zorder=0)
        ax.set_xscale("symlog", linthresh=1)
        ax.minorticks_off()
        ax.set_xlim(-0.15, 60)
        ax.set_ylim(0, 1)
        ax.set_xticks(WINDOWS)
        ax.set_xticklabels([str(w) for w in WINDOWS], fontsize=FS["tick"])
        ax.set_title(TITLE[sp], fontsize=FS["title"], fontweight="bold", pad=3)
        ax.tick_params(labelsize=FS["tick"], length=2, pad=1.5)
        ax.set_xlabel("Matching window (bp)", fontsize=FS["label"])
    axes[0].set_ylabel(ylabel, fontsize=FS["label"])
    
    #: by the caller in the service strip under the C row.  A 13-entry key inside
    #: the axes overruns the panel width and covers the curves, so the window rows
    #: keep their keys outside the axes.
    return handles, axes


# ------------------------------------------------------------- checks -------
def assert_min_fontsize(fig, min_pt: float = MIN_PT) -> float:
    sizes = [t.get_fontsize() for t in fig.findobj(Text)
             if t.get_text().strip() and t.get_visible()]
    smallest = min(sizes) if sizes else float("nan")
    assert smallest >= min_pt, f"smallest text {smallest:.2f} pt < {min_pt} pt"
    return smallest


def assert_label_fit(fig, a_axes, b_axes, cd_axes, logger) -> None:
    """Layout guards that a re-render must keep passing.

    * the widest A-row tool name (plus tick pad) must fit inside the reserved
      label gap and inside the left margin -- the failure that made the first
      version overlap the neighbouring heatmap;
    * every B / C / D panel must carry its own x axis title (user 2026-09-19);
    * the 45 deg bin labels of the B row must stay separable: the perpendicular
      distance between neighbouring labels has to exceed the text height.
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    widest, widest_txt = 0.0, ""
    for ax in a_axes:
        for t in ax.get_yticklabels():
            if not t.get_text().strip():
                continue
            w = t.get_window_extent(renderer=rend).width / fig.dpi
            if w > widest:
                widest, widest_txt = w, t.get_text()
    need = widest + 0.02                      # tick pad + breathing room
    assert need <= A_LABEL_GAP, (
        f"A-row tool label {widest_txt!r} = {widest:.3f} in (+pad {need:.3f} in) "
        f"does not fit the reserved gap {A_LABEL_GAP:.2f} in")
    assert need <= A_LEFT, (
        f"A-row first-panel label needs {need:.3f} in but A_LEFT = {A_LEFT:.2f} in")
    logger.info("A-row widest tool label %r = %.3f in (+pad %.3f in) fits gap %.2f in",
                widest_txt, widest, need, A_LABEL_GAP)

    for name, axes in (("B", b_axes), ("C/D", cd_axes)):
        missing = [i for i, ax in enumerate(axes) if not ax.get_xlabel().strip()]
        assert not missing, f"{name} row: panels {missing} have no x axis title"
    logger.info("x axis titles present on all %d B and %d C/D panels",
                len(b_axes), len(cd_axes))

    pitch = (b_axes[0].get_position().width * CANVAS_W) / len(BINS)
    clear = pitch * float(np.sin(np.deg2rad(45)))
    text_h = FS["tick"] / 72.0
    assert clear > text_h * 1.05, (
        f"B-row 45 deg bin labels would collide: clearance {clear:.3f} in "
        f"vs text height {text_h:.3f} in")
    logger.info("B-row 45 deg label clearance %.3f in > text height %.3f in",
                clear, text_h)


def assert_data_integrity(mats, bins, curve, logger) -> None:
    for sp, mat in mats.items():
        assert list(mat.columns) == METRICS
        assert mat.notna().all().all(), f"{sp}: missing mean rank"
        # ordering key took effect: SORT_METRIC rank must not increase downwards
        col = mat[SORT_METRIC].to_numpy(dtype=float)
        assert (np.diff(col) >= -1e-9).all(), f"{sp}: rows not {SORT_METRIC}-ordered"
    n_bins = bins.groupby(["species", "tool"])["hit_rate_mean"].apply(
        lambda s: int(s.notna().sum()))
    assert n_bins.min() >= 5, f"a B-panel curve has only {n_bins.min()} bins"
    assert set(curve["metric"]) == {"hit_rate", "localization_accuracy"}
    assert sorted(curve["window"].unique()) == sorted(WINDOWS)
    logger.info("integrity: %d species matrices %s-ordered, %d B curves >= 5 bins,"
                " windows %s", len(mats), SORT_METRIC, len(n_bins), WINDOWS)


# ---------------------------------------------------------------- main ------
def main() -> None:
    logger = setup_logger("39_fig5_assembled")
    apply_style()

    fig_h = (TOP + 4 * ROW_TITLE + sum(ROWS_H.values())
             + sum(ROW_GAP.values()) + LEGEND_H)
    assert fig_h <= MAX_H, f"canvas {fig_h:.2f} in exceeds {MAX_H} in"

    mats = load_rank_matrix()
    bins = load_bins()
    curve = load_window()
    assert_data_integrity(mats, bins, curve, logger)

    fig = plt.figure(figsize=(CANVAS_W, fig_h))
    rects = build_layout(fig_h)

    # panel letters sit in the gap above each row title
    letter_y = {}
    y = fig_h - TOP
    for i, row in enumerate(("A", "B", "C", "D")):
        letter_y[row] = y - 0.02
        y -= ROW_TITLE + ROWS_H[row]
        if i < 3:
            y -= ROW_GAP[row]
    for row, y_in in letter_y.items():
        fig.text(0.012, y_in / fig_h, row, fontsize=FS["letter"],
                 fontweight="bold", va="top", ha="left")

    axes_a = draw_row_a(fig, rects, mats, CANVAS_W, fig_h, logger)
    axes_b = draw_row_b(fig, rects, bins, CANVAS_W, fig_h)
    handles, axes_c = draw_window_row(fig, rects, curve, "hit_rate",
                                      f"PPV vs. GLORI ({C.PRIMARY_WINDOW} bp)", "C",
                                      CANVAS_W, fig_h)
    
    #: the C row (its own row, so C is not left without a key like before).
    fig.legend(handles=handles, labels=[DISPLAY.get(t, t) for t in TOOL_ORDER],
               loc="upper center",
               bbox_to_anchor=(0.5, (rects["C"][0] - 0.30) / fig_h), ncol=7,
               frameon=False, fontsize=FS["legend"], handlelength=1.5,
               columnspacing=1.1, handletextpad=0.45, labelspacing=0.35)
    _, axes_d = draw_window_row(fig, rects, curve, "localization_accuracy",
                                "Exact localisation", "D", CANVAS_W, fig_h)
    fig.legend(handles=handles, labels=[DISPLAY.get(t, t) for t in TOOL_ORDER],
               loc="lower center", bbox_to_anchor=(0.5, 0.004), ncol=7,
               frameon=False, fontsize=FS["legend"], handlelength=1.5,
               columnspacing=1.1, handletextpad=0.45, labelspacing=0.35)

    smallest = assert_min_fontsize(fig)
    assert_label_fit(fig, axes_a, axes_b, axes_c + axes_d, logger)
    for row in ("A", "B", "C", "D"):
        y0, y1 = rects[row]
        logger.info("row %s: data height %.2f in", row, y1 - y0)
    logger.info("canvas %.2f x %.2f in (0.95\\textwidth = %.2f in, \\textheight = 9.21 in)",
                CANVAS_W, fig_h, CANVAS_W)
    logger.info("smallest text element: %.2f pt (limit %.1f pt)", smallest, MIN_PT)

    FIG.mkdir(parents=True, exist_ok=True)
    stem = FIG / "Figure5_rev"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches=None, pad_inches=0)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches=None, pad_inches=0)
    plt.close(fig)
    logger.info("wrote %s.pdf + .png", stem)
    logger.info("output: %s", FIG)


if __name__ == "__main__":
    main()
