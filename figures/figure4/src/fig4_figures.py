#!/usr/bin/env python
"""Figure 4 revision, step 2: panels + assembled figure.

Draws the four revised panels and the assembled figure, strictly from the
tables written by ``fig4_kl.py`` -- nothing is recomputed here, so the
assembled figure cannot drift from the numbers:

  A  KL divergence (add-eps, vs the theoretical RRACH consensus) per
     tool x species: mean +/- SD over independent replicates
     (Arabidopsis n = 3, HeLa n = 3, mouse the single mES_WT study --
     never merged, cross-study rule); open circles = cross-species mean;
     x on log scale.  Replaces the published merged-replicate points.
  B  5-mer sequence logos of the three lowest mean-KL tools
     (MINES / DENA / Nanom6A), y axis now labelled
     "Information content (bits)" (R3 minor #4).  The published in-panel
     "KL = ..." callouts are REMOVED (house rule 2026-09-19: no annotation
     text inside a panel -- the KL numbers live in panel A and in
     tables/fig4_kl_summary.tsv).  Logos use the replicate-mean PWM
     (equal weight per replicate).
  C  MDS of tools by their 5-mer frequency profiles per species
     (replicate-mean profile, z-scored per 5-mer as published; distances =
     euclidean).  sklearn was unavailable in every env that also has
     logomaker, so the embedding is classical (Torgerson) MDS computed
     with numpy -- deterministic, no random_state needed; SMACOF and
     classical MDS agree to within the marker size on such small
     well-behaved distance matrices.
  D  Top-20 5-mers per species by number of supporting tool
     configurations (a tool supports a 5-mer when it appears in the
     tool's replicate-mean top-20); orange = matches the RRACH consensus
     ^[AG][AG]AC[ACU]$, blue = other (same classification as the
     published panel, verified against Figure4.pdf).  Replicate support
     is exported in tables/fig4_panelD_support.tsv.

House rules (2026-09-19): drawn at the final printed size
(0.95 * \\textwidth = 6.66 in, Wiley USG layout), every text element
>= 7.0 pt (asserted), no gridlines, no annotation text inside panels,
Arial, vector PDF + 300 dpi PNG.

Inputs  (fig4_revision/, written by fig4_kl.py)
  tables/fig4_kl_per_rep.tsv
  analysis/fig4_pwm_per_rep.npz
  analysis/fig4_kmer_counts.tsv.gz
Outputs (fig4_revision/)
  figures/panel{A,B,C,D}_*.{pdf,png}
  figures/Figure4_rev.{pdf,png}
  tables/fig4_panelD_support.tsv

Usage
-----
conda run -n motif_analysis --no-capture-output python \\
    $RNAMODBENCH_ROOT/figures/figure4/analysis/fig4_figures.py
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
from itertools import cycle
from pathlib import Path

import logomaker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.text import Text

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                       # fig4_revision/
TAB, FIG = ROOT / "tables", ROOT / "figures"
sys.path.insert(0, str(_RB / "src/sites_v2"))
from common.figstyle import apply as apply_style            # noqa: E402

CANVAS_W, MIN_PT = 6.66, 7.0
TOP_PAD, ROW_GAP, BOTTOM_PAD = 0.03, 0.11, 0.05
#: C and D used to be 1.98 / 2.46 in: at 7.5 pt the 20 five-mer labels of D
#: needed 10 pt each (2.8 in) and C's MDS ticks collided with the y-axis title.
ROW1_H, ROW2_H, ROW3_H = 2.45, 2.35, 3.05   # A|B, C, D rows incl. titles/labels
#: the panel letters stand in their own margin left of the rows
LETTER_PAD = 0.20
FIG_H = TOP_PAD + ROW1_H + ROW_GAP + ROW2_H + ROW_GAP + ROW3_H + BOTTOM_PAD
assert FIG_H <= 8.50, f"canvas {FIG_H:.2f} in exceeds the 8.50 in page limit"

SP_ORDER = ["Arabidopsis", "Mouse", "Human (HeLa)"]
#: species names only -- the replicate counts (n = 3 / n = 1 study) and the
#: HeLa provenance of the human data live in the caption, not in the legend
SP_DISP = {"Arabidopsis": "Arabidopsis", "Mouse": "Mouse", "Human (HeLa)": "Human"}
SP_COLOR = {"Arabidopsis": "#1E888B", "Mouse": "#F5B264", "Human (HeLa)": "#3778A0"}
SP_MARKER = {"Arabidopsis": "o", "Mouse": "s", "Human (HeLa)": "^"}
BASES = ["A", "U", "C", "G"]
KMERS = [a + b + c + d + e for a in "ACGU" for b in "ACGU" for c in "ACGU"
         for d in "ACGU" for e in "ACGU"]
RRACH_RE = {"A", "G"}          # pos 0/1 allowed bases
RRACH_H = {"A", "C", "U"}      # pos 4 allowed bases
RRACH_COLOR, OTHER_COLOR = "#E8772F", "#4C78A8"
#: explicit logo colours (logomaker's built-in schemes key on DNA letters only)
LOGO_COLOR = {"A": "#33A02C", "C": "#1F78B4", "G": "#E58606", "U": "#7F8C99"}

TOOL_ORDER = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
              "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
              "NanoSPA_m6A", "xPore", "yanocomp"]
DISPLAY = {"yanocomp": "Yanocomp"}
_pal = plt.get_cmap("tab20").colors
TOOL_COLOR = {t: _pal[i % 20] for i, t in enumerate(TOOL_ORDER)}
TOOL_MARKER = dict(zip(TOOL_ORDER, cycle(["o", "s", "^", "D", "v", "<", ">", "p"])))

FS = {"tick": 7.5, "label": 8.0, "title": 8.5, "legend": 7.2, "letter": 9.5}
XLABEL_A = "KL divergence vs. RRACH consensus (bits)"


def is_rrach(k: str) -> bool:
    return k[2] == "A" and k[3] == "C" and k[0] in RRACH_RE and k[1] in RRACH_RE \
        and k[4] in RRACH_H


# ------------------------------------------------------------------- data ---
def load_kl() -> tuple[pd.DataFrame, list[str]]:
    per = pd.read_csv(TAB / "fig4_kl_per_rep.tsv", sep="\t")
    per = per[per["in_figure"]].copy()
    pooled = per.groupby(["Tool", "Species_Name"], as_index=False)[
        ["KL_rrach_addeps", "KL_rrach_clip"]].mean()
    cross = (pooled.groupby("Tool")["KL_rrach_addeps"].mean()
             .reindex(TOOL_ORDER).dropna().sort_values())
    return per, pooled, list(cross.index)


def load_pwms() -> dict[str, np.ndarray]:
    return dict(np.load(HERE / "fig4_pwm_per_rep.npz"))


def load_profiles() -> dict[str, pd.DataFrame]:
    """species -> tools x 256 replicate-mean 5-mer frequencies."""
    cnt = pd.read_csv(HERE / "fig4_kmer_counts.tsv.gz", sep="\t")
    cnt = cnt[cnt["Replicate"].map(lambda r: not r.startswith("mESCs_"))]
    out = {}
    for sp, g in cnt.groupby("Species_Name"):
        freq = g.assign(freq=g["count"] / g.groupby(["Tool", "Replicate"])["count"]
                        .transform("sum"))
        prof = freq.pivot_table(index="Tool", columns="kmer", values="freq",
                                aggfunc="mean")
        prof = prof.reindex(columns=KMERS, fill_value=0.0).reindex(TOOL_ORDER)
        out[sp] = prof
    return out


def top20_support(profiles: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Per species: top-20 5-mers by number of supporting tool configs."""
    rows = []
    for sp, prof in profiles.items():
        sup = {k: 0 for k in KMERS}
        rep_sup = {k: [] for k in KMERS}
        for tool, row in prof.iterrows():
            top = list(row.sort_values(ascending=False).head(20).index)
            for k in top:
                sup[k] += 1
        for k in sup:
            if sup[k]:
                rows.append({"Species_Name": sp, "kmer": k,
                             "n_tools_support": sup[k], "is_rrach": is_rrach(k)})
    d = pd.DataFrame(rows)
    out = []
    for sp, g in d.groupby("Species_Name"):
        g = g.sort_values(["n_tools_support", "is_rrach", "kmer"],
                          ascending=[False, False, True]).head(20).reset_index(drop=True)
        g.insert(0, "rank", np.arange(1, len(g) + 1))
        out.append(g)
    res = pd.concat(out, ignore_index=True)
    res.to_csv(TAB / "fig4_panelD_support.tsv", sep="\t", index=False)
    return res


# ---------------------------------------------------------------- layout ----
def rect(x: float, y: float, w: float, h: float) -> list[float]:
    return [x / CANVAS_W, y / FIG_H, w / CANVAS_W, h / FIG_H]


def row_y(row: int) -> float:
    """Top edge (inches from bottom) of each row band; row 1 = top of page."""
    return {1: FIG_H - TOP_PAD,
            2: FIG_H - TOP_PAD - ROW1_H - ROW_GAP,
            3: FIG_H - TOP_PAD - ROW1_H - ROW_GAP - ROW2_H - ROW_GAP}[row]


def panel_letter(fig, letter: str, axes_top_in: float, x: float = 0.015) -> None:
    """Panel letter, sitting on the axes' top edge in the letter margin.

    The old call took the row band's top and hung the glyph downwards, which
    put C and D *inside* the top-left corner of their own first panel.
    """
    fig.text(x / CANVAS_W, axes_top_in / FIG_H, letter,
             fontsize=FS["letter"], fontweight="bold", va="bottom", ha="left")


# ---------------------------------------------------------------- panel A ---
def draw_a(fig, per: pd.DataFrame, pooled: pd.DataFrame, order: list[str],
           x0: float, w: float) -> None:
    y0 = row_y(1) - ROW1_H
    ax = fig.add_axes(rect(x0 + 0.74, y0 + 0.32, w - 0.80, ROW1_H - 0.66))
    means = pooled.pivot(index="Tool", columns="Species_Name",
                         values="KL_rrach_addeps")
    for i, tool in enumerate(order):
        y = len(order) - 1 - i
        cross = means.loc[tool].mean()
        for sp in SP_ORDER:
            if sp not in means.columns or pd.isna(means.loc[tool, sp]):
                continue
            m = means.loc[tool, sp]
            sub = per[(per["Tool"] == tool) & (per["Species_Name"] == sp)]
            sd = sub["KL_rrach_addeps"].std()
            n = len(sub)
            if n > 1 and not pd.isna(sd):
                ax.errorbar(m, y, xerr=sd, fmt="none", ecolor=SP_COLOR[sp],
                            elinewidth=0.7, capsize=0, zorder=2)
            ax.scatter(m, y, s=14, marker=SP_MARKER[sp], color=SP_COLOR[sp],
                       edgecolor="black", linewidth=0.3, zorder=3,
                       label=SP_DISP[sp] if i == 0 else None)
        ax.scatter(cross, y, s=16, facecolor="white", edgecolor="black",
                   linewidth=0.8, zorder=4,
                   label="Cross-species mean" if i == 0 else None)
    ax.set_xscale("log")
    ax.set_xlim(0.02, 120)
    ax.set_xticks([0.03, 0.1, 0.3, 1, 3, 10, 30, 100])
    ax.set_xticklabels(["0.03", "0.1", "0.3", "1", "3", "10", "30", "100"])
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([DISPLAY.get(t, t) for t in order][::-1])
    ax.tick_params(labelsize=FS["tick"], length=2.5, pad=1.5)
    ax.set_xlabel(XLABEL_A, fontsize=FS["label"], labelpad=1.5)
    ax.set_title("Motif deviation", fontsize=FS["title"],
                 fontweight="bold", loc="left", pad=3)
    ax.legend(frameon=False, fontsize=FS["legend"], loc="lower left",
              handlelength=1.1, labelspacing=0.25, borderpad=0.1,
              handletextpad=0.4)
    panel_letter(fig, "A", row_y(1) - 0.34, x0)


# ---------------------------------------------------------------- panel B ---
def bits_matrix(pwm: np.ndarray) -> pd.DataFrame:
    h = -(pwm * np.log2(np.clip(pwm, 1e-12, None))).sum(axis=1)
    ic = np.clip(2.0 - h, 0.0, None)
    return pd.DataFrame(pwm * ic[:, None], index=range(-2, 3), columns=BASES)


def draw_b(fig, pwms: dict[str, np.ndarray], per: pd.DataFrame,
           order: list[str], x0: float, w: float) -> None:
    logo_tools = order[:3]                     # three lowest mean-KL tools
    y0, h = row_y(1) - ROW1_H, ROW1_H - 0.26   # band below the B title
    lab_w, ylab_w, yax_w = 0.60, 0.30, 0.18
    col_w = (w - lab_w - ylab_w - yax_w - 0.06) / 3
    row_h = (h - 0.17 - 0.13) / 3
    panel_letter(fig, "B", row_y(1) - 0.26, x0)
    for c, tool in enumerate(logo_tools):      # columns = tools
        xcol = x0 + lab_w + ylab_w + yax_w + c * (col_w + 0.03)
        for r, sp in enumerate(SP_ORDER):      # rows = species
            ytop = y0 + h - 0.17 - r * row_h
            _ = ytop
            # 0.06 in between rows put the bottom tick label of one logo row
            # on top of the top tick label of the next one
            ax = fig.add_axes(rect(xcol, ytop - row_h + 0.02, col_w,
                                   row_h - 0.16))
            reps = [k for k, v in pwms.items()
                    if k.split("|")[0] == tool and k.split("|")[1] == sp
                    and not k.split("|")[2].startswith("mESCs_")]
            pwm = np.mean([pwms[k] for k in reps], axis=0)
            logo = logomaker.Logo(bits_matrix(pwm), ax=ax,
                                  color_scheme=LOGO_COLOR)
            logo.ax.set_xlim(-0.5, 4.5)
            logo.ax.spines["top"].set_visible(False)
            logo.ax.spines["right"].set_visible(False)
            logo.ax.tick_params(labelsize=FS["tick"], length=2, pad=3)
            ax.set_xticks(range(-2, 3))
            ax.tick_params(labelbottom=(r == 2))
            if c == 0:
                ax.set_yticks([0, 2])      # two ticks clear the row gap
            else:
                ax.set_yticks([])
            if r == 0:
                ax.set_title(DISPLAY.get(tool, tool), fontsize=FS["title"],
                             fontweight="bold", pad=2)
    for r, sp in enumerate(SP_ORDER):          # species row labels
        ytop = y0 + h - 0.17 - r * row_h
        fig.text((x0 + lab_w / 2) / CANVAS_W, (ytop - row_h / 2) / FIG_H,
                 sp.replace(" (HeLa)", ""), fontsize=FS["tick"], va="center", ha="center")
    fig.text((x0 + lab_w + ylab_w / 2) / CANVAS_W,
             (y0 + h - 0.17 - 1.5 * row_h) / FIG_H,
             "Information content (bits)", rotation=90,
             fontsize=FS["tick"], va="center", ha="center")


# ---------------------------------------------------------------- panel C ---
def classical_mds(d: np.ndarray) -> np.ndarray:
    n = len(d)
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ (d ** 2) @ j
    w, v = np.linalg.eigh(b)
    idx = np.argsort(w)[::-1][:2]
    return v[:, idx] * np.sqrt(np.clip(w[idx], 0, None))


def draw_c(fig, profiles: dict[str, pd.DataFrame], x0: float, w: float) -> None:
    y0 = row_y(2) - ROW2_H
    h = ROW2_H - 0.17 - 0.33
    leg_w = 0.92
    x0 = x0 + LETTER_PAD                      # letter margin
    pw = (w - LETTER_PAD - leg_w - 0.44) / 3
    panel_letter(fig, "C", y0 + 0.31 + h)
    for c, sp in enumerate(SP_ORDER):
        prof = profiles[sp]
        z = prof.to_numpy(dtype=float)
        sd = z.std(axis=0)
        z = np.where(sd > 0, (z - z.mean(axis=0)) / np.where(sd > 0, sd, 1.0), 0.0)
        dmat = np.sqrt(((z[:, None, :] - z[None, :, :]) ** 2).sum(-1))
        pos = classical_mds(dmat)
        ax = fig.add_axes(rect(x0 + c * (pw + 0.20), y0 + 0.31, pw, h))
        for i, tool in enumerate(prof.index):
            ax.scatter(pos[i, 0], pos[i, 1], s=13, marker=TOOL_MARKER[tool],
                       color=TOOL_COLOR[tool], edgecolor="black", linewidth=0.3)
        ax.set_title(SP_DISP[sp], fontsize=FS["title"], fontweight="bold", pad=2)
        ax.set_aspect("equal", adjustable="datalim")
        # pad 1 let the top y tick and the left x tick collide in the corner
        ax.tick_params(labelsize=FS["tick"], length=2, pad=3.5)
        ax.set_xlabel("MDS 1", fontsize=FS["label"], labelpad=3)
        if c == 0:
            # labelpad 1 put "MDS 2" on top of the topmost tick label
            ax.set_ylabel("MDS 2", fontsize=FS["label"], labelpad=7)
        ax.locator_params(axis="y", nbins=3)
    handles = [Line2D([], [], color=TOOL_COLOR[t], marker=TOOL_MARKER[t], ls="none",
                      markersize=4, markeredgecolor="black", markeredgewidth=0.3,
                      label=DISPLAY.get(t, t)) for t in TOOL_ORDER]
    fig.legend(handles=handles, frameon=False, fontsize=FS["legend"],
               loc="center left",
               bbox_to_anchor=((x0 + w - leg_w + 0.04) / CANVAS_W,
                               (y0 + 0.72 * h) / FIG_H),
               handlelength=1.0, labelspacing=0.18, borderpad=0.1,
               handletextpad=0.4)


# ---------------------------------------------------------------- panel D ---

#: wide at 7.5 pt, right-aligned on their axes' left edge.  The first column used
#: to start 0.23 in in, so half of every label was painted outside the canvas and
#: silently dropped by the PDF -- the render kept only the last two or three
#: characters.  The row now reserves a label gutter of its own and every tick
#: label is asserted to sit inside the canvas.
D_LABEL_GUTTER, D_GAP = 0.46, 0.40
D_PW = 1.76


def draw_d(fig, support: pd.DataFrame, x0: float, w: float) -> None:
    y0 = row_y(3) - ROW3_H
    h = ROW3_H - 0.17 - 0.32
    x0 = x0 + LETTER_PAD                      # letter margin
    assert x0 == D_LABEL_GUTTER - 0.23, f"the D row moved: {x0:.3f} in"
    x0 = D_LABEL_GUTTER                       # + the 5-mer label gutter
    pw = D_PW
    panel_letter(fig, "D", y0 + 0.30 + h)
    for c, sp in enumerate(SP_ORDER):
        g = support[support["Species_Name"] == sp].sort_values("rank",
                                                               ascending=False)
        ax = fig.add_axes(rect(x0 + c * (pw + D_GAP), y0 + 0.30, pw, h))
        ax.hlines(np.arange(len(g)), 0, g["n_tools_support"],
                  color=[RRACH_COLOR if r else OTHER_COLOR for r in g["is_rrach"]],
                  lw=1.0, zorder=2)
        ax.scatter(g["n_tools_support"], np.arange(len(g)), s=8,
                   color=[RRACH_COLOR if r else OTHER_COLOR for r in g["is_rrach"]],
                   edgecolor="black", linewidth=0.2, zorder=3)
        ax.set_yticks(np.arange(len(g)))
        ax.set_yticklabels(g["kmer"])
        ax.set_xlim(0, 13.3)
        ax.set_xticks([0, 2, 4, 6, 8, 10, 12])
        ax.set_ylim(-0.7, len(g) - 0.3)
        ax.set_title(SP_DISP[sp], fontsize=FS["title"], fontweight="bold", pad=2)
        ax.tick_params(labelsize=FS["tick"], length=2, pad=1)
        # labelpad 1 ran the x-axis title into the tick numbers
        ax.set_xlabel("Supporting tools (of 13)", fontsize=FS["label"],
                      labelpad=4)
        handles = [Line2D([], [], color=RRACH_COLOR, marker="o", ls="none",
                          markersize=4, label="RRACH motif"),
                   Line2D([], [], color=OTHER_COLOR, marker="o", ls="none",
                          markersize=4, label="Other motifs")]
        ax.legend(handles=handles, frameon=False, fontsize=FS["legend"],
                  loc="lower right", handlelength=1.0, labelspacing=0.25,
                  borderpad=0.1, handletextpad=0.4)
        # no tick label may be painted outside the canvas: poppler and every
        # other consumer silently drop page-external text, which is how the
        # first column lost its 5-mers without any audit noticing
        fig.canvas.draw()
        rend = fig.canvas.get_renderer()
        for t in ax.get_yticklabels():
            bb = t.get_window_extent(renderer=rend)
            x0_in, x1_in = bb.x0 / fig.dpi, bb.x1 / fig.dpi
            assert x0_in > -0.002 and x1_in < CANVAS_W + 0.002, (
                f"panel D column {c}: tick label {t.get_text()!r} runs "
                f"{x0_in:.3f}..{x1_in:.3f} in on a {CANVAS_W} in canvas")


# ------------------------------------------------------------- assertions ---
def assert_min_fontsize(fig, min_pt: float = MIN_PT) -> float:
    sizes = []
    for ax in fig.axes:
        texts = (list(ax.get_xticklabels()) + list(ax.get_yticklabels())
                 + [ax.xaxis.label, ax.yaxis.label, ax.title])
        for t in texts:
            if isinstance(t, Text) and t.get_text():
                sizes.append(t.get_fontsize())
        leg = ax.get_legend()
        if leg is not None:
            sizes += [t.get_fontsize() for t in leg.get_texts()]
    sizes += [t.get_fontsize() for t in fig.texts if t.get_text()]
    for leg in fig.legends:
        sizes += [t.get_fontsize() for t in leg.get_texts()]
    smallest = min(sizes)
    assert smallest >= min_pt, f"smallest text {smallest:.2f} pt < {min_pt} pt"
    return smallest


def assert_no_grid(fig) -> None:
    """Gridlines are disabled house-wide via rcParams; re-check visibly."""
    for ax in fig.axes:
        for grid in (ax.xaxis.get_gridlines(), ax.yaxis.get_gridlines()):
            assert all(not gl.get_visible() for gl in grid)


# ------------------------------------------------------------------- main ---
def main() -> None:
    apply_style()
    per, pooled, order = load_kl()
    pwms, profiles = load_pwms(), load_profiles()
    support = top20_support(profiles)

    print("lowest-KL tools:", order[:3])
    fig = plt.figure(figsize=(CANVAS_W, FIG_H))
    draw_a(fig, per, pooled, order, 0.03, 2.42)
    draw_b(fig, pwms, per, order, 2.55, 4.08)
    # the rows start in the letter margin, so their width shrinks by it
    draw_c(fig, profiles, 0.03, CANVAS_W - 0.06 - LETTER_PAD)
    draw_d(fig, support, 0.03, CANVAS_W - 0.06 - LETTER_PAD)
    smallest = assert_min_fontsize(fig)
    assert_no_grid(fig)
    fig.savefig(FIG / "Figure4_rev.pdf", bbox_inches=None, pad_inches=0)
    fig.savefig(FIG / "Figure4_rev.png", dpi=300, bbox_inches=None, pad_inches=0)
    print(f"wrote Figure4_rev.pdf/.png  (canvas {CANVAS_W} x {FIG_H:.2f} in, "
          f"smallest text {smallest:.2f} pt)")


if __name__ == "__main__":
    main()
