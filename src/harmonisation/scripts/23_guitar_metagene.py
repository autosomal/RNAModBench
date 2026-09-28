#!/usr/bin/env python
"""Replicate-aware GUITAR metagene, stage 3: figures in the manuscript's style.

The panels reproduce the look of the published Fig. 3C/D (Guitar's
``GuitarPlot`` output): one panel per species, a smooth share-weighted density
over the concatenated axis ``1kb | 5'UTR | CDS | 3'UTR | 1kb``, dotted region
boundaries, the transcript schematic bar under the axis, filled areas, and the
key below each panel labelled with the sample-group names.

What the revision adds on top of that:

* the curve is the **majority consensus** of the biological replicates (the
  published curve came from a single replicate / an undocumented union);
* the thin same-hue lines are the **individual replicates**, so the replicate
  spread is readable inside the panel;
* the pooled curve concatenates every tool's sites, exactly as the published
  panel pooled the tools' BED files.

Figures in ``$RNAMODBENCH_LOCAL/guitar_metagene/figures``:

  FigR1_metagene_condition_{mrna,ncrna}   C-style: libraries pooled over tools
  FigR2_metagene_per_tool                 D-style: per tool x library
  FigR3_region_shares / FigR3b_...        region occupancy, with replicate SD
  FigR4_merge_rules                       union / majority / intersection
  FigR5_strand_mode                       strand-aware vs strand forced to "+"
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

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from common.config import SITES_ROOT                          # noqa: E402
from common.figstyle import SEGMENT_LABELS, apply as apply_style   # noqa: E402
from common.figstyle import guitar_panel, save                     # noqa: E402
from common.regionmodel import KINDS, NCR_ORDER, SEGMENT_ORDER  # noqa: E402

OUT = (_XB / "harmonisation/guitar_metagene")
TAB, FIG = OUT / "tables", OUT / "figures"
apply_style()

GRID = 200
POOLED = "ALL (pooled)"
SPECIES_ORDER = ["Arabidopsis", "Mouse", "Human"]
COND = {  # species -> (WT group, comparison group, comparison label)
    "Arabidopsis": ("Arabidopsis_WT", "Arabidopsis_KD", "fip37 KD"),
    "Mouse": ("Mouse_WT", "Mouse_KO", "Mettl3 KO"),
    "Human": ("HeLa_WT", "HeLa_IVT", "IVT"),
}
WT_COLOR, CMP_COLOR = "#4b81b8", "#e8a76b"
TOOL_ORDER = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
              "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
              "NanoSPA_m6A", "xPore", "yanocomp"]


# --------------------------------------------------------------------------- #
def load_density(tx_class: str, strand_mode: str) -> pd.DataFrame:
    d = pd.read_csv(TAB / "metagene_density.tsv.gz", sep="\t")
    d = d[(d.tx_class == tx_class) & (d.strand_mode == strand_mode)].copy()
    d["curve"] = [np.fromiter(map(float, s.split("|")), float, GRID)
                  for s in d.density]
    return d


def spine(d: pd.DataFrame, group: str, tool: str, merge: str,
          segs: list[str]) -> np.ndarray | None:
    """Per-segment densities concatenated onto the fixed 0..len(segs) axis."""
    out = np.zeros(len(segs) * GRID)
    seen = False
    for i, name in enumerate(segs):
        r = d[(d["group"] == group) & (d["tool"] == tool)
              & (d["merge"] == merge) & (d["kind"] == name)]
        if len(r):
            out[i * GRID:(i + 1) * GRID] = r.curve.iloc[0]
            seen = True
    return out if seen else None


def units_of(d: pd.DataFrame, group: str) -> list[str]:
    return sorted({m.split(":", 1)[1] for m in d.loc[d["group"] == group, "merge"]
                   if m.startswith("unit:")})


def consensus_word(d: pd.DataFrame, group: str) -> str:
    basis = d.loc[d["group"] == group, "merge_basis"].unique()
    if len(basis) == 1 and basis[0] == "studies":
        return "both studies"
    if len(basis) == 1 and basis[0] == "sequencing_runs":
        return "both runs"
    return "majority of replicates"


def tool_order(df: pd.DataFrame) -> list[str]:
    have = [t for t in TOOL_ORDER if t in set(df.tool)]
    return have + sorted(set(df.tool) - set(TOOL_ORDER))


# --------------------------------------------------------------------------- #
def _curve_with_replicates(ax, d, group, tool, segs, color, label):
    """Consensus curve + light fill + one thin line per replicate."""
    main = spine(d, group, tool, "majority", segs)
    if main is None:
        return False
    x = (np.arange(len(main)) + 0.5) / GRID
    ax.fill_between(x, 0, main, color=color, alpha=0.16, lw=0, zorder=2)
    for u in units_of(d, group):
        c = spine(d, group, tool, f"unit:{u}", segs)
        if c is not None:
            ax.plot(x, c, color=color, lw=0.7, alpha=0.55, zorder=2)
    ax.plot(x, main, color=color, lw=2.0, zorder=4, label=label)
    return True


def _key(ax, labels, colors):
    handles = [plt.Line2D([], [], color=c, lw=2.4, marker="s", ms=9,
                          markerfacecolor=c, markeredgecolor=c)
               for c in colors]
    ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.16),
              ncol=len(labels), frameon=False, fontsize=13)


# --------------------------------------------------------------------------- #
def fig_condition(tx_class: str, tag: str) -> None:
    """C-style panels: one curve per library, pooled over all tools."""
    segs = [KINDS[i] for i in (SEGMENT_ORDER if tx_class == "mrna" else NCR_ORDER)]
    d = load_density(tx_class, "aware")
    fig, axes = plt.subplots(1, 3, figsize=(5.6 * 3, 5.2),
                             constrained_layout=True)
    for ax, sp in zip(axes, SPECIES_ORDER):
        labels, colors = [], []
        for which, color in ((0, WT_COLOR), (1, CMP_COLOR)):
            group = COND[sp][which]
            if _curve_with_replicates(ax, d, group, POOLED, segs, color, group):
                labels.append(group)
                colors.append(color)
        guitar_panel(ax, segs, sp if tx_class == "mrna" else f"{sp} (ncRNA)")
        if labels:
            _key(ax, labels, colors)
    stem = FIG / f"FigR1_metagene_condition_{tag}"
    save(fig, str(stem))
    print("wrote", stem.name)


def fig_per_tool(tag: str, min_sites: int = 200) -> None:
    """D-style panels: one curve per tool x library."""
    segs = [KINDS[i] for i in SEGMENT_ORDER]
    d = load_density("mrna", "aware")
    counts = pd.read_csv(TAB / "merge_site_counts.tsv", sep="\t")
    counts = counts[(counts.tx_class == "mrna") & (counts.strand_mode == "aware")
                    & (counts["merge"] == "majority")
                    & (counts.tool != POOLED)]
    n_sites = {(r.group, r.tool): r.n_sites for r in counts.itertuples()}
    fig, axes = plt.subplots(1, 3, figsize=(5.6 * 3, 5.0),
                             constrained_layout=True)
    cmap = plt.get_cmap("Set1")
    handles, labels = [], []
    for ax, sp in zip(axes, SPECIES_ORDER):
        tools = [t for t in tool_order(d[(d["group"] == COND[sp][0])
                                         & (d.tool != POOLED)])
                 if n_sites.get((COND[sp][0], t), 0) >= min_sites
                 or n_sites.get((COND[sp][1], t), 0) >= min_sites]
        for i, t in enumerate(tools):
            col = cmap(i % 9)
            for which, ls in ((0, "-"), (1, "--")):
                group = COND[sp][which]
                c = spine(d, group, t, "majority", segs)
                if c is None:
                    continue
                x = (np.arange(len(c)) + 0.5) / GRID
                ax.plot(x, c, color=col, lw=1.5, ls=ls, zorder=4)
                ax.fill_between(x, 0, c, color=col, alpha=0.10, lw=0, zorder=2)
                if sp == "Arabidopsis":
                    handles.append(plt.Line2D([], [], color=col, lw=2.4, ls=ls,
                                              marker="s", ms=8, markerfacecolor=col,
                                              markeredgecolor=col))
                    labels.append(f"{t}-{'WT' if which == 0 else COND[sp][2].split()[0]}")
        guitar_panel(ax, segs, sp)
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.22),
               ncol=4, frameon=False, fontsize=12)
    stem = FIG / f"FigR2_metagene_per_tool_{tag}"
    save(fig, str(stem))
    print("wrote", stem.name)


# --------------------------------------------------------------------------- #
def fig_region_shares() -> None:
    p = pd.read_csv(TAB / "region_proportions.tsv", sep="\t")
    p = p[(p.tx_class == "mrna") & (p.strand_mode == "aware")]
    segs = [KINDS[i] for i in SEGMENT_ORDER]
    tools = tool_order(p[p.tool != POOLED])
    panels = [(COND[s][i], "WT" if i == 0 else COND[s][2])
              for s in SPECIES_ORDER for i in (0, 1)]
    fig, axes = plt.subplots(1, 6, figsize=(3.2 * 6, 5.4),
                             constrained_layout=True)
    for ax, (g, lab) in zip(axes, panels):
        sub = p[(p.group == g) & (p["merge"] == "majority") & (p.tool != POOLED)]
        mat = np.full((len(tools), len(segs)), np.nan)
        for i, t in enumerate(tools):
            for j, s in enumerate(segs):
                r = sub[(sub.tool == t) & (sub.kind == s)]
                if len(r):
                    mat[i, j] = 100 * float(r.share.iloc[0])
        ax.imshow(np.ma.masked_invalid(mat), aspect="auto", cmap="viridis",
                  vmin=0, vmax=100)
        ax.set_xticks(range(len(segs)))
        ax.set_xticklabels([SEGMENT_LABELS[s] for s in segs], fontsize=10)
        ax.set_yticks(range(len(tools)))
        ax.set_yticklabels([t if len(sub[sub.tool == t]) else f"{t}  (n/a)"
                            for t in tools], fontsize=10)
        ax.set_title(f"{g}\n({lab})", fontweight="bold", fontsize=12)
        ax.tick_params(length=0)
        for i in range(len(tools)):
            for j in range(len(segs)):
                if not np.isnan(mat[i, j]):
                    ax.text(j, i, f"{mat[i, j]:.0f}", ha="center", va="center",
                            fontsize=10, color="white" if mat[i, j] < 55 else "black")
    axes[0].set_ylabel("Tool")
    sm = plt.cm.ScalarMappable(cmap="viridis",
                               norm=mcolors.Normalize(vmin=0, vmax=100))
    fig.colorbar(sm, ax=axes[-1], fraction=0.046, pad=0.02,
                 label="% of majority sites")
    stem = FIG / "FigR3_region_shares"
    save(fig, str(stem))
    print("wrote", stem.name)


def fig_region_share_bars() -> None:
    p = pd.read_csv(TAB / "region_proportions.tsv", sep="\t")
    p = p[(p.tx_class == "mrna") & (p.strand_mode == "aware")]
    fig, axes = plt.subplots(2, 3, figsize=(5.2 * 3, 4.6 * 2),
                             constrained_layout=True)
    for j, sp in enumerate(SPECIES_ORDER):
        tools = tool_order(p[(p.group == COND[sp][0]) & (p.tool != POOLED)])
        for k, kind in enumerate(["CDS", "three_prime_UTR"]):
            name = "CDS" if kind == "CDS" else "3'UTR"
            ax = axes[k, j]
            for off, (which, col) in enumerate(((0, WT_COLOR), (1, CMP_COLOR))):
                g = COND[sp][which]
                share, sd = [], []
                for t in tools:
                    r = p[(p.group == g) & (p.tool == t) & (p.kind == kind)
                          & (p["merge"] == "majority")]
                    q = p[(p.group == g) & (p.tool == t) & (p.kind == kind)
                          & (p["merge"] == "unit_mean")]
                    share.append(100 * float(r.share.iloc[0]) if len(r) else np.nan)
                    sd.append(100 * float(q.share_sd.iloc[0])
                              if len(q) and pd.notna(q.share_sd.iloc[0]) else 0.0)
                x = np.arange(len(tools)) + (off - 0.5) * 0.38
                ax.bar(x, share, width=0.36, color=col, lw=0,
                       label="WT" if off == 0 else "perturbed")
                ax.errorbar(x, share, yerr=sd, fmt="none", ecolor="0.25",
                            elinewidth=1.0, capsize=2)
            ax.set_xticks(range(len(tools)))
            ax.set_xticklabels(tools, rotation=50, ha="right", fontsize=10)
            ax.set_ylabel(f"% of majority sites in {name}")
            ax.set_ylim(bottom=0)
            ax.set_title(f"{sp} — {name} occupancy", fontweight="bold", fontsize=13)
            if j == 0 and k == 0:
                ax.legend(frameon=False, fontsize=10)
    stem = FIG / "FigR3b_region_share_bars"
    save(fig, str(stem))
    print("wrote", stem.name)


# --------------------------------------------------------------------------- #
def fig_merge_rules(tag: str) -> None:
    segs = [KINDS[i] for i in SEGMENT_ORDER]
    d = load_density("mrna", "aware")
    fig, axes = plt.subplots(1, 3, figsize=(5.6 * 3, 4.8),
                             constrained_layout=True)
    colour = {"union": "#7a7a7a", "majority": "#1b6b3a", "intersection": "#8c3b8c"}
    for ax, sp in zip(axes, SPECIES_ORDER):
        for merge in ("union", "majority", "intersection"):
            c = spine(d, COND[sp][0], POOLED, merge, segs)
            if c is None:
                continue
            x = (np.arange(len(c)) + 0.5) / GRID
            ax.plot(x, c, lw=2.0, color=colour[merge],
                    label=f"{merge} ({consensus_word(d, COND[sp][0])})")
        guitar_panel(ax, segs, f"{sp} WT")
    axes[0].legend(frameon=False, fontsize=11)
    stem = FIG / f"FigR4_merge_rules_{tag}"
    save(fig, str(stem))
    print("wrote", stem.name)


def fig_strand_mode(tag: str) -> None:
    segs = [KINDS[i] for i in SEGMENT_ORDER]
    fig, axes = plt.subplots(1, 3, figsize=(5.6 * 3, 4.8),
                             constrained_layout=True)
    for ax, sp in zip(axes, SPECIES_ORDER):
        for mode, col, lab in (("aware", "#1b1b1b", "strand-aware (revised)"),
                               ("legacy", "#a33f3f", "strand forced to '+' (as published)")):
            c = spine(load_density("mrna", mode), COND[sp][0], POOLED,
                      "majority", segs)
            if c is None:
                continue
            x = (np.arange(len(c)) + 0.5) / GRID
            ax.plot(x, c, lw=2.0, color=col, label=lab)
        guitar_panel(ax, segs, f"{sp} WT")
    axes[0].legend(frameon=False, fontsize=10)
    stem = FIG / f"FigR5_strand_mode_{tag}"
    save(fig, str(stem))
    print("wrote", stem.name)


# --------------------------------------------------------------------------- #
def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig_condition("mrna", "mrna")
    fig_condition("ncrna", "ncrna")
    fig_per_tool("mrna")
    fig_region_shares()
    fig_region_share_bars()
    fig_merge_rules("mrna")
    fig_strand_mode("mrna")


if __name__ == "__main__":
    main()
