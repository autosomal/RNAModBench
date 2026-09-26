#!/usr/bin/env python
"""R3-6: is the metagene shape a biological pattern or a motif-detection artefact?

The reviewer concern is that a tool's apparent preference for the 3'UTR / stop
codon region may simply be its DRACH preference showing through, since DRACH
density is itself non-uniform along the transcript.  The metagene site table
from stage 1 carries, for every site, the DRACH status of its *strand-aware*
5-mer, so the two can be separated directly.

Figures in ``$RNAMODBENCH_LOCAL/motif_metagene_bias/figures``:

  FigR9_drach_stratified_metagene  majority-consensus metagene split into DRACH
                                   and non-DRACH sites, per species and condition
  FigR10_motif_vs_shape            per-tool DRACH share vs 3'UTR share, and the
                                   GLORI hit rate of DRACH vs non-DRACH sites

Both use the majority-consensus call set (a site must be seen by strictly more
than half of the replicates), so a site that only one replicate happened to
call cannot drive the shape.
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
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from common.config import SITES_ROOT                        # noqa: E402
from common.figstyle import apply as apply_style            # noqa: E402
from common.figstyle import guitar_panel, save               # noqa: E402
from common.regionmodel import KINDS, SEGMENT_ORDER          # noqa: E402

GUITAR = (_XB / "harmonisation/guitar_metagene/tables")
OUT = (_RB / "analysis/motif_metagene_bias")
TAB, FIG = OUT / "tables", OUT / "figures"
apply_style()

GRID, FINE = 200, 800
BW = 0.05
SPECIES_ORDER = ["Arabidopsis", "Mouse", "Human"]
COND = {"Arabidopsis": ("Arabidopsis_WT", "Arabidopsis_KD", "fip37 KD"),
        "Mouse": ("Mouse_WT", "Mouse_KO", "Mettl3 KO"),
        "Human": ("HeLa_WT", "HeLa_IVT", "IVT")}
SEGS = [KINDS[i] for i in SEGMENT_ORDER]
DRACH_COLOR, NON_COLOR = "#1f5c8b", "#a33f3f"
TOOL_ORDER = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
              "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
              "NanoSPA_m6A", "xPore", "yanocomp"]


def density(x: np.ndarray) -> np.ndarray:
    from scipy.ndimage import gaussian_filter1d
    if x.size == 0:
        return np.zeros(GRID)
    h, _ = np.histogram(x, bins=FINE, range=(0.0, 1.0))
    s = gaussian_filter1d(h.astype(float), BW * FINE, mode="reflect")
    s = s.reshape(GRID, FINE // GRID).sum(axis=1)
    return s / s.sum() * GRID if s.sum() else np.zeros(GRID)


def load() -> pd.DataFrame:
    s = pd.read_csv(GUITAR / "metagene_sites.tsv.gz", sep="\t")
    s["kind_name"] = s["kind"].map(KINDS)
    s = s[(s.tx_class == "mrna") & (s.strand_mode == "aware")
          & s.kind_name.isin(SEGS)]
    s["majority"] = s.n_units >= (s.n_units_group // 2 + 1)   # common.consensus.quorum
    return s[s.majority & s.is_drach.notna()]


def spine(sub: pd.DataFrame) -> np.ndarray:
    """Share-weighted metagene of one stratum, pooled over tools."""
    tot = len(sub)
    out = np.zeros(len(SEGS) * GRID)
    if not tot:
        return out
    for i, name in enumerate(SEGS):
        v = sub.loc[sub.kind_name == name, "norm_pos"].to_numpy()
        if v.size:
            out[i * GRID:(i + 1) * GRID] = density(v) * (v.size / tot)
    return out


# --------------------------------------------------------------------------- #
def fig_stratified(sites: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(5.6 * 3, 5.0 * 2),
                             constrained_layout=True)
    for j, sp in enumerate(SPECIES_ORDER):
        for k, (which, lab) in enumerate(((0, "WT"), (1, COND[sp][2]))):
            group = COND[sp][which]
            sub = sites[sites.group == group]
            ax = axes[k, j]
            x = (np.arange(len(SEGS) * GRID) + 0.5) / GRID
            for mask, col, name in ((sub.is_drach == True, DRACH_COLOR, "DRACH"),
                                    (sub.is_drach == False, NON_COLOR, "non-DRACH")):
                s = sub[mask]
                if not len(s):
                    continue
                c = spine(s)
                ax.fill_between(x, 0, c, color=col, alpha=0.14, lw=0, zorder=2)
                ax.plot(x, c, color=col, lw=2.0, zorder=4,
                        label=f"{name} ({100 * len(s) / len(sub):.0f}% of sites, "
                              f"n = {len(s):,})")
            guitar_panel(ax, SEGS, f"{sp} — {lab}")
            if k == 0 and j == 0:
                ax.legend(frameon=False, fontsize=10, loc="upper left")
    stem = FIG / "FigR9_drach_stratified_metagene"
    save(fig, str(stem))
    print("wrote", stem.name)


# --------------------------------------------------------------------------- #
def fig_shape_vs_motif(sites: pd.DataFrame) -> None:
    rows = []
    for (group, tool), sub in sites.groupby(["group", "tool"]):
        tot = len(sub)
        if tot < 50:
            continue
        utr3 = (sub.kind_name == "three_prime_UTR").sum() / tot
        cds = (sub.kind_name == "CDS").sum() / tot
        drach = float(sub.is_drach.mean())
        hit2_d = float((sub[sub.is_drach == True].dist_glori <= 2).mean())
        hit2_n = float((sub[sub.is_drach == False].dist_glori <= 2).mean())
        rows.append(dict(group=group, species=sub.species.iloc[0], tool=tool,
                         n_sites=tot, drach_share=round(drach, 4),
                         utr3_share=round(utr3, 4), cds_share=round(cds, 4),
                         glori_hit_w2=round(float((sub.dist_glori <= 2).mean()), 4),
                         glori_hit_w2_drach=round(hit2_d, 4),
                         glori_hit_w2_non_drach=round(hit2_n, 4),
                         ratio_drach_over_non=(round(hit2_d / hit2_n, 2)
                                               if hit2_n else np.nan)))
    tab = pd.DataFrame(rows).sort_values(["species", "group", "tool"])
    tab.to_csv(TAB / "shape_vs_motif.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(4.9 * 4, 4.4 * 2), constrained_layout=True)
    gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 1.15])
    for j, sp in enumerate(SPECIES_ORDER):
        ax = fig.add_subplot(gs[0, j])
        _share_bars(ax, tab, COND[sp][0], f"{sp} — WT", key=(j == 0))
        ax = fig.add_subplot(gs[1, j])
        _share_bars(ax, tab, COND[sp][1], f"{sp} — {COND[sp][2]}")
    ax = fig.add_subplot(gs[:, 3])
    _glori_support(ax, tab)
    stem = FIG / "FigR10_motif_vs_shape"
    save(fig, str(stem))
    print("wrote", stem.name)


def _share_bars(ax, tab: pd.DataFrame, group: str, title: str,
                key: bool = False) -> None:
    s = tab[tab.group == group].sort_values("drach_share")
    if s.empty:
        ax.axis("off")
        return
    y = np.arange(len(s))
    ax.barh(y - 0.2, s.drach_share, height=0.38, color=DRACH_COLOR, lw=0,
            label="DRACH share")
    ax.barh(y + 0.2, s.utr3_share, height=0.38, color=NON_COLOR, lw=0,
            label="3'UTR share")
    if key:
        ax.legend(frameon=False, fontsize=10, loc="lower right")
    rho, p = stats.spearmanr(s.drach_share, s.utr3_share)
    ax.set_yticks(y)
    ax.set_yticklabels(s.tool, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("share of majority sites")
    ax.set_title(f"{title}\nSpearman $\\rho$ = {rho:.2f} (p = {p:.1g})",
                 fontweight="bold", fontsize=12)


def _glori_support(ax, tab: pd.DataFrame) -> None:
    ok = tab[(tab.glori_hit_w2_non_drach > 0)
             & tab.group.isin([COND[s][0] for s in SPECIES_ORDER])]
    ok = ok.sort_values("ratio_drach_over_non", ascending=False)
    y = np.arange(len(ok))
    ax.barh(y - 0.2, ok.glori_hit_w2_drach, height=0.38, color=DRACH_COLOR,
            label="DRACH sites")
    ax.barh(y + 0.2, ok.glori_hit_w2_non_drach, height=0.38, color=NON_COLOR,
            label="non-DRACH sites")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.group.split('_')[0]} {r.tool}" for r in ok.itertuples()],
                       fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel("fraction of majority sites within 2 bp of GLORI")
    ax.set_title("Reference support by motif class (WT libraries)",
                 fontweight="bold", fontsize=12)
    ax.legend(frameon=False, fontsize=10, loc="lower right")


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    sites = load()
    print(f"{len(sites):,} majority DRACH-known sites")
    fig_stratified(sites)
    fig_shape_vs_motif(sites)


if __name__ == "__main__":
    main()
