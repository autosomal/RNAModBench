#!/usr/bin/env python
"""R3-2 / E6: what the replicate structure does to the benchmark.

Figures in ``$RNAMODBENCH_LOCAL/replicate_structure/figures``:

  FigR6_replicate_merge        sites kept by each merge rule, per tool and group
  FigR6b_merge_summaries       unit-to-unit Jaccard, and the recall / precision
                               shift from replicate mean to consensus, both
                               recomputed on ONE shared measurable universe
                               (common/consensus_eval.py) -- the evaluation
                               table scores replicates on their own universes
                               and merged sets on the common one, so its two
                               rows are not comparable
  FigR7_per_replicate          per-replicate precision / recall with the
                               across-replicate SD (forest plot)
  FigR8_merge_vs_replicate_mean  merged call set vs the mean of its replicates:
                               does the merge rule move the numbers or only the
                               noise?

Mouse groups are two independent studies, not replicates: they are drawn as
paired points and never summarised as mean +/- SD.
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
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from common.config import SITES_ROOT                      # noqa: E402
from common import consensus_eval as ev                   # noqa: E402
from common.figstyle import apply as apply_style          # noqa: E402
from common.figstyle import save                          # noqa: E402

EVAL = (_RB / "data/evaluation/tables")
GUITAR = (_XB / "harmonisation/guitar_metagene/tables")
OUT = (_RB / "analysis/replicate_structure")
TAB, FIG = OUT / "tables", OUT / "figures"
apply_style()

GROUPS = ["Arabidopsis_WT", "Arabidopsis_KD", "Mouse_WT", "Mouse_KO",
          "HeLa_WT", "HeLa_IVT"]
SPECIES_GROUPS = {"Arabidopsis": ("Arabidopsis_WT", "Arabidopsis_KD"),
                  "Mouse": ("Mouse_WT", "Mouse_KO"),
                  "Human": ("HeLa_WT", "HeLa_IVT")}
GROUP_SPECIES = {"Arabidopsis_WT": "Arabidopsis", "Arabidopsis_KD": "Arabidopsis",
                 "Mouse_WT": "Mouse", "Mouse_KO": "Mouse",
                 "HeLa_WT": "Human", "HeLa_IVT": "Human"}
STUDY_GROUPS = {"Mouse_WT", "Mouse_KO"}
MERGES = ("union", "majority", "intersection")
MERGE_COLOR = {"union": "#8c8c8c", "majority": "#1f5c8b", "intersection": "#e08a2e"}
TOOL_ORDER = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
              "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
              "NanoSPA_m6A", "xPore", "yanocomp"]


def tool_order(df: pd.DataFrame) -> list[str]:
    have = [t for t in TOOL_ORDER if t in set(df.tool)]
    return have + sorted(set(df.tool) - set(TOOL_ORDER))


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def _split(s):
    return [v for v in (_num(p) for p in str(s).split("|")) if not np.isnan(v)] \
        if isinstance(s, str) and s.strip() else []


# --------------------------------------------------------------------------- #
def fig_merge_effect() -> None:
    counts = pd.read_csv(GUITAR / "merge_site_counts.tsv", sep="\t")
    counts = counts[(counts.tx_class == "mrna") & (counts.strand_mode == "aware")]

    fig, axes = plt.subplots(3, 2, figsize=(7.2 * 2, 4.4 * 3),
                             constrained_layout=True)
    for ax, g in zip(axes.ravel(), GROUPS):
        sub = counts[counts.group == g]
        tools = tool_order(sub)
        width = 0.26
        for i, merge in enumerate(MERGES):
            y = [_first(sub[(sub.tool == t) & (sub["merge"] == merge)].n_sites)
                 for t in tools]
            ax.bar(np.arange(len(tools)) + (i - 1) * width, y, width=width,
                   color=MERGE_COLOR[merge], lw=0)
        ax.set_xticks(range(len(tools)))
        ax.set_xticklabels(tools, rotation=40, ha="right", fontsize=10)
        ax.set_ylabel("Annotated sites")
        ax.set_ylim(bottom=0)
        basis = sub.merge_basis.iloc[0] if len(sub) else "replicates"
        ax.set_title(f"{g} — sites kept by merge rule ({basis})",
                     fontweight="bold", fontsize=12)
    handles = [Patch(color=MERGE_COLOR[m], label=f"{m} consensus")
               for m in MERGES]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=11)
    stem = FIG / "FigR6_replicate_merge"
    save(fig, str(stem))
    print("wrote", stem.name)

    fig, axes = plt.subplots(1, 3, figsize=(6.2 * 3, 5.4),
                             constrained_layout=True)
    _jaccard_heat(axes[0], counts)
    uni = common_universe_table()
    _delta_panel(axes[1], uni, "recall",
                 "Recall on one shared universe")
    _delta_panel(axes[2], uni, "precision",
                 "Precision on one shared universe")
    stem = FIG / "FigR6b_merge_summaries"
    save(fig, str(stem))
    print("wrote", stem.name)


def _first(s):
    return float(s.iloc[0]) if len(s) else np.nan


def _jaccard_heat(ax, counts: pd.DataFrame) -> None:
    tools = tool_order(counts)
    mat = np.full((len(tools), len(GROUPS)), np.nan)
    for i, t in enumerate(tools):
        for j, g in enumerate(GROUPS):
            r = counts[(counts.tool == t) & (counts.group == g)
                       & (counts["merge"] == "union")]
            if len(r) and pd.notna(r.mean_pairwise_jaccard.iloc[0]):
                mat[i, j] = float(r.mean_pairwise_jaccard.iloc[0])
    im = ax.imshow(np.ma.masked_invalid(mat), aspect="auto", cmap="magma",
                   vmin=0, vmax=1)
    ax.set_xticks(range(len(GROUPS)))
    ax.set_xticklabels([g.replace("_", "\n") for g in GROUPS], fontsize=10)
    ax.set_yticks(range(len(tools)))
    ax.set_yticklabels(tools, fontsize=10)
    ax.tick_params(length=0)
    for i in range(len(tools)):
        for j in range(len(GROUPS)):
            if not np.isnan(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                        fontsize=10, color="k" if mat[i, j] > 0.5 else "w")
    ax.set_title("Mean pairwise Jaccard between units", fontweight="bold", fontsize=11)
    fig = ax.figure
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="Jaccard")


def common_universe_table() -> pd.DataFrame:
    """Replicates and majority consensus, all scored on the group's shared universe."""
    cache = TAB / "scoring_on_common_universe.tsv"
    if cache.exists() and "--recompute" not in sys.argv:
        return pd.read_csv(cache, sep="\t")
    rows = []
    for group, species in GROUP_SPECIES.items():
        ref = ev.reference(species)
        universe = ev.common_universe("RNA002", species, group)
        if not universe:
            continue
        ref_u = ev.reference_positions_in_universe(ref, universe)
        for tool, d in ev.consensus_sets("RNA002", species, group).items():
            for tag, sites in ev.replicate_sets("RNA002", species, group, tool).items():
                rows.append(dict(group=group, tool=tool, unit=tag,
                                 **ev.score(sites, universe, ref, ref_u)))
            rows.append(dict(group=group, tool=tool, unit="majority",
                             **ev.score(d["sites"], universe, ref, ref_u)))
        print(f"  {group}: universe {len(universe):,} positions, "
              f"{len(ref_u):,} reference sites", flush=True)
    tab = pd.DataFrame(rows)
    TAB.mkdir(parents=True, exist_ok=True)
    tab.to_csv(cache, sep="\t", index=False)
    return tab


def _delta_panel(ax, tab: pd.DataFrame, metric: str, title: str) -> None:
    per = tab[tab.unit != "majority"].groupby(["group", "tool"])[metric].mean()
    maj = tab[tab.unit == "majority"].set_index(["group", "tool"])[metric]
    both = pd.concat([per.rename("replicate_mean"), maj.rename("majority")], axis=1).dropna()
    for i, g in enumerate(GROUPS):
        s = both.xs(g, level="group") if g in both.index.get_level_values(0) else None
        if s is None or not len(s):
            continue
        d = s.majority - s.replicate_mean
        ax.plot([i] * len(s), d, ".", color="0.35", ms=5, zorder=3)
        ax.bar(i, d.max() - d.min(), bottom=d.min(), width=0.55,
               color="#1f5c8b", alpha=.3, lw=0)
    ax.axhline(0, color="k", lw=1)
    ax.set_xticks(range(len(GROUPS)))
    ax.set_xticklabels([g.replace("_", "\n") for g in GROUPS], fontsize=10)
    ax.set_ylabel(f"Δ {metric} (majority − replicate mean)")
    ax.set_title(title, fontweight="bold", fontsize=12)


# --------------------------------------------------------------------------- #
def fig_per_replicate() -> None:
    fr = pd.read_csv(EVAL / "figure_ready_replicates.tsv", sep="\t")
    fr = fr[(fr.window == 2) & (fr.min_cov == 10)]
    fig, axes = plt.subplots(2, 3, figsize=(5.4 * 3, 4.6 * 2),
                             constrained_layout=True)
    for k, (metric, cap) in enumerate([("recall", "Recall"),
                                       ("precision", "Precision")]):
        for j, (sp, pair) in enumerate(SPECIES_GROUPS.items()):
            ax = axes[k, j]
            tools = tool_order(fr[fr.dataset_group.isin(pair)])
            for off, g in enumerate(pair):
                sub = fr[fr.dataset_group == g]
                col = ("#a33f3f" if g in STUDY_GROUPS
                       else ("#1f5c8b" if g.endswith("WT") else "#e08a2e"))
                for i, t in enumerate(tools):
                    r = sub[sub.tool == t]
                    if not len(r):
                        continue
                    vals = _split(r[f"{metric}_each"].iloc[0])
                    mean = _num(r[f"{metric}_mean"].iloc[0])
                    sd = _num(r[f"{metric}_sd"].iloc[0])
                    x = i + (off - 0.5) * 0.32
                    if vals:
                        ax.plot([x] * len(vals), vals, "o", ms=4.5, color=col,
                                alpha=.85, zorder=3)
                    if g in STUDY_GROUPS and len(vals) == 2:
                        ax.plot([x, x], vals, "-", color="0.4", lw=1.0)
                    if not np.isnan(mean):
                        ax.plot([x - .11, x + .11], [mean, mean], color=col,
                                lw=2.2, zorder=4)
                        if not np.isnan(sd) and g not in STUDY_GROUPS:
                            ax.plot([x, x], [mean - sd, mean + sd], color=col,
                                    lw=1.1)
            ax.set_xticks(range(len(tools)))
            ax.set_xticklabels(tools, rotation=55, ha="right", fontsize=10)
            ax.set_ylabel(f"{cap} (w = 2 bp)")
            ax.set_ylim(bottom=0)
            ax.set_title(f"{sp} — " + " vs ".join(g.split("_", 1)[1] for g in pair)
                         + ("  (two studies)" if sp == "Mouse" else ""),
                         fontweight="bold", fontsize=11)
    fig.legend(handles=[
        plt.Line2D([], [], marker="o", ls="", ms=6, color="#1f5c8b"),
        plt.Line2D([], [], color="#1f5c8b", lw=2.2),
        plt.Line2D([], [], color="#1f5c8b", lw=1.1),
        plt.Line2D([], [], marker="o", ls="", ms=6, color="#a33f3f")],
        labels=["individual replicates", "mean", "+/- SD",
                "mouse: two studies, never pooled"],
        loc="center left", bbox_to_anchor=(1.005, 0.5), frameon=False, fontsize=10)
    stem = FIG / "FigR7_per_replicate"
    save(fig, str(stem))
    print("wrote", stem.name)


# --------------------------------------------------------------------------- #
def fig_merge_vs_mean() -> None:
    uni = common_universe_table()
    per = (uni[uni.unit != "majority"]
           .groupby(["group", "tool"])[["recall", "precision"]].mean()
           .add_suffix("_replicate_mean"))
    maj = (uni[uni.unit == "majority"]
           .set_index(["group", "tool"])[["recall", "precision"]]
           .add_suffix("_majority"))
    tab = pd.concat([per, maj], axis=1).dropna().reset_index()
    tab.to_csv(TAB / "merge_vs_replicate_mean.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(6.2 * 2, 5.0),
                             constrained_layout=True)
    for ax, metric in zip(axes, ["recall", "precision"]):
        for g in GROUPS:
            s = tab[tab.group == g]
            if not len(s):
                continue
            ax.scatter(s[f"{metric}_replicate_mean"], s[f"{metric}_majority"],
                       s=34, alpha=.85, label=g)
        rho = stats.spearmanr(tab[f"{metric}_replicate_mean"],
                              tab[f"{metric}_majority"])
        lim = float(max(ax.get_xlim()[1], ax.get_ylim()[1]))
        ax.plot([0, lim], [0, lim], ls=(0, (3, 3)), color="0.4", lw=1)
        ax.set_xlabel("mean of the individual replicates")
        ax.set_ylabel("majority consensus, same universe")
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_title(f"{metric}: Spearman $\\rho$ = {rho.statistic:.2f} "
                     f"(n = {len(tab)})", fontweight="bold", fontsize=12)
    axes[0].legend(frameon=False, fontsize=10)
    stem = FIG / "FigR8_merge_vs_replicate_mean"
    save(fig, str(stem))
    print("wrote", stem.name)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    fig_merge_effect()
    fig_per_replicate()
    fig_merge_vs_mean()


if __name__ == "__main__":
    main()
