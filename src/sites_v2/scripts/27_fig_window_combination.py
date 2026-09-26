#!/usr/bin/env python
"""NA-1 (window inflation) and NA-7 (multi-tool combination), replicate-aware.

R1-3 / R3-3 / E2 objected that a +/-50 bp matching window mechanically inflates
the GLORI overlap.  Panel (a) and (b) quantify exactly how much of the gain is
localisation tolerance rather than detection.

R3-4 / R3-m2 asked whether combining tools really helps.  For that the
combination must be scored on the SAME measurable universe as the single tools,
so panel (c) rebuilds the call sets from ``sites_v2``:

  universe        positions measurable in every replicate of the group
                  (exonic, reference-base compatible, coverage >= 10)
  tool call set   majority consensus over the replicates (strictly more than half)
  combination     union / intersection of those consensus sets
  scoring         distance to the nearest GLORI site, window w = 2 bp

Figures in ``04_revision_analysis/window_and_combination/figures``.
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

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from common import consensus_eval as ev                   # noqa: E402
from common.config import SITES_ROOT                      # noqa: E402
from common.figstyle import apply as apply_style          # noqa: E402
from common.figstyle import save                          # noqa: E402

EVAL = (_RB / "data/evaluation/tables")
OUT = (_RB / "analysis/window_and_combination")
TAB, FIG = OUT / "tables", OUT / "figures"
apply_style()

WINDOW = 2
GROUPS = {"Arabidopsis": ["Arabidopsis_WT", "Arabidopsis_KD"],
          "Mouse": ["Mouse_WT", "Mouse_KO"],
          "Human": ["HeLa_WT", "HeLa_IVT"]}
TOOL_ORDER = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
              "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
              "NanoSPA_m6A", "xPore", "yanocomp"]


# --------------------------------------------------------------------------- #
def fig_window(curve: pd.DataFrame | None = None) -> pd.DataFrame:
    d = curve if curve is not None else pd.read_csv(
        EVAL / "m6a_localization_curve.tsv", sep="\t")
    d = d[d.platform == "RNA002"]
    per = (d.groupby(["species", "tool", "window"])
             .agg(hit_rate=("hit_rate", "mean"),
                  exact=("localization_accuracy", "mean"),
                  recall=("recall", "mean")).reset_index())
    per.to_csv(TAB / "window_sweep_mean_per_tool.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(5.0 * 3, 4.2), constrained_layout=True)
    for ax, sp in zip(axes, GROUPS):
        s = per[per.species == sp]
        w = np.sort(s.window.unique())
        hit = [s[s.window == x].hit_rate.mean() for x in w]
        lo = [s[s.window == x].hit_rate.min() for x in w]
        hi = [s[s.window == x].hit_rate.max() for x in w]
        ex = [s[s.window == x].exact.mean() for x in w]
        ax.fill_between(w, lo, hi, color="#1f5c8b", alpha=.18, lw=0,
                        label="range across tools")
        ax.plot(w, hit, "-o", color="#1f5c8b", ms=4, lw=1.8,
                label="mean GLORI hit rate")
        ax.plot(w, ex, "--", color="#a33f3f", lw=1.8,
                label="exact single-nucleotide localisation")
        ax.axvline(WINDOW, color="0.35", lw=1.0, ls=(0, (3, 2)))
        ax.set_xticks([int(x) for x in w if x in (0, 2, 5, 10, 20, 50)])
        ax.set_xlabel("matching window (bp)")
        ax.set_ylabel("fraction")
        ax.set_title(f"{sp}  (revised window = {WINDOW} bp)", fontweight="bold")
        ax.set_ylim(0, 1)
    axes[0].legend(frameon=False, fontsize=10)
    stem = FIG / "FigR13_window_inflation"
    save(fig, str(stem))
    print("wrote", stem.name)
    return per


# --------------------------------------------------------------------------- #
def fig_combination(recompute: bool = False) -> None:
    cached = TAB / "combination_eval.tsv"
    if not recompute and cached.exists():
        tab = pd.read_csv(cached, sep="\t")
        print(f"reusing {cached.name} ({len(tab)} rows); "
              f"pass --recompute to rebuild it")
    else:
        tab = _score_combinations()
        tab.to_csv(cached, sep="\t", index=False)
    return _plot_combinations(tab)


def _score_combinations() -> pd.DataFrame:
    rows = []
    for species, groups in GROUPS.items():
        ref = ev.reference(species)
        for group in groups:
            sets = ev.consensus_sets("RNA002", species, group)
            consensus = {t: d["sites"] for t, d in sets.items()}
            n_units = max([d["n_units"] for d in sets.values()] or [0])
            universe = ev.common_universe("RNA002", species, group)
            if not universe:
                continue
            ref_in_u = ev.reference_positions_in_universe(ref, universe)
            tools = [t for t in TOOL_ORDER if t in consensus and consensus[t]]
            single = {t: ev.score(consensus[t], universe, ref, ref_in_u, WINDOW)
                      for t in tools}
            for t, s in single.items():
                rows.append(dict(species=species, group=group, n_units=n_units,
                                 combination=t, k=1, strategy="single", **s))
            order = sorted(tools, key=lambda t: -single[t]["tp"])
            for k in range(2, len(order) + 1):
                sel = order[:k]
                rows.append(dict(species=species, group=group, n_units=n_units,
                                 combination="+".join(sel), k=k,
                                 strategy="greedy union (by TP)",
                                 **ev.score(set().union(*[consensus[t] for t in sel]),
                                            universe, ref, ref_in_u, WINDOW)))
                rows.append(dict(species=species, group=group, n_units=n_units,
                                 combination="+".join(sel), k=k,
                                 strategy="intersection",
                                 **ev.score(set.intersection(*[consensus[t] for t in sel]),
                                            universe, ref, ref_in_u, WINDOW)))
    return pd.DataFrame(rows)


def _plot_combinations(tab: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(5.2 * 3, 4.4 * 2),
                             constrained_layout=True)
    strategies = (("single", "0.45", "o"), ("greedy union (by TP)", "#1f5c8b", "s"),
                  ("intersection", "#e08a2e", "^"))
    for j, (species, groups) in enumerate(GROUPS.items()):
        for i, metric in enumerate(["recall", "precision"]):
            ax = axes[i, j]
            for gi, group in enumerate(groups):
                s = tab[tab.group == group]
                if s.empty:
                    continue
                ls = "-" if gi == 0 else ":"
                best = s[s.k == 1].sort_values("tp").iloc[-1]
                for strat, col, mk in strategies:
                    q = s[s.strategy == strat].sort_values("k")
                    ax.plot(q.k, q[metric], marker=mk, ms=5, lw=1.6, ls=ls,
                            color=col)
                ax.axhline(float(best[metric]), color="0.2", lw=1.0,
                           ls=(0, (3, 2)))
            ax.set_xlabel("number of tools combined")
            ax.set_ylabel(f"{metric} (w = {WINDOW} bp)")
            ax.set_ylim(0, 1.02)
            ax.set_xlim(0.5, 13.5)
            ax.set_title(f"{species} — {metric}", fontweight="bold", fontsize=11)
    handles = [plt.Line2D([], [], color=c, marker=m, ms=5, lw=1.6, label=n)
               for n, c, m in strategies]
    handles += [plt.Line2D([], [], color="0.2", lw=1.0, ls=(0, (3, 2)),
                           label="best single tool"),
                plt.Line2D([], [], color="0.2", lw=1.6, ls="-", label="WT library"),
                plt.Line2D([], [], color="0.2", lw=1.6, ls=":",
                           label="KD / KO / IVT library")]
    fig.legend(handles=handles, frameon=False, fontsize=10,
               loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=6)
    stem = FIG / "FigR14_tool_combination"
    save(fig, str(stem))
    print("wrote", stem.name)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    fig_window()
    fig_combination(recompute="--recompute" in sys.argv)


if __name__ == "__main__":
    main()
