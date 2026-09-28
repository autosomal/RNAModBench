#!/usr/bin/env python
"""36 -- main-text Figure 6, revision rebuild (R3-m2 / R3-4 / R1-2 / R3-2).

What the published Figure 6 was
-------------------------------
Drawn by ``the original submission's figure code/Figure6/NGS_tool_combinations.ipynb``
(three separately saved PDFs, assembled by hand) from the legacy
``output/<group>/Tools.txt`` aggregate of ONE sample per species:

* panel A ``greedy_tool_selection()`` -- greedy max-coverage over the *raw GLORI
  list*, no candidate universe, no precision constraint;
* panel B intersection precision ``|I n GLORI| / |I|`` -- again no universe, so
  FP/TN had no boundary and calls outside the measurable set counted nowhere;
* panel C -- 78 naive two-tool precisions;
* the same notebook produced supplementary S5 ("GLORI Coverage vs Tool Count").

What this script does instead
-----------------------------
Everything is scored on the harmonisation evaluation layer (``README.md`` section 3):

    universe U   exonic, reference-base-compatible, coverage >= 10 in EVERY
                 independent unit of the group (``consensus_eval.common_universe``)
    call sets    one per independent sequencing unit (never a merged aggregate)
    scoring      TP = call within ``WINDOW`` = 2 bp of a GLORI site inside U,
                 precision = TP / n_calls, recall = TP / |reference in U|

Selection criterion (stated in the Methods and the figure caption):
for k = 1..5 ALL C(n, k) combinations are enumerated; the reported combination
maximises the **mean recall over the group's independent units subject to mean
precision >= the chance level** ``p0 = |G n U| / |U|`` (the precision a random
call set would reach; the negative-control FP density is reported next to it,
panel C1).  The greedy forward choice (add the tool with the largest marginal
mean-recall gain) is written next to it and compared.  Replicate structure is
never collapsed: Arabidopsis / HeLa carry n = 3 biological replicates
(mean +- SD), the two mouse studies are drawn as two separate curves and are
never averaged.

Figures
-------
The revision **Figure 6** is produced by ``65_fig6_combination_page.py`` and
**S5** by ``62_figS5_figure.py``; both read the frozen tables below and print at
their final size (1:1, no scaling).  The superseded 4-row main layout
(``fig_main``) and the superseded 3-row S5 layout (``fig_supp``) remain
reproducible with ``--fig6-legacy`` / ``--s5-legacy`` but are **off by default**:
they were drawn on a 380 x 367 mm canvas whose 8-14 pt text printed at 3.7-6.2 pt
at the 0.95\textwidth placement of the manuscript, which is why they were
replaced (2026-09-21).

Evidence tables -> ``figures/figure6/tables/``
    fig6_combination_selected.tsv        selected combination per k and group
    fig6_per_unit_selected.tsv           per-unit values of the selected k = 1 / 5
    fig6_per_unit_allk.tsv               per-unit values of the selected k = 1..5
    fig6_selected_members.tsv            which of the 13 tools each selected k uses
    fig6_single_tool_metrics.tsv         every single tool (k = 1) in the plane
    fig6_negative_control_fp.tsv         FP/10 kb of the selected k = 1 / 2 / 5
    fig6_negative_control_fp_bytool.tsv  FP/10 kb of every single tool (controls)
    figS5_search_space.tsv               **every** enumerated k = 1..5 combination
                                         with its mean recall/precision and the
                                         ``feasible`` flag (mean union PPV >=
                                         chance level ``p0``) that *is* the
                                         selection criterion, plus the per-k
                                         ``selected`` optimum
    figS5_greedy_1to13.tsv               greedy forward path, 1 -> 13 tools
    figS5_two_tool_pairs.tsv             every C(13,2) pair

Outputs -> ``figures/figure6/{tables,figures,logs}``

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    figures/figure6/src/40_fig6_combination.py
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
import itertools
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = _RB / "src/harmonisation"
sys.path.insert(0, str(HERE))

from common import consensus_eval as ev                       # noqa: E402
from common.config import CALLSET_ROOT, SITES_ROOT, TABLE_DIR  # noqa: E402
from common.figstyle import apply as apply_style               # noqa: E402
from common.figstyle import save                               # noqa: E402
from common.io_utils import write_table                        # noqa: E402
from common.manifest import Inventory, setup_logger            # noqa: E402

OUT = (_RB / "figures/figure6")
TAB, FIG, LOG = OUT / "tables", OUT / "figures", OUT / "logs"

WINDOW = 2
KMAX = 5
KMAX_SUPP = 13
PLATFORM = "RNA002"
MOD = "m6A"

TOOLS = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
         "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
         "NanoSPA_m6A", "xPore", "yanocomp"]

#: species -> (dataset group, how the units are drawn)
#: "replicates" = mean +- SD over independent units; "studies" = one curve per
#: cross-study unit, never averaged (harmonisation/README.md section 5).
GROUPS: dict[str, tuple[str, str]] = {
    "Arabidopsis": ("Arabidopsis_WT", "replicates"),
    "Mouse": ("Mouse_WT", "studies"),
    "Human": ("HeLa_WT", "replicates"),
}
COLOUR = {"Arabidopsis": "#1E888B", "Mouse": "#F5B264", "Human": "#3778A0"}
UNION_C, ISECT_C = "#1f5c8b", "#e08a2e"
#: unmodified controls for panel C1: label -> (species, dataset group)
CONTROLS = {"Curlcake IVT": ("Curlcake", "Curlcake_IVT"),
            "HeLa IVT": ("Human", "HeLa_IVT")}
#: column order of ``figS5_search_space.tsv`` (every k = 1..5 combination)
SEARCH_COLUMNS = ["species", "group", "k", "combination", "feasible", "selected",
                  "union_recall_mean", "union_recall_sd",
                  "union_precision_mean", "union_precision_sd",
                  "isect_recall_mean", "isect_recall_sd",
                  "isect_precision_mean", "isect_precision_sd",
                  "union_n_mean", "union_tp_mean",
                  "isect_n_mean", "isect_tp_mean",
                  "n_units", "p0_chance_precision"]
#: ``fig6_single_tool_metrics.tsv`` -- the k = 1 plane (all 13 configurations)
SINGLE_COLUMNS = ["species", "group", "tool", "n_units",
                  "union_recall_mean", "union_recall_sd",
                  "union_precision_mean", "union_precision_sd",
                  "union_n_mean", "union_tp_mean",
                  "p0_chance_precision", "selected"]
#: ``fig6_selected_members.tsv`` -- membership matrix of the selected combinations
MEMBER_COLUMNS = ["species", "group", "k", "tool", "in_combination"]
#: ``fig6_per_unit_allk.tsv`` -- per-unit values of the selected k = 1..5
#: (``*_recall_all`` = the same numerator against every GLORI site of the group,
#: i.e. without the measurable-universe restriction; see ``_combo_one_unit``)
ALLK_COLUMNS = ["species", "group", "unit", "combination", "k", "union_n",
                "union_tp", "union_precision", "union_recall",
                "union_recall_all", "isect_n", "isect_tp", "isect_precision",
                "isect_recall", "isect_recall_all"]
#: ``fig6_negative_control_fp_bytool.tsv`` -- single-tool FP on the controls
CTL_TOOL_COLUMNS = ["control", "sample", "sequencing_unit", "tool", "n_calls",
                    "region_bp", "fp_per_10kb", "tools_available"]


# --------------------------------------------------------------------------- #
# scoring primitives
# --------------------------------------------------------------------------- #
def _hit_mask(pos: np.ndarray, ref_arr: np.ndarray, window: int) -> np.ndarray:
    """True where a call sits within ``window`` bp of a reference position."""
    if ref_arr.size == 0 or pos.size == 0:
        return np.zeros(pos.size, dtype=bool)
    i = np.searchsorted(ref_arr, pos)
    lo = np.clip(i - 1, 0, ref_arr.size - 1)
    hi = np.clip(i, 0, ref_arr.size - 1)
    dist = np.minimum(np.abs(ref_arr[lo] - pos), np.abs(ref_arr[hi] - pos))
    return dist <= window


def load_unit_sets(species: str, group: str, tools: list[str],
                   universe: set, ref: dict, logger) -> dict:
    """{replicate_tag: {tool: {"calls": set, "hits": set}}} on the common universe."""
    out: dict = {}
    for r in ev.units(PLATFORM, group).to_dict("records"):
        per: dict = {}
        for t in tools:
            f = CALLSET_ROOT / PLATFORM / species / group / MOD / t / f"{r['sample']}.tsv"
            if not f.exists():
                continue
            d = ev.site_keys(f)
            calls = {(c, int(p)) for c, p in zip(d["chrom"], d["pos"])
                     if (c, int(p)) in universe}
            by_chrom: dict[str, list[int]] = defaultdict(list)
            for c, p in calls:
                by_chrom[c].append(p)
            hits: set = set()
            for c, plist in by_chrom.items():
                arr = ref.get(c)
                if arr is None or arr.size == 0:
                    continue
                p = np.asarray(plist, dtype=np.int64)
                m = _hit_mask(p, arr, WINDOW)
                hits.update((c, int(x)) for x in p[m])
            per[t] = {"calls": calls, "hits": hits}
        out[r["replicate_tag"]] = per
        logger.info("[%s/%s] unit loaded: %s (%d tools)", group, r["replicate_tag"],
                    r["sample"], len(per))
    return out


def _combo_one_unit(per: dict, combo: tuple[str, ...], ref_n: int,
                    ref_n_all: int = 0) -> dict:
    """One unit's union / intersection metrics.

    ``ref_n``   = |GLORI inside the shared measurable universe| -- the recall
                  denominator of the revision (what is measurable at all);
    ``ref_n_all`` = |GLORI| of the whole group (every annotated site, no
                  coverage filter) -- the *second* denominator, added 2026-09-21
                  so the figure can state the absolute GLORI coverage next to
                  the in-universe one (reviewer R3-4: "the absolute coverage
                  remains relatively low").  The numerator is the same in both
                  cases (in-universe calls within 2 bp of any GLORI site), so the
                  two recalls differ by the denominator only.
    """
    present = [t for t in combo if t in per]
    res: dict = {}
    if not present:
        return res
    calls = [per[t]["calls"] for t in present]
    hits = [per[t]["hits"] for t in present]
    u_calls = set().union(*calls)
    u_hits = set().union(*hits)
    res["union_n"] = len(u_calls)
    res["union_tp"] = len(u_hits)
    res["union_precision"] = len(u_hits) / len(u_calls) if u_calls else np.nan
    res["union_recall"] = len(u_hits) / ref_n if ref_n else np.nan
    if ref_n_all:                       # only computed when the caller asks for it
        res["union_recall_all"] = len(u_hits) / ref_n_all
    if len(present) == len(combo):
        i_calls = set.intersection(*calls)
        res["isect_n"] = len(i_calls)
        res["isect_tp"] = len(i_calls & u_hits)
        res["isect_precision"] = (res["isect_tp"] / len(i_calls)) if i_calls else np.nan
        res["isect_recall"] = res["isect_tp"] / ref_n if ref_n else np.nan
        if ref_n_all:
            res["isect_recall_all"] = res["isect_tp"] / ref_n_all
    return res


def score_combination(units_per: dict, combo: tuple[str, ...], ref_n: int,
                      ref_n_all: int = 0) -> dict:
    """Mean/SD of one combination over the units scored together."""
    df = pd.DataFrame([_combo_one_unit(per, combo, ref_n, ref_n_all)
                       for per in units_per.values()])
    out = {"combination": "+".join(combo), "k": len(combo),
           "n_units": len(units_per)}
    for key in ("union_recall", "union_precision", "union_n", "union_tp",
                "union_recall_all", "isect_recall", "isect_precision",
                "isect_n", "isect_tp", "isect_recall_all"):
        col = df[key] if key in df else pd.Series(dtype=float)
        out[f"{key}_mean"] = float(np.nanmean(col)) if len(col) else np.nan
        out[f"{key}_sd"] = (float(np.nanstd(col, ddof=1))
                            if len(col) > 1 else np.nan)
    return out


def per_unit_selected(units_per: dict, combos: list[tuple], ref_n: int,
                      label: str, ref_n_all: int = 0) -> pd.DataFrame:
    rows = []
    for tag, per in units_per.items():
        for combo in combos:
            m = _combo_one_unit(per, combo, ref_n, ref_n_all)
            if not m:
                continue
            rows.append({"group": label, "unit": tag,
                         "combination": "+".join(combo), "k": len(combo),
                         **m})
    return pd.DataFrame(rows)


def greedy_sequence(units_per: dict, ref_n: int, kmax: int) -> list[dict]:
    """Forward selection: add the tool with the largest marginal mean-recall gain."""
    order, seq = [], []
    for _ in range(kmax):
        best, best_r = None, -np.inf
        for t in TOOLS:
            if t in order:
                continue
            m = score_combination(units_per, tuple(order + [t]), ref_n)
            if np.isnan(m["union_recall_mean"]):
                continue
            if m["union_recall_mean"] > best_r:
                best, best_r = t, m["union_recall_mean"]
        if best is None:
            break
        order.append(best)
        m = score_combination(units_per, tuple(order), ref_n)
        seq.append({"step": len(order), "tool_added": best,
                    "combination": "+".join(order),
                    "k": len(order),
                    "union_recall_mean": m["union_recall_mean"],
                    "union_precision_mean": m["union_precision_mean"],
                    "isect_precision_mean": m["isect_precision_mean"],
                    "isect_recall_mean": m["isect_recall_mean"]})
    return seq


# --------------------------------------------------------------------------- #
# negative-control FP burden
# --------------------------------------------------------------------------- #
def _control_region_bp(species: str, group: str) -> float:
    f = TABLE_DIR / "controls_ivt_fpr.tsv"
    if f.exists():
        d = pd.read_csv(f, sep="\t")
        d = d[(d.species == species) & (d.dataset_group == group) & (d.min_cov == 10)]
        if len(d):
            return float(d.region_bp.median())
    return np.nan


def control_fp(combo: tuple[str, ...]) -> pd.DataFrame:
    """FP/10 kb of a combination's union on unmodified controls, per unit."""
    rows = []
    for cname, (cspecies, cgroup) in CONTROLS.items():
        region_bp = _control_region_bp(cspecies, cgroup)
        for r in ev.units(PLATFORM, cgroup).to_dict("records"):
            n, present = 0, 0
            for t in combo:
                f = (CALLSET_ROOT / PLATFORM / cspecies / cgroup / MOD / t
                     / f"{r['sample']}.tsv")
                if not f.exists():
                    continue
                present += 1
                n += len(ev.site_keys(f))
            if not present:
                continue
            rows.append({"control": cname, "sample": r["sample"],
                         "sequencing_unit": r["sequencing_unit"],
                         "combination": "+".join(combo), "k": len(combo),
                         "tools_available": present, "n_calls": n,
                         "region_bp": region_bp,
                         "fp_per_10kb": (n / region_bp * 10_000
                                         if region_bp else np.nan)})
    return pd.DataFrame(rows)


def tool_effects(search: pd.DataFrame) -> pd.DataFrame:
    """Main effect of every configuration at fixed k, from the enumeration itself.

    For each group, k and tool: the mean ``union_recall_mean`` of the enumerated
    combinations of size k that contain the tool, minus the mean over **all**
    combinations of that k.  It answers "which configurations actually raise the
    criterion's objective once the number of tools is fixed" without ever looking
    at the selected optimum alone (the selected sets are marked separately by
    ``fig6_selected_members.tsv``).
    """
    rows: list[dict] = []
    for (gname, k), grp in search.groupby(["group", "k"]):
        ref_r = float(grp.union_recall_mean.mean())
        ref_p = float(grp.union_precision_mean.mean())
        parts = grp.combination.str.split("+")
        for t in TOOLS:
            m = grp[parts.apply(lambda c: t in c)]
            if not len(m):
                continue
            rows.append({"species": grp.species.iloc[0], "group": gname,
                         "k": int(k), "tool": t, "n_combos": int(len(m)),
                         "effect_recall": float(m.union_recall_mean.mean()) - ref_r,
                         "effect_precision": (float(m.union_precision_mean.mean())
                                              - ref_p)})
    return pd.DataFrame(rows)[["species", "group", "k", "tool", "n_combos",
                               "effect_recall", "effect_precision"]]


def control_fp_by_tool() -> pd.DataFrame:
    """FP/10 kb of **every** single tool on the unmodified controls, per unit.

    The same quantity as :func:`control_fp` at k = 1, resolved per tool: it tells
    whether the control burden of a combination is inherited from one tool or
    spread over all of them.  Only read here -- drawn by 65 / 62 from the frozen
    table.
    """
    frames = []
    for t in TOOLS:
        c = control_fp((t,))
        if len(c):
            frames.append(c.assign(tool=t).drop(columns=["combination"]))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def _curve_handles() -> list:
    from matplotlib.lines import Line2D
    return [Line2D([], [], color=UNION_C, marker="o", ms=5, lw=1.7,
                   label="union (recall / coverage)"),
            Line2D([], [], color=ISECT_C, marker="s", ms=5, lw=1.7,
                   label="intersection (precision)"),
            Line2D([], [], color="0.35", lw=1.0, ls=(0, (3, 2)),
                   label="best single tool")]


def _tier_handles() -> list:
    """Markers for the three selected tiers, shared by panels C and D."""
    from matplotlib.lines import Line2D
    return [Line2D([], [], color="0.3", marker="o", ms=6, ls="",
                   label="k = 1 (best single tool)"),
            Line2D([], [], color="#b03a2e", marker="*", ms=9, ls="",
                   label="k = 2 (selected pair)"),
            Line2D([], [], color="0.3", marker="s", ms=6, ls="",
                   label=f"k = {KMAX} (selected combination)"),
            Line2D([], [], color="0.78", lw=0.9, ls="-",
                   label="one line = one unit")]


def _panel_letter(ax, s: str) -> None:
    """Bold row letter in the title band, left of the centred title."""
    ax.annotate(s, xy=(0.0, 1.0), xycoords="axes fraction",
                xytext=(-52, 16), textcoords="offset points",
                fontsize=13, fontweight="bold", ha="center", va="bottom",
                annotation_clip=False)


def fig_main(sel: dict[str, pd.DataFrame], units_sel: dict[str, pd.DataFrame],
             ctl: pd.DataFrame, pairs: dict[str, pd.DataFrame],
             n_text: dict[str, str]) -> None:
    apply_style()
    species_list = list(GROUPS)
    fig, axes = plt.subplots(4, 3, figsize=(5.0 * 3, 3.4 * 4),
                             constrained_layout=True,
                             gridspec_kw={"height_ratios": [1.2, 1.2, 1.0, 1.0]})
    for j, sp in enumerate(species_list):
        colour = COLOUR[sp]
        style = GROUPS[sp][1]
        sel_df = sel[sp].sort_values(["group", "k"])

        # A -- union recall vs k
        ax = axes[0, j]
        for gi, (gname, sub) in enumerate(sel_df.groupby("group")):
            ls = "-" if gi == 0 else (0, (4, 2))
            if style == "studies":
                ax.plot(sub.k, 100 * sub.union_recall_mean, marker="o", ms=5,
                        lw=1.7, ls=ls, color=colour)
            else:
                ax.errorbar(sub.k, 100 * sub.union_recall_mean,
                            yerr=100 * sub.union_recall_sd, marker="o", ms=5,
                            lw=1.7, ls=ls, capsize=3, color=colour)
            # per-group best-single-tool reference (Mouse: one per study)
            ax.axhline(100 * sub[sub.k == 1].union_recall_mean.iloc[0],
                       color="0.35", lw=1.0, ls=ls)
        ax.set_xticks(range(1, KMAX + 1))
        ax.set_xlabel("tools combined (union)")
        ax.set_ylabel("GLORI coverage / recall (%)")
        ax.set_title(f"{sp}\n{n_text[sp]}", fontweight="bold", fontsize=11)

        # B -- precision cost
        ax = axes[1, j]
        for gi, (gname, sub) in enumerate(sel_df.groupby("group")):
            ls = "-" if gi == 0 else (0, (4, 2))
            ax.errorbar(sub.k, 100 * sub.union_precision_mean,
                        yerr=100 * sub.union_precision_sd, marker="o", ms=5,
                        lw=1.7, ls=ls, capsize=3, color=UNION_C)
            ax.errorbar(sub.k, 100 * sub.isect_precision_mean,
                        yerr=100 * sub.isect_precision_sd, marker="s", ms=5,
                        lw=1.7, ls=ls, capsize=3, color=ISECT_C)
        ax.set_xticks(range(1, KMAX + 1))
        ax.set_xlabel("tools combined")
        ax.set_ylabel("PPV vs. GLORI (2 bp)")
        ax.set_title(f"{sp}\n{n_text[sp]}", fontweight="bold", fontsize=11)

        # C -- two-tool landscape (restored from the published panel C)
        ax = axes[2, j]
        p = pairs[sp]
        for gi, (gname, sub) in enumerate(p.groupby("group")):
            alpha = 0.85 if gi == 0 else 0.4
            ax.scatter(100 * sub.union_recall_mean, 100 * sub.isect_precision_mean,
                       s=22, color=colour, alpha=alpha, edgecolor="white",
                       linewidth=0.3, zorder=2)
        chosen = sel_df[sel_df.k.isin((1, 2, KMAX))]
        tier_mk = {1: ("o", "0.25", 70), 2: ("*", "#b03a2e", 110),
                   KMAX: ("s", "0.25", 70)}
        for k, (mk, mcol, ms) in tier_mk.items():
            sub = chosen[chosen.k == k]
            if len(sub):
                ax.scatter(100 * sub.union_recall_mean,
                           100 * sub.isect_precision_mean, marker=mk, s=ms,
                           color=mcol, zorder=3, edgecolor="white", linewidth=0.6)
        per_g = int(p.groupby("group").size().iloc[0])
        note = f"{per_g} two-tool pairs" + (" per study" if style == "studies" else "")
        ax.set_xlabel("union recall (%)")
        ax.set_ylabel("intersection PPV vs. GLORI (2 bp)")
        ax.set_title(f"{sp}\n{note}", fontweight="bold", fontsize=11)

        # D -- negative-control FP burden, k = 1 / 2 / 5 trajectories
        ax = axes[3, j]
        sub = ctl[ctl.species == sp]
        names = list(CONTROLS)
        off = {1: -0.2, 2: 0.0, KMAX: 0.2}
        for i, name in enumerate(names):
            s2 = sub[sub.control == name]
            for unit, s3 in s2.groupby("sequencing_unit"):
                xs, ys, ks = [], [], []
                for k in (1, 2, KMAX):
                    r = s3[s3.k == k]
                    if len(r):
                        xs.append(i + off[k])
                        ys.append(float(r.fp_per_10kb.iloc[0]))
                        ks.append(k)
                if len(xs) > 1:
                    ax.plot(xs, ys, color="0.78", lw=0.9, zorder=1)
                for x, y, k in zip(xs, ys, ks):
                    if k == 1:
                        ax.scatter([x], [y], marker="o", s=30, facecolor="white",
                                   edgecolor=colour, linewidth=1.4, zorder=2)
                    elif k == 2:
                        ax.scatter([x], [y], marker="*", s=80, color="#b03a2e",
                                   zorder=2, edgecolor="white", linewidth=0.4)
                    else:
                        ax.scatter([x], [y], marker="s", s=30, color=colour,
                                   zorder=2)
        ax.set_yscale("log")
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, fontsize=9)
        ax.set_ylabel("FP / 10 kb (control)")
        ax.set_xlabel("unmodified control")
        ax.set_title(f"{sp}\n{_control_note(sub)}", fontweight="bold", fontsize=11)
    for ax_row, letter in zip(axes[:, 0], ("A", "B", "C", "D")):
        _panel_letter(ax_row, letter)
    fig.legend(handles=_curve_handles() + _tier_handles(), loc="upper center",
               ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.02))
    save(fig, str(FIG / "Fig6_combination"))


def _control_note(sub: pd.DataFrame) -> str:
    if not len(sub):
        return "no control callsets"
    parts = [f"{name}: {int(s.sequencing_unit.nunique())} units"
             for name, s in sub.groupby("control")]
    return " | ".join(parts)


def fig_supp(greedy13: dict[str, pd.DataFrame],
             units_sel: dict[str, pd.DataFrame],
             n_text: dict[str, str]) -> None:
    apply_style()
    species_list = list(GROUPS)
    fig, axes = plt.subplots(3, 3, figsize=(5.0 * 3, 3.6 * 3),
                             constrained_layout=True)
    for j, sp in enumerate(species_list):
        colour = COLOUR[sp]

        # row 1 -- union recall over the full 1..13 path
        ax = axes[0, j]
        for gi, (gname, sub) in enumerate(greedy13[sp].groupby("group")):
            ls = "-" if gi == 0 else (0, (4, 2))
            ax.plot(sub.k, 100 * sub.union_recall_mean, marker="o", ms=4,
                    lw=1.6, ls=ls, color=colour)
            # per-group best-single-tool reference (Mouse: one per study)
            ax.axhline(100 * sub[sub.k == 1].union_recall_mean.iloc[0],
                       color="0.35", lw=1.0, ls=ls)
        ax.set_xticks(range(1, KMAX_SUPP + 1, 2))
        ax.set_xlabel("tools combined (greedy order)")
        ax.set_ylabel("GLORI coverage / recall (%)")
        ax.set_title(f"{sp}\n{n_text[sp]}", fontweight="bold", fontsize=11)

        # row 2 -- union vs intersection precision over the same path
        ax = axes[1, j]
        for gi, (gname, sub) in enumerate(greedy13[sp].groupby("group")):
            ls = "-" if gi == 0 else (0, (4, 2))
            ax.plot(sub.k, 100 * sub.union_precision_mean, marker="o", ms=4,
                    lw=1.6, ls=ls, color=UNION_C)
            ax.plot(sub.k, 100 * sub.isect_precision_mean, marker="s", ms=4,
                    lw=1.6, ls=ls, color=ISECT_C)
        ax.set_xticks(range(1, KMAX_SUPP + 1, 2))
        ax.set_xlabel("tools combined (greedy order)")
        ax.set_ylabel("PPV vs. GLORI (2 bp)")
        ax.set_title(f"{sp}\n{n_text[sp]}", fontweight="bold", fontsize=11)

        # row 3 -- per-unit stability (moved here from main panel D)
        ax = axes[2, j]
        us = units_sel[sp]
        tags = list(us["unit"].drop_duplicates())
        for i, tag in enumerate(tags):
            one = us[(us.unit == tag) & (us.k == 1)]
            five = us[(us.unit == tag) & (us.k == KMAX)]
            if len(one) and len(five):
                y1 = 100 * one.union_recall.iloc[0]
                y5 = 100 * five.union_recall.iloc[0]
                ax.plot([i, i], [y1, y5], color="0.75", lw=1.0, zorder=1)
                ax.scatter([i], [y1], marker="o", s=30, facecolor="white",
                           edgecolor=colour, linewidth=1.4, zorder=2)
                ax.scatter([i], [y5], marker="s", s=30, color=colour, zorder=2)
        ax.set_xticks(range(len(tags)))
        ax.set_xticklabels(tags, rotation=20, ha="right", fontsize=9)
        ax.set_ylabel("GLORI coverage / recall (%)")
        ax.set_xlabel("independent unit")
        ax.set_title(f"{sp}\nk = 1 vs selected k = {KMAX}", fontweight="bold",
                     fontsize=11)
    for ax_row, letter in zip(axes[:, 0], ("A", "B", "C")):
        _panel_letter(ax_row, letter)
    tier = {h.get_label(): h for h in _tier_handles()}
    fig.legend(handles=_curve_handles()
               + [tier["k = 1 (best single tool)"],
                  tier[f"k = {KMAX} (selected combination)"],
                  tier["one line = one unit"]],
               loc="upper center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    # superseded 3-row S5 layout (off by default): written under an explicitly
    # legacy name so that the canonical ``FigureS5_rev`` has one producer only
    save(fig, str(FIG / "FigS5_legacy_3row"))


# --------------------------------------------------------------------------- #
def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--species", action="append", default=None,
                    help="restrict to a subset (smoke tests); default = all three")
    ap.add_argument("--no-figures", action="store_true",
                    help="tables only (debug / partial runs)")
    ap.add_argument("--fig6-legacy", action="store_true",
                    help="also draw the superseded oversized-canvas main figure; "
                         "the revision Figure 6 is produced by "
                         "65_fig6_combination_page.py")
    ap.add_argument("--s5-legacy", action="store_true",
                    help="also draw the superseded three-row S5 layout; the "
                         "revision figure is produced by 62_figS5_figure.py")
    args = ap.parse_args()
    if args.species:
        unknown = [s for s in args.species if s not in GROUPS]
        if unknown:
            sys.exit(f"unknown species: {unknown}")
        for s in [s for s in GROUPS if s not in args.species]:
            GROUPS.pop(s)
    _ = LOG  # logs are written by setup_logger
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    logger = setup_logger("40_fig6_combination")
    inv = Inventory("40_fig6_combination")

    sel_by_species, units_sel, ctl_rows = {}, {}, []
    greedy13, pairs, unit_meta = {}, {}, {}
    space_frames: list[pd.DataFrame] = []
    units_allk: list[pd.DataFrame] = []
    member_frames: list[pd.DataFrame] = []
    for sp, (group, style) in GROUPS.items():
        universe = ev.common_universe(PLATFORM, sp, group)
        ref = ev.reference(sp)
        ref_u = ev.reference_positions_in_universe(ref, universe)
        ref_n = len(ref_u)
        # every GLORI site of the group, without the measurable-universe filter:
        # the second recall denominator (reviewer R3-4, "absolute coverage")
        ref_n_all = int(sum(len(a) for a in ref.values()))
        p0 = ref_n / len(universe) if universe else np.nan
        per_unit = load_unit_sets(sp, group, TOOLS, universe, ref, logger)
        unit_meta[sp] = (per_unit, ref_n, style, group, p0)
        logger.info("[%s] common universe=%d, GLORI-in-universe=%d, GLORI-total=%d,"
                    " chance p0=%.4f", group, len(universe), ref_n, ref_n_all, p0)

        # units scored together: the whole group, or one mouse study at a time
        if style == "studies":
            unit_groups = {tag: {tag: v} for tag, v in per_unit.items()}
        else:
            unit_groups = {group: per_unit}

        sel_frames, greedy_frames = [], []
        for gname, units_per in unit_groups.items():
            combos = [c for k in range(1, KMAX + 1)
                      for c in itertools.combinations(TOOLS, k)]
            df = pd.DataFrame([score_combination(units_per, c, ref_n, ref_n_all)
                               for c in combos])
            df = df[df.union_n_mean.notna()].copy()
            df["species"], df["group"] = sp, gname
            df["p0_chance_precision"] = p0
            ok = df[df.union_precision_mean >= p0]
            best = (ok.sort_values(["k", "union_recall_mean", "union_precision_mean"],
                                   ascending=[True, False, False])
                    .groupby("k").head(1))
            sel_frames.append(best)
            # every enumerated combination of this group, with the feasibility
            # flag that defines the criterion (mean union PPV >= chance level
            # p0) and the per-k selected optimum -> figS5_search_space.tsv
            space = df.copy()
            space["feasible"] = space.union_precision_mean >= p0
            space["selected"] = space.combination.isin(set(best.combination))
            space_frames.append(space)
            g = pd.DataFrame(greedy_sequence(units_per, ref_n, KMAX))
            g["species"], g["group"] = sp, gname
            greedy_frames.append(g)
            # greedy vs exhaustive agreement per k
            agree = []
            for r in best.itertuples():
                gr = g[g.k == r.k]
                agree.append(bool(len(gr) and set(gr.iloc[0]["combination"].split("+"))
                                  == set(r.combination.split("+"))))
            best = best.assign(greedy_agrees=agree)
            sel_frames[-1] = best
            logger.info("[%s/%s] k<=%d exhaustive done; selected: %s", sp, gname,
                        KMAX, "; ".join(f"k={r.k}:{r.combination}"
                                        for r in best.itertuples()))
            # ---- frozen evidence added 2026-09-21 (Figure 6 / S5 rebuild) ---- #
            # membership matrix: which of the 13 tools each selected k uses
            member_frames.append(pd.DataFrame(
                [{"species": sp, "group": gname, "k": int(r.k), "tool": t,
                  "in_combination": bool(t in r.combination.split("+"))}
                 for r in best.itertuples() for t in TOOLS])[MEMBER_COLUMNS])
            # per-unit values of EVERY selected k = 1..5 (not only k = 1 / 5)
            combos_all = [tuple(best[best.k == k].iloc[0]["combination"].split("+"))
                          for k in range(1, KMAX + 1)]
            ak = per_unit_selected(units_per, combos_all, ref_n, gname, ref_n_all)
            ak.insert(0, "species", sp)
            units_allk.append(ak[ALLK_COLUMNS])
            # full greedy path to 13 tools (S5 row 1/2)
            g13 = pd.DataFrame(greedy_sequence(units_per, ref_n, KMAX_SUPP))
            g13["species"], g13["group"] = sp, gname
            # pad to a metric-compatible frame
            greedy13.setdefault(sp, []).append(g13)
            # two-tool pairs (S5 row 3)
            pr = pd.DataFrame([score_combination(units_per, c, ref_n, ref_n_all)
                               for c in itertools.combinations(TOOLS, 2)])
            pr["species"], pr["group"] = sp, gname
            pairs.setdefault(sp, []).append(pr)
            # per-unit table for the chosen k=1 / k=KMAX sets (panel C2)
            base = best[best.k == 1].iloc[0]["combination"].split("+")
            top = best[best.k == KMAX].iloc[0]["combination"].split("+")
            units_sel.setdefault(sp, []).append(
                per_unit_selected(units_per, [tuple(base), tuple(top)], ref_n, gname,
                                  ref_n_all))

        sel_by_species[sp] = pd.concat(sel_frames, ignore_index=True)
        # negative-control FP burden for the selected k=1 / k=2 / k=KMAX sets
        for gname, sub in sel_by_species[sp].groupby("group"):
            for k in range(1, KMAX + 1):
                combos_k = sub[sub.k == k]
                if not len(combos_k):
                    continue
                combo = tuple(combos_k.iloc[0]["combination"].split("+"))
                c = control_fp(combo)
                if len(c):
                    ctl_rows.append(c.assign(species=sp, chosen_for=gname))
        greedy13[sp] = pd.concat(greedy13[sp], ignore_index=True)
        pairs[sp] = pd.concat(pairs[sp], ignore_index=True)
        units_sel[sp] = pd.concat(units_sel[sp], ignore_index=True)

    ctl_df = pd.concat(ctl_rows, ignore_index=True) if ctl_rows else pd.DataFrame()
    ctl_tool = control_fp_by_tool()
    n_text = {}
    for sp, (per_unit, ref_n, style, group, p0) in unit_meta.items():
        n_text[sp] = (f"{len(per_unit)} biological replicates" if style == "replicates"
                      else f"{len(per_unit)} independent studies (never averaged)")

    write_table(pd.concat(sel_by_species.values(), ignore_index=True),
                TAB / "fig6_combination_selected.tsv")
    if len(ctl_df):
        write_table(ctl_df, TAB / "fig6_negative_control_fp.tsv")
    write_table(pd.concat(units_sel.values(), ignore_index=True),
                TAB / "fig6_per_unit_selected.tsv")
    write_table(pd.concat(units_allk, ignore_index=True),
                TAB / "fig6_per_unit_allk.tsv")
    write_table(pd.concat(member_frames, ignore_index=True),
                TAB / "fig6_selected_members.tsv")
    write_table(pd.concat(greedy13.values(), ignore_index=True),
                TAB / "figS5_greedy_1to13.tsv")
    write_table(pd.concat(pairs.values(), ignore_index=True),
                TAB / "figS5_two_tool_pairs.tsv")
    search = pd.concat(space_frames, ignore_index=True)[SEARCH_COLUMNS]
    n_sel = int(search.selected.sum())
    write_table(search, TAB / "figS5_search_space.tsv")
    # single-tool plane (k = 1, all 13 configurations) and the per-tool control
    # burden -- both added 2026-09-21 for the rebuilt Figure 6 / S5
    single = (search[search.k == 1].rename(columns={"combination": "tool"})
              .loc[:, SINGLE_COLUMNS].reset_index(drop=True))
    write_table(single, TAB / "fig6_single_tool_metrics.tsv")
    if len(ctl_tool):
        write_table(ctl_tool[CTL_TOOL_COLUMNS],
                    TAB / "fig6_negative_control_fp_bytool.tsv")
    effects = tool_effects(search)
    write_table(effects, TAB / "fig6_tool_effects.tsv")
    logger.info("tool effects: %d rows (%d tools x %d k x %d groups)",
                len(effects), effects.tool.nunique(), effects.k.nunique(),
                effects.group.nunique())
    logger.info("search space: %d combinations (%d selected) -> figS5_search_space.tsv",
                len(search), n_sel)
    logger.info("single-tool plane: %d configurations; per-tool control FP: %d rows",
                len(single), len(ctl_tool))

    if not args.no_figures:
        if args.fig6_legacy:
            fig_main(sel_by_species, units_sel, ctl_df, pairs, n_text)
        if args.s5_legacy:
            fig_supp(greedy13, units_sel, n_text)
        if not (args.fig6_legacy or args.s5_legacy):
            logger.info("legacy figures off: Figure 6 -> "
                        "65_fig6_combination_page.py, S5 -> 62_figS5_figure.py")
    inv.flush()
    logger.info("figures -> %s", FIG)


if __name__ == "__main__":
    main()
