#!/usr/bin/env python
"""Figure 5A, revision rebuild (R1-2 / R1-8 / R2-1 / R3-2 / R3-7).

The published panel A was a rank heatmap of Precision / Recall / F1 / MCC per
tool and species, computed by ``code/next_postprocessing/Total_new1/NGS.ipynb``
from the legacy ``output/<group>/Tools.txt`` aggregates with

    TP = interval-overlap hits (CLOSED GLORI interval, not a nucleotide)
    FN = n_GLORI - TP          (every GLORI site treated as measurable)
    TN = genome_total - n_GLORI - FP   (a hard-coded constant per species)

so precision/recall ignored coverage, MCC/specificity were built on an invented
TN, and the values carried no replicate information (one legacy aggregate per
species; mouse = one study).

This script recomputes the very same four metrics from the harmonisation evaluation
layer, where every quantity has an explicit definition
(``docs/pipeline.md``): an explicit candidate universe
U(sample, c = 10) and TP/FP/FN/TN inside it, at the primary window w = 2 bp,
one row per **independent sequencing unit**.

Ranking rule (author decision 2026-09-19): ranks are computed *within each
independent unit* (``method="min"``, high value = rank 1, NaN never ranks) and
the panel shows the mean rank across the units of that species; the rank range,
n and k_of_n live in the tables.

Replicate structure (never merged):
  * Arabidopsis_WT  -> 3 biological replicates (rep1/rep2/rep3);
  * HeLa_WT         -> 3 biological replicates (HeLa_WT1/2/3);
  * Mouse_WT        -> TWO independent studies (SRP166020 / SRP357195).  The
    panel draws the single unit ``mES_WT`` -- the same unit the published mouse
    panel was built from (its 13 m6A tools are all filled, while the other
    study's Nanocompore callset is a real zero).  The second study is written to
    ``fig5a_mouse_cross_study_check.tsv`` (rank contrast + Spearman rho) and is
    NEVER averaged into the drawn column.

``purified`` (WT-call-minus-KO/KD) is deliberately absent: the published panel A
was WT-only, and a purified set defined by KO/KD absence is circular as an
accuracy argument (R3-7).

Outputs -> figures/figure5/{tables,figures}

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/figures/figure5/src/35_fig5a_metric_ranks.py
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
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
import seaborn as sns

HERE = _RB / "src/harmonisation"
sys.path.insert(0, str(HERE))

from common import config as C                                        # noqa: E402
from common.figstyle import apply as apply_style                     # noqa: E402
from common.figstyle import save                                     # noqa: E402
from common.io_utils import write_table                              # noqa: E402
from common.manifest import setup_logger                             # noqa: E402

CONFUSION = C.TABLE_DIR / "m6a_glori_confusion.tsv"
OUT = (_RB / "figures/figure5")
TAB, FIG = OUT / "tables", OUT / "figures"

METRICS = ["precision", "recall", "f1", "mcc"]
METRIC_LABEL = {"precision": "Precision", "recall": "Recall",
                "f1": "F1", "mcc": "MCC"}
#: row ordering key of every species panel (author decision 2026-09-19: PPV, like
#: the published panel which sorted by Precision; F1 is then not monotonic).
SORT_METRIC = "precision"

#: publication tool list and order (same set as 24/28)
TOOL_ORDER = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
              "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
              "NanoSPA_m6A", "xPore", "yanocomp"]
DISPLAY = {"yanocomp": "Yanocomp"}

#: (species, dataset_group, [independent sequencing units drawn])
PANELS: list[tuple[str, str, list[str]]] = [
    ("Arabidopsis", "Arabidopsis_WT",
     ["Arabidopsis_WT_rep1", "Arabidopsis_WT_rep2", "Arabidopsis_WT_rep3"]),
    # cross-study group: exactly ONE unit is drawn (see module docstring)
    ("Mouse", "Mouse_WT", ["mES_WT"]),
    ("Human", "HeLa_WT", ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"]),
]
#: second mouse study, tabulated only
MOUSE_STUDY_A = "mESCs_Mettl3_WT"


# --------------------------------------------------------------------------- #
def load_metrics(window: int = C.PRIMARY_WINDOW,
                 min_cov: int = C.C_MIN_DEFAULT) -> pd.DataFrame:
    """One row per (group, unit sample, tool) at the primary window/coverage."""
    df = pd.read_csv(CONFUSION, sep="\t")
    df = df[(df["window"] == window) & (df["min_cov"] == min_cov)]
    df = df[df["platform"] == C.PLATFORM_RNA002]
    keep = ["species", "dataset_group", "sample", "replicate_tag", "tool",
            *METRICS, "precision_ci_lo", "precision_ci_hi",
            "recall_ci_lo", "recall_ci_hi",
            "n_calls_in_universe", "n_reference_in_universe", "universe",
            "window", "min_cov"]
    df = df[[c for c in keep if c in df.columns]].copy()
    for m in METRICS:
        df[m] = pd.to_numeric(df[m], errors="coerce")
    return df


def rank_within_unit(per_unit: pd.DataFrame) -> pd.DataFrame:
    """Add ``<metric>_rank`` **within each independent unit** (1 = highest metric).

    The grouping is ``(dataset_group, sample)``: one sample is one sequencing
    unit, so its 13 tools are ranked against each other (NaN metric -> no rank).
    """
    ranked = per_unit.copy()
    for m in METRICS:
        ranked[f"{m}_rank"] = ranked.groupby(["dataset_group", "sample"])[m].rank(
            method="min", ascending=False)                   # NaN stays NaN
    return ranked


def summarise_ranks(ranked: pd.DataFrame) -> pd.DataFrame:
    """Long table: species x group x tool x metric, rank and value statistics."""
    rows = []
    for species, group, units in PANELS:
        n_total = len(units)
        # ONLY the drawn units: the Mouse group also holds the second study, which
        # must never enter the summary (it is tabulated separately).
        sub = ranked[(ranked["dataset_group"] == group)
                     & (ranked["sample"].isin(units))]
        for tool in TOOL_ORDER:
            t = sub[sub["tool"] == tool]
            for m in METRICS:
                rk = t[f"{m}_rank"].dropna().astype(float)
                val = t[m].dropna().astype(float)
                rows.append({
                    "species": species, "dataset_group": group, "tool": tool,
                    "metric": m,
                    "n_units": int(rk.size), "n_units_total": n_total,
                    "k_of_n": f"{int(rk.size)}/{n_total}",
                    "rank_mean": float(rk.mean()) if rk.size else np.nan,
                    "rank_min": float(rk.min()) if rk.size else np.nan,
                    "rank_max": float(rk.max()) if rk.size else np.nan,
                    "metric_mean": float(val.mean()) if val.size else np.nan,
                    "metric_sd": float(val.std(ddof=1)) if val.size > 1 else np.nan,
                    "metric_min": float(val.min()) if val.size else np.nan,
                    "metric_max": float(val.max()) if val.size else np.nan,
                    "metric_ci_lo": (float(pd.to_numeric(
                        t["precision_ci_lo"], errors="coerce").mean())
                        if m == "precision" else np.nan),
                    "metric_ci_hi": (float(pd.to_numeric(
                        t["precision_ci_hi"], errors="coerce").mean())
                        if m == "precision" else np.nan),
                })
    return pd.DataFrame(rows)


def mouse_cross_study_check(ranked: pd.DataFrame) -> pd.DataFrame:
    """mES_WT (drawn) vs mESCs_Mettl3_WT (tabulated) rank contrast, never pooled."""
    rows = []
    a = ranked[ranked["sample"] == MOUSE_STUDY_A].set_index("tool")
    b = ranked[ranked["sample"] == "mES_WT"].set_index("tool")
    for m in METRICS:
        col = f"{m}_rank"
        ra, rb = a[col].astype(float), b[col].astype(float)
        both = pd.DataFrame({"studyA": ra, "studyB": rb}).dropna()
        rho = (both["studyA"].corr(both["studyB"], method="spearman")
               if len(both) > 2 else np.nan)
        rows.append({
            "metric": m, "n_tools_both": int(len(both)),
            "rank_spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
            "mean_abs_rank_diff": (float((both["studyA"] - both["studyB"]).abs().mean())
                                   if len(both) else np.nan),
            "max_abs_rank_diff": (float((both["studyA"] - both["studyB"]).abs().max())
                                  if len(both) else np.nan),
            "studyA": "mESCs_Mettl3_WT (SRP166020, tabulated only)",
            "studyB": "mES_WT (SRP357195, drawn in the panel)",
        })
    detail = pd.DataFrame({"tool": pd.Index(a.index).union(b.index)}).set_index("tool")
    for m in METRICS:
        detail[f"{m}_rank_studyA"] = a[f"{m}_rank"].astype(float)
        detail[f"{m}_rank_studyB"] = b[f"{m}_rank"].astype(float)
    detail = detail.reset_index()
    return pd.DataFrame(rows), detail


# --------------------------------------------------------------------------- #
def _rank_matrix(mean_rank: pd.DataFrame, species: str) -> pd.DataFrame:
    """tools x metrics matrix of mean ranks, rows ordered by mean **PPV** rank.

    Sort key (author decision 2026-09-19): each species panel is ordered by the
    mean rank of ``precision`` (= PPV against the GLORI reference) ascending —
    the same key the published panel used (it sorted by Precision descending).
    Consequence: the F1 column is *not* monotonic down the rows, which is
    expected, not an error.
    """
    sp = mean_rank[mean_rank["species"] == species]
    mat = sp.pivot_table(index="tool", columns="metric", values="rank_mean")
    mat = mat.reindex(columns=METRICS)
    order = mat[SORT_METRIC].sort_values(kind="stable").index
    mat = mat.reindex(order)
    missing = [t for t in TOOL_ORDER if t not in mat.index]
    if missing:
        fill = pd.DataFrame(index=missing, columns=METRICS, dtype=float)
        mat = pd.concat([mat, fill])
    return mat


def plot_panel_a(mean_rank: pd.DataFrame, ranked: pd.DataFrame) -> None:
    apply_style()
    n_units = {sp: len(units) for sp, _, units in PANELS}

    fig, axes = plt.subplots(1, 3, figsize=(13.8, 6.2), constrained_layout=True)
    cmap = sns.diverging_palette(10, 240, as_cmap=True)   # rank 1 = red
    cmap.set_bad("0.86")

    for ax, (species, group, units) in zip(axes, PANELS):
        mat = _rank_matrix(mean_rank, species)
        n = len(units)
        # annotate: mean rank, one decimal only when it is not (near) integer
        labels = mat.astype(float).apply(
            lambda col: col.map(lambda v: "" if pd.isna(v) else
                                (f"{v:.0f}" if abs(v - round(v)) < 1e-9 else f"{v:.1f}")))
        sns.heatmap(mat.astype(float), ax=ax, cmap=cmap, vmin=1, vmax=len(mat),
                    center=(1 + len(mat)) / 2,
                    annot=labels, fmt="", annot_kws={"size": 11, "color": "black"},
                    linewidths=0.0, linecolor="none", cbar=False,
                    xticklabels=[METRIC_LABEL[m] for m in METRICS],
                    yticklabels=[DISPLAY.get(t, t) for t in mat.index])
        ax.set_title(f"{species} ({'n = 1 study' if n == 1 else f'n = {n}'})",
                     fontweight="bold", fontsize=15, pad=10)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", labelsize=12, length=0)
        ax.tick_params(axis="y", labelsize=11, length=0)
        for s in ax.spines.values():
            s.set_visible(True)
            s.set_linewidth(1.1)

    norm = Normalize(vmin=1, vmax=len(_rank_matrix(mean_rank, "Arabidopsis")))
    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=axes,
                        fraction=0.022, pad=0.015)
    cbar.set_label("Mean rank  (1 = best)", fontsize=12)
    cbar.ax.tick_params(labelsize=11)
    cbar.outline.set_linewidth(1.0)

    FIG.mkdir(parents=True, exist_ok=True)
    stem = FIG / "Fig5A_metric_ranks"
    save(fig, str(stem))
    print("wrote", stem.with_suffix(".pdf").name, "+ .png")


# --------------------------------------------------------------------------- #
def main() -> None:
    logger = setup_logger("35_fig5a_metric_ranks")
    TAB.mkdir(parents=True, exist_ok=True)

    per_unit = load_metrics()
    panels = {g for _, g, _ in PANELS}
    per_unit = per_unit[per_unit["dataset_group"].isin(panels | {"Mouse_WT"})]
    per_unit = per_unit[per_unit["tool"].isin(TOOL_ORDER)]
    missing_tools = sorted(set(TOOL_ORDER) - set(per_unit["tool"]))
    if missing_tools:
        logger.warning("tools absent from the confusion table: %s",
                       ", ".join(missing_tools))

    ranked = rank_within_unit(per_unit)
    write_table(per_unit.sort_values(["species", "dataset_group", "sample", "tool"]),
                TAB / "fig5a_metrics_per_replicate.tsv")
    write_table(ranked[["species", "dataset_group", "sample", "replicate_tag", "tool",
                        *[f"{m}_rank" for m in METRICS]]]
                .sort_values(["species", "dataset_group", "sample", "tool"]),
                TAB / "fig5a_rank_per_replicate.tsv")

    mean_rank = summarise_ranks(ranked)
    write_table(mean_rank.sort_values(["species", "tool", "metric"]),
                TAB / "fig5a_mean_rank_by_group_tool.tsv")

    check, detail = mouse_cross_study_check(ranked)
    write_table(check, TAB / "fig5a_mouse_cross_study_check.tsv")
    write_table(detail, TAB / "fig5a_mouse_cross_study_rank_detail.tsv")

    # ---- integrity checks ------------------------------------------------ #
    # (a) ranks must live in 1..n of the *ranked* tools, and rank must run in the
    #     right direction: a smaller rank number means a larger metric value.
    #     Ties are allowed (method="min"), so the column need not be a strict
    #     1..n permutation -- only its direction and range are asserted.
    for species, group, units in PANELS:
        for unit in units:
            one = ranked[ranked["sample"] == unit]
            for metric in METRICS:
                rk = pd.to_numeric(one[f"{metric}_rank"], errors="coerce")
                val = pd.to_numeric(one[metric], errors="coerce")
                n_ranked = int(rk.notna().sum())
                if n_ranked == 0:
                    continue
                if rk.dropna().min() < 1 or rk.dropna().max() > n_ranked:
                    raise AssertionError(
                        f"{unit}/{metric}: rank outside 1..{n_ranked}")
                s = pd.DataFrame({"rank": rk, "value": val}).dropna()
                s = s.sort_values("rank")
                if not s["value"].is_monotonic_decreasing:
                    raise AssertionError(
                        f"{unit}/{metric}: rank runs opposite to the metric")
                if int(rk.notna().sum()) != int(val.notna().sum()):
                    raise AssertionError(
                        f"{unit}/{metric}: NaN metric still carries a rank")
    written = per_unit[["dataset_group", "sample", "tool", *METRICS]].set_index(
        ["dataset_group", "sample", "tool"])
    src = pd.read_csv(CONFUSION, sep="\t")
    src = src[(src["window"] == C.PRIMARY_WINDOW) & (src["min_cov"] == C.C_MIN_DEFAULT)]
    src = src.set_index(["dataset_group", "sample", "tool"])
    common = written.index.intersection(src.index)
    drift = 0
    for idx in common:
        for m in METRICS:
            a, b = written.loc[idx, m], src.loc[idx, m]
            if not (pd.isna(a) and pd.isna(b)) and abs(float(a) - float(b)) > 1e-12:
                drift += 1
    if drift:
        raise AssertionError(f"{drift} metric cells differ from m6a_glori_confusion.tsv")
    logger.info("integrity: metrics identical to the confusion table on %d rows; "
                "rank columns are 1..n permutations", len(common))

    for _, r in mean_rank[mean_rank["metric"] == "f1"].iterrows():
        logger.info("%-12s %-14s n=%d mean F1 rank=%.2f range=%s",
                    r["species"], r["tool"], r["n_units"],
                    r["rank_mean"], f"{r['rank_min']:.0f}-{r['rank_max']:.0f}")
    logger.info("mouse cross-study rank agreement: %s",
                "; ".join(f"{r.metric} rho={r.rank_spearman_rho:.2f}"
                          for r in check.itertuples()))

    plot_panel_a(mean_rank, ranked)
    logger.info("output: %s", OUT)


if __name__ == "__main__":
    main()
