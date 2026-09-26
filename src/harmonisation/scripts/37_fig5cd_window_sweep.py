#!/usr/bin/env python
"""Figure 5C / 5D, revision rebuild (R1-3 / R3-3 / E2, see response letter NA-1).

The published panel C ("flanking-sequence extension") plotted, per tool, the
fraction of predictions that fall within a widening window of a GLORI site.  The
reviewers pointed out that widening the window mechanically inflates that
overlap (``output/<group>/Tools.txt``, closed-interval matching, one legacy
sample per species).

The revision therefore reports **both** quantities as a function of the window
size, per independent sequencing unit, from ``m6a_localization_curve.tsv`` (the
same table ``27_fig_window_combination.py`` / FigR13 uses, so no second set of
numbers exists for the same data):

  C  ``hit_rate``               -- GLORI hit rate (precision-like), i.e. the
                                  fraction of the tool's calls within +/- w of a
                                  reference site;
  D  ``localization_accuracy``  -- among the calls that match at +/- w, the
                                  fraction that sits on the exact reference
                                  nucleotide (``tp(w = 0) / tp(w)``).  Its
                                  decline with w is the trade-off the reviewer
                                  asked to be quantified.

Windows 0, 1, 2, 5, 10, 20, 50 bp; the primary window (2 bp) is marked with a
guide line.  Drawn units are the same as Figure 5A/5B: Arabidopsis 3 biological
replicates, HeLa 3, mouse the single study ``mES_WT``; curves are the mean across
the units of a species (the per-unit values stay in the tables).

Outputs -> figures/figure5/{tables,figures}

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/src/harmonisation/scripts/37_fig5cd_window_sweep.py
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

from common import config as C                                        # noqa: E402
from common.figstyle import apply as apply_style                     # noqa: E402
from common.figstyle import save                                     # noqa: E402
from common.io_utils import write_table                              # noqa: E402
from common.manifest import setup_logger                             # noqa: E402

CURVE = C.TABLE_DIR / "m6a_localization_curve.tsv"
OUT = (_RB / "figures/figure5")
TAB, FIG = OUT / "tables", OUT / "figures"

TOOL_ORDER = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
              "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
              "NanoSPA_m6A", "xPore", "yanocomp"]
DISPLAY = {"yanocomp": "Yanocomp"}
#: 13 distinguishable colours (matplotlib tab20 pairs), fixed across C and D
_pal = plt.get_cmap("tab20").colors
TOOL_COLOR = {t: _pal[i % 20] for i, t in enumerate(TOOL_ORDER)}

PANELS: list[tuple[str, str, list[str]]] = [
    ("Arabidopsis", "Arabidopsis_WT",
     ["Arabidopsis_WT_rep1", "Arabidopsis_WT_rep2", "Arabidopsis_WT_rep3"]),
    ("Mouse", "Mouse_WT", ["mES_WT"]),
    ("Human", "HeLa_WT", ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"]),
]
WINDOWS = C.WINDOWS


def load_curve() -> pd.DataFrame:
    d = pd.read_csv(CURVE, sep="\t")
    d = d[d["platform"] == C.PLATFORM_RNA002]
    drawn = {s for _, _, units in PANELS for s in units}
    d = d[d["sample"].isin(drawn)].copy()
    # the curve table carries no dataset_group column; take it from the registry
    d["dataset_group"] = [
        C.SAMPLES_BY_NAME[s].dataset_group if s in C.SAMPLES_BY_NAME else ""
        for s in d["sample"]]
    for c in ("hit_rate", "localization_accuracy", "exact_rate", "recall",
              "window", "n_calls_in_universe"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d


def _bootstrap_ci(values: np.ndarray, *, n_boot: int = 4000, alpha: float = 0.05,
                  seed: int = 20260921) -> tuple[float, float]:
    """Percentile bootstrap 95% CI of the mean over the independent units.

    E7 / R1-8 ask for confidence intervals rather than point estimates, so the
    window-sweep curves carry a CI band next to the unit-range band.  The
    resampling unit is the sequencing unit (never a pooled call set); mouse
    columns have one unit per study, so no CI is drawn there (``nan``).
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    draws = rng.choice(values, size=(n_boot, values.size), replace=True).mean(axis=1)
    lo, hi = np.percentile(draws, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi)


def summarise(curve: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for species, group, units in PANELS:
        for tool in TOOL_ORDER:
            for w in WINDOWS:
                sub = curve[(curve["species"] == species) & (curve["tool"] == tool)
                            & (curve["window"] == w)]
                sub = sub[sub["sample"].isin(units)]
                for metric in ("hit_rate", "localization_accuracy"):
                    v = sub[metric].dropna()
                    ci_lo, ci_hi = _bootstrap_ci(v)
                    rows.append({
                        "species": species, "dataset_group": group, "tool": tool,
                        "window": w, "metric": metric,
                        "n_units": int(v.size), "n_units_total": len(units),
                        "value_mean": float(v.mean()) if v.size else np.nan,
                        "value_min": float(v.min()) if v.size else np.nan,
                        "value_max": float(v.max()) if v.size else np.nan,
                        "value_ci_lo": ci_lo, "value_ci_hi": ci_hi,
                        "value_each": "|".join(f"{x:.4f}" for x in v),
                    })
    return pd.DataFrame(rows)


def _draw(ax, summary: pd.DataFrame, species: str, metric: str,
          ylabel: str) -> None:
    for tool in TOOL_ORDER:
        s = (summary[(summary["species"] == species)
                     & (summary["tool"] == tool) & (summary["metric"] == metric)]
             .set_index("window").reindex(WINDOWS))
        y = s["value_mean"].to_numpy(dtype=float)
        lo = s["value_min"].to_numpy(dtype=float)
        hi = s["value_max"].to_numpy(dtype=float)
        ci_lo = s["value_ci_lo"].to_numpy(dtype=float)
        ci_hi = s["value_ci_hi"].to_numpy(dtype=float)
        if not np.isfinite(y).any():
            continue
        # outer band: the range of the independent units; inner band: the
        # bootstrap 95% CI of their mean (E7 / R1-8: intervals, not point
        # estimates).  Single-unit mouse columns get neither.
        if np.isfinite(lo).sum() > 1:
            ax.fill_between(WINDOWS, lo, hi, color=TOOL_COLOR[tool], alpha=.10, lw=0)
        if np.isfinite(ci_lo).sum() > 1:
            ax.fill_between(WINDOWS, ci_lo, ci_hi, color=TOOL_COLOR[tool],
                            alpha=.28, lw=0)
        ax.plot(WINDOWS, y, "-o", color=TOOL_COLOR[tool], ms=4, lw=1.6,
                label=DISPLAY.get(tool, tool))
    ax.axvline(C.PRIMARY_WINDOW, color="0.35", lw=1.0, ls=(0, (3, 2)), zorder=1)
    # the scale must be set BEFORE the ticks, otherwise the symlog formatter
    # relabels "10" as 10^1 (and re-adds a minor label)
    ax.set_xscale("symlog", linthresh=1)
    ax.minorticks_off()
    ax.set_xlim(-0.15, 60)
    ax.set_xticks(WINDOWS)
    ax.set_xticklabels([str(w) for w in WINDOWS], fontsize=9)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Matching window (bp)")
    ax.set_ylabel(ylabel)


def figures(summary: pd.DataFrame) -> None:
    apply_style()
    for metric, ylab, stem_name in (
            # naming (user decision 2026-09-19): precision against the reference,
            # i.e. PPV -- the old "GLORI hit rate" label is retired
            ("hit_rate", f"PPV vs. GLORI ({C.PRIMARY_WINDOW} bp)",
             "Fig5C_window_hitrate"),
            ("localization_accuracy", "Exact single-nucleotide localisation",
             "Fig5D_exact_localization")):
        fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.9), constrained_layout=True)
        for ax, (species, group, units) in zip(axes, PANELS):
            _draw(ax, summary, species, metric, ylab)
            n = len(units)
            ax.set_title(f"{species} ({'n = 1 study' if n == 1 else f'n = {n}'})",
                         fontweight="bold", fontsize=15, pad=10)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="outside lower center", ncol=7,
                   frameon=False, fontsize=10)
        FIG.mkdir(parents=True, exist_ok=True)
        save(fig, str(FIG / stem_name))
        print("wrote", stem_name + ".pdf + .png")


def main() -> None:
    logger = setup_logger("37_fig5cd_window_sweep")
    TAB.mkdir(parents=True, exist_ok=True)

    curve = load_curve()
    cols = ["species", "dataset_group", "sample", "tool", "window",
            "hit_rate", "localization_accuracy", "exact_rate", "recall",
            "n_calls_in_universe"]
    write_table(curve[cols].sort_values(["species", "sample", "tool", "window"]),
                TAB / "fig5cd_window_sweep_per_replicate.tsv")
    summary = summarise(curve)
    write_table(summary.sort_values(["species", "tool", "metric", "window"]),
                TAB / "fig5cd_window_sweep_summary.tsv")

    # ---- integrity checks ------------------------------------------------ #
    missing = sorted(set(TOOL_ORDER) - set(curve["tool"]))
    if missing:
        logger.warning("tools absent from the localisation curve: %s", ", ".join(missing))
    for species, group, units in PANELS:
        for tool in TOOL_ORDER:
            s = summary[(summary["species"] == species) & (summary["tool"] == tool)
                        & (summary["metric"] == "localization_accuracy")]
            v = s.set_index("window")["value_mean"].reindex(WINDOWS)
            if v.notna().sum() < 2:
                continue
            # exact localisation can only fall as the window widens (tp(w) grows),
            # so any *increase* beyond float noise is a real inconsistency
            if (np.diff(v.dropna().to_numpy(dtype=float)) > 1e-9).any():
                raise AssertionError(
                    f"{species}/{tool}: localization_accuracy rises with the window")
        logger.info("%-12s drawn units=%s", species, ",".join(units))

    # cross-check against the mean-of-tools table 27 uses (same source table)
    for species in [p[0] for p in PANELS]:
        for w in (2, 50):
            s = summary[(summary["species"] == species) & (summary["window"] == w)
                        & (summary["metric"] == "hit_rate")]
            logger.info("%-12s w=%-2d mean hit rate over tools = %.4f (n_tools=%d)",
                        species, w, s["value_mean"].mean(),
                        int(s["value_mean"].notna().sum()))

    figures(summary)
    logger.info("output: %s", OUT)


if __name__ == "__main__":
    main()
