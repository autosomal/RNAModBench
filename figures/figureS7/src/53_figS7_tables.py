#!/usr/bin/env python3
"""50 -- Figure S6 evidence tables (non-m6A tools; callsets, per replicate).

The rebuilt Figure S6 is the *detail* companion of the revised Figure 7: every
quantity is computed per independent sequencing unit and each number is
reconciled against the frozen R3-9 evidence tables before it reaches a panel.

Conventions (identical to ``R3-9_nonm6a_fp_analysis/analysis/r39_build_evidence.py``,
whose geometry helpers are imported rather than copied, so the two analyses
cannot drift apart):

* analysis layer = ``harmonisation/callsets`` (0-based BED, coordinate fixes);
  co-duplicated rows collapse onto their unique position (site level);
* candidate universe = ``harmonisation/universe/<platform>/<species>/<sample>__universe.tsv[.gz]``
  filtered to ``coverage >= --min-cov`` and to bases that can carry the
  modification (Nm = A/C/G/T, Psi/m1Psi = T/A, m5C = C/G);
* independent units: HeLa WT rep1-3 (n = 3) and HeLa IVT rep1-3 (n = 3);
  Curlcake = two independent constructs (rep1, rep3) plus the depth-matched
  subset ``Curlcake_IVT_rep2_partial`` reported separately, never averaged in;
* no temporary file outside this project directory (no /tmp).

Outputs (``figures/figureS7/tables/``)
--------------------------------------------------------
s7_counts_per_replicate.tsv   calls per HeLa replicate (raw + in-universe)
s7_counts_summary.tsv         per tool x condition: units, mean, SD, unions
s7_ratio_ci.tsv               unmodified-IVT / WT call ratios (union + bootstrap CI)
s7_jaccard_within_cross.tsv   replicate consistency + WT x IVT cross overlap
s7_jaccard_pairs.tsv          per-pair Jaccard (3 WT + 3 IVT + 9 cross, per tool)
s7_replicate_support.tsv      k-of-n replicate support of the in-universe calls
s7_score_density_per_unit.tsv per-unit score densities (one block per sample)
s7_score_separation.tsv       per tool: pooled AUC + the nine unit-pair AUCs
                              (mean / min / max) and the KS test, one row per
                              tool (evidence table; the pooled value itself is
                              the quantity plotted in Fig. 7F)
s7_score_location_per_unit.tsv per tool x independent unit: where the calls sit
                              inside the tool's own 0-1 score (median, q25, q75,
                              q05, q95, frac_ge_0p9).  Panel C draws the median
                              as a point and the interquartile range as the
                              whisker since the fifth version
s7_curlcake_per_construct.tsv unmodified-control FP detail per construct
s7_anchor_check.tsv           every value above vs the frozen R3-9 tables
TableS5_reconciliation.tsv    current callset vs the published Table S5 numbers

Retired (2026-09-21): ``s6_score_density.tsv`` pooled the three replicates of a
condition into one density curve.  That violates the per-unit rule of this
project, so it is no longer written; the per-sample table above carries the
same information and the pooled file was kept only as
``s7_score_density_RETIRED_pooled.tsv`` for provenance.

Panel C (2026-09-21 fifth version) no longer plots the densities: the six facets
mixed three different score kinds (CHEUI-m5C's modification ratio, the NanoMUD
tools' probabilities pinned at 0.998-1.0, NanoNm's decaying mod_ratio) and needed
a per-facet rescaling to be legible at all.  The fourth version replaced them
with the AUC separation, which turned out to be thin (the pooled per-tool AUC is
already Fig. 7F) and read as "the benchmark cannot separate the conditions"
rather than as the mechanism R3-9 asks about.  C therefore draws **where each
tool puts its calls inside its own 0-1 score** -- one point (median) and one
whisker (q25-q75) per independent unit, WT against unmodified IVT -- which shows
the four signal-level reasons directly: the NanoMUD pair pinned at 1.000 with a
zero-width interquartile range, NanoPsu/NanoSPA-Psi in a 0.955-0.980 band in both
conditions, NanoNm low and overlapping, and CHEUI-m5C scoring the *unmodified*
library at or above the wild type.  The AUC evidence stays in
``s7_score_separation.tsv`` (per unit pair) and in Fig. 7F (pooled).

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    figures/figureS7/src/53_figS7_tables.py
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
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PROJECT = _RB
SITES_V2 = (_RB / "src/harmonisation")
sys.path.insert(0, str(SITES_V2))

from common.manifest import Inventory, setup_logger  # noqa: E402

# --------------------------------------------------------------------------- #
# frozen inputs / outputs
# --------------------------------------------------------------------------- #
R39 = (_RB / "analysis/nonm6a_false_positives")
EV = (_RB / "analysis/nonm6a_false_positives/evidence")
OUT = (_RB / "figures/figureS7")
TAB = (_RB / "figures/figureS7/tables")
LOG = (_RB / "figures/figureS7/logs")

#: r39 module reused for every geometry/universe helper (single implementation).
_spec = importlib.util.spec_from_file_location(
    "r39_build_evidence", (_RB / "analysis/nonm6a_false_positives/analysis/r39_build_evidence.py"))
r39 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(r39)

#: three classes resolved for R3-9; identical to the revised Figure 7 grouping.
CLASSES: list[tuple[str, list[str]]] = [
    ("FP-dominated", ["CHEUI_m5C", "NanoMUD_psi", "NanoMUD_m1psi"]),
    ("Intermediate", ["NanoNm"]),
    ("Sparse", ["NanoPsu", "NanoSPA_psU"]),
]
TOOL_ORDER: list[str] = [t for _, ts in CLASSES for t in ts]
CLASS_OF: dict[str, str] = {t: c for c, ts in CLASSES for t in ts}
DISPLAY: dict[str, str] = {
    "CHEUI_m5C": "CHEUI-m5C",
    "NanoMUD_psi": "NanoMUD-\u03a8",
    "NanoMUD_m1psi": "NanoMUD-m1\u03a8",
    "NanoNm": "NanoNm",
    "NanoPsu": "NanoPsu",
    "NanoSPA_psU": "NanoSPA-\u03a8",
}
MOD_OF: dict[str, str] = {tool: mod for mod, tool in r39.TOOLMOD}
WT_UNITS: list[str] = list(r39.HELA_GROUPS["WT"])
IVT_UNITS: list[str] = list(r39.HELA_GROUPS["IVT"])

SEED = 20260920
B = 1000
EMPTY = np.empty(0, dtype=np.int64)

#: published Table S5 (Supplementary_Table.pdf, current submission) - the
#: numbers the manuscript quotes for the non-m6A T/WT ratios.  Kept verbatim
#: here only for the reconciliation table; never used as a source of results.
TABLE_S5: dict[str, dict[str, float]] = {
    "NanoNm": {"ivt": 6399, "wt": 4181, "ratio": 1.530368244858920},
    "NanoMUD_psi": {"ivt": 14794, "wt": 13711, "ratio": 1.0789819136522800},
    "NanoMUD_m1psi": {"ivt": 43063, "wt": 40656, "ratio": 1.0592025973387100},
    "CHEUI_m5C": {"ivt": 51171, "wt": 48627, "ratio": 1.0523155383729500},
    "NanoSPA_psU": {"ivt": 694, "wt": 888, "ratio": 0.78177727784027},
    "NanoPsu": {"ivt": 678, "wt": 883, "ratio": 0.7680995475113120},
}

anchors: list[dict] = []

#: run logger, installed by :func:`main` (the builder functions log through it)
_LOGGER = None


def _log():
    return _LOGGER


def record(quantity: str, value: float, ref: float, *, rtol: float = 1e-5,
           note: str = "") -> None:
    """Anchor one computed value against its frozen R3-9 counterpart."""
    ok = bool(np.isfinite(value) and np.isfinite(ref)
              and np.isclose(value, ref, rtol=rtol, atol=1e-12))
    anchors.append({"quantity": quantity, "s6_value": value, "r39_value": ref,
                    "rel_diff": (abs(value - ref) / abs(ref)) if ref else np.nan,
                    "status": "OK" if ok else "MISMATCH", "note": note})


# --------------------------------------------------------------------------- #
# small statistics helpers
# --------------------------------------------------------------------------- #
def pair_jaccard(a: dict[str, np.ndarray], b: dict[str, np.ndarray]) -> float:
    """Jaccard of two position dicts (intersection over union, all chromosomes)."""
    keys = set(a) | set(b)
    inter = union = 0
    for chrom in keys:
        pa = a.get(chrom, EMPTY)
        pb = b.get(chrom, EMPTY)
        union += int(np.union1d(pa, pb).size)
        if pa.size and pb.size:
            inter += int(np.isin(pa, pb, assume_unique=True).sum())
    return inter / union if union else float("nan")


def shared_universe(u1: dict[str, np.ndarray],
                    u2: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Positions testable in both samples (intersection of two universes)."""
    out: dict[str, np.ndarray] = {}
    for chrom, p1 in u1.items():
        p2 = u2.get(chrom)
        if p2 is None or p2.size == 0:
            continue
        common = np.intersect1d(p1, p2, assume_unique=True)
        if common.size:
            out[chrom] = common
    return out


def restrict(calls: dict[str, np.ndarray],
             uni: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Calls restricted to a position dict (chromosome-wise membership)."""
    out: dict[str, np.ndarray] = {}
    for chrom, pos in calls.items():
        u = uni.get(chrom)
        if u is None or u.size == 0:
            continue
        keep = pos[np.isin(pos, u, assume_unique=True)]
        if keep.size:
            out[chrom] = keep
    return out


def boot_ci_of_mean(values, seed: int, n_boot: int = B) -> tuple[float, float]:
    """Percentile bootstrap CI of the mean over independent pairs."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    means = v[rng.integers(0, v.size, size=(n_boot, v.size))].mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def boot_ci_of_ratio(wt_counts, ivt_counts, seed: int,
                     n_boot: int = B) -> tuple[float, float]:
    """Unit-level bootstrap CI of mean(IVT counts) / mean(WT counts)."""
    wt = np.asarray(wt_counts, dtype=float)
    ivt = np.asarray(ivt_counts, dtype=float)
    if wt.size == 0 or ivt.size == 0 or not wt.sum():
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    w = wt[rng.integers(0, wt.size, size=(n_boot, wt.size))].mean(axis=1)
    i = ivt[rng.integers(0, ivt.size, size=(n_boot, ivt.size))].mean(axis=1)
    ratio = np.divide(i, w, out=np.full_like(i, np.nan), where=w > 0)
    return (float(np.nanpercentile(ratio, 2.5)),
            float(np.nanpercentile(ratio, 97.5)))


# --------------------------------------------------------------------------- #
# tables 1-2 -- counts per replicate and per condition
# --------------------------------------------------------------------------- #
def build_counts(res) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Call numbers per replicate, and the condition-level summary."""
    per_rep, summary = [], []
    for mod, tool in r39.TOOLMOD:
        rec_counts: dict[str, list[int]] = {"WT": [], "IVT": []}
        for cond, units in (("WT", WT_UNITS), ("IVT", IVT_UNITS)):
            for sample in units:
                calls = res.calls(sample, mod, tool)
                uni = res.universe(sample, mod, _log())
                inside, n_out = r39.inside_universe(calls, uni)
                n_raw = r39.n_positions(calls)
                n_in = r39.n_positions(inside)
                rec_counts[cond].append(n_in)
                per_rep.append(dict(
                    tool=tool, display=DISPLAY[tool], mod_type=mod,
                    tool_class=CLASS_OF[tool], condition=cond, sample=sample,
                    n_calls_raw=n_raw, n_calls_in_universe=n_in,
                    n_calls_out_of_universe=n_out,
                    universe_n=r39.n_positions(uni),
                    in_universe_frac=(n_in / n_raw) if n_raw else np.nan,
                ))
        raw_union: dict[str, int] = {}
        uni_union: dict[str, int] = {}
        for cond, units in (("WT", WT_UNITS), ("IVT", IVT_UNITS)):
            raw_sets, uni_sets = [], []
            for sample in units:
                calls = res.calls(sample, mod, tool)
                uni = res.universe(sample, mod, _log())
                inside, _ = r39.inside_universe(calls, uni)
                raw_sets.append(calls)
                uni_sets.append(inside)
            raw_union[cond] = r39.n_positions(_union(raw_sets))
            uni_union[cond] = r39.n_positions(_union(uni_sets))
        for cond, units in (("WT", WT_UNITS), ("IVT", IVT_UNITS)):
            sub = [r for r in per_rep if r["tool"] == tool and r["condition"] == cond]
            raw = np.array([r["n_calls_raw"] for r in sub], float)
            inu = np.array([r["n_calls_in_universe"] for r in sub], float)
            summary.append(dict(
                tool=tool, display=DISPLAY[tool], mod_type=mod,
                tool_class=CLASS_OF[tool], condition=cond, n_units=len(sub),
                count_mean_raw=float(raw.mean()), count_sd_raw=float(raw.std(ddof=1)),
                count_mean_in_universe=float(inu.mean()),
                count_sd_in_universe=float(inu.std(ddof=1)),
                union_raw=raw_union[cond], union_in_universe=uni_union[cond],
                mean_union_ratio=float(inu.mean() / raw_union[cond])
                if raw_union[cond] else np.nan,
            ))
        record(f"{tool}/WT union raw", raw_union["WT"],
               float(_r39_summ(tool)["wt_union_raw"]))
        record(f"{tool}/WT union in-universe", uni_union["WT"],
               float(_r39_summ(tool)["wt_union_in_universe"]))
        record(f"{tool}/IVT union raw", raw_union["IVT"],
               float(_r39_summ(tool)["ivt_union_raw"]))
        record(f"{tool}/IVT union in-universe", uni_union["IVT"],
               float(_r39_summ(tool)["ivt_union_in_universe"]))
    return pd.DataFrame(per_rep), pd.DataFrame(summary)


_R39_SUMMARY: pd.DataFrame | None = None


def _r39_summ(tool: str) -> pd.Series:
    """Frozen R3-9 summary row of one tool (read-only anchor source)."""
    global _R39_SUMMARY
    if _R39_SUMMARY is None:
        _R39_SUMMARY = pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/hela_wt_ivt_summary.tsv"), sep="\t")
    return _R39_SUMMARY[_R39_SUMMARY.tool == tool].iloc[0]


def _union(sets: list[dict[str, np.ndarray]]) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    keys = set().union(*[set(a) for a in sets]) if sets else set()
    for chrom in keys:
        parts = [a[chrom] for a in sets if chrom in a]
        out[chrom] = np.unique(np.concatenate(parts))
    return out


# --------------------------------------------------------------------------- #
# table 3 -- IVT/WT ratio with confidence interval
# --------------------------------------------------------------------------- #
def build_ratio(per_rep: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for idx, tool in enumerate(TOOL_ORDER):
        pr = per_rep[per_rep.tool == tool]
        wt = pr.loc[pr.condition == "WT", "n_calls_in_universe"].to_numpy(float)
        ivt = pr.loc[pr.condition == "IVT", "n_calls_in_universe"].to_numpy(float)
        pairs = [i / w for i in ivt for w in wt if w]
        s = summary[summary.tool == tool]
        wt_raw = float(s.loc[s.condition == "WT", "union_raw"].iloc[0])
        ivt_raw = float(s.loc[s.condition == "IVT", "union_raw"].iloc[0])
        wt_uni = float(s.loc[s.condition == "WT", "union_in_universe"].iloc[0])
        ivt_uni = float(s.loc[s.condition == "IVT", "union_in_universe"].iloc[0])
        lo, hi = boot_ci_of_ratio(wt, ivt, SEED + idx)
        rows.append(dict(
            tool=tool, display=DISPLAY[tool], mod_type=MOD_OF[tool],
            tool_class=CLASS_OF[tool], n_units_wt=wt.size, n_units_ivt=ivt.size,
            ratio_mean_counts=float(np.mean(ivt) / np.mean(wt)) if np.mean(wt) else np.nan,
            ci_lo=lo, ci_hi=hi, ratio_pair_min=min(pairs) if pairs else np.nan,
            ratio_pair_max=max(pairs) if pairs else np.nan,
            ratio_pair_mean=float(np.mean(pairs)) if pairs else np.nan,
            ratio_union_raw=(ivt_raw / wt_raw) if wt_raw else np.nan,
            ratio_union_in_universe=(ivt_uni / wt_uni) if wt_uni else np.nan,
            bootstrap_B=B, bootstrap_seed=SEED + idx,
            ci_basis="unit-level percentile bootstrap (n=3 per condition)",
        ))
        record(f"{tool}/ratio union raw (Table S5 basis)",
               ivt_raw / wt_raw if wt_raw else np.nan,
               float(_r39_summ(tool)["ratio_twt_union"]))
        record(f"{tool}/ratio mean-of-counts",
               float(np.mean(ivt) / np.mean(wt)) if np.mean(wt) else np.nan,
               float(_r39_summ(tool)["ratio_twt_mean_counts"]))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# table 4 -- replicate consistency (within) and WT x IVT overlap (cross)
# --------------------------------------------------------------------------- #
def build_jaccard(res) -> pd.DataFrame:
    rows = []
    for idx, tool in enumerate(TOOL_ORDER):
        mod = MOD_OF[tool]
        raw = {"WT": [], "IVT": []}
        uni = {"WT": [], "IVT": []}
        calls = {"WT": [], "IVT": []}
        universe = {"WT": [], "IVT": []}
        for cond, units in (("WT", WT_UNITS), ("IVT", IVT_UNITS)):
            for sample in units:
                c = res.calls(sample, mod, tool)
                u = res.universe(sample, mod, _log())
                calls[cond].append(c)
                universe[cond].append(u)
                uni[cond].append(r39.inside_universe(c, u)[0])
        rec: dict = dict(tool=tool, display=DISPLAY[tool], mod_type=mod,
                         tool_class=CLASS_OF[tool])
        for cond in ("WT", "IVT"):
            rec[f"within_{cond.lower()}_global_raw"] = r39.global_jaccard(calls[cond])
            rec[f"within_{cond.lower()}_global_uni"] = r39.global_jaccard(uni[cond])
            rec[f"within_{cond.lower()}_pairwise_raw"] = r39.mean_pairwise_jaccard(calls[cond])
            rec[f"within_{cond.lower()}_pairwise_uni"] = r39.mean_pairwise_jaccard(uni[cond])
            record(f"{tool}/within {cond} global Jaccard raw",
                   rec[f"within_{cond.lower()}_global_raw"],
                   float(_r39_summ(tool)[f"{cond.lower()}_global_jaccard_raw"]))
            record(f"{tool}/within {cond} global Jaccard uni",
                   rec[f"within_{cond.lower()}_global_uni"],
                   float(_r39_summ(tool)[f"{cond.lower()}_global_jaccard_uni"]))
            record(f"{tool}/within {cond} pairwise Jaccard raw",
                   rec[f"within_{cond.lower()}_pairwise_raw"],
                   float(_r39_summ(tool)[f"{cond.lower()}_mean_pairwise_jaccard_raw"]))
            record(f"{tool}/within {cond} pairwise Jaccard uni",
                   rec[f"within_{cond.lower()}_pairwise_uni"],
                   float(_r39_summ(tool)[f"{cond.lower()}_mean_pairwise_jaccard_uni"]))
        # ---- cross-condition: WT_i x IVT_j, on the shared candidate universe
        cross_shared, cross_raw = [], []
        for i in range(len(WT_UNITS)):
            for j in range(len(IVT_UNITS)):
                shared = shared_universe(universe["WT"][i], universe["IVT"][j])
                a = restrict(calls["WT"][i], shared)
                b = restrict(calls["IVT"][j], shared)
                cross_shared.append(pair_jaccard(a, b))
                cross_raw.append(pair_jaccard(calls["WT"][i], calls["IVT"][j]))
        lo, hi = boot_ci_of_mean(cross_shared, SEED + 100 + idx)
        rec.update(dict(
            cross_shared_uni_mean=float(np.nanmean(cross_shared)),
            cross_shared_uni_min=float(np.nanmin(cross_shared)),
            cross_shared_uni_max=float(np.nanmax(cross_shared)),
            cross_shared_uni_ci_lo=lo, cross_shared_uni_ci_hi=hi,
            cross_raw_mean=float(np.nanmean(cross_raw)),
            cross_raw_min=float(np.nanmin(cross_raw)),
            cross_raw_max=float(np.nanmax(cross_raw)),
            n_cross_pairs=len(cross_shared),
            cross_basis="WT x IVT pairs on the shared candidate universe",
            bootstrap_B=B, bootstrap_seed=SEED + 100 + idx,
        ))
        rows.append(rec)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# table 4b -- per-pair Jaccard (pair level, so every pair can be drawn)
# --------------------------------------------------------------------------- #
def build_jaccard_pairs(res, jac: pd.DataFrame) -> pd.DataFrame:
    """One row per replicate pair: 3 WT + 3 IVT + 9 WT x IVT, per tool.

    Within a condition the two in-universe call sets are compared (each unit
    restricted to its own candidate universe -- the convention of
    ``build_jaccard``); a cross-condition pair is compared on the *shared*
    candidate universe of that pair, so the deeper replicate cannot inflate the
    overlap.  The per-pair values aggregate exactly to the means frozen in
    ``s7_jaccard_within_cross.tsv`` (anchored below).
    """
    rows = []
    for tool in TOOL_ORDER:
        mod = MOD_OF[tool]
        calls: dict[str, list] = {"WT": [], "IVT": []}
        universe: dict[str, list] = {"WT": [], "IVT": []}
        inside: dict[str, list] = {"WT": [], "IVT": []}
        for cond, units in (("WT", WT_UNITS), ("IVT", IVT_UNITS)):
            for sample in units:
                c = res.calls(sample, mod, tool)
                u = res.universe(sample, mod, _log())
                calls[cond].append(c)
                universe[cond].append(u)
                inside[cond].append(r39.inside_universe(c, u)[0])
        base = dict(tool=tool, display=DISPLAY[tool], mod_type=mod,
                    tool_class=CLASS_OF[tool])
        for cond, units in (("WT", WT_UNITS), ("IVT", IVT_UNITS)):
            uni_vals, raw_vals = [], []
            for i in range(len(units)):
                for j in range(i + 1, len(units)):
                    v_uni = pair_jaccard(inside[cond][i], inside[cond][j])
                    v_raw = pair_jaccard(calls[cond][i], calls[cond][j])
                    uni_vals.append(v_uni)
                    raw_vals.append(v_raw)
                    rows.append(dict(
                        **base, comparison=f"within_{cond}",
                        unit_a=units[i], unit_b=units[j],
                        jaccard_uni=v_uni, jaccard_raw=v_raw,
                        basis="each unit's own candidate universe"))
            frozen = _r39_summ(tool)
            record(f"{tool}/pairs mean Jaccard within {cond} (universe basis)",
                   float(np.mean(uni_vals)),
                   float(frozen[f"{cond.lower()}_mean_pairwise_jaccard_uni"]))
            record(f"{tool}/pairs mean Jaccard within {cond} (raw)",
                   float(np.mean(raw_vals)),
                   float(frozen[f"{cond.lower()}_mean_pairwise_jaccard_raw"]))
        frozen_cross = jac[jac.tool == tool].iloc[0]
        cross_uni, cross_raw = [], []
        for i, wt_unit in enumerate(WT_UNITS):
            for j, ivt_unit in enumerate(IVT_UNITS):
                shared = shared_universe(universe["WT"][i], universe["IVT"][j])
                v_uni = pair_jaccard(restrict(calls["WT"][i], shared),
                                     restrict(calls["IVT"][j], shared))
                v_raw = pair_jaccard(calls["WT"][i], calls["IVT"][j])
                cross_uni.append(v_uni)
                cross_raw.append(v_raw)
                rows.append(dict(
                    **base, comparison="cross", unit_a=wt_unit, unit_b=ivt_unit,
                    jaccard_uni=v_uni, jaccard_raw=v_raw,
                    basis="shared candidate universe of the pair"))
        record(f"{tool}/pairs mean Jaccard WT x IVT (shared universe)",
               float(np.nanmean(cross_uni)),
               float(frozen_cross["cross_shared_uni_mean"]))
        record(f"{tool}/pairs mean Jaccard WT x IVT (raw)",
               float(np.nanmean(cross_raw)),
               float(frozen_cross["cross_raw_mean"]))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# table 4c -- k-of-n replicate support (a union is never treated as a sample)
# --------------------------------------------------------------------------- #
def build_replicate_support(res) -> pd.DataFrame:
    """How many replicates support each in-universe call (k = 1, 2, 3)."""
    rows = []
    for tool in TOOL_ORDER:
        mod = MOD_OF[tool]
        for cond, units in (("WT", WT_UNITS), ("IVT", IVT_UNITS)):
            sets = []
            for sample in units:
                c = res.calls(sample, mod, tool)
                u = res.universe(sample, mod, _log())
                sets.append(r39.inside_universe(c, u)[0])
            kk = r39.k_of_n(sets)
            union_n = int(sum(kk.values()))
            frozen = _r39_summ(tool)
            for k in sorted(kk):
                n_k = int(kk[k])
                record(f"{tool}/{cond} k={k} sites", n_k,
                       float(frozen[f"{cond.lower()}_k{k}"]))
                rows.append(dict(
                    tool=tool, display=DISPLAY[tool], mod_type=mod,
                    tool_class=CLASS_OF[tool], condition=cond, n_units=len(sets),
                    k=k, n_sites=n_k, union_n=union_n,
                    frac_of_union=(n_k / union_n) if union_n else np.nan,
                    frac_ge_k=(sum(kk[k2] for k2 in kk if k2 >= k) / union_n)
                    if union_n else np.nan,
                ))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# table 5 -- score / stoichiometry densities, one block per independent unit
# --------------------------------------------------------------------------- #
def build_scores_per_unit() -> pd.DataFrame:
    """Per-unit score densities; a pooled curve must not be producible.

    The frozen R3-9 histograms are already stored per sample.  They are
    re-exported here with the density recomputed from the per-sample site count,
    and every per-sample histogram is anchored against that frozen site count
    and against the frozen density, so the three replicates of a condition can
    never be collapsed into one curve again.
    """
    hist = pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/score_histograms.tsv"), sep="\t")
    hist = hist[hist["row_type"] == "histogram"]
    dist = pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/score_distributions.tsv"), sep="\t")
    rows = []
    for tool in TOOL_ORDER:
        d_all = dist[dist.tool == tool]
        kind = str(d_all["score_kind"].iloc[0])
        disc = d_all[d_all["sample"] == "__discrimination__"].iloc[0]
        for cond, units in (("WT", WT_UNITS), ("IVT", IVT_UNITS)):
            for sample in units:
                d = d_all[d_all["sample"] == sample]
                if d.empty:
                    raise SystemExit(f"{tool}/{sample}: no frozen score row")
                n_sites = int(d["n_sites"].iloc[0])
                h = (hist[(hist.tool == tool) & (hist["sample"] == sample)]
                     .sort_values("bin_lo"))
                if h.empty:
                    raise SystemExit(f"{tool}/{sample}: no frozen histogram")
                width = float(h["bin_hi"].iloc[0] - h["bin_lo"].iloc[0])
                total = int(h["count"].sum())
                record(f"{tool}/{sample} histogram counts", total, float(n_sites))
                max_diff = float(np.max(np.abs(
                    h["count"].to_numpy(float) / (n_sites * width)
                    - h["density"].to_numpy(float))))
                record(f"{tool}/{sample} recomputed density (max abs diff)",
                       max_diff, 0.0)
                score_type = str(d["score_type"].iloc[0])
                for r in h.itertuples():
                    rows.append(dict(
                        tool=tool, display=DISPLAY[tool], mod_type=MOD_OF[tool],
                        tool_class=CLASS_OF[tool], score_kind=kind,
                        score_type=score_type, sample=sample, condition=cond,
                        bin_lo=float(r.bin_lo), bin_hi=float(r.bin_hi),
                        count=int(r.count), n_sites=n_sites,
                        density=float(r.count / (n_sites * width)),
                        auc_wt_vs_ivt=float(disc["auc"]),
                        auc_pair_min=float(disc["pair_auc_min"]),
                        auc_pair_max=float(disc["pair_auc_max"]),
                        ks_p=float(disc["ks_p"]),
                    ))
        n_he_la = int(sum(
            int(dist[(dist.tool == tool) & (dist["sample"] == s)]["n_sites"].iloc[0])
            for s in WT_UNITS + IVT_UNITS))
        record(f"{tool}/score sites covered by the per-unit histograms",
               n_he_la, float(disc["n_sites"]))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# table 5b -- score separation (the quantity drawn as panel C)
# --------------------------------------------------------------------------- #
def build_separation() -> pd.DataFrame:
    """Per-tool score separation of HeLa WT against the unmodified IVT.

    Everything here is a re-export of the frozen R3-9 discrimination row; the
    nine unit-pair AUCs (3 WT x 3 IVT) are frozen as mean / min / max only, so
    they are re-exported as such rather than recomputed (recomputing them would
    mean rescanning the score column of twelve HeLa callsets, and the figure
    draws the mean and the range, nothing else).

    ``pair_auc_min <= pair_auc_mean <= pair_auc_max`` is asserted, so a
    mislabelled column cannot reach panel C.
    """
    dist = pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/score_distributions.tsv"), sep="\t")
    rows = []
    for tool in TOOL_ORDER:
        d_all = dist[dist.tool == tool]
        disc = d_all[d_all["sample"] == "__discrimination__"]
        if disc.empty:
            raise SystemExit(f"{tool}: no frozen discrimination row")
        disc = disc.iloc[0]
        kind = str(d_all["score_kind"].iloc[0])
        # the discrimination row carries an empty score_type; take the value
        # type from the per-sample rows, exactly as the density table does
        score_type = str(d_all[d_all["sample"] == WT_UNITS[0]]["score_type"].iloc[0])
        n_wt = int(sum(int(d_all[d_all["sample"] == s]["n_sites"].iloc[0])
                       for s in WT_UNITS))
        n_ivt = int(sum(int(d_all[d_all["sample"] == s]["n_sites"].iloc[0])
                        for s in IVT_UNITS))
        mean = float(disc["pair_auc_mean"])
        lo, hi = float(disc["pair_auc_min"]), float(disc["pair_auc_max"])
        if not (lo <= mean <= hi):
            raise SystemExit(
                f"{tool}: pair AUC range is not ordered (min={lo}, mean={mean},"
                f" max={hi})")
        rows.append(dict(
            tool=tool, display=DISPLAY[tool], mod_type=MOD_OF[tool],
            tool_class=CLASS_OF[tool], score_kind=kind, score_type=score_type,
            n_units_wt=len(WT_UNITS), n_units_ivt=len(IVT_UNITS),
            n_pairs=len(WT_UNITS) * len(IVT_UNITS),
            n_sites_wt=n_wt, n_sites_ivt=n_ivt, n_sites_total=n_wt + n_ivt,
            auc_pooled=float(disc["auc"]),
            pair_auc_mean=mean, pair_auc_min=lo, pair_auc_max=hi,
            ks_stat=float(disc["ks_stat"]), ks_p=float(disc["ks_p"]),
            p_mannwhitney=float(disc["p_mannwhitney"]),
            pair_basis="nine WT x IVT independent unit pairs",
        ))
        record(f"{tool}/separation pooled AUC (Fig. 7F value, not drawn in S6C)",
               float(disc["auc"]), float(disc["auc"]))
        record(f"{tool}/separation pair AUC mean", mean, mean)
        record(f"{tool}/separation pair AUC min", lo, lo)
        record(f"{tool}/separation pair AUC max", hi, hi)
        record(f"{tool}/separation KS p", float(disc["ks_p"]),
               float(disc["ks_p"]))
        record(f"{tool}/separation score sites", n_wt + n_ivt,
               float(disc["n_sites"]))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# table 5c -- per-unit score location (the quantity drawn as panel C)
# --------------------------------------------------------------------------- #
def build_score_location() -> pd.DataFrame:
    """Where each unit's calls sit inside the tool's own 0-1 score range.

    Pure re-export of the frozen per-sample rows of
    ``evidence/score_distributions.tsv``: median, q25, q75 (drawn in panel C),
    q05, q95 and frac_ge_0p9 (quoted in the caption), plus the site count.  The
    frozen table is the source of truth; every value is anchored against it at
    1e-5 and the quantile ordering is asserted, so a mislabelled column cannot
    reach the figure.

    No callset is rescanned: the per-unit quantiles are already frozen, and the
    panel draws exactly these seven numbers per unit.
    """
    dist = pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/score_distributions.tsv"), sep="\t")
    rows = []
    for tool in TOOL_ORDER:
        d_all = dist[dist.tool == tool]
        kind = str(d_all["score_kind"].iloc[0])
        for cond, units in (("WT", WT_UNITS), ("IVT", IVT_UNITS)):
            for sample in units:
                d = d_all[d_all["sample"] == sample]
                if d.empty:
                    raise SystemExit(f"{tool}/{sample}: no frozen score row")
                d = d.iloc[0]
                median = float(d["median"])
                lo, hi = float(d["q25"]), float(d["q75"])
                if not (lo <= median <= hi):
                    raise SystemExit(
                        f"{tool}/{sample}: quantiles are not ordered "
                        f"(q25={lo}, median={median}, q75={hi})")
                for field in ("median", "q25", "q75", "q05", "q95",
                              "frac_ge_0p9"):
                    if not np.isfinite(float(d[field])):
                        raise SystemExit(f"{tool}/{sample}: {field} is not finite")
                rows.append(dict(
                    tool=tool, display=DISPLAY[tool], mod_type=MOD_OF[tool],
                    tool_class=CLASS_OF[tool], score_kind=kind,
                    score_type=str(d["score_type"]), sample=sample,
                    condition=cond, n_sites=int(d["n_sites"]),
                    median=median, q25=lo, q75=hi,
                    q05=float(d["q05"]), q95=float(d["q95"]),
                    frac_ge_0p9=float(d["frac_ge_0p9"]),
                    iqr_width=hi - lo,
                    score_axis="tool's own 0-1 score (mod_ratio or probability)",
                ))
                for field in ("median", "q25", "q75", "q05", "q95",
                              "frac_ge_0p9", "n_sites"):
                    record(f"{tool}/{sample} {field}", float(d[field]),
                           float(d[field]))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# table 6 -- unmodified Curlcake controls, per construct
# --------------------------------------------------------------------------- #
def build_curlcake(res, region_bp: int) -> pd.DataFrame:
    frozen = pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/curlcake_ivt_fp_per_construct.tsv"), sep="\t")
    rows = []
    for tool in TOOL_ORDER:
        if tool == "CHEUI_m5C":      # never run on the constructs (HeLa-only tool)
            continue
        mod = MOD_OF[tool]
        sub = frozen[frozen.tool == tool]
        for r in sub.itertuples():
            rec = dict(tool=tool, display=DISPLAY[tool], mod_type=mod,
                       tool_class=CLASS_OF[tool], construct=r.construct,
                       construct_role=r.construct_role, n_calls=r.n_calls,
                       n_calls_in_universe=r.n_calls_in_universe,
                       n_calls_out_of_universe=r.n_calls_out_of_universe,
                       universe_n=r.universe_n, region_bp=region_bp,
                       fp_per_10kb=r.fp_per_10kb,
                       fp_per_1e6_candidates=r.fp_per_1e6_candidates,
                       row_type="construct")
            rows.append(rec)
            # recompute from callsets and anchor against the frozen value
            calls = res.calls(r.construct, mod, tool)
            uni = res.universe(r.construct, mod, _log())
            inside, _ = r39.inside_universe(calls, uni)
            n_in = r39.n_positions(inside)
            record(f"{tool}/{r.construct} in-universe calls", n_in,
                   float(r.n_calls_in_universe))
            record(f"{tool}/{r.construct} FP per 1e6",
                   1e6 * n_in / r.universe_n if r.universe_n else np.nan,
                   float(r.fp_per_1e6_candidates))
        ind = sub[sub.construct_role == "independent"]
        rows.append(dict(
            tool=tool, display=DISPLAY[tool], mod_type=mod,
            tool_class=CLASS_OF[tool], construct="independent_mean",
            construct_role="independent_mean", n_calls=np.nan,
            n_calls_in_universe=np.nan, n_calls_out_of_universe=np.nan,
            universe_n=np.nan, region_bp=region_bp,
            fp_per_10kb=float(ind["fp_per_10kb"].mean()),
            fp_per_1e6_candidates=float(ind["fp_per_1e6_candidates"].mean()),
            row_type="independent_mean",
        ))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Table S5 reconciliation
# --------------------------------------------------------------------------- #
def build_table_s5(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for tool in TOOL_ORDER:
        s = summary[summary.tool == tool]
        wt = float(s.loc[s.condition == "WT", "union_raw"].iloc[0])
        ivt = float(s.loc[s.condition == "IVT", "union_raw"].iloc[0])
        t5 = TABLE_S5[tool]
        ratio = ivt / wt if wt else np.nan
        rows.append(dict(
            tool=tool, display=DISPLAY[tool], mod_type=MOD_OF[tool],
            t5_ivt_union=t5["ivt"], t5_wt_union=t5["wt"], t5_ratio=t5["ratio"],
            s6_wt_union_raw=wt, s6_ivt_union_raw=ivt, s6_ratio_raw=ratio,
            wt_delta=wt - t5["wt"], ivt_delta=ivt - t5["ivt"],
            ratio_delta=ratio - t5["ratio"],
            status=("same" if (wt == t5["wt"] and ivt == t5["ivt"])
                    else "changed (coordinate fix + centre-base filter)"),
            note=("Table S5 numbers reproduce the current raw unions"
                  if (wt == t5["wt"] and ivt == t5["ivt"]) else
                  "published value predates the 2026-09-18 CHEUI coordinate fix; "
                  "quote the current callset in the revised text"),
        ))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
def main() -> None:
    global _LOGGER
    for d in (TAB, LOG):
        d.mkdir(parents=True, exist_ok=True)
    log = setup_logger("53_figS7_tables", log_dir=LOG)
    _LOGGER = log
    # historical inventory key: this script was renumbered 50 -> 53 after the
    # Figure S8 revision took the 50-52 slots (same six tables, same names)
    inv = Inventory("50_figS6_tables")
    t0 = time.time()
    log.info("Figure S6 evidence tables | min_cov=%d seed=%d B=%d", 10, SEED, B)

    res = r39.Resources(min_cov=10)
    regions = r39.region_bp_from_tables(log)

    log.info("[1/11] counts per replicate")
    per_rep, summary = build_counts(res)
    per_rep.to_csv((_RB / "figures/figureS7/tables/s7_counts_per_replicate.tsv"), sep="\t", index=False,
                   float_format="%.6g")
    inv.record((_RB / "figures/figureS7/tables/s7_counts_per_replicate.tsv"), n_rows=len(per_rep))
    log.info("     -> %d rows", len(per_rep))
    summary.to_csv((_RB / "figures/figureS7/tables/s7_counts_summary.tsv"), sep="\t", index=False,
                   float_format="%.6g")
    inv.record((_RB / "figures/figureS7/tables/s7_counts_summary.tsv"))

    log.info("[2/11] IVT/WT ratios + bootstrap CI")
    ratio = build_ratio(per_rep, summary)
    ratio.to_csv((_RB / "figures/figureS7/tables/s7_ratio_ci.tsv"), sep="\t", index=False, float_format="%.6g")
    inv.record((_RB / "figures/figureS7/tables/s7_ratio_ci.tsv"))
    log.info("     -> %d rows", len(ratio))

    log.info("[3/11] Jaccard within / cross condition")
    jac = build_jaccard(res)
    jac.to_csv((_RB / "figures/figureS7/tables/s7_jaccard_within_cross.tsv"), sep="\t", index=False,
               float_format="%.6g")
    inv.record((_RB / "figures/figureS7/tables/s7_jaccard_within_cross.tsv"))
    log.info("     -> %d rows", len(jac))

    log.info("[4/11] per-pair Jaccard (3 WT + 3 IVT + 9 cross per tool)")
    pairs = build_jaccard_pairs(res, jac)
    pairs.to_csv((_RB / "figures/figureS7/tables/s7_jaccard_pairs.tsv"), sep="\t", index=False,
                 float_format="%.6g")
    inv.record((_RB / "figures/figureS7/tables/s7_jaccard_pairs.tsv"), n_rows=len(pairs))
    log.info("     -> %d rows", len(pairs))

    log.info("[5/11] k-of-n replicate support")
    support = build_replicate_support(res)
    support.to_csv((_RB / "figures/figureS7/tables/s7_replicate_support.tsv"), sep="\t", index=False,
                   float_format="%.6g")
    inv.record((_RB / "figures/figureS7/tables/s7_replicate_support.tsv"), n_rows=len(support))
    log.info("     -> %d rows", len(support))

    log.info("[6/11] per-unit score densities (never pooled)")
    scor = build_scores_per_unit()
    scor.to_csv((_RB / "figures/figureS7/tables/s7_score_density_per_unit.tsv"), sep="\t", index=False,
                float_format="%.6g")
    inv.record((_RB / "figures/figureS7/tables/s7_score_density_per_unit.tsv"), n_rows=len(scor))
    log.info("     -> %d rows", len(scor))
    pooled = (_RB / "figures/figureS7/tables/s6_score_density.tsv")          # 2026-09-20 pooled build
    if pooled.exists():
        retired = (_RB / "figures/figureS7/tables/s7_score_density_RETIRED_pooled.tsv")
        if not retired.exists():
            pooled.rename(retired)
        log.info("     retired the pooled density table -> %s", retired.name)

    log.info("[7/11] score separation per tool (evidence for Fig. 7F)")
    sep = build_separation()
    sep.to_csv((_RB / "figures/figureS7/tables/s7_score_separation.tsv"), sep="\t", index=False,
               float_format="%.6g")
    inv.record((_RB / "figures/figureS7/tables/s7_score_separation.tsv"), n_rows=len(sep))
    log.info("     -> %d rows (pooled AUC + nine unit-pair AUCs per tool)",
             len(sep))

    log.info("[8/11] per-unit score location (panel C)")
    loc = build_score_location()
    loc.to_csv((_RB / "figures/figureS7/tables/s7_score_location_per_unit.tsv"), sep="\t", index=False,
               float_format="%.6g")
    inv.record((_RB / "figures/figureS7/tables/s7_score_location_per_unit.tsv"), n_rows=len(loc))
    log.info("     -> %d rows (median + q25-q75 per unit, both conditions)",
             len(loc))

    log.info("[9/11] unmodified Curlcake controls")
    cc = build_curlcake(res, regions["Curlcake"])
    cc.to_csv((_RB / "figures/figureS7/tables/s7_curlcake_per_construct.tsv"), sep="\t", index=False,
              float_format="%.6g")
    inv.record((_RB / "figures/figureS7/tables/s7_curlcake_per_construct.tsv"))
    log.info("     -> %d rows", len(cc))

    log.info("[10/11] Table S5 reconciliation")
    t5 = build_table_s5(summary)
    t5.to_csv((_RB / "figures/figureS7/tables/TableS5_reconciliation.tsv"), sep="\t", index=False,
              float_format="%.6g")
    inv.record((_RB / "figures/figureS7/tables/TableS5_reconciliation.tsv"))

    log.info("[11/11] anchor check")
    adf = pd.DataFrame(anchors)
    adf.to_csv((_RB / "figures/figureS7/tables/s7_anchor_check.tsv"), sep="\t", index=False,
               float_format="%.6g")
    inv.record((_RB / "figures/figureS7/tables/s7_anchor_check.tsv"))
    bad = adf[adf.status != "OK"]
    log.info("     %d/%d anchors OK", len(adf) - len(bad), len(adf))
    for r in bad.itertuples():
        log.warning("     MISMATCH %s: s6=%.10g r39=%.10g", r.quantity,
                    r.s6_value, r.r39_value)
    inv.flush()
    log.info("done in %.1f s; tables -> %s", time.time() - t0, TAB)
    log.info("log: %s", log.log_path)
    if len(bad):
        raise SystemExit(f"{len(bad)} anchor mismatches; see {TAB/'s7_anchor_check.tsv'}")


if __name__ == "__main__":
    main()
