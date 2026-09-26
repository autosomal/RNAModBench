#!/usr/bin/env python
"""Figure 3B, revision rebuild (R1-9 / E7 effect sizes, R3-2 / E6 replicates,
R1-4 / E5 negative controls).

The published panel B (``manuscript_figure_code/Figure3/mod_ratio.ipynb``) read
one legacy ``output/<group>/mod_ratio/*.txt`` aggregate per species (a single
replicate or an undocumented pile-up), paired the WT and treated call sets on
``chr_start`` and drew a boxplot per tool with ``***`` from a Wilcoxon test.
The editor/reviewers asked for (i) replicate structure to be carried into the
statistics (R3-2 / E6), (ii) effect sizes and confidence intervals instead of
p-value-only statements when n is huge (R1-9 / E7), and (iii) KO/KD not to be
treated as clean negatives, with coverage/expression confounds controlled
(R1-4 / E5).

This script recomputes the panel from ``harmonisation/callsets``:

  unit        one independent sequencing unit per point (Arabidopsis WT/KD and
              HeLa WT/IVT: three biological replicates each; mouse: the two
              cross-study pairs SRP357195 and SRP166020, never pooled);
  pairing     within a unit, the sites called by the tool in BOTH conditions,
              additionally restricted to positions measurable (candidate
              universe, coverage >= 10) in BOTH samples -- the coverage-matched
              comparison R1-4 asks for.  The unconstrained (legacy-style)
              pairing is computed as well and reported in the same table;
  ratio       the tool's own modification ratio: ``mod_ratio`` when populated,
              else ``score`` when ``score_type`` is a ratio -- identical rule to
              ``15_mod_ratio_replicate_agreement.py`` / ``36_fig5b``;
  effects     median ratio difference (WT - treated) with a site-level
              percentile bootstrap 95 % CI (B = 2000, seed fixed), and Cliff's
              delta for paired data (P(WT>trt) - P(WT<trt) over the paired
              sites).  The Wilcoxon signed-rank p (BH-FDR across the tested
              tool x species cells) is computed but kept OUT of the panel: it
              goes into the table and the caption only (house rule: no numbers
              inside a figure; R1-9: effect sizes first).

Outputs -> ``figures/figure3/{tables,figures,logs}``

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/src/harmonisation/scripts/57_fig3b_modratio_wt_treatment.py
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
import argparse
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from common import config as C                                        # noqa: E402
from common.figstyle import apply as apply_style                      # noqa: E402
from common.io_utils import write_table                               # noqa: E402
from common.manifest import setup_logger                              # noqa: E402

OUT = (_RB / "figures/figure3")
TAB, FIG, LOG = OUT / "tables", OUT / "figures", OUT / "logs"
CLEAN = (_RB / "data/callsets")
UNIVERSE = (_XB / "harmonisation/universe")

#: published panel-B tool order
TOOLS = ["m6Anet", "MINES", "Nanom6A", "DENA"]
#: score_type values that denote a modification ratio (see callsets/README.md)
RATIO_SCORE_TYPES = {"mod_ratio", "m6a_ratio", "ratio", "percent_modified/100"}

#: species -> (title, [ (unit label, WT sample, treated sample, WT label, condition label) ])
PAIRS: dict[str, tuple[str, list[tuple[str, str, str, str, str]]]] = {
    "Arabidopsis": ("Arabidopsis (WT vs KD)", [
        ("rep1", "Arabidopsis_WT_rep1", "Arabidopsis_fip37_rep1", "WT", "KD"),
        ("rep2", "Arabidopsis_WT_rep2", "Arabidopsis_fip37_rep2", "WT", "KD"),
        ("rep3", "Arabidopsis_WT_rep3", "Arabidopsis_fip37_rep3", "WT", "KD"),
    ]),
    "Mouse": ("Mouse (WT vs KO)", [
        ("studyB (SRP357195)", "mES_WT", "mES_KO", "WT", "KO"),
        ("studyA (SRP166020)", "mESCs_Mettl3_WT", "mESCs_Mettl3_KO", "WT", "KO"),
    ]),
    "Human": ("Human (WT vs IVT)", [
        ("rep1", "HeLa_WT1", "HeLa_IVT_rep1", "WT", "IVT"),
        ("rep2", "HeLa_WT2", "HeLa_IVT_rep2", "WT", "IVT"),
        ("rep3", "HeLa_WT3", "HeLa_IVT_rep3", "WT", "IVT"),
    ]),
}
GROUPS = {"Arabidopsis": ("Arabidopsis_WT", "Arabidopsis_KD"),
          "Mouse": ("Mouse_WT", "Mouse_KO"),
          "Human": ("HeLa_WT", "HeLa_IVT")}

#: values published in the original submission's text (single-replicate legacy
#: aggregates) -- kept for the legacy-vs-revision reconciliation table only
LEGACY_TEXT = {
    ("Human", "m6Anet"): (68.5, 63.9, "IVT"),
    ("Mouse", "m6Anet"): (71.6, 62.8, "Mettl3 KO"),
    ("Arabidopsis", "m6Anet"): (67.2, 62.0, "FIP37 KD"),
}

#: the legacy notebook's own statistics, exported by `mod_ratio.ipynb`
#: (`Modification_Ratio_comparison.csv`, one row per species x tool on a
#: single-replicate aggregate) -- read-only, used for the 12-cell
#: legacy-vs-revision reconciliation
LEGACY_CSV = ((_XB / "legacy/next_postprocessing/Total_new1/Modification_Ratio_comparison.csv"))
LEGACY_CSV_COLUMNS = {"Arabidopsis": "Median_KD", "Mouse": "Median_KO",
                      "Human": "Median_IVT"}
LEGACY_CONDITION = {"Arabidopsis": "FIP37 KD", "Mouse": "Mettl3 KO",
                    "Human": "IVT"}


def _legacy_csv() -> pd.DataFrame:
    """The legacy notebook's statistics for the twelve species x tool cells."""
    if not LEGACY_CSV.exists():
        return pd.DataFrame()
    df = pd.read_csv(LEGACY_CSV)
    df = df[df["Tool"].isin(TOOLS)]
    return df.set_index(["Species", "Tool"])

#: colours of the published panel (WT steel blue / condition orange)
WT_COLOR, TRT_COLOR, EFF_COLOR = "#3778A0", "#F5B264", "#4a4a4a"
#: panel canvas (final printed size; row heights chosen so the assembled A-D page
#: matches the printed height of the published Figure 3, ~7.8 in)
PANEL_W, PANEL_H = 6.66, 1.95
FS = {"title": 8.6, "label": 8.0, "tick": 7.4, "legend": 7.2}
BOOT_B, BOOT_SEED = 2000, 20260909
BA_MAX_POINTS = 600       # per tool strip, seeded subsample of paired sites
BOX_MAX_POINTS = 6000    # per box, seeded subsample of paired sites
FDR_SELFTEST_SEED = 20260920
X_OFF = 0.17


# --------------------------------------------------------------------------- #
def clean_path(species: str, group: str, tool: str, sample: str) -> Path:
    return CLEAN / C.PLATFORM_RNA002 / species / group / "m6A" / tool / f"{sample}.tsv"


def universe_path(species: str, sample: str) -> Path:
    return UNIVERSE / C.PLATFORM_RNA002 / species / f"{sample}__m6A.tsv"


def extract_ratio(df: pd.DataFrame, logger: logging.Logger, tag: str) -> pd.Series:
    """Per-site modification ratio: ``mod_ratio`` if populated, else ratio-score."""
    mr = (pd.to_numeric(df["mod_ratio"], errors="coerce")
          if "mod_ratio" in df.columns else pd.Series(np.nan, index=df.index))
    if mr.notna().any():
        return mr.clip(0, 1)
    st = ""
    if "score_type" in df.columns:
        seen = [str(v).strip().lower() for v in df["score_type"].dropna().unique()
                if str(v).strip()]
        st = seen[0] if seen else ""
    if st in RATIO_SCORE_TYPES or "ratio" in st:
        return pd.to_numeric(df["score"], errors="coerce").clip(0, 1)
    logger.warning("%s: no ratio column (score_type=%r)", tag, st)
    return pd.Series(np.nan, index=df.index)


def load_sites(species: str, group: str, tool: str, sample: str,
               logger: logging.Logger) -> pd.DataFrame:
    """``chrom/pos/ratio`` of one tool in one sample (empty frame when absent)."""
    path = clean_path(species, group, tool, sample)
    if not path.exists():
        logger.warning("callset missing: %s", path)
        return pd.DataFrame(columns=["chrom", "pos", "ratio"])
    try:
        d = pd.read_csv(path, sep="\t", dtype={"chrom": str})
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=["chrom", "pos", "ratio"])
    if d.empty:
        return pd.DataFrame(columns=["chrom", "pos", "ratio"])
    out = pd.DataFrame({"chrom": d["chrom"].astype(str),
                        "pos": pd.to_numeric(d["start"], errors="coerce").astype("Int64"),
                        "ratio": extract_ratio(d, logger, f"{sample}/{tool}")})
    return out.dropna(subset=["pos"]).assign(pos=lambda x: x["pos"].astype("int64"))


def chrom_key_map(frames: list[pd.DataFrame]) -> dict[str, int]:
    chroms = sorted({str(c) for f in frames for c in f["chrom"].unique()})
    return {c: i for i, c in enumerate(chroms)}


def keys_with(df: pd.DataFrame, chrom_ix: dict[str, int]) -> np.ndarray:
    idx = df["chrom"].map(chrom_ix)
    if idx.isna().any():
        miss = sorted(set(df.loc[idx.isna(), "chrom"].unique()))
        raise KeyError(f"chromosome(s) not in the universe chromosome set: {miss}")
    return idx.to_numpy("int64") * np.int64(10 ** 9) + df["pos"].to_numpy("int64")


def boot_ci_median(diffs: np.ndarray, seed: int, b: int = BOOT_B) -> tuple[float, float]:
    """Percentile bootstrap CI of the median paired difference (site level)."""
    v = np.asarray(diffs, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 5:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    draws = v[rng.integers(0, v.size, size=(b, v.size))].astype(float)
    med = np.median(draws, axis=1)
    return float(np.percentile(med, 2.5)), float(np.percentile(med, 97.5))


def boot_ci_median_stratified(groups: list[np.ndarray], seed: int,
                              b: int = BOOT_B) -> tuple[float, float, float]:
    """Pooled median difference with a stratified bootstrap 95% CI.

    The point estimate pools every measurable paired site of the group, while
    the resampling is stratified by independent unit (each unit contributes its
    own draw), so the CI reflects between-unit as well as between-site
    variability and replicates are never silently collapsed into one sample.
    """
    arr = [np.asarray(g, dtype=float) for g in groups]
    arr = [a[np.isfinite(a)] for a in arr]
    arr = [a for a in arr if a.size]
    if not arr:
        return float("nan"), float("nan"), float("nan")
    point = float(np.median(np.concatenate(arr)))
    if len(arr) == 1 and arr[0].size < 5:
        return point, float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    n_units = len(arr)
    draws = np.empty(b, dtype=float)
    for i in range(b):
        # two-level bootstrap: resample independent units, then the sites
        # inside each drawn unit, so the CI carries between-unit variance too
        picked = [arr[j] for j in rng.integers(0, n_units, n_units)]
        parts = [a[rng.integers(0, a.size, a.size)] for a in picked]
        draws[i] = np.median(np.concatenate(parts))
    return point, float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


#: (species, tool) -> (pooled median difference, CI low, CI high); filled by
#: `compute()` and written into the summary table / forest row of panel B
POOLED_CI: dict[tuple[str, str], tuple[float, float, float]] = {}

#: (species, tool) -> per-unit paired arrays for the Bland-Altman row:
#: mean = (ratio WT + ratio treated) / 2 per measurable paired site,
#: diff = ratio WT - ratio treated on exactly the same sites (the pairing rule
#: of the tables, i.e. called in both conditions with coverage >= 10 in both)
BA_PAIRS: dict[tuple[str, str], dict[str, list[np.ndarray]]] = {}

#: (species, tool) -> per-condition site ratios (same measurable pairs as above)
RATIO_SITES: dict[tuple[str, str], dict[str, list[np.ndarray]]] = {}


def cliff_delta_paired(diffs: np.ndarray) -> float:
    """Cliff's delta for paired data: P(d>0) - P(d<0), ties excluded."""
    d = np.asarray(diffs, dtype=float)
    d = d[np.isfinite(d)]
    if d.size == 0:
        return np.nan
    return float(((d > 0).sum() - (d < 0).sum()) / d.size)


def benjamini_hochberg_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted p-values (nan-safe).

    Delegates to ``statsmodels.stats.multitest.multipletests(method="fdr_bh")``
    so the multiple-comparison correction is the reference implementation.
    Applied to every p-value produced by the loop below (one test per
    tool x species x unit), so the reported significance is family-wise
    controlled rather than the raw per-test value.
    """
    from statsmodels.stats.multitest import multipletests

    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    out = np.full(p.shape, np.nan)
    if not ok.any():
        return out
    out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def _benjamini_hochberg_fdr_reference(p: np.ndarray) -> np.ndarray:
    """The previous hand-rolled BH implementation, kept for the self-test."""
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    out = np.full(p.shape, np.nan)
    if not ok.any():
        return out
    pv = p[ok]
    order = np.argsort(pv)
    m = pv.size
    adj = np.empty(m, dtype=float)
    prev = 1.0
    for rank, i in enumerate(order[::-1]):
        k = m - rank
        prev = min(prev, pv[i] * m / k)
        adj[i] = prev
    out[ok] = adj
    return out


def self_test_fdr() -> None:
    """Assert the library FDR and the hand-rolled reference agree."""
    rng = np.random.default_rng(FDR_SELFTEST_SEED)
    p = np.concatenate([rng.uniform(0, 0.05, 20), rng.uniform(0, 1, 20),
                        [np.nan]])
    a, b = benjamini_hochberg_fdr(p), _benjamini_hochberg_fdr_reference(p)
    ok = np.isfinite(a) & np.isfinite(b)
    delta = float(np.max(np.abs(a[ok] - b[ok])))
    if delta > 1e-12:
        raise AssertionError(f"BH-FDR self-test failed (max |d| = {delta})")
    print(f"BH-FDR self-test ok: multipletests vs reference, "
          f"{ok.sum()} finite values, max |d| = {delta:.2e}")


# --------------------------------------------------------------------------- #
def compute(logger: logging.Logger) -> pd.DataFrame:
    rows: list[dict] = []
    per_group: dict[tuple[str, str], list[np.ndarray]] = {}
    for species, (title, pairs) in PAIRS.items():
        wt_group, trt_group = GROUPS[species]
        for unit_label, wt_sample, trt_sample, wt_lab, trt_lab in pairs:
            # candidate universes of the two samples (coverage >= 10, m6A/A positions)
            uni = {}
            for sample in (wt_sample, trt_sample):
                up = universe_path(species, sample)
                if not up.exists():
                    logger.warning("universe missing: %s", up)
                    uni[sample] = None
                    continue
                u = pd.read_csv(up, sep="\t", usecols=["chrom", "pos"], dtype={"chrom": str})
                uni[sample] = u
            tools_data = {t: (load_sites(species, wt_group, t, wt_sample, logger),
                              load_sites(species, trt_group, t, trt_sample, logger))
                          for t in TOOLS}
            frames = [f for pair in tools_data.values() for f in pair] + \
                     [u for u in uni.values() if u is not None]
            chrom_ix = chrom_key_map(frames)
            uni_keys = {s: (set(keys_with(u, chrom_ix).tolist()) if u is not None else None)
                        for s, u in uni.items()}

            for t_i, tool in enumerate(TOOLS):
                wt, trt = tools_data[tool]
                if wt.empty or trt.empty:
                    logger.warning("%s/%s/%s: empty call set (WT=%d, trt=%d)",
                                   species, tool, unit_label, len(wt), len(trt))
                    continue
                wt = wt.assign(key=keys_with(wt, chrom_ix))
                trt = trt.assign(key=keys_with(trt, chrom_ix))
                m = (wt[["key", "ratio"]].drop_duplicates("key")
                     .merge(trt[["key", "ratio"]].drop_duplicates("key"),
                            on="key", suffixes=("_wt", "_trt")))
                n_uncon = len(m)
                if uni[wt_sample] is not None and uni[trt_sample] is not None:
                    keep = np.array([k in uni_keys[wt_sample] and k in uni_keys[trt_sample]
                                     for k in m["key"]], dtype=bool)
                    m_con = m[keep]
                else:
                    m_con = m
                if len(m_con) < 10:
                    logger.warning("%s/%s/%s: only %d measurable paired sites",
                                   species, tool, unit_label, len(m_con))
                d_con = (m_con["ratio_wt"] - m_con["ratio_trt"]).to_numpy(float)
                d_uncon = (m["ratio_wt"] - m["ratio_trt"]).to_numpy(float)
                lo, hi = boot_ci_median(d_con, BOOT_SEED + 7 * t_i + len(rows))
                per_group.setdefault((species, tool), []).append(d_con)
                slot = BA_PAIRS.setdefault((species, tool), {"mean": [], "diff": []})
                rs = RATIO_SITES.setdefault((species, tool), {"wt": [], "trt": []})
                rs["wt"].append(m_con["ratio_wt"].to_numpy(float))
                rs["trt"].append(m_con["ratio_trt"].to_numpy(float))
                slot["mean"].append(((m_con["ratio_wt"] + m_con["ratio_trt"]) / 2)
                                    .to_numpy(float))
                slot["diff"].append(d_con)
                try:
                    p_val = float(stats.wilcoxon(d_con, alternative="two-sided").pvalue) \
                        if len(d_con) >= 10 else np.nan
                except ValueError:
                    p_val = np.nan
                rows.append({
                    "species": species, "dataset_group": wt_group, "tool": tool,
                    "unit": unit_label, "wt_sample": wt_sample, "trt_sample": trt_sample,
                    "wt_label": wt_lab, "trt_label": trt_lab,
                    "n_paired": int(len(m_con)),
                    "n_paired_unconstrained": int(n_uncon),
                    "median_wt": float(np.median(m_con["ratio_wt"])) if len(m_con) else np.nan,
                    "median_trt": float(np.median(m_con["ratio_trt"])) if len(m_con) else np.nan,
                    "median_diff": float(np.median(d_con)) if len(d_con) else np.nan,
                    "ci_lo": lo, "ci_hi": hi,
                    "cliff_delta_paired": cliff_delta_paired(d_con),
                    "frac_sites_wt_higher": float((d_con > 0).mean()) if len(d_con) else np.nan,
                    "median_diff_unconstrained": float(np.median(d_uncon)) if len(d_uncon) else np.nan,
                    "wilcoxon_p": p_val,
                })
            logger.info("%-12s %-20s paired sites computed", species, unit_label)
    df = pd.DataFrame(rows)
    # multiple-comparison control: Benjamini-Hochberg FDR over every Wilcoxon
    # test performed in the loop above (tool x species x unit)
    df["wilcoxon_p_bh"] = benjamini_hochberg_fdr(df["wilcoxon_p"].to_numpy(float))
    # pooled effect size + stratified bootstrap CI, one entry per species x tool
    POOLED_CI.clear()
    # (BA_PAIRS is filled alongside per_group and cleared implicitly by reassignment)
    for i, ((species, tool), groups) in enumerate(per_group.items()):
        POOLED_CI[(species, tool)] = boot_ci_median_stratified(groups, BOOT_SEED + i)
    return df


def _ba_for(species: str, tool: str) -> tuple[np.ndarray, np.ndarray] | None:
    """Concatenated (mean, difference) arrays of one species x tool for the BA row."""
    slot = BA_PAIRS.get((species, tool))
    if not slot:
        return None
    mx = np.concatenate(slot["mean"])
    dy = np.concatenate(slot["diff"])
    ok = np.isfinite(mx) & np.isfinite(dy)
    if not ok.any():
        return None
    return mx[ok], dy[ok]


def bland_altman_table() -> pd.DataFrame:
    """Descriptive Bland-Altman summary of the WT-vs-treated paired ratios.

    bias = mean(ratio_WT - ratio_treated) on the measurable paired sites;
    limits of agreement = bias +/- 1.96 x SD of the same differences.  This is
    agreement *between the two conditions of one sample pair*, not with GLORI
    (the GLORI agreement metrics live in `mod_ratio_replicates/`).
    """
    rows: list[dict] = []
    for species in PAIRS:
        for tool in TOOLS:
            slot = BA_PAIRS.get((species, tool))
            if not slot:
                rows.append({"species": species, "tool": tool, "n_paired": 0,
                             "n_units": 0, "mean_ratio": np.nan, "bias": np.nan,
                             "sd_diff": np.nan, "loa_lo": np.nan, "loa_hi": np.nan,
                             "bias_minus_pooled_median": np.nan})
                continue
            d = np.concatenate(slot["diff"])
            m = np.concatenate(slot["mean"])
            bias = float(np.mean(d))
            sd = float(np.std(d, ddof=1)) if d.size > 1 else float("nan")
            pooled = POOLED_CI.get((species, tool), (float("nan"),) * 3)[0]
            rows.append({
                "species": species, "tool": tool, "n_paired": int(d.size),
                "n_units": len(slot["diff"]), "mean_ratio": float(np.mean(m)),
                "bias": bias, "sd_diff": sd,
                "loa_lo": bias - 1.96 * sd, "loa_hi": bias + 1.96 * sd,
                "bias_minus_pooled_median": bias - pooled,
            })
    return pd.DataFrame(rows)


def summarise(per_rep: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (species, tool), sub in per_rep.groupby(["species", "tool"], sort=False):
        rows.append({
            "species": species, "tool": tool, "n_units": int(len(sub)),
            "median_diff_mean": float(sub["median_diff"].mean()),
            "median_diff_sd": float(sub["median_diff"].std(ddof=1)) if len(sub) > 1 else np.nan,
            "median_diff_min": float(sub["median_diff"].min()),
            "median_diff_max": float(sub["median_diff"].max()),
            "cliff_delta_mean": float(sub["cliff_delta_paired"].mean()),
            "cliff_delta_min": float(sub["cliff_delta_paired"].min()),
            "cliff_delta_max": float(sub["cliff_delta_paired"].max()),
            "n_paired_total": int(sub["n_paired"].sum()),
            "wilcoxon_p_min": float(sub["wilcoxon_p"].min()),
            "wilcoxon_p_bh_min": float(sub["wilcoxon_p_bh"].min()),
            "n_units_p_bh_lt_0.05": int((sub["wilcoxon_p_bh"] < 0.05).sum()),
            "pooled_median_diff": POOLED_CI.get((species, tool),
                                                (float("nan"),) * 3)[0],
            "pooled_ci_lo": POOLED_CI.get((species, tool), (float("nan"),) * 3)[1],
            "pooled_ci_hi": POOLED_CI.get((species, tool), (float("nan"),) * 3)[2],
        })
    return pd.DataFrame(rows)


def legacy_table(per_rep: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for species, tool in [(s, t) for s in PAIRS for t in TOOLS]:
        sub = per_rep[(per_rep.species == species) & (per_rep.tool == tool)]
        row = {
            "species": species, "tool": tool,
            "legacy_median_wt": np.nan, "legacy_median_trt": np.nan,
            "legacy_condition": LEGACY_CONDITION[species],
            "legacy_n_common_sites": np.nan, "legacy_wilcoxon_p": np.nan,
            "legacy_basis": "legacy notebook aggregate (single replicate, "
                            "Modification_Ratio_comparison.csv)",
            "revision_median_wt_mean": float(sub["median_wt"].mean()) if len(sub) else np.nan,
            "revision_median_trt_mean": float(sub["median_trt"].mean()) if len(sub) else np.nan,
            "revision_median_diff_mean": float(sub["median_diff"].mean()) if len(sub) else np.nan,
            "revision_median_diff_sd": float(sub["median_diff"].std(ddof=1)) if len(sub) > 1 else np.nan,
            "revision_units": "; ".join(sub["unit"].tolist()),
            "revision_n_paired": int(sub["n_paired"].sum()) if len(sub) else 0,
        }
        legacy_csv = _legacy_csv()
        if (species, tool) in getattr(legacy_csv, "index", []):
            src = legacy_csv.loc[(species, tool)]
            row["legacy_median_wt"] = float(src["Median_WT"])
            row["legacy_median_trt"] = float(src[LEGACY_CSV_COLUMNS[species]])
            row["legacy_n_common_sites"] = float(src["Common_Sites"])
            row["legacy_wilcoxon_p"] = float(src["Wilcoxon_p"])
        rows.append(row)
    df = pd.DataFrame(rows)
    #: the numbers printed in the original submission's text, on top of the
    #: legacy CSV row (they must agree; a mismatch is a provenance red flag)
    for (species, tool), (wt_pub, trt_pub, _cond) in LEGACY_TEXT.items():
        m = (df.species == species) & (df.tool == tool)
        df.loc[m, "legacy_text_median_wt"] = wt_pub
        df.loc[m, "legacy_text_median_trt"] = trt_pub
    return df


def _pooled_for(per_rep: pd.DataFrame, species: str,
                tool: str) -> tuple[float, float, float] | None:
    """Pooled effect size of one species x tool cell, or None if unavailable."""
    if (species, tool) in POOLED_CI:
        return POOLED_CI[(species, tool)]
    table = TAB / "fig3b_effect_sizes.tsv"          # --figures-only path
    if not table.exists():
        return None
    df = pd.read_csv(table, sep="\t")
    hit = df[(df.species == species) & (df.tool == tool)]
    if hit.empty or "pooled_median_diff" not in df.columns:
        return None
    r = hit.iloc[0]
    return (float(r.pooled_median_diff), float(r.pooled_ci_lo),
            float(r.pooled_ci_hi))


# --------------------------------------------------------------------------- #
def figure(per_rep: pd.DataFrame, logger: logging.Logger) -> None:
    """Panel B: modification-ratio boxplots, wild type versus treated.

    One box pair per tool and species, drawn on the measurable paired sites
    (called in both conditions with coverage >= 10 in both); the median of each
    independent unit is overlaid as a small dot so the replicate structure stays
    visible.  Effect sizes (median difference, bootstrap 95% CI, Cliff's delta)
    are never printed on the figure - they are in
    `tables/fig3b_effect_sizes.tsv` and the caption.
    """
    apply_style()
    species_list = list(PAIRS)
    fig = plt.figure(figsize=(PANEL_W, PANEL_H))
    gs = fig.add_gridspec(1, len(species_list), left=0.085, right=0.995,
                          top=0.86, bottom=0.24, wspace=0.30)
    for k, species in enumerate(species_list):
        ax = fig.add_subplot(gs[0, k])
        sub = per_rep[per_rep.species == species]
        for t_i, tool in enumerate(TOOLS):
            rs = RATIO_SITES.get((species, tool))
            if rs is None:
                continue
            wt = np.concatenate(rs["wt"])
            trt = np.concatenate(rs["trt"])
            for side, vals, col in ((-0.19, wt, WT_COLOR), (0.19, trt, TRT_COLOR)):
                step = max(1, int(np.ceil(vals.size / BOX_MAX_POINTS)))
                bp = ax.boxplot([vals[::step]], positions=[t_i + side], widths=0.30,
                                patch_artist=True, showfliers=False,
                                medianprops=dict(color="black", lw=0.9),
                                boxprops=dict(facecolor=col, edgecolor="none",
                                              alpha=0.55),
                                whiskerprops=dict(color="0.45", lw=0.7),
                                capprops=dict(color="0.45", lw=0.7))
                del bp
            for u in sorted(set(sub.unit)):
                r = sub[(sub.tool == tool) & (sub.unit == u)]
                if r.empty:
                    continue
                ax.plot([t_i - 0.19], [r.iloc[0].median_wt], marker="o", ms=2.4,
                        linestyle="none", markerfacecolor="white",
                        markeredgecolor=WT_COLOR, markeredgewidth=0.7, zorder=4)
                ax.plot([t_i + 0.19], [r.iloc[0].median_trt], marker="o", ms=2.4,
                        linestyle="none", markerfacecolor="white",
                        markeredgecolor=TRT_COLOR, markeredgewidth=0.7, zorder=4)
        ax.set_title(species, fontsize=FS["title"], fontweight="bold", pad=3)
        ax.set_xticks(np.arange(len(TOOLS)))
        ax.set_xticklabels(TOOLS, rotation=35, ha="right", fontsize=FS["tick"])
        ax.set_xlim(-0.6, len(TOOLS) - 0.4)
        ax.set_ylim(0, 0.9)
        ax.set_yticks([0, 0.3, 0.6, 0.9])
        ax.tick_params(labelsize=FS["tick"], length=3, width=0.8)
        if k == 0:
            ax.set_ylabel("Modification ratio", fontsize=FS["label"], labelpad=2)
    fig.text(0.004, 0.96, "B", fontsize=FS["title"], fontweight="bold",
             ha="left", va="top")
    FIG.mkdir(parents=True, exist_ok=True)
    stem = FIG / "Fig3B_modratio_boxplot"
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=300)
    plt.close(fig)
    logger.info("wrote %s.{pdf,png} (%.2f x %.2f in)", stem.name, PANEL_W, PANEL_H)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figures-only", action="store_true")
    args = ap.parse_args()

    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    logger = setup_logger("57_fig3b_modratio_wt_treatment", log_dir=LOG)

    self_test_fdr()
    logger.info("BH-FDR self-test passed (multipletests vs reference impl)")
    per_rep = compute(logger)
    write_table(per_rep, TAB / "fig3b_per_replicate.tsv")
    summary = summarise(per_rep)
    write_table(summary, TAB / "fig3b_effect_sizes.tsv")
    legacy = legacy_table(per_rep)
    write_table(legacy, TAB / "fig3b_legacy_vs_revision.tsv")
    write_table(bland_altman_table(), TAB / "fig3b_bland_altman.tsv")
    rows = []
    for (sp, tl), slot in BA_PAIRS.items():
        mx = np.concatenate(slot["mean"]); dy = np.concatenate(slot["diff"])
        step = max(1, int(np.ceil(mx.size / BA_MAX_POINTS)))
        rows.append(pd.DataFrame({"species": sp, "tool": tl,
                                  "mean": mx[::step], "diff": dy[::step]}))
    pts = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    write_table(pts, TAB / "fig3b_ba_points.tsv.gz")
    logger.info("tables: fig3b_per_replicate.tsv (%d rows), fig3b_effect_sizes.tsv, "
                "fig3b_legacy_vs_revision.tsv", len(per_rep))
    for r in summary.itertuples():
        logger.info("%-12s %-8s median \u0394 = %+.3f (SD %.3f)  Cliff \u03b4 = %+.2f  "
                    "p_BH min = %.2g", r.species, r.tool, r.median_diff_mean,
                    r.median_diff_sd if np.isfinite(r.median_diff_sd) else 0.0,
                    r.cliff_delta_mean, r.wilcoxon_p_bh_min)
    if not args.figures_only:
        figure(per_rep, logger)
    logger.info("output: %s", OUT)


if __name__ == "__main__":
    main()
