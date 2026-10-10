#!/usr/bin/env python3
"""R3-9b -- is "IVT calls > WT calls" in the RNA002 non-m6A analysis a
replicate-structure artefact or a sequencing-depth artefact?

The submitted manuscript compares raw call COUNTS between HeLa WT and the
unmodified IVT libraries (CHEUI_m5C 47,747 WT vs 51,171 IVT, unions of three
replicates).  A count is not a rate: it is the product of how many candidate
sites a library exposes (depth + transcriptome composition) and how often the
tool fires per site.  This package separates those two factors per sequencing
unit (never on a merged call set):

  * sites at risk   = candidate universe (coverage >= --min-cov and reference
                      base can carry the modification), from harmonisation/universe
  * density         = in-universe calls per 10^6 sites at risk
  * coverage strata = calls and sites split by coverage band, then direct
                      standardisation onto the pooled HeLa coverage composition
                      and re-scaling to one common number of sites at risk
  * GLM             = Poisson, calls ~ IVT + offset(log sites at risk),
                      unit-clustered SE -> depth-adjusted IVT/WT rate ratio
  * depth pairing   = WT and IVT units matched by sites-at-risk rank

Outputs (into ``../evidence/``)
------------------------------
units_depth.tsv                 per-unit library size, provenance, replicate class
cov_strata.tsv                  tool x unit x coverage band: sites, calls, rate
units_calls.tsv                 tool x unit: calls, sites at risk, densities,
                                standardised density, count at common depth
wt_vs_ivt_depth_adjusted.tsv    per tool: naive count ratio vs depth-adjusted ratio
nearest_depth_pairs.tsv         WT/IVT pairs matched on sites at risk
replicate_structure_audit.tsv   within-WT vs within-IVT reproducibility

Conventions inherited from R3-9
(``analysis/nonm6a_false_positives/analysis/README.md``): analysis
layer = ``harmonisation/callsets``, denominator = ``harmonisation/universe/RNA002/``
filtered by coverage and modification-capable base, per-unit scoring, Curlcake
``*_partial`` reported but excluded from independent means.

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    analysis/nonm6a_replicate_depth/analysis/rd1_build.py
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
import math
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

HERE = Path(__file__).resolve()
PKG = HERE.parents[1]
PROJECT = HERE.parents[3]
SITES = PROJECT / "data" / "callsets"
SC = (_RB / "data/callsets")
UNI = (_XB / "harmonisation/universe")
META = PROJECT / "$RNAMODBENCH_LOCAL" / "sample_metadata" / "tables"
sys.path.insert(0, str(PROJECT / "src/harmonisation"))
from common.match import fix_chromosome  # noqa: E402

HELA: dict[str, list[str]] = {
    "WT": ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"],
    "IVT": ["HeLa_IVT_rep1", "HeLa_IVT_rep2", "HeLa_IVT_rep3"],
}
CURLCAKE_INDEP = ["Curlcake_IVT_rep1", "Curlcake_IVT_rep3"]
CURLCAKE_SUBSET = ["Curlcake_IVT_rep2_partial"]
CURLCAKE_ALL = CURLCAKE_INDEP + CURLCAKE_SUBSET

TOOLMOD = [
    ("m5C", "CHEUI_m5C"),
    ("Nm", "NanoNm"),
    ("Psi", "NanoMUD_psi"),
    ("Psi", "NanoPsu"),
    ("Psi", "NanoSPA_psU"),
    ("m1Psi", "NanoMUD_m1psi"),
]
MOD_BASES = {"Nm": set("ACGT"), "m5C": {"C", "G"}, "Psi": {"T", "A"},
             "m1Psi": {"A", "T"}}

STRATA = [10, 20, 30, 50, 100, 200, 400, 10_000_000]
LABELS = [f"{STRATA[i]}-{min(STRATA[i + 1] - 1, 99_999)}" for i in range(len(STRATA) - 1)]

SEED = 20260928
BOOT_UNITS = 20_000


# --------------------------------------------------------------------------- #
def log_setup() -> logging.Logger:
    (PKG / "logs").mkdir(exist_ok=True)
    lg = logging.getLogger("rd1")
    lg.setLevel(logging.INFO)
    fh = logging.FileHandler(
        PKG / "logs" / f"rd1_build_{pd.Timestamp.now():%Y%m%d_%H%M%S}.log")
    fh.setFormatter(logging.Formatter("%(asctime)s %(message)s", "%H:%M:%S"))
    lg.addHandler(fh)
    lg.addHandler(logging.StreamHandler(sys.stdout))
    return lg


def species_of(sample: str) -> str:
    return "Curlcake" if sample.startswith("Curlcake") else "Human"


def group_of(sample: str) -> str:
    return "Curlcake_IVT" if sample.startswith("Curlcake") else (
        "HeLa_WT" if "WT" in sample else "HeLa_IVT")


def universe_path(sample: str) -> Path:
    p = (_XB / "harmonisation/universe/RNA002") / species_of(sample) / f"{sample}__universe.tsv"
    return p if p.exists() else p.with_suffix(".tsv.gz")


def callsets_path(sample: str, mod: str, tool: str) -> Path:
    return ((_RB / "data/callsets/RNA002") / species_of(sample) / group_of(sample) / mod
            / tool / f"{sample}.tsv")


# --------------------------------------------------------------------------- #
def load_universe(path: Path, min_cov: int, log: logging.Logger):
    """Return (per-mod {chrom: (pos, cov)}, all-base stats)."""
    df = pd.read_csv(path, sep="\t", usecols=["chrom", "pos", "base", "coverage"],
                     dtype={"chrom": "category", "pos": np.int64,
                            "base": "category", "coverage": np.int32})
    cov = df["coverage"].to_numpy(np.float64)
    stats = {
        "n_rows": int(len(df)),
        "n_cov5": int((cov >= 5).sum()), "n_cov10": int((cov >= 10).sum()),
        "n_cov20": int((cov >= 20).sum()), "n_cov50": int((cov >= 50).sum()),
        "n_cov100": int((cov >= 100).sum()),
        "sum_cov": float(cov.sum()), "sum_cov_ge10": float(cov[cov >= 10].sum()),
        "median_cov_ge10": float(np.median(cov[cov >= 10])) if (cov >= 10).any() else np.nan,
    }
    codes = df["chrom"].cat.codes.to_numpy()
    pos_all = df["pos"].to_numpy(np.int64)
    cov_all = df["coverage"].to_numpy(np.int64)
    base = df["base"].astype(str).to_numpy()
    cats = [str(c) for c in df["chrom"].cat.categories]
    out: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]] = {}
    for mod, bases in MOD_BASES.items():
        mask = np.isin(base, list(bases)) & (cov_all >= min_cov)
        c, p, v = codes[mask], pos_all[mask], cov_all[mask]
        order = np.lexsort((p, c))
        c, p, v = c[order], p[order], v[order]
        edges = np.flatnonzero(np.r_[True, c[1:] != c[:-1]])
        d: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for i, st in enumerate(edges):
            en = edges[i + 1] if i + 1 < len(edges) else len(c)
            name = fix_chromosome(cats[c[st]])
            pp, vv = p[st:en], v[st:en]
            if name in d:                                   # merged chrom names
                op, ov = d[name]
                j = np.lexsort((np.r_[pp, ov], np.r_[op, pp]))
                d[name] = (np.r_[op, pp][j], np.r_[ov, vv][j])
            else:
                d[name] = (pp, vv)
        out[mod] = d
        log.info(f"  {path.name} {mod}: sites at risk (cov>={min_cov}) = "
                 f"{sum(len(a) for a, _ in d.values()):,}")
    return out, stats


def load_calls(sample: str, mod: str, tool: str):
    path = callsets_path(sample, mod, tool)
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t", usecols=["chrom", "start"],
                     dtype={"chrom": str, "start": np.int64})
    df["chrom"] = df["chrom"].map(fix_chromosome)
    if df.empty:
        return {c: np.zeros(0, np.int64) for c in ["x"]}
    return {c: np.unique(g["start"].to_numpy(np.int64))
            for c, g in df.groupby("chrom", sort=False)}


def cov_of_calls(calls: dict[str, np.ndarray], uni: dict[str, tuple[np.ndarray, np.ndarray]]):
    """Coverage of the calls that fall inside the universe + counts per stratum."""
    picked: list[np.ndarray] = []
    n_out = 0
    for chrom, pos in calls.items():
        pair = uni.get(chrom)
        if pair is None or pair[0].size == 0:
            n_out += int(pos.size)
            continue
        upos, ucov = pair
        j = np.clip(np.searchsorted(upos, pos), 0, max(upos.size - 1, 0))
        ok = upos[j] == pos
        n_out += int((~ok).sum())
        if ok.any():
            picked.append(ucov[j[ok]])
    cov = np.concatenate(picked) if picked else np.zeros(0, np.int64)
    return cov, n_out


def hist_sites(uni: dict[str, tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
    h = np.zeros(len(LABELS), dtype=np.int64)
    for _, (_, ucov) in uni.items():
        idx = np.clip(np.searchsorted(STRATA, ucov, side="right") - 1, 0, len(LABELS) - 1)
        h += np.bincount(idx, minlength=len(LABELS))
    return h


def jaccard(a: dict[str, np.ndarray], b: dict[str, np.ndarray]) -> float:
    inter = union = 0
    for chrom in set(a) | set(b):
        pa, pb = a.get(chrom), b.get(chrom)
        if pa is None or pb is None:
            union += int(pa.size) if pa is not None else int(pb.size)
            continue
        union += int(np.union1d(pa, pb).size)
        inter += int(np.intersect1d(pa, pb).size)
    return inter / union if union else float("nan")


def global_jaccard(sets: list[dict[str, np.ndarray]]) -> float:
    inter = union = 0
    keys = set().union(*[set(s) for s in sets]) if sets else set()
    for chrom in keys:
        parts = [s[chrom] for s in sets if chrom in s]
        if not parts:
            continue
        union += int(np.unique(np.concatenate(parts)).size)
        if len(parts) == len(sets):
            common = parts[0]
            for p in parts[1:]:
                common = np.intersect1d(common, p)
            inter += int(common.size)
    return inter / union if union else float("nan")


# --------------------------------------------------------------------------- #
def build_units_table(log: logging.Logger) -> pd.DataFrame:
    master = pd.read_csv((_XB / "sample_metadata/tables/samples_master.csv"))
    keep = HELA["WT"] + HELA["IVT"] + CURLCAKE_ALL
    cols = ["canonical_name", "dataset_group", "condition", "condition_class",
            "sample_type", "replicate_structure", "chemistry", "study_accession",
            "sra_sample_accession", "sra_run_accessions", "sra_read_count_sum",
            "local_read_count", "source_publication"]
    df = master[master["canonical_name"].isin(keep)][cols].copy()
    df = df.rename(columns={"canonical_name": "unit"})
    df["dataset"] = np.where(df["unit"].str.startswith("HeLa"), "HeLa", "Curlcake")

    def verdict(u: str) -> tuple[str, str]:
        if u.endswith("_partial"):
            return ("subset", "depth-matched subset of another unit -> NOT independent")
        if u.startswith("HeLa_IVT"):
            return ("unverified",
                    "runs filed under ONE BioSample (SAMN22863220); experiment "
                    "titles say rep.1/2/3 but no separate library record -> "
                    "independence cannot be verified from the accession")
        if u.startswith("HeLa_WT"):
            return ("independent",
                    "distinct BioSample per unit (HeLa_M11A/M11B/M11C WT) -> "
                    "independent libraries")
        return ("independent", "distinct BioSample / construct -> independent")

    df[["replicate_class", "independence_note"]] = df["unit"].apply(
        lambda u: pd.Series(verdict(u)))
    df["biosample_shared_by_n_units"] = df.groupby("sra_sample_accession")[
        "unit"].transform("count")
    log.info(f"units table: {len(df)} rows")
    return df


def analyse(min_cov: int, log: logging.Logger):
    samples = HELA["WT"] + HELA["IVT"] + CURLCAKE_ALL
    uni, stats_cache, site_hist = {}, {}, {}
    for s in samples:
        p = universe_path(s)
        if not p.exists():
            log.warning(f"no universe for {s}")
            continue
        log.info(f"universe: {s}")
        uni[s], stats_cache[s] = load_universe(p, min_cov, log)
        for m in MOD_BASES:
            site_hist[(s, m)] = hist_sites(uni[s][m])

    unit_rows, strata_rows = [], []
    calls_cache: dict[tuple[str, str, str], dict[str, np.ndarray]] = {}
    for mod, tool in TOOLMOD:
        for s in samples:
            if s not in uni:
                continue
            calls = load_calls(s, mod, tool)
            if calls is None:
                log.info(f"  no callset {tool} / {s}")
                continue
            calls_cache[(mod, tool, s)] = calls
            n_at_risk = int(site_hist[(s, mod)].sum())
            if n_at_risk == 0:
                continue
            n_raw = int(sum(v.size for v in calls.values()))
            cov, n_out = cov_of_calls(calls, uni[s][mod])
            h_sites = site_hist[(s, mod)]
            h_calls = np.bincount(
                np.clip(np.searchsorted(STRATA, cov, side="right") - 1,
                        0, len(LABELS) - 1), minlength=len(LABELS)).astype(np.int64)
            for i, lab in enumerate(LABELS):
                strata_rows.append({
                    "mod_type": mod, "tool": tool, "unit": s,
                    "dataset": "HeLa" if s.startswith("HeLa") else "Curlcake",
                    "condition": "WT" if s in HELA["WT"] else "IVT",
                    "stratum": lab, "cov_lo": STRATA[i],
                    "sites_at_risk": int(h_sites[i]), "calls": int(h_calls[i]),
                    "calls_per_1e6_sites": (h_calls[i] / h_sites[i] * 1e6)
                    if h_sites[i] else np.nan,
                })
            ge20_s, ge20_c = int(h_sites[1:].sum()), int(h_calls[1:].sum())
            unit_rows.append({
                "mod_type": mod, "tool": tool, "unit": s,
                "dataset": "HeLa" if s.startswith("HeLa") else "Curlcake",
                "condition": "WT" if s in HELA["WT"] else "IVT",
                "n_calls_raw": n_raw, "n_calls_in_universe": int(cov.size),
                "n_calls_out_of_universe": n_out, "sites_at_risk": n_at_risk,
                "mean_cov_of_calls": float(cov.mean()) if cov.size else np.nan,
                "median_cov_calls": float(np.median(cov)) if cov.size else np.nan,
                "density_per_1e6": cov.size / n_at_risk * 1e6,
                "sites_at_risk_ge20": ge20_s, "calls_ge20": ge20_c,
                "density_ge20_per_1e6": (ge20_c / ge20_s * 1e6) if ge20_s else np.nan,
                "universe_rows": stats_cache[s]["n_rows"],
                "n_cov10_all_bases": stats_cache[s]["n_cov10"],
                "mean_cov_all_ge10": round(stats_cache[s]["sum_cov_ge10"]
                                           / max(stats_cache[s]["n_cov10"], 1), 3),
                "median_cov_all_ge10": stats_cache[s]["median_cov_ge10"],
                "aligned_bases_proxy": int(stats_cache[s]["sum_cov"]),
            })
        log.info(f"done {mod}/{tool}")

    units = pd.DataFrame(unit_rows)
    strata = pd.DataFrame(strata_rows)

    # ---- direct standardisation + common-depth counts ------------------------
    extra = []
    for (mod, tool), g in strata.groupby(["mod_type", "tool"]):
        ps = g.pivot_table(index="unit", columns="stratum", values="sites_at_risk",
                           aggfunc="sum", fill_value=0).reindex(columns=LABELS, fill_value=0)
        pc = g.pivot_table(index="unit", columns="stratum", values="calls",
                           aggfunc="sum", fill_value=0).reindex(columns=LABELS, fill_value=0)
        hela_u = [u for u in ps.index if u.startswith("HeLa")] or list(ps.index)
        w = ps.loc[hela_u].to_numpy().sum(axis=0).astype(float)
        w = w / w.sum()
        rate = np.full(ps.shape, np.nan)
        with np.errstate(divide="ignore", invalid="ignore"):
            rate = np.where(ps.to_numpy() > 0,
                            pc.to_numpy() / ps.to_numpy() * 1e6, np.nan)
        std = np.nansum(np.where(np.isnan(rate), 0.0, rate) * w, axis=1)
        n_common = int(ps.loc[hela_u].to_numpy().sum(axis=1).min())
        for k, u in enumerate(ps.index):
            obs_sites = int(ps.loc[u].to_numpy().sum())
            extra.append({
                "unit": u, "mod_type": mod, "tool": tool,
                "standardised_density_per_1e6": std[k],
                "count_at_common_depth": std[k] * n_common / 1e6,
                "common_depth_n": n_common,
                "standardisation_weights": "|".join(f"{x:.4f}" for x in w),
                "density_if_shallowest_library": (
                    std[k] * n_common / 1e6),
                "observed_density_check": (pc.loc[u].to_numpy().sum()
                                           / obs_sites * 1e6) if obs_sites else np.nan,
            })
    units = units.merge(pd.DataFrame(extra), on=["unit", "mod_type", "tool"],
                        how="left")
    return units, strata, calls_cache


def nearest_depth_pairs(units: pd.DataFrame, log: logging.Logger) -> pd.DataFrame:
    rows = []
    g = units[units["dataset"] == "HeLa"]
    for (mod, tool), sub in g.groupby(["mod_type", "tool"]):
        wt = sub[sub["condition"] == "WT"].sort_values("sites_at_risk").reset_index(drop=True)
        iv = sub[sub["condition"] == "IVT"].sort_values("sites_at_risk").reset_index(drop=True)
        for i in range(min(len(wt), len(iv))):
            a, b = wt.loc[i], iv.loc[i]
            rows.append({
                "mod_type": mod, "tool": tool, "pair_rank_by_depth": i + 1,
                "wt_unit": a["unit"], "ivt_unit": b["unit"],
                "wt_sites_at_risk": int(a["sites_at_risk"]),
                "ivt_sites_at_risk": int(b["sites_at_risk"]),
                "depth_ratio_ivt_over_wt": b["sites_at_risk"] / a["sites_at_risk"],
                "wt_calls": int(a["n_calls_in_universe"]),
                "ivt_calls": int(b["n_calls_in_universe"]),
                "count_ratio_ivt_over_wt": b["n_calls_in_universe"] / a["n_calls_in_universe"]
                if a["n_calls_in_universe"] else np.nan,
                "wt_density": a["density_per_1e6"], "ivt_density": b["density_per_1e6"],
                "density_ratio_ivt_over_wt": b["density_per_1e6"] / a["density_per_1e6"]
                if a["density_per_1e6"] else np.nan,
                "std_density_ratio": (b["standardised_density_per_1e6"]
                                      / a["standardised_density_per_1e6"])
                if a["standardised_density_per_1e6"] else np.nan,
            })
    return pd.DataFrame(rows)


def compare_wt_ivt(units: pd.DataFrame, log: logging.Logger) -> pd.DataFrame:
    rows = []
    for (mod, tool), g in units[units["dataset"] == "HeLa"].groupby(["mod_type", "tool"]):
        wt, iv = g[g["condition"] == "WT"], g[g["condition"] == "IVT"]
        if len(wt) < 2 or len(iv) < 2:
            continue
        base = {"mod_type": mod, "tool": tool}

        def blk(label, cc, cd):
            w_c, i_c = wt[cc].sum(), iv[cc].sum()
            w_d, i_d = wt[cd].sum(), iv[cd].sum()
            dd_w = w_c / w_d * 1e6 if w_d else np.nan
            dd_i = i_c / i_d * 1e6 if i_d else np.nan
            return {**base, "metric": label,
                    "wt_count": int(w_c), "ivt_count": int(i_c),
                    "wt_sites_at_risk": int(w_d), "ivt_sites_at_risk": int(i_d),
                    "sites_at_risk_ratio_ivt_over_wt": i_d / w_d if w_d else np.nan,
                    "count_ratio_ivt_over_wt": i_c / w_c if w_c else np.nan,
                    "wt_density": dd_w, "ivt_density": dd_i,
                    "density_ratio_ivt_over_wt": dd_i / dd_w if dd_w else np.nan}

        rows.append(blk("raw calls / universe rows", "n_calls_raw", "sites_at_risk"))
        rows.append(blk("in-universe calls / sites at risk",
                        "n_calls_in_universe", "sites_at_risk"))
        rows.append(blk("coverage>=20 band only", "calls_ge20", "sites_at_risk_ge20"))

        # standardised density + common-depth count
        w_s, i_s = wt["standardised_density_per_1e6"].mean(), iv["standardised_density_per_1e6"].mean()
        w_c, i_c = wt["count_at_common_depth"].mean(), iv["count_at_common_depth"].mean()
        rng = np.random.default_rng(SEED)
        a = wt["standardised_density_per_1e6"].to_numpy(float)
        b = iv["standardised_density_per_1e6"].to_numpy(float)
        draws = []
        for _ in range(BOOT_UNITS):
            x = a[rng.integers(0, len(a), len(a))].mean()
            y = b[rng.integers(0, len(b), len(b))].mean()
            if x > 0:
                draws.append(y / x)
        d = np.asarray(draws)
        rows.append({**base,
                     "metric": f"coverage-standardised, n={int(wt['common_depth_n'].iloc[0])} "
                               "sites at risk (common depth)",
                     "wt_count": round(w_c, 1), "ivt_count": round(i_c, 1),
                     "count_ratio_ivt_over_wt": i_c / w_c if w_c else np.nan,
                     "wt_density": w_s, "ivt_density": i_s,
                     "density_ratio_ivt_over_wt": i_s / w_s if w_s else np.nan,
                     "boot_ratio_mean": float(d.mean()),
                     "boot_ci_lo": float(np.percentile(d, 2.5)),
                     "boot_ci_hi": float(np.percentile(d, 97.5))})

        # GLM
        sub = g.copy()
        X = np.column_stack([np.ones(len(sub)),
                             (sub["condition"] == "IVT").astype(float)])
        try:
            res = sm.GLM(sub["n_calls_in_universe"].to_numpy(float), X,
                         family=sm.families.Poisson(),
                         offset=np.log(sub["sites_at_risk"].to_numpy(float))
                         ).fit(cov_type="cluster",
                               cov_kwds={"groups": sub["unit"].to_numpy()},
                               scale="X2")
            beta, se = res.params[1], res.bse[1]
            rows.append({**base,
                         "metric": "GLM Poisson + log(sites) offset, unit-clustered",
                         "glm_rate_ratio": math.exp(beta),
                         "glm_ci_lo": math.exp(beta - 1.96 * se),
                         "glm_ci_hi": math.exp(beta + 1.96 * se),
                         "glm_p": float(res.pvalues[1]),
                         "count_ratio_ivt_over_wt": math.exp(beta)})
        except Exception as exc:
            log.warning(f"GLM failed {tool}: {exc}")

        # scaling of counts with depth across the 6 units
        xs = np.log(g["sites_at_risk"].to_numpy(float))
        ys = np.log(g["n_calls_in_universe"].to_numpy(float) + 0.5)
        slope, icpt = np.polyfit(xs, ys, 1)
        resid = ys - (slope * xs + icpt)
        rows.append({**base,
                     "metric": "scaling log(calls) ~ log(sites at risk), 6 units",
                     "slope_calls_vs_depth": float(slope),
                     "r2_calls_vs_depth": float(np.corrcoef(xs, ys)[0, 1] ** 2),
                     "resid_mean_WT": float(resid[(g["condition"] == "WT").to_numpy()].mean()),
                     "resid_mean_IVT": float(resid[(g["condition"] == "IVT").to_numpy()].mean())})

    cols = ["mod_type", "tool", "metric", "wt_count", "ivt_count",
            "wt_sites_at_risk", "ivt_sites_at_risk",
            "sites_at_risk_ratio_ivt_over_wt", "count_ratio_ivt_over_wt",
            "wt_density", "ivt_density", "density_ratio_ivt_over_wt",
            "boot_ratio_mean", "boot_ci_lo", "boot_ci_hi", "glm_rate_ratio",
            "glm_ci_lo", "glm_ci_hi", "glm_p", "slope_calls_vs_depth",
            "r2_calls_vs_depth", "resid_mean_WT", "resid_mean_IVT"]
    df = pd.DataFrame(rows)
    for c in cols:
        if c not in df.columns:
            df[c] = np.nan
    return df[cols]


def replicate_structure(units, calls_cache, log) -> pd.DataFrame:
    rows = []
    for (mod, tool), g in units[units["dataset"] == "HeLa"].groupby(["mod_type", "tool"]):
        out = {"mod_type": mod, "tool": tool}
        for cond, lst in HELA.items():
            sets = [calls_cache[(mod, tool, s)] for s in lst
                    if (mod, tool, s) in calls_cache]
            dens = g[g["condition"] == cond]["density_per_1e6"].to_numpy(float)
            pw = [jaccard(x, y) for x, y in combinations(sets, 2)]
            out[f"n_units_{cond}"] = len(sets)
            out[f"global_jaccard_{cond}"] = global_jaccard(sets) if len(sets) >= 2 else np.nan
            out[f"mean_pairwise_jaccard_{cond}"] = float(np.mean(pw)) if pw else np.nan
            out[f"density_mean_{cond}"] = float(dens.mean()) if dens.size else np.nan
            out[f"density_cv_{cond}"] = (float(dens.std(ddof=1) / dens.mean())
                                         if len(dens) > 1 and dens.mean() else np.nan)
            out[f"callcount_cv_{cond}"] = out[f"density_cv_{cond}"]
        out["jaccard_IVT_minus_WT"] = (out.get("global_jaccard_IVT", np.nan)
                                       - out.get("global_jaccard_WT", np.nan))
        rows.append(out)
    return pd.DataFrame(rows)


def curlcake_summary(units: pd.DataFrame, log) -> pd.DataFrame:
    g = units[units["dataset"] == "Curlcake"].copy()
    g["independent"] = ~g["unit"].str.endswith("_partial")
    rows = []
    for (mod, tool), sub in g.groupby(["mod_type", "tool"]):
        ind = sub[sub["independent"]]
        rows.append({"mod_type": mod, "tool": tool,
                     "n_independent_units": len(ind),
                     "calls_by_unit": "|".join(f"{r.unit}:{int(r.n_calls_in_universe)}"
                                               for r in sub.itertuples()),
                     "sites_at_risk_by_unit": "|".join(f"{r.unit}:{int(r.sites_at_risk)}"
                                                       for r in sub.itertuples()),
                     "density_by_unit_per_1e6": "|".join(
                         f"{r.unit}:{r.density_per_1e6:.1f}" for r in sub.itertuples()),
                     "mean_density_independent": float(ind["density_per_1e6"].mean())
                     if len(ind) else np.nan,
                     "spread_independent_units": (
                         float(ind["density_per_1e6"].max() - ind["density_per_1e6"].min())
                         if len(ind) > 1 else np.nan)})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-cov", type=int, default=10)
    ap.add_argument("--out", type=Path, default=PKG / "evidence")
    args = ap.parse_args()
    log = log_setup()
    log.info(f"R3-9b start min_cov={args.min_cov} seed={SEED}")
    args.out.mkdir(parents=True, exist_ok=True)

    units, strata, calls_cache = analyse(args.min_cov, log)
    units.to_csv(args.out / "units_calls.tsv", sep="\t", index=False)
    strata.to_csv(args.out / "cov_strata.tsv", sep="\t", index=False)

    depth = build_units_table(log)
    agg = (units.groupby("unit")
           .agg(sites_at_risk_largest_mod=("sites_at_risk", "max"),
                sites_at_risk_Nm=("sites_at_risk", lambda x: np.nan),
                mean_cov_all_ge10=("mean_cov_all_ge10", "first"),
                median_cov_all_ge10=("median_cov_all_ge10", "first"),
                aligned_bases_proxy=("aligned_bases_proxy", "first"))
           .reset_index())
    nm = units[units["mod_type"] == "Nm"].set_index("unit")["sites_at_risk"]
    agg["sites_at_risk_Nm"] = agg["unit"].map(nm)
    depth = depth.merge(agg.drop(columns=["sites_at_risk_Nm"]), on="unit", how="left")
    depth.to_csv(args.out / "units_depth.tsv", sep="\t", index=False)

    compare_wt_ivt(units, log).to_csv(args.out / "wt_vs_ivt_depth_adjusted.tsv",
                                      sep="\t", index=False)
    nearest_depth_pairs(units, log).to_csv(args.out / "nearest_depth_pairs.tsv",
                                           sep="\t", index=False)
    replicate_structure(units, calls_cache, log).to_csv(
        args.out / "replicate_structure_audit.tsv", sep="\t", index=False)
    curlcake_summary(units, log).to_csv(args.out / "curlcake_units_summary.tsv",
                                        sep="\t", index=False)

    with pd.option_context("display.width", 240, "display.max_columns", 30):
        cmp_df = pd.read_csv(args.out / "wt_vs_ivt_depth_adjusted.tsv", sep="\t")
        log.info("\n== naive vs depth-adjusted (HeLa) ==\n" + cmp_df[
            cmp_df["metric"].str.startswith(("in-universe", "coverage-standardised",
                                             "GLM", "coverage>=20"))]
            [["mod_type", "tool", "metric", "count_ratio_ivt_over_wt",
              "density_ratio_ivt_over_wt", "sites_at_risk_ratio_ivt_over_wt"]]
            .to_string(index=False))
        rep = pd.read_csv(args.out / "replicate_structure_audit.tsv", sep="\t")
        log.info("\n== within-group reproducibility ==\n" + rep[
            ["mod_type", "tool", "global_jaccard_WT", "global_jaccard_IVT",
             "density_cv_WT", "density_cv_IVT"]].to_string(index=False))
    log.info("tables written")


if __name__ == "__main__":
    main()
