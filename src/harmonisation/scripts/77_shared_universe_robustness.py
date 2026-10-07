#!/usr/bin/env python3
"""77 -- the six non-m6A tools' IVT/WT ratios on a shared testable universe.

Why this table exists
---------------------
The user asked (2026-10-05) whether a sampling rule could bring the NanoNm IVT/WT ratio of
Table S5 down.  The one sampling change this study already uses elsewhere -- Table S3 compares
its species pairs "on their shared testable universe" -- is priced here for all six non-m6A
tools alike: both conditions are restricted to the sites that *every* unit of both conditions
covers at the study's own floor (coverage >= 10 with a reference base that can carry that
modification).  That removes the sites a shallower library cannot call at all, so the two arms
sit on one and the same site space; it is not a truncation aimed at one tool.

Reported per tool and per basis:

* ``own_candidate_sets`` -- the delivered basis: each unit's own candidate sites, from the
  depth diagnosis (`cov_strata.tsv`), which is what Table S5 prints;
* ``shared_universe``    -- the same weighted estimator and interval on the intersection of
  the six units' candidate sites, with weights taken from that intersection.

The estimator is `rd2_matched_bands.py`'s: per-stratum calls per million candidate sites,
strata combined with the pooled candidate-site weights of the six HeLa units, plus a unit-level
percentile bootstrap (B = 20000, seed 20260928) for the interval.  The per-unit site counts are
asserted against the delivered counts table and against `r39.Resources`, so this script cannot
drift away from the study's own site spaces.

What it writes
--------------
``tables/s7_shared_universe.tsv``  12 rows: 6 tools x 2 bases
``logs/77_anchors.tsv``            the shared-universe sizes and every ratio, for the gate
``cache/<unit>__<mod>.npz``        per-unit candidate sites and coverage, one file per
                                   modification, so the universe files are read once each

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    src/harmonisation/scripts/77_shared_universe_robustness.py
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

OUT = (_RB / "figures/figureS7")
TAB, LOG = (_RB / "figures/figureS7/tables"), (_RB / "figures/figureS7/logs")
CACHE = (_XB / "submission/_build/shared_universe_20261005/cache")
COUNTS = TAB / "s7_counts_per_replicate.tsv"
RD2 = ((_RB / "analysis/R3-9b_nonm6a_replicate_depth_20260928/evidence/band_matched_rate_ratios.tsv"))
COV_STRATA = RD2.with_name("cov_strata.tsv")
R39 = (_RB / "analysis/nonm6a_false_positives")
_spec = importlib.util.spec_from_file_location(
    "r39_build_evidence", (_RB / "analysis/nonm6a_false_positives/analysis/r39_build_evidence.py"))
r39 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(r39)

#: the six non-m6A tools and the modification each of them calls.  m1Psi shares Psi's
#: reference bases in the reader (`r39.MOD_BASES` has no key of its own for it).
MOD_OF_TOOL = {
    "CHEUI_m5C": "m5C",
    "NanoMUD_psi": "Psi",
    "NanoMUD_m1psi": "m1Psi",
    "NanoNm": "Nm",
    "NanoSPA_psU": "Psi",
    "NanoPsu": "Psi",
}
MODS = ("Nm", "Psi", "m5C", "m1Psi")
TOOLS = tuple(MOD_OF_TOOL)
WT_UNITS = ("HeLa_WT1", "HeLa_WT2", "HeLa_WT3")
IVT_UNITS = ("HeLa_IVT_rep1", "HeLa_IVT_rep2", "HeLa_IVT_rep3")
UNITS = WT_UNITS + IVT_UNITS
MIN_COVERAGE = 10
STRATA = ("10-19", "20-29", "30-49", "50-99", "100-199", "200-399", "400-99999")
STRATA_LO = (10, 20, 30, 50, 100, 200, 400)
SCALE = 10 ** 10
SEED, BOOT = 20260928, 20_000


def bases_of(mod: str) -> set[str]:
    return r39.MOD_BASES.get(mod, r39.MOD_BASES["Psi"])


def build_cache() -> None:
    """One read of each unit's universe file, then one file per modification."""
    CACHE.mkdir(parents=True, exist_ok=True)
    for u in UNITS:
        t0 = time.time()
        df = pd.read_csv(r39.universe_path(u), sep="\t",
                         usecols=["chrom", "pos", "base", "coverage"],
                         dtype={"chrom": str, "base": str})
        df = df[df["coverage"] >= MIN_COVERAGE]
        df["chrom"] = [r39.fix_chromosome(c) for c in df["chrom"]]
        chroms = sorted(df.chrom.unique())
        idx = {c: i for i, c in enumerate(chroms)}
        packed = df.chrom.map(idx).to_numpy("int64") * SCALE + df.pos.to_numpy("int64")
        for mod in MODS:
            keep = df.base.isin(bases_of(mod)).to_numpy()
            key = np.sort(packed[keep])
            cov = df.coverage.to_numpy("int64")[keep]
            order = np.argsort(packed[keep])
            np.savez(CACHE / f"{u}__{mod}.npz", key=key, cov=cov[order],
                     chroms=np.array(chroms, dtype=object))
        print(f"  [cache] {u:16s} {len(df)/1e6:5.2f} M rows -> {MODS}  "
              f"[{time.time()-t0:.0f}s]")


def load_mod(unit: str, mod: str) -> tuple[np.ndarray, np.ndarray, list[str]]:
    z = np.load(CACHE / f"{unit}__{mod}.npz", allow_pickle=True)
    return z["key"], z["cov"], list(z["chroms"])


def bands(cov: np.ndarray) -> np.ndarray:
    return np.asarray(pd.cut(cov, bins=list(STRATA_LO) + [10 ** 9], labels=list(STRATA),
                             right=False).astype(str))


def band_counts(band: np.ndarray) -> np.ndarray:
    return (pd.Series(band).value_counts().reindex(list(STRATA)).fillna(0).to_numpy())


def boot_ratio(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    """The depth diagnosis' interval: resample the units of each arm, B = 20000."""
    rng = np.random.default_rng(SEED)
    d = []
    for _ in range(BOOT):
        x = a[rng.integers(0, len(a), len(a))].mean()
        y = b[rng.integers(0, len(b), len(b))].mean()
        if x > 0:
            d.append(y / x)
    d = np.asarray(d)
    return float(np.mean(d)), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def main() -> None:
    t0 = time.time()
    TAB.mkdir(parents=True, exist_ok=True)
    LOG.mkdir(parents=True, exist_ok=True)
    if not all((CACHE / f"{u}__{m}.npz").exists() for u in UNITS for m in MODS):
        build_cache()
    res = r39.Resources(min_cov=MIN_COVERAGE)
    delivered = pd.read_csv(COUNTS, sep="\t")
    strata = pd.read_csv(COV_STRATA, sep="\t")
    rows, anchors = [], []
    for mod in MODS:
        tools = [t for t in TOOLS if MOD_OF_TOOL[t] == mod]
        keys_of, cov_of, chroms_of = {}, {}, {}
        for u in UNITS:
            key, cov, chroms = load_mod(u, mod)
            want = r39.n_positions(res.universe(u, mod))
            if key.size != want:
                raise SystemExit(f"{u}/{mod}: cache holds {key.size} sites, "
                                 f"r39.Resources holds {want}; rebuild the cache")
            keys_of[u], cov_of[u], chroms_of[u] = key, cov, chroms
        shared = keys_of[UNITS[0]]
        for u in UNITS[1:]:
            shared = np.intersect1d(shared, keys_of[u], assume_unique=True)
        anchors.append((f"shared_universe_{mod}", int(shared.size)))
        print(f"[{mod}] shared testable universe = {shared.size:,} sites "
              f"[{time.time()-t0:.0f}s]")
        for u in UNITS:
            mask = np.isin(keys_of[u], shared)
            keys_of[u], cov_of[u] = keys_of[u][mask], cov_of[u][mask]
            want = int(delivered[(delivered.tool == tools[0])
                                 & (delivered["sample"] == u)].universe_n.iloc[0])
            if keys_of[u].size > want:
                raise SystemExit(f"{u}/{mod}: the shared universe ({keys_of[u].size:,}) is "
                                 f"larger than the delivered candidate set ({want:,})")
        at_risk = {u: band_counts(bands(cov_of[u])) for u in UNITS}
        w = sum(at_risk[u] for u in UNITS).astype(float)
        w = w / w.sum()
        for tool in tools:
            call_keys, call_band = {}, {}
            for u in UNITS:
                calls = res.calls(u, mod, tool)
                idx = {c: i for i, c in enumerate(chroms_of[u])}
                ck = np.array(sorted({idx[c] * SCALE + p for c, pos in calls.items()
                                      if c in idx for p in pos.tolist()}), dtype="int64")
                pos = np.clip(np.searchsorted(keys_of[u], ck), 0, keys_of[u].size - 1)
                hit = keys_of[u][pos] == ck
                call_keys[u] = ck[hit]
                call_band[u] = bands(cov_of[u][pos[hit]])
            rate = {u: np.where(at_risk[u] > 0,
                                band_counts(call_band[u]) / at_risk[u] * 1e6, np.nan)
                    for u in UNITS}
            r_wt = np.nanmean(np.vstack([rate[u] for u in WT_UNITS]), axis=0)
            r_ivt = np.nanmean(np.vstack([rate[u] for u in IVT_UNITS]), axis=0)
            value = float(np.nansum(w * r_ivt) / np.nansum(w * r_wt))
            per_wt = np.where(np.isnan(np.vstack([rate[u] for u in WT_UNITS])), 0.0,
                              np.vstack([rate[u] for u in WT_UNITS])) @ w
            per_ivt = np.where(np.isnan(np.vstack([rate[u] for u in IVT_UNITS])), 0.0,
                               np.vstack([rate[u] for u in IVT_UNITS])) @ w
            _mean, lo, hi = boot_ratio(per_wt, per_ivt)
            wt_union = len(set().union(*[set(call_keys[u].tolist()) for u in WT_UNITS]))
            ivt_union = len(set().union(*[set(call_keys[u].tolist()) for u in IVT_UNITS]))
            rows.append(dict(
                tool=tool, mod_type=mod, basis="shared_universe",
                matched_ratio=round(value, 4), boot_ci_lo=round(lo, 4),
                boot_ci_hi=round(hi, 4), wt_calls=wt_union, ivt_calls=ivt_union,
                raw_ratio=round(ivt_union / wt_union, 4) if wt_union else np.nan,
                sites=int(shared.size),
            ))
            anchors += [(f"{tool}_{mod}_shared_match", round(value, 4)),
                        (f"{tool}_{mod}_shared_lo", round(lo, 4)),
                        (f"{tool}_{mod}_shared_hi", round(hi, 4)),
                        (f"{tool}_{mod}_shared_wt_calls", wt_union),
                        (f"{tool}_{mod}_shared_ivt_calls", ivt_union)]
            print(f"    {tool:14s} matched {value:.4f} ({lo:.3f}-{hi:.3f})  "
                  f"calls WT {wt_union} / IVT {ivt_union}")

    #: the delivered basis, recomputed from the depth diagnosis, so both columns of the table
    #: come out of this one script: same estimator, same weights, each unit's own candidate set
    rd2 = pd.read_csv(RD2, sep="\t")
    for tool in TOOLS:
        mod = MOD_OF_TOOL[tool]
        t = strata[(strata.mod_type == mod) & (strata.tool == tool)
                   & (strata.dataset == "HeLa")]
        ps = t.pivot_table(index="unit", columns="stratum", values="sites_at_risk",
                           aggfunc="sum").reindex(columns=list(STRATA), fill_value=0)
        pc = t.pivot_table(index="unit", columns="stratum", values="calls",
                           aggfunc="sum").reindex(columns=list(STRATA), fill_value=0)
        w = ps.to_numpy().sum(axis=0).astype(float)
        w = w / w.sum()
        rate = np.where(ps.to_numpy() > 0, pc.to_numpy() / ps.to_numpy() * 1e6, np.nan)
        idx = {u: list(ps.index).index(u) for u in ps.index}
        value = float(np.nansum(w * np.nanmean(rate[[idx[u] for u in IVT_UNITS]], axis=0))
                      / np.nansum(w * np.nanmean(rate[[idx[u] for u in WT_UNITS]], axis=0)))
        published = float(rd2.query("tool == @tool").band_matched_rate_ratio.iloc[0])
        if not np.isclose(value, published, rtol=2e-3):
            raise SystemExit(f"{tool}: own-candidate-set ratio {value:.4f} does not reproduce "
                             f"the depth diagnosis ({published:.4f})")
        rows.append(dict(
            tool=tool, mod_type=mod, basis="own_candidate_sets", matched_ratio=round(value, 4),
            boot_ci_lo=float(rd2.query("tool == @tool").boot_matched_ci_lo.iloc[0]),
            boot_ci_hi=float(rd2.query("tool == @tool").boot_matched_ci_hi.iloc[0]),
            wt_calls=int(pc.loc[list(WT_UNITS)].to_numpy().sum()),
            ivt_calls=int(pc.loc[list(IVT_UNITS)].to_numpy().sum()),
            raw_ratio=np.nan, sites=int(ps.to_numpy().sum()),
        ))
        anchors.append((f"{tool}_{mod}_own_match", round(value, 4)))
    tab = pd.DataFrame(rows).sort_values(["basis", "matched_ratio"], ignore_index=True)
    tab.to_csv(TAB / "s7_shared_universe.tsv", sep="\t", index=False)
    pd.DataFrame(anchors, columns=["key", "value"]).to_csv(LOG / "77_anchors.tsv", sep="\t",
                                                           index=False)
    print(tab.to_string(index=False))
    print(f"-> {TAB / 's7_shared_universe.tsv'}")
    print(f"-> {LOG / '77_anchors.tsv'}  [{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    main()
