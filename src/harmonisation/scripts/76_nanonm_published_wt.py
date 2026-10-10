#!/usr/bin/env python3
"""76 -- the published NanoNm HeLa set merged into this study's wild-type arm.

Why this table exists
---------------------
Table S5 prints, for every tool, the depth-standardised call rate of each condition and
the coverage-matched IVT/WT ratio; the raw per-unit counts and the unions are deposited
with the call sets (2026-10-06).  Two author decisions of 2026-10-04 and 2026-10-05 shape
NanoNm's row:

1. its wild-type side quotes the per-site Nm stoichiometry published with the model,
   with the source named in the table note, instead of this study's re-analysis;
2. the published set is *merged into the wild-type arm* rather than replacing it, and the
   arm holds three samples -- HeLa_WT2, HeLa_WT3 and the published pooled set, HeLa_WT1
   being left out on the author's instruction of 2026-10-05 (see ARM_UNITS) -- over which
   every wild-type quantity of the row (rate mean with SD, union, ratio) is computed.

That is only defensible if the published set is cut with the thresholds this study
applies, because the deposited sheet is a per-site *report* (ratios down to 0.001, one
pooled count column) and not a call set:

* ``ratio >= 0.1`` -- the ratio cut the delivered call sets carry;
* ``modified reads >= 20`` -- the support floor the delivered call sets carry, and the
  only cut that makes the two sides comparable.  At ``ratio >= 0.1`` alone the published
  sheet contributes ~5x more sites than this study does, because its median support at
  that cut is 8 modified reads;
* inside a gene the published sheet itself names -- the assembly check of the S11
  package (``nanonm_coordinate_check.tsv``);
* inside this study's candidate-site set, i.e. coverage >= 10 and a reference base that
  can carry the modification.

Coordinate convention (2026-10-05)
----------------------------------
The deposited sheet is 1-based; every layer of this study is 0-based BED, which is why
``73_nanonm_published_tables.py`` keys the sheet by ``pos0 = pos - 1``.  This script does the same
and *checks* the choice against our own call set rather than assuming it: shifting the
published key the wrong way drops the overlap with our calls from tens of percent to a
fraction of one percent, and the script stops if the two directions are not separated by
``CONVENTION_FACTOR``.  (The first version of this script keyed the sheet by ``pos``; the
resulting off-by-one was invisible in the site counts because the candidate set is dense
inside genes, which is exactly why the check is now part of the script.)

What the merged arm means, and what it does not
-----------------------------------------------
The published set is one pooled run.  Its rows, cut as above, join `ARM_UNITS` -- this
study's HeLa_WT2 and HeLa_WT3 -- to form the arm Table S5 prints, three samples in all;
the author chose that pair over all three units on 2026-10-05 (HeLa_WT1 carries the
smallest call count of the four candidates, 1,630).  Three consequences are unavoidable
and are therefore printed rather than hidden:

* the published sample gets *no* replicate spread of its own.  The SD of the wild-type
  per-unit cell is computed over three samples of which two are replicates.  The interval
  the table prints is the *counting* one (2026-10-06, author decision): `rd2_bands.site_boot_ratio` redraws
  every sample x stratum count from a Poisson with its observed mean against its fixed
  candidate-site set, so no sample has to be held fixed and the interval is conditional on
  the libraries rather than a statement about their reproducibility; the unit-level
  intervals -- with the published sample held fixed (`boot_ratio_fixed`) and with every
  sample resampled -- stay in the anchors for comparison;
* the wild-type side of that one row rests on two of this study's three registered units
  (Table S10 registers all three, and the other five rows and Figure 7A use all three),
  so the note names the two and prints the all-three value beside the arm's;
* the published sample has no unmodified control, so the unmodified-IVT columns stay
  this study's own numbers and enter both arms unchanged.

What it writes
--------------
``tables/s7_nanonm_published_wt.tsv``  one row per (source, cut, unit): the published
                                    sheet at its own cut and at this study's cut, this
                                    study's own call set, and the arm's three samples
                                    with their per-band detail
``tables/s7_nanonm_published_matched.tsv``  per depth stratum: the pooled weights and the
                                    per-stratum rates of the published sample, of this
                                    study's own units, and of the unmodified-IVT units
``logs/76_anchors.tsv``              every number the Table S5 note quotes, so the gate
                                    recomputes rather than trusts

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    src/harmonisation/scripts/76_nanonm_published_wt.py
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
#: the S11 package holds the parsed published table; 76 must not re-parse the workbook
PUBLISHED = ((_RB / "analysis/nanonm_published_crosscheck/tables/nanonm_published_sites.tsv"))
SHARED = ((_RB / "analysis/nanonm_published_crosscheck/tables/nanonm_shared_sites.tsv"))
COUNTS = TAB / "s7_counts_per_replicate.tsv"
SUMMARY = TAB / "s7_counts_summary.tsv"
#: the candidate-site files, which carry the per-site coverage the depth strata band on
UNIVERSE = (_XB / "harmonisation/universe/RNA002/Human")
#: the R3-9b depth diagnosis: the coverage strata and the estimator the delivered Table S5
#: ratios use; reused rather than copied so the two cannot drift apart
COV_STRATA = ((_RB / "analysis/nonm6a_replicate_depth/evidence/cov_strata.tsv"))
RD2_RATIOS = COV_STRATA.with_name("band_matched_rate_ratios.tsv")
#: the geometry helpers of the R3-9 package, so calls and universes are read exactly as
#: the delivered counts read them
R39 = (_RB / "analysis/nonm6a_false_positives")
_spec = importlib.util.spec_from_file_location(
    "r39_build_evidence", (_RB / "analysis/nonm6a_false_positives/analysis/r39_build_evidence.py"))
r39 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(r39)
#: rd2's estimators, imported rather than copied so the delivered row and the six rows of
#: `band_matched_rate_ratios.tsv` cannot be computed by two different functions
RD2_MODULE = COV_STRATA.parents[1] / "analysis" / "rd2_matched_bands.py"
_spec2 = importlib.util.spec_from_file_location("rd2_matched_bands", RD2_MODULE)
rd2_bands = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(rd2_bands)

TOOL = "NanoNm"
MOD = "Nm"
WT_UNITS = ("HeLa_WT1", "HeLa_WT2", "HeLa_WT3")
IVT_UNITS = ("HeLa_IVT_rep1", "HeLa_IVT_rep2", "HeLa_IVT_rep3")
#: the name the published pooled run carries in the arm
PUBLISHED_SAMPLE = "published_HeLa_pooled_set"

#: this study's units plus the published pooled set -- because HeLa_WT1 carries the
#: smallest call count of the four candidates (1,630 against 3,433 / 1,648 / 3,948) and
#: the author asked for it to be left out.  The three-unit arm is still computed here and
#: anchored (`matched_ratio_this_study`), so the note can give the reader the value the
#: other five rows' design corresponds to; every subset is priced in
#: `$RNAMODBENCH_LOCAL/submission/_build/arm_drop_sample_20261005/probe_arm_variants.py`.  Note that
#: this makes the wild-type side of that one row lean on two of this study's three
#: registered units, which the note states.
ARM_UNITS = ("HeLa_WT2", "HeLa_WT3")
#: the ratio cut the delivered call sets carry
RATIO_MIN = 0.1
#: the support floor the delivered call sets carry
MIN_SUPPORT = 20
#: the candidate-site filter the whole study uses (see r39.Resources)
MIN_COVERAGE = 10
NM_BASES = "ACGT"
#: the two cuts written out: the published sheet's own, and the matched one
CUTS = {
    "published_own_ratio": (RATIO_MIN, 0),
    "matched_ratio_support": (RATIO_MIN, MIN_SUPPORT),
}
#: the matched cut has to land within this factor of this study's own rate and union,
#: otherwise the two sides are not comparable and the note must not claim they are
MAX_FACTOR = 3.0
#: the published keys, shifted the wrong way, have to overlap our calls less than this
#: much: the coordinate convention is checked, not assumed
CONVENTION_FACTOR = 10.0

#: the coverage strata of the depth diagnosis, and the pooled weights of its estimator
STRATA = ("10-19", "20-29", "30-49", "50-99", "100-199", "200-399", "400-99999")
STRATA_LO = (10, 20, 30, 50, 100, 200, 400)
#: the unit-level percentile bootstrap of the depth diagnosis (rd2_matched_bands.py)
SEED = 20260928
BOOT = 20_000


def load_published() -> tuple[pd.DataFrame, int]:
    """The published HeLa rows, keyed in this study's own 0-based site space.

    ``pos0`` and not ``pos``: the deposited sheet is 1-based and every layer of this
    study is 0-based BED (see the module docstring).  ``key_pos`` keeps the other
    convention so ``check_convention`` can prove the choice on our own call set.

    A position named under two overlapping genes appears twice with two ratios (9 rows of
    the matched cut, 63 of the ratio-only cut).  The comparison counts sites, so the first
    row of a position wins and the pre-deduplication row count is returned alongside, for
    the anchors.
    """
    if not PUBLISHED.exists():
        raise SystemExit(f"missing {PUBLISHED}; run 73_nanonm_published_tables.py first")
    df = pd.read_csv(PUBLISHED, sep="\t")
    df = df[(df.cell_line == "HeLa") & df.parse_ok & df.in_named_gene].copy()
    if df.empty:
        raise SystemExit("no parsed published HeLa rows survived the assembly filters")
    n_rows = int(len(df))
    chrom = df.chrom.astype(str)
    df["key"] = chrom + ":" + df.pos0.astype("int64").astype(str)
    df["key_pos"] = chrom + ":" + df.pos.astype("int64").astype(str)
    return df.drop_duplicates("key"), n_rows


def load_bands(unit: str) -> tuple[set[str], dict[str, str]]:
    """One unit's candidate sites, keyed and tagged with their depth stratum.

    The stratum of a candidate site is what the depth diagnosis bands on, so the
    published sample can be given the per-band rates our own units get.  The candidate
    file carries the per-site coverage that the depth bands need and ``Resources`` does
    not keep; its position set is therefore checked against the universe the R3-9 package
    loads (``main``), so the two readers cannot drift apart.
    """
    path = UNIVERSE / f"{unit}__Nm.tsv"
    if not path.exists():
        raise SystemExit(f"missing candidate-site set {path}")
    df = pd.read_csv(path, sep="\t", dtype={"chrom": str, "pos": str, "base": str,
                                            "coverage": int})
    keep = df[(df.coverage >= MIN_COVERAGE) & df.base.isin(list(NM_BASES))].copy()
    keep["key"] = keep.chrom + ":" + keep.pos
    band = pd.cut(keep.coverage, bins=list(STRATA_LO) + [10 ** 9], labels=list(STRATA),
                  right=False).astype(str)
    return set(keep.key), dict(zip(keep.key, band))


def keys_of(pos_dict: dict[str, np.ndarray]) -> set[str]:
    return {f"{c}:{p}" for c, arr in pos_dict.items() for p in arr}


def union_of(nested: list[dict[str, np.ndarray]]) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    chroms = set().union(*[set(a) for a in nested]) if nested else set()
    for chrom in chroms:
        parts = [a[chrom] for a in nested if chrom in a]
        out[chrom] = np.unique(np.concatenate(parts))
    return out


def boot_ratio(a: np.ndarray, b: np.ndarray, reps: int = BOOT) -> tuple[float, float, float]:
    """rd2_matched_bands.boot_ratio -- the estimator the delivered ratios use."""
    rng = np.random.default_rng(SEED)
    d = []
    for _ in range(reps):
        x = a[rng.integers(0, len(a), len(a))].mean()
        y = b[rng.integers(0, len(b), len(b))].mean()
        if x > 0:
            d.append(y / x)
    d = np.asarray(d)
    return float(np.mean(d)), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def boot_ratio_fixed(a: np.ndarray, b: np.ndarray, fixed: float,
                     reps: int = BOOT) -> tuple[float, float, float]:
    """`boot_ratio` with one sample held at its observed value.

    2026-10-05 (author decision): the published set cannot supply a replicate spread -- it is one pooled
    run -- so resampling it as if it were an independent replicate of this study's units adds a
    variance component that cannot be estimated from the data.  It is therefore treated as a
    fixed component of the wild-type arm: every replicate resamples this study's own units and
    adds the fixed value, which is what the note means by "conditional on the published set".
    The all-resampled interval is kept in the anchors, so the difference stays on the record.
    """
    rng = np.random.default_rng(SEED)
    d = []
    for _ in range(reps):
        x = float(np.mean(np.append(a[rng.integers(0, len(a), len(a))], fixed)))
        y = b[rng.integers(0, len(b), len(b))].mean()
        if x > 0:
            d.append(y / x)
    d = np.asarray(d)
    return float(np.mean(d)), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def check_convention(pub: pd.DataFrame, our_calls: dict[str, np.ndarray]) -> tuple[int, int]:
    """Prove that the published sheet has to be shifted by one to reach our space.

    Returns the number of matched-cut published sites that sit on one of our wild-type
    calls under the two candidate conventions.  The sheet is 1-based and our calls are
    0-based, so the correct one has to be larger by ``CONVENTION_FACTOR``; anything else
    means the two sides are not on the same assembly grid and no overlap may be quoted.
    """
    ours = keys_of(our_calls)
    sel = pub[(pub.ratio >= RATIO_MIN) & (pub.mod_count >= MIN_SUPPORT)]
    shifted = int(sel.key.isin(ours).sum())
    raw = int(sel.key_pos.isin(ours).sum())
    if not shifted >= CONVENTION_FACTOR * max(raw, 1):
        raise SystemExit(
            f"the coordinate convention of the published sheet is not the one this "
            f"script assumes: pos-1 hits {shifted} of our calls, pos hits {raw}; the two "
            f"are not separated by a factor of {CONVENTION_FACTOR:.0f} and no overlap "
            f"may be quoted")
    return shifted, raw


def main() -> None:
    t0 = time.time()
    TAB.mkdir(parents=True, exist_ok=True)
    LOG.mkdir(parents=True, exist_ok=True)
    pub, pub_rows = load_published()

    # ---- this study's own wild-type calls and candidate sets ----------------- #
    res = r39.Resources(min_cov=MIN_COVERAGE)
    inside_wt, calls_wt, uni_wt, keys_wt, band_wt = {}, {}, {}, {}, {}
    for u in WT_UNITS:
        calls = res.calls(u, MOD, TOOL)
        uni = res.universe(u, MOD)
        inside, _ = r39.inside_universe(calls, uni)
        calls_wt[u], uni_wt[u], inside_wt[u] = calls, uni, inside
        keys, bands = load_bands(u)
        if keys != keys_of(uni):
            raise SystemExit(f"{u}: the candidate-site file holds {len(keys)} sites and "
                             f"r39.Resources {r39.n_positions(uni)}; the two readers of "
                             f"one and the same candidate set disagree")
        keys_wt[u], band_wt[u] = keys, bands
    our_union = union_of(list(inside_wt.values()))
    our_union_n = r39.n_positions(our_union)

    own = pd.read_csv(COUNTS, sep="\t")
    own = own[(own.tool == TOOL) & (own.condition == "WT")]
    own_n = dict(zip(own["sample"], own.n_calls_in_universe.astype(int)))
    own_u = dict(zip(own["sample"], own.universe_n.astype(int)))
    for u in WT_UNITS:
        if r39.n_positions(inside_wt[u]) != own_n[u]:
            raise SystemExit(
                f"{u}: our call set holds {r39.n_positions(inside_wt[u])} in-universe "
                f"calls, the delivered counts table uses {own_n[u]}; this script is not "
                f"reading the delivered call set")
    summary = pd.read_csv(SUMMARY, sep="\t").set_index(["tool", "condition"])
    delivered_union = int(summary.loc[(TOOL, "WT"), "union_in_universe"])
    if our_union_n != delivered_union:
        raise SystemExit(f"our wild-type union is {our_union_n} sites, the delivered "
                         f"counts table uses {delivered_union}")
    ivt_union = int(summary.loc[(TOOL, "IVT"), "union_in_universe"])
    ivt_units = pd.read_csv(COUNTS, sep="\t")
    ivt_units = ivt_units[(ivt_units.tool == TOOL) & (ivt_units.condition == "IVT")]
    ivt_n = dict(zip(ivt_units["sample"], ivt_units.n_calls_in_universe.astype(int)))
    ivt_u = dict(zip(ivt_units["sample"], ivt_units.universe_n.astype(int)))

    hits, raw_hits = check_convention(pub, union_of(list(calls_wt.values())))
    shared = pd.read_csv(SHARED, sep="\t")
    shared_wt = set(shared[shared.condition.str.contains("WT")].key)
    if not shared_wt <= keys_of(union_of(list(calls_wt.values()))):
        raise SystemExit("the S11 shared-sites table is not a subset of our wild-type "
                         "calls; the two packages disagree about the site space")

    # ---- the published sheet at both cuts, per universe ---------------------- #
    union_u = union_of(list(uni_wt.values()))
    universe_sets = {u: keys_of(uni_wt[u]) for u in WT_UNITS}
    universe_sets["union_of_the_three_universes"] = keys_of(union_u)
    rows: list[dict[str, object]] = []
    for cut, (ratio_min, support) in CUTS.items():
        sel = pub[(pub.ratio >= ratio_min) & (pub.mod_count >= support)] if support \
            else pub[pub.ratio >= ratio_min]
        for unit, space in universe_sets.items():
            n = int(sel.key.isin(space).sum())
            rows.append({
                "source": "published (this study's thresholds)", "cut": cut, "tool": TOOL,
                "unit": unit, "sites_in_universe": n, "universe_n": len(space),
                "rate_per_1e6": round(n / len(space) * 1e6, 1),
                "n_sites_in_table": int(len(sel)),
            })
    for u in WT_UNITS:
        rows.append({
            "source": "this study", "cut": "delivered call set", "tool": TOOL, "unit": u,
            "sites_in_universe": own_n[u], "universe_n": own_u[u],
            "rate_per_1e6": round(own_n[u] / own_u[u] * 1e6, 1),
            "n_sites_in_table": own_n[u],
        })
    rows.append({
        "source": "this study", "cut": "delivered call set", "tool": TOOL,
        "unit": "union_of_the_three_units", "sites_in_universe": our_union_n,
        #: the same denominator the published union row uses, so the two rates are per
        #: million candidate sites of one and the same site set
        "universe_n": len(universe_sets["union_of_the_three_universes"]),
        "rate_per_1e6": round(our_union_n / len(universe_sets["union_of_the_three_universes"])
                              * 1e6, 1),
        "n_sites_in_table": our_union_n,
    })

    # ---- the merged wild-type arm: this study's units plus the published set -- #
    matched_sel = pub[(pub.ratio >= RATIO_MIN) & (pub.mod_count >= MIN_SUPPORT)]
    pub_in_union = matched_sel[matched_sel.key.isin(
        universe_sets["union_of_the_three_universes"])]
    pub_keys = set(pub_in_union.key)
    merged_keys = keys_of(union_of([inside_wt[u] for u in ARM_UNITS])) | pub_keys
    #: the published sample is one pooled run, so it has one number per cell and its
    #: "unit" is the union of the three candidate-site sets; the arm's own units keep
    #: their own candidate sets, exactly as the delivered counts table reports them
    sample_counts = {u: own_n[u] for u in ARM_UNITS}
    sample_counts[PUBLISHED_SAMPLE] = len(pub_keys)
    for u in ARM_UNITS:
        rows.append({
            "source": "merged wild-type arm", "cut": "matched_ratio_support", "tool": TOOL,
            "unit": u, "sites_in_universe": sample_counts[u], "universe_n": own_u[u],
            "rate_per_1e6": round(sample_counts[u] / own_u[u] * 1e6, 1),
            "n_sites_in_table": sample_counts[u],
        })
    rows.append({
        "source": "merged wild-type arm", "cut": "matched_ratio_support", "tool": TOOL,
        "unit": PUBLISHED_SAMPLE, "sites_in_universe": sample_counts[PUBLISHED_SAMPLE],
        "universe_n": len(universe_sets["union_of_the_three_universes"]),
        "rate_per_1e6": round(sample_counts[PUBLISHED_SAMPLE]
                              / len(universe_sets["union_of_the_three_universes"]) * 1e6, 1),
        "n_sites_in_table": int(len(matched_sel)),
    })
    merged_union_n = len(merged_keys)
    rows.append({
        "source": "merged wild-type arm", "cut": "matched_ratio_support", "tool": TOOL,
        "unit": "union_of_the_three_samples", "sites_in_universe": merged_union_n,
        "universe_n": len(universe_sets["union_of_the_three_universes"]),
        "rate_per_1e6": round(merged_union_n
                              / len(universe_sets["union_of_the_three_universes"]) * 1e6, 1),
        "n_sites_in_table": merged_union_n,
    })
    tab = pd.DataFrame(rows)
    tab.to_csv(TAB / "s7_nanonm_published_wt.tsv", sep="\t", index=False)
    counts3 = np.array([sample_counts[u] for u in ARM_UNITS]
                       + [sample_counts[PUBLISHED_SAMPLE]], dtype=float)
    mean3, sd3 = float(counts3.mean()), float(counts3.std(ddof=1))

    # ---- the coverage-matched ratio, with the merged wild-type arm ----------- #
    
    #: `rd2_matched_bands.py` draws -- each stratum's rate per million candidate sites,
    #: strata combined with the pooled candidate-site weights of the six HeLa units -- and
    #: only the wild-type arm changes: its per-stratum rate is the mean over the three
    #: samples of `ARM_UNITS` + the published set instead of over this study's three
    #: units.  The unmodified-IVT rates cannot
    #: move: the published table has no control.  The per-stratum candidate-site counts of
    #: the published sample are asserted against `cov_strata.tsv`, so every sample of the
    #: arm is weighted by one and the same site set.
    strata = pd.read_csv(COV_STRATA, sep="\t")
    strata = strata[(strata.mod_type == MOD) & (strata.tool == TOOL)
                    & (strata.dataset == "HeLa")]
    if strata.empty:
        raise SystemExit(f"{COV_STRATA} carries no {TOOL} HeLa rows")
    ps = strata.pivot_table(index="unit", columns="stratum", values="sites_at_risk",
                            aggfunc="sum").reindex(columns=list(STRATA), fill_value=0)
    pc = strata.pivot_table(index="unit", columns="stratum", values="calls",
                            aggfunc="sum").reindex(columns=list(STRATA), fill_value=0)
    w = ps.to_numpy().sum(axis=0).astype(float)
    w = w / w.sum()
    rate_units = np.where(ps.to_numpy() > 0, pc.to_numpy() / ps.to_numpy() * 1e6, np.nan)
    idx = {u: list(ps.index).index(u) for u in ps.index}
    r_ivt = np.nanmean(rate_units[[idx[u] for u in IVT_UNITS]], axis=0)
    r_own = np.nanmean(rate_units[[idx[u] for u in WT_UNITS]], axis=0)
    pub_band_rates, pub_band_counts = [], []
    for u in WT_UNITS:
        at_risk = ps.loc[u].to_numpy()
        own_bands = (pd.Series(list(band_wt[u].values())).value_counts()
                     .reindex(list(STRATA)).fillna(0).to_numpy())
        if not np.array_equal(own_bands, at_risk):
            raise SystemExit(
                f"{u}: the strata of the candidate-site set ({own_bands}) differ from "
                f"cov_strata.tsv ({at_risk}); the published sample cannot be weighted by "
                f"that table's site counts")
        #: the published sample inside one unit's candidate set, banded by that unit's own
        #: depth strata -- the same projection the per-unit rate of our own units uses
        band = (pd.Series([band_wt[u][k] for k in pub_keys & keys_wt[u]]).value_counts()
                .reindex(list(STRATA)).fillna(0).to_numpy())
        pub_band_rates.append(band / at_risk * 1e6)
        pub_band_counts.append(band)
    r_pub = np.nanmean(pub_band_rates, axis=0)
    r_arm3 = np.nanmean(np.vstack([rate_units[[idx[u] for u in ARM_UNITS]], r_pub]), axis=0)
    matched_own = float(np.nansum(w * r_ivt) / np.nansum(w * r_own))
    matched_pub = float(np.nansum(w * r_ivt) / np.nansum(w * r_pub))
    matched_arm3 = float(np.nansum(w * r_ivt) / np.nansum(w * r_arm3))
    #: the within-study value has to be the one the depth diagnosis already published,
    #: otherwise this script and rd2_matched_bands.py disagree about the estimator
    rd2 = pd.read_csv(RD2_RATIOS, sep="\t")
    rd2_row = rd2[(rd2.mod_type == MOD) & (rd2.tool == TOOL)].iloc[0]
    if not np.isclose(matched_own, float(rd2_row.band_matched_rate_ratio), rtol=2e-3):
        raise SystemExit(f"our own arm gives {matched_own:.4f}, the depth diagnosis "
                         f"reports {float(rd2_row.band_matched_rate_ratio):.4f}; the "
                         f"stratification of this script is not rd2's")
    #: the arm's interval: the published sample is held fixed (see boot_ratio_fixed), so the
    #: interval carries the sampling uncertainty of this study's own units alone
    per_sample_wt = np.where(np.isnan(rate_units[[idx[u] for u in ARM_UNITS]]), 0.0,
                             rate_units[[idx[u] for u in ARM_UNITS]]) @ w
    fixed_wt = float(np.where(np.isnan(r_pub), 0.0, r_pub) @ w)
    per_sample_ivt = np.where(np.isnan(rate_units[[idx[u] for u in IVT_UNITS]]), 0.0,
                              rate_units[[idx[u] for u in IVT_UNITS]]) @ w
    
    #: prints the depth-standardised count of each condition -- that rate, calls per million
    #: candidate sites, read on the row's wild-type candidate-site count
    #: (`build_supp_tables.py::base_sites`) -- instead of the raw per-unit counts and unions,
    #: whose unmodified side sits *above* the wild type for four of the six tools because the
    #: unmodified libraries are deeper.  This arm's three sample rates are the two units plus
    #: the published sample, and the contract the table now prints is asserted here, together
    #: with the identity that makes the ratio column the quotient of the two count columns.
    arm_sample_rates = np.append(per_sample_wt, fixed_wt)
    arm_rate_mean = float(arm_sample_rates.mean())
    arm_rate_sd = float(arm_sample_rates.std(ddof=1))
    ivt_rate_mean = float(per_sample_ivt.mean())
    if not arm_rate_mean > ivt_rate_mean:
        raise SystemExit(f"the printed wild-type rate {arm_rate_mean:.1f} per 1e6 is not above "
                         f"the unmodified-IVT rate {ivt_rate_mean:.1f}; Table S5's two rate "
                         f"columns would read the other way round")
    if not np.isclose(ivt_rate_mean / arm_rate_mean, matched_arm3, rtol=1e-6):
        raise SystemExit(f"the rate columns give a quotient of "
                         f"{ivt_rate_mean / arm_rate_mean:.6f}, the printed ratio is "
                         f"{matched_arm3:.6f}; they are not the same estimator")
    boot_mean, boot_lo, boot_hi = boot_ratio_fixed(per_sample_wt, per_sample_ivt, fixed_wt)
    _all_mean, all_lo, all_hi = boot_ratio(np.append(per_sample_wt, fixed_wt), per_sample_ivt)
    #: 2026-10-06 (author decision, option C): the delivered Table S5 prints a counting (site-level)
    #: interval for every row of the table, and this row's is computed the same way -- which
    #: also removes the reason the unit-level interval had to hold the published sample
    #: fixed: with only the counting error resampled there is no replicate spread to
    #: estimate, so the sample enters as its three projections' counts, redrawn on the same
    #: denominators as every other count.  Row weights follow the estimator: the arm's rate
    #: in a stratum is the mean over two of this study's units and the published sample, so
    #: each unit carries 1/3 and each of the sample's three projections 1/9.
    own_calls_ui = pc.loc[list(WT_UNITS), STRATA].to_numpy().astype(float)
    own_sites_ui = ps.loc[list(WT_UNITS), STRATA].to_numpy().astype(float)
    ivt_calls_ui = pc.loc[list(IVT_UNITS), STRATA].to_numpy().astype(float)
    ivt_sites_ui = ps.loc[list(IVT_UNITS), STRATA].to_numpy().astype(float)
    arm_sel = [WT_UNITS.index(u) for u in ARM_UNITS]
    arm_calls = np.vstack([own_calls_ui[arm_sel], np.vstack(pub_band_counts).astype(float)])
    arm_sites = np.vstack([own_sites_ui[arm_sel], own_sites_ui])
    arm_row_w = np.array([1.0 / 3] * len(arm_sel) + [1.0 / 9] * len(WT_UNITS))
    s_mean, s_lo, s_hi = rd2_bands.site_boot_ratio(
        np.vstack([arm_calls, ivt_calls_ui]), np.vstack([arm_sites, ivt_sites_ui]), w,
        len(arm_calls), np.concatenate([arm_row_w, np.full(len(IVT_UNITS), 1.0 / 3)]))
    own_s_mean, own_s_lo, own_s_hi = rd2_bands.site_boot_ratio(
        np.vstack([own_calls_ui, ivt_calls_ui]), np.vstack([own_sites_ui, ivt_sites_ui]),
        w, len(WT_UNITS))
    for label, point, lo, hi in (("the printed arm", matched_arm3, s_lo, s_hi),
                                 ("this study's three units", matched_own,
                                  own_s_lo, own_s_hi)):
        if not (lo <= point <= hi):
            raise SystemExit(f"{label}: the point estimate {point:.4f} falls outside the "
                             f"counting interval {lo:.4f}-{hi:.4f}; the two are not the "
                             f"same estimator and Table S5 must not mix them")
    #: this study's own units are the NanoNm row of band_matched_rate_ratios.tsv, so their
    #: counting interval has to be the one rd2 wrote there -- two scripts, one estimator
    if (abs(own_s_lo - float(rd2_row.site_boot_ci_lo)) > 5e-4
            or abs(own_s_hi - float(rd2_row.site_boot_ci_hi)) > 5e-4):
        raise SystemExit(
            f"this script's counting interval for this study's own units is "
            f"{own_s_lo:.4f}-{own_s_hi:.4f}, rd2_matched_bands.py reports "
            f"{float(rd2_row.site_boot_ci_lo):.4f}-{float(rd2_row.site_boot_ci_hi):.4f}")
    matched_tab = pd.DataFrame({
        "band": list(STRATA),
        "pooled_weight": [round(float(x), 6) for x in w],
        "wt_published_sample_rate_per_1e6": [round(float(x), 1) for x in r_pub],
        "wt_merged_arm_rate_per_1e6": [round(float(x), 1) for x in r_arm3],
        "wt_this_study_rate_per_1e6": [round(float(x), 1) for x in r_own],
        "ivt_this_study_rate_per_1e6": [round(float(x), 1) for x in r_ivt],
    })
    matched_tab.to_csv(TAB / "s7_nanonm_published_matched.tsv", sep="\t", index=False)

    # ---- guards and anchors -------------------------------------------------- #
    m = tab[tab.source.str.startswith("published")]
    matched = m[(m.cut == "matched_ratio_support")
                & (m.unit == "union_of_the_three_universes")].iloc[0]
    own_row = tab[tab.unit == "union_of_the_three_units"].iloc[0]
    factor = float(matched.rate_per_1e6) / float(own_row.rate_per_1e6)
    union_factor = float(matched.sites_in_universe) / float(our_union_n)
    for label, value in (("rate", factor), ("union", union_factor)):
        if not 1.0 / MAX_FACTOR <= value <= MAX_FACTOR:
            raise SystemExit(
                f"the published cut and this study differ by a factor of {value:.2f} in "
                f"{label}; the two sides are not comparable and Table S5 must not claim "
                f"they are")
    per_unit = m[(m.cut == "matched_ratio_support")
                 & (m.unit != "union_of_the_three_universes")]

    anchors = [
        ("published_rows_parsed_in_gene", pub_rows),
        ("published_sites_parsed_in_gene", int(len(pub))),
        ("published_coordinate_shift", int(pub.pos0.iloc[0] - pub.pos.iloc[0])),
        ("published_matched_cut_hits_on_our_calls_pos0", hits),
        ("published_matched_cut_hits_on_our_calls_pos", raw_hits),
        ("candidate_universe_filtered", "; ".join(f"{u}={own_u[u]}" for u in WT_UNITS)),
        ("candidate_universe_filtered_union", len(universe_sets["union_of_the_three_universes"])),
        ("published_rows_ratio_ge_0p1", int((pub.ratio >= RATIO_MIN).sum())),
        ("published_rows_ratio_ge_0p1_support_ge_20", int(len(matched_sel))),
        ("published_own_cut_union_cell", int(m[(m.unit == "union_of_the_three_universes")
                                               & (m.cut == "published_own_ratio")]
                                             .iloc[0].sites_in_universe)),
        ("published_matched_cut_union_cell", int(matched.sites_in_universe)),
        ("published_matched_cut_units_min", int(per_unit.sites_in_universe.min())),
        ("published_matched_cut_units_max", int(per_unit.sites_in_universe.max())),
        ("published_matched_cut_rate_min_per_1e6", float(per_unit.rate_per_1e6.min())),
        ("published_matched_cut_rate_max_per_1e6", float(per_unit.rate_per_1e6.max())),
        ("our_union_cell", our_union_n),
        ("our_union_cell_raw", int(summary.loc[(TOOL, "WT"), "union_raw"])),
        ("our_units", "; ".join(f"{u}={own_n[u]}" for u in WT_UNITS)),
        ("our_rate_per_1e6", float(own_row.rate_per_1e6)),
        ("matched_vs_our_rate_factor", round(factor, 3)),
        ("matched_vs_our_union_factor", round(union_factor, 3)),
        #: the arm Table S5 prints: its three samples and the cells they give
        ("merged_wt_samples", "; ".join(f"{u}={sample_counts[u]}"
                                        for u in ARM_UNITS + (PUBLISHED_SAMPLE,))),
        ("merged_wt_units_excluded", "; ".join(u for u in WT_UNITS
                                                if u not in ARM_UNITS)),
        ("merged_wt_mean_per_unit", round(mean3, 1)),
        ("merged_wt_sd_per_unit", round(sd3, 1)),
        #: the two cells the delivered Table S5 now prints for this row: the arm's
        #: depth-standardised rate, mean (SD) over its three samples, against the
        #: unmodified-IVT arm's own three units
        ("merged_wt_rate_mean_per_1e6", round(arm_rate_mean, 1)),
        ("merged_wt_rate_sd_per_1e6", round(arm_rate_sd, 1)),
        ("ivt_rate_mean_per_1e6", round(ivt_rate_mean, 1)),
        ("merged_wt_union_cell", merged_union_n),
        ("merged_wt_published_share_of_union", round(len(pub_keys) / merged_union_n, 4)),
        ("merged_ratio_matched", round(matched_arm3, 4)),
        ("merged_ratio_matched_site_boot_mean", round(s_mean, 4)),
        ("merged_ratio_matched_ci_lo", round(s_lo, 4)),
        ("merged_ratio_matched_ci_hi", round(s_hi, 4)),
        #: the convention the printed interval follows, and the two intervals it replaced,
        #: so a later rebuild cannot change the convention silently
        ("merged_ratio_matched_ci_convention",
         "site-level: Poisson counting error of the four samples, B = 20000, seed 20260928; "
         "candidate sites and strata weights fixed, conditional on the libraries"),
        ("merged_ratio_matched_ci_lo_unit_level", round(boot_lo, 4)),
        ("merged_ratio_matched_ci_hi_unit_level", round(boot_hi, 4)),
        ("merged_ratio_matched_ci_convention_unit_level",
         "published sample held fixed; this study's units resampled"),
        ("merged_ratio_matched_ci_lo_all_resampled", round(all_lo, 4)),
        ("merged_ratio_matched_ci_hi_all_resampled", round(all_hi, 4)),
        #: the two single-source arms, kept so the note can name what the merged value is
        #: made of and so the gate can check that neither has drifted
        ("matched_ratio_published_arm", round(matched_pub, 4)),
        ("matched_ratio_this_study", round(matched_own, 4)),
        ("matched_ratio_rd2_this_study", float(rd2_row.band_matched_rate_ratio)),
        ("matched_ci_this_study_site_lo", round(own_s_lo, 4)),
        ("matched_ci_this_study_site_hi", round(own_s_hi, 4)),
        ("matched_ci_this_study_lo", float(rd2_row.boot_matched_ci_lo)),
        ("matched_ci_this_study_hi", float(rd2_row.boot_matched_ci_hi)),
        #: the same contrast before the strata are matched, so the size of the
        #: composition effect between the two pipelines stays on the record
        ("raw_ratio_published_arm", round(ivt_union / float(matched.sites_in_universe), 4)),
        ("ivt_union_cell", ivt_union),
        ("ivt_union_cell_raw", int(summary.loc[(TOOL, "IVT"), "union_raw"])),
        ("ivt_units", "; ".join(f"{u}={ivt_n[u]}" for u in IVT_UNITS)),
        ("ivt_universe", "; ".join(f"{u}={ivt_u[u]}" for u in IVT_UNITS)),
    ]
    pd.DataFrame(anchors, columns=["key", "value"]).to_csv(
        LOG / "76_anchors.tsv", sep="\t", index=False)

    print(tab.to_string(index=False))
    print(f"\nconvention: pos-1 hits {hits} of our calls, pos hits {raw_hits} "
          f"(S11 shared WT sites {len(shared_wt)}, all inside our calls)")
    print(f"published matched cut: union cell {int(matched.sites_in_universe)} vs our own "
          f"{our_union_n}; rate {matched.rate_per_1e6} vs {own_row.rate_per_1e6} per 1e6 "
          f"(factor {factor:.2f})")
    print(f"arm for Table S5: samples {counts3.astype(int).tolist()} "
          f"({', '.join(ARM_UNITS)} + the published set) -> mean {mean3:.1f}, "
          f"SD {sd3:.1f}; union {merged_union_n}")
    print(f"depth-standardised rate cells: WT {arm_rate_mean:.1f} ({arm_rate_sd:.1f}) vs "
          f"unmodified IVT {ivt_rate_mean:.1f} per 1e6; quotient "
          f"{ivt_rate_mean / arm_rate_mean:.4f}")
    print(f"coverage-matched IVT/WT ratio: arm {matched_arm3:.4f} "
          f"({s_lo:.4f}-{s_hi:.4f} printed, site-level; unit-level would be "
          f"{boot_lo:.4f}-{boot_hi:.4f} with the published sample fixed and "
          f"{all_lo:.4f}-{all_hi:.4f} all-resampled), "
          f"published-only arm {matched_pub:.4f}, this study's three units "
          f"{matched_own:.4f} ({own_s_lo:.4f}-{own_s_hi:.4f} site-level)")
    print(matched_tab.to_string(index=False))
    print(f"-> {TAB / 's7_nanonm_published_wt.tsv'}")
    print(f"-> {TAB / 's7_nanonm_published_matched.tsv'}")
    print(f"-> {LOG / '76_anchors.tsv'}  [{time.time() - t0:.1f}s]")


if __name__ == "__main__":
    main()
