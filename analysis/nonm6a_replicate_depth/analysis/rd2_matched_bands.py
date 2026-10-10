#!/usr/bin/env python3
"""R3-9b step 2 -- coverage-band-matched WT vs IVT rate ratios.

rd1_build.py showed that the HeLa WT and IVT libraries expose an almost
identical NUMBER of candidate sites (ratio 1.005-1.010) but a very different
coverage COMPOSITION (IVT carries 2-3x more sites above coverage 100), while
every tool's per-site call rate rises by orders of magnitude with coverage.
This script separates the two effects:

  rate ratio (raw)          = (calls_IVT / sites_IVT) / (calls_WT / sites_WT)
  rate ratio (matched)      = sum_k w_k rate_IVT,k / sum_k w_k rate_WT,k
                              with w_k = pooled sites at risk in band k
  depth-composition factor  = raw / matched  (how much of the excess is only
                              the two groups sitting at different depths)

It also reports the comparability of the Curlcake synthetic null with the HeLa
cellular libraries (their candidate sites sit at coverage 10^3-10^5, so a
"FP per 10^6 candidate sites" number is not on the same footing), and rebuilds
the per-unit depth table including the Curlcake units, whose canonical names in
``samples_master.csv`` are still the SRA-era names.

Outputs: band_matched_rate_ratios.tsv, curlcake_vs_hela_depth.tsv,
         units_depth_all.tsv, per_band_rates.tsv,
         band_matched_unit_rates.tsv, band_matched_pair_ratios.tsv

Each ratio carries two intervals of the *same* point estimate (2026-10-06): ``boot_*``
resamples whole sequencing units and therefore carries the spread between libraries, while
``site_boot_*`` resamples only the counting error inside the two arms against their fixed
candidate-site sets.  The delivered Table S5 prints the counting one, the figure and every
other table keep the unit-level one.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parents[1]
EV = PKG / "evidence"
META = PKG.parents[0] / "sample_metadata" / "tables"
LABELS = ["10-19", "20-29", "30-49", "50-99", "100-199", "200-399", "400-99999"]
SEED = 20260928
BOOT = 20_000

CURLCAKE_ALIAS = {
    "Curlcake_IVT_rep1": "RNAAB089716",
    "Curlcake_IVT_rep3": "RNA081120181",
    "Curlcake_IVT_rep2_partial": "RNA081120181part",
}


def boot_ratio(a: np.ndarray, b: np.ndarray, reps: int = BOOT) -> tuple[float, float, float]:
    rng = np.random.default_rng(SEED)
    d = []
    for _ in range(reps):
        x = a[rng.integers(0, len(a), len(a))].mean()
        y = b[rng.integers(0, len(b), len(b))].mean()
        if x > 0:
            d.append(y / x)
    d = np.asarray(d)
    return float(np.mean(d)), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def site_boot_ratio(calls: np.ndarray, sites: np.ndarray, w: np.ndarray, n_wt: int,
                    row_w: np.ndarray | None = None, reps: int = BOOT,
                    seed: int = SEED) -> tuple[float, float, float]:
    """The same estimator with only the *counting* error of the two arms resampled.

    2026-10-06 (author decision, option C): Table S5 prints this interval.  Rows are samples (the wild
    type first), columns are the depth strata.  Each sample x stratum count is redrawn from
    a Poisson with its observed mean while its candidate-site count and the pooled strata
    weights stay at their observed values, so the percentile interval below is the precision
    with which the two call rates are *counted* inside these libraries.  That is a different
    quantity from ``boot_ratio`` above, which resamples whole units and so also carries the
    spread between libraries; both are kept, the table prints this one and the figure and
    the other tables keep the unit-level one.

    ``row_w`` is the weight each row carries inside its arm (equal weights by default, which
    is what this package's own ratio uses); ``76_nanonm_published_wt.py`` passes weights for
    an arm whose rows are not on one footing.
    """
    if (sites <= 0).any():
        raise SystemExit("site_boot_ratio needs a positive candidate-site count in every "
                         "sample x stratum cell; the estimator it resamples would then "
                         "differ from the printed ratio")
    if row_w is None:
        row_w = np.ones(calls.shape[0]) / calls.shape[0]
    row_w = np.asarray(row_w, dtype=float)
    rng = np.random.default_rng(seed)
    draws = rng.poisson(np.maximum(calls, 0.0), size=(reps,) + calls.shape)
    rate = draws / sites * 1e6
    per_row = rate @ w
    wt = (per_row[:, :n_wt] * row_w[:n_wt]).sum(axis=1) / row_w[:n_wt].sum()
    ivt = (per_row[:, n_wt:] * row_w[n_wt:]).sum(axis=1) / row_w[n_wt:].sum()
    d = ivt / wt
    return float(d.mean()), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def units_depth_all() -> pd.DataFrame:
    master = pd.read_csv(META / "samples_master.csv")
    rows = []
    for unit, alias in CURLCAKE_ALIAS.items():
        r = master[master["canonical_name"] == alias]
        if r.empty:
            continue
        d = r.iloc[0].to_dict()
        d["canonical_name"] = unit
        rows.append(d)
    cc = pd.DataFrame(rows)
    hela = master[master["canonical_name"].isin(
        ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3", "HeLa_IVT_rep1", "HeLa_IVT_rep2",
         "HeLa_IVT_rep3"])]
    allm = pd.concat([hela, cc], ignore_index=True).rename(columns={"canonical_name": "unit"})
    allm["dataset"] = np.where(allm["unit"].str.startswith("HeLa"), "HeLa", "Curlcake")
    uc = pd.read_csv(EV / "units_calls.tsv", sep="\t")
    agg = (uc.groupby("unit")
           .agg(sites_at_risk_max=("sites_at_risk", "max"),
                mean_cov_all_ge10=("mean_cov_all_ge10", "first"),
                median_cov_all_ge10=("median_cov_all_ge10", "first"),
                universe_rows_min=("universe_rows", "min"),
                aligned_bases_proxy=("aligned_bases_proxy", "first"))
           .reset_index())
    return allm.merge(agg, on="unit", how="left")


def main() -> None:
    s = pd.read_csv(EV / "cov_strata.tsv", sep="\t")
    u = pd.read_csv(EV / "units_calls.tsv", sep="\t")
    out_rows, band_rows = [], []
    
    #: pooled ratio above; written after the loop so the files exist even if a
    #: tool is added with no HeLa unit.
    unit_rows, pair_rows = [], []

    for (mod, tool), g in s[s.dataset == "HeLa"].groupby(["mod_type", "tool"]):
        ps = g.pivot_table(index="unit", columns="stratum", values="sites_at_risk",
                           aggfunc="sum", fill_value=0).reindex(columns=LABELS, fill_value=0)
        pc = g.pivot_table(index="unit", columns="stratum", values="calls",
                           aggfunc="sum", fill_value=0).reindex(columns=LABELS, fill_value=0)
        wt = [x for x in ps.index if x in ("HeLa_WT1", "HeLa_WT2", "HeLa_WT3")]
        iv = [x for x in ps.index if x.startswith("HeLa_IVT")]
        rate = np.where(ps.to_numpy() > 0, pc.to_numpy() / ps.to_numpy() * 1e6, np.nan)
        w = ps.to_numpy().sum(axis=0).astype(float)          # pooled sites = weights
        w = w / w.sum()
        r_wt = np.nanmean(rate[[list(ps.index).index(x) for x in wt]], axis=0)
        r_iv = np.nanmean(rate[[list(ps.index).index(x) for x in iv]], axis=0)
        matched = (np.nansum(w * r_iv) / np.nansum(w * r_wt)) if np.nansum(w * r_wt) else np.nan

        # per-unit matched ratio (keeps the replicate spread for bootstrapping)
        idx_iv = [list(ps.index).index(x) for x in iv]
        idx_wt = [list(ps.index).index(x) for x in wt]
        per_unit_rate_iv = np.where(np.isnan(rate[idx_iv]), 0.0, rate[idx_iv]) @ w
        per_unit_rate_wt = np.where(np.isnan(rate[idx_wt]), 0.0, rate[idx_wt]) @ w
        b_mean, b_lo, b_hi = boot_ratio(per_unit_rate_wt, per_unit_rate_iv)

        #: 2026-10-06 (author decision, option C): the interval the delivered Table S5 prints.  Same
        #: estimator, same strata weights, but only the counting error of the two arms is
        #: resampled (see site_boot_ratio), so it is narrower than the unit-level one by
        #: exactly the between-library spread.  The point estimate does not move, and the
        #: unit-level interval stays in this table for every consumer that wants the
        #: library-level quantity (the R3-9b figure, the robustness checks).
        calls_ui = pc.loc[wt + iv, LABELS].to_numpy().astype(float)
        sites_ui = ps.loc[wt + iv, LABELS].to_numpy().astype(float)
        s_mean, s_lo, s_hi = site_boot_ratio(calls_ui, sites_ui, w, len(wt))
        if not (s_lo <= matched <= s_hi):
            raise SystemExit(
                f"{mod}/{tool}: the printed point estimate {matched:.4f} falls outside the "
                f"counting interval {s_lo:.4f}-{s_hi:.4f}; the two are not the same "
                f"estimator and Table S5 must not mix them")

        
        # pooled ratio (identical pooled-band weights ``w``), so they are written out
        # instead of being used only for the bootstrap.  They give Figure S8D one
        # mark per depth-matched WT/IVT replicate pair; the earlier
        # ``nearest_depth_pairs.tsv::std_density_ratio`` is NOT that quantity (it
        # divides each unit's own ``standardised_density_per_1e6``, i.e. it uses the
        # per-unit standardisation weights), so reusing it would put a second
        # estimator on the same page.  Pairing follows rd1_build.py: both conditions
        # are sorted by sites at risk and paired index by index, rank 1 = shallowest.
        unit_rate = dict(zip(iv, per_unit_rate_iv), **dict(zip(wt, per_unit_rate_wt)))
        for cond_name, units_list, rates in (("WT", wt, per_unit_rate_wt),
                                             ("IVT", iv, per_unit_rate_iv)):
            for unit_name, rate_u in zip(units_list, rates):
                unit_rows.append({
                    "mod_type": mod, "tool": tool, "unit": unit_name,
                    "condition": cond_name,
                    "sites_at_risk": int(ps.loc[unit_name].to_numpy().sum()),
                    "n_calls_in_universe": int(pc.loc[unit_name].to_numpy().sum()),
                    "band_matched_rate_per_1e6": float(rate_u),
                })
        wt_by_depth = sorted(((float(ps.loc[x].to_numpy().sum()), x) for x in wt))
        iv_by_depth = sorted(((float(ps.loc[x].to_numpy().sum()), x) for x in iv))
        for k, ((s_wt, u_wt), (s_iv, u_iv)) in enumerate(zip(wt_by_depth, iv_by_depth), 1):
            pair_rows.append({
                "mod_type": mod, "tool": tool, "pair_rank_by_depth": k,
                "wt_unit": u_wt, "ivt_unit": u_iv,
                "wt_sites_at_risk": int(s_wt), "ivt_sites_at_risk": int(s_iv),
                "depth_ratio_ivt_over_wt": s_iv / s_wt if s_wt else np.nan,
                "wt_rate_per_1e6": float(unit_rate[u_wt]),
                "ivt_rate_per_1e6": float(unit_rate[u_iv]),
                "matched_rate_ratio_ivt_over_wt": (float(unit_rate[u_iv])
                                                   / float(unit_rate[u_wt]))
                if unit_rate[u_wt] else np.nan,
            })

        raw_w = pc.loc[wt].to_numpy().sum() / ps.loc[wt].to_numpy().sum() * 1e6
        raw_i = pc.loc[iv].to_numpy().sum() / ps.loc[iv].to_numpy().sum() * 1e6
        raw_ratio = raw_i / raw_w
        out_rows.append({
            "mod_type": mod, "tool": tool,
            "wt_calls": int(pc.loc[wt].to_numpy().sum()),
            "ivt_calls": int(pc.loc[iv].to_numpy().sum()),
            "count_ratio_ivt_over_wt": pc.loc[iv].to_numpy().sum() / pc.loc[wt].to_numpy().sum(),
            "raw_density_ratio": raw_ratio,
            "band_matched_rate_ratio": matched,
            "depth_composition_factor": raw_ratio / matched if matched else np.nan,
            "boot_matched_mean": b_mean, "boot_matched_ci_lo": b_lo, "boot_matched_ci_hi": b_hi,
            "site_boot_mean": s_mean, "site_boot_ci_lo": s_lo, "site_boot_ci_hi": s_hi,
            "share_sites_cov_ge100_WT": float(ps.loc[wt, ["100-199", "200-399",
                                                          "400-99999"]].to_numpy().sum()
                                              / ps.loc[wt].to_numpy().sum()),
            "share_sites_cov_ge100_IVT": float(ps.loc[iv, ["100-199", "200-399",
                                                           "400-99999"]].to_numpy().sum()
                                               / ps.loc[iv].to_numpy().sum()),
            "rate_slope_log10_cov_per_decade": None,
        })
        # band-level detail
        for k, lab in enumerate(LABELS):
            band_rows.append({
                "mod_type": mod, "tool": tool, "band": lab,
                "sites_WT": int(ps.loc[wt, lab].to_numpy().sum()),
                "sites_IVT": int(ps.loc[iv, lab].to_numpy().sum()),
                "calls_WT": int(pc.loc[wt, lab].to_numpy().sum()),
                "calls_IVT": int(pc.loc[iv, lab].to_numpy().sum()),
                "rate_WT_per_1e6": r_wt[k], "rate_IVT_per_1e6": r_iv[k],
                "band_rate_ratio": (r_iv[k] / r_wt[k]) if r_wt[k] else np.nan,
                "pooled_weight": w[k],
            })

    band = pd.DataFrame(out_rows)
    # how steeply does the call rate track coverage?  geometric mid-point of each
    # band vs per-site rate, pooled over the six HeLa units -> log-log slope
    BAND_MID = np.sqrt(np.array([10 * 20, 20 * 30, 30 * 50, 50 * 100, 100 * 200,
                                 200 * 400, 400 * 1000], dtype=float))
    tools_order = [(m, t) for m, t in s[s.dataset == "HeLa"][["mod_type", "tool"]]
                   .drop_duplicates().itertuples(index=False, name=None)]
    slopes = {}
    for mod, tool in tools_order:
        g = s[(s.dataset == "HeLa") & (s.mod_type == mod) & (s.tool == tool)]
        calls = g.pivot_table(index="unit", columns="stratum", values="calls",
                              aggfunc="sum", fill_value=0).reindex(columns=LABELS, fill_value=0)
        sites = g.pivot_table(index="unit", columns="stratum", values="sites_at_risk",
                              aggfunc="sum", fill_value=0).reindex(columns=LABELS, fill_value=0)
        rate = calls.to_numpy().sum(axis=0) / np.maximum(sites.to_numpy().sum(axis=0), 1) * 1e6
        ok = sites.to_numpy().sum(axis=0) > 1000
        if ok.sum() > 2:
            x = np.log10(BAND_MID[ok])
            y = np.log10(rate[ok] + 1.0)
            slopes[(mod, tool)] = float(np.polyfit(x, y, 1)[0])
    band["rate_slope_log10_cov_per_decade"] = [
        slopes.get((m, t), np.nan) for m, t in
        zip(band["mod_type"], band["tool"])]
    band.to_csv(EV / "band_matched_rate_ratios.tsv", sep="\t", index=False)
    pd.DataFrame(band_rows).to_csv(EV / "per_band_rates.tsv", sep="\t", index=False)
    unit_tab = pd.DataFrame(unit_rows)
    pairs = pd.DataFrame(pair_rows)
    unit_tab.to_csv(EV / "band_matched_unit_rates.tsv", sep="\t", index=False)
    pairs.to_csv(EV / "band_matched_pair_ratios.tsv", sep="\t", index=False)

    # ---- Curlcake null vs HeLa libraries: are the densities comparable? -----
    rows = []
    for (mod, tool), g in u.groupby(["mod_type", "tool"]):
        for ds, sub in g.groupby("dataset"):
            indep = sub[~sub["unit"].str.endswith("_partial")]
            rows.append({
                "mod_type": mod, "tool": tool, "dataset": ds,
                "n_independent_units": len(indep),
                "sites_at_risk_range": f"{int(indep['sites_at_risk'].min()):,}-"
                                       f"{int(indep['sites_at_risk'].max()):,}",
                "median_cov_of_sites_at_risk": float(indep["median_cov_all_ge10"].median()),
                "mean_cov_of_calls": float(indep["mean_cov_of_calls"].mean()),
                "median_cov_of_calls": float(indep["median_cov_calls"].median()),
                "density_per_1e6_median": float(indep["density_per_1e6"].median()),
                "band_matched_density_per_1e6": float(
                    indep["standardised_density_per_1e6"].median()),
            })
    pd.DataFrame(rows).to_csv(EV / "curlcake_vs_hela_depth.tsv", sep="\t", index=False)
    units_depth_all().to_csv(EV / "units_depth_all.tsv", sep="\t", index=False)

    pd.set_option("display.width", 240)
    print("\n== band-matched rate ratios (IVT/WT) ==")
    print(band[["mod_type", "tool", "count_ratio_ivt_over_wt", "raw_density_ratio",
                "band_matched_rate_ratio", "depth_composition_factor",
                "boot_matched_ci_lo", "boot_matched_ci_hi",
                "site_boot_ci_lo", "site_boot_ci_hi",
                "share_sites_cov_ge100_WT", "share_sites_cov_ge100_IVT",
                "rate_slope_log10_cov_per_decade"]].round(3).to_string(index=False))
    print("\n== per-band detail ==")
    print(pd.DataFrame(band_rows).pivot_table(index="tool", columns="band",
          values="band_rate_ratio").round(2).to_string())
    
    # mean of the three pair ratios is printed next to the pooled bootstrap mean.
    #: the two are different functionals (ratio of the two resampled condition means
    #: versus the mean of three single-pair ratios), so a small gap is expected and
    #: is bounded by the figure build, not here.
    chk = (pairs.groupby("tool")["matched_rate_ratio_ivt_over_wt"]
           .agg(["mean", "min", "max"])
           .join(band.set_index("tool")["boot_matched_mean"]))
    chk["mean_minus_bootstrap"] = chk["mean"] - chk["boot_matched_mean"]
    print("\n== pair ratios vs pooled bootstrap mean ==")
    print(chk.round(4).to_string())


if __name__ == "__main__":
    main()
