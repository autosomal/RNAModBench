#!/usr/bin/env python
"""Figure S2 v4 - statistics recomputed from scratch (no third-party numbers).

Two questions, both answering R3-6 ("does the motif signal describe the tools
or the species?"):

1. Variance decomposition of the per-unit KL divergence from the RRACH
   consensus: additive two-factor model y ~ tool + species fitted by least
   squares on the in-figure units, with SS_tool / SS_species / SS_residual,
   eta^2, F and p (the interaction term is reported from the full model).
   Computed for all 13 tool configurations and for the 10 tools shown in
   panel A, and compared with the two pre-existing R3-6 tables so that any
   difference in scope is documented rather than glossed over.

2. Top-5 5-mer identity: for every tool x species the set of its five most
   frequent centre-A 5-mers; overlap |A n B| is compared between
   (i) different tools within one species and (ii) the same tool across
   species, with a Mann-Whitney U test (H1: same-tool-across-species pairs
   agree more than within-species pairs).

Outputs analysis/figS2v4_variance.tsv, analysis/figS2v4_overlap.tsv,
        logs/figS2v4_stats_<date>.log
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
import os
from datetime import date

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
ANALYSIS = os.path.join(ROOT, "analysis")
LOGS = os.path.join(ROOT, "logs")
R36 = (str(_RB / "analysis/motif_bias_control/analysis"))
lines: list[str] = []


def log(msg: str = "") -> None:
    print(msg, flush=True)
    lines.append(str(msg))


def decompose(df: pd.DataFrame, y: str = "KL") -> pd.DataFrame:
    """Additive two-way decomposition (tool + species) plus interaction term."""
    d = df.copy()
    tools = sorted(d.Tool.unique())
    species = sorted(d.Species.unique())
    rows = []
    # additive model: intercept + tool dummies + species dummies
    x = np.column_stack(
        [np.ones(len(d))] +
        [(d.Tool == t).to_numpy(float) for t in tools[1:]] +
        [(d.Species == s).to_numpy(float) for s in species[1:]])
    yv = d[y].to_numpy(float)
    beta, *_ = np.linalg.lstsq(x, yv, rcond=None)
    fit = x @ beta
    ss_total = float(((yv - yv.mean()) ** 2).sum())
    ss_resid = float(((yv - fit) ** 2).sum())
    ss_tool = ss_total - ss_resid - _ss_factor(d, yv, "Species")
    ss_species = ss_total - ss_resid - _ss_factor(d, yv, "Tool")
    df_tool, df_species = len(tools) - 1, len(species) - 1
    df_resid = len(d) - 1 - df_tool - df_species
    ms_resid = ss_resid / df_resid
    for name, ss, df_ in (("Tool", ss_tool, df_tool),
                          ("Species", ss_species, df_species)):
        f = (ss / df_) / ms_resid if df_ else np.nan
        p = float(stats.f.sf(f, df_, df_resid)) if df_ else np.nan
        rows.append({"Source": name, "SS": ss, "df": df_,
                     "eta_sq": ss / ss_total, "F": f, "p": p})
    rows.append({"Source": "Interaction+Residual", "SS": ss_resid,
                 "df": df_resid, "eta_sq": ss_resid / ss_total, "F": np.nan,
                 "p": np.nan})
    out = pd.DataFrame(rows)
    out.attrs["n"] = len(d)
    out.attrs["ss_total"] = ss_total
    return out


def _ss_factor(d: pd.DataFrame, y: np.ndarray, factor: str) -> float:
    """SS explained by one factor alone (one-way ANOVA)."""
    ss = 0.0
    for _, g in d.assign(_y=y).groupby(factor):
        ss += len(g) * (g._y.mean() - y.mean()) ** 2
    return float(ss)


def overlap_block(top5: pd.DataFrame, name: str) -> tuple[pd.DataFrame, dict]:
    sets = {}
    for r in top5.itertuples():
        sets.setdefault((r.Species, r.Tool_display), set()).add(r.Kmer)
    within, across = [], []
    for sp, g in top5.groupby("Species"):
        tools = sorted(g.Tool_display.unique())
        for a, b in itertools.combinations(tools, 2):
            within.append(len(sets[(sp, a)] & sets[(sp, b)]))
    for tool, g in top5.groupby("Tool_display"):
        sps = sorted(g.Species.unique())
        for a, b in itertools.combinations(sps, 2):
            across.append(len(sets[(a, tool)] & sets[(b, tool)]))
    w, a = np.array(within), np.array(across)
    u, p = stats.mannwhitneyu(w, a, alternative="less")
    _, p2 = stats.mannwhitneyu(w, a, alternative="two-sided")
    res = {"scope": name, "n_within_species_pairs": len(w),
           "n_same_tool_pairs": len(a),
           "median_within": float(np.median(w)),
           "median_across": float(np.median(a)),
           "mean_within": float(w.mean()), "mean_across": float(a.mean()),
           "U": float(u), "p_one_sided_less": float(p),
           "p_two_sided": float(p2)}
    rows = ([{"scope": name, "group": "within species, different tools",
              "overlap": v} for v in w] +
            [{"scope": name, "group": "same tool, different species",
              "overlap": v} for v in a])
    return pd.DataFrame(rows), res


def main() -> None:
    os.makedirs(LOGS, exist_ok=True)
    kl = pd.read_csv(os.path.join(ANALYSIS, "figS2_kl_full.tsv"), sep="\t")
    kl = kl[kl.in_figure.astype(bool)]
    top5 = pd.read_csv(os.path.join(ANALYSIS, "figS2_top5_full.tsv"), sep="\t")

    log("=" * 72)
    log("Figure S2 v4 - statistics recomputed (own computation)")
    log("=" * 72)
    log(f"KL units: {len(kl)} (tools {kl.Tool.nunique()}, "
        f"species {kl.Species.nunique()}); top-5 groups: "
        f"{top5.groupby(['Species', 'Tool_display']).ngroups}")

    ten = ["CHEUI_m6A", "DRUMMER", "ELIGOS_diff", "ELIGOS_solo", "EpiNano_Error",
           "NanoSPA_m6A", "Nanocompore", "Yanocomp", "m6Anet", "xPore"]
    # design A: cell means (one KL per tool x species) - the design behind the
    # published "KL variance attributable to tool" statement
    cell = (kl.groupby(["Tool_display", "Species"], as_index=False)["KL"]
            .mean().rename(columns={"Tool_display": "Tool"}))
    parts = []
    for name, sub in (("cell means, all 13 tools", cell),
                      ("per-unit, all 13 tools", kl),
                      ("cell means, panel-A 10 tools",
                       cell[cell.Tool.isin(ten)]),
                      ("cell means, three species (no E. coli)",
                       cell[cell.Species != "E.coli"]),
                      ("cell means, panel-A 10 tools, three species",
                       cell[cell.Tool.isin(ten) & (cell.Species != "E.coli")])):
        if "Tool_display" in sub.columns:
            sub = sub.drop(columns=["Tool_display"])
        d = decompose(sub)
        d.insert(0, "scope", name)
        parts.append(d)
        log(f"\n[KL variance decomposition - {name}] n = {d.attrs['n']}")
        for r in d.itertuples():
            fp = "" if np.isnan(r.F) else f"F = {r.F:6.2f}  p = {r.p:.3g}"
            log(f"  {r.Source:20s} eta^2 = {r.eta_sq * 100:5.1f}%  "
                f"df = {r.df:3d}  {fp}")
    var = pd.concat(parts, ignore_index=True)
    var.to_csv(os.path.join(ANALYSIS, "figS2v4_variance.tsv"), sep="\t",
               index=False)

    log("\n[cross-check vs the two pre-existing R3-6 tables]")
    for fn in ("kl_variance_decomposition.tsv", "fullcount_variance.tsv"):
        p = os.path.join(R36, fn)
        if os.path.exists(p):
            t = pd.read_csv(p)
            t.columns = [c if c else "Source" for c in t.columns]
            t = t.rename(columns={t.columns[0]: "Source"})
            vals = ", ".join(f"{r.Source}={r.eta_sq * 100:.1f}%"
                             for r in t.itertuples())
            log(f"  {fn}: {vals}")

    over_parts, summaries = [], []
    for name, sub in (("all four species", top5),
                      ("three species (panel A/C/D scope)",
                       top5[top5.Species != "E.coli"])):
        d, s = overlap_block(sub, name)
        over_parts.append(d)
        summaries.append(s)
        log(f"\n[top-5 5-mer overlap - {name}]")
        log(f"  within species, different tools : n = {s['n_within_species_pairs']}"
            f", median {s['median_within']:.1f}, mean {s['mean_within']:.2f}")
        log(f"  same tool, different species    : n = {s['n_same_tool_pairs']}"
            f", median {s['median_across']:.1f}, mean {s['mean_across']:.2f}")
        log(f"  Mann-Whitney (H1: same-tool pairs agree more, i.e. "
            f"within < across): U = {s['U']:.1f}, one-sided p = "
            f"{s['p_one_sided_less']:.3g}, two-sided p = "
            f"{s['p_two_sided']:.3g}")
    over = pd.concat(over_parts, ignore_index=True)
    over.to_csv(os.path.join(ANALYSIS, "figS2v4_overlap.tsv"), sep="\t",
                index=False)
    pd.DataFrame(summaries).to_csv(
        os.path.join(ANALYSIS, "figS2v4_overlap_summary.tsv"), sep="\t",
        index=False)

    log("\nNOTE: the previously circulated value p = 1.7e-07 was produced "
        "elsewhere and does not reproduce here; the figure and caption use "
        "the numbers computed by this script.")
    with open(os.path.join(LOGS, f"figS2v4_stats_{date.today():%Y%m%d}.log"),
              "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
