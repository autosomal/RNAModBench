#!/usr/bin/env python
"""Figure 4 revision, step 1: replicate-aware KL divergence recompute.

Rebuilds the numbers behind Figure 4A/B from the replicate-aware, site-level
``callsets`` export (13 m6A tool configurations x 3 species WT groups),
replacing the legacy per-condition merged ``output/<cond>/5mer/*_5mer.txt``
input that had no replicate structure at all.

Method (identical to the published Figure-4 pipeline
``Position-wise_Frequency_Differences.ipynb`` and to the R3-6 recompute
``R3-6_motif_algorithmic_bias/analysis/r3_6_callsets_kl.py``):
  * 5-mers are strand-normalised so the centre base is A (U -> revcomp);
  * per-position base frequencies form the observed PWM P;
  * KL = sum over positions of sum over bases b of
        p*log2(p/q) + (1-p)*log2((1-p)/(1-q)).

What changes vs the legacy code (reviewer R3 minor #1 asked for exactly
these three things):
  1. zero-frequency handling: the legacy silent ``np.clip(1e-10)`` is
     replaced by explicit add-epsilon (Laplace) smoothing of BOTH
     distributions, ``x'_b = (x_b + eps) / (1 + 4*eps)``,
     primary eps = 1e-4 (sensitivity table for 1e-5 / 1e-3);
  2. the background is stated explicitly:
       primary Q     = theoretical RRACH consensus PWM (published method);
       robustness Q  = per-replicate empirical candidate-site background
                       (exonic transcript-A positions with coverage >= 10,
                       ``universe/RNA002/<sp>/<rep>__m6A.tsv``);
  3. every (tool, species) value is computed per independent sequencing
     unit and summarised as mean +/- SD.  The two mouse WT studies are
     NEVER merged: only ``mES_WT`` enters the figure;
     ``mESCs_Mettl3_WT`` is computed for cross-checking but flagged
     ``in_figure = False`` (cross-study rule, config.py CROSS_STUDY_GROUPS).

A legacy-method column (clip 1e-10, RRACH background) is computed on the
same PWMs and asserted equal to the R3-6 table ``fullcount_kl.tsv``, which
is itself the bridge to the published numbers (pooled r = 0.991).

Inputs (read-only)
  harmonisation/callsets/RNA002/<Sp>/<group>_WT/m6A/<tool>/<rep>.tsv
      columns used: five_mer_raw
  harmonisation/universe/RNA002/<Sp>/<rep>__m6A.tsv
      columns used: five_mer, coverage
  R3-6_motif_algorithmic_bias/analysis/fullcount_kl.tsv   (validation)

Outputs
  ../tables/fig4_kl_per_rep.tsv       per replicate KL (3 variants) + N sites/universe
  ../tables/fig4_kl_summary.tsv       per tool x species mean +/- SD (figure rows)
  ../tables/fig4_legacy_vs_revision.tsv  clip vs add-eps pooled KL + tool ranks
  ../tables/fig4_eps_sensitivity.tsv  rank stability across eps
  (this directory) fig4_kmer_counts.tsv.gz / fig4_pwm_per_rep.npz  intermediates
                   consumed by fig4_figures.py

Env: conda ``motif_analysis`` (numpy/pandas/scipy only; no sklearn needed).
"""

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
import glob
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
TABLES = os.path.join(os.path.dirname(HERE), "tables")
ROOT = str(_RB)
CLEAN = f"{_RB}/data/callsets/RNA002"
UNIVERSE = f"{_XB}/harmonisation/universe/RNA002"
R36_KL = f"{ROOT}/analysis/motif_bias_control/analysis/fullcount_kl.tsv"

# locked analysis constants (harmonisation/common/config.py)
PRIMARY_C_MIN = 10
EPS_PRIMARY = 1e-4
EPS_SCAN = [1e-5, 1e-4, 1e-3]

SP = {"Arabidopsis": "Arabidopsis", "Mouse": "Mouse", "Human (HeLa)": "Human"}
BASES = ["A", "U", "C", "G"]
#: theoretical RRACH consensus PWM, rows = positions -2..+2, cols = BASES
RRACH = np.array([
    [0.5, 0.0, 0.0, 0.5],   # R
    [0.5, 0.0, 0.0, 0.5],   # R
    [1.0, 0.0, 0.0, 0.0],   # A
    [0.0, 0.0, 1.0, 0.0],   # C
    [0.33, 0.34, 0.33, 0.0]  # H  (A .33 / U .34 / C .33, as published)
])
COMP = str.maketrans("ACUG", "UGAC")


def norm5(s: str) -> str:
    """Strand-normalise a U-RNA 5-mer so the centre base is A."""
    if s[2] == "U":
        return s.translate(COMP)[::-1]
    return s


def counts_from_kmers(series: pd.Series) -> Counter:
    """Normalised centre-A 5-mer Counter from a raw k-mer column."""
    km = series.astype(str).str.replace("T", "U", regex=False)
    km = km[km.str.len() == 5]
    km = km.map(norm5)
    km = km[km.str[2] == "A"]
    return Counter(km)


def pwm_from_counts(cnt: Counter) -> np.ndarray:
    """5x4 per-position base-frequency matrix (column order = BASES)."""
    pos = np.zeros((5, 4))
    for km, c in cnt.items():
        for i, b in enumerate(km):
            pos[i, BASES.index(b)] += c
    return pos / pos.sum(1, keepdims=True)


def smooth(x: np.ndarray, eps: float) -> np.ndarray:
    """Add-epsilon (Laplace) smoothing of a probability row: (x+eps)/(1+4eps)."""
    return (x + eps) / (1.0 + 4.0 * eps)


def kl_binary(p: np.ndarray, q: np.ndarray) -> float:
    """Published KL form: per position, sum over bases of the two binary terms."""
    return float(np.sum(p * np.log2(p / q) + (1 - p) * np.log2((1 - p) / (1 - q))))


def kl_rrach_clip(pwm: np.ndarray) -> float:
    """Legacy method, kept only to bridge the published / R3-6 numbers."""
    tot = 0.0
    for i in range(5):
        p = np.clip(pwm[i], 1e-10, 1 - 1e-10)
        q = np.clip(RRACH[i], 1e-10, 1 - 1e-10)
        tot += kl_binary(p, q)
    return tot


def kl_rrach_addeps(pwm: np.ndarray, eps: float) -> float:
    return sum(kl_binary(smooth(pwm[i], eps), smooth(RRACH[i], eps)) for i in range(5))


def kl_emp_addeps(pwm: np.ndarray, q_pwm: np.ndarray, eps: float) -> float:
    return sum(kl_binary(smooth(pwm[i], eps), smooth(q_pwm[i], eps)) for i in range(5))


def universe_pwm(rep: str, spdir: str) -> tuple[np.ndarray, int]:
    """Empirical background Q: candidate transcript-A sites, coverage >= c_min."""
    path = f"{UNIVERSE}/{spdir}/{rep}__m6A.tsv"
    df = pd.read_csv(path, sep="\t", usecols=["five_mer", "coverage"])
    df = df[df["coverage"] >= PRIMARY_C_MIN]
    cnt = counts_from_kmers(df["five_mer"])
    return pwm_from_counts(cnt), int(sum(cnt.values()))


def main() -> None:
    per_rep, pwm_store = [], {}
    kmer_rows = []
    for spname, spdir in SP.items():
        for path in sorted(glob.glob(f"{CLEAN}/{spdir}/*_WT/m6A/*/*.tsv")):
            tool = path.split("/m6A/")[1].split("/")[0]
            rep = os.path.basename(path)[:-4]
            in_figure = not (spname == "Mouse" and rep.startswith("mESCs_"))
            site = pd.read_csv(path, sep="\t", usecols=["five_mer_raw"])
            if len(site) == 0:
                print(f"[skip] empty callset: {spname} {tool} {rep}", flush=True)
                continue
            cnt = counts_from_kmers(site["five_mer_raw"])
            if not cnt:
                print(f"[skip] no centre-A 5-mer: {spname} {tool} {rep}", flush=True)
                continue
            n_sites = int(sum(cnt.values()))
            pwm = pwm_from_counts(cnt)
            q_emp, n_uni = universe_pwm(rep, spdir)
            row = {
                "Tool": tool, "Species_Name": spname, "Replicate": rep,
                "in_figure": in_figure, "N_Sites": n_sites, "N_universe_c10": n_uni,
                "KL_rrach_clip": kl_rrach_clip(pwm),
                "KL_rrach_addeps": kl_rrach_addeps(pwm, EPS_PRIMARY),
                "KL_emp_addeps": kl_emp_addeps(pwm, q_emp, EPS_PRIMARY),
            }
            per_rep.append(row)
            pwm_store[f"{tool}|{spname}|{rep}"] = pwm
            for km, c in sorted(cnt.items()):
                kmer_rows.append((tool, spname, rep, km, c))
            print(f"[ok] {spname:14s} {tool:14s} {rep:22s} "
                  f"n={n_sites:>6d} uni={n_uni:>7d} "
                  f"clip={row['KL_rrach_clip']:8.3f} "
                  f"addeps={row['KL_rrach_addeps']:7.3f} "
                  f"emp={row['KL_emp_addeps']:7.3f}", flush=True)

    per = pd.DataFrame(per_rep)
    per.to_csv(f"{TABLES}/fig4_kl_per_rep.tsv", sep="\t", index=False)
    pd.DataFrame(kmer_rows, columns=["Tool", "Species_Name", "Replicate", "kmer", "count"]
                 ).to_csv(f"{HERE}/fig4_kmer_counts.tsv.gz", sep="\t", index=False,
                          compression="gzip")
    np.savez_compressed(f"{HERE}/fig4_pwm_per_rep.npz", **pwm_store)

    # ---------------- summary (figure rows only: replicate principle) -------
    fig = per[per["in_figure"]]
    order = (fig.groupby("Tool")["KL_rrach_addeps"].mean().sort_values().index.tolist())
    summ = (fig.groupby(["Tool", "Species_Name"])
            .agg(n_rep=("Replicate", "count"),
                 KL_rrach_mean=("KL_rrach_addeps", "mean"),
                 KL_rrach_sd=("KL_rrach_addeps", "std"),
                 KL_emp_mean=("KL_emp_addeps", "mean"),
                 KL_emp_sd=("KL_emp_addeps", "std"),
                 KL_clip_mean=("KL_rrach_clip", "mean"))
            .reset_index())
    summ["mean_rank_overall"] = summ["Tool"].map(
        {t: i + 1 for i, t in enumerate(order)})
    summ = summ.sort_values(["mean_rank_overall", "Species_Name"])
    summ.to_csv(f"{TABLES}/fig4_kl_summary.tsv", sep="\t", index=False)

    # ---------------- legacy (clip) vs revision (add-eps) -------------------
    pooled_clip = fig.groupby(["Tool", "Species_Name"])["KL_rrach_clip"].mean()
    pooled_eps = fig.groupby(["Tool", "Species_Name"])["KL_rrach_addeps"].mean()
    pooled_emp = fig.groupby(["Tool", "Species_Name"])["KL_emp_addeps"].mean()
    lvr = pd.DataFrame({"KL_clip_pooled": pooled_clip, "KL_addeps_pooled": pooled_eps,
                        "KL_empbg_pooled": pooled_emp}).reset_index()
    lvr["rank_clip"] = lvr.groupby("Species_Name")["KL_clip_pooled"].rank()
    lvr["rank_addeps"] = lvr.groupby("Species_Name")["KL_addeps_pooled"].rank()
    lvr["rank_empbg"] = lvr.groupby("Species_Name")["KL_empbg_pooled"].rank()
    lvr.to_csv(f"{TABLES}/fig4_legacy_vs_revision.tsv", sep="\t", index=False)

    # ---------------- eps sensitivity ---------------------------------------
    eps_rows = []
    keys = sorted({(r.Tool, r.Species_Name) for r in fig.itertuples()})
    for eps in EPS_SCAN:
        acc = {k: [] for k in keys}
        for r in fig.itertuples():
            pwm = pwm_store[f"{r.Tool}|{r.Species_Name}|{r.Replicate}"]
            acc[(r.Tool, r.Species_Name)].append(kl_rrach_addeps(pwm, eps))
        vals = {k: float(np.mean(v)) for k, v in acc.items()}
        eps_rows.append({"eps": eps,
                         "spearman_vs_primary": stats.spearmanr(
                             [vals[k] for k in keys],
                             [pooled_eps[k] for k in keys]).statistic,
                         "min": min(vals.values()), "max": max(vals.values())})
    pd.DataFrame(eps_rows).to_csv(f"{TABLES}/fig4_eps_sensitivity.tsv", sep="\t", index=False)

    # ---------------- validation vs the R3-6 bridge table -------------------
    r36 = pd.read_csv(R36_KL, sep="\t")
    merged = per.merge(r36[["Tool", "Species_Name", "Replicate", "KL_Divergence"]],
                       on=["Tool", "Species_Name", "Replicate"], how="left")
    missing = merged["KL_Divergence"].isna().sum()
    diff = (merged["KL_rrach_clip"] - merged["KL_Divergence"]).abs().max()
    print("\n=== validation vs R3-6 fullcount_kl.tsv (legacy clip method) ===")
    print(f"rows matched: {len(merged) - missing}/{len(merged)}; "
          f"max|diff| = {diff:.3e}  ->  {'OK' if diff < 1e-6 and missing == 0 else 'MISMATCH'}")

    # ---------------- headline rank-stability checks -------------------------
    print("\n=== rank stability ===")
    print(f"tool rank, clip vs add-eps (pooled, all species): "
          f"rho = {stats.spearmanr(lvr['KL_clip_pooled'], lvr['KL_addeps_pooled']).statistic:.4f}")
    print(f"tool rank, RRACH bg vs empirical bg (add-eps):     "
          f"rho = {stats.spearmanr(lvr['KL_addeps_pooled'], lvr['KL_empbg_pooled']).statistic:.4f}")
    print("\nlowest-KL tools (add-eps, in-figure mean over species):")
    print(fig.groupby("Tool")["KL_rrach_addeps"].mean().sort_values().head(5).round(4))
    print("\nKL range (add-eps, pooled): "
          f"{lvr['KL_addeps_pooled'].min():.3f} .. {lvr['KL_addeps_pooled'].max():.3f}")
    print(f"KL range (empirical bg, pooled): "
          f"{lvr['KL_empbg_pooled'].min():.3f} .. {lvr['KL_empbg_pooled'].max():.3f}")
    print("\neps sensitivity:"); print(pd.DataFrame(eps_rows).round(6))


if __name__ == "__main__":
    sys.exit(main())
