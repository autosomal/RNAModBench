#!/usr/bin/env python
"""R3-6 v2: recompute the tool-vs-species motif-bias statistics from the
replicate-aware, site-level `sites_clean` export (full 5-mer counts, no
Top-5 truncation). PWM/KL follow the original Figure-4 method
(Position-wise_Frequency_Differences.ipynb): strand-normalize so the center
base is A, per-position base frequencies, KL = sum over positions of
sum p*log2(p/q)+(1-p)*log2((1-p)/(1-q)) against the theoretical RRACH PWM.

Inputs : ../../sites_v2/sites_clean/RNA002/<sp>/<sp>_WT/m6A/<tool>/*.tsv
Outputs(analysis/):
  fullcount_kl.tsv            per tool x species x replicate KL (+pooled)
  fullcount_motif_stats.tsv   pooled per tool x species: n, DRACH frac,
                              AGAC/GGAC freqs, KL
  fullcount_variance.tsv      eta^2 / F / p for tool vs species (pooled)
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
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
CLEAN = str(_RB / "data/sites_clean/RNA002")
SP = {"Arabidopsis": "Arabidopsis", "Mouse": "Mouse", "Human (HeLa)": "Human"}
BASES = ["A", "U", "C", "G"]
RRACH = np.array([
    [0.5, 0.0, 0.0, 0.5],   # A,U,C,G order: R = .5A .5G
    [0.5, 0.0, 0.0, 0.5],
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.33, 0.34, 0.33, 0.0],
])
COMP = str.maketrans("ACUG", "UGAC")


def kl_div(p, q, eps=1e-10):
    p = np.clip(p, eps, 1 - eps)
    q = np.clip(q, eps, 1 - eps)
    return np.sum(p * np.log2(p / q) + (1 - p) * np.log2((1 - p) / (1 - q)))


def norm5(s):
    """strand-normalize a U-RNA 5-mer so the center base is A"""
    if s[2] == "U":
        return s.translate(COMP)[::-1]
    return s


def pwm_and_kl(counts):
    """counts: Counter of normalized 5-mers -> (pwm, kl, total)"""
    total = sum(counts.values())
    pos = np.zeros((5, 4))
    for km, c in counts.items():
        for i, b in enumerate(km):
            pos[i, BASES.index(b)] += c
    pwm = pos / pos.sum(1, keepdims=True)
    kl = sum(kl_div(pwm[i], RRACH[i]) for i in range(5))
    return pwm, kl, total


def main():
    kl_rows, stat_rows = [], []
    for spname, spdir in SP.items():
        for path in sorted(glob.glob(f"{CLEAN}/{spdir}/*_WT/m6A/*/*.tsv")):
            tool = path.split("/m6A/")[1].split("/")[0]
            rep = os.path.basename(path)[:-4]
            df = pd.read_csv(path, sep="\t", usecols=["strand", "five_mer_raw", "drach"])
            if len(df) == 0:
                continue
            kmers = df["five_mer_raw"].str.replace("T", "U", regex=False)
            ok = kmers.str.len() == 5
            kmers = kmers[ok].map(norm5)
            kmers = kmers[kmers.str[2] == "A"]
            cnt = defaultdict(int)
            for k in kmers:
                cnt[k] += 1
            if not cnt:
                continue
            _, kl, tot = pwm_and_kl(cnt)
            drach = df.loc[ok, "drach"].mean()
            agac = sum(c for k, c in cnt.items() if k[:4] == "AGAC") / tot
            ggac = sum(c for k, c in cnt.items() if k[:4] == "GGAC") / tot
            kl_rows.append({"Tool": tool, "Species_Name": spname, "Replicate": rep,
                            "KL_Divergence": kl, "N_Sites": len(kmers)})
            stat_rows.append({"Tool": tool, "Species_Name": spname, "N_Sites": len(kmers),
                              "KL_Divergence": kl, "DRACH_frac": drach,
                              "AGAC_frac": agac, "GGAC_frac": ggac,
                              "AGAC_minus_GGAC": agac - ggac})
    kl = pd.DataFrame(kl_rows)
    kl.to_csv(os.path.join(HERE, "fullcount_kl.tsv"), sep="\t", index=False)
    st = pd.DataFrame(stat_rows).sort_values(["Tool", "Species_Name"])
    st.to_csv(os.path.join(HERE, "fullcount_motif_stats.tsv"), sep="\t", index=False)

    pooled = kl.groupby(["Tool", "Species_Name"], as_index=False)["KL_Divergence"].mean()
    grand = pooled["KL_Divergence"].mean()

    def ss_factor(f):
        gm = pooled.groupby(f)["KL_Divergence"].transform("mean")
        return ((gm - grand) ** 2).sum()

    ss_tot = ((pooled["KL_Divergence"] - grand) ** 2).sum()
    ss_t, ss_s = ss_factor("Tool"), ss_factor("Species_Name")
    ss_r = ss_tot - ss_t - ss_s
    k_t, k_s = pooled["Tool"].nunique(), pooled["Species_Name"].nunique()
    df_t, df_s, df_r = k_t - 1, k_s - 1, (k_t - 1) * (k_s - 1)
    f_t = (ss_t / df_t) / (ss_r / df_r)
    f_s = (ss_s / df_s) / (ss_r / df_r)
    var = pd.DataFrame({
        "SS": [ss_t, ss_s, ss_r], "df": [df_t, df_s, df_r],
        "eta_sq": [ss_t / ss_tot, ss_s / ss_tot, ss_r / ss_tot],
        "F": [f_t, f_s, np.nan],
        "p": [stats.f.sf(f_t, df_t, df_r), stats.f.sf(f_s, df_s, df_r), np.nan],
    }, index=["Tool", "Species", "Interaction+Residual"])
    var.to_csv(os.path.join(HERE, "fullcount_variance.tsv"))
    print(var.round(4))
    print(pooled.pivot(index="Tool", columns="Species_Name", values="KL_Divergence").round(2))
    print(st[st.Species_Name == "Arabidopsis"][["Tool", "DRACH_frac", "AGAC_minus_GGAC"]].round(3))
    # replicate stability: within tool-species CV of KL
    cv = kl.groupby(["Tool", "Species_Name"])["KL_Divergence"].agg(["mean", "std", "count"])
    cv["cv"] = cv["std"] / cv["mean"]
    print("median within-tool-species KL CV across replicates:",
          round(cv["cv"].median(), 4), " (n groups:", len(cv), ")")


if __name__ == "__main__":
    main()
