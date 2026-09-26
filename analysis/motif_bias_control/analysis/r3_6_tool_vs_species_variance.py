#!/usr/bin/env python
"""R3-6: how much of the observed motif variation is explained by TOOL
versus SPECIES? Uses existing project outputs only (no re-calling).

Inputs
  $RNAMODBENCH_LOCAL/code/code/next_postprocessing/Total_new/combined_kl_data.csv
      per-tool per-species KL divergence of the 5-mer PWM vs RRACH
  $RNAMODBENCH_LOCAL/code/code/next_postprocessing/Total_new/All_Species_Top5_5mer_Data.csv
      per-tool per-species top-20 strand-aware 5-mers with frequencies

Outputs (into this folder's analysis/ dir)
  kl_variance_decomposition.tsv   two-way ANOVA + eta^2 for tool/species/interaction
  kl_tool_across_species.tsv      per-tool KL range across species
  top1_motif_concordance.tsv      within-species (tool-vs-tool) vs
                                  between-species (same-tool) top-5 overlap
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
import os
from math import nan
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TOTAL = str(_XB / "legacy/next_postprocessing/Total_new")


def anova_and_eta(df):
    """Balanced two-way main-effects decomposition (one obs per cell, so
    residual SS = tool x species interaction)."""
    import numpy as np
    from scipy import stats
    grand = df["KL_Divergence"].mean()
    ss_tot = ((df["KL_Divergence"] - grand) ** 2).sum()
    def ss_factor(f):
        gm = df.groupby(f)["KL_Divergence"].transform("mean")
        return ((gm - grand) ** 2).sum()
    ss_tool, ss_species = ss_factor("Tool"), ss_factor("Species_Name")
    ss_res = ss_tot - ss_tool - ss_species
    k_t = df["Tool"].nunique(); k_s = df["Species_Name"].nunique()
    df_t, df_s, df_r = k_t - 1, k_s - 1, (k_t - 1) * (k_s - 1)
    f_t = (ss_tool / df_t) / (ss_res / df_r)
    f_s = (ss_species / df_s) / (ss_res / df_r)
    out = pd.DataFrame({
        "SS": [ss_tool, ss_species, ss_res],
        "df": [df_t, df_s, df_r],
        "eta_sq": [ss_tool / ss_tot, ss_species / ss_tot, ss_res / ss_tot],
        "F": [f_t, f_s, nan],
        "p": [stats.f.sf(f_t, df_t, df_r), stats.f.sf(f_s, df_s, df_r), nan],
    }, index=["Tool", "Species", "Interaction+Residual"])
    return out


def main():
    kl = pd.read_csv(os.path.join(TOTAL, "combined_kl_data.csv"))
    tab = anova_and_eta(kl)
    tab.to_csv(os.path.join(HERE, "kl_variance_decomposition.tsv"))

    wide = kl.pivot(index="Tool", columns="Species_Name", values="KL_Divergence")
    med = wide.median(axis=1)
    wide["range"] = wide.max(axis=1) - wide.min(axis=1)
    wide["across_species_range/median"] = wide["range"] / med
    wide.sort_values("range").to_csv(os.path.join(HERE, "kl_tool_across_species.tsv"))

    top = pd.read_csv(os.path.join(TOTAL, "All_Species_Top5_5mer_Data.csv"))
    t5 = (top.sort_values(["Species", "Tool", "Rank"])
             .groupby(["Species", "Tool"])["5mer"]
             .apply(lambda s: set(s.head(5))))
    rows = []
    for (s1, t1), a in t5.items():
        for (s2, t2), b in t5.items():
            if (s1, t1) >= (s2, t2):
                continue
            kind = "within_species_diff_tool" if s1 == s2 else (
                "same_tool_diff_species" if t1 == t2 else "diff_both")
            inter = len(a & b)
            rows.append({"kind": kind, "a": f"{s1}|{t1}", "b": f"{s2}|{t2}", "overlap5": inter})
    ov = pd.DataFrame(rows)
    res = ov.groupby("kind")["overlap5"].agg(["mean", "count"])
    res.to_csv(os.path.join(HERE, "top1_motif_concordance.tsv"))

    gg = (top[top["5mer"].isin(["GGACU", "GGACA", "AGACA", "AGACU"])]
          .assign(kind=lambda x: x["5mer"].str[:2].map({"GG": "GGAC", "AG": "AGAC"}))
          .pivot_table(index=["Tool"], columns=["Species", "kind"],
                       values="Relative_Frequency", aggfunc="sum"))
    gg.columns = [f"{c1}_{c2}" for c1, c2 in gg.columns]
    gg = gg.fillna(0.0)
    for sp in ["Arabidopsis", "Mouse", "Human (HeLa)"]:
        a, b = gg[f"{sp}_AGAC"], gg[f"{sp}_GGAC"]
        gg[f"{sp}_AGAC_minus_GGAC"] = a - b
    gg.to_csv(os.path.join(HERE, "ggac_agac_within_tool.tsv"))
    print(res)
    print(wide.sort_values("range").round(2))


if __name__ == "__main__":
    main()
