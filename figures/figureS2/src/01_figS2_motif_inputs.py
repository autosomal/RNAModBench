#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Supplementary Figure S2 (revision) - full-count motif inputs from sites_clean.

Computes, per species x tool x replicate (sites_clean RNA002 WT groups):
  * strand-normalized full-count 5-mer spectrum (center base A after norm)
  * KL(position-wise PWM || theoretical RRACH PWM), two zero-frequency
    treatments:
      - KL        : add-eps smoothing x=(x+eps)/(1+4*eps), eps=1e-4
                    (manuscript Fig. 4 revision convention; PRIMARY)
      - KL_clip   : symmetric clip at 1e-10 (R3-6 legacy-clip convention;
                    reconciliation target of R3-6 fullcount_kl.tsv)
  * Top-5 motifs ranked by per-replicate equal-weight mean frequency
    (Fig. 4 revision convention), with pooled-count frequencies kept for
    cross-checking
  * per-replicate 5x4 PWM (npz) and the per tool x species
    replicate-equal-weight mean PWM (long TSV) -- the exact matrices the
    panel-A sequence logos of the redrawn Figure S2 are drawn from
    (added 2026-09-20, v2 published-layout redraw).

Mouse cross-study rule: figure uses mES_WT only (in_figure=False for
mESCs_Mettl3_WT); both studies are kept in the tables.

Reconciliations (logged to analysis/reconciliation_20260920.log):
  1. per-rep KL_clip vs R3-6 fullcount_kl.tsv           (rel tol 1e-9, must pass)
  2. pooled KL_clip (legacy 10-tool subset, 3 species) vs
     Total_new/combined_kl_data.csv                      (Pearson r, expect ~0.991)
  3. Top-5 sets vs R3-6 fullcount_top5.tsv and legacy
     All_Species_Top5_5mer_Data.csv                      (overlap report)
  4. per-rep KL (add-eps + clip) and per-rep PWM vs the Figure-4 revision
     products fig4_kl_per_rep.tsv / fig4_pwm_per_rep.npz (rel tol 1e-9,
     atol 1e-12, must pass): the 10 logos of Figure S2 and the 3 logos of
     Figure 4B must come from one and the same computation.

Outputs (analysis/):
  figS2_kl_full.tsv, figS2_kl_pooled.tsv, figS2_top5_full.tsv,
  figS2_ecoli_coverage.tsv, figS2_pwm_per_rep.npz, figS2_pwm_mean.tsv
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
import re
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

CLEAN = str(_RB / "data/sites_clean/RNA002")
HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.normpath(os.path.join(HERE, "..", "analysis"))
R36 = str(_RB / "analysis/motif_bias_control/analysis")
LEGACY_ROOT = str(_XB / "legacy/next_postprocessing")
FIG4 = str(_RB / "figures/figure4")
FIG4_PER_REP = os.path.join(FIG4, "tables", "fig4_kl_per_rep.tsv")
FIG4_PWM_NPZ = os.path.join(FIG4, "analysis", "fig4_pwm_per_rep.npz")

# species -> WT group used for KL (Panel A) and Top-5 (Panel B)
GROUPS = [
    ("Arabidopsis", "Arabidopsis_WT"),
    ("Mouse", "Mouse_WT"),
    ("Human", "HeLa_WT"),
    ("E.coli", "E.coli_WT"),
]
IVT_COVERAGE = [("E.coli", "E.coli_IVT")]  # coverage reference only
KL_SPECIES = ["Arabidopsis", "Mouse", "Human"]  # Panel A species (E.coli = Panel B only)

# Mouse cross-study rule: Panel display uses mES_WT (studyB) only
NOT_IN_FIGURE_REPS = {"mESCs_Mettl3_WT"}

# sites_clean dir name -> published display name (legacy sup2/Fig4 names)
DISPLAY = {"ELIGOS2_diff": "ELIGOS_diff", "ELIGOS2_solo": "ELIGOS_solo", "yanocomp": "Yanocomp"}
REVERSE_DISPLAY = {v: k for k, v in DISPLAY.items()}

# S2 species label -> Fig.4 revision species label (reconciliation 4 only)
FIG4_SPNAME = {"Arabidopsis": "Arabidopsis", "Mouse": "Mouse", "Human": "Human (HeLa)"}

BASES = ["A", "U", "C", "G"]
RRACH = np.array([
    [0.5, 0.0, 0.0, 0.5],   # R = .5A/.5G
    [0.5, 0.0, 0.0, 0.5],   # R
    [1.0, 0.0, 0.0, 0.0],   # A
    [0.0, 0.0, 1.0, 0.0],   # C
    [0.33, 0.34, 0.33, 0.0],  # H  (A .33 / U .34 / C .33, as published;
                              # identical to fig4_revision & R3-6 constants)
])
COMP = str.maketrans("ACUG", "UGAC")
RRACH_RE = re.compile(r"^[AG][AG]AC[ACU]$")
EPS_CLIP = 1e-10
EPS_ADDEPS = 1e-4

LOG_PATH = os.path.join(ANALYSIS, "reconciliation_20260920.log")
_log_lines = []


def log(msg=""):
    print(msg, flush=True)
    _log_lines.append(str(msg))


def norm5(s):
    """strand-normalize a U-RNA 5-mer so the center base is A"""
    return s.translate(COMP)[::-1] if s[2] == "U" else s


def pwm_from_counts(counts):
    pos = np.zeros((5, 4))
    for km, c in counts.items():
        for i, b in enumerate(km):
            pos[i, BASES.index(b)] += c
    return pos / pos.sum(axis=1, keepdims=True)


def kl_clip(p, q, eps=EPS_CLIP):
    p = np.clip(p, eps, 1 - eps)
    q = np.clip(q, eps, 1 - eps)
    return float(np.sum(p * np.log2(p / q) + (1 - p) * np.log2((1 - p) / (1 - q))))


def kl_addeps(p, q, eps=EPS_ADDEPS):
    p = (p + eps) / (1 + 4 * eps)
    q = (q + eps) / (1 + 4 * eps)
    return float(np.sum(p * np.log2(p / q) + (1 - p) * np.log2((1 - p) / (1 - q))))


def read_kmers(path):
    """return (n_rows, Counter of normalized center-A 5-mers)"""
    df = pd.read_csv(path, sep="\t", usecols=["five_mer_raw"], dtype=str)
    s = df["five_mer_raw"].dropna().str.replace("T", "U", regex=False)
    s = s[s.str.len() == 5].map(norm5)
    s = s[s.str[2] == "A"]
    return len(df), Counter(s)


def load_legacy_kl():
    for sub in ("Total_new", "Total_new1"):
        p = os.path.join(LEGACY_ROOT, sub, "combined_kl_data.csv")
        if os.path.exists(p):
            df = pd.read_csv(p)
            log(f"[legacy-KL] source: {p}  rows={len(df)}")
            return df, sub
    raise FileNotFoundError("combined_kl_data.csv not found under Total_new/Total_new1")


def norm_species(s):
    s = str(s).strip().lower()
    if "hela" in s or s == "human":
        return "Human"
    if "mouse" in s:
        return "Mouse"
    if "arab" in s:
        return "Arabidopsis"
    if "coli" in s:
        return "E.coli"
    return s


def main():
    os.makedirs(ANALYSIS, exist_ok=True)
    log("=" * 72)
    log("FigS2 revision - full-count motif inputs  (2026-09-20)")
    log(f"sites_clean root: {CLEAN}")
    log("=" * 72)

    full_rows, cov_rows = [], []
    spectra = defaultdict(dict)  # (sp, tool) -> {rep: Counter}  (WT groups only)
    pwm_store = {}               # "tool|species|rep" -> 5x4 PWM (all species)

    for sp, cond in GROUPS + IVT_COVERAGE:
        paths = sorted(glob.glob(os.path.join(CLEAN, sp, cond, "m6A", "*", "*.tsv")))
        tools = sorted({p.split("/m6A/")[1].split("/")[0] for p in paths})
        log(f"\n[{sp}/{cond}] {len(paths)} files, {len(tools)} tools: {', '.join(tools)}")
        for p in paths:
            tool = p.split("/m6A/")[1].split("/")[0]
            rep = os.path.basename(p)[:-4]
            n_rows, cnt = read_kmers(p)
            if sp == "E.coli":  # coverage reference (WT + IVT)
                cov_rows.append({"Species": sp, "Condition": cond, "Tool": tool,
                                 "Tool_display": DISPLAY.get(tool, tool),
                                 "Replicate": rep, "n_rows": n_rows,
                                 "N_Sites_centerA": sum(cnt.values())})
            if cond == "E.coli_IVT":
                continue
            in_figure = rep not in NOT_IN_FIGURE_REPS
            if sum(cnt.values()) == 0:
                log(f"  WARNING empty spectrum: {sp}/{cond}/{tool}/{rep}")
                continue
            mat = pwm_from_counts(cnt)
            full_rows.append({
                "Species": sp, "Condition": cond, "Tool": tool,
                "Tool_display": DISPLAY.get(tool, tool), "Replicate": rep,
                "KL": kl_addeps(mat, RRACH), "KL_clip": kl_clip(mat, RRACH),
                "N_Sites": sum(cnt.values()), "n_rows": n_rows,
                "in_figure": in_figure,
            })
            spectra[(sp, tool)][rep] = cnt
            pwm_store[f"{tool}|{sp}|{rep}"] = mat

    kl = pd.DataFrame(full_rows)
    kl.to_csv(os.path.join(ANALYSIS, "figS2_kl_full.tsv"), sep="\t", index=False)
    log(f"\nfigS2_kl_full.tsv: {len(kl)} rows "
        f"({kl['Species'].nunique()} species x {kl['Tool'].nunique()} tools)")

    # ---------------- pooled table ----------------
    legacy_kl, legacy_sub = load_legacy_kl()
    legacy_kl["Species_norm"] = legacy_kl["Species_Name"].map(norm_species)
    legacy_tools = sorted(legacy_kl["Tool"].unique())  # published display names
    log(f"[legacy-KL] tools ({len(legacy_tools)}): {', '.join(legacy_tools)}")
    log(f"[legacy-KL] species: {sorted(legacy_kl['Species_norm'].unique())}")

    pooled_rows = []
    for (sp, tool), reps in sorted(spectra.items()):
        reps_infig = [r for r in reps if r not in NOT_IN_FIGURE_REPS]
        kl_rows_infig = kl[(kl.Species == sp) & (kl.Tool == tool) & kl.in_figure]
        kl_rows_all = kl[(kl.Species == sp) & (kl.Tool == tool)]
        n = len(kl_rows_infig)
        mean = kl_rows_infig["KL"].mean() if n else np.nan
        sd = kl_rows_infig["KL"].std(ddof=1) if n > 1 else np.nan
        pooled_rows.append({
            "Species": sp, "Tool": tool, "Tool_display": DISPLAY.get(tool, tool),
            "N_reps_in_figure": n, "Replicates_in_figure": ";".join(sorted(reps_infig)),
            "KL_mean": mean, "KL_sd": sd,
            "KL_clip_mean_in_figure": kl_rows_infig["KL_clip"].mean() if n else np.nan,
            "KL_clip_mean_all_reps": kl_rows_all["KL_clip"].mean(),
            "N_Sites_in_figure": int(kl_rows_infig["N_Sites"].sum()),
            "in_panelA": bool(sp in KL_SPECIES and DISPLAY.get(tool, tool) in legacy_tools),
        })
    pooled = pd.DataFrame(pooled_rows)
    pooled.to_csv(os.path.join(ANALYSIS, "figS2_kl_pooled.tsv"), sep="\t", index=False)
    log(f"figS2_kl_pooled.tsv: {len(pooled)} rows; Panel A subset (10-tool legacy set, "
        f"3 species): {int(pooled.in_panelA.sum())}")

    # ---------------- Top-5 (full count, per-rep equal-weight mean frequency) ----
    top_rows = []
    for (sp, tool), reps in sorted(spectra.items()):
        reps_infig = [r for r in reps if r not in NOT_IN_FIGURE_REPS]
        if not reps_infig:
            continue
        totals = {r: sum(reps[r].values()) for r in reps_infig}
        repmix, pooled_cnt = defaultdict(float), defaultdict(int)
        for r in reps_infig:
            t = totals[r]
            for km, c in reps[r].items():
                repmix[km] += c / t
                pooled_cnt[km] += c
        for km in repmix:
            repmix[km] /= len(reps_infig)
        top = sorted(repmix.items(), key=lambda x: (-x[1], -pooled_cnt[x[0]], x[0]))[:5]
        tot_repmix = sum(repmix.values())
        for rank, (km, f) in enumerate(top, 1):
            top_rows.append({
                "Species": sp, "Tool": tool, "Tool_display": DISPLAY.get(tool, tool),
                "Rank": rank, "Kmer": km,
                "RelFreq_repmix": f, "RelFreq_repmix_of_centerA": f / tot_repmix,
                "Count_pooled": pooled_cnt[km],
                "RelFreq_pooled": pooled_cnt[km] / sum(pooled_cnt.values()),
                "Is_RRACH": bool(RRACH_RE.match(km)),
            })
    top5 = pd.DataFrame(top_rows)
    top5.to_csv(os.path.join(ANALYSIS, "figS2_top5_full.tsv"), sep="\t", index=False)
    log(f"figS2_top5_full.tsv: {len(top5)} rows over "
        f"{top5.groupby(['Species', 'Tool']).ngroups} species x tool groups")

    # ---------------- PWM export (panel-A sequence logos) ----------------
    np.savez_compressed(os.path.join(ANALYSIS, "figS2_pwm_per_rep.npz"), **pwm_store)
    pwm_rows = []
    for (sp, tool), reps in sorted(spectra.items()):
        reps_infig = sorted(r for r in reps if r not in NOT_IN_FIGURE_REPS)
        if not reps_infig:
            continue
        mean = np.mean([pwm_store[f"{tool}|{sp}|{r}"] for r in reps_infig], axis=0)
        for i, pos in enumerate(range(-2, 3)):
            row = {"Species": sp, "Tool": tool, "Tool_display": DISPLAY.get(tool, tool),
                   "N_reps_in_figure": len(reps_infig),
                   "Replicates_in_figure": ";".join(reps_infig), "Pos": pos}
            row.update({b: float(mean[i, j]) for j, b in enumerate(BASES)})
            pwm_rows.append(row)
    pwm_mean = pd.DataFrame(pwm_rows)
    pwm_mean.to_csv(os.path.join(ANALYSIS, "figS2_pwm_mean.tsv"), sep="\t", index=False)
    log(f"figS2_pwm_per_rep.npz: {len(pwm_store)} replicate PWMs; "
        f"figS2_pwm_mean.tsv: {len(pwm_mean)} rows "
        f"({pwm_mean.groupby(['Species', 'Tool']).ngroups} species x tool groups)")

    # ---------------- E.coli coverage ----------------
    cov = pd.DataFrame(cov_rows).sort_values(["Condition", "Tool"])
    cov.to_csv(os.path.join(ANALYSIS, "figS2_ecoli_coverage.tsv"), sep="\t", index=False)
    if len(cov):
        log(f"figS2_ecoli_coverage.tsv: {len(cov)} rows "
            f"(conditions: {sorted(cov.Condition.unique())})")

    # ================= reconciliation 1: vs R3-6 fullcount_kl.tsv ==============
    log("\n" + "=" * 72)
    log("RECONCILIATION 1: per-rep KL_clip vs R3-6 fullcount_kl.tsv (rel tol 1e-9)")
    r36 = pd.read_csv(os.path.join(R36, "fullcount_kl.tsv"), sep="\t")
    r36["Species_Name"] = r36["Species_Name"].map(norm_species)  # "Human (HeLa)" -> "Human"
    ours = kl[["Tool", "Species", "Replicate", "KL_clip", "N_Sites"]].rename(
        columns={"Species": "Species_Name"})
    m = ours.merge(r36, on=["Tool", "Species_Name", "Replicate"],
                   how="inner", suffixes=("_ours", "_r36"))
    rel = (m["KL_clip"] - m["KL_Divergence"]).abs() / m["KL_Divergence"].abs().clip(lower=1e-300)
    n_pass = int((rel < 1e-9).sum())
    ns_mismatch = int((m["N_Sites_ours"] != m["N_Sites_r36"]).sum())
    only_ours = len(ours) - len(m)
    only_r36 = len(r36) - len(m)
    log(f"  matched rows      : {len(m)}")
    log(f"  PASS (rel < 1e-9) : {n_pass}/{len(m)}   max_rel_diff = {rel.max():.3e}")
    log(f"  N_Sites mismatches: {ns_mismatch}")
    log(f"  rows only in ours : {only_ours}   rows only in R3-6: {only_r36}")
    rec1_ok = (len(m) > 0) and (n_pass == len(m))
    log(f"  => RECONCILIATION 1 {'PASS' if rec1_ok else 'FAIL'}")

    # ================= reconciliation 2: vs legacy combined_kl_data.csv ========
    log("\nRECONCILIATION 2: pooled KL_clip vs legacy combined_kl_data.csv "
        f"({legacy_sub}, expect r ~ 0.991)")
    lg = legacy_kl[["Tool", "Species_norm", "KL_Divergence"]].copy()
    pl = pooled[["Species", "Tool", "KL_clip_mean_in_figure", "KL_clip_mean_all_reps"]].copy()
    pl["Tool"] = pl["Tool"].map(lambda t: DISPLAY.get(t, t))  # -> published names
    pl = pl.rename(columns={"Species": "Species_norm"})
    m2 = pl.merge(lg, on=["Tool", "Species_norm"], how="inner")
    m2 = m2[m2.Species_norm.isin(KL_SPECIES)]
    if len(m2) >= 3:
        from scipy import stats as sps
        # primary: in-figure pooled (Mouse = mES_WT only, published scope) -> ~0.991
        r_pear = float(np.corrcoef(m2.KL_clip_mean_in_figure, m2.KL_Divergence)[0, 1])
        r_all = float(np.corrcoef(m2.KL_clip_mean_all_reps, m2.KL_Divergence)[0, 1])
        rho = float(sps.spearmanr(m2.KL_clip_mean_in_figure, m2.KL_Divergence).statistic)
        log(f"  pairs={len(m2)}  Pearson r (in-figure pooled) = {r_pear:.4f}  "
            f"[all-reps pooled = {r_all:.4f}]  Spearman rho = {rho:.4f}")
        miss = lg[lg.Species_norm.isin(KL_SPECIES)].merge(
            pl, on=["Tool", "Species_norm"], how="left", indicator=True)
        miss = miss[miss._merge != "both"]
        if len(miss):
            log(f"  legacy rows without sites_clean match: "
                f"{miss[['Tool', 'Species_norm']].to_dict('records')}")
    else:
        r_pear = np.nan
        log("  WARNING: too few matched pairs")

    # ================= reconciliation 3: Top-5 sets ============================
    log("\nRECONCILIATION 3: Top-5 sets")
    r36_top = pd.read_csv(os.path.join(R36, "fullcount_top5.tsv"), sep="\t")
    r36_sets = r36_top.groupby(["Species", "Tool"])["5mer"].apply(set).to_dict()
    our_sets = top5.groupby(["Species", "Tool"])["Kmer"].apply(set).to_dict()
    common = sorted(set(r36_sets) & set(our_sets))
    ovl = [len(r36_sets[k] & our_sets[k]) for k in common]
    log(f"  vs R3-6 fullcount_top5 (pooled-count ranking): groups={len(common)}, "
        f"mean overlap={np.mean(ovl):.2f}/5, identical={sum(o == 5 for o in ovl)}")
    diffs = [k for k in common if r36_sets[k] != our_sets[k]]
    for k in diffs[:12]:
        log(f"    diff {k}: ours={sorted(our_sets[k])} r36={sorted(r36_sets[k])}")

    leg_top_path = os.path.join(LEGACY_ROOT, "Total_new1", "All_Species_Top5_5mer_Data.csv")
    if os.path.exists(leg_top_path):
        lt = pd.read_csv(leg_top_path)
        lt["Species_norm"] = lt["Species"].map(norm_species)
        lt["Tool_sc"] = lt["Tool"].map(lambda t: REVERSE_DISPLAY.get(t, t))
        lt_sets = lt.groupby(["Species_norm", "Tool_sc"])["5mer"].apply(set).to_dict()
        common2 = sorted(set(lt_sets) & set(our_sets))
        ovl2 = [len(lt_sets[k] & our_sets[k]) for k in common2]
        log(f"  vs legacy All_Species_Top5 (Top-20-truncated table): groups={len(common2)}, "
            f"mean overlap={np.mean(ovl2):.2f}/5 (differences expected: "
            f"full count vs Top-20 truncation)")

    # ================= reconciliation 4: vs Figure-4 revision products =========
    # Both figures must be drawn from the same computation: S2 panel A carries
    # the 10 tools that Figure 4B does not show.
    log("\n" + "=" * 72)
    log("RECONCILIATION 4: parity with the Figure-4 revision products")
    f4 = pd.read_csv(FIG4_PER_REP, sep="\t")
    ours4 = kl[["Tool", "Species", "Replicate", "KL", "KL_clip", "N_Sites"]].copy()
    ours4["Species_Name"] = ours4["Species"].map(FIG4_SPNAME)   # S2 label -> Fig.4 label
    n_dropped4 = int(ours4["Species_Name"].isna().sum())        # E.coli only
    ours4 = ours4.dropna(subset=["Species_Name"])
    m4 = ours4.merge(f4, on=["Tool", "Species_Name", "Replicate"], how="inner",
                     suffixes=("_ours", "_f4"))
    rec4_ok = len(m4) > 0
    if rec4_ok:
        rel_add = ((m4["KL"] - m4["KL_rrach_addeps"]).abs()
                   / m4["KL_rrach_addeps"].abs().clip(lower=1.0))
        rel_clip = ((m4["KL_clip"] - m4["KL_rrach_clip"]).abs()
                    / m4["KL_rrach_clip"].abs().clip(lower=1.0))
        ns_mismatch4 = int((m4["N_Sites_ours"] != m4["N_Sites_f4"]).sum())
        n_pass4 = int(((rel_add < 1e-9) & (rel_clip < 1e-9)).sum())
        log(f"  matched rows      : {len(m4)} (ours {len(ours4)}, fig4 {len(f4)}; "
            f"non-Fig.4 species rows skipped: {n_dropped4})")
        log(f"  PASS (rel < 1e-9) : {n_pass4}/{len(m4)}   max rel diff: "
            f"add-eps {rel_add.max():.3e}, clip {rel_clip.max():.3e}")
        log(f"  N_Sites mismatches: {ns_mismatch4}")
        rec4_ok = (n_pass4 == len(m4)) and (ns_mismatch4 == 0) \
            and (len(ours4) == len(f4))
    else:
        log("  FAIL: no overlapping rows with fig4_kl_per_rep.tsv")

    # per-replicate PWM equality (relabelled to the Fig.4 species names)
    f4p = np.load(FIG4_PWM_NPZ)
    ours_pwm = {}
    for key, arr in pwm_store.items():
        tool, sp, rep = key.split("|")
        if sp in FIG4_SPNAME:
            ours_pwm[f"{tool}|{FIG4_SPNAME[sp]}|{rep}"] = arr
    common_keys = sorted(set(f4p.files) & set(ours_pwm))
    only_f4 = sorted(set(f4p.files) - set(ours_pwm))
    only_ours = sorted(set(ours_pwm) - set(f4p.files))
    maxabs = max((float(np.abs(f4p[k] - ours_pwm[k]).max()) for k in common_keys),
                 default=float("nan"))
    log(f"  PWM keys          : ours {len(ours_pwm)}, fig4 {len(f4p.files)}, "
        f"common {len(common_keys)}")
    log(f"  PWM max |diff|    : {maxabs:.3e}   (only in fig4: {len(only_f4)}, "
        f"only in ours: {len(only_ours)})")
    if only_f4 or only_ours or not (maxabs < 1e-12) \
            or len(common_keys) != len(f4p.files):
        rec4_ok = False
    log(f"  => RECONCILIATION 4 {'PASS' if rec4_ok else 'FAIL'}")

    # ---------------- summary ----------------
    log("\n" + "=" * 72)
    log("SUMMARY")
    log(f"  rec1 (R3-6 exactness)    : {'PASS' if rec1_ok else 'FAIL'}")
    log(f"  rec2 (legacy r)          : r = {r_pear:.4f} (n={len(m2)})")
    log(f"  rec4 (Fig.4 rev parity)  : {'PASS' if rec4_ok else 'FAIL'}")
    log(f"  Panel A subset           : {int(pooled.in_panelA.sum())} tool x species cells")
    log(f"  Panel B groups           : {top5.groupby(['Species', 'Tool']).ngroups}")
    log("=" * 72)

    with open(LOG_PATH, "w") as fh:
        fh.write("\n".join(_log_lines) + "\n")
    log(f"log written: {LOG_PATH}")

    if not (rec1_ok and rec4_ok):
        sys.exit(1)


if __name__ == "__main__":
    main()
