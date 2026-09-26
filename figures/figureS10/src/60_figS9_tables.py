#!/usr/bin/env python3
"""60 -- Figure S9 tables: region shares and WT-vs-IVT contrast.

Supplementary Figure S9 shows the metagene distribution of the Dorado
*other-modification* models (pseU, m5C) on the RNA004 HeLa libraries.  The
published figure (``05_submission_work/02_AS_working_copy_and_revisions/sup/sup9.pdf``) contained
only curves; this script produces the numbers that carry the claim, so the
figure cannot be read as "m5C/pseU modification distribution" any more:

* every call comes from the analysis-ready layer ``sites_v2/sites_clean``
  (site-level de-duplication, own strand, delivered operating point
  ``percent_modified >= 90 %``);
* sites are placed on the **Ensembl** GRCh38p14 release-112 protein-coding
  transcript model (GENCODE is forbidden project-wide, guarded below);
* the unmodified IVT library is the negative control: a position preference
  that is reproduced on IVT is not modification-specific (reviewer R2-2/E1,
  R3-9, and the moderated wording demanded by R3-8/E8).

Outputs (``../tables/``):

* ``figS9_region_shares.tsv``    per model x condition segment counts/shares
* ``figS9_wt_ivt_contrast.tsv``  per model WT-minus-IVT differences/ratios
* ``figS9_key_numbers.md``       every number quoted in the legend/text

Usage
-----
conda run -n benchmark-revision --no-capture-output \
    python $RNAMODBENCH_ROOT/figures/figureS10/scripts/60_figS9_tables.py
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
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

PROJECT = Path(str(_RB))
SITES = (_RB / "data")
CLEAN = (_RB / "data/sites_clean")
OUT = (_RB / "figures/figureS10")
TABLES = (_RB / "figures/figureS10/tables")

sys.path.insert(0, str((_RB / "src/sites_v2")))
from common.match import fix_chromosome                                  # noqa: E402
from common.regionmodel import (CDS, KINDS, UTR3, UTR5, RegionIndex,     # noqa: E402
                                assign_sites)

#: the segments of the transcript body (the denominator of every ``share_*``)
BODY_KINDS = (UTR5, CDS, UTR3)

#: the six model panels of block A, in the order they are drawn (left to
#: right, top to bottom -> SLOTS; the page prints no per-panel letter)
MODELS = [
    ("Psi", "Dorado_hac@v5.0.0_pseU@v1_otherMod", "hac@v5.0.0_pseU"),
    ("Psi", "Dorado_hac@v5.1.0_pseU@v1_otherMod", "hac@v5.1.0_pseU"),
    ("Psi", "Dorado_sup@v5.0.0_pseU@v1_otherMod", "sup@v5.0.0_pseU"),
    ("Psi", "Dorado_sup@v5.1.0_pseU@v1_otherMod", "sup@v5.1.0_pseU"),
    ("m5C", "Dorado_hac@v5.1.0_m5C@v1_otherMod", "hac@v5.1.0_m5C"),
    ("m5C", "Dorado_sup@v5.1.0_m5C@v1_otherMod", "sup@v5.1.0_m5C"),
]
#: the page letters BLOCKS, not panels (user decision 2026-09-21): A = the six
#: Guitar density panels (one block), B = the Curlcake threshold scan (bottom
#: left), C = the HeLa score validity (bottom right).  Inside block A a panel is
#: addressed by its grid slot R1C1 .. R2C3 (left to right, top to bottom);
#: no per-panel letter is printed anywhere.
BLOCK_A = "A"
SLOTS = [f"R{r}C{c}" for r in (1, 2) for c in (1, 2, 3)]
assert len(SLOTS) == len(MODELS) == 6

#: models that exist in the data layer but are NOT drawn (reason in the legend,
#: reviewer R1-6 asks for explicit inclusion/exclusion reasons)
EXCLUDED = [
    ("inosine", "Dorado_hac@v5.1.0_inosine_m6A_otherMod", 15, 0),
    ("inosine", "Dorado_sup@v5.1.0_inosine_m6A_otherMod", 32, 6),
]
SAMPLES = {"WT": ("RNA004_HeLa_WT", "HeLa_RNA004_WT"),
           "IVT": ("RNA004_HeLa_IVT", "HeLa_RNA004_IVT")}
#: raw call rows of every input file (frozen 2026-09-20; any change upstream
#: must fail loudly instead of silently re-scaling the figure)
EXPECTED_RAW = {
    ("Psi", "Dorado_hac@v5.0.0_pseU@v1_otherMod", "WT"): 279,
    ("Psi", "Dorado_hac@v5.0.0_pseU@v1_otherMod", "IVT"): 138,
    ("Psi", "Dorado_hac@v5.1.0_pseU@v1_otherMod", "WT"): 235,
    ("Psi", "Dorado_hac@v5.1.0_pseU@v1_otherMod", "IVT"): 70,
    ("Psi", "Dorado_sup@v5.0.0_pseU@v1_otherMod", "WT"): 549,
    ("Psi", "Dorado_sup@v5.0.0_pseU@v1_otherMod", "IVT"): 953,
    ("Psi", "Dorado_sup@v5.1.0_pseU@v1_otherMod", "WT"): 610,
    ("Psi", "Dorado_sup@v5.1.0_pseU@v1_otherMod", "IVT"): 233,
    ("m5C", "Dorado_hac@v5.1.0_m5C@v1_otherMod", "WT"): 725,
    ("m5C", "Dorado_hac@v5.1.0_m5C@v1_otherMod", "IVT"): 261,
    ("m5C", "Dorado_sup@v5.1.0_m5C@v1_otherMod", "WT"): 218,
    ("m5C", "Dorado_sup@v5.1.0_m5C@v1_otherMod", "IVT"): 97,
}
#: assigned-site counts of the 2026-09-20 dry run (same convention as below)
EXPECTED_ASSIGNED = {
    ("Psi", "Dorado_hac@v5.0.0_pseU@v1_otherMod", "WT"): 94,
    ("Psi", "Dorado_hac@v5.0.0_pseU@v1_otherMod", "IVT"): 56,
    ("Psi", "Dorado_hac@v5.1.0_pseU@v1_otherMod", "WT"): 62,
    ("Psi", "Dorado_hac@v5.1.0_pseU@v1_otherMod", "IVT"): 27,
    ("Psi", "Dorado_sup@v5.0.0_pseU@v1_otherMod", "WT"): 221,
    ("Psi", "Dorado_sup@v5.0.0_pseU@v1_otherMod", "IVT"): 175,
    ("Psi", "Dorado_sup@v5.1.0_pseU@v1_otherMod", "WT"): 201,
    ("Psi", "Dorado_sup@v5.1.0_pseU@v1_otherMod", "IVT"): 96,
    ("m5C", "Dorado_hac@v5.1.0_m5C@v1_otherMod", "WT"): 396,
    ("m5C", "Dorado_hac@v5.1.0_m5C@v1_otherMod", "IVT"): 136,
    ("m5C", "Dorado_sup@v5.1.0_m5C@v1_otherMod", "WT"): 118,
    ("m5C", "Dorado_sup@v5.1.0_m5C@v1_otherMod", "IVT"): 52,
}
SEG = ["five_prime_UTR", "CDS", "three_prime_UTR"]
FLANKS = ["five_prime_flank", "three_prime_flank"]
COV_MIN = 10
SHARE_EQUIV_PP = 5.0     # |delta 3'UTR share| below this = "reproduced on IVT"
BOOTSTRAP_B = 1000       # site-level resamples per curve
RNG_SEED = 20260920      # same seed convention as the S6 unit-level bootstrap

# --------------------------------------------------------------------------- #
# Strict false-positive rate of the same six models on the unmodified control
# (reviewer R3-8/E8: "clearly defined FPR/specificity metrics").  The numbers
# are NOT recomputed here: they are moved from the frozen evaluation table used
# by Fig. 8B/D, with hard cross-checks so that any upstream drift fails loudly.
# --------------------------------------------------------------------------- #
FROZEN_FPR = (_RB / "data/evaluation/tables/controls_ivt_fpr.tsv")
FPR_SAMPLE = "HeLa_RNA004_IVT"

# --------------------------------------------------------------------------- #
# Block B (bottom left): threshold scan of the false-positive density on the
# unmodified Curlcake control (reviewer R3-8/E8: the FP metric in its most
# informative, threshold-resolved form).  Unlike the HeLa Dorado call sets,
# the Curlcake call sets are NOT pre-filtered at 90 % modified read, so the
# threshold can be swept.  Naming follows the Curlcake run of the same model
# families (the synthetic control has no plain pseU@v1 model).
# --------------------------------------------------------------------------- #
CURLC_GROUP = "RNA004_Curlcake_IVT"    # dataset-group directory name
CURLC_SAMPLE = "Curlcake_RNA004_IVT"   # sample name in the frozen table
CURLC_BP = 10135                      # 4 constructs, 10,135 bp
CURLC_UNIVERSE = {"Psi": 4928, "m5C": 3963}   # coverage >= 10, base-compatible
THRESHOLDS = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90]
CURLC_MODELS = [   # (mod, model, HeLa label, slot in block A, basecaller, version)
    ("Psi", "Dorado_hac@v5.0.0_pseU_m6A_Psi", "hac@v5.0.0_pseU", SLOTS[0], "hac", "v5.0.0"),
    ("Psi", "Dorado_hac@v5.1.0_all_Psi", "hac@v5.1.0_pseU", SLOTS[1], "hac", "v5.1.0"),
    ("Psi", "Dorado_sup@v5.0.0_pseU_m6A_Psi", "sup@v5.0.0_pseU", SLOTS[2], "sup", "v5.0.0"),
    ("Psi", "Dorado_sup@v5.1.0_all_Psi", "sup@v5.1.0_pseU", SLOTS[3], "sup", "v5.1.0"),
    ("m5C", "Dorado_hac@v5.1.0_all_m5C", "hac@v5.1.0_m5C", SLOTS[4], "hac", "v5.1.0"),
    ("m5C", "Dorado_sup@v5.1.0_all_m5C", "sup@v5.1.0_m5C", SLOTS[5], "sup", "v5.1.0"),
]
#: frozen anchors: (raw calls, calls >= 5 %, calls >= 50 %, calls >= 90 %)
EXPECTED_SCAN = {
    "Dorado_hac@v5.0.0_pseU_m6A_Psi": (2422, 89, 2, 0),
    "Dorado_hac@v5.1.0_all_Psi": (2400, 107, 5, 0),
    "Dorado_sup@v5.0.0_pseU_m6A_Psi": (2406, 70, 3, 0),
    "Dorado_sup@v5.1.0_all_Psi": (2368, 111, 7, 0),
    "Dorado_hac@v5.1.0_all_m5C": (2488, 246, 4, 0),
    "Dorado_sup@v5.1.0_all_m5C": (2435, 106, 2, 0),
}

# --------------------------------------------------------------------------- #
# Block C (bottom right): does the model's own reported modified fraction separate
# wild-type from unmodified RNA?  (reviewer R3-9: signal-level reasons for the
# non-m6A failures).  Scores are the delivered percent_modified values, i.e.
# all >= 0.90; the question is whether the tail above 0.99 differs.
# --------------------------------------------------------------------------- #
SCORE_SATURATION = 0.99
#: frozen AUC anchors (WT vs IVT, Mann-Whitney)
EXPECTED_AUC = {
    "Dorado_hac@v5.0.0_pseU@v1_otherMod": 0.4767,
    "Dorado_hac@v5.1.0_pseU@v1_otherMod": 0.5799,
    "Dorado_sup@v5.0.0_pseU@v1_otherMod": 0.2394,
    "Dorado_sup@v5.1.0_pseU@v1_otherMod": 0.6074,
    "Dorado_hac@v5.1.0_m5C@v1_otherMod": 0.5486,
    "Dorado_sup@v5.1.0_m5C@v1_otherMod": 0.5036,
}
EXPECTED_N_UNIVERSE = {"Psi": 15817652,      # base-compatible, coverage >= 10
                       "m5C": 14155581}
EXPECTED_REGION_BP = 157105601               # HeLa annotated exons
FPR_TOL = 1e-5                               # relative; the frozen table is %.6g
#: (mod, tool) -> frozen (fp_per_10kb, fp_per_1e6_candidates)
EXPECTED_FPR = {
    ("Psi", "Dorado_hac@v5.0.0_pseU@v1_otherMod"): (0.00878390, 8.724430),
    ("Psi", "Dorado_hac@v5.1.0_pseU@v1_otherMod"): (0.00445560, 4.425440),
    ("Psi", "Dorado_sup@v5.0.0_pseU@v1_otherMod"): (0.06065980, 60.249100),
    ("Psi", "Dorado_sup@v5.1.0_pseU@v1_otherMod"): (0.01483080, 14.730400),
    ("m5C", "Dorado_hac@v5.1.0_m5C@v1_otherMod"): (0.01661300, 18.438000),
    ("m5C", "Dorado_sup@v5.1.0_m5C@v1_otherMod"): (0.00617419, 6.852420),
}


def clean_path(mod: str, tool: str, condition: str) -> Path:
    group, sample = SAMPLES[condition]
    return (_RB / "data/sites_clean/RNA004/Human") / group / mod / tool / f"{sample}.tsv"


def load_sites(path: Path) -> pd.DataFrame:
    """One call per genomic position (highest score kept), own strand."""
    d = pd.read_csv(path, sep="\t",
                    usecols=["chrom", "start", "strand", "score", "coverage"],
                    low_memory=False)
    d = d.assign(pos=pd.to_numeric(d["start"], errors="coerce"))
    d = d.dropna(subset=["pos"]).copy()
    d["pos"] = d["pos"].astype("int64")
    d["score"] = pd.to_numeric(d["score"], errors="coerce")
    d["coverage"] = pd.to_numeric(d["coverage"], errors="coerce")
    d["strand"] = np.where(d["strand"].astype(str).isin(["-", "-1"]), "-", "+")
    d["chrom"] = d["chrom"].astype(str).map(fix_chromosome)
    d = d.sort_values("score", ascending=False, kind="stable")
    return d.drop_duplicates(["chrom", "pos"]).reset_index(drop=True)


def assign_kinds(model: RegionIndex, sites: pd.DataFrame) -> np.ndarray:
    """Segment kind per site (``-1`` = outside the mRNA region model)."""
    kind = np.full(len(sites), -1, dtype=np.int8)
    for strand, sub in sites.groupby("strand", sort=False):
        for chrom, ss in sub.groupby("chrom", sort=False):
            si, ki, _nz, _tl = assign_sites(
                model, chrom, strand, ss["pos"].to_numpy(np.int64))
            if len(si) == 0:
                continue
            kind[ss.index.to_numpy()[si]] = ki
    return kind


def describe(n_assigned: int, kind: np.ndarray, idx: np.ndarray) -> dict:
    """Counts and shares of one subset of the assigned sites."""
    sub = kind[idx]
    n_body = int(np.isin(sub, [0, 1, 2]).sum())
    out = {"n_assigned": int(idx.size), "n_body": n_body}
    for k, name in KINDS.items():
        n = int((sub == k).sum())
        out[f"n_{name}"] = n
        out[f"share_{name}"] = (n / n_body) if name in SEG and n_body else np.nan
        out[f"assigned_share_{name}"] = (n / idx.size) if idx.size else np.nan
    return out


def bootstrap_ci(masks: dict[str, np.ndarray], rng: np.random.Generator,
                 b: int = BOOTSTRAP_B) -> dict[str, tuple[float, float, np.ndarray]]:
    """Percentile bootstrap of the transcript-body segment shares.

    ``masks`` maps a segment name to a boolean membership vector over the
    library's transcript-body sites (the ``share_*`` denominator).  ``b``
    resample sets of the same size are drawn **with replacement** from those
    sites, one index matrix is shared by all segments, and the 2.5/97.5
    percentiles give the 95 % interval.  The resampling unit is a called site
    and not a sequencing unit: each condition is one library (n = 1), so the
    interval quantifies call sampling, never replicate variability.
    """
    segs = list(masks)
    n = int(masks[segs[0]].size) if segs else 0
    if n == 0:
        return {s: (np.nan, np.nan, np.zeros(b)) for s in segs}
    idx = rng.integers(0, n, size=(b, n))
    out: dict[str, tuple[float, float, np.ndarray]] = {}
    for s in segs:
        vals = masks[s][idx].mean(axis=1)
        lo, hi = np.percentile(vals, [2.5, 97.5])
        out[s] = (float(lo), float(hi), vals)
    return out


def ivt_fpr_table(shares: pd.DataFrame) -> pd.DataFrame:
    """Strict FPR of the six models on the unmodified HeLa IVT control.

    FPR = calls on a genuinely unmodified library, normalised as in Fig. 8B/D:
    per 10 kb of mappable sequence and per 10^6 base-compatible candidate sites
    (coverage >= 10; Psi: A/T, m5C: C/G - hence two different denominators).
    Every number is cross-checked against the frozen ``controls_ivt_fpr.tsv``
    row and against this figure's own IVT call count.
    """
    frozen = pd.read_csv(FROZEN_FPR, sep="\t")
    frozen = frozen[frozen["sample"] == FPR_SAMPLE]
    rows = []
    for i, (mod, tool, label) in enumerate(MODELS, start=1):
        sub = frozen[(frozen["tool"] == tool) & (frozen["mod_type"] == mod)]
        assert len(sub) == 1, f"expected one frozen FPR row for {tool}, got {len(sub)}"
        r = sub.iloc[0]
        n_ivt = int(shares[(shares.model == tool)
                           & (shares.condition == "IVT")]["n_calls_raw"].iloc[0])
        assert int(r["n_calls"]) == n_ivt, \
            f"IVT call count drift: {label} {int(r['n_calls'])} != {n_ivt}"
        assert int(r["n_universe"]) == EXPECTED_N_UNIVERSE[mod], \
            f"candidate-universe drift: {label} {int(r['n_universe'])}"
        assert int(r["region_bp"]) == EXPECTED_REGION_BP, \
            f"mappable-length drift: {label} {int(r['region_bp'])}"
        want = EXPECTED_FPR[(mod, tool)]
        got = (float(r["fp_per_10kb"]), float(r["fp_per_1e6_candidates"]))
        assert all(abs(g - w) <= FPR_TOL * abs(w) for g, w in zip(got, want)), \
            f"frozen FPR drift: {label} {got} != {want}"
        # the two columns must be the definitions of the two normalisations
        assert abs(got[0] - 1e4 * n_ivt / EXPECTED_REGION_BP) <= FPR_TOL * got[0]
        assert abs(got[1] - 1e6 * n_ivt / EXPECTED_N_UNIVERSE[mod]) <= FPR_TOL * got[1]
        rows.append({
            "block": BLOCK_A, "slot": SLOTS[i - 1], "mod": mod, "model": tool,
            "label": label,
            "n_calls": n_ivt, "n_calls_in_universe": int(r["n_calls_in_universe"]),
            "n_universe": int(r["n_universe"]), "region_bp": int(r["region_bp"]),
            "fp_per_10kb": got[0], "fp_per_1e6_candidates": got[1],
            "min_cov": int(r["min_cov"]), "source": str(FROZEN_FPR.relative_to(PROJECT)),
        })
    return pd.DataFrame(rows)


def curlcake_scan() -> pd.DataFrame:
    """False-positive density on the unmodified Curlcake control vs threshold.

    Every call on an unmodified construct is a false positive; the density is
    reported per 10 kb of the 10,135-bp construct sequence and per 10^6
    base-compatible candidate sites (coverage >= 10).  The frozen evaluation
    table pins the raw call number and the two denominators; the threshold
    profile must be non-increasing and hits the frozen 5 % / 50 % / 90 % anchors.
    """
    frozen = pd.read_csv(FROZEN_FPR, sep="\t")
    frozen = frozen[frozen["sample"] == CURLC_SAMPLE]
    rows = []
    for mod, model, label, slot, bc, ver in CURLC_MODELS:
        path = ((_RB / "data/sites_clean/RNA004/Curlcake") / CURLC_GROUP / mod / model
                / "Curlcake_RNA004_IVT.tsv")
        d = pd.read_csv(path, sep="\t", usecols=["score", "coverage"],
                        low_memory=False)
        s = pd.to_numeric(d["score"], errors="coerce").dropna().to_numpy(float) * 100
        raw, n5, n50, n90 = EXPECTED_SCAN[model]
        assert len(d) == raw, f"Curlcake row drift: {model} {len(d)} != {raw}"
        f = frozen[frozen["tool"] == model]
        assert len(f) == 1, f"expected one frozen Curlcake row for {model}, got {len(f)}"
        assert int(f["n_calls"].iloc[0]) == raw, "frozen Curlcake call drift"
        assert int(f["n_universe"].iloc[0]) == CURLC_UNIVERSE[mod], \
            f"candidate-universe drift on Curlcake: {model}"
        assert int(f["region_bp"].iloc[0]) == CURLC_BP, "construct length drift"
        counts = {t: int((s >= t).sum()) for t in THRESHOLDS}
        assert (counts[5], counts[50], counts[90]) == (n5, n50, n90), (
            f"threshold anchor drift: {label} "
            f"{counts[5]}/{counts[50]}/{counts[90]} != {n5}/{n50}/{n90}")
        prev = None
        for t in THRESHOLDS:
            n = counts[t]
            if prev is not None:
                assert n <= prev, f"scan not monotone: {label} at {t} %"
            prev = n
            rows.append({
                "block": BLOCK_A, "slot": slot, "mod": mod, "model": model,
                "label": label,
                "basecaller": bc, "version": ver, "threshold_pct": t,
                "n_calls": n, "fp_per_10kb": 1e4 * n / CURLC_BP,
                "fp_per_1e6_candidates": 1e6 * n / CURLC_UNIVERSE[mod],
                "n_candidates": CURLC_UNIVERSE[mod],
                "source": str(path.relative_to(PROJECT)),
            })
    return pd.DataFrame(rows)


def score_validity() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-call scores (long table) and the per-model WT-vs-IVT AUC summary.

    The delivered calls all satisfy >= 90 % modified read, so the question is
    whether the model's own confidence separates the wild-type library from the
    unmodified control (AUC = 0.5 means no discrimination; < 0.4 means the
    unmodified library scores higher).
    """
    long_rows, summ = [], []
    for i, (mod, tool, label) in enumerate(MODELS, start=1):
        slot = SLOTS[i - 1]
        per: dict[str, pd.DataFrame] = {}
        for cond in ("WT", "IVT"):
            sites = load_sites(clean_path(mod, tool, cond))
            per[cond] = sites
            for score, cov in zip(sites["score"].to_numpy(float),
                                  sites["coverage"].to_numpy(float)):
                long_rows.append({"block": BLOCK_A, "slot": slot, "mod": mod,
                                  "model": tool, "label": label,
                                  "condition": cond,
                                  "score": score, "coverage": cov})
        w = per["WT"]["score"].to_numpy(float)
        v = per["IVT"]["score"].to_numpy(float)
        u = mannwhitneyu(w, v, alternative="two-sided")
        auc = float(u.statistic) / (len(w) * len(v))
        assert abs(auc - EXPECTED_AUC[tool]) < 5e-4, \
            f"AUC drift: {label} {auc:.4f} != {EXPECTED_AUC[tool]:.4f}"
        verdict = ("inverted (unmodified control scores higher)" if auc < 0.4 else
                   "no discrimination" if auc <= 0.6 else
                   "wild type scores higher")
        summ.append({
            "block": BLOCK_A, "slot": slot, "mod": mod, "model": tool,
            "label": label,
            "n_WT": len(w), "n_IVT": len(v),
            "median_score_WT": float(np.median(w)),
            "median_score_IVT": float(np.median(v)),
            "auc_WT_vs_IVT": round(auc, 4), "p_mannwhitney": float(u.pvalue),
            "frac_ge_0p99_WT": float((w >= SCORE_SATURATION).mean()),
            "frac_ge_0p99_IVT": float((v >= SCORE_SATURATION).mean()),
            "median_coverage_WT": float(per["WT"]["coverage"].median()),
            "median_coverage_IVT": float(per["IVT"]["coverage"].median()),
            "verdict": verdict,
        })
    return pd.DataFrame(long_rows), pd.DataFrame(summ)


def main() -> None:
    t0 = time.time()
    model_path = [p for p in sorted(((_XB / "reference/regionmodels")).glob(
        "Human.*.mrna.regionmodel.pkl"))
        if "gencode" not in p.name.lower()]
    assert len(model_path) == 1, f"expected exactly one Ensembl human mRNA model: {model_path}"
    print(f"[60] region model: {model_path[0].name}", flush=True)
    region = RegionIndex.load(model_path[0])

    rows = []
    boot_vals: dict[tuple[str, str], np.ndarray] = {}
    rng = np.random.default_rng(RNG_SEED)   # one stream, MODELS x (WT, IVT) order
    for mod, tool, label in MODELS:
        for condition in ("WT", "IVT"):
            path = clean_path(mod, tool, condition)
            raw = pd.read_csv(path, sep="\t", usecols=["chrom"], low_memory=False)
            n_raw = len(raw)
            assert n_raw == EXPECTED_RAW[(mod, tool, condition)], \
                f"raw row count drift: {path} {n_raw} != " \
                f"{EXPECTED_RAW[(mod, tool, condition)]}"
            sites = load_sites(path)
            kind = assign_kinds(region, sites)
            inside = np.flatnonzero(kind >= 0)
            assert inside.size == EXPECTED_ASSIGNED[(mod, tool, condition)], \
                (f"assigned-site drift: {label} {condition} {inside.size} != "
                 f"{EXPECTED_ASSIGNED[(mod, tool, condition)]}")
            d = describe(len(sites), kind, inside)
            # site-level bootstrap of the body shares (resampling unit = a site)
            body_kind = kind[inside][np.isin(kind[inside], BODY_KINDS)]
            assert body_kind.size == d["n_body"], "body-site count drift"
            ci = bootstrap_ci({KINDS[k]: (body_kind == k) for k in BODY_KINDS}, rng)
            for s, (lo, hi, _v) in ci.items():
                point = d[f"share_{s}"]
                assert np.isfinite(point) and lo - 1e-9 <= point <= hi + 1e-9, (
                    f"point estimate outside its bootstrap CI: {label}/{condition}/{s}")
            boot_vals[(tool, condition)] = ci["three_prime_UTR"][2]
            cov10 = np.flatnonzero((kind >= 0) & (sites["coverage"].to_numpy() >= COV_MIN))
            d10 = describe(len(sites), kind, cov10)
            print(f"[60] {label:<22} {condition:<3} raw={n_raw:>4} sites={len(sites):>4} "
                  f"assigned={d['n_assigned']:>4} 3'UTR={100*d['share_three_prime_UTR']:.1f}% "
                  f"[{100*ci['three_prime_UTR'][0]:.1f}, {100*ci['three_prime_UTR'][1]:.1f}] "
                  f"(cov>=10: {100*d10['share_three_prime_UTR']:.1f}%, n={d10['n_assigned']})",
                  flush=True)
            share_min = float(sites["score"].min())
            assert share_min >= 0.90, f"delivered threshold drift: {share_min}"
            rows.append({
                "mod": mod, "model": tool, "label": label, "condition": condition,
                "n_units": 1, "n_calls_raw": n_raw, "n_calls_site": len(sites),
                **{k: v for k, v in d.items()},
                **{f"share_{s}_ci_lo": ci[s][0] for s in SEG},
                **{f"share_{s}_ci_hi": ci[s][1] for s in SEG},
                "share_cov10_three_prime_UTR": d10["share_three_prime_UTR"],
                "n_cov10_body": d10["n_body"],
                "score_min": round(share_min, 4),
                "source": str(path.relative_to(PROJECT)),
            })
    shares = pd.DataFrame(rows)

    # ------------------------------------------------------------------ contrast
    body_sum = shares[[f"n_{s}" for s in SEG]].sum(axis=1)
    assert (body_sum == shares["n_body"]).all(), "body counts do not add up"
    assert np.allclose(shares[[f"share_{s}" for s in SEG]].sum(axis=1), 1.0,
                       atol=1e-9), "body shares do not sum to 1"

    con = []
    for mod, tool, label in MODELS:
        w = shares[(shares.model == tool) & (shares.condition == "WT")].iloc[0]
        v = shares[(shares.model == tool) & (shares.condition == "IVT")].iloc[0]
        rec = {"mod": mod, "model": tool, "label": label,
               "n_assigned_WT": int(w.n_assigned), "n_assigned_IVT": int(v.n_assigned)}
        for s in SEG:
            rec[f"share_{s}_WT"] = w[f"share_{s}"]
            rec[f"share_{s}_IVT"] = v[f"share_{s}"]
            rec[f"delta_pp_{s}"] = 100 * (w[f"share_{s}"] - v[f"share_{s}"])
            rec[f"ratio_WT_over_IVT_{s}"] = (w[f"share_{s}"] / v[f"share_{s}"]
                                             if v[f"share_{s}"] > 0 else np.nan)
        for ratio_name, num, den in (("three_prime_UTR_to_CDS", "three_prime_UTR", "CDS"),):
            rec[f"{ratio_name}_WT"] = (w[f"share_{num}"] / w[f"share_{den}"]
                                       if w[f"share_{den}"] > 0 else np.nan)
            rec[f"{ratio_name}_IVT"] = (v[f"share_{num}"] / v[f"share_{den}"]
                                        if v[f"share_{den}"] > 0 else np.nan)
        delta = rec["delta_pp_three_prime_UTR"]
        rec["three_prime_UTR_delta_pp"] = delta
        rec["verdict"] = ("reproduced on unmodified IVT"
                          if abs(delta) <= SHARE_EQUIV_PP else
                          ("higher in WT" if delta > 0 else "lower in WT"))
        # WT-minus-IVT difference from the two independent site-level bootstrap
        # distributions of the 3'UTR share (same seed stream, library-wise)
        d_vals = 100 * (boot_vals[(tool, "WT")] - boot_vals[(tool, "IVT")])
        d_lo, d_hi = np.percentile(d_vals, [2.5, 97.5])
        rec["delta_pp_ci_lo"] = float(d_lo)
        rec["delta_pp_ci_hi"] = float(d_hi)
        rec["delta_pp_ci_width"] = float(d_hi - d_lo)
        assert d_lo - 1e-9 <= delta <= d_hi + 1e-9, \
            f"delta outside its bootstrap CI: {label}"
        con.append(rec)
    contrast = pd.DataFrame(con)

    TABLES.mkdir(parents=True, exist_ok=True)
    shares.to_csv((_RB / "figures/figureS10/tables/figS9_region_shares.tsv"), sep="\t", index=False)
    contrast.to_csv((_RB / "figures/figureS10/tables/figS9_wt_ivt_contrast.tsv"), sep="\t", index=False)
    print(f"[60] wrote {TABLES/'figS9_region_shares.tsv'} and "
          f"{TABLES/'figS9_wt_ivt_contrast.tsv'}", flush=True)

    # --------------------------------- reported values: strict FPR (R3-8/E8)
    fpr = ivt_fpr_table(shares)
    fpr.to_csv((_RB / "figures/figureS10/tables/figS9_ivt_fpr.tsv"), sep="\t", index=False)
    print(f"[60] wrote {TABLES/'figS9_ivt_fpr.tsv'} (reported values, frozen source "
          f"cross-checked)", flush=True)
    for _, r in fpr.iterrows():
        print(f"[60]   {r.slot} {r.label:<22} FP/10kb={r.fp_per_10kb:.5f} "
              f"FP/1e6={r.fp_per_1e6_candidates:.2f} (IVT calls {r.n_calls}, "
              f"universe {r.n_universe})", flush=True)

    # ------------------------------ block B (bottom left): Curlcake threshold scan
    scan = curlcake_scan()
    scan.to_csv((_RB / "figures/figureS10/tables/figS9_curlcake_scan.tsv"), sep="\t", index=False)
    print(f"[60] wrote {TABLES/'figS9_curlcake_scan.tsv'} (block B, bottom left)", flush=True)
    for label in scan["label"].unique():
        sub = scan[scan["label"] == label]
        calls = dict(zip(sub["threshold_pct"].astype(int), sub["n_calls"].astype(int)))
        rate50 = float(sub[sub["threshold_pct"] == 50]["fp_per_10kb"].iloc[0])
        print(f"[60]   {label:<22} calls @5/50/90 % = "
              f"{calls[5]:>3}/{calls[50]}/{calls[90]}  FP/10kb @50 % = {rate50:.2f}",
              flush=True)

    # ---------------------------- block C (bottom right): score validity in HeLa
    validity, vsum = score_validity()
    validity.to_csv((_RB / "figures/figureS10/tables/figS9_score_validity.tsv"), sep="\t", index=False)
    vsum.to_csv((_RB / "figures/figureS10/tables/figS9_score_validity_summary.tsv"), sep="\t", index=False)
    print(f"[60] wrote {TABLES/'figS9_score_validity.tsv'} and "
          f"{TABLES/'figS9_score_validity_summary.tsv'} (block C, bottom right)", flush=True)
    for _, r in vsum.iterrows():
        print(f"[60]   {r.slot} {r.label:<22} AUC={r.auc_WT_vs_IVT:.3f} "
              f"med {r.median_score_WT:.3f}/{r.median_score_IVT:.3f} "
              f"cov {r.median_coverage_WT:.0f}/{r.median_coverage_IVT:.0f} "
              f"-- {r.verdict}", flush=True)

    # ------------------------------------------------------------- key numbers
    def fmt_share(df: pd.DataFrame, s: str) -> str:
        return " / ".join(f"{100*x:.1f}%" for x in df[f"share_{s}"])

    lines = []
    lines.append("# Figure S9 -- key numbers (single source of truth)\n")
    lines.append(f"Generated by `60_figS9_tables.py` on "
                 f"{pd.Timestamp.now():%Y-%m-%d %H:%M}. Region model: "
                 f"`{model_path[0].name}` (Ensembl GRCh38p14 release 112, "
                 "protein-coding transcripts; GENCODE forbidden).\n")
    lines.append("## Inputs and definitions\n")
    lines.append("- Data layer: `sites_v2/sites_clean` (one call per genomic "
                 "position, highest score kept; the sample's own strand).")
    lines.append("- Delivered operating point: every Dorado call set in "
                 f"`sites_clean` has `percent_modified >= 0.90` "
                 f"(min score observed = {shares.score_min.min():.2f}).")
    lines.append("- Sequencing units: RNA004 HeLa WT n = 1, unmodified IVT n = 1 "
                 "(one library each, no biological replication; R3-2/E6).")
    lines.append("- Sites outside the mRNA region model (intronic, intergenic, "
                 "unannotated contigs) are reported as `n_calls_site - n_assigned` "
                 "and are not part of any share.")
    lines.append("- `share_*` uses the 5'UTR + CDS + 3'UTR body as denominator; "
                 "the 1 kb flanks are tabulated separately.")
    lines.append("- **Model names** are written `<mode>@<basecaller version>_"
                 "<modification>` (for example `hac@v5.1.0_pseU`); the modification "
                 "model version is `@v1` in every case, so the six full Dorado model "
                 "ids are `hac@v5.0.0_pseU@v1`, `hac@v5.1.0_pseU@v1`, "
                 "`sup@v5.0.0_pseU@v1`, `sup@v5.1.0_pseU@v1`, `hac@v5.1.0_m5C@v1` "
                 "and `sup@v5.1.0_m5C@v1`.")
    lines.append("- **Figure letters**: **A** = the six Guitar density panels as ONE "
                 "block (slots R1C1-R2C3 in the order hac@v5.0.0_pseU, "
                 "hac@v5.1.0_pseU, sup@v5.0.0_pseU, sup@v5.1.0_pseU, "
                 "hac@v5.1.0_m5C, sup@v5.1.0_m5C), **B** = the Curkcake threshold "
                 "scan (bottom left), **C** = the HeLa score validity (bottom "
                 "right); no per-panel letter is printed.")
    lines.append("- Uncertainty: **site-level percentile bootstrap**, "
                 f"B = {BOOTSTRAP_B}, seed = {RNG_SEED}, resampling the "
                 "transcript-body sites of one library; with n = 1 sequencing "
                 "unit per condition these intervals describe call sampling "
                 "within a library, not replicate uncertainty.\n")
    lines.append("## Calls and shares (per model x condition)\n")
    lines.append("| model | Dorado model id | condition | raw rows | sites | on mRNA | "
                 "5'UTR | CDS | 3'UTR [95% CI] | 3'UTR/CDS | "
                 "3'UTR share, coverage >= 10 |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for _, r in shares.iterrows():
        lines.append(
            f"| {r.label} | {r.model} | {r.condition} | {r.n_calls_raw} | "
            f"{r.n_calls_site} | {r.n_assigned} | "
            f"{100*r.share_five_prime_UTR:.1f}% ({r.n_five_prime_UTR}) | "
            f"{100*r.share_CDS:.1f}% ({r.n_CDS}) | "
            f"{100*r.share_three_prime_UTR:.1f}% [{100*r.share_three_prime_UTR_ci_lo:.1f}, "
            f"{100*r.share_three_prime_UTR_ci_hi:.1f}] ({r.n_three_prime_UTR}) | "
            f"{r.share_three_prime_UTR/r.share_CDS:.2f} | "
            f"{100*r.share_cov10_three_prime_UTR:.1f}% (n={r.n_cov10_body}) |")
    lines.append("\n## WT vs unmodified IVT\n")
    lines.append("| model | 3'UTR WT | 3'UTR IVT | delta (pp) | delta 95% CI (pp) | "
                 "WT/IVT ratio | verdict (|delta| <= 5 pp) |")
    lines.append("|---|---|---|---|---|---|---|")
    for _, r in contrast.iterrows():
        lines.append(f"| {r.label} | {100*r.share_three_prime_UTR_WT:.1f}% | "
                     f"{100*r.share_three_prime_UTR_IVT:.1f}% | "
                     f"{r.three_prime_UTR_delta_pp:+.1f} | "
                     f"[{r.delta_pp_ci_lo:+.1f}, {r.delta_pp_ci_hi:+.1f}] | "
                     f"{r.ratio_WT_over_IVT_three_prime_UTR:.2f} | {r.verdict} |")
    lines.append("\n## Excluded from this figure (R1-6 inclusion/exclusion reasons)\n")
    lines.append("| modification | model | WT calls | IVT calls | reason |")
    lines.append("|---|---|---|---|---|")
    for mod, tool, n_wt, n_ivt in EXCLUDED:
        lines.append(f"| {mod} | {tool} | {n_wt} | {n_ivt} | too few calls for a "
                     "density curve; the inosine+m6A models are counted in "
                     "Fig. S8A/C |")
    lines.append("\n## Headline\n")
    lines.append(
        f"- The 3'UTR-dominated profile of the other-modification models is "
        f"reproduced on the **unmodified IVT negative control** in every model "
        f"(IVT 3'UTR share {100*contrast.share_three_prime_UTR_IVT.min():.1f}--"
        f"{100*contrast.share_three_prime_UTR_IVT.max():.1f}% versus WT "
        f"{100*contrast.share_three_prime_UTR_WT.min():.1f}--"
        f"{100*contrast.share_three_prime_UTR_WT.max():.1f}%).")
    lines.append(
        f"- The largest WT-minus-IVT difference in the 3'UTR share is "
        f"{contrast.three_prime_UTR_delta_pp.abs().max():.1f} pp "
        f"({contrast.loc[contrast.three_prime_UTR_delta_pp.abs().idxmax(), 'label']}); "
        f"in {(contrast.three_prime_UTR_delta_pp.abs() <= SHARE_EQUIV_PP).sum()} of "
        f"{len(contrast)} models it stays within +/-{SHARE_EQUIV_PP:.0f} pp, and in "
        f"two pseU models the IVT share is even higher than the WT share.")
    half = np.maximum(
        (shares.share_three_prime_UTR_ci_hi - shares.share_three_prime_UTR),
        (shares.share_three_prime_UTR - shares.share_three_prime_UTR_ci_lo)) * 100
    fmax = contrast.loc[contrast.three_prime_UTR_delta_pp.abs().idxmax()]
    lines.append(
        f"- Uncertainty is site-level, not replicate-level: the 3'UTR shares carry "
        f"95 % bootstrap intervals of at most +/-{half.max():.1f} pp "
        f"(B = {BOOTSTRAP_B}, seed = {RNG_SEED}, resampling the transcript-body "
        f"sites of one library). The largest WT-minus-IVT gap ({fmax['label']}, "
        f"{fmax.three_prime_UTR_delta_pp:+.1f} pp) has a difference interval of "
        f"[{fmax.delta_pp_ci_lo:+.1f}, {fmax.delta_pp_ci_hi:+.1f}] pp.")
    lines.append(
        f"- The same six models stay at {fpr.fp_per_10kb.min():.5g}--"
        f"{fpr.fp_per_10kb.max():.5g} false positives per 10 kb and "
        f"{fpr.fp_per_1e6_candidates.min():.1f}--"
        f"{fpr.fp_per_1e6_candidates.max():.1f} per 10^6 candidate sites on the "
        f"unmodified HeLa IVT control (reported values, definitions identical to "
        f"Fig. 8B/D).")
    lines.append("\n## Uncertainty (site-level bootstrap, B = "
                 f"{BOOTSTRAP_B}, seed = {RNG_SEED})\n")
    lines.append("| model | condition | 3'UTR share | 95% CI (%) | CI width (pp) | body sites |")
    lines.append("|---|---|---|---|---|---|")
    for _, r in shares.iterrows():
        lines.append(
            f"| {r.label} | {r.condition} | {100*r.share_three_prime_UTR:.1f}% | "
            f"[{100*r.share_three_prime_UTR_ci_lo:.1f}, "
            f"{100*r.share_three_prime_UTR_ci_hi:.1f}] | "
            f"{100*(r.share_three_prime_UTR_ci_hi - r.share_three_prime_UTR_ci_lo):.1f} | "
            f"{r.n_body} |")
    lines.append("\nNote: with one sequencing unit per condition the interval "
                 "describes call sampling inside a library and must not be read as "
                 "replicate uncertainty; it answers the R1-8/R1-9/E7 request for an "
                 "uncertainty statement on the reported shares.\n")
    lines.append("\n## Block B (bottom left) -- threshold scan on the unmodified "
                 "Curkcake control (R3-8/E8)\n")
    lines.append("| slot in A | model (Curkcake run) | basecaller | version | "
                 "calls @5 % | @50 % | @90 % | FP per 10 kb @5 % | @50 % |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for label in scan["label"].unique():
        sub = scan[scan["label"] == label]
        r = sub.iloc[0]
        calls = dict(zip(sub["threshold_pct"].astype(int), sub["n_calls"].astype(int)))
        rate = dict(zip(sub["threshold_pct"].astype(int),
                        sub["fp_per_10kb"].astype(float)))
        lines.append(f"| {r['slot']} | {r['model']} | {r['basecaller']} | "
                     f"{r['version']} | {calls[5]} | {calls[50]} | {calls[90]} | "
                     f"{rate[5]:.2f} | {rate[50]:.2f} |")
    lines.append(
        f"\nDefinitions: every call on the unmodified synthetic construct "
        f"({CURLC_BP:,} bp, 4 constructs) is a false positive; the density is "
        f"per 10 kb of construct sequence (1 call = {1e4/CURLC_BP:.3f} per 10 kb) "
        f"and per 10^6 base-compatible candidate sites (Psi "
        f"{CURLC_UNIVERSE['Psi']:,}, m5C {CURLC_UNIVERSE['m5C']:,}; coverage >= 10). "
        f"The Curkcake run names the same model families differently from HeLa "
        f"(v5.0.0: `pseU_m6A_Psi`; v5.1.0: `all_Psi` / `all_m5C`; no plain "
        f"`pseU@v1`), so the rows carry the HeLa model names and their slot in block "
        f"A. Values at "
        f"every threshold are in `{scan['source'].iloc[0].rsplit('/', 1)[0]}` "
        f"(`figS9_curlcake_scan.tsv`).\n")
    lines.append("\n## Block C (bottom right) -- score validity in HeLa (R3-9)\n")
    lines.append("| slot in A | model | WT/IVT calls | median score WT / IVT | "
                 ">= 0.99 WT / IVT | median coverage WT / IVT | AUC (WT vs IVT) | p | verdict |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for _, r in vsum.iterrows():
        lines.append(
            f"| {r['slot']} | {r['label']} | {r['n_WT']} / {r['n_IVT']} | "
            f"{r['median_score_WT']:.3f} / {r['median_score_IVT']:.3f} | "
            f"{100*r['frac_ge_0p99_WT']:.1f} % / {100*r['frac_ge_0p99_IVT']:.1f} % | "
            f"{r['median_coverage_WT']:.0f} / {r['median_coverage_IVT']:.0f} | "
            f"{r['auc_WT_vs_IVT']:.3f} | {r['p_mannwhitney']:.2g} | {r['verdict']} |")
    lines.append(
        "\nDefinitions: score = the model's own reported modified fraction "
        f"(delivered call sets already satisfy >= 90 %, so the axis spans "
        f"0.90-1.00); AUC = Mann-Whitney common-language effect size between the "
        f"wild-type and unmodified-IVT calls (0.5 = no discrimination, < 0.4 = the "
        f"unmodified control scores higher). Verdict rule: < 0.4 inverted, "
        f"0.4-0.6 no discrimination, > 0.6 wild type higher.\n")
    lines.append("\n## Reported values -- strict false-positive rate on the "
                 "unmodified HeLa IVT control (delivered operating point; "
                 "reported values, not drawn)\n")
    lines.append("| slot in A | model | IVT calls (in universe) | FP per 10 kb | "
                 "FP per 10^6 candidates |")
    lines.append("|---|---|---|---|---|")
    for _, r in fpr.iterrows():
        lines.append(f"| {r.slot} | {r.label} | {r.n_calls} "
                     f"({r.n_calls_in_universe}) | {r.fp_per_10kb:.5g} | "
                     f"{r.fp_per_1e6_candidates:.3f} |")
    lines.append(
        f"\nDefinitions (identical to Fig. 8B/D): false positives are the calls on the "
        f"unmodified HeLa IVT library (delivered at >= 90 % modified-read fraction, "
        f"coverage >= {int(fpr.min_cov.iloc[0])}); normalised per 10 kb of mappable "
        f"sequence ({int(fpr.region_bp.iloc[0]):,} bp of annotated exons) and per "
        f"10^6 base-compatible candidate sites (Psi "
        f"{EXPECTED_N_UNIVERSE['Psi']:,}, m5C {EXPECTED_N_UNIVERSE['m5C']:,}; positions "
        f"with coverage >= 10). Source: `{fpr.source.iloc[0]}` (frozen evaluation "
        f"table used by Fig. 8B/D, cross-checked row by row here at "
        f"{FPR_TOL:g} relative tolerance).\n")
    lines.append("\n## Model-level statements\n")
    for _, r in contrast.iterrows():
        lines.append(f"- **{r['label']}** ({r['mod']}): 3'UTR "
                     f"{100*r.share_three_prime_UTR_WT:.1f}% in WT vs "
                     f"{100*r.share_three_prime_UTR_IVT:.1f}% on unmodified IVT "
                     f"({r.three_prime_UTR_delta_pp:+.1f} pp) -> {r.verdict}.")
    ((_RB / "figures/figureS10/tables/figS9_key_numbers.md")).write_text("\n".join(lines) + "\n")
    print(f"[60] wrote {TABLES/'figS9_key_numbers.md'} "
          f"({time.time() - t0:.1f} s)", flush=True)


if __name__ == "__main__":
    main()
