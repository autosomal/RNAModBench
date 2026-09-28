#!/usr/bin/env python
"""41 -- Figure S4 revision tables (R1-3/R3-3 window + AUPRC, R3-7 purified sites).

Rebuilds every number behind the revised Supplementary Figure S4 from the clean
per-replicate callsets (``harmonisation/callsets/``), replacing the legacy
assembly (``code/next_postprocessing/Total_new1/{AB_ratio,PRAUC}.ipynb``) which
used one union sample per species, no candidate universe and the pre-unify GLORI
beds.

Panels fed
----------
A  purified-site counts per tool, per independent WT/deficient pair
   (Arabidopsis WT x fip37 KD n=3, Mouse two cross-studies n=1 each, HeLa
   WT x IVT n=3 -- the HeLa pairs are NEW here; the frozen pipeline table
   ``evaluation/tables/purified_sites.tsv`` only covers KO/KD).
B  wider-window PPV -- read-only reuse of ``m6a_localization_curve.tsv``
   (drawn by 42; nothing computed here).
C  precision on the purified subset vs the full WT-common call set, both
   against GLORI at the primary 2-bp window (DESCRIPTIVE ONLY, R3-7).
D  AUPRC vs matching window {0,1,2,5,10,20,50} bp: calls ranked by a
   direction-normalised score inside the sample's own candidate universe,
   TP = call within +/- w of a testable GLORI site, recall denominator =
   number of testable GLORI sites; one PR curve per independent sequencing
   unit, then mean +- SD per species group (mouse studies never merged).
E/F/G independent validation of the purified subsets (R3-7 (iii)):
   E  GLORI(2 bp) hit rate of purified vs shared vs def-only site groups;
   F  DRACH fraction vs the per-sample universe background;
   G  score (stoichiometry-semantic) quantiles of purified vs shared.

Outputs -> figures/figureS4/tables/
  figS4_purified_sites.tsv
  figS4_precision_on_purified.tsv
  figS4_validation_groups.tsv
  figS4_universe_drach.tsv
  figS4_auprc_per_unit.tsv
  figS4_auprc_summary.tsv
  figS4_reconciliation.tsv          (anchor checks vs the frozen pipeline)

Anchor checks (script exits 2 on any mismatch)
  * KO/KD pairs: n_purified must equal ``purified_sites.tsv`` row by row;
  * HeLa/AUPRC side: precision@2bp on the sample's own universe must equal
    ``m6a_glori_confusion.tsv`` (window == 2) within 1e-9.

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/figures/figureS4/src/41_figS4_tables.py
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
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(_RB / "src/harmonisation"))

from common.config import (GLORI, PRIMARY_WINDOW, SITES_ROOT, TABLE_DIR,  # noqa: E402
                           UNIVERSE_ROOT, WINDOWS)
from common.evaluation import (load_reference, load_universe,  # noqa: E402
                               reference_in_universe)
from common.io_utils import read_table, write_table  # noqa: E402
from common.manifest import Inventory, log_time, setup_logger  # noqa: E402
from common.match import fix_chromosome, window_hit_mask  # noqa: E402

OUT = (_RB / "figures/figureS4")
TAB = OUT / "tables"
LOG = OUT / "logs"
SITES_CLEAN = (_RB / "data/callsets")

PLATFORM = "RNA002"
MOD = "m6A"
MIN_COV = 10
#: tools with fewer in-universe calls get PR-AUC = NA (curve unstable)
AP_MIN_CALLS = 100

#: fixed tool order / display (same as 37_fig5cd_window_sweep.py)
TOOL_ORDER = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
              "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
              "NanoSPA_m6A", "xPore", "yanocomp"]

#: score_type values that must be inverted to -log10 (larger = more modified)
SCORE_P_TYPES = {"Pvalue", "Padj", "FDR"}

#: (species, WT sample, deficient sample) -- mouse pairs stay study-separated
PAIRS: list[tuple[str, str, str]] = [
    ("Arabidopsis", "Arabidopsis_WT_rep1", "Arabidopsis_fip37_rep1"),
    ("Arabidopsis", "Arabidopsis_WT_rep2", "Arabidopsis_fip37_rep2"),
    ("Arabidopsis", "Arabidopsis_WT_rep3", "Arabidopsis_fip37_rep3"),
    ("Mouse", "mESCs_Mettl3_WT", "mESCs_Mettl3_KO"),
    ("Mouse", "mES_WT", "mES_KO"),
    ("Human", "HeLa_WT1", "HeLa_IVT_rep1"),
    ("Human", "HeLa_WT2", "HeLa_IVT_rep2"),
    ("Human", "HeLa_WT3", "HeLa_IVT_rep3"),
]
PAIR_ID = {
    ("Arabidopsis", "Arabidopsis_WT_rep1", "Arabidopsis_fip37_rep1"): "Ath_rep1",
    ("Arabidopsis", "Arabidopsis_WT_rep2", "Arabidopsis_fip37_rep2"): "Ath_rep2",
    ("Arabidopsis", "Arabidopsis_WT_rep3", "Arabidopsis_fip37_rep3"): "Ath_rep3",
    ("Mouse", "mESCs_Mettl3_WT", "mESCs_Mettl3_KO"): "Mouse_studyA",
    ("Mouse", "mES_WT", "mES_KO"): "Mouse_studyB",
    ("Human", "HeLa_WT1", "HeLa_IVT_rep1"): "HeLa_rep1",
    ("Human", "HeLa_WT2", "HeLa_IVT_rep2"): "HeLa_rep2",
    ("Human", "HeLa_WT3", "HeLa_IVT_rep3"): "HeLa_rep3",
}

#: AUPRC drawing units = the WT reference units of Fig. 5 (37 script PANELS)
AUPRC_GROUPS: dict[str, tuple[str, list[str]]] = {
    "Arabidopsis": ("Arabidopsis",
                    ["Arabidopsis_WT_rep1", "Arabidopsis_WT_rep2", "Arabidopsis_WT_rep3"]),
    "Mouse_studyA": ("Mouse", ["mESCs_Mettl3_WT"]),
    "Mouse_studyB": ("Mouse", ["mES_WT"]),
    "HeLa": ("Human", ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"]),
}

CRITERION_KOKD = ("called in WT and absent from the KO/KD callset, both restricted "
                  "to the common testable universe; circularity caveat per R3-7")
CRITERION_IVT = ("called in WT and absent from the matched IVT negative-control "
                 "callset; IVT is an unmodified transcriptome (true negative "
                 "control), NOT a knockdown of the same RNA -- descriptive "
                 "filter only, circularity caveat per R3-7")

logger = setup_logger("41_figS4_tables", log_dir=LOG)

# --------------------------------------------------------------------------- #
# loading helpers
# --------------------------------------------------------------------------- #
_universe_cache: dict[str, dict] = {}
_drach_cache: dict[str, tuple[int, int]] = {}
_rows_cache: dict[tuple[str, str], pd.DataFrame] = {}
_glori_cache: dict[str, dict] = {}
_refu_cache: dict[tuple[str, int], dict] = {}


def _universe_path(sample_name: str) -> Path:
    plain = UNIVERSE_ROOT / PLATFORM / species_of(sample_name) / f"{sample_name}__universe.tsv"
    gz = plain.with_suffix(".tsv.gz")
    return plain if plain.exists() else gz


def species_of(sample_name: str) -> str:
    s = SAMPLE_BY_NAME.get(sample_name)
    if s is None:
        raise KeyError(sample_name)
    return s.species


def universe_of(sample_name: str) -> dict[str, np.ndarray]:
    """Candidate-universe positions (m6A bases, coverage >= 10), cached."""
    if sample_name not in _universe_cache:
        up = _universe_path(sample_name)
        if not up.exists():
            raise FileNotFoundError(up)
        _universe_cache[sample_name] = load_universe(up, MOD, MIN_COV)
    return _universe_cache[sample_name]


def universe_drach_background(sample_name: str) -> tuple[int, int]:
    """(n_drach, n_universe) inside the sample's testable universe.

    The universe file marks DRACH compatibility with the strand of the match
    (``+``/``-``) or ``.`` when the position is not a DRACH context.
    """
    if sample_name in _drach_cache:
        return _drach_cache[sample_name]
    df = pd.read_csv(_universe_path(sample_name), sep="\t",
                     usecols=["base", "coverage", "drach"],
                     dtype={"base": str, "drach": str})
    want = {"A", "T"}  # MOD_GENOME_BASES["m6A"]
    m = df["coverage"].astype(int) >= MIN_COV
    m &= df["base"].isin(want)
    sub = df[m]
    n = int(len(sub))
    k = int(sub["drach"].isin(["+", "-"]).sum())
    _drach_cache[sample_name] = (k, n)
    return _drach_cache[sample_name]


def glori_of(species: str) -> dict[str, np.ndarray]:
    if species not in _glori_cache:
        _glori_cache[species] = load_reference(GLORI[species])
    return _glori_cache[species]


def ref_in_universe(species: str, universe: dict[str, np.ndarray],
                    key: str) -> dict[str, np.ndarray]:
    """Testable GLORI sites for one universe (cached by string key)."""
    if key not in _refu_cache:
        _refu_cache[key] = reference_in_universe(glori_of(species), universe)
    return _refu_cache[key]


def _callset_path(sample_name: str, tool: str) -> Path:
    return (SITES_CLEAN / PLATFORM / species_of(sample_name)
            / SAMPLE_BY_NAME[sample_name].dataset_group / MOD / tool
            / f"{sample_name}.tsv")


def unit_rows(sample_name: str, tool: str) -> pd.DataFrame:
    """callsets rows of one (sample, tool): chrom, pos, score, drach.

    ``score`` is the direction-normalised ranking score (larger = more likely
    modified): -log10 for Pvalue/Padj/FDR tools, the mod_ratio column for
    CHEUI_m6A (its Prob is constant ~1), otherwise the raw score.
    """
    key = (sample_name, tool)
    if key in _rows_cache:
        return _rows_cache[key]
    path = _callset_path(sample_name, tool)
    if not path.exists():
        _rows_cache[key] = pd.DataFrame(
            {"chrom": pd.Series(dtype=object), "pos": pd.Series(dtype=np.int64),
             "score": pd.Series(dtype=float), "drach": pd.Series(dtype=object)})
        return _rows_cache[key]
    df = read_table(path)
    if df.empty or "start" not in df.columns:
        _rows_cache[key] = pd.DataFrame(
            {"chrom": pd.Series(dtype=object), "pos": pd.Series(dtype=np.int64),
             "score": pd.Series(dtype=float), "drach": pd.Series(dtype=object)})
        return _rows_cache[key]
    out = pd.DataFrame({
        "chrom": [fix_chromosome(c) for c in df["chrom"]],
        "pos": pd.to_numeric(df["start"], errors="coerce").astype("Int64"),
    })
    if tool == "CHEUI_m6A" and "mod_ratio" in df.columns:
        sc = pd.to_numeric(df["mod_ratio"], errors="coerce")
    else:
        sc = pd.to_numeric(df["score"], errors="coerce")
        if "score_type" in df.columns:
            pm = df["score_type"].isin(SCORE_P_TYPES)
            if pm.any():
                neg = -np.log10(sc[pm].clip(lower=1e-300, upper=1.0))
                sc = sc.copy()
                sc.loc[pm] = neg
    out["score"] = sc.astype(float)
    if "drach" in df.columns:
        out["drach"] = (df["drach"].map({"True": True, "False": False})
                        .astype("boolean").fillna(False))
    else:
        out["drach"] = pd.Series([False] * len(out), dtype="boolean")
    n0 = len(out)
    out = out.dropna(subset=["pos", "score"])
    # site-level semantics (same as 06's np.unique dedup): overlapping genes can
    # make a tool emit the same genomic position twice with different scores;
    # keep the strongest evidence (highest score) per position.
    out = (out.sort_values("score", ascending=False, kind="mergesort")
              .drop_duplicates(subset=["chrom", "pos"], keep="first"))
    if len(out) != n0:
        logger.warning("  %s/%s: collapsed %d duplicate positions "
                       "(kept best score)", sample_name, tool, n0 - len(out))
    out["pos"] = out["pos"].astype(np.int64)
    _rows_cache[key] = out
    return _rows_cache[key]


def restrict(rows: pd.DataFrame, universe: dict[str, np.ndarray]) -> pd.DataFrame:
    """Keep only rows whose (chrom, pos) is an exact universe member."""
    parts = []
    for chrom, sub in rows.groupby("chrom", sort=False):
        u = universe.get(chrom)
        if u is None or u.size == 0:
            continue
        pos = sub["pos"].to_numpy(np.int64)
        idx = np.clip(np.searchsorted(u, pos), 0, u.size - 1)
        inside = u[idx] == pos
        if inside.any():
            parts.append(sub[inside])
    if not parts:
        return rows.iloc[0:0]
    return pd.concat(parts, ignore_index=True)


def common_universe(a: str, b: str) -> dict[str, np.ndarray]:
    """Intersection of two samples' universes (position-level, per chrom)."""
    ua, ub = universe_of(a), universe_of(b)
    out: dict[str, np.ndarray] = {}
    for chrom in set(ua) & set(ub):
        inter = np.intersect1d(ua[chrom], ub[chrom])
        if inter.size:
            out[chrom] = inter
    return out


def hit_rate_w2(rows: pd.DataFrame, ref_u: dict[str, np.ndarray]) -> tuple[int, int]:
    """(#rows within +/-2 of a testable GLORI site, #rows)."""
    n_hit = n_tot = 0
    for chrom, sub in rows.groupby("chrom", sort=False):
        refs = ref_u.get(chrom)
        pos = sub["pos"].to_numpy(np.int64)
        n_tot += pos.size
        if refs is None or not refs.size or pos.size == 0:
            continue
        n_hit += int(window_hit_mask(pos, refs, PRIMARY_WINDOW).sum())
    return n_hit, n_tot


def _arrays(rows: pd.DataFrame) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """{chrom: (pos, score, drach_bool)} for vectorised set operations."""
    out: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chrom, sub in rows.groupby("chrom", sort=False):
        dr = sub["drach"].astype("boolean").fillna(False).to_numpy(dtype=bool)
        out[chrom] = (sub["pos"].to_numpy(np.int64),
                      sub["score"].to_numpy(float), dr)
    return out


# --------------------------------------------------------------------------- #
# average precision (sklearn-compatible, distinct-threshold groups)
# --------------------------------------------------------------------------- #
def average_precision(hit: np.ndarray, score: np.ndarray, n_ref: int) -> float:
    """AP of a ranked call list against ``n_ref`` reference positives.

    Calls are sorted by descending score; ties share one threshold, whose
    precision/recall are evaluated at the END of the tie group (sklearn
    ``average_precision_score`` semantics).  ``hit[i]`` = call i is a TP.
    """
    if n_ref <= 0 or hit.size == 0:
        return float("nan")
    order = np.argsort(-score, kind="mergesort")
    h = hit[order].astype(np.int64)
    s = score[order]
    cum_tp = np.cumsum(h)
    ranks = np.arange(1, h.size + 1, dtype=float)
    prec = cum_tp / ranks
    rec = cum_tp / float(n_ref)
    ends = np.flatnonzero(np.diff(s) != 0)
    idx = np.append(ends, s.size - 1)
    r = rec[idx]
    p = prec[idx]
    dr = np.diff(np.concatenate(([0.0], r)))
    return float(np.sum(np.clip(dr, 0.0, None) * p))


# --------------------------------------------------------------------------- #
# part 1 -- pairs: purified counts + precision + validation groups
# --------------------------------------------------------------------------- #
def run_pairs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows_counts, rows_prec, rows_valid = [], [], []
    for species, wt, dfn in PAIRS:
        pair_id = PAIR_ID[(species, wt, dfn)]
        def_type = SAMPLE_BY_NAME[dfn].condition_class
        study_wt = SAMPLE_BY_NAME[wt].study
        study_def = SAMPLE_BY_NAME[dfn].study
        uni = common_universe(wt, dfn)
        n_common = int(sum(v.size for v in uni.values()))
        ref_u = ref_in_universe(species, uni, f"pair|{wt}|{dfn}")
        logger.info("[%s] %s vs %s (%s): common universe=%d, testable GLORI=%d",
                    pair_id, wt, dfn, def_type, n_common,
                    int(sum(v.size for v in ref_u.values())))
        for tool in TOOL_ORDER:
            wt_rows = restrict(unit_rows(wt, tool), uni)
            df_rows = restrict(unit_rows(dfn, tool), uni)
            wt_a, df_a = _arrays(wt_rows), _arrays(df_rows)

            pur_a, sha_a, don_a = {}, {}, {}
            for chrom in sorted(set(wt_a) | set(df_a)):
                w = wt_a.get(chrom, (np.zeros(0, np.int64), np.zeros(0), np.zeros(0, bool)))
                d = df_a.get(chrom, (np.zeros(0, np.int64), np.zeros(0), np.zeros(0, bool)))
                wpos, wsc, wdr = w
                dpos, dsc, ddr = d
                if wpos.size and dpos.size:
                    in_def = np.isin(wpos, dpos)
                    in_wt = np.isin(dpos, wpos)
                else:
                    in_def = np.zeros(wpos.size, bool)
                    in_wt = np.zeros(dpos.size, bool)
                pur_a[chrom] = (wpos[~in_def], wsc[~in_def], wdr[~in_def])
                sha_a[chrom] = (wpos[in_def], wsc[in_def], wdr[in_def])
                don_a[chrom] = (dpos[~in_wt], dsc[~in_wt], ddr[~in_wt])

            n_wt = sum(v[0].size for v in wt_a.values())
            n_def = sum(v[0].size for v in df_a.values())
            n_pur = sum(v[0].size for v in pur_a.values())
            n_sha = sum(v[0].size for v in sha_a.values())
            n_don = sum(v[0].size for v in don_a.values())

            crit = CRITERION_IVT if def_type == "IVT" else CRITERION_KOKD
            rows_counts.append({
                "species": species, "pair_id": pair_id, "def_type": def_type,
                "study_wt": study_wt, "study_def": study_def,
                "same_study": int(study_wt == study_def),
                "sample_wt": wt, "sample_def": dfn, "tool": tool,
                "n_common_universe": n_common,
                "n_wt_common": n_wt, "n_def_common": n_def,
                "n_shared": n_sha, "n_purified": n_pur, "n_def_only": n_don,
                "criterion": crit,
            })

            # -- C: precision on purified vs full WT-common (descriptive) ----
            pur_hit, pur_n = hit_rate_w2(_from_arrays(pur_a), ref_u)
            wt_hit, wtc_n = hit_rate_w2(_from_arrays(wt_a), ref_u)
            rows_prec.append({
                "species": species, "pair_id": pair_id, "tool": tool,
                "n_purified": pur_n, "n_purified_hit_glori_w2": pur_hit,
                "precision_purified_w2": (pur_hit / pur_n) if pur_n else np.nan,
                "n_wt_common": wtc_n, "n_wt_hit_glori_w2": wt_hit,
                "precision_wt_common_w2": (wt_hit / wtc_n) if wtc_n else np.nan,
                "delta_precision": ((pur_hit / pur_n - wt_hit / wtc_n)
                                    if pur_n and wtc_n else np.nan),
            })

            # -- E/F/G: three site groups -----------------------------------
            for gname, ga in (("purified", pur_a), ("shared", sha_a),
                              ("def_only", don_a)):
                g_rows = _from_arrays(ga)
                g_hit, g_n = hit_rate_w2(g_rows, ref_u)
                n_drach = sum(int(v[2].sum()) for v in ga.values())
                sc = (np.concatenate([v[1] for v in ga.values()])
                      if ga and any(v[1].size for v in ga.values())
                      else np.zeros(0))
                rows_valid.append({
                    "species": species, "pair_id": pair_id, "tool": tool,
                    "group": gname, "n_sites": g_n,
                    "n_hit_glori_w2": g_hit,
                    "glori_hit_rate_w2": (g_hit / g_n) if g_n else np.nan,
                    "n_drach": n_drach,
                    "drach_rate": (n_drach / g_n) if g_n else np.nan,
                    "score_q25": (float(np.percentile(sc, 25)) if sc.size else np.nan),
                    "score_median": (float(np.median(sc)) if sc.size else np.nan),
                    "score_q75": (float(np.percentile(sc, 75)) if sc.size else np.nan),
                    "score_semantic": SCORE_SEMANTIC[tool],
                })
    cols = ["species", "pair_id", "def_type", "study_wt", "study_def",
            "same_study", "sample_wt", "sample_def", "tool",
            "n_common_universe", "n_wt_common", "n_def_common", "n_shared",
            "n_purified", "n_def_only", "criterion"]
    return (pd.DataFrame(rows_counts)[cols],
            pd.DataFrame(rows_prec),
            pd.DataFrame(rows_valid))


def _from_arrays(a: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]
                 ) -> pd.DataFrame:
    """Rebuild a chrom/pos frame (for hit_rate_w2) from per-chrom arrays."""
    parts = []
    for chrom, (pos, _score, _dr) in a.items():
        if pos.size:
            parts.append(pd.DataFrame({"chrom": chrom, "pos": pos}))
    if not parts:
        return pd.DataFrame({"chrom": pd.Series(dtype=object),
                             "pos": pd.Series(dtype=np.int64)})
    return pd.concat(parts, ignore_index=True)


# --------------------------------------------------------------------------- #
# part 2 -- AUPRC per independent unit x tool x window
# --------------------------------------------------------------------------- #
def run_auprc() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, per_unit_w2 = [], []
    for grp, (species, units) in AUPRC_GROUPS.items():
        for unit in units:
            uni = universe_of(unit)
            ref_u = ref_in_universe(species, uni, f"unit|{unit}")
            n_ref = int(sum(v.size for v in ref_u.values()))
            for tool in TOOL_ORDER:
                u_rows = restrict(unit_rows(unit, tool), uni)
                n_calls = int(len(u_rows))
                rec = {"species_group": grp, "species": species, "unit": unit,
                       "study": SAMPLE_BY_NAME[unit].study, "tool": tool,
                       "n_calls_in_universe": n_calls,
                       "n_reference_in_universe": n_ref}
                if n_calls and n_ref:
                    pos_all, sc_all = [], []
                    hits = {w: [] for w in WINDOWS}
                    for chrom, sub in u_rows.groupby("chrom", sort=False):
                        pos = sub["pos"].to_numpy(np.int64)
                        refs = ref_u.get(chrom)
                        pos_all.append(pos)
                        sc_all.append(sub["score"].to_numpy(float))
                        for w in WINDOWS:
                            if refs is None or not refs.size:
                                hits[w].append(np.zeros(pos.size, bool))
                            else:
                                hits[w].append(window_hit_mask(pos, refs, w))
                    pos_v = np.concatenate(pos_all)
                    sc_v = np.concatenate(sc_all)
                    for w in WINDOWS:
                        hit_v = np.concatenate(hits[w])
                        rec[f"hit_rate_w{w}"] = float(hit_v.mean()) if n_calls else np.nan
                        rec[f"pr_auc_w{w}"] = (average_precision(hit_v, sc_v, n_ref)
                                               if n_calls >= AP_MIN_CALLS else np.nan)
                    prec2 = rec["hit_rate_w2"]
                    rows.append(rec)
                    per_unit_w2.append({**rec, "precision_w2": prec2})
                    logger.info("[AUPRC] %s/%s n_calls=%d n_ref=%d AP: %s",
                                unit, tool, n_calls, n_ref,
                                ", ".join(f"w{w}={rec[f'pr_auc_w{w}']:.4f}"
                                          for w in WINDOWS
                                          if rec[f'pr_auc_w{w}'] == rec[f'pr_auc_w{w}']))
                else:
                    for w in WINDOWS:
                        rec[f"hit_rate_w{w}"] = np.nan
                        rec[f"pr_auc_w{w}"] = np.nan
                    rows.append(rec)
                    per_unit_w2.append({**rec, "precision_w2": np.nan})
    per_unit = pd.DataFrame(rows)
    cols = (["species_group", "species", "unit", "study", "tool",
             "n_calls_in_universe", "n_reference_in_universe"]
            + [c for w in WINDOWS for c in (f"hit_rate_w{w}", f"pr_auc_w{w}")])
    per_unit = per_unit[cols]
    # mean +- SD across the units of a species group (per tool x window)
    summ_rows = []
    for grp, (_species, _units) in AUPRC_GROUPS.items():
        sub = per_unit[per_unit["species_group"] == grp]
        for tool in TOOL_ORDER:
            t = sub[sub["tool"] == tool]
            for w in WINDOWS:
                v = pd.to_numeric(t[f"pr_auc_w{w}"], errors="coerce").dropna()
                summ_rows.append({
                    "species_group": grp, "tool": tool, "window": w,
                    "n_units": int(len(v)),
                    "pr_auc_mean": float(v.mean()) if len(v) else np.nan,
                    "pr_auc_sd": float(v.std(ddof=1)) if len(v) > 1 else np.nan,
                })
    return per_unit, pd.DataFrame(summ_rows)


# --------------------------------------------------------------------------- #
# part 3 -- anchor checks against the frozen pipeline tables
# --------------------------------------------------------------------------- #
def reconcile(df_counts: pd.DataFrame, per_unit: pd.DataFrame) -> pd.DataFrame:
    recs = []
    frozen_pur = read_table(TABLE_DIR / "purified_sites.tsv")
    if not frozen_pur.empty:
        fp = {(r["sample_wt"], r["sample_ko"], r["tool"]): int(r["n_purified_sites"])
              for _, r in frozen_pur.iterrows()}
        for _, r in df_counts.iterrows():
            k = (r["sample_wt"], r["sample_def"], r["tool"])
            if r["def_type"] == "IVT" or k not in fp:
                continue
            mine = int(r["n_purified"])
            recs.append({"check": "purified_counts_vs_frozen",
                         "key": "|".join(k),
                         "frozen": fp[k], "mine": mine,
                         "ok": bool(fp[k] == mine)})
    frozen_conf = read_table(TABLE_DIR / "m6a_glori_confusion.tsv")
    if not frozen_conf.empty and not per_unit.empty:
        fc = frozen_conf[frozen_conf["window"] == "2"]
        fcm: dict[tuple[str, str], float] = {}
        for _, r in fc.iterrows():
            try:
                fcm[(r["sample"], r["tool"])] = float(r["precision"])
            except (TypeError, ValueError):
                continue  # empty precision (zero-call rows)
        for _, r in per_unit.iterrows():
            k = (r["unit"], r["tool"])
            if k not in fcm:
                continue
            mine = r["hit_rate_w2"]
            if mine == mine:
                # the frozen table stores %.6g floats -> relative tolerance
                ok = abs(float(mine) - fcm[k]) <= 1e-5 * max(1.0, abs(fcm[k]))
                recs.append({"check": "precision_w2_vs_confusion",
                             "key": "|".join(k), "frozen": fcm[k],
                             "mine": float(mine), "ok": bool(ok)})
    return pd.DataFrame(recs, columns=["check", "key", "frozen", "mine", "ok"])


SAMPLE_BY_NAME = {}
SCORE_SEMANTIC = {}


def main() -> None:
    global SAMPLE_BY_NAME, SCORE_SEMANTIC
    from common.config import SAMPLES_BY_NAME  # local import keeps module import cheap
    SAMPLE_BY_NAME = SAMPLES_BY_NAME
    SCORE_SEMANTIC = {
        "CHEUI_m6A": "mod_ratio", "m6Anet": "probability", "NanoSPA_m6A": "probability",
        "Nanom6A": "mod_ratio", "MINES": "mod_ratio", "DENA": "mod_ratio",
        "EpiNano_Error": "delta_sum_err",
        "DRUMMER": "neg_log10_p", "Nanocompore": "neg_log10_p",
        "ELIGOS2_diff": "neg_log10_padj", "ELIGOS2_solo": "neg_log10_padj",
        "xPore": "neg_log10_fdr", "yanocomp": "neg_log10_fdr",
    }
    TAB.mkdir(parents=True, exist_ok=True)
    inv = Inventory("41_figS4_tables")

    with log_time(logger, "pairs (purified + precision + validation)"):
        df_counts, df_prec, df_valid = run_pairs()
    with log_time(logger, "AUPRC per unit"):
        per_unit, per_summ = run_auprc()

    write_table(df_counts, TAB / "figS4_purified_sites.tsv")
    write_table(df_prec, TAB / "figS4_precision_on_purified.tsv")
    write_table(df_valid, TAB / "figS4_validation_groups.tsv")
    dr = [{"sample": s, "species": species_of(s),
           "n_universe": bg[1], "n_drach": bg[0],
           "drach_rate": (bg[0] / bg[1]) if bg[1] else np.nan}
          for s in sorted({w for _, w, _ in PAIRS} | {d for _, _, d in PAIRS}
                          | {u for _sp, (_s, us) in AUPRC_GROUPS.items() for u in us})
          for bg in [universe_drach_background(s)]]
    write_table(pd.DataFrame(dr), TAB / "figS4_universe_drach.tsv")
    write_table(per_unit, TAB / "figS4_auprc_per_unit.tsv")
    write_table(per_summ, TAB / "figS4_auprc_summary.tsv")

    with log_time(logger, "reconciliation"):
        rec = reconcile(df_counts, per_unit)
    write_table(rec, TAB / "figS4_reconciliation.tsv")

    for t in (TAB / "figS4_purified_sites.tsv", TAB / "figS4_precision_on_purified.tsv",
              TAB / "figS4_validation_groups.tsv", TAB / "figS4_universe_drach.tsv",
              TAB / "figS4_auprc_per_unit.tsv", TAB / "figS4_auprc_summary.tsv",
              TAB / "figS4_reconciliation.tsv"):
        try:
            n = sum(1 for _ in open(t)) - 1
            inv.record(t, n_rows=max(n, 0))
        except FileNotFoundError:
            pass
    inv.flush()

    bad = rec[~rec["ok"]] if not rec.empty else rec
    if bad is not None and len(bad):
        logger.error("RECONCILIATION FAILED: %d mismatches, e.g.\n%s",
                     len(bad), bad.head(10).to_string(index=False))
        sys.exit(2)
    logger.info("all anchor checks passed (%s rows)",
                0 if rec is None else len(rec))
    logger.info("tables -> %s", TAB)


if __name__ == "__main__":
    main()
