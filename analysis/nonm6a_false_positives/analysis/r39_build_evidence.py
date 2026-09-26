#!/usr/bin/env python3
"""R3-9 -- rebuild the non-m6A false-positive evidence tables (2026-09-19).

Replaces the superseded union-level shell scripts (archived in
``_superseded_20260919/``).  Everything here is **per replicate**: unions are
never treated as samples; the Curlcake depth-matched subset library is reported
but excluded from the independent mean.

Outputs (into ``../evidence/``)
------------------------------
curlcake_ivt_fp_per_construct.tsv    FP calls on the unmodified synthetic RNA, per construct
hela_wt_ivt_per_replicate.tsv        per-replicate HeLa WT/IVT calls, in-universe calls, FP density
hela_wt_ivt_summary.tsv              per tool-mod: unions, Jaccard (raw + universe), k-of-n, T/WT
truth_precision_per_replicate.tsv    precision vs RMBase+DirectRMDB / orthogonal NGS + permutation null
score_distributions.tsv              per-sample score summary + WT-vs-IVT discrimination
score_histograms.tsv                 score histograms used by the figure (fixed 0-1, 50 bins)
ecoli_gse271571_prob.tsv             third-party GSE271571 (E. coli) WT/IVT/rlm probability stats

Conventions
-----------
* analysis layer = ``sites_v2/sites_clean`` (0-based BED, coordinate fixes applied);
* scores / stoichiometry are read from ``sites_v2/callsets`` (holds ``mod_ratio``; the
  two layers are row-level equivalent after the centre-base filter);
* candidate universe = ``sites_v2/universe/<platform>/<species>/<sample>__universe.tsv[.gz]``
  filtered to ``coverage >= --min-cov`` and to bases that can carry the modification
  (identical to ``common.evaluation.load_universe``, i.e. the denominator behind the
  published FPR tables);
* permutation null = chromosome-stratified uniform resampling of the calls inside the
  sample's own candidate universe (as many draws as there were in-universe calls on
  that chromosome), R = ``--perm``; expected = mean over permutations; enrichment =
  observed / expected; empirical p = (1 + #{perm >= observed}) / (R + 1);
* no temporary file is written outside this project directory (no /tmp).

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    04_revision_analysis/R3-9_nonm6a_fp_analysis/analysis/r39_build_evidence.py
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
import argparse
import logging
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve()
PKG = HERE.parents[1]                     # R3-9_nonm6a_fp_analysis/
PROJECT = HERE.parents[3]                 # $RNAMODBENCH_ROOT
SITES = PROJECT / "04_revision_analysis" / "sites_v2"
sys.path.insert(0, str(PROJECT / "src/sites_v2"))

from common.config import SAMPLES_BY_NAME, TABLE_DIR        # noqa: E402
from common.match import fix_chromosome                     # noqa: E402
from common.manifest import setup_logger                    # noqa: E402

SC = (_RB / "data/sites_clean")
CS = (_XB / "sites_v2/callsets")
UNI = (_XB / "sites_v2/universe")

HELA_GROUPS: dict[str, list[str]] = {
    "WT": ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"],
    "IVT": ["HeLa_IVT_rep1", "HeLa_IVT_rep2", "HeLa_IVT_rep3"],
}
CURLCAKE = ["Curlcake_IVT_rep1", "Curlcake_IVT_rep2_partial", "Curlcake_IVT_rep3"]
CURLCAKE_SUBSET = "Curlcake_IVT_rep2_partial"
TOOLMOD: list[tuple[str, str]] = [
    ("Nm", "NanoNm"),
    ("Psi", "NanoMUD_psi"),
    ("Psi", "NanoPsu"),
    ("Psi", "NanoSPA_psU"),
    ("m1Psi", "NanoMUD_m1psi"),
    ("m5C", "CHEUI_m5C"),
]
#: RMBase + DirectRMDB database reference (human, hg38); note the source file
#: lives in a directory named ``orca_annotation`` -- that is only a folder
#: name, ORCA itself is an unrelated modification-calling tool.
ORCA_MOD = {"Psi": "pseudoU", "m5C": "m5C", "Nm": "Nm"}
ORCA_PATH = Path(str(_XB / "nanopore/reference/orca_annotation/Answer_from_RMBase_and_DirectRMDB_human.csv"))
NGS_REF = {
    "m5C": ("UBS-seq_m5C_HeLa", PROJECT / "07_third_party" / "NGS" / "UBS-seq" /
            "UBS-seq_m5C_HeLa.txt"),
    "Nm": ("Nm-Mut-seq_Nm_HeLa", PROJECT / "07_third_party" / "NGS" / "Nm-Mut-seq" /
           "Nm-Mut-seq_Nm_HeLa.bed"),
}
GLORI_HELA = PROJECT / "07_third_party" / "NGS" / "GLORI" / "Hela_GLORI.bed"
ECOLI = PROJECT / "07_third_party" / "GEO" / "GSE271571_Ecoli_epitranscriptome"

#: score column that actually carries information for each tool (CHEUI's Prob is
#: constant ~1; its stoichiometry lives in ``mod_ratio``)
SCORE_KIND = {"CHEUI_m5C": "mod_ratio"}
PROB_BINS = np.linspace(0.0, 1.0, 51)
#: genome-orientation bases that can carry each modification (mirrors
#: ``common.evaluation.MOD_GENOME_BASES``)
MOD_BASES = {"Nm": set("ACGT"), "m5C": {"C", "G"}, "Psi": {"T", "A"},
             "m6A": {"A", "T"}}


# --------------------------------------------------------------------------- #
# paths
# --------------------------------------------------------------------------- #
def _layer_path(root: Path, sample: str, mod: str, tool: str) -> Path:
    s = SAMPLES_BY_NAME[sample]
    return root / s.platform / s.species / s.dataset_group / mod / tool / f"{sample}.tsv"


def universe_path(sample: str) -> Path:
    s = SAMPLES_BY_NAME[sample]
    plain = UNI / s.platform / s.species / f"{s.canonical}__universe.tsv"
    return plain if plain.exists() else plain.with_suffix(".tsv.gz")


# --------------------------------------------------------------------------- #
# geometry helpers (0-based single-nucleotide positions, sorted int64 arrays)
# --------------------------------------------------------------------------- #
def _to_pos_dict(chroms, positions) -> dict[str, np.ndarray]:
    df = pd.DataFrame({"chrom": [fix_chromosome(c) for c in chroms],
                       "pos": pd.to_numeric(pd.Series(positions)).to_numpy(np.int64)})
    return {c: np.unique(sub["pos"].to_numpy(np.int64))
            for c, sub in df.groupby("chrom", sort=False)}


def load_positions(path: Path, chrom_col: str = "chrom", pos_col: str = "start"
                   ) -> dict[str, np.ndarray]:
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if df.empty:
        return {}
    return _to_pos_dict(df[chrom_col], df[pos_col])


def load_orca(mod_label: str) -> dict[str, np.ndarray]:
    df = pd.read_csv(ORCA_PATH, header=None, usecols=[0, 1, 3],
                     names=["chrom", "pos", "mod"], dtype={0: str})
    sub = df[df["mod"] == mod_label]
    return _to_pos_dict(sub["chrom"], sub["pos"])


def load_bed(path: Path) -> dict[str, np.ndarray]:
    df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1],
                     names=["chrom", "start"], dtype={0: str})
    return _to_pos_dict(df["chrom"], df["start"])


def count_within(call_pos: np.ndarray, ref: np.ndarray, w: int) -> int:
    """Number of call positions within +-w bp of any reference position."""
    if call_pos.size == 0 or ref.size == 0:
        return 0
    j = np.searchsorted(ref, call_pos, side="left")
    hit = np.zeros(call_pos.size, dtype=bool)
    m = j > 0
    if m.any():
        hit[m] = ref[j[m] - 1] >= call_pos[m] - w
    k = j < ref.size
    if k.any():
        hit[k] |= ref[j[k]] <= call_pos[k] + w
    return int(np.count_nonzero(hit))


def inside_universe(calls: dict[str, np.ndarray], uni: dict[str, np.ndarray]
                    ) -> tuple[dict[str, np.ndarray], int]:
    out: dict[str, np.ndarray] = {}
    n_out = 0
    for chrom, pos in calls.items():
        u = uni.get(chrom)
        if u is None or u.size == 0:
            n_out += int(pos.size)
            continue
        j = np.clip(np.searchsorted(u, pos), 0, u.size - 1)
        ok = u[j] == pos
        n_out += int((~ok).sum())
        if ok.any():
            out[chrom] = pos[ok]
    return out, n_out


def n_positions(pos_dict: dict[str, np.ndarray]) -> int:
    return int(sum(v.size for v in pos_dict.values()))


def global_jaccard(sets: list[dict[str, np.ndarray]]) -> float:
    arrs = [s for s in sets if n_positions(s)]
    if len(arrs) < 2:
        return float("nan")
    keys = set().union(*[set(a) for a in arrs])
    inter = union = 0
    for chrom in keys:
        parts = [a[chrom] for a in arrs if chrom in a]
        u = np.unique(np.concatenate(parts))
        union += u.size
        if len(parts) == len(arrs):
            common = parts[0]
            for other in parts[1:]:
                common = common[np.isin(common, other)]
            inter += common.size
    return inter / union if union else float("nan")


def mean_pairwise_jaccard(sets: list[dict[str, np.ndarray]]) -> float:
    vals = []
    for a, b in combinations(sets, 2):
        keys = set(a) | set(b)
        inter = union = 0
        for chrom in keys:
            pa = a.get(chrom, np.empty(0, np.int64))
            pb = b.get(chrom, np.empty(0, np.int64))
            union += np.union1d(pa, pb).size
            if pa.size and pb.size:
                inter += int(np.isin(pa, pb, assume_unique=True).sum())
        if union:
            vals.append(inter / union)
    return float(np.mean(vals)) if vals else float("nan")


def k_of_n(sets: list[dict[str, np.ndarray]]) -> dict[int, int]:
    """{k: number of union sites detected in exactly k replicates}."""
    if not sets:
        return {}
    keys = set().union(*[set(a) for a in sets])
    counts: dict[int, int] = {k: 0 for k in range(1, len(sets) + 1)}
    for chrom in keys:
        parts = [a[chrom] for a in sets if chrom in a]
        u = np.unique(np.concatenate(parts))
        kk = np.zeros(u.size, dtype=int)
        for p in parts:
            kk += np.isin(u, p, assume_unique=True).astype(int)
        for k, c in zip(*np.unique(kk, return_counts=True)):
            counts[int(k)] += int(c)
    return counts


def permutation_null(call_pos: dict[str, np.ndarray], uni: dict[str, np.ndarray],
                     ref: dict[str, np.ndarray], w: int, n_perm: int,
                     rng: np.random.Generator, chunk: int = 100) -> np.ndarray:
    """Chromosome-stratified resampling of the calls inside the candidate universe."""
    tot = np.zeros(n_perm, dtype=np.int64)
    for chrom, pos in call_pos.items():
        u, r = uni.get(chrom), ref.get(chrom)
        if r is None or u is None or u.size == 0 or r.size == 0 or pos.size == 0:
            continue
        n_c = int(pos.size)
        for s0 in range(0, n_perm, chunk):
            k = min(chunk, n_perm - s0)
            sim = u[rng.integers(0, u.size, size=(k, n_c))].ravel()
            j = np.searchsorted(r, sim, side="left")
            hit = np.zeros(sim.size, dtype=bool)
            m = j > 0
            if m.any():
                hit[m] = r[j[m] - 1] >= sim[m] - w
            kk = j < r.size
            if kk.any():
                hit[kk] |= r[j[kk]] <= sim[kk] + w
            tot[s0:s0 + k] += hit.reshape(k, n_c).sum(axis=1)
    return tot


# --------------------------------------------------------------------------- #
# cached inputs
# --------------------------------------------------------------------------- #
class Resources:
    def __init__(self, min_cov: int) -> None:
        self.min_cov = min_cov
        self._uni: dict[tuple[str, str], dict[str, np.ndarray]] = {}
        self._sites: dict[tuple[str, str, str], dict[str, np.ndarray]] = {}
        self._scores: dict[tuple[str, str, str], pd.DataFrame | None] = {}

    def _ensure_sample(self, sample: str, log: logging.Logger | None = None) -> None:
        """Read one candidate-universe file once and derive every mod subset.

        The per-sample universe file backs all modifications (``coverage >=
        min_cov`` + the reference-base column); deriving the base subsets in one
        pass avoids re-reading a ~40 MB gzip six times per sample and yields the
        same arrays as ``common.evaluation.load_universe``.
        """
        if any((sample, m) in self._uni for m in ("Nm", "m5C", "Psi", "m6A")):
            return
        t0 = time.time()
        df = pd.read_csv(universe_path(sample), sep="\t",
                         usecols=["chrom", "pos", "base", "coverage"],
                         dtype={"chrom": str, "base": str})
        df = df[df["coverage"] >= self.min_cov]
        df["chrom"] = [fix_chromosome(c) for c in df["chrom"]]
        for mod, bases in MOD_BASES.items():
            sub = df[df["base"].isin(bases)]
            self._uni[(sample, mod)] = {
                c: np.unique(s["pos"].to_numpy(np.int64))
                for c, s in sub.groupby("chrom", sort=False)}
        # Psi and m1Psi share the base constraint: reuse the same arrays
        self._uni[(sample, "m1Psi")] = self._uni[(sample, "Psi")]
        if log is not None:
            log.info("    universe %-22s Nm=%-10d m5C=%-9d Psi=%-9d m6A=%-9d (%.1fs)",
                     sample, n_positions(self._uni[(sample, "Nm")]),
                     n_positions(self._uni[(sample, "m5C")]),
                     n_positions(self._uni[(sample, "Psi")]),
                     n_positions(self._uni[(sample, "m6A")]), time.time() - t0)
        del df

    def universe(self, sample: str, mod: str,
                 log: logging.Logger | None = None) -> dict[str, np.ndarray]:
        self._ensure_sample(sample, log)
        return self._uni.get((sample, mod), {})

    def calls(self, sample: str, mod: str, tool: str) -> dict[str, np.ndarray]:
        key = (sample, mod, tool)
        if key not in self._sites:
            self._sites[key] = load_positions(_layer_path(SC, sample, mod, tool))
        return self._sites[key]

    def scored(self, sample: str, mod: str, tool: str) -> pd.DataFrame | None:
        """Callsets rows restricted to the sample's candidate universe."""
        key = (sample, mod, tool)
        if key in self._scores:
            return self._scores[key]
        path = _layer_path(CS, sample, mod, tool)
        if not path.exists():
            self._scores[key] = None
            return None
        head = pd.read_csv(path, sep="\t", nrows=0)
        wanted = ["chrom", "pos_raw", "score", "score_type", "mod_ratio", "coverage"]
        cols = [c for c in wanted if c in head.columns]
        df = pd.read_csv(path, sep="\t", low_memory=False, usecols=cols)
        for c in wanted:
            if c not in df.columns:
                df[c] = np.nan
        if df.empty:
            self._scores[key] = None
            return None
        df["chrom"] = [fix_chromosome(c) for c in df["chrom"]]
        keep = []
        uni = self.universe(sample, mod)
        for chrom, sub in df.groupby("chrom", sort=False):
            u = None if uni is None else uni.get(chrom)
            if u is None or u.size == 0:
                continue
            pos = sub["pos_raw"].to_numpy(np.int64)
            j = np.clip(np.searchsorted(u, pos), 0, u.size - 1)
            keep.append(sub[u[j] == pos])
        df = pd.concat(keep) if keep else df.iloc[0:0]
        self._scores[key] = df
        return df


# --------------------------------------------------------------------------- #
# layer 1 -- unmodified Curlcake controls
# --------------------------------------------------------------------------- #
def build_curlcake(res: Resources, region_bp: int, log: logging.Logger) -> pd.DataFrame:
    rows = []
    for mod, tool in TOOLMOD:
        for sample in CURLCAKE:
            if not _layer_path(SC, sample, mod, tool).exists():
                continue
            calls = res.calls(sample, mod, tool)
            uni = res.universe(sample, mod, log)
            uni_n = n_positions(uni)
            inside, n_out = inside_universe(calls, uni)
            n_in = n_positions(inside)
            rows.append(dict(
                mod_type=mod, tool=tool, construct=sample,
                construct_role=("depth_matched_subset_of_rep3"
                                if sample == CURLCAKE_SUBSET else "independent"),
                n_calls=n_positions(calls), n_calls_in_universe=n_in,
                n_calls_out_of_universe=n_out, universe_n=uni_n, region_bp=region_bp,
                fp_per_10kb=1e4 * n_in / region_bp if region_bp else np.nan,
                fp_per_1e6_candidates=1e6 * n_in / uni_n if uni_n else np.nan,
            ))
    df = pd.DataFrame(rows)
    return df.sort_values(["mod_type", "tool", "construct"]).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# layer 2 -- HeLa WT vs IVT, per replicate
# --------------------------------------------------------------------------- #
def build_hela(res: Resources, region_bp: int, log: logging.Logger
               ) -> tuple[pd.DataFrame, pd.DataFrame]:
    per_rep, summary = [], []
    for mod, tool in TOOLMOD:
        if not _layer_path(SC, HELA_GROUPS["WT"][0], mod, tool).exists():
            log.warning("  skip (no HeLa callset): %s / %s", mod, tool)
            continue
        log.info("  %s / %s", mod, tool)
        grp_sets: dict[str, list[dict[str, np.ndarray]]] = {}
        grp_sets_uni: dict[str, list[dict[str, np.ndarray]]] = {}
        counts: dict[str, list[int]] = {}
        union_len: dict[str, int] = {}
        for grp, samples in HELA_GROUPS.items():
            raw, unis, cnt = [], [], []
            for sample in samples:
                calls = res.calls(sample, mod, tool)
                uni = res.universe(sample, mod, log)
                uni_n = n_positions(uni)
                inside, n_out = inside_universe(calls, uni)
                n_in = n_positions(inside)
                raw.append(calls)
                unis.append(inside)
                cnt.append(n_in)
                per_rep.append(dict(
                    mod_type=mod, tool=tool, group=grp, sample=sample,
                    n_calls=n_positions(calls), n_calls_in_universe=n_in,
                    n_calls_out_of_universe=n_out, universe_n=uni_n,
                    region_bp=region_bp,
                    fp_per_10kb=1e4 * n_in / region_bp if region_bp else np.nan,
                    fp_per_1e6_candidates=1e6 * n_in / uni_n if uni_n else np.nan,
                ))
            grp_sets[grp], grp_sets_uni[grp], counts[grp] = raw, unis, cnt

        def _union(sets):
            out: dict[str, np.ndarray] = {}
            keys = set().union(*[set(a) for a in sets]) if sets else set()
            for chrom in keys:
                parts = [a[chrom] for a in sets if chrom in a]
                out[chrom] = np.unique(np.concatenate(parts))
            return out

        rec = dict(mod_type=mod, tool=tool)
        for grp in HELA_GROUPS:
            sets_raw, sets_uni = grp_sets[grp], grp_sets_uni[grp]
            kk = k_of_n(sets_uni)
            rec.update({
                f"{grp.lower()}_reps": len(sets_raw),
                f"{grp.lower()}_union_raw": n_positions(_union(sets_raw)),
                f"{grp.lower()}_union_in_universe": n_positions(_union(sets_uni)),
                f"{grp.lower()}_count_mean": float(np.mean(counts[grp])),
                f"{grp.lower()}_count_sd": float(np.std(counts[grp], ddof=1))
                if len(counts[grp]) > 1 else np.nan,
                f"{grp.lower()}_global_jaccard_raw": global_jaccard(sets_raw),
                f"{grp.lower()}_mean_pairwise_jaccard_raw": mean_pairwise_jaccard(sets_raw),
                f"{grp.lower()}_global_jaccard_uni": global_jaccard(sets_uni),
                f"{grp.lower()}_mean_pairwise_jaccard_uni": mean_pairwise_jaccard(sets_uni),
                f"{grp.lower()}_k1": kk.get(1, 0), f"{grp.lower()}_k2": kk.get(2, 0),
                f"{grp.lower()}_k3": kk.get(3, 0),
            })
            if sum(kk.values()):
                rec[f"{grp.lower()}_frac_ge2"] = (kk.get(2, 0) + kk.get(3, 0)) / sum(kk.values())
            else:
                rec[f"{grp.lower()}_frac_ge2"] = np.nan
        wt, ivt = counts["WT"], counts["IVT"]
        pairs = [i / w for i in ivt for w in wt if w]
        rec["ratio_twt_union"] = (rec["ivt_union_raw"] / rec["wt_union_raw"]
                                  if rec["wt_union_raw"] else np.nan)
        rec["ratio_twt_mean_counts"] = (float(np.mean(ivt)) / float(np.mean(wt))
                                        if np.mean(wt) else np.nan)
        rec["ratio_twt_pair_min"] = min(pairs) if pairs else np.nan
        rec["ratio_twt_pair_max"] = max(pairs) if pairs else np.nan
        summary.append(rec)
    return (pd.DataFrame(per_rep).sort_values(["mod_type", "tool", "group", "sample"])
            .reset_index(drop=True),
            pd.DataFrame(summary))


# --------------------------------------------------------------------------- #
# layer 3 -- truth-anchored precision, per replicate
# --------------------------------------------------------------------------- #
def build_truth(res: Resources, refs: dict[str, dict[str, dict[str, np.ndarray]]],
                n_perm: int, rng: np.random.Generator, log: logging.Logger) -> pd.DataFrame:
    rows = []
    samples = HELA_GROUPS["WT"] + HELA_GROUPS["IVT"]
    for sample in samples:
        grp = "WT" if sample in HELA_GROUPS["WT"] else "IVT"
        for mod, tool in TOOLMOD:
            if not _layer_path(SC, sample, mod, tool).exists():
                continue
            uni = res.universe(sample, mod, log)
            uni_n = n_positions(uni)
            called, _ = inside_universe(res.calls(sample, mod, tool), uni)
            n_in = n_positions(called)
            for ref_name, ref_by_mod in refs.items():
                ref = ref_by_mod.get(mod)
                if ref is None:
                    continue
                ref_in, _ = inside_universe(ref, uni)
                n_ref_in = n_positions(ref_in)
                obs0 = sum(count_within(p, ref_in.get(c, np.empty(0, np.int64)), 0)
                           for c, p in called.items())
                obs1 = sum(count_within(p, ref_in.get(c, np.empty(0, np.int64)), 1)
                           for c, p in called.items())
                perm0 = permutation_null(called, uni, ref_in, 0, n_perm, rng)
                perm1 = permutation_null(called, uni, ref_in, 1, n_perm, rng)
                for w, obs, perm in ((0, obs0, perm0), (1, obs1, perm1)):
                    exp = float(perm.mean()) if perm.size else np.nan
                    pval = ((1 + int(np.count_nonzero(perm >= obs))) / (perm.size + 1)
                            if perm.size else np.nan)
                    if exp == exp and exp > 0:
                        enrich = obs / exp
                    elif exp == 0:
                        enrich = 0.0 if obs == 0 else np.nan
                    else:
                        enrich = np.nan
                    rows.append(dict(
                        cohort="HeLa", sample=sample, group=grp, mod_type=mod, tool=tool,
                        reference=ref_name, window_bp=w,
                        n_calls=n_positions(res.calls(sample, mod, tool)),
                        n_calls_in_universe=n_in, universe_n=uni_n,
                        n_reference_sites_in_universe=n_ref_in,
                        n_overlap=int(obs),
                        precision=(obs / n_in) if n_in else np.nan,
                        chance_precision=(exp / n_in) if (n_in and exp == exp) else np.nan,
                        expected_overlap=exp,
                        enrichment=enrich,
                        p_empirical=pval, n_perm=n_perm,
                    ))
        log.info("  truth layer done: %s", sample)
    # GLORI sanity row: m6A pipeline check on the same HeLa libraries
    for sample in samples:
        grp = "WT" if sample in HELA_GROUPS["WT"] else "IVT"
        tool = "CHEUI_m6A"
        if not _layer_path(SC, sample, "m6A", tool).exists():
            continue
        uni = res.universe(sample, "m6A", log)
        called, _ = inside_universe(res.calls(sample, "m6A", tool), uni)
        ref_in, _ = inside_universe(refs["GLORI"]["m6A"], uni)
        obs = sum(count_within(p, ref_in.get(c, np.empty(0, np.int64)), 1)
                  for c, p in called.items())
        n_in = n_positions(called)
        rows.append(dict(
            cohort="validation", sample=sample, group=grp, mod_type="m6A", tool=tool,
            reference="GLORI_m6A", window_bp=1,
            n_calls=n_positions(res.calls(sample, "m6A", tool)),
            n_calls_in_universe=n_in, universe_n=n_positions(uni),
            n_reference_sites_in_universe=n_positions(ref_in), n_overlap=int(obs),
            precision=(obs / n_in) if n_in else np.nan,
            chance_precision=np.nan, expected_overlap=np.nan, enrichment=np.nan,
            p_empirical=np.nan, n_perm=0,
        ))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# layer 4 -- score distributions
# --------------------------------------------------------------------------- #
def build_scores(res: Resources, log: logging.Logger
                 ) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries, hist_rows = [], []
    samples = (HELA_GROUPS["WT"] + HELA_GROUPS["IVT"]
               + [s for s in CURLCAKE if s != CURLCAKE_SUBSET])
    for mod, tool in TOOLMOD:
        kind = SCORE_KIND.get(tool, "score")
        col_vals: dict[tuple[str, str], np.ndarray] = {}
        by_sample: dict[tuple[str, str], np.ndarray] = {}
        for sample in samples:
            if not _layer_path(CS, sample, mod, tool).exists():
                continue
            res.universe(sample, mod, log)          # populate the candidate universe first
            df = res.scored(sample, mod, tool)
            if df is None or df.empty:
                continue
            vals = pd.to_numeric(df[kind], errors="coerce").dropna().to_numpy(float)
            if vals.size == 0:
                continue
            grp_short = ("WT" if sample in HELA_GROUPS["WT"] else
                         ("IVT" if sample in HELA_GROUPS["IVT"] else "IVT_Curlcake"))
            col_vals[(sample, grp_short)] = vals
            by_sample[(sample, grp_short)] = vals
            summaries.append(dict(
                mod_type=mod, tool=tool, score_kind=kind,
                sample=sample,
                group=("WT" if sample in HELA_GROUPS["WT"] else
                       ("IVT" if sample in HELA_GROUPS["IVT"] else "Curlcake_IVT")),
                score_type=str(df["score_type"].iloc[0]), n_sites=int(vals.size),
                mean=float(vals.mean()), sd=float(vals.std(ddof=1)) if vals.size > 1 else np.nan,
                median=float(np.median(vals)), q05=float(np.quantile(vals, .05)),
                q25=float(np.quantile(vals, .25)), q75=float(np.quantile(vals, .75)),
                q95=float(np.quantile(vals, .95)),
                min=float(vals.min()), max=float(vals.max()),
                frac_ge_0p9=float(np.mean(vals >= 0.9)),
            ))
        # histogram edges: quantised scores (Psi tools, all ~0.95-1.0) need a zoomed
        # range, otherwise every call lands in the last 0-1 bin
        hela_pool = [v for (s, g), v in by_sample.items()
                     if g in ("WT", "IVT") and v.size]
        if hela_pool:
            pooled = np.concatenate(hela_pool)
            if np.nanmin(pooled) >= 0.9:
                edges = np.linspace(0.9, 1.0, 51)
            elif np.nanmin(pooled) >= 0.0:
                edges = np.linspace(0.0, 1.0, 51)
            else:
                lo, hi = float(np.nanmin(pooled)), float(np.nanmax(pooled))
                edges = np.linspace(lo, hi, 51)
        else:
            edges = PROB_BINS
        for (sample, grp), vals in by_sample.items():
            counts, _ = np.histogram(vals, bins=edges)
            for i, c in enumerate(counts):
                hist_rows.append(dict(
                    row_type="histogram", mod_type=mod, tool=tool, score_kind=kind,
                    sample=sample, group=grp,
                    bin_lo=float(edges[i]), bin_hi=float(edges[i + 1]),
                    count=int(c),
                    density=float(c / (vals.size * (edges[1] - edges[0]))),
                ))
        # WT vs IVT discrimination (HeLa only, pooled + per replicate pair)
        wt = [(s, v) for (s, g), v in col_vals.items() if g == "WT"]
        ivt = [(s, v) for (s, g), v in col_vals.items() if g == "IVT"]
        if wt and ivt:
            a = np.concatenate([v for _, v in wt])
            b = np.concatenate([v for _, v in ivt])
            u = stats.mannwhitneyu(a, b, alternative="two-sided")
            ks = stats.ks_2samp(a, b)
            pair_aucs = []
            for _, va in wt:
                for _, vb in ivt:
                    uu = stats.mannwhitneyu(va, vb, alternative="two-sided").statistic
                    pair_aucs.append(uu / (va.size * vb.size))
            summaries.append(dict(
                mod_type=mod, tool=tool, score_kind=kind, sample="__discrimination__",
                group="WT_vs_IVT", score_type="", n_sites=int(a.size + b.size),
                mean=np.nan, sd=np.nan, median=np.nan, q05=np.nan, q25=np.nan,
                q75=np.nan, q95=np.nan, min=np.nan, max=np.nan, frac_ge_0p9=np.nan,
                auc=float(u.statistic / (a.size * b.size)), p_mannwhitney=float(u.pvalue),
                ks_stat=float(ks.statistic), ks_p=float(ks.pvalue),
                pair_auc_mean=float(np.mean(pair_aucs)),
                pair_auc_min=float(np.min(pair_aucs)), pair_auc_max=float(np.max(pair_aucs)),
            ))
        log.info("  score layer done: %s / %s (%s)", mod, tool, kind)
    return pd.DataFrame(summaries), pd.DataFrame(hist_rows)


# --------------------------------------------------------------------------- #
# layer 5 -- third-party GSE271571 (E. coli)
# --------------------------------------------------------------------------- #
def build_ecoli(log: logging.Logger) -> pd.DataFrame:
    rows = []
    sl = (_XB / "third_party/GEO/GSE271571_Ecoli_epitranscriptome/CHEUI/site_level")
    for cond in ("WT", "IVT", "rlm"):
        for mod in ("m5C", "m6A"):
            path = sl / f"GSE271571_{cond}.cheui.{mod}.site.level.predictions.txt"
            if not path.exists():
                continue
            df = pd.read_csv(path, sep="\t", usecols=["coverage", "probability"])
            cov = pd.to_numeric(df["coverage"], errors="coerce")
            for tag, thr in (("all", 0), ("cov_ge_20", 20), ("cov_ge_50", 50)):
                sub = df if thr == 0 else df.loc[cov >= thr]
                p = pd.to_numeric(sub["probability"], errors="coerce").dropna().to_numpy(float)
                if p.size == 0:
                    continue
                rows.append(dict(
                    row_type="site_summary", dataset="GSE271571", species="E.coli",
                    condition=cond, mod_type=mod, coverage_filter=tag,
                    n_sites=int(p.size), n_prob_ge_0p9=int(np.count_nonzero(p >= 0.9)),
                    frac_prob_ge_0p9=float(np.mean(p >= 0.9)),
                    median_prob=float(np.median(p)), mean_prob=float(p.mean()),
                    q90_prob=float(np.quantile(p, .90)),
                ))
                if tag == "cov_ge_20":
                    counts, _ = np.histogram(p, bins=PROB_BINS)
                    for i, c in enumerate(counts):
                        rows.append(dict(
                            row_type="prob_bin", dataset="GSE271571", species="E.coli",
                            condition=cond, mod_type=mod, coverage_filter=tag,
                            bin_lo=PROB_BINS[i], bin_hi=PROB_BINS[i + 1], count=int(c),
                            density=float(c / (p.size * (PROB_BINS[1] - PROB_BINS[0]))),
                        ))
    diff_dir = (_XB / "third_party/GEO/GSE271571_Ecoli_epitranscriptome/CHEUI/differential")
    for comp in ("WT_vs_IVT", "WT_vs_rlm"):
        for mod in ("m5C", "m6A"):
            path = diff_dir / f"GSE271571_{comp}_cheui_differential_{mod}_sites.txt"
            if not path.exists():
                continue
            df = pd.read_csv(path, sep="\t")
            pv = pd.to_numeric(df.get("pval_U"), errors="coerce")
            d = pd.to_numeric(df.get("stoichiometry_diff"), errors="coerce")
            rows.append(dict(
                row_type="differential", dataset="GSE271571", species="E.coli",
                condition=comp, mod_type=mod, coverage_filter="",
                n_sites=len(df), n_pval_lt_0p05=int(np.count_nonzero(pv < 0.05)),
                n_diff_pos=int(np.count_nonzero(d > 0)),
                n_diff_neg=int(np.count_nonzero(d < 0)),
            ))
    log.info("  E.coli layer done (%d rows)", len(rows))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
def region_bp_from_tables(log: logging.Logger) -> dict[str, int]:
    out = {"Human": 153_613_673, "Curlcake": 10_135}
    path = TABLE_DIR / "controls_ivt_fpr.tsv"
    if path.exists():
        df = pd.read_csv(path, sep="\t")
        for species, sub in df.groupby("species"):
            vals = pd.to_numeric(sub["region_bp"], errors="coerce").dropna().unique()
            if len(vals) == 1:
                out[str(species)] = int(vals[0])
        log.info("  region_bp from %s: %s", path.name, out)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="R3-9 non-m6A evidence builder")
    ap.add_argument("--min-cov", type=int, default=10)
    ap.add_argument("--perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20260919)
    ap.add_argument("--outdir", type=Path, default=PKG / "evidence")
    ap.add_argument("--logdir", type=Path, default=PKG / "logs")
    ap.add_argument("--skip-scores", action="store_true")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.logdir.mkdir(parents=True, exist_ok=True)
    log = setup_logger("r39_build_evidence", log_dir=args.logdir)
    log.info("R3-9 evidence rebuild | min_cov=%d perm=%d seed=%d",
             args.min_cov, args.perm, args.seed)
    t_start = time.time()
    rng = np.random.default_rng(args.seed)
    res = Resources(args.min_cov)
    regions = region_bp_from_tables(log)

    # ---------------------------------------------------------------- layer 1
    log.info("[1/5] Curlcake unmodified controls")
    cc = build_curlcake(res, regions["Curlcake"], log)
    cc.to_csv(args.outdir / "curlcake_ivt_fp_per_construct.tsv", sep="\t", index=False)
    log.info("  -> %d rows", len(cc))

    # ---------------------------------------------------------------- layer 2
    log.info("[2/5] HeLa WT vs IVT per replicate")
    hela_rep, hela_sum = build_hela(res, regions["Human"], log)
    hela_rep.to_csv(args.outdir / "hela_wt_ivt_per_replicate.tsv", sep="\t", index=False)
    hela_sum.to_csv(args.outdir / "hela_wt_ivt_summary.tsv", sep="\t", index=False)
    log.info("  -> %d per-replicate rows", len(hela_rep))

    # ---------------------------------------------------------------- layer 3
    log.info("[3/5] truth-anchored precision (RMBase+DirectRMDB + NGS, permutation null)")
    refs: dict[str, dict[str, dict[str, np.ndarray]]] = {"RMBase+DirectRMDB": {}, "NGS": {}, "GLORI": {}}
    for mod, label in ORCA_MOD.items():
        refs["RMBase+DirectRMDB"][mod] = load_orca(label)
    for mod, (name, path) in NGS_REF.items():
        refs["NGS"][mod] = load_bed(path)
    refs["GLORI"]["m6A"] = load_bed(GLORI_HELA)
    truth = build_truth(res, refs, args.perm, rng, log)
    truth.to_csv(args.outdir / "truth_precision_per_replicate.tsv", sep="\t", index=False)
    log.info("  -> %d rows", len(truth))

    # ---------------------------------------------------------------- layer 4
    if not args.skip_scores:
        log.info("[4/5] score distributions")
        sc_sum, sc_hist = build_scores(res, log)
        sc_sum.to_csv(args.outdir / "score_distributions.tsv", sep="\t", index=False)
        sc_hist.to_csv(args.outdir / "score_histograms.tsv", sep="\t", index=False)
        log.info("  -> %d summary / %d histogram rows", len(sc_sum), len(sc_hist))

    # ---------------------------------------------------------------- layer 5
    log.info("[5/5] third-party GSE271571 (E. coli)")
    ec = build_ecoli(log)
    ec.to_csv(args.outdir / "ecoli_gse271571_prob.tsv", sep="\t", index=False)
    log.info("  -> %d rows", len(ec))

    log.info("done in %.1f min; tables -> %s", (time.time() - t_start) / 60, args.outdir)


if __name__ == "__main__":
    main()
