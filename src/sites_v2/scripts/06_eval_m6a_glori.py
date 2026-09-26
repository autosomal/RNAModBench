#!/usr/bin/env python3
"""06 -- m6A evaluation against GLORI, per sample x tool x window.

For every m6A callset on a GLORI species (Arabidopsis / Mouse / HeLa, RNA002 and
RNA004) the script computes TP/FP/FN/TN inside the sample's candidate universe
(``04_build_universe.py``, c >= 10) for the window sweep ``0,1,2,5,10,20,50``,
plus precision/recall/F1/MCC/specificity with site-level bootstrap CIs at the
primary windows (0 and 2).

It also produces the replicate analysis (R3-2 / E6) **without merging
biological replicates**: every replicate is evaluated in its own candidate
universe, between-replicate agreement is reported as pairwise Jaccard plus a
mean/SD summary, and sites are additionally broken down by how many of the
independent replicates detected them (``k_of_n`` tiers, k >= 2).  The union of
replicates is not evaluated as if it were a sample.

Outputs
-------
sites_v2/evaluation/tables/m6a_glori_confusion.tsv      long: sample x tool x window
sites_v2/evaluation/tables/m6a_localization_curve.tsv   extension vs exact accuracy
sites_v2/evaluation/tables/reproducibility.tsv   per tool: per-replicate rows,
                                                  k_of_n tiers and group summary

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/sites_v2/scripts/06_eval_m6a_glori.py [--species S] [--sample S]
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import (CALLSET_ROOT, GENOMES, GLORI, PRIMARY_WINDOW, SAMPLES,
                           SAMPLES_BY_NAME, TABLE_DIR, UNIVERSE_ROOT, WINDOWS,
                           independent_units, replicate_class)
from common.evaluation import (confusion_windows, load_reference, load_universe,
                               positions_in_universe, reference_in_universe,
                               strand_base_consistency)
from common.io_utils import read_table, write_table
from common.manifest import Inventory, log_time, setup_logger
from common.metrics import precision_recall_ci


def universe_path(sample) -> Path:
    plain = UNIVERSE_ROOT / sample.platform / sample.species / f"{sample.canonical}__universe.tsv"
    gz = plain.with_suffix(".tsv.gz")
    return plain if plain.exists() else gz


def _restrict_to(calls: dict[str, np.ndarray],
                 universe: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Keep only the positions that are inside ``universe`` (exact membership)."""
    out: dict[str, np.ndarray] = {}
    for chrom, pos in calls.items():
        u = universe.get(chrom)
        if u is None or u.size == 0:
            continue
        idx = np.clip(np.searchsorted(u, pos), 0, u.size - 1)
        keep = u[idx] == pos
        if keep.any():
            out[chrom] = pos[keep]
    return out


def genomic_calls(df: pd.DataFrame, genome: Path) -> tuple[pd.DataFrame, int]:
    """Rows whose chromosome exists in the genome FASTA (transcript-space rows dropped)."""
    from common.io_utils import FASTA, normalize_chrom_for_fasta

    fa_keys = FASTA.keys(genome)
    keep = []
    for c in df["chrom"].astype(str).unique():
        if normalize_chrom_for_fasta(fa_keys, c) is not None:
            keep.append(c)
    keep_set = set(keep)
    mask = df["chrom"].astype(str).isin(keep_set)
    return df[mask], int((~mask).sum())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", action="append")
    ap.add_argument("--sample", action="append")
    ap.add_argument("--min-cov", type=int, default=10)
    args = ap.parse_args()

    logger = setup_logger("06_eval_m6a_glori")
    inv = Inventory("06_eval_m6a_glori")
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    rows, curve_rows, repro_rows = [], [], []
    #: (platform, species, dataset_group) -> tool -> sample -> calls-in-universe
    group_calls: dict = defaultdict(lambda: defaultdict(dict))

    index = []
    for f in sorted(CALLSET_ROOT.rglob("*.tsv")):
        platform, species, group, mod_type, tool, fname = f.relative_to(CALLSET_ROOT).parts
        index.append({"platform": platform, "species": species, "group": group,
                      "mod_type": mod_type, "tool": tool, "sample": fname[:-4]})
    index_df = pd.DataFrame(index)
    logger.info("callset index: %d files, %d samples", len(index_df),
                index_df["sample"].nunique())

    with log_time(logger, "m6A GLORI evaluation"):
        for sample_name in sorted(index_df["sample"].unique()):
            sample = SAMPLES_BY_NAME.get(sample_name)
            if sample is None or not GLORI.get(sample.species):
                continue
            if args.species and sample.species not in args.species:
                continue
            if args.sample and sample.canonical not in args.sample:
                continue
            upath = universe_path(sample)
            if not upath.exists():
                logger.warning("[%s] no universe file, skipped", sample_name)
                continue
            genome = GENOMES[sample.species]
            universe = load_universe(upath, "m6A", args.min_cov)
            n_u = int(sum(v.size for v in universe.values()))
            glori = load_reference(GLORI[sample.species])
            ref_u = reference_in_universe(glori, universe)
            n_ref_u = int(sum(v.size for v in ref_u.values()))
            logger.info("[%s] universe=%d, testable GLORI=%d", sample_name, n_u, n_ref_u)

            # ---------- per-tool evaluation ---------------------------------
            tool_dir = CALLSET_ROOT / sample.platform / sample.species / \
                sample.dataset_group / "m6A"
            for tool_path in sorted(tool_dir.glob("*/*.tsv")) if tool_dir.exists() else []:
                tool = tool_path.parent.name
                if tool_path.stem != sample.canonical:
                    continue
                df = read_table(tool_path)
                # A header-only callset (parser == ``zero_sites``: raw present but
                # Total_Detected == 0 after the tool's own threshold, e.g. DRUMMER
                # Arabidopsis_fip37_rep2) still flows through the confusion math:
                # an empty call set gives tp = fp = 0, fn = n_reference, so the
                # replicate is represented as a real 0-detection result rather than
                # being silently dropped.  Only a callset that is genuinely
                # unreadable (no ``pos_raw`` column at all) is skipped.
                if "pos_raw" not in getattr(df, "columns", []):
                    continue
                df["pos_raw"] = pd.to_numeric(df["pos_raw"]).astype(np.int64)
                df_gen, n_nonge = genomic_calls(df, genome)
                calls_u, n_out = positions_in_universe(df_gen, universe)
                conf = confusion_windows(universe, ref_u, calls_u, WINDOWS)
                ok_strand, n_checked = strand_base_consistency(df_gen, "m6A", genome)
                if not calls_u:
                    logger.warning("  %-28s no calls inside the universe", tool)
                for w in WINDOWS:
                    c = conf[w]
                    rec = {"platform": sample.platform, "species": sample.species,
                           "dataset_group": sample.dataset_group, "sample": sample_name,
                           "replicate_tag": sample.replicate_tag, "tool": tool,
                           "window": w, **{k: c[k] for k in
                                           ("tp", "fp", "fn", "tn", "universe",
                                            "n_calls_in_universe", "n_reference_in_universe",
                                            "precision", "recall", "f1", "specificity", "mcc")},
                           "n_calls_total": len(df), "n_calls_out_of_universe": n_out,
                           "n_calls_transcript_space": n_nonge,
                           "pct_strand_base_ok": 100 * ok_strand / n_checked if n_checked else np.nan,
                           "min_cov": args.min_cov}
                    if w in (0, PRIMARY_WINDOW):
                        rec.update(precision_recall_ci(c["tp"], c["fp"], c["fn"]))
                    rows.append(rec)
                # localization curve (extension distance vs exact accuracy)
                base = conf[0]
                for w in WINDOWS:
                    c = conf[w]
                    exact_rate = (c["tp"] / c["n_calls_in_universe"]
                                  if c["n_calls_in_universe"] else np.nan)
                    loc_acc = (base["tp"] / c["tp"]) if c["tp"] else np.nan
                    curve_rows.append({
                        "platform": sample.platform, "species": sample.species,
                        "sample": sample_name, "tool": tool, "window": w,
                        "hit_rate": (c["tp"] / c["n_calls_in_universe"]
                                     if c["n_calls_in_universe"] else np.nan),
                        "exact_rate": exact_rate, "localization_accuracy": loc_acc,
                        "recall": c["recall"], "n_calls_in_universe": c["n_calls_in_universe"],
                    })
                group_calls[(sample.platform, sample.species,
                             sample.dataset_group)][tool][sample_name] = calls_u

    # ------------------------------------------------------------------ #
    # Replicate analysis (R3-2 / E6).  The unit of analysis is ONE
    # independent sequencing unit -- biological replicates are NOT merged into
    # a pseudo-sample, because a union callset has no sampling distribution of
    # its own and would answer the editor's independence concern with
    # pseudoreplication.  What is reported per dataset group:
    #   * one row per replicate, in that replicate's OWN candidate universe
    #     (the headline metrics), plus the between-replicate mean/SD;
    #   * pairwise Jaccard on the common testable universe;
    #   * reproducibility TIERS ``k_of_n`` (a site detected in >= k of the n
    #     independent replicates) -- a support-level breakdown, not a new
    #     "sample".  k = 1 (the union) is deliberately not emitted.
    # Nested subsets / same-run splits collapse to one unit and two studies are
    # concordance, not replication (see ``config.replicate_class``).
    # ------------------------------------------------------------------ #
    with log_time(logger, "replicate-level reproducibility"):
        for (platform, species, group), tool_calls in sorted(group_calls.items()):
            members = [s for s in SAMPLES if s.platform == platform
                       and s.species == species and s.dataset_group == group]
            if not members:
                continue
            klass = replicate_class(group)
            units = independent_units(members)
            if len(units) < 2:
                continue
            rep_member = {u: (u if u in SAMPLES_BY_NAME else ms[0])
                          for u, ms in units.items()}
            genome = GENOMES[species]
            glori = load_reference(GLORI[species])
            rep_universe: dict[str, dict[str, np.ndarray]] = {}
            for unit, member in rep_member.items():
                up = universe_path(SAMPLES_BY_NAME[member])
                if up.exists():
                    rep_universe[unit] = load_universe(up, "m6A", args.min_cov)
            if len(rep_universe) < 2:
                continue
            logger.info("[%s/%s] n_independent=%d class=%s", platform, group,
                        len(rep_universe), klass)

            # common testable universe is built per tool below (same for all tools
            # of a group, so it is cached after the first computation)
            for tool, by_sample in sorted(tool_calls.items()):
                avail = {u: by_sample[m] for u, m in rep_member.items()
                         if m in by_sample and u in rep_universe}
                own = [r for r in rows if r["tool"] == tool
                       and r["dataset_group"] == group
                       and r["window"] == PRIMARY_WINDOW]
                # (a) per-replicate rows, each in its own universe
                for rec in own:
                    repro_rows.append({**{k: rec[k] for k in
                                          ("platform", "species", "dataset_group",
                                           "tool", "tp", "fp", "fn", "tn",
                                           "precision", "recall", "f1", "mcc",
                                           "specificity")},
                                       "sample": rec["sample"],
                                       "replicate_tag": rec["replicate_tag"],
                                       "universe": rec["universe"],
                                       "n_sites": rec["n_calls_in_universe"],
                                       "n_replicates": len(rep_universe),
                                       "replicate_class": klass,
                                       "support": "per_replicate",
                                       "window": PRIMARY_WINDOW})
                # (b) between-replicate summary + tiers on the common universe
                if len(avail) < 2:
                    continue
                common = {}
                for chrom in set.intersection(*[set(u) for u in rep_universe.values()]):
                    arrays = [u[chrom] for u in rep_universe.values()]
                    acc = arrays[0]
                    for a in arrays[1:]:
                        acc = np.intersect1d(acc, a)
                    if acc.size:
                        common[chrom] = acc
                n_common = int(sum(v.size for v in common.values()))
                if not n_common:
                    continue
                ref_c = reference_in_universe(glori, common)
                sets = {u: _restrict_to(calls, common) for u, calls in avail.items()}
                keys = sorted(sets)
                jac = []
                for i in range(len(keys)):
                    for j in range(i + 1, len(keys)):
                        a, b = sets[keys[i]], sets[keys[j]]
                        inter = sum(np.intersect1d(a[c], b[c]).size
                                    for c in set(a) & set(b))
                        z = np.zeros(0, dtype=np.int64)
                        uni = sum(np.union1d(a.get(c, z), b.get(c, z)).size
                                  for c in set(a) | set(b))
                        jac.append(inter / uni if uni else np.nan)
                n_rep = len(sets)
                cnt: dict[str, dict[int, int]] = {}
                for chrom in {c for s in sets.values() for c in s}:
                    tally: dict[int, int] = {}
                    for s in sets.values():
                        for p in s.get(chrom, np.zeros(0, dtype=np.int64)):
                            tally[int(p)] = tally.get(int(p), 0) + 1
                    cnt[chrom] = tally
                tiers = {"all": n_rep, "majority": math.ceil(n_rep / 2)}
                for level, need in tiers.items():
                    calls_l = {c: np.array(sorted(p for p, k in t.items() if k >= need),
                                            dtype=np.int64)
                               for c, t in cnt.items()}
                    calls_l = {c: v for c, v in calls_l.items() if v.size}
                    conf = confusion_windows(common, ref_c, calls_l,
                                             [0, PRIMARY_WINDOW])
                    for w, c in conf.items():
                        repro_rows.append({
                            "platform": platform, "species": species,
                            "dataset_group": group, "tool": tool, "sample": "",
                            "replicate_tag": "", "window": w,
                            "n_sites": int(sum(v.size for v in calls_l.values())),
                            "n_common_universe": n_common,
                            "n_replicates": n_rep, "replicate_class": klass,
                            "support": f"{need}_of_{n_rep}_{level}",
                            **{k: c[k] for k in
                               ("tp", "fp", "fn", "tn", "precision", "recall",
                                "f1", "mcc", "specificity")}})
                pre = [r["precision"] for r in own]
                rec = [r["recall"] for r in own]
                repro_rows.append({
                    "platform": platform, "species": species, "dataset_group": group,
                    "tool": tool, "sample": "", "replicate_tag": "",
                    "window": PRIMARY_WINDOW, "support": "group_summary",
                    "n_replicates": 1 if klass == "cross_study" else n_rep,
                    "n_independent_datasets": n_rep, "replicate_class": klass,
                    "n_common_universe": n_common,
                    "replicates": "|".join(keys),
                    "mean_pairwise_jaccard": (float(np.nanmean(jac)) if jac else np.nan),
                    "min_pairwise_jaccard": (float(np.nanmin(jac)) if jac else np.nan),
                    "precision_mean": float(np.nanmean(pre)) if pre else np.nan,
                    "precision_sd": float(np.nanstd(pre, ddof=1)) if len(pre) > 1 else np.nan,
                    "recall_mean": float(np.nanmean(rec)) if rec else np.nan,
                    "recall_sd": float(np.nanstd(rec, ddof=1)) if len(rec) > 1 else np.nan,
                })
                logger.info("  repro %-28s n_indep=%d jaccard=%.3f", tool, n_rep,
                            float(np.nanmean(jac)) if jac else float("nan"))

    for name, data in (("m6a_glori_confusion.tsv", rows),
                       ("m6a_localization_curve.tsv", curve_rows),
                       ("reproducibility.tsv", repro_rows)):
        if data:
            write_table(pd.DataFrame(data), TABLE_DIR / name)
            inv.record(TABLE_DIR / name, n_rows=len(data))
            logger.info("%s: %d rows", name, len(data))
    inv.flush()
    logger.info("tables -> %s", TABLE_DIR)
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
