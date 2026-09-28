#!/usr/bin/env python3
"""12 -- Fig. 7 aligned non-m6A summary (HeLa + unmodified Curlcake).

The manuscript's non-m6A panel (Fig. 7A-C) is built on exactly two libraries:
the three HeLa WT / three HeLa IVT replicates and the unmodified Curlcake
constructs.  This script recomputes that panel from ``harmonisation/callsets`` and
puts it side by side with the numbers actually printed in the paper, so the
response letter can quote a single, reproducible table.

Definitions used (all stated in the output as well)
---------------------------------------------------
* **Fig. 7A "T/WT"** - numerically ``IVT union / WT union`` per tool (verified
  against the published counts: CHEUI_m5C 48,627; NanoPsu 883; NanoSPA_psU 888
  all match the **WT** union, and only NanoPsu / NanoSPA_psU fall below 1).
  The manuscript calls this "treatment/wild-type"; in this panel the treatment
  arm *is* the unmodified IVT library - flagged for the response letter.
* **Fig. 7C "Global Jaccard"** - ``|intersection over replicates| / |union over
  replicates|`` on the three replicates of one group (WT or IVT).
* **Fig. 7B false positives** - calls made on the unmodified Curlcake
  constructs.  ``Curlcake_IVT_rep2_partial`` is a depth-matched subset of
  ``Curlcake_IVT_rep3`` and is therefore excluded (reported separately).
  (Its pre-2026-09-16 working name is listed in
  ``harmonisation/manifest/curlcake_semantic_rename_map.csv``.)

Outputs
-------
harmonisation/evaluation/tables/nonm6a_fig7_summary.tsv
harmonisation/evaluation/tables/nonm6a_curlcake_detail.tsv

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/12_nonm6a_summary.py
"""

from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import (CALLSET_ROOT, CURLCAKE_FASTA, LEGACY_OUTPUT,
                           SAMPLES_BY_NAME, TABLE_DIR, UNIVERSE_ROOT)
from common.evaluation import load_universe
from common.io_utils import read_table, write_table
from common.manifest import Inventory, log_time, setup_logger
from common.match import fix_chromosome

#: HeLa replicates that form the non-m6A panel.
HELA_REPS = {
    "WT": ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"],
    "IVT": ["HeLa_IVT_rep1", "HeLa_IVT_rep2", "HeLa_IVT_rep3"],
}
#: Curlcake IVT constructs; ``Curlcake_IVT_rep2_partial`` (= the old
#: ``rep3part``) is a depth-matched subset of ``Curlcake_IVT_rep3`` and is
#: excluded from the union, reported separately (``construct_role``).
CURLCAKE_CONSTRUCTS = ["Curlcake_IVT_rep3", "Curlcake_IVT_rep1"]
CURLCAKE_SUBSET = ["Curlcake_IVT_rep2_partial"]

#: numbers actually printed in the manuscript (verbatim) for cross-checking.
PUBLISHED_FIG7A_SITES = {"CHEUI_m5C": 48627, "NanoPsu": 883, "NanoSPA_psU": 888}
PUBLISHED_GLOBAL_JACCARD = {"NanoNm": 0.233, "NanoMUD_m1psi": 0.170,
                            "CHEUI_m5C": 0.0002}

LEGACY_GROUP_DIR = {"WT": "HeLa_WT", "IVT": "HeLa_IVT"}


def hela_callset(mod_type: str, tool: str, sample: str) -> Path:
    return (CALLSET_ROOT / "RNA002" / "Human" / SAMPLES_BY_NAME[sample].dataset_group
            / mod_type / tool / f"{sample}.tsv")


def curlcake_callset(mod_type: str, tool: str, sample: str) -> Path:
    return CALLSET_ROOT / "RNA002" / "Curlcake" / "Curlcake_IVT" / mod_type / tool / f"{sample}.tsv"


def universe_path(sample_name: str) -> Path:
    s = SAMPLES_BY_NAME[sample_name]
    plain = (UNIVERSE_ROOT / s.platform / s.species /
             f"{s.canonical}__universe.tsv")
    gz = plain.with_suffix(".tsv.gz")
    return plain if plain.exists() else gz


def pos_set(path: Path) -> set[tuple[str, int]]:
    """{(chrom, 0-based pos)} of a callsets file."""
    df = read_table(path)
    if df.empty:
        return set()
    return {(fix_chromosome(c), int(p))
            for c, p in zip(df["chrom"], pd.to_numeric(df["pos_raw"]))}


def pos_set_universe(path: Path, mod_type: str, min_cov: int) -> set[tuple[str, int]]:
    """Same, restricted to the sample's own candidate universe."""
    df = read_table(path)
    if df.empty:
        return set()
    up = universe_path(path.stem)
    if not up.exists():
        return set()
    uni = load_universe(up, mod_type, min_cov)
    out = set()
    for c, p in zip(df["chrom"], pd.to_numeric(df["pos_raw"])):
        chrom = fix_chromosome(c)
        arr = uni.get(chrom)
        if arr is None or arr.size == 0:
            continue
        i = int(np.searchsorted(arr, p))
        if i < arr.size and arr[i] == p:
            out.add((chrom, int(p)))
    return out


def global_jaccard(sets: list[set]) -> float:
    """|intersection of all sets| / |union of all sets|."""
    sets = [s for s in sets if s is not None]
    if not sets:
        return np.nan
    inter = set.intersection(*sets)
    union = set.union(*sets)
    return len(inter) / len(union) if union else np.nan


def mean_pairwise_jaccard(sets: list[set]) -> float:
    vals = []
    for a, b in combinations(sets, 2):
        union = a | b
        if union:
            vals.append(len(a & b) / len(union))
    return float(np.mean(vals)) if vals else np.nan


def legacy_sets(group: str, tool: str) -> list[set]:
    """Replicate sets from the legacy ``output/`` tree (0-based Start column)."""
    d = LEGACY_OUTPUT / LEGACY_GROUP_DIR[group] / "tools" / "others" / tool
    if not d.exists():
        return []
    out = []
    for f in sorted(d.glob("*rep[123]*.txt")):
        df = read_table(f)
        if df.empty:
            continue
        out.append({(fix_chromosome(c), int(p))
                    for c, p in zip(df["Chr"], pd.to_numeric(df["Start"]))})
    return out


def discover_tools() -> list[tuple[str, str]]:
    """(mod_type, tool) pairs present for HeLa, sorted for stable output."""
    base = CALLSET_ROOT / "RNA002" / "Human" / "HeLa_WT"
    pairs = []
    for mod_dir in sorted(p for p in base.iterdir() if p.is_dir()):
        if mod_dir.name == "m6A":
            continue
        for tool_dir in sorted(p for p in mod_dir.iterdir() if p.is_dir()):
            pairs.append((mod_dir.name, tool_dir.name))
    return pairs


def region_bp_curlcake() -> int:
    fai = Path(str(CURLCAKE_FASTA) + ".fai")
    total = 0
    with fai.open() as fh:
        for line in fh:
            total += int(line.split("\t")[1])
    return total


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-cov", type=int, default=10)
    args = ap.parse_args()

    logger = setup_logger("12_nonm6a_summary")
    inv = Inventory("12_nonm6a_summary")
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    cc_bp = region_bp_curlcake()
    rows, cc_rows = [], []

    with log_time(logger, "Fig. 7 non-m6A summary"):
        for mod_type, tool in discover_tools():
            if not hela_callset(mod_type, tool, HELA_REPS["WT"][0]).exists():
                continue
            rec: dict = {"mod_type": mod_type, "tool": tool}
            sets_raw, sets_uni, legacy = {}, {}, {}
            for group, samples in HELA_REPS.items():
                raw, uni = [], []
                for s in samples:
                    p = hela_callset(mod_type, tool, s)
                    if not p.exists():
                        continue
                    raw.append(pos_set(p))
                    uni.append(pos_set_universe(p, mod_type, args.min_cov))
                sets_raw[group], sets_uni[group] = raw, uni
                gset = set.union(*raw) if raw else set()
                gset_u = set.union(*uni) if uni else set()
                rec[f"{group.lower()}_reps"] = len(raw)
                rec[f"{group.lower()}_union_sites"] = len(gset)
                rec[f"{group.lower()}_union_sites_in_universe"] = len(gset_u)
                rec[f"global_jaccard_{group.lower()}"] = global_jaccard(raw)
                rec[f"mean_pairwise_jaccard_{group.lower()}"] = mean_pairwise_jaccard(raw)
                legacy[group] = legacy_sets(group, tool)

            wt, ivt = rec.get("wt_union_sites", np.nan), rec.get("ivt_union_sites", np.nan)
            rec["ratio_twt"] = ivt / wt if wt else np.nan
            rec["ratio_wt_over_ivt"] = wt / ivt if ivt else np.nan
            lwt = len(set.union(*legacy["WT"])) if legacy.get("WT") else np.nan
            livt = len(set.union(*legacy["IVT"])) if legacy.get("IVT") else np.nan
            rec["legacy_wt_union_sites"] = lwt
            rec["legacy_ivt_union_sites"] = livt
            rec["legacy_ratio_twt"] = (livt / lwt) if lwt == lwt and lwt else np.nan
            rec["legacy_global_jaccard_wt"] = global_jaccard(legacy.get("WT", []))
            rec["legacy_global_jaccard_ivt"] = global_jaccard(legacy.get("IVT", []))

            # -------------------------------------------------- Fig. 7B --------
            cc_sets, n_uni = [], 0
            for s in CURLCAKE_CONSTRUCTS:
                p = curlcake_callset(mod_type, tool, s)
                if not p.exists():
                    continue
                st = pos_set(p)
                st_u = pos_set_universe(p, mod_type, args.min_cov)
                cc_sets.append(st)
                up = universe_path(s)
                if up.exists():
                    n_uni = max(n_uni, int(sum(v.size for v in
                                              load_universe(up, mod_type, args.min_cov).values())))
                cc_rows.append({"mod_type": mod_type, "tool": tool, "construct": s,
                                "construct_role": "independent",
                                "n_sites": len(st),
                                "n_sites_in_universe": len(st_u)})
            for s in CURLCAKE_SUBSET:
                p = curlcake_callset(mod_type, tool, s)
                if p.exists():
                    cc_rows.append({"mod_type": mod_type, "tool": tool, "construct": s,
                                    "construct_role": "subset_of_Curlcake_IVT_rep3",
                                    "n_sites": len(pos_set(p)),
                                    "n_sites_in_universe": len(pos_set_universe(p, mod_type, args.min_cov))})
            cc_union = set.union(*cc_sets) if cc_sets else set()
            rec["curlcake_ivt_constructs"] = ",".join(
                s for s in CURLCAKE_CONSTRUCTS if curlcake_callset(mod_type, tool, s).exists())
            rec["curlcake_ivt_union_sites"] = len(cc_union)
            rec["curlcake_ivt_per_10kb"] = (1e4 * len(cc_union) / cc_bp) if cc_bp else np.nan
            rec["curlcake_ivt_per_1e6_candidates"] = (1e6 * len(cc_union) / n_uni) if n_uni else np.nan
            rec["curlcake_ivt_universe"] = n_uni

            # -------------------------------------------------- published ------
            rec["published_fig7a_sites"] = PUBLISHED_FIG7A_SITES.get(tool, "")
            rec["published_global_jaccard"] = PUBLISHED_GLOBAL_JACCARD.get(tool, "")
            rec["matches_published_count"] = (
                "" if tool not in PUBLISHED_FIG7A_SITES
                else int(rec["wt_union_sites"]) == PUBLISHED_FIG7A_SITES[tool])
            rec["matches_published_jaccard"] = (
                "" if tool not in PUBLISHED_GLOBAL_JACCARD
                else bool(abs(rec["global_jaccard_wt"]
                              - PUBLISHED_GLOBAL_JACCARD[tool]) < 5e-4))
            rec["fig7a_twt_definition"] = "IVT union / WT union (III)"
            rows.append(rec)
            logger.info("  %-8s %-16s WT=%-7d IVT=%-7d T/WT=%.3f  gJ(WT)=%.4f  Curlcake=%d",
                        mod_type, tool, rec["wt_union_sites"], rec["ivt_union_sites"],
                        rec["ratio_twt"], rec["global_jaccard_wt"], len(cc_union))

    summary = pd.DataFrame(rows)
    if len(summary):
        write_table(summary, TABLE_DIR / "nonm6a_fig7_summary.tsv")
        inv.record(TABLE_DIR / "nonm6a_fig7_summary.tsv", n_rows=len(summary))
    detail = pd.DataFrame(cc_rows)
    if len(detail):
        write_table(detail, TABLE_DIR / "nonm6a_curlcake_detail.tsv")
        inv.record(TABLE_DIR / "nonm6a_curlcake_detail.tsv", n_rows=len(detail))
    inv.flush()

    if len(summary):
        ok = summary[summary["matches_published_count"] == True]  # noqa: E712
        logger.info("published Fig. 7A counts reproduced for %d/%d tools that quote one",
                    len(ok), int((summary["matches_published_count"] != "").sum()))
        jok = summary[summary["matches_published_jaccard"] == True]  # noqa: E712
        logger.info("published Global Jaccard reproduced for %d/%d tools that quote one",
                    len(jok), int((summary["matches_published_jaccard"] != "").sum()))
    logger.info("tables -> %s", TABLE_DIR)
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
