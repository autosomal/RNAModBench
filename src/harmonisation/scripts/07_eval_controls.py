#!/usr/bin/env python3
"""07 -- negative controls (IVT), partial negative controls (KO/KD) and purified sites.

Panels
------
1. **IVT false positives** (R3-8, R3-9): in an unmodified library every call is a
   false positive by construction.  Because absolute counts depend on library
   size, the FP burden is normalised three ways:
   ``fp_per_10kb`` (per 10 kb of annotated exon sequence),
   ``fp_per_1e6_candidates`` (per 10^6 universe positions) and
   ``fp_fraction_of_universe``.
   Panels: Curlcake IVT (RNA002 + RNA004), HeLa IVT (RNA002 + RNA004),
   E. coli IVT (RNA002).

2. **KO / KD partial negatives** (R1-4): the KO/KD libraries are NOT clean
   negatives.  For each tool the WT and the matched KO/KD sample are compared on
   the *intersection* of their testable universes, so coverage/expression
   differences cannot masquerade as specificity:
   ``n_wt_common``, ``n_ko_common``, ``n_shared``, ``wt_only``, ``ko_only`` and
   ``ko_fraction`` (share of common-universe calls that are KO-only).

3. **Purified sites** (R3-7): sites called in WT on the common universe and
   absent from the KO/KD callset, with coverage of both sides recorded.  The
   circularity caveat is stated in the README.

Outputs
-------
harmonisation/evaluation/tables/controls_ivt_fpr.tsv
harmonisation/evaluation/tables/ko_kd_metrics.tsv
harmonisation/evaluation/tables/purified_sites.tsv

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/07_eval_controls.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import (CALLSET_ROOT, CURLCAKE_FASTA, GENOMES, SAMPLES,
                           SAMPLES_BY_NAME, TABLE_DIR, UNIVERSE_ROOT,
                           independence_class, sequencing_unit)
from common.evaluation import (MOD_GENOME_BASES, load_universe,
                               positions_in_universe)
from common.io_utils import read_table, write_table
from common.manifest import Inventory, log_time, setup_logger
from common.refs import build_exon_bed

#: (WT sample, KO/KD sample) pairs that belong to the same study
KO_PAIRS = [
    ("Arabidopsis_WT_rep1", "Arabidopsis_fip37_rep1"),
    ("Arabidopsis_WT_rep2", "Arabidopsis_fip37_rep2"),
    ("Arabidopsis_WT_rep3", "Arabidopsis_fip37_rep3"),
    ("mES_WT", "mES_KO"),
    ("mESCs_Mettl3_WT", "mESCs_Mettl3_KO"),
    ("HeLa_WT1", "HeLa_IVT_rep1"),   # not a KO - excluded below by condition class
    ("HeLa_WT2", "HeLa_IVT_rep2"),
    ("HeLa_WT3", "HeLa_IVT_rep3"),
]


def universe_path(sample) -> Path:
    plain = UNIVERSE_ROOT / sample.platform / sample.species / f"{sample.canonical}__universe.tsv"
    gz = plain.with_suffix(".tsv.gz")
    return plain if plain.exists() else gz


def region_sizes(species: str) -> dict[str, int]:
    """Annotated region length per chromosome (universe denominator for FP/10kb)."""
    if species == "Curlcake":
        sizes = {}
        with open(str(CURLCAKE_FASTA) + ".fai") as fh:
            for line in fh:
                f = line.split("\t")
                sizes[f[0]] = int(f[1])
        return sizes
    bed = build_exon_bed(species)
    sizes: dict[str, int] = {}
    with open(bed) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sizes[f[0]] = sizes.get(f[0], 0) + (int(f[2]) - int(f[1]))
    return sizes


def listsets(base: Path) -> list[tuple[str, str, Path, str]]:
    """[(sample, tool, path, mod_type)] for every callset under CALLSET_ROOT."""
    out = []
    for f in sorted(base.rglob("*.tsv")):
        platform, species, group, mod_type, tool, fname = f.relative_to(base).parts
        out.append((fname[:-4], tool, f, mod_type))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-cov", type=int, default=10)
    args = ap.parse_args()

    logger = setup_logger("07_eval_controls")
    inv = Inventory("07_eval_controls")
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    all_sets = listsets(CALLSET_ROOT)

    # ------------------------------------------------------------------ IVT --
    fpr_rows = []
    for sample_name, tool, path, mod_type in all_sets:
        sample = SAMPLES_BY_NAME.get(sample_name)
        if sample is None or sample.condition_class != "IVT":
            continue
        up = universe_path(sample)
        df = read_table(path)
        if df.empty:
            continue
        df["pos_raw"] = pd.to_numeric(df["pos_raw"]).astype(np.int64)
        n_calls = len(df)
        universe = load_universe(up, mod_type, args.min_cov) if up.exists() else {}
        n_u = int(sum(v.size for v in universe.values()))
        calls_u, n_out = positions_in_universe(df, universe) if universe else ({}, n_calls)
        sizes = region_sizes(sample.species)
        region_bp = int(sum(sizes.values()))
        fpr_rows.append({
            "platform": sample.platform, "species": sample.species,
            "sample": sample_name, "dataset_group": sample.dataset_group,
            "sequencing_unit": sequencing_unit(sample),
            "independence_class": independence_class(sample),
            "mod_type": mod_type, "tool": tool, "n_calls": n_calls,
            "n_calls_in_universe": int(sum(v.size for v in calls_u.values())),
            "n_calls_out_of_universe": n_out,
            "n_universe": n_u,
            "region_bp": region_bp,
            "fp_per_10kb": 1e4 * n_calls / region_bp if region_bp else np.nan,
            "fp_per_1e6_candidates": 1e6 * n_calls / n_u if n_u else np.nan,
            "fp_fraction_of_universe": n_calls / n_u if n_u else np.nan,
            "min_cov": args.min_cov,
        })
    if fpr_rows:
        df_fpr = pd.DataFrame(fpr_rows)
        write_table(df_fpr, TABLE_DIR / "controls_ivt_fpr.tsv")
        inv.record(TABLE_DIR / "controls_ivt_fpr.tsv", n_rows=len(df_fpr))
        logger.info("IVT panel: %d rows", len(df_fpr))

    # ------------------------------------------------------------- KO / KD --
    ko_rows = []
    for wt_name, ko_name in KO_PAIRS:
        wt = SAMPLES_BY_NAME.get(wt_name)
        ko = SAMPLES_BY_NAME.get(ko_name)
        if wt is None or ko is None or ko.condition_class not in ("KO", "KD"):
            continue
        wt_up, ko_up = universe_path(wt), universe_path(ko)
        if not wt_up.exists() or not ko_up.exists():
            continue
        tools = [t for s, t, _, _ in all_sets if s == wt_name]
        for tool in sorted(set(tools)):
            wt_files = [p for s, t, p, _ in all_sets if s == wt_name and t == tool]
            ko_files = [p for s, t, p, _ in all_sets if s == ko_name and t == tool]
            if not wt_files or not ko_files:
                continue
            mod_type = [m for s, t, _, m in all_sets if s == wt_name and t == tool][0]
            u_wt = load_universe(wt_up, mod_type, args.min_cov)
            u_ko = load_universe(ko_up, mod_type, args.min_cov)
            common = {}
            for chrom in set(u_wt) & set(u_ko):
                common[chrom] = np.intersect1d(u_wt[chrom], u_ko[chrom])
            n_common = int(sum(v.size for v in common.values()))

            def _u(df: pd.DataFrame) -> dict[str, np.ndarray]:
                df = df.copy()
                df["pos_raw"] = pd.to_numeric(df["pos_raw"]).astype(np.int64)
                calls, _ = positions_in_universe(df, common)
                return calls

            c_wt = _u(read_table(wt_files[0]))
            c_ko = _u(read_table(ko_files[0]))
            shared = sum(np.intersect1d(c_wt.get(ch, np.zeros(0, np.int64)),
                                        c_ko.get(ch, np.zeros(0, np.int64))).size
                         for ch in sorted(set(c_wt) | set(c_ko)))
            n_wt = int(sum(v.size for v in c_wt.values()))
            n_ko = int(sum(v.size for v in c_ko.values()))
            ko_rows.append({
                "species": wt.species, "dataset_group_wt": wt.dataset_group,
                "sample_wt": wt_name, "sample_ko": ko_name,
                "study_wt": wt.study, "study_ko": ko.study,
                "same_study": int(wt.study == ko.study),
                "mod_type": mod_type, "tool": tool,
                "n_common_universe": n_common,
                "n_wt_common": n_wt, "n_ko_common": n_ko, "n_shared": int(shared),
                "wt_only": n_wt - int(shared), "ko_only": n_ko - int(shared),
                "ko_fraction": n_ko / n_common if n_common else np.nan,
                "wt_fraction": n_wt / n_common if n_common else np.nan,
                "ko_wt_ratio": (n_ko / n_wt) if n_wt else np.nan,
            })
    if ko_rows:
        df_ko = pd.DataFrame(ko_rows)
        write_table(df_ko, TABLE_DIR / "ko_kd_metrics.tsv")
        inv.record(TABLE_DIR / "ko_kd_metrics.tsv", n_rows=len(df_ko))
        logger.info("KO/KD panel: %d rows", len(df_ko))
        purified = df_ko[df_ko["ko_only"] >= 0].copy()
        purified = purified.rename(columns={"wt_only": "n_purified_sites"})
        purified = purified[["species", "sample_wt", "sample_ko", "tool", "mod_type",
                             "n_wt_common", "n_ko_common", "n_shared",
                             "n_purified_sites", "ko_only", "n_common_universe"]]
        purified["criterion"] = ("called in WT and absent from the KO/KD callset, "
                                 "both restricted to the common testable universe; "
                                 "circularity caveat per R3-7 (see README)")
        write_table(purified, TABLE_DIR / "purified_sites.tsv")
        inv.record(TABLE_DIR / "purified_sites.tsv", n_rows=len(purified))
        logger.info("purified sites table: %d rows", len(purified))

    inv.flush()
    logger.info("tables -> %s", TABLE_DIR)
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
