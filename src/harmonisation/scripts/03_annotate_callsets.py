#!/usr/bin/env python3
"""03 -- annotate every callset in place.

Adds, per site (see ``common/annotate.py`` for the exact conventions):

``ref_base five_mer_raw five_mer_center center_base_expected is_drach_raw
is_drach_center pos_center dist_center dist_drach_a coverage cov_source
dist_glori glori_exact offset_flag``

Coverage is computed once per sample from the sample's genomic BAM
(``resolve_coverage_bam`` / ``rna004_coverage_bam``) over the union of the
sample's called positions and merged back onto every tool of that sample.
Curlcake RNA004 has no local BAM; it falls back to the modkit pileup
``valid_coverage`` carried in the source rows.

The script is idempotent: previously added annotation columns are dropped
before re-annotating, so ``01 -> 03`` can be re-run at any time.

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/03_annotate_callsets.py [--sample S] [--tool T]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.annotate import ANNOTATION_COLUMNS, annotate_frame, load_glori
from common.bamcov import coverage_at_positions
from common.config import (CALLSET_ROOT, GENOMES, GLORI, MANIFEST_DIR,
                           SAMPLES_BY_NAME, resolve_coverage_bam, rna004_coverage_bam)
from common.io_utils import read_table, rel, write_table
from common.manifest import Inventory, log_time, setup_logger


def list_callsets() -> pd.DataFrame:
    rows = []
    for f in sorted(CALLSET_ROOT.rglob("*.tsv")):
        parts = f.relative_to(CALLSET_ROOT).parts
        # <platform>/<species>/<dataset_group>/<mod_type>/<tool>/<sample>.tsv
        platform, species, group, mod_type, tool, fname = parts
        rows.append({"platform": platform, "species": species, "dataset_group": group,
                     "mod_type": mod_type, "tool": tool,
                     "sample": fname[:-4], "path": f})
    return pd.DataFrame(rows)


def sample_bam(sample):
    """(path, source-label) of the genomic BAM used for coverage."""
    if sample.platform == "RNA004":
        p = rna004_coverage_bam(sample)
        return (p, "rna004_minimap2_bam") if p else (None, "")
    return resolve_coverage_bam(sample)


def build_coverage_table(sample, frames: list[pd.DataFrame], logger):
    """Union-position coverage table for one sample: DataFrame(chrom,pos_raw,coverage)."""
    path, src = sample_bam(sample)
    if path is None:
        return None, "no_bam"
    allpos = pd.concat([f[["chrom", "pos_raw"]] for f in frames], ignore_index=True)
    allpos["pos_raw"] = pd.to_numeric(allpos["pos_raw"]).astype(np.int64)
    allpos = allpos.drop_duplicates().reset_index(drop=True)
    chrom_positions = {c: sub["pos_raw"].to_numpy() for c, sub in allpos.groupby("chrom")}
    logger.info("[%s] coverage: %s (%d unique positions)", sample.canonical,
                Path(path).name, len(allpos))
    cov_map = coverage_at_positions(path, chrom_positions)

    parts = []
    for chrom, pos in chrom_positions.items():
        parts.append(pd.DataFrame({"chrom": chrom, "pos_raw": pos,
                                   "coverage": cov_map[chrom]}))
    table = pd.concat(parts, ignore_index=True)
    return table, src


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", action="append")
    ap.add_argument("--tool", action="append")
    ap.add_argument("--no-coverage", action="store_true",
                    help="skip the samtools-depth step and keep the coverage "
                         "column already stored in the callset (2026-09-18: "
                         "the strand / centre-base correction only needs the "
                         "reference base, and a full depth pass over the 10 GB "
                         "RNA004 BAMs costs hours)")
    args = ap.parse_args()

    logger = setup_logger("03_annotate_callsets")
    inv = Inventory("03_annotate_callsets")

    index = list_callsets()
    if args.sample:
        index = index[index["sample"].isin(args.sample)]
    if args.tool:
        index = index[index["tool"].isin(args.tool)]

    summaries = []
    with log_time(logger, "annotation"):
        for sample_name, grp in index.groupby("sample"):
            sample = SAMPLES_BY_NAME.get(sample_name)
            if sample is None:
                logger.warning("unknown sample in callset tree: %s", sample_name)
                continue
            frames = []
            for _, r in grp.iterrows():
                df = read_table(r["path"])
                if df.empty:
                    continue
                df["pos_raw"] = pd.to_numeric(df["pos_raw"]).astype(np.int64)
                #: ``--no-coverage`` keeps the depth already stored in the file
                #: (computed from the sample's own BAM on an earlier run); only
                #: the strand / centre-base columns are recomputed.
                keep_cov = None
                keep_src = ""
                if args.no_coverage and "coverage" in df.columns:
                    keep_cov = pd.to_numeric(df["coverage"], errors="coerce").to_numpy()
                    keep_src = (str(df["cov_source"].iloc[0])
                                if "cov_source" in df.columns else "unchanged")
                frames.append((r, df, keep_cov, keep_src))
            if not frames:
                continue

            cov_table, cov_source = (None, "no_coverage") if args.no_coverage else \
                build_coverage_table(sample, [f for _, f, _, _ in frames], logger)
            glori = None
            if GLORI.get(sample.species):
                glori = load_glori(GLORI[sample.species])

            for r, df, keep_cov, keep_src in frames:
                # idempotency: drop any annotation column from a previous run so
                # the coverage merge cannot collide with an existing 'coverage'
                df = df.drop(columns=[c for c in ANNOTATION_COLUMNS if c in df.columns])
                cov = keep_cov
                src = keep_src or "no_coverage"
                if keep_cov is not None:
                    pass          # --no-coverage: keep the depth already stored
                elif cov_table is not None:
                    merged = df.merge(cov_table, on=["chrom", "pos_raw"], how="left")
                    cov = pd.to_numeric(merged["coverage"], errors="coerce").to_numpy()
                    src = cov_source
                elif "src_coverage" in df.columns:
                    cov = pd.to_numeric(df["src_coverage"], errors="coerce").to_numpy()
                    src = f"src_coverage({r['tool']})"
                else:
                    src = "no_coverage"

                out = annotate_frame(df, genome=GENOMES[sample.species],
                                     mod_type=r["mod_type"], coverage=cov,
                                     cov_source=src, glori=glori)
                write_table(out, r["path"])
                summaries.append({
                    "sample": sample_name, "tool": r["tool"], "mod_type": r["mod_type"],
                    "rows": len(out),
                    "pct_ref_base_found": 100 * (out["ref_base"].astype(str) != "").mean(),
                    "pct_expected_base_at_raw": 100 * out["center_base_expected"].mean(),
                    "pct_is_drach_center": 100 * out["is_drach_center"].mean(),
                    "median_abs_dist_center": float(np.nanmedian(np.abs(out["dist_center"]))),
                    "median_coverage": float(np.nanmedian(out["coverage"])),
                    "pct_glori_exact": (100 * out["glori_exact"].mean()) if glori else np.nan,
                    "cov_source": src,
                })
                logger.info("  %-22s %-30s rows=%-7d expA=%5.1f%% drach=%5.1f%% medcov=%6.0f",
                            sample_name, r["tool"], len(out),
                            summaries[-1]["pct_expected_base_at_raw"],
                            summaries[-1]["pct_is_drach_center"],
                            summaries[-1]["median_coverage"])

    sdf = pd.DataFrame(summaries)
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    write_table(sdf, MANIFEST_DIR / "annotation_summary.csv")
    inv.record(MANIFEST_DIR / "annotation_summary.csv", n_rows=len(sdf))
    inv.flush()

    logger.info("annotated %d callsets (%d rows)", len(sdf), int(sdf["rows"].sum()) if len(sdf) else 0)
    logger.info("summary -> %s", MANIFEST_DIR / "annotation_summary.csv")
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
