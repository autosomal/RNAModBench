#!/usr/bin/env python3
"""05 -- per-tool positional-offset audit (the "site displacement" check).

For every (sample, tool) the signed distance ``pos_raw - nearest_reference`` is
computed against the sample's positive reference:

* m6A species with GLORI -> GLORI bed (liftover version for mouse);
* Curlcake fully-modified constructs -> all DRACH A sites of ``cc.fasta``;
* everything else -> the audit is skipped (no positives; ``n_audited = 0``).

The audit reports, per (species, tool):

``n_audited`` (calls within +/-5 bp of a reference), ``mode_offset`` (most
frequent signed distance), ``mode_share`` (its fraction among audited calls),
``pct_exact`` (share at distance 0), ``pct_within_2``, and a verdict
``offset_flag`` in ``{none, systematic_+k, systematic_-k, diffuse}``.

Verdict rule (documented in README): a tool is flagged ``systematic_<k>`` when
``k != 0`` is the mode, ``mode_share >= 0.40`` and ``pct_exact < 0.60``.
``diffuse`` means >30 % of calls sit within +/-5 bp but without a dominant
offset (signals positional noise rather than a coordinate convention).

The per-row ``offset_flag`` column in every callset is filled with
``systematic_<k>`` for rows whose distance equals the tool's mode and ``""``
otherwise, so downstream analyses can filter/annotate shifted calls.

Outputs
-------
harmonisation/evaluation/tables/offset_audit.tsv
harmonisation/evaluation/tables/offset_histogram.tsv

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/05_offset_audit.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.annotate import load_glori
from common.config import (CALLSET_ROOT, CURLCAKE_FASTA, GENOMES, GLORI,
                           MANIFEST_DIR, TABLE_DIR)
from common.io_utils import FASTA, read_table, write_table
from common.manifest import Inventory, log_time, setup_logger
from common.match import fix_chromosome, signed_nearest_distance
from common.refs import build_exon_bed, curlcake_drach_a  # noqa: F401

AUDIT_WINDOW = 5          # calls further than this are not informative
MODE_MIN_SHARE = 0.40
EXACT_MAX_FOR_FLAG = 0.60
DIFFUSE_MIN_NEARBY = 0.30

def audit_reference(sample: str, species: str, dataset_group: str, mod_type: str, logger):
    """(reference dict, label) used for the offset audit of one sample.

    GLORI is m6A truth, so it audits only m6A callsets; Curlcake DRACH-A is the
    m6A truth for Curlcake m6A samples. Other modification types (m5C, Psi, Nm,
    inosine) have no single-base truth reference here and are skipped.
    Reference keys are normalised with ``fix_chromosome`` so they line up with
    the callset ``chrom`` column (e.g. ``chrcurlcake1``).
    """
    if mod_type != "m6A":
        return None, "none"
    if GLORI.get(species):
        return load_glori(GLORI[species]), Path(GLORI[species]).name
    if species == "Curlcake" and "m6A" in dataset_group:
        raw = curlcake_drach_a(CURLCAKE_FASTA)
        return {fix_chromosome(c): v for c, v in raw.items()}, "curlcake_DRACH_A"
    return None, "none"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", action="append")
    args = ap.parse_args()

    logger = setup_logger("05_offset_audit")
    inv = Inventory("05_offset_audit")
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    audit_rows, hist_rows = [], []

    index = []
    for f in sorted(CALLSET_ROOT.rglob("*.tsv")):
        platform, species, group, mod_type, tool, fname = f.relative_to(CALLSET_ROOT).parts
        index.append({"platform": platform, "species": species, "dataset_group": group,
                      "mod_type": mod_type, "tool": tool, "sample": fname[:-4],
                      "path": f})

    with log_time(logger, "offset audit"):
        for r in index:
            if args.sample and r["sample"] not in args.sample:
                continue
            df = read_table(r["path"])
            if df.empty:
                continue
            df["pos_raw"] = pd.to_numeric(df["pos_raw"]).astype(np.int64)
            ref, ref_label = audit_reference(r["sample"], r["species"],
                                             r["dataset_group"], r["mod_type"], logger)
            dist = np.full(len(df), np.nan)
            if ref:
                for chrom, sub in df.groupby("chrom"):
                    g = ref.get(chrom)
                    if g is None or g.size == 0:
                        continue
                    dist[sub.index.to_numpy()] = signed_nearest_distance(
                        sub["pos_raw"].to_numpy(), g)
            near = np.abs(dist) <= AUDIT_WINDOW
            n_near = int(near.sum())
            mode_offset, mode_share, pct_exact, pct_within_2 = np.nan, np.nan, np.nan, np.nan
            flag = "no_reference" if not ref else ("none" if n_near == 0 else "diffuse")
            if n_near:
                d_near = dist[near].astype(int)
                pct_exact = float((d_near == 0).mean())
                pct_within_2 = float((np.abs(d_near) <= 2).mean())
                counts = pd.Series(d_near).value_counts()
                mode_offset = int(counts.index[0])
                mode_share = float(counts.iloc[0] / n_near)
                if (mode_offset != 0 and mode_share >= MODE_MIN_SHARE
                        and pct_exact < EXACT_MAX_FOR_FLAG):
                    flag = f"systematic_{mode_offset:+d}"
                elif pct_exact >= 0.5 or n_near < 50:
                    flag = "none"
                elif pct_within_2 >= 0.5:
                    flag = "diffuse"
                else:
                    flag = "none"

                for off, n in counts.items():
                    hist_rows.append({"sample": r["sample"], "tool": r["tool"],
                                      "species": r["species"], "offset": int(off),
                                      "n": int(n)})

            audit_rows.append({
                "sample": r["sample"], "species": r["species"],
                "dataset_group": r["dataset_group"], "mod_type": r["mod_type"],
                "tool": r["tool"], "n_calls": len(df), "n_audited": n_near,
                "reference": ref_label, "mode_offset": mode_offset,
                "mode_share": mode_share, "pct_exact": pct_exact,
                "pct_within_2": pct_within_2, "offset_flag": flag,
            })

            # write the per-row flag back into the callset: only rows sitting at
            # the tool's dominant offset keep the flag
            if ref:
                if flag.startswith("systematic") and not np.isnan(mode_offset):
                    df["offset_flag"] = np.where(near & (dist == mode_offset), flag, "")
                else:
                    df["offset_flag"] = ""
                write_table(df, r["path"])

            logger.info("  %-22s %-30s n=%-7d audited=%-6d mode=%s share=%.2f flag=%s",
                        r["sample"], r["tool"], len(df), n_near,
                        "-" if np.isnan(mode_offset) else f"{mode_offset:+.0f}",
                        -1 if np.isnan(mode_share) else mode_share, flag)

    adf = pd.DataFrame(audit_rows)
    hdf = pd.DataFrame(hist_rows)
    write_table(adf, TABLE_DIR / "offset_audit.tsv")
    write_table(hdf, TABLE_DIR / "offset_histogram.tsv")
    inv.record(TABLE_DIR / "offset_audit.tsv", n_rows=len(adf))
    inv.record(TABLE_DIR / "offset_histogram.tsv", n_rows=len(hdf))
    inv.flush()

    if len(adf):
        summary = (adf.groupby(["species", "tool"])["offset_flag"]
                   .agg(lambda s: s.value_counts().to_dict()))
        logger.info("offset verdicts (per tool, all samples):\n%s", summary.to_string())
    logger.info("audit -> %s", TABLE_DIR / "offset_audit.tsv")
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
