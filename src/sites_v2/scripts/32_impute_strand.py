#!/usr/bin/env python3
"""32 -- impute the transcript strand for callsets that carry none.

Why
---
``common/annotate.py`` treats any strand that is not exactly ``'-'`` as ``'+'``.
Nanom6A writes ``'*'`` for **every** row (25 callsets, ~286 k rows), so its
minus-strand calls -- where the modified adenosine shows up as a genomic ``T``
-- were annotated with a bogus ``dist_center`` (-5..+5) and a ``pos_center``
2-5 bp off the real site.  The calls themselves are correct: reverse-
complementing their 5mer lifts the DRACH share from 0.000 to 0.734
(Arabidopsis) / 0.998 (HeLa), exactly the share the plus-strand rows show.

This step fills the strand from the species exon annotation
(``common/strand.py``, same ``config.GTF_EXON`` the candidate universe is built
from) before ``03_annotate_callsets.py`` runs, so the strand-aware columns are
computed on the right strand.  It never moves a coordinate.

Columns written
---------------
``strand``      ``'+'``/``'-'`` when known or imputed, ``'*'`` when unresolved
``strand_src``  ``tool`` / ``imputed_read_bed`` / ``imputed_exon`` / ``unresolved``

Two evidence sources, in this order
-----------------------------------
1. the tool's own read-alignment BED (``config.READ_STRAND_BED``) -- the strand
   of the reads the calls came from;
2. the species exon annotation (``common/strand.py``).

Run BEFORE ``03`` (idempotent: rows that already carry a strand are left alone)
and AFTER ``02b`` (which compares the rebuilt callsets against the legacy
converted layer row-for-row and must see the extraction as it came out of 01).

Outputs
-------
sites_v2/manifest/strand_imputation.csv   per callset: rows / unknown before /
                                          imputed + / imputed - / ambiguous /
                                          no_exon / unknown after

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/sites_v2/scripts/32_impute_strand.py [--dry-run] [--sample S]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import CALLSET_ROOT, MANIFEST_DIR, READ_STRAND_BED, RESULT_RNA002  # noqa: E402
from common.io_utils import read_table, write_table  # noqa: E402
from common.manifest import Inventory, log_time, setup_logger  # noqa: E402
from common.strand import (BOTH, MINUS, NO_EXON, PLUS, impute_frame,  # noqa: E402
                           load_alignment_index)


def iter_callsets():
    """(path, platform, species, group, mod_type, tool, sample) for every callset."""
    for f in sorted(CALLSET_ROOT.rglob("*.tsv")):
        parts = f.relative_to(CALLSET_ROOT).parts
        if len(parts) != 6:
            continue
        platform, species, group, mod_type, tool, fname = parts
        yield f, platform, species, group, mod_type, tool, fname[:-4]


_ALIGN_CACHE: dict[str, dict] = {}


def alignment_index_for(tool: str, sample: str) -> tuple[dict, str]:
    """(strand index, label) from the tool's own read-alignment BED, if any."""
    spec = READ_STRAND_BED.get(tool)
    if spec is None:
        return {}, "none"
    path = RESULT_RNA002 / spec[0] / sample / spec[1]
    if not path.exists():
        return {}, "missing"
    if str(path) not in _ALIGN_CACHE:
        _ALIGN_CACHE[str(path)] = load_alignment_index(path)
    return _ALIGN_CACHE[str(path)], "read_bed"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="measure only; do not touch the callsets")
    ap.add_argument("--sample", action="append", default=None)
    args = ap.parse_args()

    logger = setup_logger("32_impute_strand")
    inv = Inventory("32_impute_strand")

    rows = []
    with log_time(logger, "strand imputation" + (" (dry-run)" if args.dry_run else "")):
        for path, platform, species, group, mod_type, tool, sample in iter_callsets():
            if args.sample and sample not in args.sample:
                continue
            df = read_table(path)
            if df.empty:
                continue
            cur = df["strand"].astype(str) if "strand" in df.columns else pd.Series([""] * len(df))
            unknown_before = int((~cur.isin(["+", "-"])).sum())
            rec = {"platform": platform, "species": species, "dataset_group": group,
                   "mod_type": mod_type, "tool": tool, "sample": sample,
                   "n_rows": len(df), "n_unknown_before": unknown_before,
                   "n_imputed_plus": 0, "n_imputed_minus": 0,
                   "n_ambiguous": 0, "n_no_exon": 0, "n_unknown_after": unknown_before,
                   "pct_unknown_after": round(100 * unknown_before / len(df), 4) if len(df) else 0.0,
                   "strand_source": "none"}
            if unknown_before:
                aidx, alabel = alignment_index_for(tool, sample)
                new_strand, status, stats = impute_frame(df, species, alignment_index=aidx)
                rec["strand_source"] = alabel if aidx else "exon"
                todo = status >= 0
                rec["n_imputed_plus"] = int((status == PLUS).sum())
                rec["n_imputed_minus"] = int((status == MINUS).sum())
                rec["n_ambiguous"] = int((status == BOTH).sum())
                rec["n_no_exon"] = int((status == NO_EXON).sum())
                rec.update({f"read_{k}": v for k, v in stats.items()})
                after = int((~pd.Series(new_strand).isin(["+", "-"])).sum())
                rec["n_unknown_after"] = after
                rec["pct_unknown_after"] = round(100 * after / len(df), 4)
                if not args.dry_run:
                    label = "imputed_read_bed" if aidx else "imputed_exon"
                    src = np.where(todo, label, "tool").astype(object)
                    src[todo & (status != PLUS) & (status != MINUS)] = "unresolved"
                    df["strand"] = new_strand
                    df["strand_src"] = src
                    write_table(df, path)
                logger.info("  %-22s %-26s unknown %7d -> +%-7d -%-7d amb%-6d noexon%-7d after %d",
                            sample, tool, unknown_before, rec["n_imputed_plus"],
                            rec["n_imputed_minus"], rec["n_ambiguous"], rec["n_no_exon"],
                            rec["n_unknown_after"])
            rows.append(rec)

    sdf = pd.DataFrame(rows)
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    out = MANIFEST_DIR / "strand_imputation.csv"
    write_table(sdf, out)
    inv.record(out, n_rows=len(sdf))
    inv.flush()

    if len(sdf):
        tot = int(sdf["n_unknown_before"].sum())
        after = int(sdf["n_unknown_after"].sum())
        logger.info("unknown strand rows: %d -> %d (%.2f%% of %d rows)",
                    tot, after, 100 * after / max(int(sdf["n_rows"].sum()), 1),
                    int(sdf["n_rows"].sum()))
        by_tool = (sdf.groupby("tool")[["n_unknown_before", "n_imputed_plus",
                                        "n_imputed_minus", "n_ambiguous", "n_no_exon",
                                        "n_unknown_after"]].sum())
        logger.info("per tool:\n%s", by_tool.to_string())
    logger.info("-> %s", out)
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
