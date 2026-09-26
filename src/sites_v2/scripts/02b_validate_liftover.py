#!/usr/bin/env python3
"""02b -- validate the R2Dtool liftover reproduction on samples where the legacy
files already exist.

For every tool with a builder in ``common/legacy_liftover.py`` this script:

1. rebuilds the transcript-space input from the *raw* result file,
2. lifts it over with ``common.liftover`` (exon arithmetic from the species GTF),
3. compares the output with the legacy ``*_liftover.txt`` file row-for-row
   (as sets, the legacy files are not sorted).

For tools whose legacy step also applied a chromosome whitelist the comparison is
repeated against ``*_remove_chr.txt`` (``filtered_verdict``), and for the
construct-space Curlcake Nanocompore runs the builder output is compared directly
against ``*_Nanocompore.txt``.

Historical note: the archived pipeline used the ``r2d liftover`` binary.  It
reproduces the archived files for Human/Mouse/E. coli but mis-maps ~46 % of the
Arabidopsis rows in this workspace; the internal mapper reproduces all of them.

A tool/sample is only allowed into the replicate-fill step
(``01b_liftover_missing.py``) if the reproduction is exact here.

Output: sites_v2/manifest/liftover_validation.csv (+ console log)

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/sites_v2/scripts/02b_validate_liftover.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import (CONVERTED_CALLSETS, MANIFEST_DIR, RAW_OVER_CONVERTED,
                           RESULT_RNA002, RNA002_TOOLS, SAMPLES_BY_NAME)
from common.io_utils import write_table
from common.legacy_liftover import (BUILDERS, EXPLICIT_DIR_SAMPLE, GENOMIC_LEGACY_GLOB,
                                    SPECIES_GTF, build, needs_liftover,
                                    post_liftover_filter)
from common.liftover import liftover_like_r2d, load_model
from common.manifest import Inventory, setup_logger
from common.registry import canonical_from_dir

WORK = MANIFEST_DIR.parent / "_liftover"


def _rows(df: pd.DataFrame, cols: list[str] | None = None) -> set:
    d = df if cols is None else df[cols]
    return set(map(tuple, d.astype(str).values))


def _verdict(ka: set, kb: set, header_same: bool = True) -> tuple[str, int, int]:
    n_only_ours, n_only_theirs = len(ka - kb), len(kb - ka)
    if not (ka ^ kb):
        verdict = "identical" if header_same else "identical_values_header_diff"
    elif max(n_only_ours, n_only_theirs) < 0.005 * max(len(kb), 1):
        verdict = "nearly_identical"
    else:
        verdict = "different"
    return verdict, n_only_ours, n_only_theirs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", action="append", default=None,
                    help="only validate these tools (repeatable)")
    ap.add_argument("--out", default="liftover_validation.csv",
                    help="manifest file name to write")
    args = ap.parse_args()

    logger = setup_logger("02b_validate_liftover")
    inv = Inventory("02b_validate_liftover")
    WORK.mkdir(parents=True, exist_ok=True)

    rows = []
    by_tool = {t.tool: t for t in RNA002_TOOLS}
    for tool_id, builder in BUILDERS.items():
        spec = by_tool.get(tool_id)
        if spec is None:
            continue
        if args.tool and tool_id not in args.tool:
            continue
        tool_dir = RESULT_RNA002 / spec.result_subdir
        if not tool_dir.is_dir():
            continue
        for d in sorted(p for p in tool_dir.iterdir() if p.is_dir()):
            legacy = sorted(d.glob(f"*_{tool_id}_liftover.txt"))
            # for a genomic (no-liftover) tool the archived final file differs per
            # tool (Nanocompore ``*_Nanocompore.txt``, NanoSPA ``*_5mer.txt``,
            # EpiNano_Error ``*_remove_chr.txt``).
            genomic_glob = GENOMIC_LEGACY_GLOB.get(tool_id, "*_Nanocompore.txt")
            legacy_txt = sorted(d.glob(genomic_glob))
            sample, _ = canonical_from_dir(d.name)
            if sample is None:
                sample = EXPLICIT_DIR_SAMPLE.get(tool_id, {}).get(d.name)
            if sample is None or sample not in SAMPLES_BY_NAME:
                continue
            superseded = (sample, tool_id) in RAW_OVER_CONVERTED
            #: the 2026-09-15 reorganization promoted the archival copies out of
            #: ``result/`` into the converted layer, so fall back to it for *every*
            #: pair -- without this the comparison silently does not happen at all
            #: (02b produced no rows between 2026-09-15 and 2026-09-18).
            if not legacy or not legacy_txt:
                alt = CONVERTED_CALLSETS / spec.result_subdir / d.name
                if alt.is_dir():
                    if not legacy:
                        legacy = sorted(alt.glob(f"*_{tool_id}_liftover.txt"))
                    if not legacy_txt:
                        legacy_txt = sorted(alt.glob(genomic_glob))
            species = SAMPLES_BY_NAME[sample].species
            gtf = SPECIES_GTF.get(species)
            if needs_liftover(tool_id, species):
                if not legacy or gtf is None or not gtf.exists():
                    continue
            elif not legacy_txt:
                continue
            try:
                inp = build(tool_id, d, spec.mod_type, species)
            except Exception as exc:  # noqa: BLE001
                logger.warning("builder failed: %s %s: %s", sample, tool_id, exc)
                continue
            if inp is None or inp.empty:
                continue
            rec = {"sample": sample, "species": species, "tool": tool_id,
                   "gtf": gtf.name if gtf else "",
                   "legacy_file": (legacy or legacy_txt)[0].name}
            in_path = WORK / f"{sample}__{tool_id}__input.txt"
            inp.to_csv(in_path, sep="\t", index=False)

            if not needs_liftover(tool_id, species):
                # construct-space rebuild: compare the builder output directly
                theirs = pd.read_csv(legacy_txt[0], sep="\t", dtype=str, keep_default_na=False)
                ka, kb = _rows(inp), _rows(theirs)
                verdict, n_ours, n_theirs = _verdict(
                    ka, kb, list(inp.columns) == list(theirs.columns))
                rec.update({"n_legacy": len(theirs), "n_rebuilt": len(inp),
                            "only_rebuilt": n_ours, "only_legacy": n_theirs,
                            "verdict": verdict, "filtered_verdict": "",
                            "n_filtered": "", "n_legacy_remove_chr": "",
                            "legacy_status": "superseded"
                            if (sample, tool_id) in RAW_OVER_CONVERTED else ""})
                rows.append(rec)
                if (sample, tool_id) in RAW_OVER_CONVERTED:
                    #: archival file known-defective and deliberately superseded
                    #: (2026-09-18: the Arabidopsis EpiNano per-site tables were
                    #: built against a human FASTA -- see code/epinano_refix/);
                    #: the rebuild is the truth, so this is not a regression.
                    logger.warning("%-22s %-14s legacy file SUPERSEDED (wrong-reference "
                                   "per-site table): legacy=%-7d rebuilt=%-7d "
                                   "(verdict=%s recorded for the record)",
                                   sample, tool_id, len(theirs), len(inp), verdict)
                else:
                    logger.info("%-22s %-14s legacy=%-7d rebuilt=%-7d verdict=%s (no liftover)",
                                sample, tool_id, len(theirs), len(inp), verdict)
                continue

            out_path = WORK / f"{sample}__{tool_id}__liftover.txt"
            ours = liftover_like_r2d(inp, load_model(gtf))
            ours.to_csv(out_path, sep="\t", index=False)
            theirs = pd.read_csv(legacy[0], sep="\t", dtype=str, keep_default_na=False)
            ours = ours.astype(str)
            # the two frames may use different score-column names; compare the
            # genomic block by position and the shared columns by name
            shared = [c for c in theirs.columns if c in ours.columns]
            ours = ours[shared]
            theirs = theirs[shared]
            if ours.shape[1] != theirs.shape[1]:
                verdict, n_only_ours, n_only_theirs = "shape_mismatch", None, None
            else:
                verdict, n_only_ours, n_only_theirs = _verdict(
                    _rows(ours), _rows(theirs),
                    list(ours.columns) == list(theirs.columns))
            rec.update({"n_legacy": len(theirs), "n_rebuilt": len(ours),
                        "only_rebuilt": n_only_ours, "only_legacy": n_only_theirs,
                        "verdict": verdict, "filtered_verdict": "",
                        "n_filtered": "", "n_legacy_remove_chr": ""})

            # end-to-end check of the chromosome whitelist against *_remove_chr.txt
            rc = sorted(d.glob("*_remove_chr.txt"))
            if rc:
                kept = post_liftover_filter(tool_id, species, ours)
                cols = ["chromosome", "start", "end", "Status", "Pvalue", "strand"]
                cols = [c for c in cols if c in kept.columns]
                f_theirs = pd.read_csv(rc[0], sep="\t", dtype=str, keep_default_na=False)
                fverdict, _, _ = _verdict(_rows(kept, cols), _rows(f_theirs))
                rec.update({"filtered_verdict": fverdict, "n_filtered": len(kept),
                            "n_legacy_remove_chr": len(f_theirs)})
            rows.append(rec)
            logger.info("%-22s %-14s legacy=%-7d rebuilt=%-7d verdict=%-18s filtered=%s",
                        sample, tool_id, len(theirs), len(ours), verdict,
                        rec["filtered_verdict"] or "-")

    if rows:
        df = pd.DataFrame(rows)
        out_csv = MANIFEST_DIR / args.out
        write_table(df, out_csv)
        inv.record(out_csv, n_rows=len(df))
        logger.info("verdict counts: %s", df["verdict"].value_counts().to_dict())
        if "filtered_verdict" in df.columns:
            logger.info("filtered verdict counts: %s",
                        df["filtered_verdict"].replace("", pd.NA).dropna()
                          .value_counts().to_dict())
        ok = {t for t, g in df.groupby("tool") if (g["verdict"] == "identical").all()}
        logger.info("tools with 100%% identical reproduction: %s", sorted(ok))
    inv.flush()
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
