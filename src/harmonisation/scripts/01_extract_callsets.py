#!/usr/bin/env python3
"""01 -- extract per-sample (== per-replicate) callsets from result/ + result_RNA004/.

Reads the registry produced by ``00_build_registry.py`` and writes one TSV per
(sample, tool, mod_type) pair under::

    harmonisation/callsets/<platform>/<species>/<dataset_group>/<mod_type>/<tool>/<sample>.tsv

Every row keeps the tool's raw coordinate (``pos_raw``, 0-based, base
conversion only).  Annotation columns (5mer / DRACH / coverage / GLORI
distance) are added later by ``03_annotate_callsets.py``; the offset audit
(``05_offset_audit.py``) reports systematic shifts per tool.

A summary of every written callset (rows in/out, parser, filters, source file
fingerprint) goes to ``harmonisation/manifest/callsets_summary.csv``.

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/01_extract_callsets.py [--sample S] [--tool T]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import (MANIFEST_DIR, SAMPLES_BY_NAME, callset_dir,
                           canonical_mod_type, tool_in_scope)
from common.io_utils import ensure_dirs, file_fingerprint, rel, write_table
from common.legacy_liftover import dedupe_epinano_sites
from common.manifest import Inventory, log_time, setup_logger
from common.parsers import CANONICAL_COLUMNS, parse_callset
from common.registry import build_registry, parse_dorado_name

META_COLUMNS = ["platform", "sample", "replicate_tag", "species", "dataset_group",
                "condition_class", "mod_type", "tool", "source_file", "parser",
                "model_set"]


def extract_one(row: pd.Series, logger) -> dict:
    sample = SAMPLES_BY_NAME[row["sample"]]
    src = Path(row["source_file"])
    parser = row["parser"] or "standard"
    #: ASCII-normalise the modification id so ``Psi`` never splits into a second
    #: ``\u03a8`` directory (see ``config.canonical_mod_type``).
    mod_type = canonical_mod_type(row["mod_type"])
    out_path = (callset_dir(sample, mod_type, row["tool"]) /
                f"{sample.canonical}.tsv")

    rec = {
        "platform": row["platform"], "sample": row["sample"], "tool": row["tool"],
        "mod_type": mod_type, "species": row["species"],
        "dataset_group": row["dataset_group"],
        "source_file": rel(src) if row["source_file"] else "",
        "source_fingerprint": file_fingerprint(src) if row["source_file"] else "",
        "parser": parser, "is_raw": int(row["is_raw"]), "status": row["status"],
        "rows_in": "", "rows_out": "", "note": row["note"],
        "out_file": rel(out_path),
    }
    # ``empty`` sources (a 0-byte or header-only converted callset) are parsed too:
    # they are genuine 0-detection results and must leave a header-only file.
    if row["status"] not in ("ok", "empty"):
        rec["rows_in"] = row["n_rows"]
        rec["rows_out"] = ""
        return rec

    try:
        kwargs: dict = {}
        if parser == "m6anet_csv":
            kwargs["platform"] = row["platform"]
        #: ``parser_arg`` (JSON) carries the per-callset parser options the registry
        #: derived from the source itself -- e.g. ``{"mod_code": "a"}`` for a Dorado
        #: pileup that mixes several modifications in one file.
        parser_arg = str(row.get("parser_arg", "") or "").strip()
        if parser_arg:
            try:
                kwargs.update(json.loads(parser_arg))
            except ValueError:
                logger.warning("unreadable parser_arg for %s x %s: %r",
                               row["sample"], row["tool"], parser_arg)
            else:
                rec["note"] = "|".join(x for x in (rec["note"],
                                                   f"parser_arg={parser_arg}") if x)
        try:
            df = parse_callset(src, parser, **kwargs)
        except pd.errors.EmptyDataError:
            # a 0-byte / newline-only source is a real 0-detection result (e.g.
            # the Curlcake ELIGOS2 combine files), not a broken parse
            df = pd.DataFrame(columns=CANONICAL_COLUMNS)
    except Exception as exc:  # noqa: BLE001 - recorded, not swallowed
        logger.error("parse failed: %s x %s (%s): %s", row["sample"], row["tool"],
                     src, exc)
        rec["status"] = "parse_error"
        rec["note"] = f"{row['note']}|{exc}"
        return rec

    if row["tool"] == "EpiNano_Error" and len(df) and "chrom" in df.columns:
        #: the legacy DiffErr prediction files were appended to by repeated R runs,
        #: so one site can appear up to 5x (identical values).  A callset is a *set*
        #: of sites -> collapse the copies and record how many rows that removed
        #: (2026-09-18; also applied inside build_epinano_error for the rebuilt path,
        #: keyed on the same species because the mouse chains have no per-site
        #: tables left to re-run from).
        n_before = len(df)
        df = dedupe_epinano_sites(df, subset=("chrom", "pos_raw", "strand"))
        if len(df) != n_before:
            rec["note"] = "|".join(x for x in (
                rec["note"], f"deduped {n_before - len(df)} repeated site rows") if x)

    if df.empty:
        # A real 0-detection result (the tool ran and its own threshold kept
        # nothing) still gets a header-only callset, so the (sample, tool) cell is
        # present and auditable instead of silently absent -- the same convention
        # ``01b`` uses for ``zero_sites``.
        rec["status"] = "empty_parsed"
        if not len(df.columns):
            df = pd.DataFrame(columns=CANONICAL_COLUMNS)
        values = {c: row.get(c, "") for c in META_COLUMNS}
        meta = pd.DataFrame({c: [values[c]] * 0 for c in META_COLUMNS})
        write_table(pd.concat([meta, df.reset_index(drop=True)], axis=1), out_path)
        rec["rows_out"] = 0
        return rec

    model_set = ""
    if parser == "dorado_pileup":
        _, _, model = parse_dorado_name(src.stem)
        model_set = model
    values = {c: row.get(c, "") for c in META_COLUMNS}
    values["model_set"] = model_set
    meta = pd.DataFrame({c: [values[c]] * len(df) for c in META_COLUMNS})
    out = pd.concat([meta, df.reset_index(drop=True)], axis=1)
    n = write_table(out, out_path)
    rec["rows_in"] = row["n_rows"]
    rec["rows_out"] = n
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", action="append", default=None)
    ap.add_argument("--tool", action="append", default=None)
    args = ap.parse_args()

    logger = setup_logger("01_extract_callsets")
    inv = Inventory("01_extract_callsets")

    _, pairs, _ = build_registry()
    if args.sample:
        pairs = pairs[pairs["sample"].isin(args.sample)]
    if args.tool:
        pairs = pairs[pairs["tool"].isin(args.tool)]

    logger.info("extracting %d (sample, tool) pairs", len(pairs))

    recs = []
    with log_time(logger, "extraction"):
        for _, row in pairs.iterrows():
            rec = extract_one(row, logger)
            recs.append(rec)
            if rec.get("status") not in ("ok", "empty"):
                logger.warning("  %-22s %-16s %s (%s)", rec["sample"], rec["tool"],
                               rec.get("status"), rec.get("note", ""))

    summary = pd.DataFrame(recs)
    ensure_dirs(MANIFEST_DIR)
    ledger = MANIFEST_DIR / "callsets_summary.csv"
    if (args.sample or args.tool) and ledger.exists():
        #: A filtered run must not shrink the ledger to its own subset -- it did,
        #: and the "full" callsets_summary.csv silently became 85 Curlcake rows.
        #: Keep every other (sample, tool) record untouched and replace only the
        #: pairs this run actually processed.
        prev = pd.read_csv(ledger, sep="\t", dtype=str, engine="python")
        done = set(zip(summary["sample"], summary["tool"]))
        keep = [(s, t) not in done
                for s, t in zip(prev["sample"].fillna(""), prev["tool"].fillna(""))]
        summary = pd.concat([prev[keep], summary], ignore_index=True)
        logger.info("filtered run: merged into the existing ledger -> %d rows total",
                    len(summary))
    write_table(summary, ledger)
    inv.record(ledger, n_rows=len(summary))
    inv.flush()

    ok = summary[summary["status"] == "ok"]
    logger.info("extracted callsets: %d ok / %d rows total", len(ok),
                int(pd.to_numeric(ok["rows_out"], errors="coerce").fillna(0).sum()))
    logger.info("by species:\n%s", ok.groupby(["platform", "species"]).size().to_string())
    bad = summary[summary["status"] != "ok"]
    logger.info("non-ok pairs: %d", len(bad))
    if len(bad):
        logger.info("by status: %s", bad["status"].value_counts().to_dict())
    logger.info("summary -> %s", MANIFEST_DIR / "callsets_summary.csv")
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
