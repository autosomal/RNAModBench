#!/usr/bin/env python3
"""01b -- fill the missing per-replicate callsets with R2Dtool liftover.

Many tools only had rep3 (or rep1) converted by the legacy pipeline even though
their **transcript-space raw output exists for every replicate** (e.g. CHEUI
``site_level_m6A_predictions.txt`` for rep1/rep2, m6Anet ``data.site_proba.csv``,
DENA ``*.tsv``, MINES ``*.bed``, DRUMMER ``summary.txt``).  This script rebuilds
the legacy transcript-space input for those samples (byte-identical: verified by
``02b_validate_liftover.py``) and runs ``r2d liftover`` to obtain genomic
callsets, so that biological replicates can actually be compared.

The liftover itself is done by ``common.liftover`` (exon arithmetic from the
GTF) instead of the ``r2d`` binary: r2d reproduces the archived files for
Human/Mouse/E. coli but mis-maps ~46 % of the Arabidopsis rows, while the
internal mapper is row-for-row identical to the legacy files for every species
(``02b_validate_liftover.py``).

Written callsets carry ``parser = gtf_liftover`` and the source directory is kept
in ``source_file``; they land in the same per-replicate tree as the callsets
extracted by ``01_extract_callsets.py`` (which always takes precedence).

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/01b_liftover_missing.py [--sample S] [--tool T] [--force]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import (CALLSET_ROOT_EXTENDED, MANIFEST_DIR, RESULT_RNA002,
                           RNA002_TOOLS, SAMPLES, SAMPLES_BY_NAME, callset_dir,
                           in_scope)
from common.io_utils import write_table
from common.legacy_liftover import (BUILDERS, RAW_PROBE, SPECIES_GTF, build,
                                    needs_liftover, post_liftover_filter)
from common.liftover import liftover_like_r2d, load_model
from common.manifest import Inventory, log_time, setup_logger
from common.parsers import CANONICAL_COLUMNS, parse_callset, parse_liftover
from common.registry import sample_of_dir

WORK = MANIFEST_DIR.parent / "_liftover"

META_COLUMNS = ["platform", "sample", "replicate_tag", "species", "dataset_group",
                "condition_class", "mod_type", "tool", "source_file", "parser",
                "model_set"]


def _empty_canonical() -> pd.DataFrame:
    """0-row frame with the canonical callset columns (a genuine Total_Detected=0)."""
    return pd.DataFrame(columns=CANONICAL_COLUMNS)


def _assemble(sample, spec, source_dir: Path, df: pd.DataFrame,
              parser_id: str) -> pd.DataFrame:
    """Prepend the metadata columns to a (possibly empty) canonical frame."""
    values = {c: "" for c in META_COLUMNS}
    values.update({
        "platform": sample.platform, "sample": sample.canonical,
        "replicate_tag": sample.replicate_tag, "species": sample.species,
        "dataset_group": sample.dataset_group,
        "condition_class": sample.condition_class,
        "mod_type": spec.mod_type, "tool": spec.tool,
        "source_file": str(source_dir.relative_to(RESULT_RNA002)),
        "parser": parser_id,
    })
    meta = pd.DataFrame({c: [values[c]] * len(df) for c in META_COLUMNS})
    return pd.concat([meta, df.reset_index(drop=True)], axis=1)


#: parsers written by earlier versions of this script.  Callsets carrying one of
#: them are regenerated even without ``--force`` (the older engine, the ``r2d``
#: binary, mis-maps part of the Arabidopsis rows).
SUPERSEDED_PARSERS = {"r2d_liftover"}


def engine_of(path: Path) -> str:
    """``parser`` value of an existing callset ('' when unreadable)."""
    try:
        with path.open() as fh:
            header = fh.readline().rstrip("\n").split("\t")
            if "parser" not in header:
                return ""
            idx = header.index("parser")
            first = fh.readline().rstrip("\n").split("\t")
        return first[idx] if len(first) > idx else ""
    except OSError:
        return ""


def sample_dirs_for(tool: str, tool_subdir: str, sample_name: str) -> list[Path]:
    """Raw directories that belong to ``sample_name``.

    Most directories canonicalise by name; the group-level comparison directories
    (Curlcake Nanocompore ``nanocompore_resultN``, the E.coli IVT
    ``E_IVT_neg_result`` / ``E_IVT_neg``) are mapped through
    ``registry.sample_of_dir`` so 01/01b/13/14 always agree.
    """
    tool_dir = RESULT_RNA002 / tool_subdir
    if not tool_dir.is_dir():
        return []
    out = []
    for d in sorted(p for p in tool_dir.iterdir() if p.is_dir()):
        canonical, _ = sample_of_dir(tool, d.name)
        if canonical == sample_name:
            out.append(d)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", action="append")
    ap.add_argument("--tool", action="append")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    logger = setup_logger("01b_liftover_missing")
    inv = Inventory("01b_liftover_missing")
    WORK.mkdir(parents=True, exist_ok=True)

    rows = []
    with log_time(logger, "liftover fill"):
        for spec in RNA002_TOOLS:
            builder = BUILDERS.get(spec.tool)
            if builder is None:
                continue
            if args.tool and spec.tool not in args.tool:
                continue
            for sample in SAMPLES:
                if sample.platform != "RNA002":
                    continue
                if args.sample and sample.canonical not in args.sample:
                    continue
                if not in_scope(sample.platform, sample.species, sample.dataset_group,
                                spec.mod_type):
                    # archived combination (no figure, no reference): replicates
                    # are deliberately not filled - see 11_scope_split.py
                    continue
                out_path = callset_dir(sample, spec.mod_type, spec.tool) / \
                    f"{sample.canonical}.tsv"
                if out_path.exists() and not args.force:
                    if engine_of(out_path) not in SUPERSEDED_PARSERS:
                        # the original r2d pipeline mis-mapped part of the
                        # Arabidopsis rows (see common/liftover.py), so any
                        # Arabidopsis callset derived from a legacy r2d liftover is
                        # regenerated from raw with the internal mapper.
                        if sample.species != "Arabidopsis":
                            continue
                        logger.info("  regenerating %-20s %-14s (Arabidopsis: "
                                    "re-derive from raw)", sample.canonical, spec.tool)
                    logger.info("  refilling %-20s %-14s (superseded engine)",
                                sample.canonical, spec.tool)
                use_liftover = needs_liftover(spec.tool, sample.species)
                gtf = SPECIES_GTF.get(sample.species)
                if use_liftover and (gtf is None or not gtf.exists()):
                    continue
                dirs = sample_dirs_for(spec.tool, spec.result_subdir, sample.canonical)
                if not dirs:
                    continue
                made = False
                for d in dirs:
                    # Does the tool's raw input actually exist in this directory?
                    # (independent of how many sites survive the tool threshold)
                    probe = RAW_PROBE.get(spec.tool)
                    raw = None
                    if probe is not None:
                        try:
                            raw = probe(d, spec.mod_type)
                        except Exception:  # noqa: BLE001
                            raw = None
                    try:
                        inp = build(spec.tool, d, spec.mod_type, sample.species)
                    except Exception as exc:  # noqa: BLE001
                        logger.warning("builder failed %s %s: %s", sample.canonical,
                                       spec.tool, exc)
                        continue

                    # ---- genuine 0-site: raw present, nothing passes the tool
                    # threshold (e.g. DRUMMER Arabidopsis_fip37_rep2, max
                    # |frac_diff|=0.075 < 0.1).  Write a header-only callset so the
                    # replicate is represented as Total_Detected=0, never silently
                    # dropped nor mislabelled "no raw input".
                    if inp is None or inp.empty:
                        if raw is None:
                            continue          # raw truly absent -> try another dir
                        out = _assemble(sample, spec, d, _empty_canonical(), "zero_sites")
                        n = write_table(out, out_path)
                        rows.append({
                            "sample": sample.canonical, "tool": spec.tool,
                            "mod_type": spec.mod_type, "dir": d.name,
                            "gtf": gtf.name if gtf else "", "n_rows": 0,
                            "status": "ok_zero",
                            "out_file": str(out_path.relative_to(MANIFEST_DIR.parent)),
                        })
                        logger.info("  zero   %-22s %-14s Total_Detected=0 (raw present, "
                                    "0 sites after filter; dir=%s)", sample.canonical,
                                    spec.tool, d.name)
                        inv.record(out_path, n_rows=0, source=str(d))
                        made = True
                        break

                    in_path = WORK / f"{sample.canonical}__{spec.tool}__input.txt"
                    out_txt = WORK / f"{sample.canonical}__{spec.tool}__liftover.txt"
                    inp.to_csv(in_path, sep="\t", index=False)
                    if use_liftover:
                        # internal exon-arithmetic liftover (``common.liftover``):
                        # reproduces the legacy ``*_liftover.txt`` files exactly,
                        # while the r2d binary mis-maps ~46% of the Arabidopsis
                        # rows (see the module docstring).
                        model = load_model(gtf)
                        lifted = liftover_like_r2d(inp, model)
                        lifted.to_csv(out_txt, sep="\t", index=False)
                        df = post_liftover_filter(spec.tool, sample.species,
                                                  parse_liftover(out_txt))
                        parser_id = "gtf_liftover"
                    else:
                        # already in genome (construct) coordinates - e.g. the
                        # Curlcake Nanocompore runs on ``cc.fasta``, or the genomic
                        # tools (NanoSPA / EpiNano_Error)
                        df = parse_callset(in_path, "standard")
                        parser_id = "legacy_rebuild"
                    if df.empty:
                        # liftover/parse produced nothing although raw existed:
                        # treat as a genuine 0-site result too.
                        df = _empty_canonical()
                        parser_id = "zero_sites"
                    out = _assemble(sample, spec, d, df, parser_id)
                    n = write_table(out, out_path)
                    status = "ok_zero" if n == 0 else "filled"
                    rows.append({
                        "sample": sample.canonical, "tool": spec.tool,
                        "mod_type": spec.mod_type, "dir": d.name,
                        "gtf": gtf.name if gtf else "", "n_rows": n,
                        "status": status,
                        "out_file": str(out_path.relative_to(MANIFEST_DIR.parent)),
                    })
                    logger.info("  %-6s %-22s %-14s rows=%-6d (dir=%s)", status,
                                sample.canonical, spec.tool, n, d.name)
                    inv.record(out_path, n_rows=n, source=str(d))
                    made = True
                    break
                if not made:
                    logger.info("  absent %-22s %-14s (no raw input in %s)",
                                sample.canonical, spec.tool,
                                [d.name for d in dirs] or "no dir")

    if rows:
        df = pd.DataFrame(rows)
        write_table(df, MANIFEST_DIR / "liftover_fill.csv")
        inv.record(MANIFEST_DIR / "liftover_fill.csv", n_rows=len(df))
        logger.info("filled %d callsets", len(df))
    inv.flush()
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
