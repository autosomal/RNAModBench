#!/usr/bin/env python3
"""13 -- completeness audit: for every (sample, tool) say what happened.

The registry (``00``) discovers *converted* files under ``result/``; ``01b``
fills missing replicates from *raw* tool output and records genuine
``Total_Detected = 0`` results as header-only callsets; ``11`` deletes
out-of-scope non-m6A for Arabidopsis / Mouse / E. coli.  This script overlays
those three views onto what is actually on disk and classifies every pair:

``filled``               callset present, >0 sites (provenance in ``how``)
``ok_zero``              header-only callset: raw present, tool threshold gave 0
                         sites (a real result, NOT a missing file)
``raw_present_unfilled`` NO callset but the tool's raw input exists and the tool
                         has a builder -> ``01b`` should have filled it (ANOMALY)
``raw_missing``          no callset, no raw, no converted file -> the tool was
                         simply never run on this sample (expected)
``raw_available_legacy_unused``
                         a converted/raw result EXISTS but the legacy assembly
                         deliberately never used it (see
                         ``legacy_liftover.LEGACY_UNUSED_RAW``) -- recorded so it
                         is never mistaken for a silent gap
``out_of_scope_deleted`` non-m6A on Arabidopsis/Mouse/E.coli -> physically
                         removed by ``11`` (verify it is really gone)
``out_of_scope_combination``
                         a manuscript tool on a (sample, modification) combination
                         the manuscript never used (DENA / MINES on E. coli, the
                         non-m6A panel outside HeLa & Curlcake, ...)
                         -- the tool IS part of the paper, the combination is not
``out_of_scope_tool``    tool outside the manuscript's 15-tool set (differr /
                         EpiNano_SVM / Tombo_com / CHEUI-diff / mAFiA / CHEUI on
                         Curlcake) -> physically removed by ``11``
``out_of_scope_archived`` non-m6A parked in ``callsets_extended/`` (HeLa/Curlcake
                         RNA004 Psi etc.)

Output: harmonisation/manifest/completeness_audit.csv (+ console summary).

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/13_completeness_audit.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import (ARTICLE_M6A_TOOLS, ARTICLE_NONM6A_TOOLS, CALLSET_ROOT,
                           CALLSET_ROOT_EXTENDED, DORADO_TOOL_PREFIX, MANIFEST_DIR,
                           PROJECT, RESULT_RNA002, RNA002_TOOLS, SAMPLES,
                           callset_dir, canonical_mod_type, delete_nonm6a,
                           in_scope, tool_in_scope)
from common.io_utils import read_table, write_table
from common.legacy_liftover import LEGACY_UNUSED_RAW, RAW_PROBE
from common.manifest import Inventory, log_time, setup_logger
from common.registry import sample_of_dir

REGISTRY = MANIFEST_DIR / "sample_tool_registry.csv"

#: parsers written by 01b (rebuilt from raw) vs 01_extract (converted file).
_RAW_PARSERS = {"gtf_liftover", "legacy_rebuild", "zero_sites", "r2d_liftover"}


def _row_count(path: Path) -> int:
    """Data-row count of a callset (header excluded); callsets are small."""
    try:
        with path.open("rb") as fh:
            return max(fh.read().count(b"\n") - 1, 0)
    except OSError:
        return -1


def _parser_of(path: Path) -> str:
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


def _raw_probe(spec, sample) -> Path | None:
    """The tool's raw input for this sample, via RAW_PROBE (builder tools only)."""
    probe = RAW_PROBE.get(spec.tool)
    if probe is None:
        return None
    tdir = RESULT_RNA002 / spec.result_subdir
    if not tdir.is_dir():
        return None
    for d in sorted(p for p in tdir.iterdir() if p.is_dir()):
        canonical, _ = sample_of_dir(spec.tool, d.name)
        if canonical == sample.canonical:
            try:
                r = probe(d, spec.mod_type)
            except Exception:  # noqa: BLE001
                r = None
            if r is not None:
                return r
    return None


def _archived(spec, sample) -> Path | None:
    mod = canonical_mod_type(spec.mod_type)
    p = (CALLSET_ROOT_EXTENDED / sample.platform / sample.species /
         sample.dataset_group / mod / spec.tool / f"{sample.canonical}.tsv")
    return p if p.is_file() else None


def _sync_manifests(df: pd.DataFrame, logger, inv) -> None:
    """Reconcile ``callsets_summary.csv`` / ``sample_tool_registry.csv`` with disk.

    ``00``/``01`` write these ledgers *before* ``01b`` raw-fills the missing
    replicates, so pairs that 01b created (e.g. the header-only DRUMMER
    ``Arabidopsis_fip37_rep2``) still read ``status=missing`` there.  Running
    last, this step stamps the on-disk truth back: status ok / ok_zero, the
    real row count, and (for ok_zero) *why* the total detection count is 0.
    """
    disk = {(r["sample"], r["tool"], r["mod_type"]): r for r in df.to_dict("records")
            if r["status"] in ("filled", "ok_zero")}

    # ---- callsets_summary.csv ----
    SUM = MANIFEST_DIR / "callsets_summary.csv"
    if SUM.exists():
        s = read_table(SUM)
        if "n_sites" not in s.columns:
            s["n_sites"] = pd.NA
        # pandas 3 str-dtype columns reject int setitem -> coerce up front
        s["rows_out"] = pd.to_numeric(s["rows_out"], errors="coerce").astype("Int64")
        s["n_sites"] = pd.to_numeric(s["n_sites"], errors="coerce").astype("Int64")
        s["is_raw"] = pd.to_numeric(s["is_raw"], errors="coerce").astype("Int64")
        hit = s["sample"].isin({k[0] for k in disk}) & s["tool"].isin({k[1] for k in disk})
        changed = 0
        for i, r in s[hit].iterrows():
            d = disk.get((r["sample"], r["tool"], r["mod_type"]))
            if d is None:
                continue
            new_status = "ok" if d["status"] == "filled" else "ok_zero"
            if r["status"] == new_status and int(r["rows_out"]) == d["n_sites"]:
                continue  # already consistent (01_extract converted-file row)
            note = ("raw present, 0 sites after tool threshold "
                    "(Total_Detected=0, filled by 01b)" if new_status == "ok_zero"
                    else "filled by 01b from raw")
            changed += 1
            s.at[i, "status"] = new_status
            s.at[i, "rows_out"] = d["n_sites"]
            s.at[i, "n_sites"] = d["n_sites"]
            s.at[i, "parser"] = "zero_sites" if new_status == "ok_zero" else "legacy_rebuild"
            s.at[i, "is_raw"] = 1
            s.at[i, "note"] = note
        write_table(s, SUM)
        inv.record(SUM, n_rows=len(s))
        logger.info("callsets_summary.csv reconciled: %d rows stamped from disk", changed)

    # ---- sample_tool_registry.csv ----
    REG = MANIFEST_DIR / "sample_tool_registry.csv"
    if REG.exists():
        g = read_table(REG)
        g["n_rows"] = pd.to_numeric(g["n_rows"], errors="coerce").astype("Int64")
        hit = g["sample"].isin({k[0] for k in disk}) & g["tool"].isin({k[1] for k in disk})
        changed = 0
        for i, r in g[hit].iterrows():
            d = disk.get((r["sample"], r["tool"], r["mod_type"]))
            if d is None:
                continue
            new_status = "ok" if d["status"] == "filled" else "ok_zero"
            if r["status"] != new_status:
                changed += 1
            g.at[i, "status"] = new_status
            g.at[i, "n_rows"] = d["n_sites"]
            if new_status == "ok_zero":
                g.at[i, "note"] = "raw present, Total_Detected=0 after tool threshold (01b)"
        write_table(g, REG)
        inv.record(REG, n_rows=len(g))
        logger.info("sample_tool_registry.csv reconciled: %d status flips", changed)


def main() -> None:
    logger = setup_logger("13_completeness_audit")
    inv = Inventory("13_completeness_audit")

    reg = read_table(REGISTRY) if REGISTRY.exists() else pd.DataFrame()
    reg_key = {}
    if len(reg):
        for _, r in reg.iterrows():
            reg_key[(r["sample"], r["tool"])] = r

    subdirs_present = {d.name for d in RESULT_RNA002.iterdir() if d.is_dir()}
    rows = []
    with log_time(logger, "completeness audit"):
        for sample in SAMPLES:
            if sample.platform != "RNA002":
                continue
            for spec in RNA002_TOOLS:
                mod = canonical_mod_type(spec.mod_type)
                cpath = callset_dir(sample, mod, spec.tool) / f"{sample.canonical}.tsv"
                exists = cpath.is_file()
                n = _row_count(cpath) if exists else 0
                parser = _parser_of(cpath) if exists else ""
                rr = reg_key.get((sample.canonical, spec.tool))
                reg_status = rr["status"] if rr is not None else ""
                src = (rr["source_file"] if rr is not None else "") or ""
                scoped = in_scope(sample.platform, sample.species,
                                  sample.dataset_group, mod)
                tool_ok = tool_in_scope(sample.platform, sample.species,
                                        sample.dataset_group, mod, spec.tool)

                if exists:
                    status = "filled" if n > 0 else "ok_zero"
                    how = "01b_from_raw" if parser in _RAW_PARSERS else "01_converted"
                elif delete_nonm6a(sample.species, mod):
                    arch = _archived(spec, sample)
                    status = "out_of_scope_deleted" if arch is None else "out_of_scope_archived"
                    how = "removed_by_11" if status == "out_of_scope_deleted" else "archived"
                    cpath = arch or cpath
                elif not tool_ok:
                    #: Two genuinely different facts, previously collapsed into one
                    #: label (which made the console line list DENA / DRUMMER /
                    #: MINES as "not a manuscript tool" -- they are, they just were
                    #: never run on *this* combination).
                    article_tool = (spec.tool in ARTICLE_M6A_TOOLS
                                    or spec.tool in ARTICLE_NONM6A_TOOLS
                                    or str(spec.tool).startswith(DORADO_TOOL_PREFIX))
                    status = "out_of_scope_combination" if article_tool else "out_of_scope_tool"
                    how = ("combination_not_in_manuscript" if article_tool
                           else "removed_by_11")
                elif not scoped:
                    arch = _archived(spec, sample)
                    status = "out_of_scope_archived" if arch is not None else "raw_missing"
                    how = "archived_by_11" if arch is not None else "no_raw"
                else:
                    raw = _raw_probe(spec, sample)
                    unused = LEGACY_UNUSED_RAW.get((sample.canonical, spec.tool))
                    if raw is not None and spec.tool in RAW_PROBE:
                        status, how = "raw_present_unfilled", str(raw)
                    elif reg_status in ("empty", "pending", "ok_zero"):
                        # ``empty``/``pending``: 00 discovered a converted file
                        # with 0 rows; ``ok_zero``: _sync_manifests (this same
                        # script, earlier run) already stamped the pair.  Either
                        # way the converted result IS the genuine
                        # Total_Detected = 0 -- never "missing".
                        status, how = "ok_zero", "converted_but_empty"
                    elif unused is not None and (PROJECT / unused).exists():
                        # a real result the legacy assembly deliberately skipped:
                        # recorded, never a silent gap
                        status, how = "raw_available_legacy_unused", unused
                    else:
                        status, how = "raw_missing", ""

                rows.append({
                    "platform": sample.platform, "species": sample.species,
                    "dataset_group": sample.dataset_group, "sample": sample.canonical,
                    "replicate_tag": sample.replicate_tag, "condition_class": sample.condition_class,
                    "mod_type": mod, "tool": spec.tool,
                    "in_scope": int(scoped and tool_ok),
                    "status": status, "n_sites": (n if exists else 0),
                    "provenance": how,
                    "callset": (str(cpath.relative_to(CALLSET_ROOT.parent))
                                if (exists or status == "out_of_scope_archived") else ""),
                    "result_subdir": spec.result_subdir, "registry_status": reg_status,
                })

    df = pd.DataFrame(rows)
    write_table(df, MANIFEST_DIR / "completeness_audit.csv")
    inv.record(MANIFEST_DIR / "completeness_audit.csv", n_rows=len(df))
    inv.flush()

    _sync_manifests(df, logger, inv)

    # ---------------- summary ----------------
    logger.info("completeness audit: %d (sample,tool) pairs", len(df))
    piv = (df.pivot_table(index=["species", "mod_type"], columns="status",
                          values="sample", aggfunc="size", fill_value=0, observed=True))
    logger.info("counts per species x mod_type x status:\n%s", piv.to_string())

    m6a = df[(df["mod_type"] == "m6A") & (df["in_scope"] == 1)]
    gap = m6a[m6a["status"] == "raw_present_unfilled"]
    zero = m6a[m6a["status"] == "ok_zero"]
    logger.info("m6A in-scope: filled=%d ok_zero=%d raw_missing=%d UNFILLED-GAP=%d",
                int((m6a['status'] == 'filled').sum()), len(zero),
                int((m6a['status'] == 'raw_missing').sum()), len(gap))
    if len(gap):
        for _, r in gap.iterrows():
            logger.warning("UNFILLED m6A gap: %s / %s (raw=%s)", r["sample"],
                           r["tool"], r["provenance"])
    if len(zero):
        logger.info("genuine 0-site (Total_Detected=0) samples:")
        for _, r in zero.iterrows():
            logger.info("   %-22s %-14s (%s)", r["sample"], r["tool"], r["provenance"])
    unused = df[df["status"] == "raw_available_legacy_unused"]
    if len(unused):
        logger.info("results the legacy assembly never used (%d):", len(unused))
        for _, r in unused.iterrows():
            logger.info("   %-22s %-14s (%s)", r["sample"], r["tool"], r["provenance"])
    oos = df[df["status"] == "out_of_scope_tool"]
    if len(oos):
        logger.info("tools the manuscript never used (removed by 11): "
                    "%d (sample, tool) pairs -> %s", len(oos),
                    sorted(oos["tool"].unique()))
    oosc = df[df["status"] == "out_of_scope_combination"]
    if len(oosc):
        logger.info("manuscript tools on a combination the paper did not use "
                    "(expected absence, NOT a missing tool): %d pairs -> %s",
                    len(oosc),
                    sorted({f"{r.dataset_group}/{r.mod_type}:{r.tool}"
                            for _, r in oosc.groupby(
                                ["dataset_group", "mod_type", "tool"]).size().reset_index()
                            .iterrows()}))
    logger.info("manifest -> %s", MANIFEST_DIR / "completeness_audit.csv")
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
