#!/usr/bin/env python3
"""11 -- split the callsets into the manuscript scope vs the archived remainder.

Why
---
The non-m6A tools (m5C / Psi / m1Psi / Nm) carry the manuscript's **Fig. 7**,
which was produced on **HeLa (WT + IVT) and the unmodified Curlcake constructs
only** - the reviewer questions about non-m6A data (R2-2: "why an m6A-centred
GLORI benchmark can validly support performance claims for non-m6A tools";
review comment about "substantial false-positive detections in unmodified
Curlcake controls and poor reproducibility") are all anchored on those samples.

The tools were nevertheless run on other species too (CHEUI_m5C on Arabidopsis
/ mouse / E. coli, NanoMUD/NanoSPA/NanoPsu on Arabidopsis, ...).  Those callsets
are *not* part of any figure and have no reference, so on revision request
(2026-09-15) they are **physically deleted** (``shutil.rmtree``) and only
*recorded* in the deletion manifests.  Nothing is lost: every deleted callset is
reproducible from ``$RNAMODBENCH_LOCAL/raw/result/`` + the recipes in
``common/legacy_liftover.BUILDERS``, which is what R1-10 ("deposit all processed
callsets") actually needs.

What it does
------------
1. **Normalises the modification directory name** (``\u03a8`` -> ``Psi``), so the
   same modification never splits into two directories again.
2. **Deletes** every out-of-scope callset, in two independent passes:
   (a) the whole non-m6A ``mod_type`` directory of Arabidopsis / Mouse / E. coli;
   (b) any tool directory -- **including under ``m6A``** -- whose tool is not one
   of the 15 tools the manuscript evaluates (differr / EpiNano_SVM / Tombo_com /
   CHEUI-diff / mAFiA / CHEUI on Curlcake).  Emptied directories are pruned.
   ``--undo`` restores them from the deletion manifests.
3. Writes ``manifest/nonm6a_scope.csv`` (one row per callset file: scope, reason,
   action) + the two deletion ledgers, and adds a ``scope`` column to
   ``manifest/callsets_summary.csv``.

Scope rule
----------
* ``m6A``                                  -> in scope (all species, Fig. 2-6 + RNA004)
* non-m6A on ``HeLa_WT`` / ``HeLa_IVT`` /
  ``Curlcake_IVT`` (RNA002)                -> in scope (Fig. 7 panel)
* non-m6A on the three RNA004 groups in
  ``config.IN_SCOPE_NONM6A``               -> in scope (Fig. S9A)
* every other non-m6A combination          -> deleted (Arabidopsis/Mouse/E.coli)
* any tool outside the 15-tool set         -> deleted (whatever its mod_type)

Dorado is one of the 15 tools and its RNA004 pileups carry **four** modifications
(m6A / m5C / Psi / inosine, read from the modkit ``name`` code -- see
``docs/pipeline.md``), so ``config._DORADO_GROUPS`` lists all four per
sample.  Listing only m6A there made this script delete the correctly mod-typed
callsets that ``01`` had just written (2026-09-17).

The rule itself lives in ``common/config.py`` (``IN_SCOPE_NONM6A``,
``delete_nonm6a()``, ``in_scope()``, ``tool_in_scope()``) -- this script only
applies it to the tree.

Where it belongs in the pipeline
--------------------------------
**After** ``01b``/``03`` (``01`` does not know about scope and re-creates the
out-of-scope files on every clean rebuild, so 11 must run after it) and, since
2026-09-16, **before** ``05``/``04``/``06-09``: 05, 06 and 09 index the callset
tree with an unfiltered ``rglob``, so with the old order (11 after 09) their
tables kept rows for callsets that were deleted a moment later.  Running 11
first makes "no evaluation table can contain a deleted callset" structural
instead of a convention.

Outputs
-------
harmonisation/manifest/nonm6a_scope.csv
harmonisation/manifest/callsets_summary.csv   (``scope`` column added/refreshed)
harmonisation/manifest/nonm6a_deleted.csv
harmonisation/manifest/out_of_scope_tools_deleted.csv

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/11_scope_split.py --dry-run
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/11_scope_split.py
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/11_scope_split.py --undo
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import (CALLSET_ROOT, CALLSET_ROOT_EXTENDED, IN_SCOPE_NONM6A,
                           MANIFEST_DIR, MOD_LABEL, MOD_TYPE_ORDER, PROJECT,
                           SITES_ROOT, TABLE_DIR, canonical_mod_type, delete_nonm6a,
                           tool_in_scope)
from common.io_utils import ensure_dirs, read_table, write_table
from common.manifest import Inventory, log_time, setup_logger

MANIFEST = MANIFEST_DIR / "nonm6a_scope.csv"
DELETED_MANIFEST = MANIFEST_DIR / "nonm6a_deleted.csv"
DELETED_TOOLS_MANIFEST = MANIFEST_DIR / "out_of_scope_tools_deleted.csv"
SUMMARY = MANIFEST_DIR / "callsets_summary.csv"

#: The in-scope (platform, species, dataset_group) triples live in
#: ``config.IN_SCOPE_NONM6A`` so that ``01b_liftover_missing.py`` honours the
#: same rule when it creates new per-replicate callsets.

REASON_M6A = "m6A panel: all species / samples are part of the manuscript"
REASON_NONM6A_IN = ("non-m6A panel of the manuscript (Fig. 7 + R2-2): "
                    "HeLa (WT/IVT) and the unmodified Curlcake constructs")
REASON_NONM6A_OUT = ("non-m6A tool run outside the manuscript scope: no figure, "
                     "no reference; archived for traceability (R1-10)")
REASON_NONM6A_DELETE = ("non-m6A on Arabidopsis/Mouse/E.coli: out of the manuscript "
                        "scope (Fig. 7 is HeLa + Curlcake only), no reference; "
                        "DELETED on revision request (2026-09-14), listed in "
                        "manifest/nonm6a_deleted.csv for traceability")
REASON_TOOL_OUT = ("tool outside the manuscript's tool set: the paper evaluates 15 "
                   "dRNA-seq tools (manuscript.tex, Experimental Section) and this "
                   "tool is not among them (differr / EpiNano_SVM / Tombo_com / "
                   "CHEUI-diff / mAFiA / CHEUI on Curlcake); DELETED on revision request "
                   "(2026-09-15), listed in manifest/out_of_scope_tools_deleted.csv")


def iter_mod_dirs(root: Path):
    """Yield every ``<platform>/<species>/<group>/<mod_type>`` directory."""
    for platform in sorted(p for p in root.iterdir() if p.is_dir()):
        for species in sorted(p for p in platform.iterdir() if p.is_dir()):
            for group in sorted(p for p in species.iterdir() if p.is_dir()):
                for mod in sorted(p for p in group.iterdir() if p.is_dir()):
                    yield platform.name, species.name, group.name, mod


def normalize_mod_dirs(logger, dry_run: bool) -> int:
    """Rename non-ASCII modification directories (``\u03a8``) onto the ASCII id."""
    merged = 0
    touched: list[Path] = []
    for platform, species, group, mod in iter_mod_dirs(CALLSET_ROOT):
        canonical = canonical_mod_type(mod.name)
        if canonical == mod.name:
            continue
        target = mod.parent / canonical
        logger.warning("mod_type normalisation: %s -> %s",
                       mod.relative_to(SITES_ROOT), target.relative_to(SITES_ROOT))
        if dry_run:
            merged += 1
            continue
        target.mkdir(parents=True, exist_ok=True)
        for tool_dir in sorted(p for p in mod.iterdir() if p.is_dir()):
            dest = target / tool_dir.name
            if dest.exists():
                for f in sorted(tool_dir.glob("*.tsv")):
                    (dest / f.name).unlink(missing_ok=True)
                    shutil.move(str(f), str(dest / f.name))
                for leftover in sorted(tool_dir.glob("*")):
                    shutil.move(str(leftover), str(dest / leftover.name))
                tool_dir.rmdir()
            else:
                shutil.move(str(tool_dir), str(dest))
        for leftover in sorted(mod.iterdir()):
            shutil.move(str(leftover), str(target / leftover.name))
        mod.rmdir()
        touched.extend(sorted(target.rglob("*.tsv")))
        merged += 1
    if touched and not dry_run:
        fix_mod_column(touched, logger)
    return merged


def fix_mod_column(files: list[Path], logger) -> int:
    """Rewrite the in-file ``mod_type`` column so it matches the directory name."""
    n_fixed = 0
    for f in files:
        df = read_table(f)
        if "mod_type" not in df.columns:
            continue
        canonical = canonical_mod_type(df["mod_type"].iloc[0]) if len(df) else ""
        if not canonical or df["mod_type"].astype(str).eq(canonical).all():
            continue
        df["mod_type"] = canonical
        write_table(df, f)
        n_fixed += 1
    if n_fixed:
        logger.warning("in-file mod_type column normalised in %d callsets", n_fixed)
    return n_fixed


def callset_files(mod_dir: Path) -> list[Path]:
    return sorted(mod_dir.rglob("*.tsv"))


def normalize_mod_columns(logger, dry_run: bool) -> int:
    """Normalise the ``mod_type`` column of every callset whose dir was fixed.

    Probes only the header + first data row, so the multi-hundred-MB callsets are
    never loaded unless their label actually disagrees with their directory.
    """
    fixed = 0
    for root in (CALLSET_ROOT, CALLSET_ROOT_EXTENDED):
        if not root.exists():
            continue
        for f in sorted(root.rglob("*.tsv")):
            with f.open() as fh:
                header = fh.readline().rstrip("\n").split("\t")
                if "mod_type" not in header:
                    continue
                idx = header.index("mod_type")
                first = fh.readline().rstrip("\n").split("\t")
            if len(first) <= idx:
                continue
            canonical = canonical_mod_type(first[idx])
            if first[idx] == canonical:      # already canonical, nothing to do
                continue
            if dry_run:
                logger.warning("would normalise mod_type column: %s",
                               f.relative_to(SITES_ROOT))
                fixed += 1
                continue
            df = read_table(f)
            df["mod_type"] = canonical
            write_table(df, f)
            fixed += 1
    if fixed:
        logger.warning("in-file mod_type column normalised in %d callsets", fixed)
    return fixed


def scope_of(platform: str, species: str, group: str, mod_type: str) -> tuple[str, str]:
    if canonical_mod_type(mod_type) == "m6A":
        return "in_scope", REASON_M6A
    if (platform, species, group) in IN_SCOPE_NONM6A:
        return "in_scope", REASON_NONM6A_IN
    return "out_of_scope", REASON_NONM6A_OUT


def archived_rows() -> list[dict]:
    """Rows for everything already parked in ``callsets_extended/`` (idempotency)."""
    out: list[dict] = []
    if not CALLSET_ROOT_EXTENDED.exists():
        return out
    for f in sorted(CALLSET_ROOT_EXTENDED.rglob("*.tsv")):
        platform, species, group, mod_type, tool = f.relative_to(
            CALLSET_ROOT_EXTENDED).parts[:5]
        scope, reason = scope_of(platform, species, group,
                                 canonical_mod_type(mod_type))
        out.append({
            "scope": scope, "reason": reason,
            "platform": platform, "species": species, "dataset_group": group,
            "mod_type": canonical_mod_type(mod_type),
            "mod_display": MOD_LABEL.get(canonical_mod_type(mod_type), mod_type),
            "tool": tool, "sample": f.stem,
            "n_rows": sum(1 for _ in f.open()) - 1,
            "size_bytes": f.stat().st_size,
            "callsets_path": str(f.relative_to(CALLSET_ROOT_EXTENDED)),
            "archived_path": str(f.relative_to(SITES_ROOT)),
            "action": "moved",
        })
    return out


def move_file(src: Path, dest: Path, dry_run: bool) -> None:
    if dry_run:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        dest.unlink()
    shutil.move(str(src), str(dest))


def prune_empty(root: Path) -> None:
    for p in sorted((d for d in root.rglob("*") if d.is_dir()),
                    key=lambda p: -len(p.parts)):
        try:
            p.rmdir()
        except OSError:
            pass


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="report the planned moves without touching the filesystem")
    ap.add_argument("--undo", action="store_true",
                    help="move everything in callsets_extended/ back into callsets/")
    args = ap.parse_args()

    logger = setup_logger("11_scope_split")
    inv = Inventory("11_scope_split")
    ensure_dirs(MANIFEST_DIR, TABLE_DIR)

    # ---------------------------------------------------------------- undo ----
    if args.undo:
        if not MANIFEST.exists():
            raise SystemExit(f"no manifest to undo from: {MANIFEST}")
        man = read_table(MANIFEST)
        moved = man[man["action"] == "moved"]
        with log_time(logger, "undo scope split"):
            for _, r in moved.iterrows():
                src = SITES_ROOT / r["archived_path"]
                dest = SITES_ROOT / r["callsets_path"]
                if src.exists():
                    move_file(src, dest, args.dry_run)
            prune_empty(CALLSET_ROOT_EXTENDED)
        logger.info("restored %d callsets from %s", len(moved),
                    CALLSET_ROOT_EXTENDED.relative_to(SITES_ROOT))
        return

    rows: list[dict] = []
    deleted: list[dict] = []
    deleted_tools: list[dict] = []
    n_deleted_tool_files = 0

    def _del_row(f: Path, root: Path, platform, species, group, mod_type,
                 reason: str = REASON_NONM6A_DELETE) -> dict:
        return {
            "tree": root.name, "platform": platform, "species": species,
            "dataset_group": group, "mod_type": mod_type,
            "mod_display": MOD_LABEL.get(mod_type, mod_type),
            "tool": f.parent.name, "sample": f.stem,
            "n_rows": sum(1 for _ in f.open()) - 1,
            "size_bytes": f.stat().st_size,
            "deleted_path": str(f.relative_to(SITES_ROOT)),
            "reason": reason, "action": "deleted",
        }

    def _delete_mod_dir(mod: Path, root: Path, platform, species, group, mod_type) -> int:
        """Physically remove one out-of-scope non-m6A ``<mod_type>`` directory."""
        files = callset_files(mod)
        for f in files:
            deleted.append(_del_row(f, root, platform, species, group, mod_type))
        if not args.dry_run:
            shutil.rmtree(mod, ignore_errors=True)
        logger.warning("DELETE non-m6A %s/%s/%s/%s/%s (%d files)", root.name, platform,
                       species, group, mod_type, len(files))
        return len(files)

    def _delete_out_of_scope_tools(mod: Path, root: Path, platform, species, group,
                                   mod_type) -> int:
        """Physically remove the tool directories the manuscript never used.

        The paper evaluates 15 dRNA-seq tools (``config.ARTICLE_TOOL_SCOPE``);
        anything else found in ``result/`` (differr, EpiNano_SVM, Tombo_com,
        CHEUI-diff, mAFiA, CHEUI on Curlcake, non-m6A tools on
        Arabidopsis/Mouse/E. coli) belongs to no figure and is removed outright,
        empty directories included.
        """
        n = 0
        for tool_dir in sorted(p for p in mod.iterdir() if p.is_dir()):
            if tool_in_scope(platform, species, group, mod_type, tool_dir.name):
                continue
            files = callset_files(tool_dir)
            for f in files:
                deleted_tools.append(_del_row(f, root, platform, species, group,
                                              mod_type, REASON_TOOL_OUT))
            if not args.dry_run:
                shutil.rmtree(tool_dir, ignore_errors=True)
            logger.warning("DELETE out-of-scope tool %s/%s/%s/%s/%s/%s (%d files)",
                           root.name, platform, species, group, mod_type,
                           tool_dir.name, len(files))
            n += len(files)
        return n

    with log_time(logger, "scope split"):
        normalize_mod_dirs(logger, args.dry_run)
        normalize_mod_columns(logger, args.dry_run)

        for platform, species, group, mod in iter_mod_dirs(CALLSET_ROOT):
            mod_type = canonical_mod_type(mod.name)
            # (1) DELETE: Arabidopsis/Mouse/E.coli non-m6A -- removed outright
            #     (whole ``<mod_type>`` dir, empty tool subdirs included).
            if delete_nonm6a(species, mod_type):
                _delete_mod_dir(mod, CALLSET_ROOT, platform, species, group, mod_type)
                continue
            # (1b) DELETE: tools outside the manuscript's 15-tool set
            #      (config.tool_in_scope; differr / EpiNano_SVM / Tombo_com /
            #      CHEUI-diff / mAFiA / CHEUI on Curlcake).
            n_tool = _delete_out_of_scope_tools(mod, CALLSET_ROOT, platform, species,
                                                group, mod_type)
            if n_tool:
                n_deleted_tool_files += n_tool
            # (2) archive / keep everything else
            scope, reason = scope_of(platform, species, group, mod_type)
            files = callset_files(mod)
            if not files:
                continue
            for f in files:
                rel = f.relative_to(CALLSET_ROOT)
                tool = f.parent.name
                sample = f.stem
                archived = CALLSET_ROOT_EXTENDED / rel
                row = {
                    "scope": scope, "reason": reason,
                    "platform": platform, "species": species,
                    "dataset_group": group, "mod_type": mod_type,
                    "mod_display": MOD_LABEL.get(mod_type, mod_type),
                    "tool": tool, "sample": sample,
                    "n_rows": sum(1 for _ in f.open()) - 1,
                    "size_bytes": f.stat().st_size,
                    "callsets_path": str(rel),
                    "archived_path": str(rel) if scope == "out_of_scope" else "",
                    "action": "kept" if scope == "in_scope" else "moved",
                }
                rows.append(row)
                if scope == "out_of_scope":
                    move_file(f, archived, args.dry_run)

        # (3) also purge anything a previous run had ARCHIVED for the deleted
        #     species, so callsets_extended/ keeps only HeLa/Curlcake non-m6A.
        if CALLSET_ROOT_EXTENDED.exists():
            for platform, species, group, mod in iter_mod_dirs(CALLSET_ROOT_EXTENDED):
                mod_type = canonical_mod_type(mod.name)
                if delete_nonm6a(species, mod_type):
                    _delete_mod_dir(mod, CALLSET_ROOT_EXTENDED, platform, species,
                                    group, mod_type)
                    continue
                n_deleted_tool_files += _delete_out_of_scope_tools(
                    mod, CALLSET_ROOT_EXTENDED, platform, species, group, mod_type)

        if not args.dry_run:
            prune_empty(CALLSET_ROOT_EXTENDED)
            prune_empty(CALLSET_ROOT)

        # everything parked by a previous run: keep the manifest complete so that
        # ``--undo`` can always restore the full set.
        known = {r["callsets_path"] for r in rows}
        rows.extend(r for r in archived_rows() if r["callsets_path"] not in known)

    man = pd.DataFrame(rows)
    man["mod_type"] = pd.Categorical(man["mod_type"], MOD_TYPE_ORDER, ordered=True)
    man = man.sort_values(["scope", "species", "dataset_group", "mod_type", "tool",
                           "sample"]).reset_index(drop=True)
    if not args.dry_run:
        write_table(man, MANIFEST)
        inv.record(MANIFEST, n_rows=len(man))

    # keep a cumulative deletion record across re-runs (traceability): once the
    # directories are gone a later run finds nothing to delete, so merge with
    # whatever an earlier run already recorded.
    prior_deleted: list[dict] = []
    if DELETED_MANIFEST.exists():
        try:
            prior_deleted = read_table(DELETED_MANIFEST).to_dict("records")
        except Exception:  # noqa: BLE001 - a corrupt manifest must not block the run
            prior_deleted = []
    seen_del = {d["deleted_path"] for d in deleted}
    all_deleted = deleted + [d for d in prior_deleted
                             if d.get("deleted_path") not in seen_del]
    del_man = pd.DataFrame(all_deleted)
    if len(del_man):
        del_man["mod_type"] = pd.Categorical(del_man["mod_type"], MOD_TYPE_ORDER,
                                             ordered=True)
        del_man = del_man.sort_values(["tree", "species", "dataset_group", "mod_type",
                                       "tool", "sample"]).reset_index(drop=True)
    if not args.dry_run:
        write_table(del_man, DELETED_MANIFEST)
        inv.record(DELETED_MANIFEST, n_rows=len(del_man))

    # tools outside the manuscript's 15-tool set get their own cumulative ledger
    # (same merge rule: once the directories are gone a later run finds nothing).
    prior_tools: list[dict] = []
    if DELETED_TOOLS_MANIFEST.exists():
        try:
            prior_tools = read_table(DELETED_TOOLS_MANIFEST).to_dict("records")
        except Exception:  # noqa: BLE001 - a corrupt manifest must not block the run
            prior_tools = []
    seen_tools = {d["deleted_path"] for d in deleted_tools}
    all_deleted_tools = deleted_tools + [d for d in prior_tools
                                         if d.get("deleted_path") not in seen_tools]
    del_tools_man = pd.DataFrame(all_deleted_tools)
    if len(del_tools_man):
        del_tools_man["mod_type"] = pd.Categorical(del_tools_man["mod_type"],
                                                   MOD_TYPE_ORDER, ordered=True)
        del_tools_man = del_tools_man.sort_values(
            ["tree", "species", "dataset_group", "mod_type", "tool", "sample"]
        ).reset_index(drop=True)
    if not args.dry_run:
        write_table(del_tools_man, DELETED_TOOLS_MANIFEST)
        inv.record(DELETED_TOOLS_MANIFEST, n_rows=len(del_tools_man))

    # ------------------------------------------------- callsets_summary join --
    if SUMMARY.exists() and not args.dry_run:
        summary = read_table(SUMMARY)
        # scope follows from the row identity, never from a path string; rows
        # removed by this step (out-of-scope species OR out-of-scope tool) are
        # marked ``deleted`` (no path).
        summary["scope"] = [
            "deleted" if (
                delete_nonm6a(str(r["species"]), canonical_mod_type(r["mod_type"]))
                or not tool_in_scope(str(r["platform"]), str(r["species"]),
                                     str(r["dataset_group"]),
                                     canonical_mod_type(r["mod_type"]), str(r["tool"])))
            else scope_of(str(r["platform"]), str(r["species"]), str(r["dataset_group"]),
                          canonical_mod_type(r["mod_type"]))[0]
            for _, r in summary.iterrows()]
        summary["mod_type"] = [canonical_mod_type(m) for m in summary["mod_type"]]
        # rebuild ``out_file`` from the row's identity so a half-written value
        # from an earlier (buggy) pass can never survive.
        roots = {True: CALLSET_ROOT, False: CALLSET_ROOT_EXTENDED}
        new_out, n_exist = [], 0
        for _, r in summary.iterrows():
            if str(r["scope"]) == "deleted":
                new_out.append("")          # physically removed -> no path
                continue
            rel = Path(str(r["platform"]), str(r["species"]), str(r["dataset_group"]),
                       canonical_mod_type(r["mod_type"]), str(r["tool"]),
                       f"{r['sample']}.tsv")
            dest = roots[str(r["scope"]) == "in_scope"] / rel
            if dest.exists():
                n_exist += 1
            new_out.append(str(dest.relative_to(PROJECT)))
        summary["out_file"] = new_out
        write_table(summary, SUMMARY)
        inv.record(SUMMARY, n_rows=len(summary))
        logger.info("callsets_summary.csv: scope + out_file refreshed "
                    "(%d/%d files present on disk)", n_exist, len(summary))

    inv.flush()

    if len(man):
        logger.info("scope split summary:\n%s",
                    man.groupby(["scope", "mod_type"], observed=True)
                       .agg(files=("n_rows", "size"), rows=("n_rows", "sum"))
                       .to_string())
        out = man[man["scope"] == "out_of_scope"]
        logger.info("archived: %d files, %d rows -> %s", len(out),
                    int(out["n_rows"].sum()), CALLSET_ROOT_EXTENDED.relative_to(SITES_ROOT))
        kept = man[man["scope"] == "in_scope"]
        logger.info("kept: %d files, %d rows", len(kept), int(kept["n_rows"].sum()))
    if len(del_man):
        logger.info("DELETED (Arabidopsis/Mouse/E.coli non-m6A): %d files, %d rows "
                    "-> manifest/nonm6a_deleted.csv", len(del_man),
                    int(pd.to_numeric(del_man["n_rows"], errors="coerce").fillna(0).sum()))
    else:
        logger.info("DELETED: 0 files (no out-of-scope non-m6A present)")
    if len(del_tools_man):
        logger.info("DELETED (tools outside the manuscript's 15-tool set): %d files, "
                    "%d rows -> manifest/out_of_scope_tools_deleted.csv",
                    len(del_tools_man),
                    int(pd.to_numeric(del_tools_man["n_rows"],
                                      errors="coerce").fillna(0).sum()))
    else:
        logger.info("DELETED tools: 0 files (no out-of-scope tool present)")
    if args.dry_run:
        logger.info("DRY RUN - nothing was moved or deleted")
    logger.info("manifest -> %s", MANIFEST)
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
