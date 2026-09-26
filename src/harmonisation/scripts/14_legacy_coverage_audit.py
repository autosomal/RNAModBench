#!/usr/bin/env python3
"""14 -- legacy coverage audit: keep EVERY file the legacy assembly used.

Independent, reproducible reconciliation against two legacy evidence sources
(nothing here is derived from the harmonisation registry itself):

1. the seven ``cp.sh`` copy maps (``$RNAMODBENCH_LOCAL/code/code/{next,python}_postprocessing``;
   the pre-2026-09-15 location was ``tool_scripts/``) -- the exact file list the
   published ``output/`` tree was assembled from (``cp <src> <dst>``).  Those
   lines record the **legacy** ``result/`` paths, which
   :func:`legacy_src_current` replays onto today's tree through the
   ``result_tidy/rename_map_step{1,2,3}`` maps; and
2. the assembled ``output/`` tree (one directory per tool per group; now at
   ``config.LEGACY_OUTPUT``, the archived copy); the derived ``*_purified``
   groups (WT-called / KO-absent site lists) and non-callset helpers are
   skipped, since they are analysis products, not callsets.

Each entry is resolved to a canonical ``(sample, tool)``:

* group-level comparison directories through
  ``legacy_liftover.EXPLICIT_DIR_SAMPLE`` (the tool's own "test side"), then
* directory-name aliases, then
* the Curlcake construct library recovered from the **file stem**
  (``Curlcake_m6A_result1_xPore.txt`` -> ``Curlcake_m6A_rep1``).

and matched against ``manifest/completeness_audit.csv``.  Verdicts:

``covered``            callset present with >= 1 site
``covered_zero``       header-only callset (genuine ``Total_Detected = 0``)
``legacy_unused``      result exists, the legacy assembly skipped it on purpose
                       (``legacy_liftover.LEGACY_UNUSED_RAW``)
``out_of_scope``       non-m6A outside the manuscript scope (deleted/archived)
``group_covered``      cp.sh row is group-level; the group has data for the tool
``legacy_src_missing`` the legacy source file itself is absent (typo, ...)
``legacy_empty_dir``   the legacy output directory is empty (never used)
``unmapped_tool``      the legacy label could not be mapped onto a tool id
``MISMATCH``           the legacy assembly used a file we have NO callset for

``MISMATCH`` is the only blocking verdict: it is what silently dropped
``result/xPore/result{1..4}``, ``result/yanocomp/result{1..4}``,
``result/DRUMMER/result{3,4}`` and ``result/*/E_IVT_neg_result`` before.

Output: harmonisation/manifest/legacy_coverage_audit.csv (+ console summary).

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/14_legacy_coverage_audit.py
"""

from __future__ import annotations

# --- RNAModBench path bootstrap (added when this file was deposited) ----------
import os as _rb_os, pathlib as _rb_pl


def _rb_find(start):
    for p in (start, *start.parents):
        if (p / "RNAMOD_BENCH_ROOT").exists():
            return p
    return start


_RB = _rb_pl.Path(_rb_os.environ.get("RNAMODBENCH_ROOT") or _rb_find(_rb_pl.Path(__file__).resolve().parent))
_XB = _rb_pl.Path(_rb_os.environ.get("RNAMODBENCH_LOCAL") or (_RB / "_local"))
# --------------------------------------------------------------------------- #
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import (CONVERTED_CALLSETS, LEGACY_OUTPUT, MANIFEST_DIR,
                           PROJECT, RESULT_RNA002, RNA002_TOOLS)
from common.io_utils import ensure_dirs, read_table, write_table
from common.legacy_liftover import LEGACY_UNUSED_RAW
from common.manifest import Inventory, log_time, setup_logger
from common.registry import (curlcake_sample_from_stem,
                             infer_tool_from_filename, legacy_tool_map,
                             sample_of_dir)

#: the copy maps of the legacy assembly (verbatim ``cp <src> <dst>`` lines).
#: Each entry lists candidate locations, first existing wins: the 2026-09-15
#: reorganization split the tree between ``$RNAMODBENCH_LOCAL/tool_scripts/'' (verbatim copy the
#: reconciliation evidence was extracted from) and ``$RNAMODBENCH_LOCAL/code/code/'' (kept after
#: the dedup).  Their content still records the *legacy* ``result/`` paths --
#: :func:`legacy_src_current` resolves those onto the current tree.
CP_SH_FILES: tuple[tuple[str, ...], ...] = (
    ("$RNAMODBENCH_LOCAL/tool_scripts/next_postprocessing/Arabidopsis/m6A/cp.sh",
     "$RNAMODBENCH_LOCAL/code/code/next_postprocessing/Arabidopsis/m6A/cp.sh"),
    ("$RNAMODBENCH_LOCAL/tool_scripts/next_postprocessing/Mouse/m6A/cp.sh",
     "$RNAMODBENCH_LOCAL/code/code/next_postprocessing/Mouse/m6A/cp.sh"),
    ("$RNAMODBENCH_LOCAL/tool_scripts/next_postprocessing/Hela/m6A/cp.sh",
     "$RNAMODBENCH_LOCAL/code/code/next_postprocessing/Hela/m6A/cp.sh"),
    ("$RNAMODBENCH_LOCAL/code/code/next_postprocessing/Hela/other/cp.sh",),
    ("$RNAMODBENCH_LOCAL/code/code/next_postprocessing/Curlcake/cp.sh",),
    ("$RNAMODBENCH_LOCAL/code/code/next_postprocessing/E.coli/cp.sh",),
    ("$RNAMODBENCH_LOCAL/code/code/python_postprocessing/HeLa/cp.sh",),
)

#: output/ groups that are analysis products, not callsets.
_DERIVED_GROUP_SUFFIXES = ("_purified",)
_DERIVED_GROUP_NAMES = {"RRACH", "total"}

#: file suffixes that never carry callsites inside the legacy output tree.
_NON_CALLSET_SUFFIXES = (".log", ".md", ".json", ".gzip", ".ipynb", ".sh")

_KNOWN_TOOLS: set[str] = {t.tool for t in RNA002_TOOLS}
#: completeness-audit statuses that mean "the rebuilt tree has this data".
_OK_STATUSES = {"filled", "ok_zero"}
#: completeness-audit statuses that mean "outside the manuscript's scope by
#: design": deleted/archived files **and** tool x sample combinations the
#: manuscript never used (``13`` splits those into four labels; all are
#: non-blocking here, otherwise e.g. the legacy ELIGOS2_diff run on the E.coli
#: IVT control -- a 33-byte empty callset the paper never used -- reads as a
#: blocking MISMATCH).
_OUT_OF_SCOPE_STATUSES = {"out_of_scope_deleted", "out_of_scope_archived",
                          "out_of_scope_combination", "out_of_scope_tool"}


def _stem(name: str) -> str:
    return Path(name).stem


def _group_from_dst(dst: str, groups: set[str]) -> str:
    """Output group of a cp.sh destination path (``/output/<group>/...``)."""
    parts = Path(dst).parts
    if "output" in parts:
        idx = parts.index("output")
        if len(parts) > idx + 1 and parts[idx + 1] in groups:
            return parts[idx + 1]
    for g in sorted(groups, key=len, reverse=True):
        if f"/{g}/" in dst or dst.rstrip("/").endswith(f"/{g}"):
            return g
    return ""


def _tool_label(text: str) -> str:
    """Tool id of a destination path component or file name ('' = unknown)."""
    last = Path(text.rstrip("/")).name
    mapped = legacy_tool_map().get(last)
    if mapped:
        return mapped
    mapped = legacy_tool_map().get(last.split(".")[0])
    if mapped:
        return mapped
    return ""


def _resolve_sample(tool: str, parent: str, stem: str) -> str:
    """Canonical sample of a legacy source file ('' when group-level)."""
    sample, _ = sample_of_dir(tool, parent)
    if sample is not None:
        return sample
    return curlcake_sample_from_stem(stem) or ""


# --------------------------------------------------------------------------- #
# legacy path resolution: $RNAMODBENCH_LOCAL/raw/result -> today's tree
# --------------------------------------------------------------------------- #
#: The cp.sh maps were written against the pre-2026-09-15 tree
#: (``$RNAMODBENCH_LOCAL/raw/result<tool>/<sample_dir>/<file>``).  The reorganization moved
#: the tree to ``$RNAMODBENCH_LOCAL/raw/result`` and the tidy pass renamed every
#: tool/sample directory (sometimes twice), so checking existence at the recorded
#: path marks everything as ``legacy_src_missing``.  The relocations were
#: recorded -- :func:`legacy_src_current` replays them instead.
_LEGACY_RESULT_ROOT = str(_XB / "raw/result")
_TIDY_DIR = (_RB / "analysis/result_tidy")
_RELOCATION_MAPS = ("rename_map_step1.csv", "rename_map_step2.csv",
                    "rename_map_step3_curlcake.csv",
                    "rename_map_step4_pair_dirs.csv")
#: roots that now hold the files the legacy assembly copied (converted callsets
#: first: they are the promoted ``*_remove_chr.txt``/``*_<tool>.txt`` files)
_RESOLVE_ROOTS = (CONVERTED_CALLSETS, RESULT_RNA002)


def _load_relocations() -> tuple[dict[str, str], dict[str, str]]:
    """(directory renames, tool renames) replayed from the tidy maps.

    ``rename_map_step{1,2,3}`` record ``<tool>/<sample_dir>`` renames (kinds
    ``sample_dir`` / ``archived_sample_dir`` / ``tool_internal_dir``) and whole
    tool-directory renames (``Epinano_DiffErr`` -> ``EpiNano_DiffErr``).
    """
    dirs: dict[str, str] = {}
    tools: dict[str, str] = {}
    for name in _RELOCATION_MAPS:
        path = _TIDY_DIR / name
        if not path.exists():
            continue
        for _, r in read_table(path, sep=",").iterrows():
            kind, old, new = str(r["kind"]), str(r["old"]), str(r["new"])
            if kind in ("sample_dir", "archived_sample_dir", "tool_internal_dir"):
                dirs[old] = new
            elif kind in ("tool_dir", "tool_dir_out_of_scope") and "/" not in old:
                tools[old] = new
    return dirs, tools


_DIR_RENAME, _TOOL_RENAME = _load_relocations()


def _rename_dir(tool: str, sample_dir: str) -> tuple[str, str]:
    """Replay the recorded renames on a ``<tool>/<sample_dir>`` pair.

    The maps are applied in a loop because one directory can be renamed more than
    once (``result1`` -> ``Curlcake_IVT_result1`` -> ``..._rep2_partial_vs_IVT_rep1``);
    every lookup is tried with the current tool spelling and with the original one,
    because the ``Epinano_DiffErr`` rows of step 1 predate the tool rename.
    """
    cur_tool = _TOOL_RENAME.get(tool, tool)
    cur_dir = sample_dir
    for _ in range(6):
        key = f"{cur_tool}/{cur_dir}"
        new = _DIR_RENAME.get(key) or _DIR_RENAME.get(f"{tool}/{cur_dir}")
        if not new or new == key:
            break
        cur_tool, _, cur_dir = new.partition("/")
    return cur_tool, cur_dir


def legacy_src_current(src: Path) -> Path | None:
    """The file a legacy ``result/...`` path points at today (``None`` if gone).

    The tidy pass renamed directories only, so the recorded file name is kept.
    Both current roots are tried; a handful of legacy rows are genuine typos
    (e.g. ``..._ELIGOS2_solo.txtt``) and stay ``legacy_src_missing`` by design.
    """
    text = str(src)
    if not text.startswith(_LEGACY_RESULT_ROOT + "/"):
        return src if src.exists() else None
    parts = text[len(_LEGACY_RESULT_ROOT) + 1:].split("/")
    if len(parts) < 3:
        return None
    tool, sample_dir = _rename_dir(parts[0], parts[1])
    name = parts[-1]
    rel_file = "/".join(parts[2:])
    for root in _RESOLVE_ROOTS:
        for cand in (root / tool / sample_dir / name,
                     root / tool / sample_dir / rel_file):
            if cand.exists():
                return cand
    return None


def _output_tree_groups() -> tuple[dict[str, dict[str, list[Path]]],
                                   list[tuple[str, str]]]:
    """``({group: {tool: [files]}}, [(group, tool) for empty tool dirs])``."""
    groups: dict[str, dict[str, list[Path]]] = {}
    empty_dirs: list[tuple[str, str]] = []
    if not LEGACY_OUTPUT.is_dir():
        return groups, empty_dirs
    for gdir in sorted(p for p in LEGACY_OUTPUT.iterdir() if p.is_dir()):
        if gdir.name.endswith(_DERIVED_GROUP_SUFFIXES) or gdir.name in _DERIVED_GROUP_NAMES:
            continue
        found: dict[str, list[Path]] = {}
        for f in gdir.rglob("*"):
            if not f.is_file() or f.stat().st_size == 0:
                continue
            if f.name.lower().endswith(_NON_CALLSET_SUFFIXES):
                continue
            tool = _tool_label(str(f)) or infer_tool_from_filename(f.name)
            if tool:
                found.setdefault(tool, []).append(f)
        for tdir in sorted(p for p in gdir.rglob("*") if p.is_dir()):
            tool = _tool_label(str(tdir))
            if tool and tool not in found:
                empty_dirs.append((gdir.name, tool))
        if found:
            groups[gdir.name] = found
    return groups, empty_dirs


def main() -> None:
    logger = setup_logger("14_legacy_coverage_audit")
    inv = Inventory("14_legacy_coverage_audit")
    ensure_dirs(MANIFEST_DIR)

    audit = read_table(MANIFEST_DIR / "completeness_audit.csv")
    status_by = {(r["sample"], r["tool"]): r["status"] for _, r in audit.iterrows()}
    group_tools = {(r["dataset_group"], r["tool"]) for _, r in audit.iterrows()
                   if r["status"] in _OK_STATUSES}
    out_groups, empty_dirs = _output_tree_groups()
    known_groups = set(out_groups) | set(audit["dataset_group"].astype(str))
    if LEGACY_OUTPUT.is_dir():
        known_groups |= {p.name for p in LEGACY_OUTPUT.iterdir() if p.is_dir()}

    rows: list[dict] = []
    with log_time(logger, "legacy coverage audit"):
        # ---------------------------------------------------- 1. cp.sh maps --
        for candidates in CP_SH_FILES:
            path = next((PROJECT / c for c in candidates
                         if (PROJECT / c).exists()), None)
            if path is None:
                logger.warning("cp.sh not found: %s", " | ".join(candidates))
                continue
            rel = str(path.relative_to(PROJECT))
            for lineno, line in enumerate(path.read_text().splitlines(), start=1):
                line = line.strip()
                if not line.startswith("cp "):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                src, dst = Path(parts[1]), parts[2]
                tool = _tool_label(dst) or infer_tool_from_filename(src.name) or ""
                group = _group_from_dst(dst, known_groups)
                sample = _resolve_sample(tool, src.parent.name, _stem(src.name))
                current = legacy_src_current(src)
                exists = current is not None
                declared_unused = LEGACY_UNUSED_RAW.get((sample, tool), "")
                declared_path = PROJECT / declared_unused if declared_unused else None

                if not exists:
                    if declared_path is not None and declared_path.exists():
                        verdict, note = ("legacy_unused",
                                         "declared in LEGACY_UNUSED_RAW "
                                         "(file renamed by the tidy pass)")
                    else:
                        verdict, note = ("legacy_src_missing",
                                         "legacy source file absent on disk")
                elif not tool:
                    verdict, note = "unmapped_tool", "legacy label not mapped to a tool id"
                elif declared_path is not None and current == declared_path:
                    verdict, note = "legacy_unused", "result exists, legacy assembly skipped it"
                elif (sample, tool) in status_by:
                    st = status_by[(sample, tool)]
                    if st in _OK_STATUSES:
                        verdict, note = ("covered" if st == "filled" else "covered_zero"), ""
                    elif st == "raw_available_legacy_unused":
                        verdict, note = "legacy_unused", "result exists, legacy assembly skipped it"
                    elif st in _OUT_OF_SCOPE_STATUSES:
                        verdict, note = "out_of_scope", st
                    else:
                        verdict, note = "MISMATCH", f"audit status={st}"
                elif (group, tool) in group_tools:
                    verdict, note = "group_covered", "group-level cp.sh row"
                else:
                    verdict, note = "MISMATCH", "no callset for this legacy source"
                rows.append({"source": f"cp.sh:{rel}", "line_no": lineno,
                             "group": group, "sample": sample, "tool": tool,
                             "verdict": verdict, "src_file": str(src),
                             "src_resolved": str(current) if current else "",
                             "src_exists": exists, "dst": dst, "note": note})

        # ------------------------------------------------ 2. output/ tree --
        for group, tools in sorted(out_groups.items()):
            for tool, files in sorted(tools.items()):
                if (group, tool) in group_tools:
                    verdict, note = "group_covered", f"{len(files)} file(s) in output/"
                else:
                    verdict, note = "MISMATCH", "output/ has data, no rebuilt callset"
                rows.append({"source": "output/", "line_no": 0, "group": group,
                             "sample": "", "tool": tool, "verdict": verdict,
                             "src_file": str(files[0]),
                             "src_resolved": str(files[0]), "src_exists": True,
                             "dst": f"output/{group}", "note": note})
        for group, tool in sorted(set(empty_dirs)):
            rows.append({"source": "output/", "line_no": 0, "group": group,
                         "sample": "", "tool": tool, "verdict": "legacy_empty_dir",
                         "src_file": f"output/{group}/{tool}", "src_resolved": "",
                         "src_exists": True, "dst": f"output/{group}",
                         "note": "legacy output directory is empty (never used)"})

        # --------------------------- 3. declared legacy-unused results ------
        for (sample, tool), rel in sorted(LEGACY_UNUSED_RAW.items()):
            p = PROJECT / rel
            rows.append({"source": "LEGACY_UNUSED_RAW", "line_no": 0,
                         "group": "", "sample": sample, "tool": tool,
                         "verdict": "legacy_unused" if p.exists() else "legacy_src_missing",
                         "src_file": str(p),
                         "src_resolved": str(p) if p.exists() else "",
                         "src_exists": p.exists(),
                         "dst": "", "note": "declared legacy-unused result"})

    df = pd.DataFrame(rows)
    write_table(df, MANIFEST_DIR / "legacy_coverage_audit.csv")
    inv.record(MANIFEST_DIR / "legacy_coverage_audit.csv", n_rows=len(df))
    #: old -> resolved path ledger: the audit trail of ``legacy_src_current``
    cp_rows = df[df["source"].str.startswith("cp.sh")].copy()
    ledger = cp_rows[["sample", "tool", "group", "verdict", "src_file",
                      "src_resolved", "src_exists", "note"]]
    write_table(ledger, MANIFEST_DIR / "legacy_src_resolution.csv")
    inv.record(MANIFEST_DIR / "legacy_src_resolution.csv", n_rows=len(ledger))
    inv.flush()

    logger.info("legacy coverage audit: %d entries", len(df))
    logger.info("verdicts: %s", df["verdict"].value_counts().to_dict())
    bad = df[df["verdict"] == "MISMATCH"]
    if len(bad):
        logger.warning("MISMATCH: %d legacy-used file(s) have NO rebuilt callset",
                       len(bad))
        for _, r in bad.iterrows():
            logger.warning("   %-18s %-16s %s", r["group"] or r["sample"], r["tool"],
                           r["src_file"])
    else:
        logger.info("every file the legacy assembly used is covered by the rebuild")
    for verdict in ("legacy_src_missing", "legacy_empty_dir", "legacy_unused",
                    "unmapped_tool"):
        sub = df[df["verdict"] == verdict]
        if len(sub):
            logger.info("%s (%d):", verdict, len(sub))
            for _, r in sub.iterrows():
                logger.info("   %-18s %-16s %s", r["group"] or r["sample"], r["tool"],
                            r["note"])
    logger.info("manifest -> %s", MANIFEST_DIR / "legacy_coverage_audit.csv")
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
