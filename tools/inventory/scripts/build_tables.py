#!/usr/bin/env python3
"""Merge the three evidence layers into the supplementary tables (R1-5, R1-10).

Inputs (all produced inside this directory)
-------------------------------------------
* ``raw/commands_raw.csv``    -- exact command lines            (static_script)
* ``raw/params_raw.csv``      -- parameters from yml/log/md     (static_*)
* ``raw/versions_raw.csv``    -- versions from conda/pip/binary (dynamic)
* ``curated/tool_inventory_curated.csv`` -- manual curation     (curated)

Merging rule: **curated > dynamic probe > static yml/md > static log > static
script**.  Nothing is invented — unresolved cells are ``not recorded`` and are
listed in ``TI5_todo_report.csv``.

Outputs (``tables/``)
---------------------
* ``TI1_per_tool_implementation.csv`` / ``.md`` -- R1-5 ten fields + evidence
* ``TI2_command_lines.csv``                     -- R1-10 exact command lines
* ``TI3_software_versions.csv``                 -- tool / env / version / method
* ``TI4_model_checkpoints.csv``                 -- models and checkpoints
* ``TI5_todo_report.csv``                       -- everything still missing

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/tools/inventory/scripts/build_tables.py
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
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ti_common as tc  # noqa: E402

#: lower number == higher priority
PRIORITY: dict[str, int] = {
    "curated_confirmed": 0,
    "curated_check": 1,
    "conda_list": 2,
    "pip_show": 3,
    "bin_version": 4,
    "static_path": 5,
    "source_git": 5,
    "source_dir": 5,
    "source_file": 5,
    "static_yml": 6,
    "model_file": 6,
    "static_md": 6,
    "static_log": 7,
    "static_script": 8,
}

MAX_VALUES = 4      # per cell, before collapsing with "..."
MAX_EVIDENCE = 3   # provenance entries quoted per cell


def _read(path: Path, cols: list[str]) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=cols)
    df = pd.read_csv(path).fillna("")
    for c in cols:
        if c not in df.columns:
            df[c] = ""
    return df


def collect_candidates() -> dict[tuple[str, str], list[dict]]:
    """(tool, field) -> list of {value, priority, evidence}."""
    cands: dict[tuple[str, str], list[dict]] = {}

    def add(tool: str, field: str, value: str, method: str, evidence: str) -> None:
        tool = str(tool).strip()
        value = str(value).strip()
        if not tool or not value or value == tc.MISSING or value.lower() == "nan":
            return
        if field not in tc.FIELDS:
            return
        cands.setdefault((tool, field), []).append(
            {"value": value, "priority": PRIORITY.get(method, 9),
             "evidence": f"{evidence} [{method}]"})

    # ---- curated --------------------------------------------------------
    cur = _read(tc.CURATED_DIR / "tool_inventory_curated.csv",
                ["tool_canonical", "field", "value", "evidence", "status"])
    for _, r in cur.iterrows():
        status = str(r["status"]).strip().upper()
        if status not in {"CONFIRMED", "NEEDS_CHECK"} or not str(r["value"]).strip():
            continue
        method = "curated_confirmed" if status == "CONFIRMED" else "curated_check"
        add(r["tool_canonical"], r["field"], r["value"], method,
            str(r["evidence"]) or "manual curation")

    # ---- versions -------------------------------------------------------
    ver = _read(tc.RAW_DIR / "versions_raw.csv",
                ["entity_type", "entity", "package", "version", "method",
                 "evidence"])
    for _, r in ver.iterrows():
        if str(r["entity_type"]) != "tool":
            continue
        pkg = str(r["package"])
        val = f"{pkg} {r['version']}".strip() if pkg not in {"", "-"} \
            else str(r["version"])
        m = re.search(r"conda list -n (\S+)", str(r["evidence"]))
        if m:
            val = f"{val} [env:{m.group(1)}]"
        add(r["entity"], "software_or_version", val, str(r["method"]),
            str(r["evidence"]))

    # ---- params ---------------------------------------------------------
    par = _read(tc.RAW_DIR / "params_raw.csv",
                ["tool_canonical", "dataset", "field", "value", "source",
                 "line_no", "method"])
    for _, r in par.iterrows():
        add(r["tool_canonical"], r["field"], r["value"], str(r["method"]),
            f"{r['source']}:{r['line_no']}")

    # ---- vendored source trees -----------------------------------------
    srcf = _read(tc.RAW_DIR / "source_raw.csv",
                 ["tool_canonical", "dataset", "field", "value", "source",
                  "line_no", "method"])
    for _, r in srcf.iterrows():
        ev = f"{r['source']}:{r['line_no']}" if str(r["line_no"]) not in {"", "0"} \
            else str(r["source"])
        add(r["tool_canonical"], r["field"], r["value"], str(r["method"]), ev)

    return cands


def resolve(cands: dict[tuple[str, str], list[dict]], tool: str, field: str
            ) -> tuple[str, str, str]:
    """Return (value, evidence, status) for one cell."""
    items = sorted(cands.get((tool, field), []), key=lambda d: (d["priority"],
                                                                d["value"]))
    if not items:
        return tc.MISSING, "", "TODO"
    best = items[0]["priority"]
    group = [i for i in items if i["priority"] == best]
    values, evids = [], []
    for i in group:
        if i["value"] not in values:
            values.append(i["value"])
        if i["evidence"] not in evids:
            evids.append(i["evidence"])
    value = " | ".join(values[:MAX_VALUES])
    if len(values) > MAX_VALUES:
        value += " | ..."
    evidence = "; ".join(evids[:MAX_EVIDENCE])
    status = {0: "curated", 1: "curated (needs check)", 2: "auto (conda)",
              3: "auto (pip)", 4: "auto (binary --version)"}.get(
        best, "auto (static artefact)")
    return value, evidence, status


def benchmark_toolset(logger) -> set[str]:
    """Tools that actually produced a callset / result directory.

    A tool counts as *included in the benchmark* when it has a directory under
    ``result/`` or ``result_RNA004/**/RNA004_result/``, or a callset file under
    ``output/<sample>/``.  This is the objective answer to R1-6 ("why is a tool
    in one figure but not another").
    """
    tools: set[str] = set()
    for d in tc.RESULT.iterdir():
        if d.is_dir():
            c = tc.canonical_known(d.name)
            if c:
                tools.add(c)
    for d in (tc.RNA004_ROOT / "HeLa" / "raw_calls" / "RNA004_result",
              tc.RNA004_ROOT / "Curlcake" / "data" / "RNA004_result"):
        if d.exists():
            for sub in d.iterdir():
                if sub.is_dir():
                    c = tc.canonical_known(sub.name)
                    if c:
                        tools.add(c)
    out_root = (_XB / "archive/output")
    if out_root.exists():
        for f in out_root.glob("*/*.txt"):
            if f.name.lower() == "tools.txt":
                continue
            c = tc.canonical_known(f.stem)
            if c:
                tools.add(c)
    tools.add("Dorado")  # RNA004 modification calling, evidenced by IMPORTANT_INFO.md
    logger.info("tools with a result directory or callset: %s", sorted(tools))
    return tools


def main() -> None:
    logger = tc.setup_logger("build_tables")
    tc.TABLE_DIR.mkdir(parents=True, exist_ok=True)
    cands = collect_candidates()
    logger.info("candidate facts: %d (tool, field) pairs", len(cands))

    cmds = _read(tc.RAW_DIR / "commands_raw.csv",
                 ["tool_canonical", "tool_raw", "dataset", "conda_env", "command",
                  "source", "line_no", "method"])
    bench = benchmark_toolset(logger)
    vers = _read(tc.RAW_DIR / "versions_raw.csv",
                 ["entity_type", "entity", "package", "version", "method",
                  "evidence"])

    # ---------------------------------------------------------------- TI1
    seen = {str(t) for t in cmds["tool_canonical"].unique()
            if str(t) and str(t) != "unclassified"}
    tools = sorted(set(tc.TOOL_META) | bench | seen)
    rows = []
    for tool in tools:
        cat, mod, role = tc.TOOL_META.get(tool, ("unclassified", "n/a", "n/a"))
        sub = cmds[cmds["tool_canonical"] == tool]
        rna004 = bool((sub["dataset"] == "RNA004").any())
        rna002 = bool((sub["dataset"] != "RNA004").any())
        ev_lines = [f"{r['source']}:{r['line_no']}" for _, r in sub.head(3).iterrows()]
        row = {"tool_canonical": tool, "category": cat, "modification": mod,
               "role": role, "in_benchmark_callsets": tool in bench,
               "ran_on_RNA002": rna002, "ran_on_RNA004": rna004,
               "n_command_lines": int(len(sub)),
               "command_line_evidence": "; ".join(ev_lines)}
        n_todo = 0
        for field in tc.FIELDS:
            value, evidence, status = resolve(cands, tool, field)
            row[field] = value
            row[f"{field}_evidence"] = evidence
            if value == tc.MISSING and field != "notes":
                n_todo += 1
        row["status"] = ("complete" if n_todo == 0 else
                         f"partial ({n_todo} field(s) missing)")
        rows.append(row)

    ti1 = pd.DataFrame(rows)
    cols = (["tool_canonical", "category", "modification", "role",
             "in_benchmark_callsets", "ran_on_RNA002", "ran_on_RNA004",
             "n_command_lines"] +
            tc.FIELDS +
            [f"{f}_evidence" for f in tc.FIELDS if f != "notes"] +
            ["command_line_evidence", "status"])
    ti1 = ti1[[c for c in cols if c in ti1.columns]]
    ti1["generated_at"] = tc.stamp()
    ti1.to_csv(tc.TABLE_DIR / "TI1_per_tool_implementation.csv", index=False)
    logger.info("TI1: %d tools", len(ti1))

    # ---------------------------------------------------------------- TI2
    ti2 = cmds.drop_duplicates(
        subset=["tool_canonical", "dataset", "command", "source", "line_no"])
    ti2 = ti2.sort_values(["tool_canonical", "dataset", "source", "line_no"])
    ti2["generated_at"] = tc.stamp()
    ti2.to_csv(tc.TABLE_DIR / "TI2_command_lines.csv", index=False)
    logger.info("TI2: %d command lines", len(ti2))

    # ---------------------------------------------------------------- TI3
    ti3 = vers.copy()
    ti3["generated_at"] = tc.stamp()
    ti3.to_csv(tc.TABLE_DIR / "TI3_software_versions.csv", index=False)
    ok = int((ti3["version"] != tc.MISSING).sum())
    logger.info("TI3: %d rows (%d with a version)", len(ti3), ok)

    # ---------------------------------------------------------------- TI4
    par = _read(tc.RAW_DIR / "params_raw.csv",
                ["tool_canonical", "dataset", "field", "value", "source",
                 "line_no", "method"])
    models = par[par["field"] == "model_checkpoint"].copy()
    cur = _read(tc.CURATED_DIR / "tool_inventory_curated.csv",
                ["tool_canonical", "field", "value", "evidence", "status"])
    cm = cur[(cur["field"] == "model_checkpoint") & (cur["value"] != "")]
    cm = cm.rename(columns={"evidence": "source"})
    cm["dataset"] = "all"
    cm["line_no"] = 0
    cm["method"] = "curated"
    srcf = _read(tc.RAW_DIR / "source_raw.csv",
                 ["tool_canonical", "dataset", "field", "value", "source",
                  "line_no", "method"])
    sm = srcf[srcf["field"] == "model_checkpoint"]
    ti4 = pd.concat([models[["tool_canonical", "dataset", "value", "source",
                             "line_no", "method"]],
                     sm[["tool_canonical", "dataset", "value", "source",
                         "line_no", "method"]],
                     cm[["tool_canonical", "dataset", "value", "source",
                         "line_no", "method"]]], ignore_index=True)
    ti4 = ti4.drop_duplicates(subset=["tool_canonical", "value"])
    ti4 = ti4.sort_values(["tool_canonical", "value"])
    ti4["generated_at"] = tc.stamp()
    ti4.to_csv(tc.TABLE_DIR / "TI4_model_checkpoints.csv", index=False)
    logger.info("TI4: %d model / checkpoint rows", len(ti4))

    # ---------------------------------------------------------------- TI5
    def attempt_note(tool: str) -> str:
        """What was already tried for this tool, so nobody repeats it."""
        parts = []
        ms = sorted({str(m) for m in
                     vers[(vers["entity_type"] == "tool")
                          & (vers["entity"] == tool)]["method"]})
        if ms:
            parts.append("probed: " + ", ".join(ms))
        sub = srcf[(srcf["tool_canonical"] == tool)
                   & (srcf["method"].isin(["source_dir", "source_git",
                                           "source_file"]))]
        dirs = sorted({Path(str(s)).name for s in sub["source"]})[:3]
        if dirs:
            parts.append("source dirs: " + ", ".join(dirs))
        n_cmd = int((cmds["tool_canonical"] == tool).sum())
        if n_cmd:
            parts.append(f"{n_cmd} command lines scanned")
        return "; ".join(parts) or "no automatic probe reached this tool"

    todo = []
    for tool in tools:
        for field in tc.FIELDS:
            if field == "notes":
                continue
            value, evidence, status = resolve(cands, tool, field)
            if value == tc.MISSING:
                todo.append({"tool_canonical": tool, "field": field,
                             "status": "TODO",
                             "how_to_obtain": tc.HOW_TO_OBTAIN.get(field, ""),
                             "attempts": attempt_note(tool)})
    ti5 = pd.DataFrame(todo, columns=["tool_canonical", "field", "status",
                                      "how_to_obtain", "attempts"])
    ti5["generated_at"] = tc.stamp()
    ti5.to_csv(tc.TABLE_DIR / "TI5_todo_report.csv", index=False)
    logger.info("TI5: %d fields still to be filled in", len(ti5))

    # ---------------------------------------------------------------- TI6
    unc = cmds[cmds["tool_canonical"] == "unclassified"]
    if len(unc):
        ti6 = (unc.groupby("tool_raw")
                  .agg(n_commands=("command", "size"),
                       datasets=("dataset", lambda s: ";".join(sorted(set(s))[:4])),
                       example_command=("command", "first"),
                       script=("source", "first"))
                  .reset_index()
                  .sort_values("n_commands", ascending=False))
        ti6["generated_at"] = tc.stamp()
        ti6.to_csv(tc.TABLE_DIR / "TI6_unclassified_commands.csv", index=False)
        logger.info("TI6: %d scripts whose tool could not be attributed "
                    "(needs manual classification)", len(ti6))
    else:
        logger.info("TI6: every command line was attributed to a known tool")

    # ---------------------------------------------------------------- MD
    md = ["# TI1 — Per-tool implementation details (Supplementary Table)\n",
          "Addresses **R1-5 / E4** (software/version or commit, model/checkpoint, "
          "required input, minimum read coverage, filtering parameters, "
          "probability/p-value thresholds, multiple-testing correction, default "
          "vs optimised parameters, coordinate harmonisation).\n",
          f"Generated: {tc.stamp()}\n",
          "| Tool | Category | Modification | In callsets | Version | "
          "Model / checkpoint | Required input | Min coverage | Calling threshold | "
          "Multiple testing | Default/optimised | Coordinates | Commands | Evidence |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in ti1.iterrows():
        md.append("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
            r["tool_canonical"], r["category"], r["modification"],
            "yes" if r["in_benchmark_callsets"] else "no",
            r["software_or_version"], r["model_checkpoint"], r["required_input"],
            r["min_read_coverage"], r["calling_threshold"],
            r["multiple_testing_correction"], r["default_vs_optimised"],
            r["coordinate_harmonisation"], r["n_command_lines"],
            (r["command_line_evidence"] or "-")[:180]))
    md.append("\nProvenance for every non-empty cell is in "
              "`TI1_per_tool_implementation.csv` (`*_evidence` columns); the exact "
              "command lines are in `TI2_command_lines.csv` (R1-10).\n")
    (tc.TABLE_DIR / "TI1_per_tool_implementation.md").write_text(
        "\n".join(md) + "\n", encoding="utf-8")

    logger.info("all tables written -> %s", tc.TABLE_DIR)


if __name__ == "__main__":
    main()
