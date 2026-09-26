#!/usr/bin/env python3
"""Apply the researched R1-5 facts to the tool inventory (reviewer comment R1-5).

What it does
------------
1. Reads ``tables/per_tool_implementation.csv`` and ``raw/research_facts.csv``
   (produced by ``parse_research.py`` from ``research/*.md``).
2. Fills every ``not recorded`` cell and, where the research evidence contradicts a
   previously auto-guessed cell, replaces it -- the displaced text is kept in
   ``tables/change_log.csv`` so nothing is silently dropped.
3. Writes the curation layer ``curated/tool_inventory_curated.csv`` (the only
   hand-editable file in this directory).  Because the legacy merge ranks curated
   facts above every auto probe, the corrected values survive a future
   ``build_tables.py`` run.  ``notes`` is deliberately *not* curated -- the legacy
   template injects an unevidenced note that would otherwise stick forever.
4. Rewrites ``tables/per_tool_implementation.csv`` (+ ``.md``) with live paths
   and the merged values, and refreshes ``tables/todo_report.csv``.

Inputs are never modified; the previous tables are copied to ``*.bak_<ts>``.

Usage
-----
    python3 apply_research.py [--dry-run] [--include-low]
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ti_paths as tp  # noqa: E402

CURATED_COLS = ["tool_canonical", "category", "modification", "role", "field",
                "value", "evidence", "status", "how_to_obtain"]
#: how many distinct values are kept per cell before collapsing with "..."
MAX_VALUES = 4
MAX_EVIDENCE = 3
CURATED_FIELDS = [f for f in tp.FIELDS if f != "notes"]
#: statements the legacy template injected without evidence -- safe to displace
TEMPLATE_STRINGS = {
    "transcript -> genomic via BED12/gtf2bed12 transcript alignment",
    "tool_scripts/detection/<species>/<sample>/<tool>.sh (see TI2)",
    "post-processing scripts in tool_scripts/*postprocessing",
}


def split_values(cell: str) -> list[str]:
    return [v.strip() for v in (cell or "").split("|") if v.strip()]


def dedupe(values: list[str]) -> list[str]:
    seen, out = set(), []
    for v in values:
        key = " ".join(v.lower().split())
        if key not in seen:
            seen.add(key)
            out.append(v)
    return out


def join_cell(values: list[str]) -> str:
    if len(values) <= MAX_VALUES:
        return " | ".join(values)
    return " | ".join(values[:MAX_VALUES]) + " | ..."


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--include-low", action="store_true",
                    help="also use facts the research notes rated confidence=low")
    args = ap.parse_args()

    ti1_path = tp.TABLE_DIR / "per_tool_implementation.csv"
    ti1 = tp.read_csv(ti1_path)
    if not ti1:
        print("TI1 missing", file=sys.stderr)
        return 1
    meta = {r["tool_canonical"]: (r["category"], r["modification"], r["role"])
            for r in ti1}
    facts = [f for f in tp.read_csv(tp.RAW_DIR / "research_facts.csv")
             if f["field"] in CURATED_FIELDS and f["tool_canonical"] in meta]
    if not args.include_low:
        dropped = [f for f in facts if f["confidence"] == "low"]
        facts = [f for f in facts if f["confidence"] != "low"]
        if dropped:
            print(f"skipped {len(dropped)} low-confidence facts")

    by_cell: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for f in facts:
        by_cell[(f["tool_canonical"], f["field"])].append(f)
    rank = {"high": 0, "medium": 1, "low": 2}

    changes: list[dict] = []
    new_cells: dict[tuple[str, str], tuple[str, str, str]] = {}
    for (tool, field), items in sorted(by_cell.items()):
        items = sorted(items, key=lambda d: (rank[d["confidence"]], d["value"]))
        values = dedupe([i["value"] for i in items])
        evidence = "; ".join(f"{i['evidence']}" for i in items[:MAX_EVIDENCE])
        prov = ("research:" + "; ".join(sorted({i["source_file"] for i in items}))
                + f" ({len(items)} facts)")
        best = items[0]["confidence"]
        # nothing is auto-CONFIRMED: a human promotes cells with promote_status.py
        status = tp.STATUS_CHECK
        old = next((r[field] for r in ti1 if r["tool_canonical"] == tool), "")
        new = join_cell(values)
        if old == tp.MISSING or not old:
            kind = "filled gap"
        elif " ".join(old.lower().split()) in {" ".join(v.lower().split())
                                               for v in values}:
            kind = "confirmed existing"
        else:
            kind = "replaced"
        new_cells[(tool, field)] = (new, f"{evidence} {prov}", status)
        if kind in ("filled gap", "replaced"):
            changes.append({"tool_canonical": tool, "field": field, "change": kind,
                            "previous_value": old,
                            "previous_evidence": next(
                                (r.get(f"{field}_evidence", "") for r in ti1
                                 if r["tool_canonical"] == tool), ""),
                            "new_value": new, "new_evidence": f"{evidence} {prov}",
                            "confidence": best,
                            "reason": " | ".join(dedupe(
                                [i["note"] for i in items if i["note"]])[:3])})

    n_gap = sum(1 for c in changes if c["change"] == "filled gap")
    n_rep = sum(1 for c in changes if c["change"] == "replaced")
    print(f"{len(by_cell)} (tool,field) cells from research: "
          f"{n_gap} gaps filled, {n_rep} cells corrected")

    # ------------------------------------------------------------------ write
    if args.dry_run:
        for c in changes[:25]:
            print(f"  [{c['change']}] {c['tool_canonical']}/{c['field']}: "
                  f"{c['previous_value'][:50]!r} -> {c['new_value'][:70]!r}")
        print("dry run: nothing written")
        return 0

    ti1_cols = list(ti1[0].keys())
    curated_path = tp.CURATED_DIR / "tool_inventory_curated.csv"
    previous = {}
    for r in tp.read_csv(curated_path):
        if r.get("value"):
            previous[(r["tool_canonical"], r["field"])] = (
                r["value"], r.get("evidence", ""), r.get("status", "NEEDS_CHECK"))
    bak = tp.backup(curated_path)
    if bak:
        print(f"previous curated -> {bak.name}")

    notes = [f for f in tp.read_csv(tp.RAW_DIR / "research_facts.csv")
             if f["field"] == "notes"]
    curated_rows = []
    for tool in sorted(meta):
        cat, mod, role = meta[tool]
        for field in CURATED_FIELDS:
            if (tool, field) in new_cells:
                value, evidence, status = new_cells[(tool, field)]
            elif (tool, field) in previous:
                value, evidence, status = previous[(tool, field)]
            else:
                row = next(r for r in ti1 if r["tool_canonical"] == tool)
                old = row[field]
                old_ev = row.get(f"{field}_evidence", "")
                if old == tp.MISSING or old_ev in TEMPLATE_STRINGS:
                    value, evidence, status = "", "", "TODO"
                else:
                    value, evidence, status = old, f"{old_ev} [auto-artefact]", \
                        "NEEDS_CHECK"
            curated_rows.append({
                "tool_canonical": tool, "category": cat, "modification": mod,
                "role": role, "field": field, "value": value, "evidence": evidence,
                "status": status, "how_to_obtain": ""})
    tp.write_csv(curated_path, curated_rows, CURATED_COLS)
    n_cur = sum(1 for r in curated_rows if r["value"])
    n_conf = sum(1 for r in curated_rows if r["status"] == tp.STATUS_CONFIRMED)
    print(f"curated: {len(curated_rows)} rows ({n_cur} carrying a value) "
          f"-> {curated_path.relative_to(tp.INV_ROOT)}")

    for row in ti1:
        tool = row["tool_canonical"]
        for field in tp.FIELDS:
            if field == "notes":
                rel = [f for f in notes if f["tool_canonical"] == tool]
                cell = " || ".join(f"{f['value']} [{f['evidence']}]" for f in rel[:8])
                row[field] = cell or tp.MISSING
                row[f"{field}_evidence"] = (
                    f"research facts: {', '.join(sorted({f['source_file'] for f in rel}))}"
                    if rel else "")
                continue
            if (tool, field) in new_cells:
                value, evidence, _ = new_cells[(tool, field)]
                row[field], row[f"{field}_evidence"] = value, evidence
            row[f"{field}_evidence"] = tp.normalise_evidence(
                row.get(f"{field}_evidence", ""))
        todo = [f for f in CURATED_FIELDS if row[f] == tp.MISSING]
        row["status"] = "complete" if not todo else \
            f"partial ({len(todo)} field(s) missing: {', '.join(todo)})"
    stamp = tp.stamp()
    for row in ti1:
        row["generated_at"] = stamp
    tp.backup(ti1_path)
    tp.write_csv(ti1_path, ti1, ti1_cols)
    print(f"TI1 rewritten -> {ti1_path.relative_to(tp.INV_ROOT)}")

    md_cols = (["tool_canonical", "category", "modification", "role"] +
               CURATED_FIELDS)
    with (tp.TABLE_DIR / "per_tool_implementation.md").open(
            "w", encoding="utf-8") as fh:
        fh.write("# TI1 -- per-tool implementation detail (R1-5)\n\n")
        fh.write(f"Generated {stamp}. Values carry provenance in the CSV's "
                 f"*{{field}}_evidence columns.\n\n")
        for row in ti1:
            fh.write(f"## {row['tool_canonical']} "
                     f"({row['category']} / {row['modification']} / {row['role']})\n\n")
            for field in md_cols[4:]:
                fh.write(f"- **{field}**: {row[field]}\n")
            for field in CURATED_FIELDS:
                if row.get(f"{field}_evidence"):
                    fh.write(f"  - evidence: {row[f'{field}_evidence']}\n")
            fh.write(f"- **status**: {row['status']}\n\n")

    tp.write_csv(tp.TABLE_DIR / "change_log.csv", changes,
                 ["tool_canonical", "field", "change", "previous_value",
                  "previous_evidence", "new_value", "new_evidence", "confidence",
                  "reason"])
    todo_rows = []
    for row in ti1:
        for field in CURATED_FIELDS:
            if row[field] == tp.MISSING:
                todo_rows.append({
                    "tool_canonical": row["tool_canonical"], "field": field,
                    "tried": "conda/pip/binary probes; yml/log/md scan; vendored "
                             "source scan (research/*.md); our own scripts "
                             "(commands_raw.csv)",
                    "note": next((n["note"] for n in notes
                                  if n["tool_canonical"] == row["tool_canonical"]),
                                 "")})
    tp.write_csv(tp.TABLE_DIR / "todo_report.csv", todo_rows,
                 ["tool_canonical", "field", "tried", "note"])
    covered = {r["tool_canonical"] for r in todo_rows}
    print(f"TI5: {len(todo_rows)} cells still unresolved across {len(covered)} tools")

    report = ["# TI8 -- R1-5 field coverage after research\n",
              f"Generated {stamp}.\n",
              "## By field\n", "| field | tools with a value | not recorded |",
              "|---|---|---|"]
    for field in CURATED_FIELDS:
        have = sum(1 for r in ti1 if r[field] != tp.MISSING)
        report.append(f"| {field} | {have} | {len(ti1) - have} |")
    report += ["\n## By tool\n", "| tool | resolved | open fields |", "|---|---|---|"]
    for row in sorted(ti1, key=lambda r: r["tool_canonical"]):
        open_fields = [f for f in CURATED_FIELDS if row[f] == tp.MISSING]
        report.append(f"| {row['tool_canonical']} | "
                      f"{len(CURATED_FIELDS) - len(open_fields)}/{len(CURATED_FIELDS)} "
                      f"| {', '.join(open_fields) or '-'} |")
    (tp.TABLE_DIR / "coverage_report.md").write_text("\n".join(report) + "\n",
                                                         encoding="utf-8")
    print("TI8 coverage report -> tables/coverage_report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
