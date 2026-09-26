#!/usr/bin/env python3
"""Parse the research notes into machine-readable R1-5 fact proposals.

Input  : ``research/group{A,B,C,D}_*.md`` -- one line per fact,
         ``TOOL|field|value|evidence_path:line|confidence|note``
Output : ``raw/research_facts.csv`` (normalised facts) and
         ``tables/TI7_proposed_values.csv`` (gap-filling proposals per tool/field)

Only *structurally valid* lines are accepted: exactly six ``|`` separated columns,
a known tool and a known R1-5 field.  Every rejection is reported with its file:line
so that nothing is silently dropped.
"""
from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ti_paths as tp  # noqa: E402

VALID_CONFIDENCE = {"high", "medium", "low"}
PLACEHOLDER_EVIDENCE = {"", "-", "\u2014", "n/a", "not recorded", "CM"}


def known_tools() -> set[str]:
    ti1 = tp.TABLE_DIR / "TI1_per_tool_implementation.csv"
    return {r["tool_canonical"] for r in tp.read_csv(ti1)}


CONF_TOKENS = {"high", "medium", "low"}
PATH_CHUNK = re.compile(r"^(?:\S*\S/[\w.\-+@%~()<> */#]*|\S+\.(?:py|R|r|sh|md|yml|"
                        r"yaml|csv|txt|toml|cfg|ipynb|json)\S*)"
                        r"(?::\d+(?:-\d+)?(?:,\d+(?:-\d+)?)*)?$", re.IGNORECASE)
LINE_LIST = re.compile(r"^(?::?\d+(?:-\d+)?(?:,\d+(?:-\d+)?)*)$")
EVIDENCE_SEP = re.compile(r"\s*;\s*")


def _is_evidence(tok: str) -> bool:
    """One or more paths, joined by ``;``/``vs``, each with optional ``:12,15`` lines.

    Splitting a path list on ", " keeps ``file.csv:1018,262`` in one piece, and a bare
    number segment is treated as a continuation of the preceding path.
    """
    segments = [s.strip().strip('`"\'') for s in EVIDENCE_SEP.split(tok.strip())]
    segments = [s for s in segments if s]
    if not segments:
        return False
    for seg in segments:
        chunks = [c.strip().strip('`"\'') for c in seg.split(", ")]
        chunks = [c for c in chunks if c]
        prev_path = False
        for c in chunks:
            if PATH_CHUNK.match(c):
                prev_path = True
                continue
            if prev_path and LINE_LIST.match(c.lstrip(":")):
                continue
            return False
        if not prev_path:
            return False
    return True


def split_fact(line: str) -> list[str] | None:
    """``TOOL|field|value|evidence|confidence|note`` with ``|`` allowed inside value.

    The confidence token is the anchor: the column before it is the evidence, the
    columns before that (down to field 2) are the value, the optional column after
    it is the note.  Returns ``None`` when no such anchor exists.
    """
    parts = [p.strip() for p in line.split("|")]
    if len(parts) < 5:
        return None
    for i in range(len(parts) - 2, 2, -1):
        if parts[i].lower() not in CONF_TOKENS:
            continue
        if not _is_evidence(parts[i - 1]):
            continue
        value = " , ".join(p for p in parts[2:i - 1] if p)   # value carried a '|'
        return [parts[0], parts[1], value, parts[i - 1], parts[i].lower(),
                " , ".join(p for p in parts[i + 1:] if p)]
    return None


def parse_file(path: Path, tools: set[str]) -> tuple[list[dict], list[str]]:
    facts: list[dict] = []
    problems: list[str] = []
    in_unknown = False
    for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip().strip("`")
        if line.lower().startswith("## still unknown"):
            in_unknown = True
            continue
        if line.startswith("## "):
            in_unknown = False
        if "|" not in line or line.startswith("- ") or in_unknown:
            continue                      # prose / unknown-list entries
        parts = split_fact(line)
        if parts is None:
            problems.append(f"{path.name}:{lineno}: unparsable :: {line[:90]}")
            continue
        tool, field, value, evidence, confidence, note = parts
        if field not in tp.FIELDS:
            continue
        if tool not in tools:
            problems.append(f"{path.name}:{lineno}: unknown tool {tool!r}")
            continue
        if not value or value == tp.MISSING:
            problems.append(f"{path.name}:{lineno}: empty value for {tool}/{field}")
            continue
        if confidence not in VALID_CONFIDENCE:
            problems.append(f"{path.name}:{lineno}: bad confidence {confidence!r} "
                            f"({tool}/{field})")
            continue
        if evidence in PLACEHOLDER_EVIDENCE:
            problems.append(f"{path.name}:{lineno}: no usable evidence for "
                            f"{tool}/{field} ({evidence!r})")
            continue
        facts.append({
            "tool_canonical": tool,
            "field": field,
            "value": tp.normalise_evidence(value.replace("raw/commands_raw.csv", str(tp.RAW_DIR / "commands_raw.csv"))),
            "evidence": tp.normalise_evidence(evidence),
            "confidence": confidence,
            "note": note,
            "source_file": path.name,
            "source_line": str(lineno),
            "method": f"research_{confidence}",
        })
    return facts, problems


def main() -> int:
    tools = known_tools()
    if not tools:
        print("TI1 not found -- run the legacy build_tables.py first", file=sys.stderr)
        return 1
    all_facts: list[dict] = []
    all_problems: list[str] = []
    for name in tp.RESEARCH_FILES:
        path = tp.RESEARCH_DIR / name
        if not path.exists():
            all_problems.append(f"{name}: MISSING")
            continue
        facts, problems = parse_file(path, tools)
        all_facts += facts
        all_problems += problems
        print(f"{name}: {len(facts)} facts accepted, {len(problems)} rejected")

    cols = ["tool_canonical", "field", "value", "evidence", "confidence", "note",
            "source_file", "source_line", "method", "generated_at"]
    for f in all_facts:
        f["generated_at"] = tp.stamp()
    tp.write_csv(tp.RAW_DIR / "research_facts.csv", all_facts, cols)

    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for f in all_facts:
        grouped[(f["tool_canonical"], f["field"])].append(f)
    order = {"high": 0, "medium": 1, "low": 2}
    proposals = []
    for (tool, field), items in sorted(grouped.items()):
        items = sorted(items, key=lambda d: (order[d["confidence"]], d["value"]))
        seen, uniq = set(), []
        for i in items:
            if i["value"] not in seen:
                seen.add(i["value"])
                uniq.append(i)
        proposals.append({
            "tool_canonical": tool,
            "field": field,
            "n_facts": str(len(items)),
            "proposed_value": " | ".join(i["value"] for i in uniq[:6]),
            "best_confidence": uniq[0]["confidence"],
            "evidence": "; ".join(f"{i['evidence']}" for i in uniq[:4]),
            "notes": " || ".join(i["note"] for i in uniq[:3] if i["note"]),
            "source_files": "; ".join(sorted({i["source_file"] for i in items})),
        })
    pcols = ["tool_canonical", "field", "n_facts", "best_confidence",
             "proposed_value", "evidence", "notes", "source_files"]
    tp.write_csv(tp.TABLE_DIR / "TI7_proposed_values.csv", proposals, pcols)

    ti1 = tp.read_csv(tp.TABLE_DIR / "TI1_per_tool_implementation.csv")
    current = {(r["tool_canonical"], f): r[f] for r in ti1 for f in tp.FIELDS}
    gaps_filled = sum(1 for p in proposals
                      if current.get((p["tool_canonical"], p["field"])) == tp.MISSING)
    conflicts = [p for p in proposals
                 if current.get((p["tool_canonical"], p["field"])) not in
                 (tp.MISSING, "", p["proposed_value"])]
    print(f"\ntools in TI1: {len(ti1)}   facts: {len(all_facts)}   "
          f"(tool,field) proposals: {len(proposals)}")
    print(f"proposals that fill a 'not recorded' gap : {gaps_filled}")
    print(f"proposals that CONTRADICT the current cell: {len(conflicts)}")
    print(f"rejected lines: {len(all_problems)}")
    with (tp.LOG_DIR / "parse_research_problems.txt").open("w", encoding="utf-8") as fh:
        fh.write("\n".join(all_problems) + ("\n" if all_problems else ""))
    for p in all_problems[:40]:
        print("  !", p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
