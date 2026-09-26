#!/usr/bin/env python3
"""Promote curated cells from NEEDS_CHECK to CONFIRMED after human verification.

The research layer proposes values with status ``NEEDS_CHECK``: the value and its
``file:line`` evidence are in the repository, but a human has not eyeballed every
citation yet.  Once a cell has been checked, promote it so that ``build_tables.py``
(and any downstream table) ranks it above every automatic probe.

Usage
-----
    python3 promote_status.py --list                     # what is still open
    python3 promote_status.py --tool xPore --field all   # promote verified cells
    python3 promote_status.py --tool xPore --field calling_threshold --revoke
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ti_paths as tp  # noqa: E402

CURATED = tp.CURATED_DIR / "tool_inventory_curated.csv"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--tool")
    ap.add_argument("--field")
    ap.add_argument("--revoke", action="store_true",
                    help="set the matched cells back to NEEDS_CHECK")
    args = ap.parse_args()

    rows = tp.read_csv(CURATED)
    if not rows:
        print("curated file missing", file=sys.stderr)
        return 1

    def matches(r: dict) -> bool:
        return ((not args.tool or r["tool_canonical"] == args.tool)
                and (not args.field or args.field == "all" or r["field"] == args.field))

    if args.list:
        counts = Counter((r["status"], bool(r["value"])) for r in rows)
        for (status, has_value), n in sorted(counts.items()):
            print(f"{status or '(blank)':14} value={has_value!s:5} {n:4}")
        open_cells = [f"{r['tool_canonical']}/{r['field']}" for r in rows
                      if r["status"] != tp.STATUS_CONFIRMED and matches(r)]
        print(f"\nnot yet CONFIRMED and matching filter: {len(open_cells)}")
        for c in open_cells[:60]:
            print("  ", c)
        return 0

    if not (args.tool or args.field):
        print("refusing to touch every cell: pass --tool and/or --field", file=sys.stderr)
        return 2
    target = tp.STATUS_CHECK if args.revoke else tp.STATUS_CONFIRMED
    changed = 0
    for r in rows:
        if matches(r) and r["value"] and r["status"] != target:
            r["status"] = target
            r["evidence"] = (r["evidence"] + f" promoted {tp.stamp()}").strip()
            changed += 1
    tp.backup(CURATED)
    tp.write_csv(CURATED, rows, list(rows[0].keys()))
    print(f"{changed} cells -> {target}; re-run apply_research.py to refresh TI1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
