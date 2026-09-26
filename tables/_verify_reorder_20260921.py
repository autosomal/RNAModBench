#!/usr/bin/env python3
"""Verify the 2026-09-21 reorder of Supplementary_Tables.pdf.

Old PDF order was alphabetical (S10, S11, S1 ... S9); the new one is numeric
(S1 ... S11).  This script asserts

* the new page order is 1, 2, 3, ... 11 by table number, and
* the normalised text of every table is byte-identical to the old build
  (page footers stripped, since page numbers legitimately shift).

Usage: conda run -n benchmark-revision --no-capture-output python _verify_reorder_20260921.py
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD = HERE / "Supplementary_Tables.pdf.bak_20260921_1046"
NEW = HERE / "Supplementary_Tables.pdf"


def page_count(pdf: Path) -> int:
    info = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True).stdout
    return int(re.search(r"Pages:\s+(\d+)", info).group(1))


def pages(pdf: Path) -> list[str]:
    out = []
    for p in range(1, page_count(pdf) + 1):
        txt = subprocess.run(["pdftotext", "-f", str(p), "-l", str(p), str(pdf), "-"],
                             capture_output=True, text=True).stdout
        # footer is a lone page number at the very end of the page text
        txt = re.sub(r"\n\s*\d+\s*\n?\s*$", "\n", txt)
        out.append(txt)
    return out


def blocks(pdf: Path) -> dict[str, str]:
    """Map 'S1' -> concatenated normalised text of every page of that table."""
    res: dict[str, list[str]] = {}
    cur: str | None = None
    for txt in pages(pdf):
        m = re.search(r"Table (S\d+):", txt)
        if m:
            cur = m.group(1)
        assert cur, "page before the first table heading"
        res.setdefault(cur, []).append(re.sub(r"\s+", " ", txt).strip())
    return {k: " ".join(v) for k, v in res.items()}


def order(pdf: Path) -> list[str]:
    seen: list[str] = []
    for txt in pages(pdf):
        m = re.search(r"Table (S\d+):", txt)
        if m and (not seen or seen[-1] != m.group(1)):
            seen.append(m.group(1))
    return seen


old_order, new_order = order(OLD), order(NEW)
expected = ["S%d" % i for i in range(1, 12)]
print("old order:", " ".join(old_order))
print("new order:", " ".join(new_order))
assert new_order == expected, "new order is not S1..S11"

old_b, new_b = blocks(OLD), blocks(NEW)
assert set(old_b) == set(new_b) == set(expected), (sorted(old_b), sorted(new_b))
bad = []
for k in expected:
    if old_b[k] != new_b[k]:
        bad.append(k)
        # show the first divergence to make the failure actionable
        a, b = old_b[k], new_b[k]
        i = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
        print("DIFF", k, "at char", i)
        print("  old:", a[max(0, i - 60):i + 60])
        print("  new:", b[max(0, i - 60):i + 60])
print("tables compared:", len(expected), "content mismatches:", len(bad))
sys.exit(1 if bad else 0)
