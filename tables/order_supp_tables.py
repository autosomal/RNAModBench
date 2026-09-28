#!/usr/bin/env python3
"""Put every supplementary table in the order the manuscript uses (2026-09-27).

The built tables inherited whatever order their source files happened to have:
S1 and S11 a curated order, S6 the figure order with later rows appended, S7
alphabetical, and the tables grouped by species started with Human in two cases
and with Arabidopsis in others.  A supplement reads as one document, so this
script rewrites the row order of every table by the two orders the manuscript
itself fixes:

* **species** -- the order of the introductory sentence and of the figure code
  (`harmonisation/common/config.py`, `SPECIES_ORDER`): Arabidopsis, mouse, HeLa,
  E. coli, Curlcake;
* **tools** -- the configuration order of the main figures (`TOOL_ORDER` of
  `harmonisation/scripts/23_guitar_metagene.py` and its siblings, mirrored by
  `ARTICLE_M6A_TOOLS`): the thirteen m6A configurations first, then the
  non-m6A tools in the order the non-m6A results present them (Fig. 7), then
  the RNA004 built-in modification models (Fig. 8), then the tools that appear
  only in the Supporting Information, alphabetically.

Sorting is **stable**: inside one species or one tool, the rows keep the order
their own source gave them (pairs of S3, ranks of S4, strata of S12, and so on).
No value is touched -- only the sequence of rows.

Tables whose first column is neither a species nor a tool are left alone: S2 is
one tool per row but sorted by detected sites as in Fig. 1B, S9 is one row per
modification type in the order the boundary is discussed.

Usage (after the builders, before `render_supp_tables.py`)
---------------------------------------------------------
conda run -n benchmark-revision --no-capture-output python order_supp_tables.py
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TAB = HERE / "tables"

#: the order of the introductory sentence and of `SPECIES_ORDER` in the figure code
SPECIES = ["arabidopsis", "mouse", "human", "escherichia", "curlcake"]

#: the thirteen m6A configurations, in the order of the main figures' `TOOL_ORDER`
M6A13 = ["cheui_m6a", "dena", "drummer", "eligos2_diff", "eligos2_solo",
         "epinano_error", "m6anet", "mines", "nanocompore", "nanom6a",
         "nanospa_m6a", "xpore", "yanocomp"]

#: the non-m6A tools, in the order the non-m6A results present them (Fig. 7)
NONM6A = ["nanopsu", "nanospa_psi", "nanomud_psi", "cheui_m5c", "nanonm",
          "nanomud_m1psi"]

#: rows that name one configuration *of* a tool above: they stay next to it
NEIGHBOURS = {"cheui_diff": "cheui_m6a", "epinano_svm": "epinano_error",
              "nanomud": "nanomud_psi",
              "eligos2": "eligos2_diff", "epinano": "epinano_error",
              "cheui": "cheui_m6a", "nanospa": "nanospa_m6a"}

SPECIES_ALIAS = {
    "arabidopsis": "arabidopsis", "arabidopsis thaliana": "arabidopsis",
    "mouse": "mouse", "mus musculus": "mouse", "mef": "mouse",
    "human": "human", "hela": "human", "homo sapiens": "human",
    "e.coli": "escherichia", "e. coli": "escherichia",
    "escherichia coli": "escherichia", "escherichia": "escherichia",
    "curlcake": "curlcake",
}


def tool_key(name: str) -> str:
    """Normalise a tool label of any table to the key used by the orders."""
    s = name.strip().lower()
    s = re.sub(r"\(.*?\)", " ", s)                     # "Dorado (RNA004 ...)"
    s = s.replace("-", "_").replace(" ", "_")
    s = re.sub(r"_+", "_", s).strip("_")
    return NEIGHBOURS.get(s, s)


def tool_rank(name: str) -> tuple:
    key = tool_key(name)
    if key in M6A13:
        return (0, M6A13.index(key), key)
    for base, near in NEIGHBOURS.items():
        if key == base and near in M6A13:              # e.g. CHEUI-diff
            return (0, M6A13.index(near) + 0.5, key)
    if key in NONM6A:
        return (1, NONM6A.index(key), key)
    if key.startswith("dorado"):
        return (2, 0, key)
    return (3, 0, key)                                 # SI-only tools, alphabetical


def species_rank(name: str) -> tuple:
    key = SPECIES_ALIAS.get(name.strip().lower(), name.strip().lower())
    return (SPECIES.index(key) if key in SPECIES else len(SPECIES), key)


def rewrite(path: Path, col: int, rank) -> int:
    """Stable-sort the data rows of one table, keeping its line layout."""
    lines = path.read_text(encoding="utf-8").split("\n")
    body = [ln for ln in lines if ln.strip() and not ln.startswith("#")]
    header, data = body[0], body[1:]
    first = lines.index(header)                       #: comments above the header
    preamble = [ln for ln in lines[:first] if ln.strip()]
    last = max(lines.index(r) for r in data)
    epilogue = [ln for ln in lines[last + 1:] if ln.strip()]
    before = [r.split("\t")[col] for r in data]
    data.sort(key=lambda r: rank(r.split("\t")[col]))  #: stable
    after = [r.split("\t")[col] for r in data]
    if before == after:
        return 0
    path.write_text("\n".join(preamble + [header] + data + epilogue) + "\n",
                    encoding="utf-8")
    return sum(1 for a, b in zip(before, after) if a != b)


#: table -> (column index of the sort key, the key function)
PLAN = {
    "TableS1_tools": (0, tool_rank),
    "TableS3_shared_universe_ratio": (0, species_rank),
    "TableS4_top5_motifs": (0, species_rank),
    "TableS5_non_m6a": (0, tool_rank),
    "TableS8_glori_matching": (0, species_rank),
    "TableS6_per_tool_implementation": (0, tool_rank),
    "TableS7_chemistry_matrix": (0, tool_rank),
    "TableS10_datasets": (0, species_rank),
    "TableS11_nomenclature": (0, tool_rank),
    "TableS12_coverage_rank_stability": (0, species_rank),
}


def main() -> None:
    moved = 0
    for stem, (col, rank) in PLAN.items():
        path = TAB / f"{stem}.tsv"
        if not path.exists():
            print(f"  {stem}: missing, skipped")
            continue
        n = rewrite(path, col, rank)
        moved += n
        print(f"  {stem}: {n} row(s) in a new position")

    #: gate -- every table now starts in the manuscript's order
    bad = []
    for stem, (col, rank) in PLAN.items():
        path = TAB / f"{stem}.tsv"
        if not path.exists():
            continue
        body = [ln for ln in path.read_text(encoding="utf-8").split("\n")
                if ln.strip() and not ln.startswith("#")]
        keys = [rank(r.split("\t")[col]) for r in body[1:]]
        if keys != sorted(keys):
            bad.append(stem)
    if bad:
        print("NOT IN ORDER:", ", ".join(bad))
        sys.exit(1)
    print(f"gate: {len(PLAN)} tables follow the manuscript's species and tool order "
          f"({moved} row(s) moved in total)")


if __name__ == "__main__":
    main()
