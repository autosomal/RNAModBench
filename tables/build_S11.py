#!/usr/bin/env python3
"""Table S11: tool nomenclature (standardized names, modes, configuration evaluated)."""
import csv
from pathlib import Path

OUT = Path(__file__).resolve().parent / "tables" / "TableS11_nomenclature.tsv"

modes = {}
with (Path(__file__).resolve().parent / "tables" / "TableS1_tools.tsv").open() as fh:
    rows = [l for l in fh if not l.startswith("#")]
    for r in csv.DictReader(rows, delimiter="\t"):
        modes[r["Tool"]] = r["Modes evaluated in this study"]


#: built by looking the tool up in `NA4_per_tool_implementation.csv` and cutting the
#: string with `([A-Za-z][^|;]{0,70})`.  The lookup key was the first segment of
#: `tool_display` while the table is keyed by `CHEUI_m6A` / `ELIGOS2_diff` /
#: `ELIGOS2_solo`, so those four rows fell through to "not recorded" and contradicted
#: Table S6, and the 71-character cut left two long cells ending mid-word (Nanom6A,
#: MINES).  The column then became a pointer to Table S6.


#: A column whose every row repeats one pointer carries no information, so it is
#: **dropped**: the table is a nomenclature table (standardized name, tool,
#: configuration evaluated) and the note keeps the single pointer to Table S6, which
#: is where versions, checkpoints and parameters live.  The manuscript names the same
#: two columns ("the standardized name and the configuration evaluated"), so text and
#: table stay aligned.

ORDER = ["xPore", "Nanom6A", "ELIGOS2_solo", "ELIGOS2_diff", "DRUMMER", "Nanocompore",
         "EpiNano_Error", "m6Anet", "MINES", "DENA", "Yanocomp", "CHEUI_m6A", "CHEUI_m5C",
         "NanoSPA_m6A", "NanoSPA_\u03a8", "NanoMUD_\u03a8", "NanoMUD_m1\u03a8", "NanoPsu",
         "NanoNm"]
BASE = {"ELIGOS2_solo": "ELIGOS2", "ELIGOS2_diff": "ELIGOS2", "CHEUI_m6A": "CHEUI",
        "CHEUI_m5C": "CHEUI", "NanoSPA_m6A": "NanoSPA", "NanoSPA_\u03a8": "NanoSPA",
        "NanoMUD_\u03a8": "NanoMUD", "NanoMUD_m1\u03a8": "NanoMUD"}

with OUT.open("w") as fh:
    fh.write("# Table S11: Tool nomenclature used throughout the manuscript, figures and "
             "Supporting Information.\n")
    fh.write("Standardized name\tTool\tConfiguration evaluated\n")
    for name in ORDER:
        base = BASE.get(name, name)
        fh.write("\t".join([name, base, modes.get(base, "")]) + "\n")
    fh.write("# note: names are standardized as in the manuscript; the reference for every tool "
             "is given in the numbered reference list of the manuscript. The software version, "
             "model checkpoint, input requirements, coverage floor, filtering parameters, "
             "thresholds, multiple-testing correction and coordinate harmonisation of every "
             "configuration are in Table S6.\n")
print("TableS11 written")
