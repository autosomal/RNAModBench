#!/usr/bin/env python3
"""Table S11: tool nomenclature (standardized names, modes, software versions)."""

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
import csv
import re
from pathlib import Path

RA = Path(str(_RB / "analysis"))
OUT = Path(__file__).resolve().parent / "tables" / "TableS11_nomenclature.tsv"

modes = {}
with (Path(__file__).resolve().parent / "tables" / "TableS1_tools.tsv").open() as fh:
    rows = [l for l in fh if not l.startswith("#")]
    for r in csv.DictReader(rows, delimiter="\t"):
        modes[r["Tool"]] = r["Modes evaluated in this study"]

ver = {}
with ((_RB / "analysis/supp_table_inputs/per_tool_implementation.csv")).open() as fh:
    for r in csv.DictReader(fh):
        name = (r["tool_display"] or r["tool"]).split(" /")[0].strip()
        s = re.sub(r"\s+", " ", r["software_or_version"]).strip()
        m = re.match(r"([A-Za-z][^|;]{0,70})", s)
        ver[name] = (m.group(1).strip() if m else s[:60]) or "not recorded"

ORDER = ["xPore", "Nanom6A", "ELIGOS2_solo", "ELIGOS2_diff", "DRUMMER", "Nanocompore",
         "EpiNano_Error", "m6Anet", "MINES", "DENA", "Yanocomp", "CHEUI_m6A", "CHEUI_m5C",
         "NanoSPA_m6A", "NanoSPA_Psi", "NanoMUD_Psi", "NanoMUD_m1Psi", "NanoPsu", "NanoNm"]
BASE = {"ELIGOS2_solo": "ELIGOS2", "ELIGOS2_diff": "ELIGOS2", "CHEUI_m6A": "CHEUI",
        "CHEUI_m5C": "CHEUI", "NanoSPA_m6A": "NanoSPA", "NanoSPA_Psi": "NanoSPA",
        "NanoMUD_Psi": "NanoMUD", "NanoMUD_m1Psi": "NanoMUD"}

with OUT.open("w") as fh:
    fh.write("# Table S11: Tool nomenclature used throughout the manuscript, figures and "
             "Supporting Information.\n")
    fh.write("Standardized name\tTool\tConfiguration evaluated\tSoftware version / model\n")
    for name in ORDER:
        base = BASE.get(name, name)
        fh.write("\t".join([name, base, modes.get(base, ""), ver.get(base, "not recorded")]) + "\n")
    fh.write("# note: names are standardized as in the manuscript; the reference for every tool "
             "is given in the numbered reference list of the manuscript. Technical parameters and "
             "the coordinate harmonisation are in Table S6.\n")
print("TableS11 written")
