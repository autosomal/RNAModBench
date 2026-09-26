#!/usr/bin/env python3
"""Shared paths + constants for the R1-5 gap-filling scripts.

The project was reorganised (see $RNAMODBENCH_ROOT/00_REORG_MAP.csv) after the
original collectors ran, so every path used by the legacy
``tool_inventory/scripts/ti_common.py`` is remapped here.  Nothing in this module
writes outside ``INV_ROOT``.
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
import csv
import re
import time
from pathlib import Path

PROJECT_ROOT = Path(str(_RB))
INV_ROOT = (_RB / "tools/inventory")
SCRIPTS_DIR = (_RB / "tools/inventory/scripts")
RAW_DIR = (_RB / "tools/inventory/raw")
CURATED_DIR = (_RB / "tools/inventory/curated")
TABLE_DIR = (_RB / "tools/inventory/tables")
LOG_DIR = (_RB / "tools/inventory/logs")
RESEARCH_DIR = (_RB / "tools/inventory/research")

#: live locations of the things the legacy collectors called code/ result / yaml
CODE_USER = (_XB / "tool_scripts")
CODE = (_XB / "code/code")
YAML_DIR = (_RB / "envs/as_run")
RESULT = (_XB / "raw/result")
RESULT_RNA004 = (_XB / "raw/result_RNA004")
SOURCE_ROOT = Path(str(_XB / "source_code/benchmark"))

#: prefix tokens the research notes use to keep lines short
PREFIXES: dict[str, str] = {
    "NM": str((_XB / "source_code/benchmark/nanom6A_2022_12_22")),
    "MA": str((_XB / "source_code/benchmark/m6anet-v-2.1.0")),
    "CH": str((_XB / "source_code/benchmark/CHEUI")),
    "DM": str((_XB / "source_code/benchmark/DRUMMER")),
    "EL": str((_XB / "source_code/benchmark/eligos2-v2.1.0")),
    "SC": str((_XB / "tool_scripts/detection")),
    "HA": str((_RB / "src/sites_v2/common/legacy_liftover.py")),
    "CM": str((_RB / "tools/inventory/raw/commands_raw.csv")),
    "CU/": str(CODE_USER),
    "CO/": str(CODE),
    "TI/": str(INV_ROOT),
    "RR": str((_XB / "raw")),
    # group E / F abbreviations
    "MK": str((_XB / "source_code/benchmark/mAFiA-0.1.0")),
    "MFI": str(_XB / "miniconda3/envs/mafia"),
    "TB": str(_XB / "miniconda3/envs/tombo/lib/python3.7/site-packages/tombo"),
    "PB": str((_XB / "tool_scripts/python_postprocessing")),
    "SV": str((_RB / "src/sites_v2")),
    "GU": str(_XB / "miniconda3/envs/guitar_asm/lib/R/library/Guitar"),
    "RS": str((_XB / "raw/result/_not_in_manuscript")),
    "CC": str((_XB / "raw/converted_callsets/_not_in_manuscript")),
    "CM": str((_RB / "tools/inventory/raw/commands_raw.csv")),
    "VW": str((_RB / "tools/inventory/raw/versions_raw.csv")),
    "R2": str(_XB / "source_code/nanopore/R2Dtool"),
    "GTR": str((_XB / "sites_v2/guitar_metagene")),
}

MISSING = "not recorded"

#: the ten reviewer-mandated R1-5 fields (same contract as ti_common.FIELDS)
FIELDS: list[str] = [
    "software_or_version",
    "model_checkpoint",
    "required_input",
    "min_read_coverage",
    "filtering_parameters",
    "calling_threshold",
    "multiple_testing_correction",
    "default_vs_optimised",
    "coordinate_harmonisation",
    "notes",
]

RESEARCH_FILES = [
    "groupA_detection_a.md",
    "groupB_detection_b.md",
    "groupC_detection_c.md",
    "groupD_upstream_and_harmonisation.md",
    "groupE_missing_detection_tools.md",
    "groupF_infrastructure.md",
]

STATUS_CONFIRMED = "CONFIRMED"
STATUS_CHECK = "NEEDS_CHECK"

#: old -> current absolute path fragments, applied to every evidence string
PATH_FIXUPS: list[tuple[str, str]] = [
    (str(_XB / "tool_scripts"), f"{CODE_USER}/"),
    (str(_XB / "raw/result_RNA004"), f"{RESULT_RNA004}/"),
    (str(_XB / "raw/result"), f"{RESULT}/"),
    (str(_RB / "envs/as_run"), f"{YAML_DIR}/"),
    (str(_XB / "code"), f"{CODE}/"),
    (str(_RB / "tools/inventory"), f"{INV_ROOT}/"),
    (str(_RB / "revision_output"),
     f"{PROJECT_ROOT}/04_revision_analysis/revision_output/"),
    (str(_XB / "archive/output"),
     f"{PROJECT_ROOT}/04_revision_analysis/output/"),
]


def normalise_evidence(text: str) -> str:
    """Expand research prefixes and repoint stale absolute paths at the live tree."""
    out = (text or "").strip()
    for tok, target in PREFIXES.items():
        if tok == "HA":
            out = re.sub(rf"\b{tok}(?=[\s:|]|$)", target, out)
        else:
            out = re.sub(rf"(?<![\w/.]){tok}/", f"{target}/", out)
            out = re.sub(rf"(?<![\w/.]){tok}(?=[\s:|]|$)", target, out)
    out = re.sub(r"\braw/(commands|params|versions|source|research)_raw\.csv",
                 lambda m: str(RAW_DIR / (m.group(1) + "_raw.csv")), out)
    for old, new in PATH_FIXUPS:
        if old != new:
            out = out.replace(old, new)
    return out


def read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", newline="") as fh:
        return [{k: (v or "") for k, v in row.items()}
                for row in csv.DictReader(fh)]


def write_csv(path: Path, rows: list[dict], cols: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def stamp() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def backup(path: Path) -> Path | None:
    """Copy ``path`` next to itself with a timestamp suffix (never overwrite)."""
    if not path.exists():
        return None
    dest = path.with_name(f"{path.stem}.bak_{time.strftime('%Y%m%d_%H%M%S')}{path.suffix}")
    dest.write_bytes(path.read_bytes())
    return dest
