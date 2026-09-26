#!/usr/bin/env python3
"""Collect run-time parameters, models and versions from the artefacts (R1-5).

Static (read-only) sources
--------------------------
* ``yaml/*.yml``                       -- conda environment definitions
* ``result/**/*.yml``                  -- xPore (readcount_min/max), CHEUI-diff
                                          (upper/lower_cutoff), m6Anet configs
* ``result/**/*.log``                  -- Nanocompore SampComp header
                                          (package_version, min_coverage,
                                          comparison_methods, ...),
                                          f5c-mode pipeline logs (--support)
* ``result_RNA004/documents/IMPORTANT_INFO.md`` -- Dorado models and version
* ``raw/commands_raw.csv``             -- parameters given on the command line

Only the *header* of big logs is read (parameters are logged at start-up).

Output: ``raw/params_raw.csv`` (RawFact columns)

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/tools/inventory/scripts/collect_params.py
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

#: artefact key (lower-cased) -> inventory field
KEY_TO_FIELD: dict[str, str] = {
    "package_version": "software_or_version",
    "version": "software_or_version",
    "software_version": "software_or_version",
    "min_coverage": "min_read_coverage",
    "mincov": "min_read_coverage",
    "readcount_min": "min_read_coverage",
    "readcount_max": "filtering_parameters",
    "coverage": "min_read_coverage",
    "support": "min_read_coverage",
    "upper_cutoff": "calling_threshold",
    "lower_cutoff": "calling_threshold",
    "cutoff": "calling_threshold",
    "threshold": "calling_threshold",
    "pval": "calling_threshold",
    "pvalue": "calling_threshold",
    "p_value": "calling_threshold",
    "qvalue": "calling_threshold",
    "comparison_methods": "multiple_testing_correction",
    "sequence_context": "filtering_parameters",
    "downsample_high_coverage": "filtering_parameters",
    "max_invalid_kmers_freq": "filtering_parameters",
    "min_ref_length": "filtering_parameters",
    "logit": "filtering_parameters",
    "model": "model_checkpoint",
    "model_path": "model_checkpoint",
    "checkpoint": "model_checkpoint",
    "nthreads": "filtering_parameters",
}

#: command-line flag -> inventory field.  Compute-only flags (threads, device)
#: are deliberately absent: they are not implementation parameters a reviewer
#: needs and would only drown the table in noise.
FLAG_TO_FIELD: dict[str, str] = {
    "--min-coverage": "min_read_coverage",
    "--coverage": "min_read_coverage",
    "--support": "min_read_coverage",
    "--min-reads": "min_read_coverage",
    "--readcount-min": "min_read_coverage",
    "--threshold": "calling_threshold",
    "--cutoff": "calling_threshold",
    "--pval": "calling_threshold",
    "--p-value": "calling_threshold",
    "--qval": "calling_threshold",
    "--fdr": "multiple_testing_correction",
    "--model": "model_checkpoint",
    "--modified-bases-models": "model_checkpoint",
}

#: fields whose value must look like a number (coverage / threshold / p-value)
NUMERIC_FIELDS = {"min_read_coverage", "calling_threshold"}

YML_KV = re.compile(r"^\s*([A-Za-z_][\w.\-]*)\s*:\s*(\S.*?)\s*$")
LOG_KV = re.compile(r"[|\s]\s*([a-z_][a-z0-9_]*)\s*:\s*([^\s|,]+)")
DORADO_MODEL = re.compile(r"rna00[24]_[0-9]+bps_(?:sup|hac)@v[\d.]+(?:_[A-Za-z0-9_]+@v\d+)?")
DORADO_BIN = re.compile(r"dorado-(\d+\.\d+(?:\.\d+)?)")
VERSION_IN_PATH = re.compile(r"/([A-Za-z][A-Za-z0-9_.]*?)[-_](\d+\.\d+(?:\.\d+)?)(?:/|$)")


def field_for_key(key: str) -> str | None:
    k = key.lower().strip()
    if k in KEY_TO_FIELD:
        return KEY_TO_FIELD[k]
    for pat, field in (("coverage", "min_read_coverage"),
                       ("cutoff", "calling_threshold"),
                       ("threshold", "calling_threshold"),
                       ("fdr", "multiple_testing_correction"),
                       ("model", "model_checkpoint"),
                       ("version", "software_or_version")):
        if pat in k:
            return field
    return None


def scan_yml(path: Path, logger) -> list[tc.RawFact]:
    facts: list[tc.RawFact] = []
    dataset = tc.dataset_from_path(path)
    tool = tc.canonical_tool(tool_hint(path))
    text = tc.read_text_head(path, logger)
    for i, line in enumerate(text.splitlines(), start=1):
        m = YML_KV.match(line)
        if not m:
            continue
        key, value = m.group(1), m.group(2).strip()
        # conda env definitions: dependency pins
        if key.lower() in {"dependencies", "channels", "name", "prefix"}:
            continue
        dep = re.match(r"^\s*-\s*([A-Za-z0-9_.\-]+)\s*([=<>!~]+)\s*([\d][^\s]*)", line)
        if dep and path.parent == tc.YAML_DIR:
            pkg, ver = dep.group(1), dep.group(3)
            cand = tc.canonical_tool(pkg)
            if cand in tc.TOOL_META or pkg.lower() in tc.ALIASES:
                facts.append(tc.RawFact(cand, dataset, "software_or_version",
                                        f"{pkg} {ver}", str(path), i, "static_yml"))
            continue
        field = field_for_key(key)
        if not field:
            continue
        facts.append(tc.RawFact(tool, dataset, field, value, str(path), i,
                                "static_yml"))
    return facts


def tool_hint(path: Path) -> str:
    """Best guess of the owning tool from the artefact path."""
    parts = path.parts
    if "result" in parts:
        idx = parts.index("result")
        if idx + 1 < len(parts):
            return parts[idx + 1]
    if path.parent == tc.YAML_DIR:
        return path.stem
    return path.stem


def scan_log(path: Path, logger) -> list[tc.RawFact]:
    """Parse the header of a log file (first 400 lines)."""
    facts: list[tc.RawFact] = []
    dataset = tc.dataset_from_path(path)
    tool = tc.canonical_tool(tool_hint(path))
    text = tc.read_text_head(path, logger)
    lines = text.splitlines()[:400]
    seen: set[str] = set()
    for i, line in enumerate(lines, start=1):
        if len(line) > 4000:
            continue
        m = LOG_KV.search(line)
        if not m:
            continue
        key, value = m.group(1), m.group(2).strip().rstrip(",")
        field = field_for_key(key)
        if not field or key in seen:
            continue
        seen.add(key)
        facts.append(tc.RawFact(tool, dataset, field, value, str(path), i,
                                "static_log"))
    return facts


def scan_md(path: Path, logger) -> list[tc.RawFact]:
    """Dorado models / version from result_RNA004/documents/IMPORTANT_INFO.md."""
    facts: list[tc.RawFact] = []
    text = tc.read_text_head(path, logger)
    models = sorted(set(DORADO_MODEL.findall(text)))
    if models:
        facts.append(tc.RawFact("Dorado", "RNA004", "model_checkpoint",
                                "; ".join(models), str(path), 1, "static_md"))
    vers = sorted(set(DORADO_BIN.findall(text)))
    if vers:
        facts.append(tc.RawFact("Dorado", "RNA004", "software_or_version",
                                "/".join(vers), str(path), 1, "static_md"))
    # explicit command blocks
    for i, line in enumerate(text.splitlines(), start=1):
        if line.strip().startswith(str(_XB / "tool")):
            facts.append(tc.RawFact("Dorado", "RNA004", "command_line",
                                    line.strip(), str(path), i, "static_md"))
    return facts


def _plausible(field: str, value: str) -> bool:
    """Reject values that are obviously not a parameter (paths, shell noise)."""
    v = value.strip()
    if not v or v in {"-", "--"}:
        return False
    if field in NUMERIC_FIELDS:
        return bool(re.fullmatch(r"\d+(\.\d+)?", v))
    if field in {"required_input", "notes"}:
        return len(v) <= 200
    if field in {"model_checkpoint", "command_line"}:
        return len(v) <= 1000
    # filtering / multiple-testing / version: short tokens only
    return len(v) <= 40 and "/" not in v


def scan_commands(path: Path, logger) -> list[tc.RawFact]:
    """Pull parameter values out of the collected command lines."""
    facts: list[tc.RawFact] = []
    if not path.exists():
        logger.warning("commands_raw.csv missing - run collect_commands.py first")
        return facts
    df = pd.read_csv(path)
    for _, r in df.iterrows():
        cmd = str(r["command"])
        for flag, field in FLAG_TO_FIELD.items():
            m = re.search(re.escape(flag) + r"[=\s]+([^\s|;&]+)", cmd)
            if not m or not _plausible(field, m.group(1)):
                continue
            facts.append(tc.RawFact(str(r["tool_canonical"]), str(r["dataset"]),
                                    field, m.group(1), str(r["source"]),
                                    int(r["line_no"]), "static_script"))
        # version encoded in a source path, e.g. .../NanoNm-1.0.0/predict_sites...
        for m in VERSION_IN_PATH.finditer(cmd):
            cand = tc.canonical_tool(m.group(1))
            if cand in tc.TOOL_META:
                facts.append(tc.RawFact(cand, str(r["dataset"]),
                                        "software_or_version", m.group(2),
                                        str(r["source"]), int(r["line_no"]),
                                        "static_script"))
    return facts


def main() -> None:
    logger = tc.setup_logger("collect_params")
    facts: list[tc.RawFact] = []

    for p in sorted(tc.YAML_DIR.glob("*.y*ml")):
        facts += scan_yml(p, logger)
    logger.info("yaml definitions: %d facts", len(facts))

    n0 = len(facts)
    for p in sorted(tc.RESULT.rglob("*.y*ml")):
        if p.suffix.lower() in tc.SKIP_SUFFIXES:
            continue
        facts += scan_yml(p, logger)
    logger.info("result/**/*.yml: %d facts", len(facts) - n0)

    n0 = len(facts)
    n_log = 0
    for p in sorted(tc.RESULT.rglob("*.log")):
        facts += scan_log(p, logger)
        n_log += 1
    logger.info("result/**/*.log: %d files -> %d facts", n_log, len(facts) - n0)

    info = tc.RNA004_ROOT / "documents" / "IMPORTANT_INFO.md"
    if info.exists():
        n0 = len(facts)
        facts += scan_md(info, logger)
        logger.info("IMPORTANT_INFO.md: %d facts", len(facts) - n0)

    n0 = len(facts)
    facts += scan_commands(tc.RAW_DIR / "commands_raw.csv", logger)
    logger.info("command-line parameters: %d facts", len(facts) - n0)

    facts = [f for f in facts
             if f.value and f.value.lower() != "none"
             and _plausible(f.field, f.value)]
    n = tc.write_facts(facts, tc.RAW_DIR / "params_raw.csv", generated_at=tc.stamp())
    logger.info("written %d facts -> %s", n, tc.RAW_DIR / "params_raw.csv")


if __name__ == "__main__":
    main()
