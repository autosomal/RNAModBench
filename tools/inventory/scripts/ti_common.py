"""Shared constants and helpers for the tool-inventory collectors.

Everything here is read-only and dependency-light: the collectors must be
re-runnable from any conda env that has ``pandas``.

Conventions
-----------
* ``RawFact`` is the single row contract produced by every collector, so that
  ``build_tables.py`` can merge the three layers with one code path.
* Tool names are normalised through :func:`canonical_tool`, which implements the
  R1-7 naming standardisation used elsewhere in the revision
  (xPore / Nanom6A / ELIGOS2_diff / NanoSPA-Psi).
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
import logging
import re
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
PROJECT_ROOT = Path(str(_RB))
INV_ROOT = (_RB / "tools/inventory")
SCRIPTS_DIR = (_RB / "tools/inventory/scripts")
RAW_DIR = (_RB / "tools/inventory/raw")
CURATED_DIR = (_RB / "tools/inventory/curated")
TABLE_DIR = (_RB / "tools/inventory/tables")
LOG_DIR = (_RB / "tools/inventory/logs")

#: script trees scanned for exact command lines (R1-10)
SCRIPT_ROOTS: list[Path] = [
    (_XB / "tool_scripts/detection"),
    (_XB / "tool_scripts/basecalling"),
    (_XB / "tool_scripts/shell_postprocess"),
    (_XB / "code/rerun"),
    (_XB / "code/f5c_mode"),
    (_XB / "raw/result_RNA004/scripts"),
]

RESULT = (_XB / "raw/result")
RNA004_ROOT = (_XB / "raw/result_RNA004")
YAML_DIR = (_RB / "envs/as_run")

#: never read these (binary / huge)
SKIP_SUFFIXES = {
    ".bam", ".bai", ".cram", ".hdf5", ".h5", ".pod5", ".blow5", ".fast5",
    ".gz", ".zip", ".npz", ".pt", ".pkl", ".pdf", ".png", ".jpg", ".pptx",
    ".gzi", ".index", ".mmi", ".so", ".pyc",
}
MAX_TEXT_BYTES = 20 * 1024 * 1024     # skip files above this
HEAD_BYTES = 200 * 1024               # for big logs: parameters live in the header

# --------------------------------------------------------------------------- #
# R1-5 : the ten reviewer-mandated fields
# --------------------------------------------------------------------------- #
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

#: how a missing value should be chased down (used by the curated template)
HOW_TO_OBTAIN: dict[str, str] = {
    "software_or_version":
        "conda activate <env> && conda list | grep -i <tool> <tool> --version"
 "GitHub release / commit hash",
    "model_checkpoint":
        " models/ ",
    "required_input":
        " TI2 and README",
    "min_read_coverage":
        " --min-coverage/--support/readcount_min min_coverage",
    "filtering_parameters":
        " and *.yml / *.toml",
    "calling_threshold":
        " --threshold/-p/--pval",
    "multiple_testing_correction":
        "/ BH Bonferroni qvalue / FDR / adj.P.Val",
    "default_vs_optimised":
        " optimised/changed",
    "coordinate_harmonisation":
        "tool_scripts/*_postprocessing transcript→genomic ",
    "notes":
        " ",
}

MISSING = "not recorded"

# --------------------------------------------------------------------------- #
# Tool catalogue
# --------------------------------------------------------------------------- #
#: canonical tool name -> (category, modification, role in the benchmark)
TOOL_META: dict[str, tuple[str, str, str]] = {
    # ---- m6A / modification detection tools (RNA002) ----------------------
    "Nanom6A": ("detection", "m6A", "de novo"),
    "m6Anet": ("detection", "m6A", "de novo"),
    "CHEUI_m6A": ("detection", "m6A", "de novo"),
    "CHEUI_m5C": ("detection", "m5C", "de novo"),
    "CHEUI-diff": ("detection", "m6A/m5C", "comparative"),
    "DENA": ("detection", "m6A", "comparative"),
    "DRUMMER": ("detection", "m6A", "comparative"),
    "ELIGOS2_diff": ("detection", "m6A (any)", "comparative"),
    "ELIGOS2_solo": ("detection", "m6A (any)", "de novo"),
    "EpiNano_Error": ("detection", "m6A", "comparative"),
    "EpiNano_SVM": ("detection", "m6A", "de novo"),
    "differr": ("detection", "m6A", "comparative"),
    "MINES": ("detection", "m6A", "de novo"),
    "Nanocompore": ("detection", "any", "comparative"),
    "NanoMUD": ("detection", "Psi/m1Psi", "de novo"),
    "NanoNm": ("detection", "Nm", "de novo"),
    "NanoPsu": ("detection", "Psi", "de novo"),
    "NanoSPA": ("detection", "m6A/Psi", "de novo"),
    "Tombo": ("detection", "any", "de novo"),
    "Tombo_com": ("detection", "any", "comparative"),
    "xPore": ("detection", "m6A", "comparative"),
    "Yanocomp": ("detection", "any", "comparative"),
    "mAFiA": ("detection", "m6A", "de novo"),
    "SingleMod": ("detection", "m6A", "de novo"),
    "TandemMod": ("detection", "m6A/m5C/...", "de novo"),
    "nanodoc2": ("detection", "any", "comparative"),
    "nanoRMS": ("detection", "any", "comparative"),
    "m1a-prediction": ("detection", "m1A", "de novo"),
    "penguin": ("detection", "Psi", "de novo"),
    "Dorado": ("detection", "m6A/pseU/m5C/inosine", "de novo"),
    # ---- upstream / infrastructure ----------------------------------------
    "guppy": ("basecalling", "n/a", "upstream"),
    "minimap2": ("alignment", "n/a", "upstream"),
    "samtools": ("file-handling", "n/a", "upstream"),
    "nanopolish": ("signal", "n/a", "upstream"),
    "f5c": ("signal", "n/a", "upstream"),
    "modkit": ("file-handling", "n/a", "upstream"),
    "slow5tools": ("file-handling", "n/a", "upstream"),
    "gffread": ("annotation", "n/a", "upstream"),
    "bedtools": ("file-handling", "n/a", "upstream"),
    "seqkit": ("file-handling", "n/a", "upstream"),
    "R2Dtool": ("annotation", "n/a", "upstream"),
    "Guitar": ("annotation", "n/a", "upstream"),
}

#: lower-cased token (script stem, source dir, executable) -> canonical name
ALIASES: dict[str, str] = {
    # detection tools
    "nanom6a": "Nanom6A", "nanom6a.sh": "Nanom6A", "nanom6a1": "Nanom6A",
    "m6anet": "m6Anet", "m6anet.sh": "m6Anet",
    "cheui": "CHEUI_m6A", "cheui_m6a": "CHEUI_m6A", "cheui_m5c": "CHEUI_m5C",
    "cheui-diff": "CHEUI-diff", "cheuidiff": "CHEUI-diff",
    "dena": "DENA", "drummer": "DRUMMER",
    "eligos2_diff": "ELIGOS2_diff", "eligos_diff": "ELIGOS2_diff",
    "eligos2_solo": "ELIGOS2_solo", "eligos_solo": "ELIGOS2_solo",
    "eligos2": "ELIGOS2_diff", "eligos": "ELIGOS2_diff",
    "epinano_differr": "EpiNano_Error", "epinano_error": "EpiNano_Error",
    "epinanoerror": "EpiNano_Error",
    "epinano_svm": "EpiNano_SVM", "epinanosvm": "EpiNano_SVM",
    "epinano": "EpiNano_Error", "singlemod": "SingleMod",
    "differr": "differr", "mines": "MINES",
    "nanocompore": "Nanocompore", "nanomud": "NanoMUD",
    "nanonm": "NanoNm", "nanonm-1.0.0": "NanoNm", "nanopsu": "NanoPsu",
    "nanospa": "NanoSPA", "tombo": "Tombo", "tombo_com": "Tombo_com",
    "xpore": "xPore", "xpor": "xPore", "yanocomp": "Yanocomp",
    "mafia": "mAFiA", "tandemmod": "TandemMod", "nanodoc2": "nanodoc2",
    "nanorms": "nanoRMS", "m1a-prediction": "m1a-prediction",
    "penguin": "penguin", "dorado": "Dorado",
    # upstream
    "guppy": "guppy", "guppy_basecaller": "guppy",
    "minimap2": "minimap2", "samtools": "samtools",
    "nanopolish": "nanopolish", "f5c": "f5c", "modkit": "modkit",
    "slow5tools": "slow5tools", "slow5tools-f2s": "slow5tools",
    "gffread": "gffread", "bedtools": "bedtools", "seqkit": "seqkit",
    "r2dtool": "R2Dtool", "guitar": "Guitar",
}

#: canonical -> conda env used by the benchmark scripts (filled by the collectors)
DEFAULT_ENVS: dict[str, str] = {
    "Nanom6A": "nanom6A", "DENA": "dena", "Tombo": "tombo", "Tombo_com": "tombo",
    "EpiNano_Error": "epinano", "EpiNano_SVM": "epinano",
    "NanoSPA": "NanoSPA", "NanoNm": "nanom6A", "NanoPsu": "nanopsu",
    "CHEUI_m6A": "cheui", "CHEUI_m5C": "cheui", "CHEUI-diff": "cheui",
    "m1a-prediction": "m1a-prediction", "TandemMod": "TandemMod",
    "nanodoc2": "nanodoc2", "nanoRMS": "nanoRMS", "penguin": "penguin",
    "Nanocompore": "nanocompore", "f5c": "f5c_env",
    "xPore": "xPore", "Yanocomp": "nanocompore", "m6Anet": "m6anet",
    "MINES": "tombo", "DRUMMER": "tombo", "ELIGOS2_diff": "eligos2",
    "ELIGOS2_solo": "eligos2",
}


def canonical_known(token: str) -> str | None:
    """Like :func:`canonical_tool` but returns ``None`` instead of guessing.

    Used by the collectors so that an unrecognised script is reported as
    ``unclassified`` (and kept in a side table) rather than silently becoming a
    fake tool name.
    """
    t = str(token).strip()
    if not t:
        return None
    key = t.lower()
    if key in ALIASES:
        return ALIASES[key]
    stem = Path(t).stem.lower()
    if stem in ALIASES:
        return ALIASES[stem]
    base = stem.split("-")[0]
    if base in ALIASES:
        return ALIASES[base]
    for name in TOOL_META:
        if key == name.lower() or stem == name.lower():
            return name
    # strip a trailing version: Epinano1.2.0 -> epinano, m6anet-v-2.1.0 -> m6anet,
    # nanom6A_2022_12_22 -> nanom6A
    stripped = re.sub(r"[-_]?\d{4}[-_]\d{2}[-_]\d{2}$", "", stem)
    stripped = re.sub(r"[-_]?v?\d+(\.\d+)+$", "", stripped).rstrip("-_v")
    if stripped:
        if stripped in ALIASES:
            return ALIASES[stripped]
        for name in sorted(TOOL_META, key=len, reverse=True):
            if stripped.startswith(name.lower()):
                return name
    # longest prefix first, e.g. Nanocompore_mESCs_f5c -> Nanocompore
    for name in sorted(TOOL_META, key=len, reverse=True):
        n = name.lower()
        if stem.startswith(n) or key.startswith(n):
            return name
    return None


def canonical_tool(token: str) -> str:
    """Normalise a raw token (script stem / source dir / executable) to a
    canonical tool name.  Unknown tokens are returned unchanged so that they
    still show up in the inventory instead of silently disappearing."""
    return canonical_known(token) or str(token)


# --------------------------------------------------------------------------- #
# Row contract
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class RawFact:
    """One observed fact about one tool, always with its provenance."""

    tool_canonical: str
    dataset: str
    field: str
    value: str
    source: str
    line_no: int
    method: str

    def as_row(self) -> dict:
        return asdict(self)


def write_facts(facts: Iterable[RawFact], path: Path, columns: list[str] | None = None,
                generated_at: str | None = None) -> int:
    """Write facts to CSV (pandas keeps the encoding sane for Chinese paths)."""
    import pandas as pd

    rows = [f.as_row() for f in facts]
    cols = columns or ["tool_canonical", "dataset", "field", "value",
                       "source", "line_no", "method"]
    df = pd.DataFrame(rows, columns=cols)
    if generated_at:
        df["generated_at"] = generated_at
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return len(df)


# --------------------------------------------------------------------------- #
# Logging / subprocess
# --------------------------------------------------------------------------- #
def setup_logger(name: str) -> logging.Logger:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(f"tool_inventory.{name}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False
    fmt = logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s",
                            "%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    fh = logging.FileHandler(LOG_DIR / f"{name}.log", mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    return logger


def run_cmd(cmd: list[str], timeout: int = 60, logger: logging.Logger | None = None
            ) -> tuple[int, str, str]:
    """Run a read-only command.  Never raises: failures are logged and returned."""
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout, p.stderr
    except FileNotFoundError as exc:
        if logger:
            logger.warning("command not found: %s (%s)", cmd, exc)
        return 127, "", str(exc)
    except subprocess.TimeoutExpired:
        if logger:
            logger.warning("timeout (%ss): %s", timeout, " ".join(cmd))
        return 124, "", f"timeout after {timeout}s"
    except Exception as exc:  # pragma: no cover - defensive
        if logger:
            logger.warning("command failed: %s (%s)", cmd, exc)
        return 1, "", str(exc)


def read_text_head(path: Path, logger: logging.Logger | None = None) -> str:
    """Read a text file safely: skip binaries/huge files, cap at ``HEAD_BYTES``."""
    try:
        if path.suffix.lower() in SKIP_SUFFIXES:
            return ""
        size = path.stat().st_size
        if size > MAX_TEXT_BYTES:
            with path.open("r", encoding="utf-8", errors="ignore") as fh:
                return fh.read(HEAD_BYTES)
        return path.read_text(encoding="utf-8", errors="ignore")
    except Exception as exc:
        if logger:
            logger.debug("cannot read %s: %s", path, exc)
        return ""


def iter_files(root: Path, suffixes: set[str]) -> Iterable[Path]:
    """Yield files under ``root`` with one of ``suffixes`` (case-insensitive)."""
    if not root.exists():
        return []
    return (p for p in root.rglob("*")
            if p.is_file() and p.suffix.lower() in suffixes)


def dataset_from_path(path: Path) -> str:
    """Infer the dataset / species a script or artefact belongs to."""
    parts = {p.lower() for p in path.parts}
    for key, label in (("arabidopsis", "Arabidopsis"), ("mouse", "Mouse"),
                       ("hela", "HeLa"), ("e. coli", "E.coli"), ("e.coli", "E.coli"),
                       ("curlcake", "Curlcake"), ("rna004", "RNA004"),
                       ("mESCs", "Mouse"), ("mES", "Mouse")):
        if key in parts:
            return label
    if "result_RNA004" in path.parts:
        return "RNA004"
    return "all"


def stamp() -> str:
    import time
    return time.strftime("%Y-%m-%d %H:%M:%S")
