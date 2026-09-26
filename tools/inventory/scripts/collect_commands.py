#!/usr/bin/env python3
"""Collect the *exact command lines* the benchmark actually executed (R1-10).

Reviewer 1 (R1-10) asks whether the exact command lines / configuration files
are deposited.  Every shell script under ``tool_scripts/detection``,
``tool_scripts/basecalling``, ``code/rerun``, ``code/f5c_mode`` and
``result_RNA004/scripts`` is therefore parsed line by line:

* ``conda activate <env>`` / ``source activate <env>`` -> conda environment
* logical commands (continuation lines joined) -> one row each, with the
  ``file:line`` provenance that can be quoted in the response letter

Output: ``raw/commands_raw.csv``
columns: tool_canonical, dataset, conda_env, command, source, line_no, method

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/tools/inventory/scripts/collect_commands.py
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ti_common as tc  # noqa: E402

NON_COMMAND = re.compile(
    r"^\s*(#|$|set\s|cd\s|echo\s|conda\s+activate|source\s+activate"
    r"|export\s+PATH|ulimit\s|mkdir\s|rm\s|cp\s|mv\s)", re.I)
ACTIVATE = re.compile(r"^\s*(?:conda|source)\s+activate\s+([^\s;#]+)", re.I)
ENV_EXPORT = re.compile(r"^\s*(?:export\s+)?(CUDA_VISIBLE_DEVICES|OMP_NUM_THREADS)\s*=")


@dataclass
class CommandRow:
    tool_canonical: str
    tool_raw: str
    dataset: str
    conda_env: str
    command: str
    source: str
    line_no: int
    method: str


def logical_lines(text: str) -> list[tuple[int, str]]:
    """Join backslash continuations; return (start_line, joined_command)."""
    out: list[tuple[int, str]] = []
    buf = ""
    start = 0
    for i, raw in enumerate(text.splitlines(), start=1):
        line = raw.rstrip()
        if not buf:
            start = i
        if line.endswith("\\"):
            buf += line[:-1].rstrip() + " "
            continue
        buf += line
        out.append((start, buf.strip()))
        buf = ""
    if buf:
        out.append((start, buf.strip()))
    return out


EXECS = ("guppy_basecaller", "dorado", "minimap2", "samtools", "nanopolish",
         "f5c", "modkit", "slow5tools", "tombo", "nanocompore", "bedtools",
         "gffread", "seqkit", "guitar", "R2Dtool")


def attribute(script: Path, cmd: str) -> tuple[str, str]:
    """Attribute a command to a tool.

    Order: script stem (``DENA.sh``) > source-code directory inside the command
    (``.../source_code/benchmark/NanoNm-1.0.0/...``) > a known executable >
    a ``*.py`` referenced by the command.  When nothing matches the row is
    reported as ``unclassified`` with the script stem kept in ``tool_raw``.
    """
    stem = script.stem
    if stem.lower() not in {"run", "all", "main", "step1", "step2", "step3",
                            "step4", "run_all", "single_tools"}:
        cand = tc.canonical_known(stem)
        if cand:
            return cand, stem
    m = re.search(r"source_code/benchmark([^/\s]+)", cmd)
    if m:
        cand = tc.canonical_known(m.group(1))
        if cand:
            return cand, stem
    for exe in EXECS:
        if re.search(rf"(^|[\s/;|&]){re.escape(exe)}(\s|$)", cmd):
            cand = tc.canonical_known(exe)
            if cand:
                return cand, stem
    m = re.search(r"([A-Za-z0-9_\-]+)\.(?:py|r|R|sh)\b", cmd)
    if m:
        cand = tc.canonical_known(m.group(1))
        if cand:
            return cand, stem
    return "unclassified", stem


def main() -> None:
    logger = tc.setup_logger("collect_commands")
    rows: list[CommandRow] = []
    n_files = 0

    for root in tc.SCRIPT_ROOTS:
        if not root.exists():
            logger.warning("script root missing: %s", root)
            continue
        for script in sorted(root.rglob("*.sh")):
            text = tc.read_text_head(script, logger)
            if not text:
                continue
            n_files += 1
            dataset = tc.dataset_from_path(script)
            env = ""
            for line_no, cmd in logical_lines(text):
                m_act = ACTIVATE.match(cmd)
                if m_act:
                    env = m_act.group(1).strip()
                    continue
                if not cmd or NON_COMMAND.match(cmd):
                    continue
                if ENV_EXPORT.match(cmd):
                    continue
                cmd = re.sub(r"\s+", " ", cmd).strip()
                if len(cmd) < 3:
                    continue
                tool, tool_raw = attribute(script, cmd)
                rows.append(CommandRow(
                    tool_canonical=tool,
                    tool_raw=tool_raw,
                    dataset=dataset,
                    conda_env=env,
                    command=cmd,
                    source=str(script),
                    line_no=line_no,
                    method="static_script",
                ))

    df = pd.DataFrame([r.__dict__ for r in rows],
                      columns=["tool_canonical", "tool_raw", "dataset", "conda_env",
                               "command", "source", "line_no", "method"])
    df["generated_at"] = tc.stamp()
    tc.RAW_DIR.mkdir(parents=True, exist_ok=True)
    out = tc.RAW_DIR / "commands_raw.csv"
    df.to_csv(out, index=False)

    logger.info("scanned %d shell scripts -> %d command rows", n_files, len(df))
    logger.info("tools seen: %s", sorted(df["tool_canonical"].unique()))
    logger.info("conda envs seen: %s",
                sorted({e for e in df["conda_env"].unique() if e}))
    logger.info("written -> %s", out)


if __name__ == "__main__":
    main()
