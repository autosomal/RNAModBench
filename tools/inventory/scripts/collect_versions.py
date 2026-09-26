#!/usr/bin/env python3
"""Probe real software versions for every tool in the benchmark (R1-5).

Three strategies, all read-only, each with its own evidence tag:

1. **conda**  -- ``conda list -n <env> --json`` for every environment that a
   benchmark script activated (env names come from ``raw/commands_raw.csv``).
2. **binary** -- ``<exe> --version`` for standalone tools (dorado, guppy,
   minimap2, samtools, nanopolish, f5c, modkit, slow5tools, ...).
3. **static** -- version encoded in an install path (``dorado-0.9.1-linux-x64``,
   ``$RNAMODBENCH_LOCAL/source_code/benchmark/NanoNm-1.0.0``), in ``yaml/*.yml`` pins or
   in a log header (``package_version:``).

conda is located automatically (``which`` -> ``conda info --base`` -> known
candidate roots).  If it cannot be located the collector still finishes and
tags every row ``conda not located`` instead of failing.

Output: ``raw/versions_raw.csv``
columns: entity_type, entity, package, version, method, evidence

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/tools/inventory/scripts/collect_versions.py
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
import json
import os
import re
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ti_common as tc  # noqa: E402

#: canonical tool -> candidate executables (name looked up on PATH first)
BIN_PROBES: dict[str, list[str]] = {
    # the benchmark used the 0.9.1 install (see IMPORTANT_INFO.md); the 2.0.0
    # install is listed afterwards so that both are on record
    "Dorado": [str(_XB / "tool/dorado-0.9.1-linux-x64/bin/dorado"),
               str(_XB / "tool/dorado-2.0.0-linux-x64/bin/dorado"),
               "dorado"],
    "guppy": ["guppy_basecaller", str(_XB / "tool/guppy/bin/guppy_basecaller")],
    "minimap2": ["minimap2"],
    "samtools": ["samtools"],
    "nanopolish": ["nanopolish"],
    "f5c": ["f5c"],
    "modkit": ["modkit"],
    "slow5tools": ["slow5tools"],
    "bedtools": ["bedtools"],
    "gffread": ["gffread"],
    "seqkit": ["seqkit"],
    "Tombo": ["tombo"],
    "Nanocompore": ["nanocompore"],
    "Nanom6A": ["nanom6A"],
}

#: keywords used to decide whether a conda package belongs to a tool
KEYWORD_TO_TOOL: dict[str, str] = {
    "nanom6a": "Nanom6A", "m6anet": "m6Anet", "cheui": "CHEUI_m6A",
    "dena": "DENA", "drummer": "DRUMMER", "eligos": "ELIGOS2_diff",
    "epinano": "EpiNano_Error", "differr": "differr", "mines": "MINES",
    "nanocompore": "Nanocompore", "nanomud": "NanoMUD", "nanonm": "NanoNm",
    "nanopsu": "NanoPsu", "nanospa": "NanoSPA", "tombo": "Tombo",
    "xpore": "xPore", "yanocomp": "Yanocomp", "mafia": "mAFiA",
    "tandemmod": "TandemMod", "nanodoc": "nanodoc2", "nanorms": "nanoRMS",
    "penguin": "penguin", "dorado": "Dorado", "minimap2": "minimap2",
    "ont-tombo": "Tombo", "eligos2": "ELIGOS2_diff", "drummer": "DRUMMER",
    "nanom6a": "Nanom6A", "epinano": "EpiNano_Error", "nano-doc": "nanodoc2",
    "samtools": "samtools", "nanopolish": "nanopolish", "f5c": "f5c",
    "modkit": "modkit", "slow5tools": "slow5tools", "ont-guppy": "guppy",
    "guppy": "guppy", "bedtools": "bedtools", "gffread": "gffread",
    "seqkit": "seqkit", "hdf5": "", "python": "", "numpy": "", "pandas": "",
    "scipy": "", "scikit-learn": "", "torch": "", "tensorflow": "",
    "xgboost": "", "h5py": "", "pysam": "", "ont-fast5-api": "",
    "ont-vbz-hdf-plugin": "", "biopython": "", "matplotlib": "",
    "seaborn": "", "statsmodels": "", "numba": "", "cython": "", "pip": "",
}

#: site-specific conda roots are supplied through $CONDA_EXTRA_ROOTS (colon-separated)
CONDA_CANDIDATES = [
    *(p for p in os.environ.get("CONDA_EXTRA_ROOTS", "").split(":") if p),
    "/opt/conda", "/opt/anaconda3", "/opt/miniconda3",
    "/usr/local/anaconda3", "/usr/local/miniconda3",
]

VERSION_RE = re.compile(r"(\d+\.\d+(?:\.\d+)?(?:[a-z0-9.\-+]*)?)")


def _base_has_envs(base: Path, envs: set[str]) -> bool:
    """True when ``base/envs`` contains at least one of the benchmark envs."""
    ed = base / "envs"
    if not ed.exists():
        return False
    have = {p.name.lower() for p in ed.iterdir() if p.is_dir()}
    return bool(have & {e.lower() for e in envs})


def locate_conda(logger, envs: set[str]) -> Path | None:
    """Find the conda installation that actually holds the benchmark envs.

    Preference: ``$CONDA_EXE`` (set when the collector itself runs inside a
    conda env) -> ``conda info --base`` -> ``which conda`` -> known roots.
    A candidate is only accepted when its ``envs/`` directory contains one of
    the environments the benchmark scripts activated.
    """
    cands: list[Path] = []
    env_exe = os.environ.get("CONDA_EXE")
    if env_exe:
        cands.append(Path(env_exe))
    exe = shutil.which("conda")
    if exe:
        rc, out, _ = tc.run_cmd(["conda", "info", "--base"], timeout=60,
                                logger=logger)
        if rc == 0 and out.strip():
            cands.append(Path(out.strip().splitlines()[-1]))
        cands.append(Path(exe).resolve().parent.parent)
    cands += [Path(c) for c in CONDA_CANDIDATES]

    seen: set[Path] = set()
    fallback: Path | None = None
    for c in cands:
        if not c.exists():
            continue
        base = c if c.is_dir() else c.parent.parent
        if base in seen:
            continue
        seen.add(base)
        if not ((base / "bin" / "conda").exists()
                or (base / "condabin" / "conda").exists()):
            continue
        fallback = fallback or base
        if _base_has_envs(base, envs):
            return base
    if fallback:
        logger.warning("no conda root holds the benchmark envs; using %s "
                       "(conda list will report 'env not found')", fallback)
        return fallback
    logger.warning("conda not located - falling back to static version inference")
    return None


def pick_version(text: str) -> str | None:
    """First plausible ``x.y`` / ``x.y.z`` token in a ``--version`` banner.

    Some tools print build numbers or CUDA versions first; requiring a major
    version below 1000 keeps those out.
    """
    for m in VERSION_RE.finditer(text):
        try:
            if float(m.group(1).split(".")[0]) < 1000:
                return m.group(1)
        except ValueError:
            continue
    return None


def conda_list(base: Path, env: str, logger) -> list[dict]:
    conda = base / "bin" / "conda"
    if not conda.exists():
        conda = base / "condabin" / "conda"
    rc, out, err = tc.run_cmd([str(conda), "list", "-n", env, "--json"],
                              timeout=180, logger=logger)
    if rc != 0:
        return []
    try:
        data = json.loads(out)
    except json.JSONDecodeError:
        # older conda prints a warning line before the JSON
        m = re.search(r"\[.*\]", out, re.S)
        if not m:
            logger.warning("conda list %s: unparsable output (%s)", env,
                           err.strip()[:120])
            return []
        try:
            data = json.loads(m.group(0))
        except json.JSONDecodeError:
            return []
    return [d for d in data if isinstance(d, dict)]


SOURCE_ROOT = Path(str(_XB / "source_code/benchmark"))
MODEL_SUFFIXES = {".h5", ".hdf5", ".pth", ".pt", ".ckpt", ".joblib", ".pkl",
                  ".model", ".bin"}


def version_from_source(d: Path, logger) -> str | None:
    """Read a version string from packaging files (setup.py/pyproject/__init__)."""
    for name in ("setup.py", "pyproject.toml", "__init__.py", "setup.cfg"):
        cands = list(d.glob(name))[:1] + list(d.glob(f"*/{name}"))[:3]
        for p in cands:
            txt = tc.read_text_head(p, logger)
            m = re.search(r"(?:__version__|version)\s*=\s*[\"']([\d][^\"']{0,20})[\"']",
                          txt)
            if m:
                return m.group(1)
    return None


def scan_source_code(logger) -> list[tc.RawFact]:
    """Version / commit / model files from the vendored tool source trees.

    Many benchmark tools were installed from source (CHEUI, DENA, DRUMMER,
    MINES, ELIGOS2, EpiNano, Nanom6A, ...), so ``conda list`` cannot see them.
    Here we take (a) the version encoded in the directory name, (b) the git
    commit when the checkout still carries ``.git``, and (c) the model files
    shipped inside the tree — which is exactly the ``model/checkpoint``
    information R1-5 asks for.
    """
    facts: list[tc.RawFact] = []
    if not SOURCE_ROOT.exists():
        logger.warning("source tree not found: %s", SOURCE_ROOT)
        return facts

    for d in sorted(SOURCE_ROOT.iterdir()):
        if not d.is_dir():
            continue
        tool = tc.canonical_known(d.name)
        if not tool:
            continue
        # (a) version in the directory name (x.y.z or a release date)
        vm = re.findall(r"(\d+\.\d+(?:\.\d+)?|\d{4}[-_]\d{2}[-_]\d{2})", d.name)
        if vm:
            facts.append(tc.RawFact(tool, "all", "software_or_version",
                                    vm[-1].replace("_", "-"), str(d), 0,
                                    "source_dir"))
        else:  # (a') version declared in the packaging files
            ver = version_from_source(d, logger)
            if ver:
                facts.append(tc.RawFact(tool, "all", "software_or_version", ver,
                                        str(d), 0, "source_file"))
        # (b) git commit
        if (d / ".git").exists():
            rc, out, _ = tc.run_cmd(["git", "-C", str(d), "log", "-1",
                                     "--format=%h|%ci"], timeout=30,
                                    logger=logger)
            if rc == 0 and out.strip():
                facts.append(tc.RawFact(tool, "all", "software_or_version",
                                        out.strip().replace("|", " "), str(d), 0,
                                        "source_git"))
        # (c) model / checkpoint files
        n = 0
        for p in d.rglob("*"):
            if not p.is_file() or p.suffix.lower() not in MODEL_SUFFIXES:
                continue
            if p.stat().st_size > 500 * 1024 * 1024:
                continue
            facts.append(tc.RawFact(tool, "all", "model_checkpoint",
                                    str(p.relative_to(SOURCE_ROOT)), str(p), 0,
                                    "model_file"))
            n += 1
            if n >= 10:
                break
        if n:
            logger.debug("%s: %d model files", d.name, n)
    logger.info("source tree: %d facts", len(facts))
    return facts


def main() -> None:
    logger = tc.setup_logger("collect_versions")
    rows: list[dict] = []

    # ---- environments actually used by the benchmark --------------------
    envs: set[str] = set(tc.DEFAULT_ENVS.values())
    cmd_csv = tc.RAW_DIR / "commands_raw.csv"
    if cmd_csv.exists():
        df = pd.read_csv(cmd_csv)
        envs |= {str(e).strip() for e in df["conda_env"].dropna().unique()
                 if str(e).strip()}
    logger.info("conda environments referenced by the benchmark: %s", sorted(envs))

    # ---- 1. conda -------------------------------------------------------
    base = locate_conda(logger, envs)
    if base:
        logger.info("conda base: %s", base)
        for env in sorted(envs):
            pkgs = conda_list(base, env, logger)
            if not pkgs:
                rows.append({"entity_type": "conda_env", "entity": env,
                             "package": "-", "version": tc.MISSING,
                             "method": "conda_list",
                             "evidence": "env not found or conda list failed"})
                continue
            for d in pkgs:
                name, ver = str(d.get("name", "")), str(d.get("version", ""))
                if not name:
                    continue
                rows.append({"entity_type": "conda_env", "entity": env,
                             "package": name, "version": ver,
                             "method": "conda_list",
                             "evidence": f"conda list -n {env}"})
                tool = KEYWORD_TO_TOOL.get(name.lower())
                if tool:
                    rows.append({"entity_type": "tool", "entity": tool,
                                 "package": name, "version": ver,
                                 "method": "conda_list",
                                 "evidence": f"conda list -n {env}"})
    else:
        for env in sorted(envs):
            rows.append({"entity_type": "conda_env", "entity": env,
                         "package": "-", "version": tc.MISSING,
                         "method": "conda_list",
                         "evidence": "conda not located"})

    # ---- 2. binaries ----------------------------------------------------
    for tool, cands in BIN_PROBES.items():
        found = False
        for cand in cands:
            exe = cand if cand.startswith("/") else shutil.which(cand)
            if not exe or not Path(exe).exists():
                continue
            rc, out, err = tc.run_cmd([str(exe), "--version"], timeout=60,
                                      logger=logger)
            # stdout first: some tools print CUDA/GPU banners on stderr
            first, ver = "", None
            for stream in ((out or "").strip(), (err or "").strip()):
                if not stream or rc != 0:
                    continue
                first = stream.splitlines()[0]
                ver = pick_version(first) or pick_version(stream)
                if ver:
                    break
            if ver:
                rows.append({"entity_type": "tool", "entity": tool,
                             "package": Path(exe).name,
                             "version": ver, "method": "bin_version",
                             "evidence": f"{exe} --version -> {first[:120]}"})
                found = True
                break
        if not found:
            rows.append({"entity_type": "tool", "entity": tool, "package": "-",
                         "version": tc.MISSING, "method": "bin_version",
                         "evidence": "executable not found on PATH / known dirs"})

    # ---- 3. static inference -------------------------------------------
    # conda environment definitions kept in the repository (yaml/*.yml)
    for p in sorted(tc.YAML_DIR.glob("*.y*ml")):
        text = tc.read_text_head(p, logger)
        in_dep = False
        for i, line in enumerate(text.splitlines(), start=1):
            if re.match(r"^\s*dependencies\s*:", line):
                in_dep = True
                continue
            if not in_dep:
                continue
            m = re.match(r"^\s*-\s*([A-Za-z0-9_.:=<>!~\-]+)\s*$", line)
            if not m:
                continue
            spec = m.group(1)
            name = spec.split("::")[-1]
            mv = re.search(r"[=<>!~]+([\d][^\s]*)", spec)
            ver = mv.group(1) if mv else "no pin"
            rows.append({"entity_type": "conda_env", "entity": p.stem,
                         "package": name, "version": ver,
                         "method": "static_yml", "evidence": f"{p}:{i}"})
            tool = KEYWORD_TO_TOOL.get(name.lower())
            if tool and ver != "no pin":
                rows.append({"entity_type": "tool", "entity": tool,
                             "package": name, "version": ver,
                             "method": "static_yml", "evidence": f"{p}:{i}"})

    tool_dir = Path(str(_XB / "tool"))
    if tool_dir.exists():
        for p in sorted(tool_dir.iterdir()):
            m = re.match(r"^([A-Za-z][A-Za-z0-9_]*)-\d", p.name)
            if not m:
                continue
            vm = VERSION_RE.search(p.name)
            if not vm:
                continue
            cand = tc.canonical_tool(m.group(1))
            if cand in tc.TOOL_META:
                rows.append({"entity_type": "tool", "entity": cand,
                             "package": p.name, "version": vm.group(1),
                             "method": "static_path",
                             "evidence": f"install directory {p}"})

    if cmd_csv.exists():
        df = pd.read_csv(cmd_csv)
        for _, r in df.iterrows():
            for m in re.finditer(
                    r"/([A-Za-z][A-Za-z0-9_.]*?)-(\d+\.\d+(?:\.\d+)?)(?:/|\s|$)",
                    str(r["command"])):
                cand = tc.canonical_tool(m.group(1))
                if cand not in tc.TOOL_META:
                    continue
                rows.append({"entity_type": "tool", "entity": cand,
                             "package": m.group(1), "version": m.group(2),
                             "method": "static_path",
                             "evidence": f"{r['source']}:{r['line_no']}"})

    out_df = pd.DataFrame(rows, columns=["entity_type", "entity", "package",
                                         "version", "method", "evidence"])
    out_df["generated_at"] = tc.stamp()
    tc.RAW_DIR.mkdir(parents=True, exist_ok=True)
    out = tc.RAW_DIR / "versions_raw.csv"
    out_df.to_csv(out, index=False)

    # ---- 4. vendored source trees (git commit / model files) -------------
    src = scan_source_code(logger)
    tc.write_facts(src, tc.RAW_DIR / "source_raw.csv", generated_at=tc.stamp())
    logger.info("source_raw.csv: %d facts", len(src))
    ok = int((out_df["version"] != tc.MISSING).sum())
    logger.info("written %d rows (%d with a real version) -> %s",
                len(out_df), ok, out)


if __name__ == "__main__":
    main()
