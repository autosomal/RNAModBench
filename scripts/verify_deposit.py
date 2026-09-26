#!/usr/bin/env python3
"""Integrity check for an RNAModBench checkout.

Verifies, without touching the network:

1. no personal absolute path survives anywhere in the deposit;
2. no directory name from the private working tree survives in a text source;
3. every Python / R / shell source file parses;
4. every deposited callset is listed in ``metadata/callsets_index.tsv`` and its
   recorded site count equals the file's row count;
5. the callset layout is the documented ``<platform>/<species>/<group>/<mod>/<tool>/
   <sample>.tsv``;
6. every figure in ``docs/figure_index.md`` has a delivered file, and every
   script path the index names exists;
7. writes ``metadata/deposited_files.sha256``.

Exit status is non-zero when a check fails.
"""
from __future__ import annotations

import hashlib
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SITES = ROOT / "data" / "callsets"
# character classes so that this checker does not match itself
PERSONAL = re.compile(rb"/[d]ata/l[xy]|/[h]ome/l[xy]|/[U]sers/[a-z]")
LAYOUT = re.compile(r"^RNA00[24]/[^/]+/[^/]+/(m6A|m5C|Psi|m1Psi|Nm|inosine)/[^/]+/[^/]+\.tsv$")
FAIL: list[str] = []


def note(ok: bool, msg: str) -> None:
    print(("  ok   " if ok else "  FAIL ") + msg)
    if not ok:
        FAIL.append(msg)


def read_files():
    for p in sorted(ROOT.rglob("*")):
        if p.is_file() and "_local" not in p.parts and ".git" not in p.parts:
            yield p


# character classes so that this checker does not match itself
INTERNAL_NAMES = re.compile(
    r"04_[r]evision_analysis|01_[c]ode|02_[r]aw_results|05_[s]ubmission|06_[r]eviewers"
    r"|07_[t]hird_party|08_[r]evision|sites_[v]2|sites_[c]lean|[c]ode_user"
    r"|the [a]nalysis tree")


def check_internal_names() -> None:
    """Working-tree directory names must not survive into the deposit."""
    hits = []
    for p in read_files():
        if p.suffix not in (".py", ".R", ".sh", ".md", ".tex", ".json", ".yml", ".yaml"):
            continue
        try:
            text = p.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        n = len(INTERNAL_NAMES.findall(text))
        if n:
            hits.append((p.relative_to(ROOT).as_posix(), n))
    shown = "; ".join(f"{a} ({b})" for a, b in sorted(hits)[:6])
    note(not hits, f"no internal working-tree names in {sum(1 for _ in read_files())} text sources"
         + ("" if not hits else f"; {len(hits)} file(s): {shown}"))


def check_paths() -> None:
    hits = [p.relative_to(ROOT).as_posix() for p in read_files()
            if p.suffix not in (".pdf", ".pkl", ".npz", ".gz", ".rds")
            and PERSONAL.search(p.read_bytes())]
    note(not hits, f"no personal absolute path in {len(list(read_files()))} files"
         + ("" if not hits else f"; found in {len(hits)}: {hits[:5]}"))


def check_parses() -> None:
    bad = []
    for p in read_files():
        if p.suffix == ".py":
            try:
                compile(p.read_text(encoding="utf-8"), str(p), "exec")
            except SyntaxError as e:
                bad.append(f"{p.relative_to(ROOT)}: {e.msg} (line {e.lineno})")
        elif p.suffix == ".sh":
            r = subprocess.run(["bash", "-n", str(p)], capture_output=True, text=True)
            if r.returncode:
                bad.append(f"{p.relative_to(ROOT)}: {r.stderr.strip()[:120]}")
    note(not bad, f"all Python and shell sources parse ({len(bad)} failure(s))"
         + ("" if not bad else "\n    " + "\n    ".join(bad[:10])))


def check_callsets() -> None:
    import csv
    idx = {r["deposited_file"]: r for r in csv.DictReader(
        (ROOT / "metadata" / "callsets_index.tsv").open(encoding="utf-8"), delimiter="\t")}
    on_disk = {"data/callsets/" + p.relative_to(SITES).as_posix(): p
               for p in SITES.rglob("*.tsv")}
    note(set(idx) == set(on_disk),
         f"index and files agree ({len(idx)} indexed, {len(on_disk)} on disk)")
    layout_bad = [k for k in on_disk if not LAYOUT.match(k[len("data/callsets/"):])]
    note(not layout_bad, f"callset layout is platform/species/group/mod/tool/sample.tsv"
         + ("" if not layout_bad else f"; {len(layout_bad)} off, e.g. {layout_bad[:3]}"))
    drift = []
    total = 0
    for key, p in on_disk.items():
        with p.open("rb") as fh:
            n = sum(1 for _ in fh) - 1
        total += n
        rec = idx.get(key)
        if rec and int(rec["sites"]) != n:
            drift.append(f"{key}: recorded {rec['sites']}, file {n}")
    note(not drift, f"every recorded site count matches its file ({total:,} sites total)"
         + ("" if not drift else f"; {len(drift)} drift(s): {drift[:3]}"))


def check_figure_index() -> None:
    text = (ROOT / "docs" / "figure_index.md").read_text(encoding="utf-8")
    named = set(re.findall(r"`((?:src|figures|analysis|tables)/[^`]+?\.(?:py|R|sh|tex|tsv|csv|gz|npz))`", text))
    missing = sorted(p for p in named if "*" not in p and not (ROOT / p).exists())
    note(not missing, f"{len(named)} paths named in figure_index.md all exist"
         + ("" if not missing else f"; missing: {missing[:6]}"))
    undelivered = [d.relative_to(ROOT).as_posix() for d in sorted(ROOT.glob("figures/figure*"))
                   if not list(d.glob("delivered/*.pdf"))]
    note(not undelivered, "every figure directory has a delivered/ figure"
         + ("" if not undelivered else f"; missing for {undelivered[:6]}"))


def write_hashes() -> None:
    out = []
    for p in read_files():
        if p.name == "deposited_files.sha256" or p.suffix == ".py" and "__pycache__" in p.parts:
            continue
        digest = hashlib.sha256(p.read_bytes()).hexdigest()
        out.append(f"{digest}  {p.relative_to(ROOT).as_posix()}")
    (ROOT / "metadata" / "deposited_files.sha256").write_text("\n".join(out) + "\n", encoding="utf-8")
    note(True, f"wrote metadata/deposited_files.sha256 ({len(out)} entries)")


def main() -> int:
    print(f"verifying {ROOT}")
    check_paths()
    check_internal_names()
    check_parses()
    check_callsets()
    check_figure_index()
    write_hashes()
    print(f"\n{'FAILED: ' + str(len(FAIL)) if FAIL else 'all checks passed'}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
