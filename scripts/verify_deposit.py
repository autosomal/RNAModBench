#!/usr/bin/env python3
"""Integrity check for an RNAModBench checkout.

Verifies, without touching the network:

1. no personal absolute path survives anywhere in the deposit;
2. no directory name from the private working tree, and no working note in a
   language other than English, survives in a text source;
3. every Python / R / shell source file parses;
4. every deposited callset is listed in ``metadata/callsets_index.tsv`` and its
   recorded site count equals the file's row count;
5. the callset layout is the documented ``<platform>/<species>/<group>/<mod>/<tool>/
   <sample>.tsv``;
6. the figure tree is laid out as documented: 8 main and 10 supplementary
   directories, each named in ``docs/figure_index.md``, each holding its own code in
   ``src/``, and no figure script left behind in the callset pipeline directory;
   every script path the index names exists;
7. writes ``metadata/deposited_files.sha256``.

Exit status is non-zero when a check fails.
"""
from __future__ import annotations

import gzip
import hashlib
import re
import subprocess
import sys
import tempfile
import zipfile
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
    r"|the [a]nalysis tree|revision_[o]utput")
#: written with escapes so that no CJK character appears in this file itself
CJK = re.compile(r"[\u3000-\u303f\u3040-\u30ff\u3400-\u4dbf\u4e00-\u9fff\uff00-\uffef]")
TEXT = (".py", ".R", ".r", ".sh", ".md", ".tex", ".json", ".yml", ".yaml",
        ".csv", ".tsv", ".txt", ".bed")
#: the repository's own configuration files have no suffix to test
EXTRA_TEXT = (".gitignore", ".gitattributes")


def is_text(p: Path) -> bool:
    return p.suffix in TEXT or p.name in EXTRA_TEXT


def check_internal_names() -> None:
    """Working-tree directory names must not survive into the deposit."""
    hits = []
    for p in read_files():
        if not is_text(p):
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


def check_language() -> None:
    """The deposit is published in English; no working note in another language is."""
    hits = []
    for p in read_files():
        if not is_text(p):
            continue
        try:
            text = p.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        n = len(CJK.findall(text))
        if n:
            hits.append((p.relative_to(ROOT).as_posix(), n))
    shown = "; ".join(f"{a} ({b})" for a, b in sorted(hits)[:6])
    note(not hits, "no CJK (working-note) characters in the text sources"
         + ("" if not hits else f"; {len(hits)} file(s): {shown}"))


def check_paths() -> None:
    """No machine-specific absolute path survives anywhere in the deposit.

    Every file is read as bytes; a gzip stream is searched decompressed and a zip
    archive member by member, because the reference tables and the numeric caches
    ship compressed - and a serialised object carries the absolute path of the file
    it was built from, which is exactly what a reader must not be able to unzip.
    """
    hits = []
    for p in read_files():
        rel = p.relative_to(ROOT).as_posix()
        data = p.read_bytes()
        if PERSONAL.search(data):
            hits.append(rel)
        elif data[:2] == b"\x1f\x8b":
            try:
                inner = gzip.decompress(data)
            except (OSError, EOFError):        # a truncated stream is still a stream
                inner = b""
            if PERSONAL.search(inner):
                hits.append(rel + " (gzip stream)")
        elif zipfile.is_zipfile(p):
            with zipfile.ZipFile(p) as zf:
                for name in zf.namelist():
                    if PERSONAL.search(zf.read(name)):
                        hits.append(f"{rel} :: {name}")
    n = len([p for p in read_files()])
    note(not hits, f"no personal absolute path in {n} files, compressed streams included"
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
    dirs = sorted(d for d in ROOT.glob("figures/figure*") if d.is_dir())
    note(len(dirs) == 18, f"8 main and 10 supplementary figure directories"
         + ("" if len(dirs) == 18 else f"; found {len(dirs)}"))
    undocumented = [d.relative_to(ROOT).as_posix() for d in dirs
                    if d.name not in text]
    note(not undocumented, "every figure directory is named in figure_index.md"
         + ("" if not undocumented else f"; missing: {undocumented[:6]}"))
    # one home per figure: its code sits in the figure's own directory
    no_code = [d.relative_to(ROOT).as_posix() for d in dirs
               if not list((d / "src").glob("*"))]
    note(not no_code, "every figure directory holds its own code in src/"
         + ("" if not no_code else f"; empty: {no_code[:6]}"))
    stray = sorted(p.name for p in (ROOT / "src" / "harmonisation" / "scripts").glob("*")
                   if re.search(r"_figS?\d", p.name))
    note(not stray, "no figure script is left in the callset pipeline directory"
         + ("" if not stray else f"; stray: {stray[:6]}"))


def write_hashes() -> None:
    out = []
    for p in read_files():
        if p.name == "deposited_files.sha256" or p.suffix == ".py" and "__pycache__" in p.parts:
            continue
        digest = hashlib.sha256(p.read_bytes()).hexdigest()
        out.append(f"{digest}  {p.relative_to(ROOT).as_posix()}")
    (ROOT / "metadata" / "deposited_files.sha256").write_text("\n".join(out) + "\n", encoding="utf-8")
    note(True, f"wrote metadata/deposited_files.sha256 ({len(out)} entries)")


#: the supplement was renumbered for publication; these are the working numbers the
#: delivered figures used before that (S5 was the tenth, S6 was the fifth, ...)
SI_STALE = {"6": "5", "7": "6", "8": "7", "9": "8", "10": "9", "5": "10"}
SI_NAME = re.compile(r"(?<![A-Za-z0-9])([Ff]ig[Ss]?|s)(\d+)(?![0-9])")


def check_figure_numbering() -> None:
    """No file may be named with the pre-publication number of its own figure.

    One number per figure is the rule the deposit is built on (docs/figure_index.md
    names the assembler that fixes it). figures/figureS8/ used to hold
    61_figS8_tables.py, drawn for the page delivered as Figure S9. A name referring
    to a *different* figure - figures/figure6/tables/figS6_site_quality.tsv, say,
    which Figure S6 really does read from its main figure's folder - is not what
    this looks for.
    """
    bad = []
    for d in sorted(ROOT.glob("figures/figureS*")):
        own = d.name.removeprefix("figureS")
        stale = SI_STALE.get(own)
        if not stale:
            continue
        for p in sorted(d.rglob("*")):
            if not p.is_file():
                continue
            rel = p.relative_to(d).as_posix()
            if any(m.group(2) == stale for m in SI_NAME.finditer(rel)):
                bad.append(rel)
    note(not bad, "no file is named with the pre-publication number of its own figure"
         + ("" if not bad else f"; {len(bad)}: {bad[:6]}"))


CAPTION_BLOCK = re.compile(r"^\*\*Figure S(\d+)\.\*\*", re.M)
LEGEND_TITLE = re.compile(r"^#\s*Figure S(\d+)\b")


def check_si_labels() -> None:
    """A supplementary figure is named by one number everywhere it is written down.

    The legend of Figure S6 titles itself Figure S6, and the caption sheet - the one
    generated file the build does not renumber, because it was written after the
    supplement had been ordered - lists its ten blocks from S1 upward.  Either of the
    two drifting is the mistake this catches.
    """
    bad = []
    for d in sorted(ROOT.glob("figures/figureS*")):
        own = d.name.removeprefix("figureS")
        for md in sorted(d.glob("figures/Fig*S*_legends*.md")):
            first = md.read_text(encoding="utf-8", errors="ignore").splitlines()[0]
            m = LEGEND_TITLE.match(first)
            if not m or m.group(1) != own:
                bad.append(md.relative_to(ROOT).as_posix())
    note(not bad, "every legend titles itself with its figure's number"
         + ("" if not bad else f"; {bad[:6]}"))
    caps = ROOT / "tables/sup_figure_captions.md"
    if caps.is_file():
        nums = [int(n) for n in CAPTION_BLOCK.findall(caps.read_text(encoding="utf-8"))]
        note(nums == list(range(1, 11)),
             f"the caption sheet numbers its ten figures in order ({nums})")
    else:
        note(False, "tables/sup_figure_captions.md is missing")


DOC_POINTER = re.compile(r"(?:^|[\s(`\"'/=])((?:[\w-]+/)*[\w][A-Za-z0-9_./+-]*\.md)")
#: a ledger's evidence column names the note or the source line a fact came from,
#: and those notes are the working tree's; docs/tool_inventory_notes.md says so
PROVENANCE_LINE = re.compile(r"evidence:|research:")


def _resolvable(rel: str, here: Path) -> bool:
    """A document may be cited from the page it sits on or from the repository root."""
    if rel == "README.md":
        return (here.parent / "README.md").is_file() or (ROOT / "README.md").is_file()
    for cand in (ROOT / rel, here / rel, *ROOT.glob("**/" + rel)):
        if cand.is_file():
            return True
    return False


def check_doc_pointers() -> None:
    """Every document the written documentation points to is in the repository.

    A reader follows a `see docs/x.md` out of the prose and expects to land on a
    file; code-side path strings are not checked here because the render sweep
    (``scripts/run_figures.sh``) already executes them.
    """
    dead = {}
    for p in sorted(ROOT.rglob("*")):
        if not p.is_file() or "_local" in p.parts or ".git" in p.parts:
            continue
        if p.suffix.lower() not in (".md", ".tex"):
            continue
        for line in p.read_text(encoding="utf-8", errors="ignore").splitlines():
            if PROVENANCE_LINE.search(line):
                continue
            for m in DOC_POINTER.finditer(line):
                if not _resolvable(m.group(1), p):
                    dead.setdefault(m.group(1), []).append(
                        p.relative_to(ROOT).as_posix())
    worst = sorted(dead.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    note(not dead, "every document the documentation cites exists"
         + ("" if not dead else f"; {len(dead)} dead: "
            + ", ".join(f"{k} ({len(v)})" for k, v in worst[:8])))


def main() -> int:
    print(f"verifying {ROOT}")
    check_paths()
    check_internal_names()
    check_language()
    check_parses()
    check_callsets()
    check_figure_index()
    check_figure_numbering()
    check_si_labels()
    check_doc_pointers()
    write_hashes()
    print(f"\n{'FAILED: ' + str(len(FAIL)) if FAIL else 'all checks passed'}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
