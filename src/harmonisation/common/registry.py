"""Sample / tool registry: resolve every (sample, tool) onto a source file.

The resolution rules were derived from a full inventory of ``result/`` and
``result_RNA004/`` (see the per-tool table in harmonisation/README.md):

* a tool's sample directories are matched to canonical samples through the
  alias lists in :mod:`config` (directory names vary per tool:
  ``HeLa_WT1`` / ``HeLa_WT_result1`` / ``HeLa_WTresult1``),
* for the Curlcake construct libraries, the tools used generic ``resultN``
  directories; the canonical library is recovered from the *file stem*
  (``Curlcake_IVT_result3_*`` -> ``Curlcake_IVT_rep3``),
* within a sample directory the tool-specific patterns from
  :data:`config.RNA002_TOOLS` are tried in order (``_remove_chr.txt`` beats the
  unfiltered ``.txt``),
* files that are header-only / zero bytes are reported as ``empty`` rather than
  silently dropped.
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
import gzip
import json
import re
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .config import (CALLSET_ROOT, CONVERTED_CALLSETS, LEGACY_OUTPUT, MODKIT_CODE,
                     PROJECT, RAW_OVER_CONVERTED, RESULT_RNA002, RESULT_RNA004,
                     RNA002_SOURCE_ROOTS, RNA002_TOOLS, RNA004_SOURCES, SAMPLES,
                     SAMPLES_BY_NAME, Sample, ToolSpec, independence_class,
                     resolve_coverage_bam, sequencing_unit)
from .io_utils import read_table
from .legacy_liftover import (EXPLICIT_DIR_SAMPLE, LEGACY_UNUSED_RAW,
                              PAIR_DIR_SAMPLE)

# --------------------------------------------------------------------------- #
# Curlcake construct libraries
# --------------------------------------------------------------------------- #
#: The legacy conversion notebooks used **two different numbering spaces** in a
#: Curlcake file stem, and conflating them hands one library's callset to another:
#:
#: * ``Curlcake_<class>_rep<N>``  -- N is the **library** index (1..3), the
#:   spelling DENA / MINES / m6Anet / ELIGOS2_solo used, e.g.
#:   ``Curlcake_IVT_rep2`` = ``RNA081120181part`` (copy map, ELIGOS2_solo rows).
#: * ``Curlcake_<class>_result<N>`` -- N is the **comparison** index within the
#:   class, i.e. the order in which that tool ran its Curlcake comparisons, NOT a
#:   library number.  The E. coli / Curlcake ``cp.sh`` lines and the
#:   ``result_tidy`` directory names agree on it: comparison IVT#1 =
#:   ``rep2_partial vs rep1``, IVT#2 = ``rep3 vs rep1``, m6A#1 = ``m6A_rep2 vs
#:   IVT_rep3``, m6A#2 = ``m6A_rep1 vs IVT_rep1``.
_STEM_REP_INDEX = re.compile(r"Curlcake_(IVT|m6A)_rep(?:licate)?(\d+)", re.IGNORECASE)
_STEM_RESULT_INDEX = re.compile(r"Curlcake_(IVT|m6A)_result(\d+)", re.IGNORECASE)

CURLCAKE_LIBRARY_INDEX: dict[tuple[str, int], str] = {
    ("IVT", 1): "Curlcake_IVT_rep1",
    ("IVT", 2): "Curlcake_IVT_rep2_partial",
    ("IVT", 3): "Curlcake_IVT_rep3",
    ("m6A", 1): "Curlcake_m6A_rep1",
    ("m6A", 2): "Curlcake_m6A_rep2",
}
CURLCAKE_COMPARISON_INDEX: dict[tuple[str, int], str] = {
    ("IVT", 1): "Curlcake_IVT_rep2_partial",
    ("IVT", 2): "Curlcake_IVT_rep3",
    ("m6A", 1): "Curlcake_m6A_rep2",
    ("m6A", 2): "Curlcake_m6A_rep1",
}


#: Tools whose Curlcake stems number the **comparison**, not the library.  The
#: remaining tools (DRUMMER, ELIGOS2_diff, xPore, and every single-library tool)
#: number the library, and the per-file row counts show it: e.g. xPore
#: ``Curlcake_IVT_result2`` (14 sites) is the ``rep2_partial`` library while
#: Nanocompore ``Curlcake_m6A_result2`` (1,734 sites) is ``m6A_rep1``.  Seeded
#: from the 2026-09-15 extraction ledger, which reproduces the published Fig.7/F1
#: site counts; re-derive from ``aggregation_copy_map.csv`` if a new comparison
#: directory appears.
CURLCAKE_COMPARISON_INDEX_TOOLS: frozenset[str] = frozenset({"Nanocompore"})


def curlcake_sample_from_stem(stem: str, tool: str = "") -> str | None:
    """Canonical library of a Curlcake file stem (see the two index spaces above)."""
    comparison = tool in CURLCAKE_COMPARISON_INDEX_TOOLS
    pairs = ((_STEM_RESULT_INDEX, CURLCAKE_COMPARISON_INDEX if comparison
              else CURLCAKE_LIBRARY_INDEX),
             (_STEM_REP_INDEX, CURLCAKE_LIBRARY_INDEX))
    for rex, table in pairs:
        m = rex.search(stem)
        if m:
            return table.get((m.group(1), int(m.group(2))))
    return None

#: ambiguous directory names: canonical -> alternative(s), flagged in output.
AMBIGUOUS_DIRS = {"E_IVT_neg": ["E_IVT_neg1", "E_IVT_neg2"]}

#: ``<base>_[rep|result]<n>[T]`` -- the replicate index must be preserved!
_VS_CTRL_RE = re.compile(r"_vs_.*$")
_REP_RE = re.compile(r"^(?P<base>.+?)_(?:rep|result)(?P<idx>\d+)(?P<tail>T?)$")
_TRAILING_SUFFIX_RE = re.compile(
    r"(_fast5|_transcripts_inference|_inference|_sort|_result)$")


def dir_base(name: str) -> tuple[str, int | None]:
    """``Arabidopsis_WT_result2`` -> ``('Arabidopsis_WT', 2)``.

    Trailing tool suffixes are removed *first* so that combined forms like
    ``<sample>_rep2_result`` keep their replicate index.  (Earlier versions
    stripped the index, silently collapsing rep1/rep2/rep3 onto one sample - the
    root cause of the "missing biological replicates".)
    """
    core = _TRAILING_SUFFIX_RE.sub("", name)
    #: ``<A>_vs_<B>`` comparison directories are deliberately NOT canonicalised
    #: from their name: which side owns the result differs per tool (DRUMMER and
    #: ELIGOS2 take their ``-t`` side, xPore the ``wt`` field of the yml), and the
    #: tidy-step directory names do not encode it consistently.  Their sample is
    #: recovered from the legacy file **stem** instead, which always carries
    #: ``Curlcake_<class>_result<N>`` (see ``_stem_source``).
    m = _REP_RE.match(core)
    if m:
        return m.group("base"), int(m.group("idx"))
    return core, None


def build_dir_index() -> dict[tuple[str, int | None], str]:
    """{(base, replicate_index): canonical sample} built from the alias lists."""
    index: dict[tuple[str, int | None], str] = {}
    for s in SAMPLES:
        for alias in (s.canonical, *s.aliases):
            base, idx = dir_base(alias)
            index.setdefault((base, idx), s.canonical)
    return index


DIR_INDEX = build_dir_index()
#: distinct canonical samples per base (used for the unambiguous fallback)
_BASE_SAMPLES: dict[str, set[str]] = {}
for (_b, _i), _c in DIR_INDEX.items():
    _BASE_SAMPLES.setdefault(_b, set()).add(_c)


def canonical_from_dir(dirname: str) -> tuple[str | None, bool]:
    """Map a result directory name to a canonical sample.  Returns (name, ambiguous)."""
    if dirname in AMBIGUOUS_DIRS:
        return AMBIGUOUS_DIRS[dirname][0], True
    #: ``<A>_vs_<B>`` comparison directories are never attributed by name: the
    #: side that owns the result differs per tool (DRUMMER/ELIGOS2 take ``-t``,
    #: xPore the yml ``wt``) and the tidy-step names do not encode it.  They are
    #: resolved from the legacy file stem by ``_stem_source`` instead; without
    #: this, a semantic sample alias is a prefix of its own comparison dir and
    #: silently hands one library's callset to another.
    if "_vs_" in dirname:
        return None, False
    base, idx = dir_base(dirname)
    if (base, idx) in DIR_INDEX:
        return DIR_INDEX[(base, idx)], False
    # a replicate-specific directory for a sample with a single replicate
    if idx is not None:
        cands = _BASE_SAMPLES.get(base, set())
        if len(cands) == 1:
            return next(iter(cands)), True
    # longest-prefix fallback (Curlcake_IVT_rep2_partial vs Curlcake_IVT_rep3)
    best: tuple[str, int] | None = None
    for (b, i), canonical in DIR_INDEX.items():
        if base.startswith(b) and (best is None or len(b) > best[1]):
            best = (canonical, len(b))
    if best is not None:
        return best[0], True
    return None, False


def sample_of_dir(tool: str, dirname: str) -> tuple[str | None, bool]:
    """Canonical sample of a tool result directory (+ ambiguous flag).

    Canonical copy shared by ``00/01`` (registry), ``01b`` (raw fill),
    ``13`` (completeness audit) and ``14`` (legacy coverage audit).  Resolution
    order:

    1. ``legacy_liftover.PAIR_DIR_SAMPLE`` -- declared ``<A>_vs_<B>`` comparison
       directories; the declaration is authoritative, hence **not** ambiguous;
    2. ``legacy_liftover.EXPLICIT_DIR_SAMPLE`` -- legacy group-level directory
       names (e.g. Nanocompore's ``E_IVT_neg``); ambiguous flag = True;
    3. the alias index (``canonical_from_dir``), which refuses ``_vs_`` names by
       design -- an undeclared pair directory therefore resolves to nothing and
       must be added to (1), never guessed.
    """
    declared = PAIR_DIR_SAMPLE.get(tool, {})
    if dirname in declared:
        return declared[dirname], False
    explicit = EXPLICIT_DIR_SAMPLE.get(tool, {})
    if dirname in explicit:
        return explicit[dirname], True
    return canonical_from_dir(dirname)


# --------------------------------------------------------------------------- #
# file resolution
# --------------------------------------------------------------------------- #
@dataclass
class Source:
    tool: str
    mod_type: str
    sample: str
    platform: str
    path: Path | None
    pattern: str
    kind: str            # final | raw | none
    status: str          # ok | empty | missing | pending
    size: int
    n_rows: int | None
    ambiguous: bool = False
    note: str = ""
    parser: str = "standard"
    is_raw: bool = False     # True => the source file is a raw format, not a callset
    #: JSON string of extra kwargs handed to the parser (e.g. ``{"mod_code": "a"}``
    #: for a Dorado pileup that mixes several modifications in one file).  Empty
    #: for every other source, so the column is purely additive.
    parser_arg: str = ""


def _file_status(path: Path) -> tuple[str, int, int | None]:
    size = path.stat().st_size
    if size == 0:
        return "empty", 0, 0
    n_rows: int | None = None
    if size < 5 * 1024 * 1024:  # only count tiny files (CephFS rule)
        try:
            with open(path, "rb") as fh:
                n_rows = max(sum(1 for _ in fh) - 1, 0)
        except OSError:
            n_rows = None
    kind = "empty" if n_rows == 0 else "ok"
    return kind, size, n_rows


def find_final_file(sample_dir: Path, patterns: tuple[str, ...]) -> tuple[Path | None, str]:
    """First matching file for the tool patterns (order = priority)."""
    for pat in patterns:
        hits = sorted(p for p in sample_dir.glob(pat) if p.is_file())
        if hits:
            return hits[0], pat
    return None, ""


# --------------------------------------------------------------------------- #
# generic ranked search (naming variants: _remove_chr / _liftover / _5mer / ...)
# --------------------------------------------------------------------------- #
_NON_CALLSET_SUFFIXES = (
    ".log", ".yml", ".yaml", ".hdf5", ".json", ".json.gzip", ".stats", ".sh",
    ".fa", ".fasta", ".fastq", ".fastq.gz", ".bam", ".bai", ".fai", ".dict",
    ".py", ".r", ".ipynb", ".out", ".tmp", ".partial",
)
_NON_CALLSET_NAMES = {"summary.txt", "sam_parse2.txt", "diffmod.table",
                      "read_level_m6A_sorted.txt"}


def derive_tokens(tool: ToolSpec) -> tuple[str, ...]:
    """Name tokens that identify a tool inside a result file name."""
    if tool.tokens:
        return tool.tokens
    toks = {tool.tool}
    base = tool.tool
    for suffix in ("_m6A", "_m5C", "_psi", "_m1psi", "_com"):
        if base.endswith(suffix):
            toks.add(base[: -len(suffix)])
    if base == "ELIGOS2_solo":
        toks.add("ELIGOS_solo")
    if base == "ELIGOS2_diff":
        toks.add("ELIGOS_diff")
    if base == "EpiNano_Error":
        toks.add("Epinano_Error")
    if base == "Tombo_com":
        toks.add("Tombo")
    if base == "NanoSPA_psU":
        toks.add("NanoSPA_Psu")
    if base == "yanocomp":
        toks.add("Yanocomp")
    return tuple(sorted(toks, key=len, reverse=True))


def _mod_marker_conflict(name_low: str, mod_type: str) -> bool:
    """Reject files that clearly belong to a *different* modification type."""
    if mod_type == "m6A" and ("_m5c" in name_low or "m5c_" in name_low):
        return True
    if mod_type == "m5C" and ("_m6a" in name_low or "m6a_" in name_low):
        return True
    if mod_type == "Psi" and "m1psi" in name_low:
        return True
    if mod_type == "m1Psi" and ("_psi_" in name_low or name_low.endswith("_psi.txt")):
        return True
    return False


def _excluded_name(name: str, tool: ToolSpec) -> bool:
    """True for known transcript-space / intermediate files of this tool."""
    low = name.lower()
    return any(low.endswith(suf.lower()) for suf in tool.exclude_suffixes)


def iter_callset_candidates(sample_dir: Path, tool: ToolSpec) -> Iterator[tuple[int, Path, str, str]]:
    """Yield ``(score, path, variant, parser)`` for every callset-like file.

    Shared by :func:`find_ranked_candidate` (best file of a *sample* directory)
    and :func:`find_stem_resolved_candidate` (files inside *generic*
    ``resultN`` directories, whose library is only recoverable from the stem).

    Intermediate/auxiliary files (logs, summaries, fastq/bam, hdf5, diffmod
    tables) are dropped.  ``*_unfiltered*`` files are penalised: for xPore the
    filtered ``*_xPore.txt`` is the callset the legacy copy map used.
    """
    tokens = derive_tokens(tool)
    stem = sample_dir.name
    for f in sample_dir.rglob("*"):
        if not f.is_file() or f.stat().st_size == 0:
            continue
        name = f.name
        low = name.lower()
        if low.endswith(_NON_CALLSET_SUFFIXES) or name in _NON_CALLSET_NAMES:
            continue
        if _excluded_name(name, tool):
            # the ``<sample>_<TOOL>.txt`` files of the *transcript* species are
            # r2d-input intermediates and stay excluded; the identically-shaped
            # Curlcake files (``Curlcake_IVT_repN_<tool>.txt`` /
            # ``Curlcake_m6A_resultN_<tool>.txt``) are the FINAL construct-space
            # callsets (no liftover stage exists there) -- never exclude those.
            if curlcake_sample_from_stem(name.rsplit(".", 1)[0], tool.tool) is None:
                continue
        if not low.endswith((".txt", ".tsv", ".bed")):
            continue
        if not any(t.lower() in low for t in tokens):
            continue
        if _mod_marker_conflict(low, tool.mod_type):
            continue
        if "cheui-diff" in low or "cheui_diff" in low:
            continue
        score = 0
        if "_remove_chr" in low:
            score += 100
        if "_liftover" in low:
            score += 50
        if "_5mer" in low:
            score += 20
        if "unfiltered" in low:
            score -= 10
        if low.endswith((".txt", ".tsv")):
            score += 5
        if stem.lower() in low:
            score += 3
        parser = "liftover" if "_liftover" in low else "standard"
        variant = ("remove_chr" if "_remove_chr" in low else
                   "liftover" if "_liftover" in low else
                   "5mer" if "_5mer" in low else "plain")
        yield score, f, variant, parser


def find_ranked_candidate(sample_dir: Path, tool: ToolSpec) -> tuple[Path | None, str, str]:
    """Best callset-like file in ``sample_dir`` matching the tool's tokens.

    Ranking (higher is better): ``_remove_chr`` > ``_liftover`` > ``_5mer`` >
    other ``.txt/.tsv/.bed``; files containing the sample stem get a bonus.

    Returns ``(path, variant, parser)``; the parser is ``liftover`` for
    ``*_liftover.txt`` files (13-column junction tables) and ``standard``
    otherwise.
    """
    best: tuple[int, Path, str, str] | None = None
    for cand in iter_callset_candidates(sample_dir, tool):
        if best is None or cand[0] > best[0]:
            best = cand
    if best is None:
        return None, "", "standard"
    return best[1], best[2], best[3]


def find_stem_resolved_candidate(tool_dir: Path, tool: ToolSpec,
                                 canonical: str) -> tuple[Path, str, str, str] | None:
    """Best callset-like file under ``tool_dir`` whose *stem* names ``canonical``.

    The Curlcake conversion notebooks wrote generic ``resultN`` directories
    (``result1``, ``result3``, ...) whose **name carries no sample**; the library
    is only recoverable from the file stem
    (``Curlcake_m6A_result1_xPore.txt`` -> ``Curlcake_m6A_rep1``).  Without this
    fallback those callsets were silently reported as ``raw_missing`` whenever a
    sample-named directory (holding only intermediate ``*.hdf5``) existed.

    Returns ``(path, variant, parser, dirname)`` or ``None``.
    """
    best: tuple[int, Path, str, str, str] | None = None
    for d in sorted(p for p in tool_dir.iterdir() if p.is_dir()):
        for score, f, variant, parser in iter_callset_candidates(d, tool):
            if curlcake_sample_from_stem(f.stem, tool.tool) != canonical:
                continue
            cand = (score, f, variant, parser, d.name)
            if best is None or cand[0] > best[0]:
                best = cand
    if best is None:
        return None
    return best[1], best[2], best[3], best[4]


def _stem_source(tool: ToolSpec, sample: Sample, tool_dir: Path,
                 ambiguous: bool) -> Source | None:
    """Curkcake fallback: recover the library from the file **stem**.

    The conversion notebooks wrote generic ``resultN`` directories (``result1``,
    ``result3``, ...) whose name carries no sample, so a Curlcake library is only
    recoverable from the stem (``Curlcake_m6A_result1_xPore.txt``).  Returns
    ``None`` for other species or when nothing matches.
    """
    if sample.species != "Curlcake":
        return None
    hit = find_stem_resolved_candidate(tool_dir, tool, sample.canonical)
    if hit is None:
        return None
    path, variant, parser, dirname = hit
    kind, size, n_rows = _file_status(path)
    return Source(tool.tool, tool.mod_type, sample.canonical, sample.platform,
                  path, path.name, "final", kind, size, n_rows,
                  ambiguous=ambiguous, parser=parser,
                  note=f"stem-resolved from {dirname} (variant={variant})")


def resolve_rna002_source(tool: ToolSpec, sample: Sample) -> Source:
    """Find the source file for one (tool, sample) pair.

    Walks :data:`~common.config.RNA002_SOURCE_ROOTS` in order -- the legacy
    converted callset (``$RNAMODBENCH_LOCAL/raw/converted_callsets/``) comes first,
    because it holds the recipe output the published numbers were computed from,
    and ``result/`` is searched afterwards.  Without the converted root almost
    every ``*_remove_chr.txt`` pattern resolves to ``missing`` (those files were
    promoted out of ``result/`` and exist nowhere else).

    Exception: pairs in :data:`~common.config.RAW_OVER_CONVERTED` skip the
    converted root, so the project's own rebuild under ``result/`` wins.  The six
    Arabidopsis EpiNano_Error chains are there because their legacy files were
    derived from a per-site table built against a human FASTA (2026-09-18).
    """
    for root in RNA002_SOURCE_ROOTS:
        if root == CONVERTED_CALLSETS and (sample.canonical, tool.tool) in RAW_OVER_CONVERTED:
            #: superseded legacy callset (see ``config.RAW_OVER_CONVERTED``): the
            #: archival file is known-defective, the rebuild under ``result/`` is
            #: the truth for this pair.
            continue
        src = _resolve_rna002_source_in(tool, sample, root / tool.result_subdir)
        if src is not None:
            return src
    return Source(tool.tool, tool.mod_type, sample.canonical, sample.platform,
                  None, "", "none", "missing", 0, None,
                  note="no source in any root", parser=tool.parser)


def _resolve_rna002_source_in(tool: ToolSpec, sample: Sample,
                              tool_dir: Path) -> Source | None:
    """One root of :func:`resolve_rna002_source`; ``None`` when nothing matched."""
    if not tool_dir.is_dir():
        return None

    subdirs = [d for d in sorted(tool_dir.iterdir()) if d.is_dir()]
    #: group-level comparison directories (no per-sample index): their sample is
    #: the tool's own "test side", see ``legacy_liftover.EXPLICIT_DIR_SAMPLE``.
    mapped: list[tuple[str, Path, bool]] = []
    for d in subdirs:
        canonical, amb = sample_of_dir(tool.tool, d.name)
        if canonical is not None and canonical == sample.canonical:
            mapped.append((d.name, d, amb))

    ambiguous = any(m[2] for m in mapped)
    if not mapped:
        return _stem_source(tool, sample, tool_dir, ambiguous)

    # multiple dirs may map to one sample (e.g. Curlcake_IVT_rep3 vs ...part):
    # prefer the exact-name directory, then the first alphabetical hit.
    mapped.sort(key=lambda t: (t[1].name != sample.canonical, t[1].name))
    for dirname, d, amb in mapped:
        path, pat = find_final_file(d, tool.patterns)
        if path is not None:
            kind, size, n_rows = _file_status(path)
            note = "" if dirname == sample.canonical else f"dir={dirname}"
            # `*_remove_chr.txt` is always the legacy standard-format callset;
            # other patterns may need the tool-specific raw parser.
            parser = "standard" if pat.endswith("_remove_chr.txt") else tool.parser
            return Source(tool.tool, tool.mod_type, sample.canonical, sample.platform,
                          path, pat, "final", kind, size, n_rows,
                          ambiguous=ambiguous, note=note, parser=parser)
        # ranked search over the tool's naming variants (_remove_chr/_liftover/...)
        alt, variant, alt_parser = find_ranked_candidate(d, tool)
        if alt is not None:
            kind, size, n_rows = _file_status(alt)
            note = f"variant={variant}"
            if dirname != sample.canonical:
                note += f"; dir={dirname}"
            return Source(tool.tool, tool.mod_type, sample.canonical, sample.platform,
                          alt, alt.name, "final", kind, size, n_rows,
                          ambiguous=ambiguous, note=note, parser=alt_parser)
        # raw fallbacks (e.g. f5c-mode Nanom6A ratio.0.5.tsv)
        for pat, parser in tool.raw_fallbacks:
            hits = sorted(p for p in d.glob(pat) if p.is_file())
            if hits:
                kind, size, n_rows = _file_status(hits[0])
                note = f"raw fallback; dir={dirname}" if dirname != sample.canonical else "raw fallback"
                return Source(tool.tool, tool.mod_type, sample.canonical, sample.platform,
                              hits[0], pat, "raw", kind, size, n_rows,
                              ambiguous=ambiguous, note=note, parser=parser, is_raw=True)
    # directories exist but no file matched: the Curlcake tools only keep
    # intermediates in the sample-named directories, while the converted callset
    # sits in the generic ``resultN`` directory named after the library.
    return _stem_source(tool, sample, tool_dir, ambiguous)


_DORADO_NAME_RE = re.compile(r"_(hac|sup)@(v[\d.]+)_(.*?)_pileup$")

#: sentinel used in ``config.RNA004_SOURCES``: do not trust the declared
#: modification, read it from the pileup's modkit ``name`` codes instead.
AUTO_MOD = "auto"


def parse_dorado_name(stem: str) -> tuple[str, str, str]:
    """``..._hac@v5.0.0_m6A_DRACH_pileup`` -> ('hac', 'v5.0.0', 'm6A_DRACH')."""
    m = _DORADO_NAME_RE.search(stem)
    if not m:
        return "NA", "", stem.split("_pileup")[0]
    basecall, version, model = m.group(1), m.group(2), m.group(3)
    model = model.split("#")[0]  # 'inosine_m6A#m6A' -> 'inosine_m6A'
    return basecall, version, model


def peek_modkit_codes(path: Path, limit: int = 200_000) -> list[str]:
    """Modkit codes present in a pileup, in file order (``[]`` if undeterminable).

    ``modkit pileup`` writes the modification into its ``name`` column
    (``a`` = m6A, ``m`` = m5C, ``17802`` = Psi, ``17596`` = inosine -- see
    ``parsers.MODKIT_CODE``).  A file may contain **several** codes: the Curlcake
    ``*_all_pileup.bed`` carries ``17802``/``a``/``m`` in one file, and the
    HeLa ``other_modification/`` files carry exactly one non-m6A code.  Only the
    first ``limit`` data rows are scanned (200 000 by default: far above the
    largest pileup here -- 52 k rows -- while still cheap, since only the ``name``
    field of each line is touched and the file is never loaded as a frame).
    """
    codes: list[str] = []
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = next((i for i, c in enumerate(header)
                        if c.strip().lstrip("#").lower() in ("name", "mod")), None)
            if idx is None:
                return []
            for i, line in enumerate(fh):
                if i >= limit:
                    break
                fields = line.rstrip("\n").split("\t")
                if len(fields) <= idx:
                    continue
                code = fields[idx].split("#")[0].strip()
                if code in MODKIT_CODE and code not in codes:
                    codes.append(code)
    except OSError:
        return []
    return codes


def resolve_rna004_sources(sample: Sample) -> list[Source]:
    """RNA004 sources are declared explicitly in config.RNA004_SOURCES.

    Dorado entries are *families*: every matching pileup file becomes its own
    callset, labelled ``Dorado_<basecall>_<model>`` (hac/sup x model x version).

    A Dorado entry declares ``mod_type="auto"``: the modification is read from the
    file (``peek_modkit_codes``), because the source directories do **not** always
    separate the modifications -- ``other_modification/`` holds files that are
    entirely 5mC / Psi / inosine, and the Curlcake ``dorado_model/`` pileups mix
    several codes in one file.  A file with several codes yields one callset per
    code; the m6A channel keeps the historical tool label so existing tables do
    not shift, the other channels get a ``_<mod_type>`` suffix.

    A declaration may carry a 5th element with extra parser kwargs (the Dorado
    pileup *call filter*, e.g. ``{"min_cov": 20, "min_pct": 0.0}``); they are
    merged into ``Source.parser_arg`` together with the modkit ``mod_code`` and
    handed to ``parse_callset`` unchanged.
    """
    specs = RNA004_SOURCES.get(sample.canonical, [])
    #: Both species keep their data under ``<species>/raw_calls/``; the Curlcake
    #: branch used to point at ``Curlcake/data`` which no longer exists after the
    #: 2026-09-15 reorganization, so all five Curlcake RNA004 sources silently
    #: resolved to ``missing`` while their pre-reorganization callset files were
    #: still on disk (and still read by ``09``).  Fixed 2026-09-17.
    root = (RESULT_RNA004 / "HeLa" if sample.species == "Human"
            else RESULT_RNA004 / "Curlcake" / "raw_calls")
    out: list[Source] = []
    for spec in specs:
        family, mod_type, pattern, parser = spec[:4]
        #: optional per-source parser kwargs (5th element); JSON-encoded into
        #: ``Source.parser_arg`` so no new column is needed.
        extra_kwargs: dict = spec[4] if len(spec) > 4 else {}

        def _arg(kwargs: dict) -> str:
            return json.dumps(kwargs) if kwargs else ""

        hits = [p for p in sorted(root.glob(pattern)) if p.is_file()]
        if not hits:
            out.append(Source(family, "m6A" if mod_type == AUTO_MOD else mod_type,
                              sample.canonical, sample.platform, None, pattern,
                              "none", "missing", 0, None, parser=parser,
                              note="declared auto, no source file"))
            continue
        is_raw = parser in ("dorado_pileup", "m6anet_csv", "nanospa_csv", "mafia_bed")
        for path in hits:
            #: (tool, mod_type, parser_arg) -- one Dorado pileup may hold several
            #: modifications, each of which becomes its own callset.
            variants: list[tuple[str, str, str]] = []
            if parser == "dorado_pileup":
                basecall, version, model = parse_dorado_name(path.stem)
                # the same model file name appears in several model-set
                # directories (m6A vs m6A_guitar vs other_modification), so the
                # directory must be part of the label.
                tag = ""
                for part in path.parts:
                    if part == "m6A_guitar":
                        tag += "_drachGuitar"
                    elif part == "other_modification":
                        tag += "_otherMod"
                tool = f"Dorado_{basecall}" + (f"@{version}" if version else "") + f"_{model}{tag}"
                codes = peek_modkit_codes(path)
                if not codes:
                    # header-only, unreadable, or a code we do not map: keep the
                    # declaration (and never invent a modification).
                    variants = [(tool, "m6A" if mod_type == AUTO_MOD else mod_type,
                                 _arg(extra_kwargs))]
                elif len(codes) == 1:
                    variants = [(tool, MODKIT_CODE[codes[0]],
                                 _arg(dict(extra_kwargs, mod_code=codes[0])))]
                else:
                    for code in codes:
                        mt = MODKIT_CODE[code]
                        label = tool if mt == "m6A" else f"{tool}_{mt}"
                        variants.append((label, mt,
                                         _arg(dict(extra_kwargs, mod_code=code))))
            else:
                variants = [(family, mod_type, _arg(extra_kwargs))]
            kind, size, n_rows = _file_status(path)
            for tool, mt, parser_arg in variants:
                out.append(Source(tool, mt, sample.canonical, sample.platform,
                                  path, pattern, "raw" if is_raw else "final",
                                  kind, size, n_rows, parser=parser, is_raw=is_raw,
                                  parser_arg=parser_arg,
                                  note=f"model={path.name}" if parser == "dorado_pileup" else ""))
    return out


# --------------------------------------------------------------------------- #
# legacy cross-check
# --------------------------------------------------------------------------- #
def load_legacy_copy_map() -> pd.DataFrame:
    
    
    
    p = ((_RB / "metadata/sample_metadata/scripts/mapping/aggregation_copy_map.csv"))
    if not p.exists():
        return pd.DataFrame(columns=["tool", "header_tool", "output_group", "src_parent_dir",
                                     "src_file", "src_exists"])
    # NOTE: the curated mapping files are CSV (comma separated), not TSV.
    return read_table(p, sep=",")


def infer_tool_from_filename(fname: str) -> str | None:
    """Recover the tool id from a callset file name (longest match wins)."""
    best: str | None = None
    for t in RNA002_TOOLS:
        for token in (t.tool, t.tool.replace("_m6A", ""), t.tool.replace("_m5C", "")):
            if token and token.lower() in fname.lower() and (best is None or len(token) > len(best)):
                best = t.tool
    return best


def legacy_expected_pairs(copy_map: pd.DataFrame) -> pd.DataFrame:
    """(canonical sample, tool) pairs that the legacy cp.sh expected to exist."""
    tmap = legacy_tool_map()
    known = {t.tool for t in RNA002_TOOLS}
    rows = []
    for _, r in copy_map.iterrows():
        tool = tmap.get(str(r["tool"]), str(r["tool"]))
        if tool not in known:
            # a few cp.sh lines put the output group in the tool column; recover
            # the tool from the source file name instead.
            tool = infer_tool_from_filename(str(r["src_file"])) or tool
        if tool not in known:
            continue
        # ``sample_of_dir`` applies the per-tool group-directory rules, so the
        # E.coli IVT comparison files resolve to the tool's tested side.
        sample, amb = sample_of_dir(tool, str(r["src_parent_dir"]))
        if sample is None:
            # Curkcake: recover from the file stem
            sample = curlcake_sample_from_stem(str(r["src_file"]), str(r["tool"]))
        if sample is None:
            continue
        rows.append({"sample": sample, "tool": tool,
                     "src_file": str(r["src_file"]),
                     "src_exists": str(r["src_exists"]),
                     "output_group": str(r["output_group"]),
                     "ambiguous": amb})
    return pd.DataFrame(rows)


def legacy_tool_map() -> dict[str, str]:
    """Map the tool label used in cp.sh onto RNA002_TOOLS ids."""
    mapping: dict[str, str] = {}
    for t in RNA002_TOOLS:
        mapping[t.tool] = t.tool
    mapping.update({
        "ELIGOS_diff": "ELIGOS2_diff", "ELIGOS_solo": "ELIGOS2_solo",
        "EpiNano_DiffErr": "EpiNano_Error", "Epinano_DiffErr": "EpiNano_Error",
        "Tombo": "Tombo_com", "Xpore": "xPore", "xPore": "xPore",
        "yanocomp": "yanocomp", "Yanocomp": "yanocomp",
        "NanoSPA": "NanoSPA_m6A",
    })
    return mapping


# --------------------------------------------------------------------------- #
# top level
# --------------------------------------------------------------------------- #
def build_registry() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return (sample_registry, sample_tool_registry, pending)."""
    sample_rows = []
    for s in SAMPLES:
        bam, bam_src = resolve_coverage_bam(s)
        sample_rows.append({
            "sample": s.canonical, "platform": s.platform, "species": s.species,
            "dataset_group": s.dataset_group, "condition_class": s.condition_class,
            "sample_type": s.sample_type, "study": s.study, "role": s.role,
            "replicate_tag": s.replicate_tag,
            "sequencing_unit": sequencing_unit(s),
            "independence_class": independence_class(s),
            "aliases": "|".join(s.aliases),
            "coverage_bam": str(bam) if bam else "",
            "coverage_bam_source": bam_src,
            "note": s.note,
        })
    sample_registry = pd.DataFrame(sample_rows)

    pairs: list[dict] = []
    for tool in RNA002_TOOLS:
        for s in SAMPLES:
            if s.platform != "RNA002":
                continue
            #: A tool may keep its callset-carrying directory only in the
            #: converted layer (``result/yanocomp`` was promoted out of
            #: ``result/`` before the rebuild): checking ``RESULT_RNA002`` alone
            #: silently skipped the whole tool, so its ~30 (sample, tool) pairs
            #: never reached ``sample_tool_registry.csv`` and showed up as
            #: "not discovered" in ``pending.csv`` although every callset exists.
            if not any(tool.result_subdir in [d.name for d in root.iterdir() if d.is_dir()]
                       for root in RNA002_SOURCE_ROOTS if root.is_dir()):
                continue
            src = resolve_rna002_source(tool, s)
            pairs.append({
                "platform": src.platform, "species": s.species,
                "dataset_group": s.dataset_group, "sample": s.canonical,
                "replicate_tag": s.replicate_tag, "condition_class": s.condition_class,
                "sequencing_unit": sequencing_unit(s),
                "independence_class": independence_class(s),
                "mod_type": src.mod_type, "tool": src.tool,
                "source_file": str(src.path) if src.path else "",
                "pattern": src.pattern, "kind": src.kind, "status": src.status,
                "parser": src.parser, "is_raw": int(src.is_raw),
                "size_bytes": src.size, "n_rows": src.n_rows if src.n_rows is not None else "",
                "ambiguous_dir": int(src.ambiguous), "note": src.note,
                "parser_arg": src.parser_arg,
            })
    for s in SAMPLES:
        if s.platform != "RNA004":
            continue
        for src in resolve_rna004_sources(s):
            pairs.append({
                "platform": src.platform, "species": s.species,
                "dataset_group": s.dataset_group, "sample": s.canonical,
                "replicate_tag": s.replicate_tag, "condition_class": s.condition_class,
                "sequencing_unit": sequencing_unit(s),
                "independence_class": independence_class(s),
                "mod_type": src.mod_type, "tool": src.tool,
                "source_file": str(src.path) if src.path else "",
                "pattern": src.pattern, "kind": src.kind, "status": src.status,
                "parser": src.parser, "is_raw": int(src.is_raw),
                "size_bytes": src.size, "n_rows": src.n_rows if src.n_rows is not None else "",
                "ambiguous_dir": int(src.ambiguous), "note": src.note,
                "parser_arg": src.parser_arg,
            })
    pair_registry = pd.DataFrame(pairs)

    # ---- cross-check with the legacy copy map -> pending / gaps -------------
    copy_map = load_legacy_copy_map()
    expected = legacy_expected_pairs(copy_map)
    tmap = legacy_tool_map()
    have = set(zip(pair_registry["sample"], pair_registry["tool"]))
    #: ``empty`` (converted file with 0 rows -> header-only callset) and
    #: ``ok_zero`` count as covered: the pair HAS a callset, it simply has no
    #: sites (``13`` labels both ``ok_zero``).  Requiring status == "ok" made the
    #: legacy gap list report "discovered but empty/unreadable" for genuine
    #: 0-site results.
    _covered = pair_registry["status"].isin(("ok", "ok_zero", "empty"))
    have_ok = set(zip(pair_registry.loc[_covered, "sample"],
                      pair_registry.loc[_covered, "tool"]))
    leg_rows = []
    if not expected.empty:
        seen: set[tuple[str, str]] = set()
        for _, r in expected.iterrows():
            tool = tmap.get(r["tool"], r["tool"])
            key = (r["sample"], tool)
            if key in seen:
                continue
            seen.add(key)
            leg_rows.append({
                "sample": r["sample"], "tool": tool,
                "legacy_src_exists": r["src_exists"],
                "legacy_output_group": r["output_group"],
                "legacy_ambiguous_dir": int(r["ambiguous"]),
                "discovered": key in have,
                "discovered_ok": key in have_ok,
            })
    legacy_check = pd.DataFrame(leg_rows)

    # pending = legacy expected but not discovered (or discovered empty)
    if not legacy_check.empty:
        pending = legacy_check[~legacy_check["discovered_ok"]].copy()
        pending["reason"] = [
            "not discovered" if not d else "discovered but empty/unreadable"
            for d in pending["discovered"]]
    else:
        pending = pd.DataFrame(columns=["sample", "tool", "reason"])

    return sample_registry, pair_registry, pending
