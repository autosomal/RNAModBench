"""Central configuration for the sites_v2 rebuild.

The rebuild replaces the legacy ``output/`` assembly chain (cp.sh copies +
ad-hoc notebooks) with a reproducible extraction from the raw tool results in
``result/`` (RNA002) and ``result_RNA004/`` (RNA004).

Locked conventions (do not re-derive, see also code/revision/common/config.py):

* GLORI bed files are BED (0-based, half-open) with ``end - start == 1``; the
  reference nucleotide coordinate is ``start``.
* Every callset site is stored as a 1-bp interval ``[pos_raw, pos_raw + 1)`` in
  the *same* coordinate space as the GLORI files.  A tool site ``p`` matches a
  GLORI site ``g`` at window ``w`` iff ``same chromosome and |p - g| <= w``;
  ``w = 0`` is an exact single-nucleotide match.
* Chromosome labels are normalised with
  :func:`code.revision.common.match.fix_chromosome` (lowercase ``chr*``) — the
  GLORI files use ``1``/``chr1``/``Chromosome`` inconsistently, so all matching
  goes through the normaliser.
* ``pos_raw`` is the coordinate exactly as reported by the tool (after
  documented 0-based conversion only).  Offsets are *never* silently applied;
  they are recorded in the per-tool offset audit and exposed as annotation
  columns (``pos_center``/``dist_center``/``dist_drach_a``/``offset_flag``).

Outputs live under ``$RNAMODBENCH_LOCAL/sites_v2``; nothing outside
``$RNAMODBENCH_ROOT`` is written.
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
from dataclasses import dataclass, field
from pathlib import Path

# --------------------------------------------------------------------------- #
# Roots
# --------------------------------------------------------------------------- #
PROJECT = _RB  #: this deposit
REF_ROOT = _XB / "reference"  #: Ensembl/TAIR/GRC annotations -- not redistributed
NANOPORE_DATA = _XB / "nanopore" / "data"  #: raw FASTQ/fast5 -- not redistributed

#: Deposited layer: the per-site callsets (``sites_clean/``) and the frozen
#: evaluation tables live in the repository itself.
SITES_ROOT = PROJECT / "data"
EVALUATION_ROOT = SITES_ROOT / "evaluation"
TABLE_DIR = EVALUATION_ROOT / "tables"
LOG_DIR = EVALUATION_ROOT / "logs"
MANIFEST_DIR = PROJECT / "metadata"

#: Rebuild layer: everything the extraction stages read or write is intermediate
#: and is NOT redistributed.  Point $RNAMODBENCH_LOCAL at a tree that holds them
#: (see ``_local/README.md``) to re-run stages 00->34 from scratch.
CALLSET_ROOT = _XB / "sites_v2" / "callsets"
#: callsets that were produced by a tool but are NOT part of the manuscript
#: analysis (e.g. non-m6A tools accidentally run on Arabidopsis/Mouse/E.coli).
CALLSET_ROOT_EXTENDED = _XB / "sites_v2" / "callsets_extended"
UNIVERSE_ROOT = _XB / "sites_v2" / "universe"

RESULT_RNA002 = _XB / "raw" / "result"
RESULT_RNA004 = _XB / "raw" / "result_RNA004"
#: The published assembly's ``output/`` tree, archived on 2026-09-16.  Steps 12/14
#: read it read-only to cross-check the rebuild.
LEGACY_OUTPUT = _XB / "archive" / "output_legacy_20260916"
#: Legacy converted-callset layer (``*_remove_chr.txt`` / ``*_5mer.txt`` /
#: ``*_output.bed``): the on-disk input the ``01_extract`` stage globs, hence
#: required to reproduce ``callsets/`` from scratch.
CONVERTED_CALLSETS = _XB / "raw" / "converted_callsets"

#: Roots searched by :func:`common.registry.resolve_rna002_source`, in priority
#: order.  The converted layer comes FIRST because it holds the recipe output the
#: published numbers were computed from; ``result/`` is searched afterwards, where
#: a tool's intermediate files can outrank a converted callset (Nanom6A's
#: ``ratio.0.5.tsv`` is an aggregate over libraries, the per-library
#: ``*_Nanom6A.txt`` is not).
RNA002_SOURCE_ROOTS: tuple[Path, ...] = (CONVERTED_CALLSETS, RESULT_RNA002)

#: (sample, tool) pairs whose converted-legacy callsets are **superseded** and must
#: therefore not shadow the project's own rebuild under ``result/``.
#: 2026-09-18: the six Arabidopsis ``EpiNano_Error`` chains (WT rep1-3 + fip37
#: rep1-3) were built on per-site tables generated with a HUMAN FASTA
#: (``Epinano_Variants.py -r``), so their ``*_remove_chr.txt`` selected
#: "positions where the human genome carries an A": the KD callsets carry no GLORI
#: enrichment (2.19 % vs 15-26 % for the other tools on the same sample) and the WT
#: callsets are a ~26 % subsample of the tool's real output.  The converted files
#: stay untouched (archival copy, read-only) -- the re-derived callsets come from
#: ``result/EpiNano_DiffErr/`` (see ``code/epinano_refix/README.md``).
RAW_OVER_CONVERTED: frozenset[tuple[str, str]] = frozenset(
    {(s, "EpiNano_Error") for s in (
        "Arabidopsis_WT_rep1", "Arabidopsis_WT_rep2", "Arabidopsis_WT_rep3",
        "Arabidopsis_fip37_rep1", "Arabidopsis_fip37_rep2", "Arabidopsis_fip37_rep3")}
)

# --------------------------------------------------------------------------- #
# Platform / modification vocabulary
# --------------------------------------------------------------------------- #
PLATFORM_RNA002 = "RNA002"
PLATFORM_RNA004 = "RNA004"

#: modification label -> reference base the modification sits on (upper case,
#: {U -> T} DNA alphabet as used in the FASTA files).  ``None`` = any base.
MOD_REF_BASE: dict[str, str | None] = {
    "m6A": "A", "m5C": "C", "Psi": "T", "m1Psi": "T",
    "inosine": "A", "Nm": None,
}
MOD_LABEL = {
    "m6A": "m6A", "m5C": "m5C", "Psi": "\u03a8", "m1Psi": "m1\u03a8",
    "inosine": "inosine", "Nm": "Nm",
}

#: display order used in reports / manifests.
MOD_TYPE_ORDER: list[str] = ["m6A", "m5C", "Psi", "m1Psi", "Nm", "inosine"]

#: modkit ``name`` code -> canonical modification id.  Lives here (not in
#: ``parsers``) so that ``registry`` can read a pileup's true modification without
#: importing the parser package; ``parsers`` re-exports it for compatibility.
#: ``a`` = m6A, ``m`` = 5mC, ``17802`` = pseudouridine, ``17596`` = inosine,
#: ``76792`` = m1Psi, ``19228`` = m1A, ``19229`` = Nm.
MODKIT_CODE: dict[str, str] = {
    "a": "m6A", "m": "m5C", "17802": "Psi", "17596": "inosine",
    "76792": "m1Psi", "19228": "m1A", "19229": "Nm",
}

#: alternative spellings seen in tool tables (Greek letters, prose names).
#: Anything funnelling through :func:`canonical_mod_type` is normalised to the
#: ASCII id; the Greek spelling used to create a second ``Psi``/``\u03a8``
#: directory for the same modification and made ``MOD_REF_BASE`` lookup miss.
MOD_TYPE_ALIASES: dict[str, str] = {
    "\u03a8": "Psi", "psi": "Psi", "psU": "Psi", "pseU": "Psi",
    "pseudouridine": "Psi", "m1\u03a8": "m1Psi", "m1psi": "m1Psi",
    "2'-O-methylation": "Nm", "Nm": "Nm",
}


def canonical_mod_type(label: str) -> str:
    """Collapse any spelling of a modification label onto the ASCII id."""
    txt = str(label).strip()
    if txt in MOD_REF_BASE:
        return txt
    return MOD_TYPE_ALIASES.get(txt, txt)


#: (platform, species, dataset_group) combinations whose **non-m6A** callsets are
#: part of the manuscript (Fig. 7 + R2-2: the non-m6A tools were only ever run on
#: HeLa and the Curlcake constructs; RNA004 Psi is part of Fig. S8A, which plots
#: NanoPsu / NanoSPA_psU next to the Dorado models).  Everything else is archived
#: under ``callsets_extended/`` by ``11_scope_split.py`` and is NOT filled per
#: replicate, because it belongs to no figure and has no reference.
IN_SCOPE_NONM6A: set[tuple[str, str, str]] = {
    ("RNA002", "Human", "HeLa_WT"),
    ("RNA002", "Human", "HeLa_IVT"),
    ("RNA002", "Curlcake", "Curlcake_IVT"),
    # RNA004 Psi (Fig. S8A): NanoPsu + NanoSPA_psU on HeLa WT/IVT and Curlcake
    ("RNA004", "Human", "RNA004_HeLa_WT"),
    ("RNA004", "Human", "RNA004_HeLa_IVT"),
    ("RNA004", "Curlcake", "RNA004_Curlcake_IVT"),
}

#: Species whose **non-m6A** callsets are physically DELETED (not archived) by
#: ``11_scope_split.py``.  The non-m6A tools were only ever meant for the HeLa
#: cell line and the Curlcake constructs (Fig. 7 + R2-2); any m5C/Psi/m1Psi/Nm
#: callset produced for these three species belongs to no figure, has no
#: reference, and the user asked for them to be removed outright (empty folders
#: included).  Every other non-m6A combination is kept **in ``callsets/``** --
#: the RNA002 Fig. 7 groups and, since Fig. S8A, the three RNA004 groups listed
#: in ``IN_SCOPE_NONM6A`` above.  Nothing is archived any more, so
#: ``callsets_extended/`` is expected to be empty.
DELETE_NONM6A_SPECIES: set[str] = {"Arabidopsis", "Mouse", "E.coli"}


def delete_nonm6a(species: str, mod_type: str) -> bool:
    """True when a callset must be physically removed (out-of-scope species)."""
    return canonical_mod_type(mod_type) != "m6A" and species in DELETE_NONM6A_SPECIES


def in_scope(platform: str, species: str, dataset_group: str, mod_type: str) -> bool:
    """True when a (sample, modification) combination belongs to the manuscript."""
    if canonical_mod_type(mod_type) == "m6A":
        return True
    return (platform, species, dataset_group) in IN_SCOPE_NONM6A


# --------------------------------------------------------------------------- #
# Manuscript TOOL scope (which tool was used for which sample / modification)
# --------------------------------------------------------------------------- #
#: The manuscript (``manuscript_benchmark/manuscript.tex``, Experimental Section
#: "RNA Modification Detection Tools") states that **15 dRNA-seq tools** were
#: included: 12 for m6A, 3 for pseudouridine (Psi) and/or
#: N1-methylpseudouridine (m1Psi), one of them also being the only m5C tool,
#: plus one dedicated Nm tool.  The m6A tools were run on every species; the
#: non-m6A tools (m5C / Psi / m1Psi / Nm) were **only ever run on the HeLa cell
#: line and the Curlcake constructs** (Fig. 7), never on Arabidopsis / Mouse /
#: E. coli.
#:
#: Everything else found in ``result/`` (differr, EpiNano_SVM, Tombo_com,
#: CHEUI-diff, mAFiA, CHEUI on Curlcake) is **not part of the manuscript**: those
#: callsets are physically removed by ``11_scope_split.py`` and are never
#: re-extracted by ``01``/``01b`` (see ``tool_in_scope``).
#:
#: Evidence (locked, do not re-derive):
#: * ``revision_output/tables/NA5_compatibility_matrix.md`` ->
#:   "Per-sample presence of callsets (``output/<sample>/``)" -- the assembled
#:   tree that produced every figure;
#: * the per-sample tool directories of ``output/`` itself;
#: * ``code/next_postprocessing/Curlcake/F1 Scores for Different Tools.ipynb``:
#:   the manuscript uses the 12-tool variant (Yanocomp 2086 / Nanocompore 1931 /
#:   ELIGOS2_diff 1634, matching the published 2,085 / 1,930 / 1,633); the
#:   14-tool variant that adds Tombo_com / EpiNano_SVM was discarded.
ARTICLE_M6A_TOOLS: frozenset[str] = frozenset({
    "CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
    "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
    "NanoSPA_m6A", "xPore", "yanocomp",
})

#: non-m6A tools of the manuscript (Fig. 7): 3 Psi tools, 1 m1Psi tool,
#: 1 m5C tool (CHEUI, shared with the m6A panel), 1 Nm tool.
ARTICLE_NONM6A_TOOLS: frozenset[str] = frozenset({
    "CHEUI_m5C", "NanoMUD_psi", "NanoMUD_m1psi", "NanoNm", "NanoPsu",
    "NanoSPA_psU",
})

#: RNA004 tool label prefix of the Dorado built-in model families
#: (``Dorado_<hac|sup>@<version>_<model>``); the exact model list is documented
#: in ``result_RNA004/documents/IMPORTANT_INFO.md`` and shown in Fig. S8A.
DORADO_TOOL_PREFIX = "Dorado_"

#: E. coli: DENA / MINES were **never run** (no ``result/<tool>/E_*`` directory at
#: all, and ``code_user/detection/E. coli/*/{DENA,mines}.sh`` produced nothing);
#: the IVT run has no ELIGOS2_diff.
#: Nanom6A **is** part of the E. coli panel (scripts exist for WT+IVT and the raw
#: results are in ``result/Nanom6A/E_*``): it was missing from the published
#: Figure 1B only because those four runs came out empty; the 2026-09 f5c fill +
#: the 2026-09-18 anchor salvage produced them, so it is in scope now

_EC_M6A_WT = ARTICLE_M6A_TOOLS - {"DENA", "MINES"}
_EC_M6A_IVT = _EC_M6A_WT - {"ELIGOS2_diff"}
#: Curlcake_IVT: the legacy pipeline ran neither DRUMMER nor ELIGOS2_diff
#: (``result/DRUMMER/result{1,2}`` / ``result/ELIGOS2_diff/result{1,2}`` are
#: declared ``legacy_unused`` in ``common/legacy_liftover.LEGACY_UNUSED_RAW``).
_CC_IVT_M6A = frozenset({"DENA", "ELIGOS2_solo", "EpiNano_Error", "m6Anet", "MINES",
                         "Nanocompore", "Nanom6A", "NanoSPA_m6A", "xPore", "yanocomp"})
_CC_M6A_M6A = _CC_IVT_M6A | {"DRUMMER", "ELIGOS2_diff"}
#: HeLa non-m6A panel (Fig. 7A-C); Curlcake_IVT has the same tools except that
#: CHEUI_m5C was never run on the constructs.
_HELA_NONM6A: dict[str, frozenset[str]] = {
    "m5C": frozenset({"CHEUI_m5C"}),
    "Psi": frozenset({"NanoMUD_psi", "NanoPsu", "NanoSPA_psU"}),
    "m1Psi": frozenset({"NanoMUD_m1psi"}),
    "Nm": frozenset({"NanoNm"}),
}
_CC_IVT_NONM6A: dict[str, frozenset[str]] = {
    "m5C": frozenset(),  # CHEUI_m5C is HeLa-only in the manuscript
    "Psi": _HELA_NONM6A["Psi"],
    "m1Psi": _HELA_NONM6A["m1Psi"],
    "Nm": _HELA_NONM6A["Nm"],
}

#: RNA004: Dorado built-in model families + the BAM-based tools that could be
#: re-run on the new chemistry (Fig. 8, Fig. S8).
_RNA004_M6A_HELA = frozenset({"m6Anet", "NanoSPA_m6A", "ELIGOS2_solo",
                              "ELIGOS2_diff", "DRUMMER"})
_RNA004_M6A_CURLCAKE = frozenset({"m6Anet", "NanoSPA_m6A"})
_RNA004_PSI = frozenset({"NanoPsu", "NanoSPA_psU"})
#: Dorado pileups exist for **every** modification the RNA004 models call, not just
#: m6A: the ``other_modification/`` split holds 5mC / Psi / inosine channels and the
#: Curlcake ``*_all`` / ``*_pseU_m6A`` / ``*_inosine_m6A`` pileups carry several
#: modkit codes in one file (2026-09-17).  Listing only the m6A keys made
#: ``tool_in_scope`` reject the correctly-mod-typed callsets and ``11_scope_split``
#: deleted them as "out-of-scope tool".
_DORADO_GROUPS = frozenset({
    ("RNA004", group, mod)
    for group in ("RNA004_HeLa_WT", "RNA004_HeLa_IVT", "RNA004_Curlcake_IVT")
    for mod in ("m6A", "m5C", "Psi", "inosine")
})

#: ``(platform, dataset_group, mod_type) -> tools the manuscript used``.  This is
#: an *allow list* (not a required list): combinations the legacy pipeline never
#: produced (e.g. DRUMMER on Curlcake_IVT) are simply absent here and are NOT
#: anomalies.
_FULL_M6A_GROUPS = ("Arabidopsis_WT", "Arabidopsis_KD", "Mouse_WT", "Mouse_KO",
                    "HeLa_WT", "HeLa_IVT")

ARTICLE_TOOL_SCOPE: dict[tuple[str, str, str], frozenset[str]] = {
    # ---------------------------- RNA002, m6A (Fig. 1-6) --------------------
    **{("RNA002", g, "m6A"): ARTICLE_M6A_TOOLS for g in _FULL_M6A_GROUPS},
    ("RNA002", "E.coli_WT", "m6A"): _EC_M6A_WT,
    ("RNA002", "E.coli_IVT", "m6A"): _EC_M6A_IVT,
    ("RNA002", "Curlcake_IVT", "m6A"): _CC_IVT_M6A,
    ("RNA002", "Curlcake_m6A", "m6A"): _CC_M6A_M6A,
    # ---------------------------- RNA002, non-m6A (Fig. 7) ------------------
    **{("RNA002", "HeLa_WT", m): t for m, t in _HELA_NONM6A.items()},
    **{("RNA002", "HeLa_IVT", m): t for m, t in _HELA_NONM6A.items()},
    **{("RNA002", "Curlcake_IVT", m): t for m, t in _CC_IVT_NONM6A.items()},
    # ---------------------------- RNA004 (Fig. 8, Fig. S8) ------------------
    ("RNA004", "RNA004_HeLa_WT", "m6A"): _RNA004_M6A_HELA,
    ("RNA004", "RNA004_HeLa_IVT", "m6A"): _RNA004_M6A_HELA,
    ("RNA004", "RNA004_Curlcake_IVT", "m6A"): _RNA004_M6A_CURLCAKE,
    ("RNA004", "RNA004_HeLa_WT", "Psi"): _RNA004_PSI,
    ("RNA004", "RNA004_HeLa_IVT", "Psi"): _RNA004_PSI,
    ("RNA004", "RNA004_Curlcake_IVT", "Psi"): _RNA004_PSI,
}


def article_tool_scope(platform: str, species: str, dataset_group: str,
                       mod_type: str) -> frozenset[str]:
    """Tools the manuscript used for one (platform, group, modification).

    Falls back to the manuscript-wide m6A / non-m6A set for combinations that
    are not enumerated explicitly (keeps an unexpected new group from being
    silently emptied).
    """
    mod = canonical_mod_type(mod_type)
    allowed = ARTICLE_TOOL_SCOPE.get((platform, dataset_group, mod))
    if allowed is not None:
        return allowed
    return ARTICLE_M6A_TOOLS if mod == "m6A" else ARTICLE_NONM6A_TOOLS


def tool_in_scope(platform: str, species: str, dataset_group: str,
                  mod_type: str, tool: str) -> bool:
    """True when the manuscript actually used ``tool`` on this sample.

    Combines the modification-level scope (:func:`in_scope`) with the tool-level
    allow list, so a tool that was only ever run on HeLa / Curlcake (or not at
    all: differr, EpiNano_SVM, Tombo_com, CHEUI-diff, mAFiA) can neither be
    re-extracted by ``01``/``01b`` nor survive in ``callsets/``.
    """
    if not in_scope(platform, species, dataset_group, mod_type):
        return False
    key = (platform, dataset_group, canonical_mod_type(mod_type))
    if str(tool).startswith(DORADO_TOOL_PREFIX) and key in _DORADO_GROUPS:
        return True
    return str(tool) in article_tool_scope(platform, species, dataset_group, mod_type)


# --------------------------------------------------------------------------- #
# Replicate independence (R3-2 / E6 and the editor's statistical-independence
# concern).  ``dataset_group`` groups samples by *condition*; it is NOT a
# statement that its members are replicates of one another.  Three different
# facts have to be recorded separately, because treating a group as "n
# replicates" when it is not is precisely pseudoreplication:
#
#   nested_subset / same_run  -> one sequencing unit counted twice
#   cross_study               -> independent experiments, but a different
#                                comparison (study effect, not replicate spread)
# --------------------------------------------------------------------------- #
#: Sample whose reads are not an independent sequencing unit, mapped onto the
#: sample that owns that unit.  ``Curlcake_IVT_rep2_partial`` is a depth-matched
#: subset of the SRR8767348 run (site-set containment in the full run = 1.00;
#: 49,538 of 646,251 eventalign reads, matched to ``Curlcake_IVT_rep1``'s 50,218);
#: ``E_IVT_neg1``/``E_IVT_neg2`` are the two halves of ONE IVT sample
#: (SRR27228854, user-confirmed 2026-09-16; local halves 882,018 / 506,384 reads)
#: -- SRP478171 contains exactly one IVT_neg MinION run.
NON_INDEPENDENT_SAMPLE_OF: dict[str, str] = {
    "Curlcake_IVT_rep2_partial": "Curlcake_IVT_rep3",
    "E_IVT_neg1": "E_IVT_neg2",
}

#: Classification of the folds above, exported as the ``independence_class``
#: column of the manifests / IVT-FPR table so that downstream consumers cannot
#: count "two halves" or "parent + subset" as n = 2:
#:   nested_subset  -- depth-matched subset of another sample's run
#:   same_run_split -- second half of one physical run
NON_INDEPENDENT_REASON: dict[str, str] = {
    "Curlcake_IVT_rep2_partial": "nested_subset",
    "E_IVT_neg1": "same_run_split",
}

#: Groups whose members come from different SRA studies.  Each member is an
#: independent dataset, but they are not biological replicates of one another,
#: so between-member agreement is reported as *cross-study concordance* and the
#: group is never assigned ``n_replicates > 1``.
CROSS_STUDY_GROUPS: frozenset[str] = frozenset({"Mouse_WT", "Mouse_KO"})


def sequencing_unit(sample: Sample) -> str:
    """Canonical name of the independent sequencing unit a sample belongs to."""
    return NON_INDEPENDENT_SAMPLE_OF.get(sample.canonical, sample.canonical)


def independent_units(samples: list[Sample]) -> dict[str, list[str]]:
    """{unit_key: [member canonicals]} -- the denominators for any n_replicates."""
    units: dict[str, list[str]] = {}
    for s in samples:
        units.setdefault(sequencing_unit(s), []).append(s.canonical)
    return {k: v for k, v in sorted(units.items()) if v}


def replicate_class(group: Sample | str) -> str:
    """How the members of a dataset group may be interpreted statistically."""
    name = group.dataset_group if isinstance(group, Sample) else group
    if name in CROSS_STUDY_GROUPS:
        return "cross_study"
    members = [s for s in SAMPLES if s.dataset_group == name]
    return "replicate" if len(independent_units(members)) > 1 else "single"


def independence_class(sample: Sample | str) -> str:
    """How one sample may be used in a replicate-statistics context.

    ``independent``     -- its own sequencing unit (one library / run set);
    ``nested_subset``   -- depth-matched subset of another sample's run;
    ``same_run_split``  -- one half of a physical run split in two;
    ``cross_study``     -- an independent dataset from a different study
                           (never collapsed into n > 1).
    """
    if isinstance(sample, Sample):
        if sample.dataset_group in CROSS_STUDY_GROUPS:
            return "cross_study"
        return NON_INDEPENDENT_REASON.get(sample.canonical, "independent")
    return "cross_study" if sample in CROSS_STUDY_GROUPS else "independent"

# --------------------------------------------------------------------------- #
# Analysis constants
# --------------------------------------------------------------------------- #
SEED = 20260914
BOOTSTRAP_B = 1000

#: matching windows (bp, +/-) for the distance-vs-accuracy sweep.
WINDOWS: list[int] = [0, 1, 2, 5, 10, 20, 50]
#: primary window used for headline numbers (matches code/revision).
PRIMARY_WINDOW = 2
#: coverage thresholds scanned when building the candidate universe.
C_MIN_SCAN: list[int] = [5, 10, 20]
#: default minimum coverage for a position to enter the universe.
C_MIN_DEFAULT = 10

#: DRACH / RRACH regexes on the DNA alphabet (U -> T), plus strand.
DRACH_REGEX = re.compile(r"[AGT][AG]AC[ACT]")
RRACH_REGEX = re.compile(r"[AG][AG]AC[ACT]")

#: maximum offset searched when locating the nearest expected base / DRACH A.
CENTER_SEARCH_MAX = 5

# --------------------------------------------------------------------------- #
# References
# --------------------------------------------------------------------------- #
GENOMES: dict[str, Path] = {
    "Arabidopsis": (_XB / "reference/arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"),
    "Mouse": (_XB / "reference/GRCm39/ensembl/Mus_musculus.GRCm39.dna.primary_assembly.fa"),
    # NOTE: the HeLa BAMs with '1'/'2'/... contigs were aligned against
    # GRCh38p14/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa, which no
    # longer exists on disk; GRCh38p13/GRCh38.primary_assembly.genome.fa is the
    # same primary assembly with 'chr' prefixes (normalised on access).
    "Human": (_XB / "reference/GRCh38p13/GRCh38.primary_assembly.genome.fa"),
    "E.coli": (_XB / "reference/K_12/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.toplevel.fa"),
    "Curlcake": (_XB / "reference/curlcakes/cc.fasta"),
}

TRANSCRIPTS: dict[str, Path] = {
    "Arabidopsis": (_XB / "reference/arabidopsis/Arabidopsis_thaliana.TAIR10.61.cdna.all.fa"),
    "Mouse": (_XB / "reference/GRCm39/ensembl/GRCm39.transcripts.fa"),
    #: Ensembl 112 (user rule 2026-09-18: GENCODE is banned project-wide)
    "Human": (_XB / "reference/GRCh38p14/ensembl112/GRCh38.transcripts.fa"),
    "E.coli": (_XB / "reference/K_12/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.61.cdna.all.fa"),
    "Curlcake": (_XB / "reference/curlcakes/cc.fasta"),
}

#: species -> precomputed exon-interval file (CSV: contig,source,type,start,end,
#: unk1,strand,unk2,meta,txid; 1-based start/end) used to build the strand-aware
#: candidate universe.  Missing ``.gtf.exon`` files fall back to parsing the GTF
#: directly (exon features only).
GTF_EXON: dict[str, Path] = {
    "Arabidopsis": (_XB / "reference/arabidopsis/Arabidopsis_thaliana.TAIR10.61.ensembl.gtf.exon"),
    "Mouse": (_XB / "reference/GRCm39/ensembl/Mus_musculus.GRCm39.114.gtf.exon"),
    #: Ensembl 112 (user rule 2026-09-18: GENCODE is banned project-wide)
    "Human": (_XB / "reference/GRCh38p14/ensembl112/Homo_sapiens.GRCh38.112.chr.gtf.exon"),
    "E.coli": (_XB / "reference/K_12/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.61.gtf.exon"),
    "Curlcake": (_XB / "reference/curlcakes/Curlcake.gtf.exon"),
}

#: tool -> (directory under ``RESULT_RNA002``, file name) of the tool's own
#: read-alignment BED.  Used by ``32_impute_strand`` to recover the transcript
#: strand for callsets that report none: Nanom6A writes ``'*'`` for **every**
#: row, and its ``extract.bed12`` carries the strand of the very reads the calls
#: came from -- on HeLa_WT1 it agrees with the exon annotation on 99.74 % of the
#: positions both can resolve and settles 70 % of the ones the exons cannot
#: (genes overlapping on both strands).
READ_STRAND_BED: dict[str, tuple[str, str]] = {
    "Nanom6A": ("Nanom6A", "extract.bed12"),
}

#: species -> GLORI reference bed (None = no GLORI available).
GLORI: dict[str, Path | None] = {
    "Arabidopsis": (_XB / "third_party/NGS/GLORI/Arabidopsis_GLORI.bed"),
    # Mouse MUST use the liftover version (41,961 sites, matched to GRCm39).
    "Mouse": (_XB / "third_party/NGS/GLORI/Mouse_GLORI_liftover.bed"),
    "Human": (_XB / "third_party/NGS/GLORI/Hela_GLORI.bed"),
    "E.coli": None,
    "Curlcake": None,
}

#: GLORI "all A" (coverage >= 15) candidate files, used only as a sanity
#: universe cross-check for Arabidopsis / Mouse.
GLORI_ALL_A: dict[str, list[Path]] = {
    "Arabidopsis": [
        (_XB / "third_party/NGS/GLORI/GSM7873784_Arabidopsis_allAs_cov15_rep1.bed"),
        (_XB / "third_party/NGS/GLORI/GSM7873785_Arabidopsis_allAs_cov15_rep2.bed"),
    ],
    "Mouse": [
        (_XB / "third_party/NGS/GLORI/GSM7873782_mESC_allAs_cov15_rep1.bed"),
        (_XB / "third_party/NGS/GLORI/GSM7873783_mESC_allAs_cov15_rep2.bed"),
    ],
}

#: Curlcake synthetic reference (4 constructs).
CURLCAKE_FASTA = (_XB / "reference/curlcakes/cc.fasta")
CURLCAKE_ALL_A = (_XB / "reference/curlcakes/curlcake_A.txt")

# --------------------------------------------------------------------------- #
# Samples
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Sample:
    """A canonical sequencing sample (one library / one sequencing run set)."""

    canonical: str
    platform: str
    species: str
    dataset_group: str
    condition_class: str  # WT / KD / KO / IVT / modified
    sample_type: str
    study: str
    role: str
    replicate_tag: str = ""
    aliases: tuple[str, ...] = ()  # directory-name aliases across result/ trees
    note: str = ""

    @property
    def is_negative_control(self) -> bool:
        return self.condition_class == "IVT"

    @property
    def is_partial_negative(self) -> bool:
        return self.condition_class in ("KD", "KO")


def _aliases(*names: str) -> tuple[str, ...]:
    return tuple(dict.fromkeys(n for n in names if n))


_AT = "Arabidopsis"
_MO = "Mouse"
_HE = "Human"
_EC = "E.coli"
_CC = "Curlcake"

SAMPLES: list[Sample] = [
    # ---------------- HeLa (human, RNA002) ----------------
    Sample("HeLa_WT1", PLATFORM_RNA002, _HE, "HeLa_WT", "WT", "cell line",
           "SRP393373", "WT reference", "rep1", _aliases("HeLa_WT1", "HeLa_WT_rep1", "HeLa_WT_result1", "HeLa_WTresult1")),
    Sample("HeLa_WT2", PLATFORM_RNA002, _HE, "HeLa_WT", "WT", "cell line",
           "SRP393373", "WT reference", "rep2", _aliases("HeLa_WT2", "HeLa_WT_rep2", "HeLa_WT_result2", "HeLa_WTresult2")),
    Sample("HeLa_WT3", PLATFORM_RNA002, _HE, "HeLa_WT", "WT", "cell line",
           "SRP393373", "WT reference", "rep3", _aliases("HeLa_WT3", "HeLa_WT_rep3", "HeLa_WT_result3", "HeLa_WTresult3")),
    Sample("HeLa_IVT_rep1", PLATFORM_RNA002, _HE, "HeLa_IVT", "IVT", "IVT control",
           "SRP428418", "negative control", "rep1",
           _aliases("HeLa_IVT_rep1", "HeLa_mRNA_IVT_rep1", "HeLa_IVT_rep1_fast5", "HeLa_IVT_result1", "HeLa_IVTresult1")),
    Sample("HeLa_IVT_rep2", PLATFORM_RNA002, _HE, "HeLa_IVT", "IVT", "IVT control",
           "SRP428418", "negative control", "rep2",
           _aliases("HeLa_IVT_rep2", "HeLa_mRNA_IVT_rep2", "HeLa_IVT_rep2_fast5", "HeLa_IVT_result2", "HeLa_IVTresult2")),
    Sample("HeLa_IVT_rep3", PLATFORM_RNA002, _HE, "HeLa_IVT", "IVT", "IVT control",
           "SRP428418", "negative control", "rep3",
           _aliases("HeLa_IVT_rep3", "HeLa_mRNA_IVT_rep3", "HeLa_IVT_rep3_fast5", "HeLa_IVT_result3", "HeLa_IVTresult3")),
    # ---------------- Arabidopsis (RNA002) ----------------
    Sample("Arabidopsis_WT_rep1", PLATFORM_RNA002, _AT, "Arabidopsis_WT", "WT", "whole seedling",
           "SRP329449", "WT reference", "rep1", _aliases("Arabidopsis_WT_rep1", "Arabidopsis_WT_result1")),
    Sample("Arabidopsis_WT_rep2", PLATFORM_RNA002, _AT, "Arabidopsis_WT", "WT", "whole seedling",
           "SRP329449", "WT reference", "rep2", _aliases("Arabidopsis_WT_rep2", "Arabidopsis_WT_result2")),
    Sample("Arabidopsis_WT_rep3", PLATFORM_RNA002, _AT, "Arabidopsis_WT", "WT", "whole seedling",
           "SRP329449", "WT reference", "rep3", _aliases("Arabidopsis_WT_rep3", "Arabidopsis_WT_result3")),
    Sample("Arabidopsis_fip37_rep1", PLATFORM_RNA002, _AT, "Arabidopsis_KD", "KD", "FIP37 knockdown",
           "SRP363914", "partial negative", "rep1", _aliases("Arabidopsis_fip37_rep1", "Arabidopsis_fip37_result1")),
    Sample("Arabidopsis_fip37_rep2", PLATFORM_RNA002, _AT, "Arabidopsis_KD", "KD", "FIP37 knockdown",
           "SRP363914", "partial negative", "rep2", _aliases("Arabidopsis_fip37_rep2", "Arabidopsis_fip37_result2")),
    Sample("Arabidopsis_fip37_rep3", PLATFORM_RNA002, _AT, "Arabidopsis_KD", "KD", "FIP37 knockdown",
           "SRP363914", "partial negative", "rep3", _aliases("Arabidopsis_fip37_rep3", "Arabidopsis_fip37_result3")),
    # ---------------- Mouse (RNA002, two independent studies) ----------------
    Sample("mES_WT", PLATFORM_RNA002, _MO, "Mouse_WT", "WT", "mESC (SRP357195)",
           "SRP357195", "WT reference", "studyB", _aliases("mES_WT", "mES_WT_result")),
    Sample("mES_KO", PLATFORM_RNA002, _MO, "Mouse_KO", "KO", "Mettl3 KO mESC (SRP357195)",
           "SRP357195", "partial negative", "studyB", _aliases("mES_KO", "mES_KO_result")),
    Sample("mESCs_Mettl3_WT", PLATFORM_RNA002, _MO, "Mouse_WT", "WT", "mESC (SRP166020)",
           "SRP166020", "WT reference", "studyA", _aliases("mESCs_Mettl3_WT", "mESCs_Mettl3_WT_result", "mESCs_Mettl3_WT_result1")),
    Sample("mESCs_Mettl3_KO", PLATFORM_RNA002, _MO, "Mouse_KO", "KO", "Mettl3 KO mESC (SRP166020)",
           "SRP166020", "partial negative", "studyA",
           _aliases("mESCs_Mettl3_KO", "mESCs_Mettl3_KO_result", "mESCs_Mettl3_KO_result1")),
    # ---------------- E.coli (RNA002) ----------------
    Sample("E_ss_rd_RNA1", PLATFORM_RNA002, _EC, "E.coli_WT", "WT", "rRNA-depleted polyA RNA",
           "SRP478171", "WT reference", "run1",
           _aliases("E_ss_rd_RNA1", "E_ss_rd_RNA1_result", "E_ss_rd_RNA_result1")),
    Sample("E_ss_rd_RNA2", PLATFORM_RNA002, _EC, "E.coli_WT", "WT", "rRNA-depleted polyA RNA",
           "SRP478171", "WT reference", "run2",
           _aliases("E_ss_rd_RNA2", "E_ss_rd_RNA2_result", "E_ss_rd_RNA_result2")),
    Sample("E_IVT_neg1", PLATFORM_RNA002, _EC, "E.coli_IVT", "IVT", "IVT control",
           "SRP478171", "negative control", "run1", _aliases("E_IVT_neg1", "E_IVT_neg1_result")),
    Sample("E_IVT_neg2", PLATFORM_RNA002, _EC, "E.coli_IVT", "IVT", "IVT control",
           "SRP478171", "negative control", "run2", _aliases("E_IVT_neg2", "E_IVT_neg2_result")),
    # ---------------- Curlcake (RNA002 synthetic constructs) ----------------
    # Semantic library names (user-confirmed 2026-09-15), identical to the
    # directory names in 02_raw_results/result/<tool>/ and
    # 04_revision_analysis/converted_callsets/<tool>/.  The SRA library accessions are kept
    # as aliases: they are what the legacy conversion notebooks named their
    # outputs, and they stay the traceability handle in sample_metadata/tables.
    Sample("Curlcake_IVT_rep1", PLATFORM_RNA002, _CC, "Curlcake_IVT", "IVT",
           "Curkcake IVT", "SRP174366", "negative control", "rep1",
           _aliases("Curlcake_IVT_rep1", "RNAAB089716")),
    Sample("Curlcake_IVT_rep2_partial", PLATFORM_RNA002, _CC, "Curlcake_IVT", "IVT",
           "Curkcake IVT (subset run)", "SRP174366", "negative control", "rep2",
           _aliases("Curlcake_IVT_rep2_partial", "RNA081120181part"),
           note="depth-matched subset of the Curlcake_IVT_rep3 run "
                "(SRR8767348, 49,538 reads) - NOT an independent replicate"),
    Sample("Curlcake_IVT_rep3", PLATFORM_RNA002, _CC, "Curlcake_IVT", "IVT",
           "Curkcake IVT (all)", "SRP174366", "negative control", "rep3",
           _aliases("Curlcake_IVT_rep3", "RNA081120181", "RNA081120181all")),
    Sample("Curlcake_m6A_rep1", PLATFORM_RNA002, _CC, "Curlcake_m6A", "modified",
           "Curkcake fully m6A-modified", "SRP174366", "positive control", "rep1",
           _aliases("Curlcake_m6A_rep1", "RNAAB090763")),
    Sample("Curlcake_m6A_rep2", PLATFORM_RNA002, _CC, "Curlcake_m6A", "modified",
           "Curkcake fully m6A-modified", "SRP174366", "positive control", "rep2",
           _aliases("Curlcake_m6A_rep2", "RNA081120182")),
    # ---------------- RNA004 ----------------
    Sample("HeLa_RNA004_WT", PLATFORM_RNA004, _HE, "RNA004_HeLa_WT", "WT", "HeLa cell line",
           "ERP164259", "WT reference (RNA004)", "rep1", _aliases("HeLa_WT", "HeLa_RNA004_WT")),
    Sample("HeLa_RNA004_IVT", PLATFORM_RNA004, _HE, "RNA004_HeLa_IVT", "IVT", "HeLa IVT control",
           "ERP164259", "negative control (RNA004)", "rep1", _aliases("HeLa_IVT", "HeLa_RNA004_IVT")),
    Sample("Curlcake_RNA004_IVT", PLATFORM_RNA004, _CC, "RNA004_Curlcake_IVT", "IVT",
           "Curkcake IVT (RNA004)", "ERP162788", "negative control (RNA004)", "rep1",
           # No bare "Curlcake" alias: it is a prefix of every RNA002 Curlcake
           # library and hijacked those directories in the longest-prefix fallback.
           _aliases("Curlcake_RNA004_IVT", "Curlcake_RNA004")),
]

SAMPLES_BY_NAME: dict[str, Sample] = {s.canonical: s for s in SAMPLES}

# --------------------------------------------------------------------------- #
# Tools
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ToolSpec:
    """One tool (or tool+modality) and how to find its final per-sample file."""

    tool: str                 # internal id == callset file prefix used by the tools
    display: str              # publication label
    mod_type: str             # m6A / m5C / Psi / m1Psi / inosine / Nm
    platform: str
    result_subdir: str        # directory under result/ (RNA002) or result_RNA004/
    patterns: tuple[str, ...]  # glob patterns, relative to the sample directory
    parser: str = "standard"   # parser key (see parsers/)
    genomic: bool = True       # False => transcript-space, excluded from genomic eval
    note: str = ""
    #: extra (glob, parser) fallbacks tried after ``patterns`` -- used when the
    #: converted file was never produced (e.g. f5c-mode Nanom6A samples).
    raw_fallbacks: tuple[tuple[str, str], ...] = ()
    #: name tokens used by the generic ranked file search (registry.py);
    #: empty => derived from ``tool``.
    tokens: tuple[str, ...] = ()
    #: file-name suffixes that are KNOWN to be transcript-space (pre-liftover)
    #: and must therefore never be used as a genomic callset.
    exclude_suffixes: tuple[str, ...] = ()


def _t(tool: str, display: str, mod: str, subdir: str, patterns: tuple[str, ...],
       parser: str = "standard", genomic: bool = True, note: str = "",
       raw_fallbacks: tuple[tuple[str, str], ...] = (),
       tokens: tuple[str, ...] = (),
       exclude_suffixes: tuple[str, ...] = ()) -> ToolSpec:
    return ToolSpec(tool, display, mod, PLATFORM_RNA002, subdir, patterns, parser,
                    genomic, note, raw_fallbacks, tokens, exclude_suffixes)


RNA002_TOOLS: list[ToolSpec] = [
    _t("CHEUI_m6A", "CHEUI_m6A", "m6A", "CHEUI",
       ("*_CHEUI_m6A_remove_chr.txt",), exclude_suffixes=("_CHEUI_m6A.txt",)),
    _t("CHEUI_m5C", "CHEUI_m5C", "m5C", "CHEUI",
       ("*_CHEUI_m5C_remove_chr.txt",), exclude_suffixes=("_CHEUI_m5C.txt",)),
    _t("CHEUI-diff_m6A", "CHEUI-diff_m6A", "m6A", "CHEUI-diff",
       ("*_CHEUI-diff_m6A.txt", "*_CHEUI_diff_m6A_*.txt"), genomic=False,
       note="transcript (ENST) coordinates only; excluded from genomic evaluation"),
    _t("CHEUI-diff_m5C", "CHEUI-diff_m5C", "m5C", "CHEUI-diff",
       ("*_CHEUI-diff_m5C.txt", "*_CHEUI_diff_m5C_*.txt"), genomic=False,
       note="transcript (ENST) coordinates only; excluded from genomic evaluation"),
    _t("DENA", "DENA", "m6A", "DENA", ("*_DENA_remove_chr.txt",),
       exclude_suffixes=("_DENA.txt",)),
    _t("differr", "differr", "m6A", "differr", ("*_differr_remove_chr.txt", "*_differr.txt")),
    _t("DRUMMER", "DRUMMER", "m6A", "DRUMMER",
       ("*_DRUMMER_remove_chr.txt",), exclude_suffixes=("_DRUMMER.txt",)),
    _t("ELIGOS2_diff", "ELIGOS2_diff", "m6A", "ELIGOS2_diff",
       ("*_ELIGOS2_diff_remove_chr.txt", "*_combine.txt"), parser="eligos2",
       exclude_suffixes=("_ELIGOS2_diff.txt",),
       note="combine.txt fallback applies the same legacy filter as ELIGOS2_solo"),
    _t("ELIGOS2_solo", "ELIGOS2_solo", "m6A", "ELIGOS2_solo",
       ("*_ELIGOS2_solo_remove_chr.txt", "*_combine.txt"), parser="eligos2",
       exclude_suffixes=("_ELIGOS2_solo.txt",),
       note="combine.txt fallback applies the legacy filter ref==A & total_reads>20 & adjPval<1e-4 & oddR>1.2"),
    #: subdir spelling must match the REAL directory (the 2026-09-15 tidy renamed
    #: ``Epinano_DiffErr`` -> ``EpiNano_DiffErr``); with the old spelling the
    #: registry skipped the tool entirely and 16 EpiNano pairs never appeared in
    #: ``sample_tool_registry.csv`` (2026-09-16 fix).
    _t("EpiNano_Error", "EpiNano_Error", "m6A", "EpiNano_DiffErr",
       ("*_EpiNano_Error_remove_chr.txt", "*_EpiNano_Error.txt")),
    _t("EpiNano_SVM", "EpiNano_SVM", "m6A", "Epinano_SVM",
       ("*_EpiNano_SVM_remove_chr.txt", "*_EpiNano_SVM.txt")),
    _t("m6Anet", "m6Anet", "m6A", "m6Anet", ("*_m6Anet_remove_chr.txt",),
       exclude_suffixes=("_m6Anet.txt",)),
    _t("mAFiA", "mAFiA", "m6A", "mAFiA", ("out_dir/mAFiA.sites.bed",), parser="mafia_bed",
       note="Bed with coverage/modRatio columns; modRatio is a COUNT not a ratio"),
    _t("MINES", "MINES", "m6A", "MINES", ("*_MINES_remove_chr.txt",),
       exclude_suffixes=("_MINES.txt",)),
    _t("Nanocompore", "Nanocompore", "m6A", "Nanocompore",
       ("*_Nanocompore_remove_chr.txt", "*_Nanocompore.txt")),
    _t("Nanom6A", "Nanom6A", "m6A", "Nanom6A", ("*_Nanom6A_remove_chr.txt", "*_Nanom6A.txt"),
       #: 2026-09-18: the four f5c-filled samples (Arabidopsis_WT_rep1/2,
       #: mESCs_Mettl3_WT/KO) had their site anchors written 7 bp upstream of the
       #: modified A by the ``--legacy-binary`` converter (kmer =
       #: seq[indx+5..8]+seq[indx+3], A at indx+7, key = indx); the salvaged
       #: ``ratio.0.5.anchorfix.tsv`` (read-frame +7 remap, same aggregation) comes
       #: first so it wins over the un-rectified ``ratio.0.5.tsv``.
       
       raw_fallbacks=(("ratio.0.5.anchorfix.tsv", "nanom6a_ratio_tsv"),
                      ("ratio.0.5.tsv", "nanom6a_ratio_tsv")),
       note="fallback parses raw ratio.0.5[.anchorfix].tsv (f5c-mode fill samples)"),
    _t("NanoMUD_psi", "NanoMUD-\u03a8", "Psi", "NanoMUD",
       ("*_NanoMUD_psi_remove_chr.txt", "*_NanoMUD_psi.txt")),
    _t("NanoMUD_m1psi", "NanoMUD-m1\u03a8", "m1Psi", "NanoMUD",
       ("*_NanoMUD_m1psi_remove_chr.txt", "*_NanoMUD_m1psi.txt")),
    _t("NanoNm", "NanoNm", "Nm", "NanoNm", ("*_NanoNm_remove_chr.txt", "*_NanoNm.txt")),
    _t("NanoPsu", "NanoPsu", "Psi", "NanoPsu",
       ("*_NanoPsu_5mer.txt", "*_NanoPsu_remove_chr.txt", "*_NanoPsu.txt")),
    _t("NanoSPA_m6A", "NanoSPA_m6A", "m6A", "NanoSPA",
       ("*_NanoSPA_m6A_5mer.txt", "*_NanoSPA_m6A_remove_chr.txt", "*_NanoSPA_m6A.txt")),
    _t("NanoSPA_psU", "NanoSPA-\u03a8", "Psi", "NanoSPA",
       ("*_NanoSPA_psU_5mer.txt", "*_NanoSPA_psU_remove_chr.txt", "*_NanoSPA_psU.txt")),
    _t("Tombo_com", "Tombo", "m6A", "Tombo_com",
       ("*_Tombo_com_remove_chr.txt",), exclude_suffixes=("_Tombo_com.txt",)),
    _t("xPore", "xPore", "m6A", "xPore",
       ("*_xPore_remove_chr.txt",), exclude_suffixes=("_xPore.txt", "_xPore_unfiltered.txt")),
    _t("yanocomp", "Yanocomp", "m6A", "yanocomp",
       ("*_yanocomp_remove_chr.txt", "*_output.bed"), parser="yanocomp_bed",
       exclude_suffixes=("_yanocomp.txt",),
       note="output.bed fallback: middle base of the 5-mer, Start=bed_start+2"),
]

#: Tools in result/ that never produced a genomic callset (documented, skipped).
SKIPPED_TOOLS = {
    "Tombo": "only wig/stats in transcript coordinates, no callsets were ever produced",
}

# --------------------------------------------------------------------------- #
# RNA004 sources (explicit, the tree is hand-curated)
# --------------------------------------------------------------------------- #
#: one declared RNA004 source: ``(family, mod_type, pattern, parser[, kwargs])``.
#: The optional 5th element is a dict of extra parser kwargs -- currently the
#: Dorado pileup call filter (``{"min_cov": …, "min_pct": …}``), see
#: ``parsers.parse_dorado_pileup``.
RNA004SourceSpec = tuple[str, str, str, str] | tuple[str, str, str, str, dict]

#: sample canonical -> list of declared RNA004 sources
RNA004_SOURCES: dict[str, list[RNA004SourceSpec]] = {
    "HeLa_RNA004_WT": [
        # Dorado entries are families: each matching pileup file becomes its own
        # callset labelled ``Dorado_<hac|sup>_<model>@<version>``.
        # ``mod_type="auto"``: the modification is NOT declared here but read from
        # the pileup itself -- modkit writes the modification code into its ``name``
        # column (``a``=m6A, ``m``=m5C, ``17802``=Psi, ``17596``=inosine, see
        # ``common/parsers.MODKIT_CODE``).  A file may even mix several codes (that
        # was the case for Curlcake), so ``registry.resolve_rna004_sources`` emits
        # one callset per code and the parser keeps only that code's rows.
        # ``m6A_guitar/`` is deliberately NOT declared (2026-09-17): those files are
        # header-less 3-column BEDs, so ``parse_dorado_pileup`` could never read
        # them and every ``*_drachGuitar`` callset was silently empty (0 rows out
        # of 2 k - 52 k source rows).  The channel is not part of the manuscript.
        ("Dorado", "auto", "raw_calls/RNA004_result/dorado_model_split/m6A/HeLa_WT_*_pileup.bed", "dorado_pileup"),
        ("Dorado", "auto", "raw_calls/RNA004_result/dorado_model_split/other_modification/HeLa_WT_*_pileup.bed", "dorado_pileup"),
        ("m6Anet", "m6A", "raw_calls/RNA004_result/m6Anet/HeLa_WT_RNA004_m6Anet_remove_chr.txt", "standard"),
        ("NanoSPA_m6A", "m6A", "raw_calls/RNA004_result/NanoSPA/HeLa_WT_NanoSPA_m6A_5mer.txt", "standard"),
        ("NanoSPA_psU", "Psi", "raw_calls/RNA004_result/NanoSPA/HeLa_WT_NanoSPA_psU_5mer.txt", "standard"),
        ("NanoPsu", "Psi", "raw_calls/RNA004_result/NanoPsu/HeLa_WT_NanoPsu_5mer.txt", "standard"),
        ("ELIGOS2_solo", "m6A", "raw_calls/RNA004_result/ELIGOS2_solo/HeLa_WT_RNA004_ELIGOS2_solo_remove_chr_processed.txt", "standard"),
        ("ELIGOS2_diff", "m6A", "raw_calls/RNA004_result/ELIGOS2_diff/HeLa_WT_RNA004_ELIGOS2_diff_remove_chr_processed.txt", "standard"),
        ("DRUMMER", "m6A", "raw_calls/RNA004_result/DRUMMER/HeLa_WT_RNA004_DRUMMER_remove_chr.txt", "standard"),
        # mAFiA never produced RNA004 sites (documented gap, R3-10).
        ("mAFiA", "m6A", "raw_calls/RNA004_result/mAFiA/*/out_dir/mAFiA.sites.bed", "mafia_bed"),
    ],
    "HeLa_RNA004_IVT": [
        ("Dorado", "auto", "raw_calls/RNA004_result/dorado_model_split/m6A/HeLa_IVT_*_pileup.bed", "dorado_pileup"),
        ("Dorado", "auto", "raw_calls/RNA004_result/dorado_model_split/other_modification/HeLa_IVT_*_pileup.bed", "dorado_pileup"),
        ("m6Anet", "m6A", "raw_calls/RNA004_result/m6Anet/HeLa_IVT_RNA004_m6Anet_remove_chr.txt", "standard"),
        ("NanoSPA_m6A", "m6A", "raw_calls/RNA004_result/NanoSPA/HeLa_IVT_NanoSPA_m6A_5mer.txt", "standard"),
        ("NanoSPA_psU", "Psi", "raw_calls/RNA004_result/NanoSPA/HeLa_IVT_NanoSPA_psU_5mer.txt", "standard"),
        ("NanoPsu", "Psi", "raw_calls/RNA004_result/NanoPsu/HeLa_IVT_NanoPsu_5mer.txt", "standard"),
        ("ELIGOS2_solo", "m6A", "raw_calls/RNA004_result/ELIGOS2_solo/HeLa_IVT_RNA004_ELIGOS2_solo_remove_chr_processed.txt", "standard"),
    ],
    "Curlcake_RNA004_IVT": [
        # available models: m6A_DRACH, pseU_m6A, inosine_m6A, all (no m5C pileup).
        # ``auto`` matters here: these pileups were never split by modification --
        # ``*_all_pileup.bed`` carries ``17802``/``a``/``m`` in one file and
        # ``*_pseU_m6A_pileup.bed`` carries ``17802``/``a`` (62 % of the rows the
        # old declaration labelled "m6A" were Psi / m5C / inosine).
        #: This is the only Dorado family whose files were never *called*-filtered:
        #: they hold one row per covered position x modification code, most of them
        #: no-call rows (``percent_modified = 0``).  Without a filter the callsets
        #: collapsed onto the construct base composition (m6A/m5C/Psi all ~58 %
        #: "expected base", i.e. random).  The HeLa family was already filtered at
        #: source (``dorado_model_split/``: every row cov >= 20 and pct > 0), so the
        #: filter is declared here only (2026-09-18).  ``percent_modified`` is kept
        #: as the score, so ``09_eval_rna004``'s 5/10/20/50 % scan is unchanged.
        ("Dorado", "auto", "RNA004_result/dorado_model/*_pileup.bed", "dorado_pileup",
         {"min_cov": 20, "min_pct": 0.0}),
        ("m6Anet", "m6A", "RNA004_result/m6Anet/data.site_proba.csv", "m6anet_csv"),
        ("NanoSPA_m6A", "m6A", "RNA004_result/NanoSPA/prediction_m6A.csv", "nanospa_csv"),
        ("NanoSPA_psU", "Psi", "RNA004_result/NanoSPA/prediction_psU.csv", "nanospa_csv"),
        ("NanoPsu", "Psi", "RNA004_result/NanoPsu/prediction.csv", "nanospa_csv"),
    ],
}

# --------------------------------------------------------------------------- #
# Coverage BAMs
# --------------------------------------------------------------------------- #
def nanom6a_bam(sample: Sample) -> Path | None:
    """Genomic BAM produced by the Nanom6A runs (same basecalls as callsets)."""
    for alias in sample.aliases:
        p = (_XB / "raw/result/Nanom6A") / alias / "extract.sort.bam"
        if p.exists():
            return p
    return None




_FASTQ_BACKUP_GROUP: dict[str, str] = {
    "HeLa_WT": "HeLa_WT",
    "HeLa_IVT": "HeLa_IVT",
    "Arabidopsis_WT": "Arabidopsis_WT",
    "Arabidopsis_KD": "Arabidopsis_fip37",
    "E.coli_WT": "E_coli",
    "E.coli_IVT": "E_coli",
}
_FASTQ_BACKUP_MOUSE: dict[str, str] = {
    "SRP357195": "mES",
    "SRP166020": "mESCs_Mettl3",
}


def fastq_backup_bam(sample: Sample) -> Path | None:
    """Legacy guppy alignment in fastq_backup (genomic, *G_sort.bam)."""
    if sample.dataset_group.startswith("Mouse"):
        group = _FASTQ_BACKUP_MOUSE.get(sample.study)
    else:
        group = _FASTQ_BACKUP_GROUP.get(sample.dataset_group)
    if group is None:
        return None
    base = (_XB / "raw/fastq_backup") / group / sample.canonical / "fastq"
    if not base.is_dir():
        return None
    hits = sorted(base.glob("*G_sort.bam")) or sorted(base.glob("*G.bam"))
    return hits[0] if hits else None


def nanopore_minimap2_bam(sample: Sample) -> Path | None:
    """Current-data minimap2 genomic BAM under $RNAMODBENCH_LOCAL/nanopore/data (read-only)."""
    roots = [NANOPORE_DATA / sample.study, (_XB / "nanopore/data/other") / sample.study]
    for root in roots:
        d = root / sample.canonical / "minimap2"
        if not d.is_dir():
            continue
        hits = sorted(d.glob("*G_sort.bam")) or sorted(d.glob("*G.bam"))
        if hits:
            return hits[0]
    return None


def resolve_coverage_bam(sample: Sample) -> tuple[Path | None, str]:
    """Best available genomic BAM for coverage annotation (first hit wins)."""
    for fn, label in ((nanom6a_bam, "nanom6a_bam"),
                      (fastq_backup_bam, "fastq_backup_bam"),
                      (nanopore_minimap2_bam, "nanopore_minimap2_bam")):
        p = fn(sample)
        if p is not None:
            return p, label
    return None, "missing"


# --------------------------------------------------------------------------- #
# RNA004 coverage BAMs (Dorado basecalled reads aligned with minimap2)
# --------------------------------------------------------------------------- #
def rna004_coverage_bam(sample: Sample) -> Path | None:
    """Genomic BAM for RNA004 samples (read-only nanopore tree).

    ERP164259 (HeLa RNA004) has minimap2 ``*_allG_sort.bam``; the Curlcake
    RNA004 library (ERP162788) is not on disk, so its coverage comes from the
    modkit pileup ``valid_coverage`` instead.
    """
    sub = {"HeLa_RNA004_WT": "HeLa_WT", "HeLa_RNA004_IVT": "HeLa_IVT"}.get(sample.canonical)
    if sub:
        d = (_XB / "nanopore/data/ERP164259") / sub / "minimap2"
        if d.is_dir():
            hits = sorted(d.glob("*_allG_sort.bam")) or sorted(d.glob("*G_sort.bam"))
            if hits:
                return hits[0]
    return None


def sample_result_dirs(sample: Sample) -> list[Path]:
    """Legacy helper kept for reconciliation scripts."""
    return [RESULT_RNA002 / sample.canonical]


def callset_dir(sample: Sample, mod_type: str, tool: str) -> Path:
    """Where a standardised callset is written (per sample == per replicate)."""
    return CALLSET_ROOT / sample.platform / sample.species / sample.dataset_group / mod_type / tool


def universe_path(sample: Sample, mod_type: str) -> Path:
    return UNIVERSE_ROOT / sample.platform / sample.species / f"{sample.canonical}__{mod_type}.tsv"
