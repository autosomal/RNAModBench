#!/usr/bin/env python
"""Shared helpers for the R2-1 coverage-stratified rank-stability analysis.

New analysis, 2026-09-25 (reviewer R2-1: "how GLORI-site distribution and
transcript coverage may affect the comparisons and reported rankings").
Nothing here writes outside
``analysis/coverage_rank_stability/``.

Locked (mirrors the manuscript's primary benchmark, see
``figures/figure5/README.md`` and ``harmonisation/README.md`` §3)
-------------------------------------------------------------------------------
* one row = one **independent sequencing unit** (never a merged call set);
* candidate universe ``U(unit, c = 10)``: exonic, reference-base compatible,
  coverage >= 10 in the unit's **own** alignment;
* the stratifier is the unit's own per-position coverage (a property of the
  sequencing library, *not* of a tool): the tool-reported coverage columns
  (NanoNm ``src_coverage``, Dorado modkit ``N_valid_cov``) are deliberately
  unused, otherwise a tool's own filtering would define its stratum;
* matching window w = 2 bp (``config.PRIMARY_WINDOW``), TP/FP/FN/TN inside the
  stratum's own universe: ``U_b`` (stratum positions), ``G_b`` (GLORI sites in
  ``U_b``), ``C_b`` (tool calls at ``U_b`` positions);
* ``pooled>=10`` (= ``U``, the Fig. 5A/5C/5D universe) must reproduce
  ``harmonisation/evaluation/tables/m6a_glori_confusion.tsv`` (w = 2, min_cov = 10)
  exactly -- enforced by ``r21_verify.py``;
* ``5-9`` is a **sensitivity stratum outside the primary universe**; it is never
  pooled with the primary strata and is flagged ``stratum_kind="sensitivity"``.

Modules are imported read-only from ``src/harmonisation`` (``common.*``);
no file under ``$RNAMODBENCH_LOCAL/code`` is modified.
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
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ANALYSIS_DIR = Path(__file__).resolve().parent
OUT = ANALYSIS_DIR.parent
TAB = OUT / "tables"
FIG = OUT / "figures"
LOG = OUT / "logs"

SITES_CODE = Path(str(_RB / "src/harmonisation"))
if str(SITES_CODE) not in sys.path:
    sys.path.insert(0, str(SITES_CODE))

from common import config as C                                       # noqa: E402
from common.evaluation import MOD_GENOME_BASES, confusion_windows    # noqa: E402
from common.evaluation import load_reference, reference_in_universe  # noqa: E402
from common.io_utils import read_table, write_table                  # noqa: E402
from common.manifest import setup_logger                             # noqa: E402
from common.match import fix_chromosome                              # noqa: E402
from common.metrics import window_hit_mask                           # noqa: E402

MOD_TYPE = "m6A"
WINDOW = C.PRIMARY_WINDOW          # 2 bp: the manuscript's primary window
MIN_COV = C.C_MIN_DEFAULT          # 10: the manuscript's candidate universe

# --------------------------------------------------------------------------- #
# strata
# --------------------------------------------------------------------------- #
#: pooled = the manuscript's candidate universe; the three primary strata
#: partition it exactly (universe / reference / call counts are additive).
POOLED = "pooled>=10"
MAIN_STRATA = ["10-19", "20-49", ">=50"]
#: sensitivity stratum OUTSIDE the primary universe (coverage 5-9 positions are
#: stored in the universe file but excluded from U by the c >= 10 rule).
SENS_STRATUM = "5-9"
SENS_COV_RANGE = (5, 10)
STRATUM_KIND = {POOLED: "pooled", SENS_STRATUM: "sensitivity",
                **{s: "primary" for s in MAIN_STRATA}}
#: [low, high) coverage of each stratum; None = open ended.
BAND_EDGES: dict[str, tuple[int, int | None]] = {
    SENS_STRATUM: (5, 10), "10-19": (10, 20), "20-49": (20, 50), ">=50": (50, None),
}
#: coverage code used while streaming the universe file (np.searchsorted order).
CODE_LABEL = {0: SENS_STRATUM, 1: "10-19", 2: "20-49", 3: ">=50"}
_CODE_EDGES = np.array([10, 20, 50], dtype=np.int64)
ALL_STRATA = [POOLED, SENS_STRATUM, *MAIN_STRATA]

#: coverage edges of the primary analysis; ``load_banded_universe`` accepts any
#: ascending edge list, which is what the depth sweep (r21b) uses.
DEFAULT_EDGES: tuple[int, ...] = (10, 20, 50)


def pooled_label(edges=DEFAULT_EDGES) -> str:
    """Label of the cov >= edges[0] universe (the manuscript's candidate universe)."""
    return f"pooled>={list(edges)[0]}"


def bin_labels(edges=DEFAULT_EDGES) -> list[str]:
    """Mutually exclusive coverage bins, one per ``searchsorted`` code.

    ``(10, 20, 50)`` -> ``['5-9', '10-19', '20-49', '>=50']`` (the primary strata);
    ``(10, 15, 20, 30, 50, 100)`` -> the seven fine bins of the depth sweep.
    """
    e = list(edges)
    return [f"5-{e[0] - 1}"] + [f"{a}-{b - 1}" for a, b in zip(e, e[1:])] + [f">={e[-1]}"]


def threshold_labels(edges=DEFAULT_EDGES) -> list[str]:
    """Cumulative universes ``cov >= e``; the first edge is the pooled universe."""
    return [f"cov>={e}" for e in list(edges)[1:]]


def strata_for(edges=DEFAULT_EDGES) -> list[str]:
    """Every stratum label produced by ``load_banded_universe(path, edges)``."""
    return [pooled_label(edges)] + bin_labels(edges) + threshold_labels(edges)


def stratum_kind(label: str) -> str:
    """Role of a stratum label (for tables / plotting / documentation)."""
    if label == POOLED or label.startswith("pooled>="):
        return "pooled"
    if label.startswith("cov>="):
        return "threshold"
    if label == SENS_STRATUM:
        return "sensitivity"
    if label in MAIN_STRATA:
        return "primary"
    return "fine"

METRICS = ["precision", "recall", "f1", "mcc"]
METRIC_LABEL = {"precision": "PPV", "recall": "Recall", "f1": "F1", "mcc": "MCC"}

#: publication tool list / order (identical to scripts 35 / 36)
TOOL_ORDER = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
              "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
              "NanoSPA_m6A", "xPore", "yanocomp"]
TOOL_DISPLAY = {"yanocomp": "Yanocomp"}

#: (species, dataset_group, [independent sequencing units drawn in the figures])
#: identical to Fig. 5A/5B/5C-D (mouse = the single study mES_WT).
PANELS: list[tuple[str, str, list[str]]] = [
    ("Arabidopsis", "Arabidopsis_WT",
     ["Arabidopsis_WT_rep1", "Arabidopsis_WT_rep2", "Arabidopsis_WT_rep3"]),
    ("Mouse", "Mouse_WT", ["mES_WT"]),
    ("Human", "HeLa_WT", ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"]),
]
#: tabulated-only units: KD / IVT / the second mouse study (never averaged into
#: the drawn panels; used by the sensitivity section).  The two mouse KO units
#: come from two independent studies -- reported as such, never as replicates.
EXTRA_GROUPS: list[tuple[str, str, list[str]]] = [
    ("Arabidopsis", "Arabidopsis_KD",
     ["Arabidopsis_fip37_rep1", "Arabidopsis_fip37_rep2", "Arabidopsis_fip37_rep3"]),
    ("Mouse", "Mouse_WT", ["mESCs_Mettl3_WT"]),
    ("Mouse", "Mouse_KO", ["mES_KO", "mESCs_Mettl3_KO"]),
    ("Human", "HeLa_IVT", ["HeLa_IVT_rep1", "HeLa_IVT_rep2", "HeLa_IVT_rep3"]),
]
#: every RNA002 m6A unit of the three GLORI species that has a universe file
PRIMARY_SPECIES = ["Arabidopsis", "Mouse", "Human"]
ALL_GROUPS: list[tuple[str, str]] = [
    ("Arabidopsis", "Arabidopsis_WT"), ("Arabidopsis", "Arabidopsis_KD"),
    ("Mouse", "Mouse_WT"), ("Mouse", "Mouse_KO"),
    ("Human", "HeLa_WT"), ("Human", "HeLa_IVT"),
]

BOOTSTRAP_B = C.BOOTSTRAP_B        # 1000
SEED = C.SEED                      # 20260914


# --------------------------------------------------------------------------- #
# paths
# --------------------------------------------------------------------------- #
def universe_file(sample) -> Path:
    """Generic per-unit universe file (plain or gzipped), exactly what script 06 read."""
    plain = (C.UNIVERSE_ROOT / sample.platform / sample.species
             / f"{sample.canonical}__universe.tsv")
    gz = plain.with_suffix(".tsv.gz")
    return plain if plain.exists() else gz


def callset_path(sample, tool: str) -> Path:
    return (C.CALLSET_ROOT / sample.platform / sample.species / sample.dataset_group
            / MOD_TYPE / tool / f"{sample.canonical}.tsv")


def units_for_groups(logger=None):
    """[(sample, species, group)] for every unit with a universe file."""
    out = []
    for species, group in ALL_GROUPS:
        for s in C.SAMPLES:
            if (s.platform == C.PLATFORM_RNA002 and s.species == species
                    and s.dataset_group == group):
                up = universe_file(s)
                if not up.exists():
                    if logger:
                        logger.warning("[%s] no universe file (%s), skipped",
                                       s.canonical, up.name)
                    continue
                out.append(s)
    return out


# --------------------------------------------------------------------------- #
# banded universe
# --------------------------------------------------------------------------- #
@dataclass
class BandedUniverse:
    """One unit's m6A universe split by its own per-position coverage.

    With the default edges ``strata["pooled>=10"]`` is the manuscript universe
    U(c = 10), the three primary bins partition it and ``strata["5-9"]`` is
    disjoint (sensitivity); ``load_banded_universe(edges=...)`` adds one key per
    fine bin plus one per cumulative threshold universe ``cov >= e``.
    """

    strata: dict[str, dict[str, np.ndarray]] = field(default_factory=dict)
    #: coverage >= 5 rows read (all A/T bases), for logging / additivity checks
    n_rows_read: int = 0
    #: in_glori counts as stored in the universe file, per stratum
    n_in_glori_file: dict[str, int] = field(default_factory=dict)
    #: label -> coverage codes it is built from (bins = one code, pooled /
    #: thresholds = several); needed to split call sets generically
    label_codes: dict[str, tuple[int, ...]] = field(default_factory=dict)
    #: coverage code -> bin label (one entry per code, i.e. per fine bin)
    code_labels: dict[int, str] = field(default_factory=dict)

    @property
    def n_pooled(self) -> int:
        return int(sum(v.size for v in self.strata[POOLED].values()))


def _norm_chroms(values: pd.Series, cache: dict[str, str]) -> pd.Series:
    """Vectorised ``fix_chromosome`` (per unique label, identical semantics)."""
    uniq = pd.unique(values)
    mapping = {u: cache.setdefault(u, fix_chromosome(u)) for u in uniq}
    return values.map(mapping)


def load_banded_universe(path: Path, chunk_rows: int = 4_000_000,
                         edges=DEFAULT_EDGES) -> BandedUniverse:
    """Stream a universe file once and split A/T positions into coverage strata.

    Semantics identical to ``common.evaluation.load_universe`` (base filter,
    ``fix_chromosome``, ``np.unique`` per chromosome) plus the coverage split:
    with the default edges ``(10, 20, 50)`` the bins are
    ``5-9 / 10-19 / 20-49 / >=50`` and ``pooled>=10`` is their cov >= 10 union --
    exactly the manuscript universe.  Any ascending ``edges`` list is accepted
    (r21b uses the seven fine bins ``10, 15, 20, 30, 50, 100``), and every
    threshold universe ``cov >= e`` is materialised as well.
    """
    edges_list = [int(e) for e in edges]
    edges_arr = np.array(edges_list, dtype=np.int64)
    labels = bin_labels(edges_list)                    # one label per code
    pooled_lab = pooled_label(edges_list)
    thr_labs = threshold_labels(edges_list)
    want = MOD_GENOME_BASES[MOD_TYPE]           # {'A', 'T'} for m6A
    acc: dict[tuple[int, str], list[np.ndarray]] = {}
    glori_counts = {s: 0 for s in labels}
    cache: dict[str, str] = {}
    n_read = 0

    for chunk in pd.read_csv(path, sep="\t",
                             usecols=["chrom", "pos", "base", "coverage", "in_glori"],
                             dtype={"chrom": str, "base": str, "in_glori": str},
                             chunksize=chunk_rows):
        d = chunk[chunk["base"].isin(want)]
        if d.empty:
            continue
        cov = pd.to_numeric(d["coverage"]).to_numpy(dtype=np.int64)
        pos = pd.to_numeric(d["pos"]).to_numpy(dtype=np.int64)
        gl = (d["in_glori"].to_numpy() == "True")
        chroms = _norm_chroms(d["chrom"], cache).to_numpy()
        code = np.searchsorted(edges_arr, cov, side="right")     # 0..len(edges)
        n_read += int(code.size)
        for c, label in enumerate(labels):
            m = code == c
            if not m.any():
                continue
            glori_counts[label] += int(gl[m].sum())
            sub_chroms = chroms[m]
            sub_pos = pos[m]
            for ch in pd.unique(sub_chroms):
                acc.setdefault((c, ch), []).append(sub_pos[sub_chroms == ch])

    by_code: dict[int, dict[str, np.ndarray]] = {c: {} for c in range(len(labels))}
    for (c, ch), parts in acc.items():
        by_code[c][ch] = np.unique(np.concatenate(parts))

    strata: dict[str, dict[str, np.ndarray]] = {label: {} for label in labels}
    for c, label in enumerate(labels):
        strata[label] = by_code[c]

    def _union(codes: range) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
        """Union of the codes in ``codes`` with the aligned coverage code array."""
        out: dict[str, np.ndarray] = {}
        code_tab: dict[str, np.ndarray] = {}
        chroms = sorted({ch for c in codes for ch in by_code[c]})
        for ch in chroms:
            parts, code_parts = [], []
            for c in codes:
                arr = by_code[c].get(ch)
                if arr is not None and arr.size:
                    parts.append(arr)
                    code_parts.append(np.full(arr.size, c, dtype=np.int8))
            pos = np.concatenate(parts)
            cod = np.concatenate(code_parts)
            order = np.argsort(pos, kind="stable")
            pos, cod = pos[order], cod[order]
            if pos.size > 1:                  # a position cannot be in two bands
                dup = pos[1:] == pos[:-1]
                if dup.any():
                    keep = np.ones(pos.size, dtype=bool)
                    keep[1:] = ~dup
                    pos, cod = pos[keep], cod[keep]
            out[ch] = pos
            code_tab[ch] = cod
        return out, code_tab

    pooled, codes = _union(range(1, len(labels)))          # cov >= edges[0]
    strata[pooled_lab] = pooled
    for k, lab in enumerate(thr_labs):                     # cov >= edges[k+1]
        strata[lab] = _union(range(k + 2, len(labels)))[0]
        glori_counts[lab] = sum(glori_counts[labels[c]]
                                for c in range(k + 2, len(labels)))
    glori_counts[pooled_lab] = sum(glori_counts[labels[c]]
                                   for c in range(1, len(labels)))

    label_codes: dict[str, tuple[int, ...]] = {lab: (c,) for c, lab in enumerate(labels)}
    label_codes[pooled_lab] = tuple(range(1, len(labels)))
    for k, lab in enumerate(thr_labs):
        label_codes[lab] = tuple(range(k + 2, len(labels)))

    #: per-stratum arrays as stored; the pooled codes are kept in a side table so
    #: call sets can be assigned to bands in one searchsorted pass.
    uni = BandedUniverse(strata=strata, n_rows_read=n_read, n_in_glori_file=glori_counts,
                         label_codes=label_codes,
                         code_labels={c: lab for c, lab in enumerate(labels)})
    uni.codes = codes                                  # type: ignore[attr-defined]
    return uni


def pooled_codes(uni: BandedUniverse) -> dict[str, np.ndarray]:
    """Coverage codes aligned to ``uni.strata[POOLED]`` (channel for call assignment)."""
    if not hasattr(uni, "codes"):
        raise RuntimeError("load_banded_universe must be used to build the codes")
    return uni.codes                                   # type: ignore[attr-defined]


# --------------------------------------------------------------------------- #
# calls
# --------------------------------------------------------------------------- #
def load_calls(sample, tool: str) -> pd.DataFrame:
    """Callset rows with a numeric ``pos_raw`` and normalised chromosome labels.

    Mirrors script ``06_eval_m6a_glori.py``: header-only / unreadable files
    (no ``pos_raw`` column) return an empty frame and are skipped by the caller;
    a header-only call set is a real zero-detection result, not a missing sample.
    """
    path = callset_path(sample, tool)
    df = read_table(path)
    if "pos_raw" not in getattr(df, "columns", []):
        return pd.DataFrame(columns=["chrom", "pos_raw"])
    df = df.copy()
    df["pos_raw"] = pd.to_numeric(df["pos_raw"]).astype(np.int64)
    df["chrom"] = _norm_chroms(df["chrom"].astype(str), {})
    return df


def calls_by_stratum(df: pd.DataFrame, uni: BandedUniverse, labels=None
                     ) -> tuple[dict[str, dict[str, np.ndarray]], int, int]:
    """Split a unit's tool calls into the coverage strata of ``labels``.

    Every label is a union of coverage codes (``uni.label_codes``), so the split is
    done once per code and the labels are then assembled by union -- which keeps
    bins, the pooled universe and the cumulative threshold universes consistent by
    construction.  ``labels`` defaults to the primary strata (identical behaviour
    to the original implementation).

    Returns ``({stratum: {chrom: sorted unique call positions}}, n_out, n_in)``
    where ``n_out`` / ``n_in`` are counted against the **pooled** universe only
    (same quantity as ``n_calls_out_of_universe`` in script 06).
    """
    labels = list(labels or ALL_STRATA)
    pooled = uni.strata[POOLED]
    codes = pooled_codes(uni)
    out: dict[str, dict[str, np.ndarray]] = {label: {} for label in labels}
    n_out = n_in = 0
    if df.empty:
        return out, n_out, n_in

    pooled_code_min = min(uni.label_codes[POOLED])
    by_code: dict[int, dict[str, np.ndarray]] = {}

    for chrom, sub in df.groupby("chrom", sort=False):
        pos = np.unique(sub["pos_raw"].to_numpy(dtype=np.int64))
        u = pooled.get(chrom)
        if u is not None and u.size:
            idx = np.clip(np.searchsorted(u, pos), 0, u.size - 1)
            inside = u[idx] == pos
            n_out += int((~inside).sum())
            if inside.any():
                n_in += int(inside.sum())
                calls_in = pos[inside]
                code_in = codes[chrom][np.clip(np.searchsorted(u, calls_in), 0, u.size - 1)]
                for c in np.unique(code_in):
                    by_code.setdefault(int(c), {})[chrom] = calls_in[code_in == c]
        else:
            n_out += int(pos.size)
        # codes below the pooled universe (coverage 5-9 by default): their own bins
        for c in range(pooled_code_min):
            lab = uni.code_labels[c]
            arr = uni.strata[lab].get(chrom)
            if arr is None or not arr.size:
                continue
            idx = np.clip(np.searchsorted(arr, pos), 0, arr.size - 1)
            keep = arr[idx] == pos
            if keep.any():
                by_code.setdefault(c, {})[chrom] = pos[keep]

    for label in labels:
        member_codes = uni.label_codes.get(label, ())
        chroms = sorted({ch for c in member_codes for ch in by_code.get(c, {})})
        for chrom in chroms:
            parts = [by_code[c][chrom] for c in member_codes if by_code.get(c, {}).get(chrom)
                     is not None]
            if parts:
                out[label][chrom] = (parts[0] if len(parts) == 1
                                     else np.unique(np.concatenate(parts)))
    return out, n_out, n_in


# --------------------------------------------------------------------------- #
# per-stratum confusion + metrics
# --------------------------------------------------------------------------- #
def strata_confusion(uni: BandedUniverse, ref: dict[str, np.ndarray],
                     calls: dict[str, dict[str, np.ndarray]], labels=None
                     ) -> dict[str, dict]:
    """TP/FP/FN/TN + metrics for every stratum of one unit x tool.

    Two leak-free cross-checks are added per stratum, where the *matching
    partner* is NOT restricted to the stratum (so they carry no cross-band
    boundary spill, see the README):

    * ``n_calls_hit_global`` -- calls inside the stratum matched at w = 2 bp
      against the **pooled** reference; ``PPV_global = hits / n_calls``;
    * ``n_glori_hit_global`` -- pooled GLORI sites whose own coverage lies in
      the stratum ("reference-side" stratum) detected by *any* call of the unit;
      ``recall_global = hits / |G_stratum|``.
    """
    out: dict[str, dict] = {}
    U_pooled = uni.strata[POOLED]
    G_pooled = reference_in_universe(ref, U_pooled)
    C_pooled = calls.get(POOLED, {})
    for label in list(labels or ALL_STRATA):
        U = uni.strata[label]
        G = reference_in_universe(ref, U)
        Cb = calls.get(label, {})
        conf = confusion_windows(U, G, Cb, [WINDOW])[WINDOW]
        n_calls_hit = 0
        n_glori_hit = 0
        for chrom in set(Cb) | set(G):
            c = Cb.get(chrom, np.zeros(0, dtype=np.int64))
            g = G.get(chrom, np.zeros(0, dtype=np.int64))
            gp = G_pooled.get(chrom, np.zeros(0, dtype=np.int64))
            cp = C_pooled.get(chrom, np.zeros(0, dtype=np.int64))
            if c.size and gp.size:
                n_calls_hit += int(window_hit_mask(c, gp, WINDOW).sum())
            if g.size and cp.size:
                n_glori_hit += int(window_hit_mask(g, cp, WINDOW).sum())
        out[label] = {**conf, "n_glori": int(sum(g.size for g in G.values())),
                      "n_calls_hit_global": n_calls_hit,
                      "n_glori_hit_global": n_glori_hit}
    return out


def metric_columns() -> list[str]:
    return ["tp", "fp", "fn", "tn", "universe", "n_calls_in_universe",
            "n_reference_in_universe", "precision", "recall", "f1", "specificity", "mcc"]


# --------------------------------------------------------------------------- #
# ranking + stability statistics
# --------------------------------------------------------------------------- #
#: two leak-free companions of the strict strata metrics (the matching partner is
#: NOT restricted to the stratum, so cross-band boundary spill cannot leak in):
#: ``precision_global`` = calls in the stratum matched against the pooled GLORI
#: reference; ``recall_global`` = pooled GLORI sites whose own coverage lies in
#: the stratum, detected by any call of the unit.
GLOBAL_METRICS = ["precision_global", "recall_global"]


def add_global_metrics(per_unit: pd.DataFrame) -> pd.DataFrame:
    """Derive the leak-free PPV / recall columns from the ``*_global`` counts.

    Denominators are the definitional ones: calls inside the stratum universe and
    the GLORI sites inside the stratum universe (``n_glori``), never the legacy
    ``in_glori`` annotation column of the universe file.
    """
    out = per_unit.copy()
    calls = pd.to_numeric(out["n_calls_in_universe"], errors="coerce")
    hits = pd.to_numeric(out["n_calls_hit_global"], errors="coerce")
    glori = pd.to_numeric(out["n_glori"], errors="coerce")
    hits_g = pd.to_numeric(out["n_glori_hit_global"], errors="coerce")
    out["precision_global"] = hits / calls.where(calls > 0)
    out["recall_global"] = hits_g / glori.where(glori > 0)
    return out


def rank_within_unit(per_unit: pd.DataFrame, metrics: list[str] | None = None
                     ) -> pd.DataFrame:
    """Add ``<metric>_rank`` within each (sample, stratum): 1 = highest metric.

    Same rule as ``scripts/35_fig5a_metric_ranks.py`` (``method="min"``,
    NaN never ranks).
    """
    ranked = per_unit.copy()
    for m in list(metrics or METRICS):
        ranked[f"{m}_rank"] = ranked.groupby(["sample", "stratum"])[m].rank(
            method="min", ascending=False)
    return ranked


def mean_rank_table(ranked: pd.DataFrame, units: list[str], species: str,
                    dataset_group: str, metrics: list[str] | None = None,
                    strata: list[str] | None = None) -> pd.DataFrame:
    """Mean rank / value per tool x stratum x metric over ``units`` (drawn only)."""
    metrics = list(metrics or METRICS)
    strata = list(strata or ALL_STRATA)
    sub = ranked[(ranked["species"] == species) & (ranked["dataset_group"] == dataset_group)
                 & (ranked["sample"].isin(units))]
    rows = []
    n_total = len(units)
    for tool in TOOL_ORDER:
        t = sub[sub["tool"] == tool]
        if t.empty:
            continue
        for stratum in strata:
            ts = t[t["stratum"] == stratum]
            for m in metrics:
                rk = ts[f"{m}_rank"].dropna().astype(float)
                val = ts[m].dropna().astype(float)
                rows.append({
                    "species": species, "dataset_group": dataset_group, "tool": tool,
                    "stratum": stratum, "stratum_kind": stratum_kind(stratum),
                    "metric": m, "n_units": int(rk.size), "n_units_total": n_total,
                    "k_of_n": f"{int(rk.size)}/{n_total}",
                    "rank_mean": float(rk.mean()) if rk.size else np.nan,
                    "rank_min": float(rk.min()) if rk.size else np.nan,
                    "rank_max": float(rk.max()) if rk.size else np.nan,
                    "value_mean": float(val.mean()) if val.size else np.nan,
                    "value_min": float(val.min()) if val.size else np.nan,
                    "value_max": float(val.max()) if val.size else np.nan,
                })
    return pd.DataFrame(rows)


def _rho(a: pd.Series, b: pd.Series, method: str = "spearman") -> tuple[float, float]:
    both = pd.DataFrame({"a": a, "b": b}).dropna()
    if len(both) < 3:
        return float("nan"), float("nan")
    if method == "kendall":
        res = stats.kendalltau(both["a"], both["b"])
    else:
        res = stats.spearmanr(both["a"], both["b"])
    return float(res.statistic), float(res.pvalue)


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg step-up correction (same algorithm as NA2).

    ``adjusted_(i) = min_{j >= i} p_(j) * m / j`` on the ascending p-values, then
    mapped **back to the original order** (``adj[order] = ...``; forgetting the
    inverse permutation silently scrambles the column).
    """
    pv = np.asarray(pvals, dtype=float)
    mask = np.isfinite(pv)
    out = np.full(pv.shape, np.nan)
    if mask.sum() > 1:
        m = int(mask.sum())
        order = np.argsort(pv[mask])                     # index into the masked array
        ranked = pv[mask][order]
        adj = np.minimum.accumulate((ranked * m / np.arange(1, m + 1))[::-1])[::-1]
        back = np.empty(m, dtype=float)
        back[order] = np.clip(adj, 0, 1)
        out[mask] = back
    return out


def rank_pivot(ranked: pd.DataFrame, units: list[str], species: str,
               dataset_group: str, metric: str) -> pd.DataFrame:
    """tools x (stratum, sample) matrix of within-unit ranks for one metric."""
    sub = ranked[(ranked["species"] == species) & (ranked["dataset_group"] == dataset_group)
                 & (ranked["sample"].isin(units))]
    return sub.pivot_table(index="tool", columns=["stratum", "sample"],
                           values=f"{metric}_rank", aggfunc="mean")


def bootstrap_rho_ci(piv: pd.DataFrame, units: list[str], stratum: str,
                     b: int = BOOTSTRAP_B, seed: int = SEED) -> tuple[float, float]:
    """Unit-level bootstrap CI of the stratum-vs-pooled Spearman rho (MCC etc.).

    ``piv`` is the tool x (stratum, sample) rank matrix of one species x metric.
    Units are resampled with replacement; for each resample the mean rank per
    tool is recomputed on both strata and the correlation between the two
    resampled orderings is taken.  Fewer than 2 units -> NaN (nothing to
    resample; mouse has a single drawn unit by the Fig. 5A convention).
    """
    if len(units) < 2:
        return float("nan"), float("nan")
    if (POOLED, units[0]) not in piv.columns or (stratum, units[0]) not in piv.columns:
        return float("nan"), float("nan")
    ref = piv[POOLED].to_numpy(dtype=float)
    alt = piv[stratum].to_numpy(dtype=float)
    n = ref.shape[1]                      # units actually present in the pivot
    if n < 2:                             # partial run (fewer units than declared)
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    draws = np.empty(b, dtype=float)
    for i in range(b):
        pick = rng.integers(0, n, n)
        with warnings.catch_warnings():          # all-NaN rows are possible and
            warnings.simplefilter("ignore", RuntimeWarning)   # expected here
            a = np.nanmean(ref[:, pick], axis=1)
            c = np.nanmean(alt[:, pick], axis=1)
        draws[i], _ = _rho(pd.Series(a), pd.Series(c))
    draws = draws[np.isfinite(draws)]
    if draws.size == 0:
        return float("nan"), float("nan")
    return float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def stability_table(ranked: pd.DataFrame, panels, metrics=("mcc", "precision", "recall", "f1"),
                    strata=None, with_bootstrap: bool = True, family_key: str = "primary"
                    ) -> pd.DataFrame:
    """Rank stability of every stratum vs the pooled ranking, per species x metric.

    Rows: one per (species, dataset_group, metric, stratum).  ``p_value`` comes
    from the 13-tool Spearman test; ``p_value_fdr_bh`` is corrected across the
    whole returned family (documented in the README).
    """
    strata = list(strata or [s for s in ALL_STRATA if s != POOLED])
    metrics = list(metrics)
    rows = []
    for species, group, units in panels:
        mean = mean_rank_table(ranked, units, species, group, metrics=metrics,
                               strata=[POOLED] + [s for s in strata if s != POOLED])
        for m in metrics:
            piv = rank_pivot(ranked, units, species, group, m)
            wide = (mean[mean["metric"] == m]
                    .pivot_table(index="tool", columns="stratum", values="rank_mean"))
            for s in strata:
                if POOLED not in wide.columns or s not in wide.columns:
                    continue
                rho, p = _rho(wide[POOLED], wide[s])
                tau, _ = _rho(wide[POOLED], wide[s], method="kendall")
                both = wide[[POOLED, s]].dropna()
                shift = (both[s] - both[POOLED]).abs()
                top1_pooled = wide[POOLED].idxmin() if wide[POOLED].notna().any() else None
                top1_s = wide[s].idxmin() if wide[s].notna().any() else None
                top3_p = set(wide[POOLED].nsmallest(3).index)
                top3_s = set(wide[s].nsmallest(3).index)
                per_unit = []
                for u in units:
                    if (POOLED, u) not in piv.columns or (s, u) not in piv.columns:
                        continue
                    r, _p = _rho(piv[(POOLED, u)], piv[(s, u)])
                    per_unit.append(r)
                per_unit = np.asarray(per_unit, dtype=float)
                per_unit = per_unit[np.isfinite(per_unit)]
                ci_lo, ci_hi = (bootstrap_rho_ci(piv, units, s)
                                if with_bootstrap else (np.nan, np.nan))
                rows.append({
                    "family": family_key,
                    "species": species, "dataset_group": group, "metric": m,
                    "stratum": s, "stratum_kind": stratum_kind(s),
                    "n_units": len(units), "n_tools": int(both.shape[0]),
                    "spearman_rho": rho, "kendall_tau": tau, "p_value": p,
                    "rho_ci_lo": ci_lo, "rho_ci_hi": ci_hi,
                    "rho_per_unit_mean": float(per_unit.mean()) if per_unit.size else np.nan,
                    "rho_per_unit_min": float(per_unit.min()) if per_unit.size else np.nan,
                    "max_abs_rank_shift": float(shift.max()) if len(shift) else np.nan,
                    "mean_abs_rank_shift": float(shift.mean()) if len(shift) else np.nan,
                    "top1_tool_pooled": top1_pooled, "top1_tool_stratum": top1_s,
                    "top1_same": bool(top1_pooled is not None and top1_pooled == top1_s),
                    "top3_overlap": (len(top3_p & top3_s) / 3.0) if len(top3_p) == 3 else np.nan,
                })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["p_value_fdr_bh"] = _bh_fdr(out["p_value"].to_numpy(dtype=float))
    return out


def pairwise_stratum_table(ranked: pd.DataFrame, panels, metrics=("mcc",),
                           strata: list[str] | None = None) -> pd.DataFrame:
    """Spearman rho between the mean-rank orderings of every stratum pair."""
    strata = list(strata or ALL_STRATA)
    rows = []
    for species, group, units in panels:
        mean = mean_rank_table(ranked, units, species, group, strata=strata)
        for m in metrics:
            wide = (mean[mean["metric"] == m]
                    .pivot_table(index="tool", columns="stratum", values="rank_mean"))
            for i, a in enumerate(strata):
                for b_ in strata[i + 1:]:
                    if a not in wide.columns or b_ not in wide.columns:
                        continue
                    rho, p = _rho(wide[a], wide[b_])
                    rows.append({"species": species, "dataset_group": group, "metric": m,
                                 "stratum_a": a, "stratum_b": b_, "spearman_rho": rho,
                                 "p_value": p,
                                 "n_tools": int(wide[[a, b_]].dropna().shape[0])})
    return pd.DataFrame(rows)


__all__ = [
    "OUT", "TAB", "FIG", "LOG", "MOD_TYPE", "WINDOW", "MIN_COV", "POOLED",
    "MAIN_STRATA", "SENS_STRATUM", "ALL_STRATA", "BAND_EDGES", "DEFAULT_EDGES",
    "pooled_label", "bin_labels", "threshold_labels", "strata_for", "stratum_kind",
    "METRICS", "METRIC_LABEL", "GLOBAL_METRICS", "add_global_metrics",
    "TOOL_ORDER", "TOOL_DISPLAY", "PANELS", "EXTRA_GROUPS", "ALL_GROUPS",
    "PRIMARY_SPECIES", "BOOTSTRAP_B", "SEED",
    "BandedUniverse", "load_banded_universe", "pooled_codes", "universe_file",
    "callset_path", "units_for_groups", "load_calls", "calls_by_stratum",
    "strata_confusion", "metric_columns", "rank_within_unit", "mean_rank_table",
    "rank_pivot", "bootstrap_rho_ci", "stability_table", "pairwise_stratum_table",
    "read_table", "write_table", "setup_logger", "C",
]
