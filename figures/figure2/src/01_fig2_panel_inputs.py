#!/usr/bin/env python3
"""01 -- Figure 2 (revision): panel input tables A-E.

Why this script exists
----------------------
The published Figure 2 was drawn from ``output/<group>/Tools.txt`` -- a single
replicate in most groups (or an undocumented pile-up), with the panels computed
on a pooled union.  The revision rebuilds every panel from the per-sample layer
(``sites_v2/sites_clean/``) plus the replicate-aware evaluation tables, so that

* every number is traceable to one sample (= one independent sequencing unit);
* mouse is never pooled across its two studies (SRP166020 / SRP357195);
* the Arabidopsis WT x fip37 KD pairing is flagged as *cross-study*
  (``sample_registry.csv``: WT = SRP329449, KD = SRP363914);
* KO / KD are partial negatives and only HeLa IVT is a clean negative control.

Universe convention (locked)
----------------------------
Exactly the one used by the evaluation tables (``07_eval_controls.py``):
the **generic** ``universe/<sample>__universe.tsv[.gz]`` read through
``common.evaluation.load_universe`` -- coverage >= 10 **and** a
modification-compatible reference base -- grouped per chromosome as sorted
unique int64 arrays.  ``sites_clean/`` supplies the call sets (position level,
one row per site x overlapping gene collapsed to unique ``(chrom, pos)``).

Panels (drawing is ``03_fig2_figure.py``; GO enrichment is
``02_fig2_go_enrichment.R``)

A  detection counts per independent unit (WT x control, log-log) and the
   control/WT ratio on the *shared testable universe* -- the same definition as
   ``ko_kd_metrics.tsv`` (HeLa is computed here because that table has no human
   rows);
B  inter-tool support: cumulative fraction of a unit's sites supported by >= k
   of the 13 m6A tools, one curve per independent unit;
C  high-confidence sites (support >= 5 inside a unit, then the majority
   consensus over the species' units) -> host genes -> foreground/background
   gene lists for the species-native GO:BP enrichment;
D  wild-type vs control overlap (Jaccard over every unit pair) and the explicit
   control-side false-positive burden (control calls not within 2 bp of GLORI,
   per 10^6 candidate positions of that control's testable universe);
E  replicate consistency per tool: mean pairwise Jaccard vs the pooled
   ("global") Jaccard |intersection| / |union| over the same units.

Usage
-----
conda run -n benchmark-revision --no-capture-output \
  python 01_fig2_panel_inputs.py            # ~2 min (GTF parse cached)

Outputs go to ``../tables/``; ``--out`` overrides.
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
import argparse
import subprocess
import sys
import time
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
REV = HERE.parents[1]
CODE_ROOT = Path(str(_RB / "src/sites_v2"))
sys.path.insert(0, str(CODE_ROOT))

from common.config import (  # noqa: E402
    BOOTSTRAP_B, GLORI, REF_ROOT, SEED, SITES_ROOT, UNIVERSE_ROOT,
)
from common.consensus import quorum  # noqa: E402
from common.evaluation import load_universe  # noqa: E402
from common.match import fix_chromosome  # noqa: E402

PLATFORM = "RNA002"
MOD = "m6A"
MIN_COV = 10            # candidate-universe coverage threshold (main convention)
SUPPORT_HIGH = 5        # published "high-confidence" support level
WINDOW = 2              # primary GLORI matching window

SITES_CLEAN = (_RB / "data/sites_clean")
TABLE_DIR = (_RB / "data/evaluation/tables")

#: species block -> dataset groups; ``pairing`` = how a WT unit is matched to a
#: control unit: by replicate tag, or inside a study (mouse: never across).
SPECIES_SPEC: dict[str, dict[str, str]] = {
    "Arabidopsis": dict(wt="Arabidopsis_WT", ctrl="Arabidopsis_KD",
                        ctrl_label="fip37 KD", pairing="replicate_tag"),
    "Mouse": dict(wt="Mouse_WT", ctrl="Mouse_KO",
                  ctrl_label="Mettl3 KO", pairing="study"),
    "Human": dict(wt="HeLa_WT", ctrl="HeLa_IVT",
                  ctrl_label="IVT", pairing="replicate_tag"),
}
BLOCK_ORDER = ["Arabidopsis", "Mouse", "Human"]
TOOLS = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
         "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
         "NanoSPA_m6A", "xPore", "yanocomp"]

#: gene-level GTF per species (Ensembl/TAIR annotation; project rule: no GENCODE)
GTF_PATH = {
    "Arabidopsis": REF_ROOT / "arabidopsis" / "Arabidopsis_thaliana.TAIR10.61.ensembl.gtf",
    "Mouse": REF_ROOT / "GRCm39" / "ensembl" / "Mus_musculus.GRCm39.114.gtf",
    "Human": REF_ROOT / "GRCh38p14" / "ensembl112" / "Homo_sapiens.GRCh38.112.chr.gtf",
}

#: {chrom: sorted unique positions} -- the shape every helper works on
Pos = dict[str, np.ndarray]
#: parsed universes, keyed by (species, sample, mod) -- see :func:`universe`
_UNI_CACHE: dict[tuple[str, str, str], Pos] = {}
log_lines: list[str] = []


def log(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    log_lines.append(line)


# --------------------------------------------------------------------------- #
# generic position-set helpers
# --------------------------------------------------------------------------- #
def n_pos(d: Pos) -> int:
    return int(sum(v.size for v in d.values()))


def intersect(a: Pos, b: Pos) -> Pos:
    return {c: np.intersect1d(a[c], b[c]) for c in set(a) & set(b)
            if np.intersect1d(a[c], b[c]).size}


def n_intersect(a: Pos, b: Pos) -> int:
    return int(sum(np.intersect1d(a[c], b[c]).size for c in set(a) & set(b)))


def jaccard(a: Pos, b: Pos) -> float:
    inter = n_intersect(a, b)
    union = n_pos(a) + n_pos(b) - inter
    return inter / union if union else np.nan


def to_tuples(d: Pos) -> set[tuple[str, int]]:
    return {(c, int(p)) for c, v in d.items() for p in v}


def from_tuples(s: set[tuple[str, int]]) -> Pos:
    out: dict[str, list[int]] = {}
    for c, p in s:
        out.setdefault(c, []).append(p)
    return {c: np.unique(np.asarray(v, dtype=np.int64)) for c, v in out.items()}


def boot_ci(values, lo=2.5, hi=97.5) -> tuple[float, float]:
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    if v.size < 2:
        return (np.nan, np.nan)
    rng = np.random.default_rng(SEED)
    idx = rng.integers(0, v.size, size=(BOOTSTRAP_B, v.size))
    return float(np.percentile(v[idx].mean(axis=1), lo)), \
        float(np.percentile(v[idx].mean(axis=1), hi))


# --------------------------------------------------------------------------- #
# loading
# --------------------------------------------------------------------------- #
def units_of(group: str) -> pd.DataFrame:
    """Independent sequencing units of a dataset group (nested subsets removed)."""
    reg = pd.read_csv((_RB / "metadata/sample_registry.csv"), sep="\t",
                      dtype=str)
    reg = reg[(reg.platform == PLATFORM) & (reg.dataset_group == group)
              & (reg.independence_class != "nested_subset")]
    return reg.drop_duplicates("sequencing_unit").sort_values(
        "replicate_tag").reset_index(drop=True)


def calls(species: str, group: str, tool: str, sample: str) -> Pos:
    """Position-level call set: unique ``(chrom, pos)`` per chromosome.

    A coordinate appears once per overlapping gene on two transcripts, so the
    site-level semantics is the unique key (same as the ``np.unique`` used by
    the confusion tables).
    """
    p = SITES_CLEAN / PLATFORM / species / group / MOD / tool / f"{sample}.tsv"
    if not p.exists():
        return {}
    d = pd.read_csv(p, sep="\t", usecols=["chrom", "start"], dtype={"chrom": str})
    if d.empty:
        return {}
    d = d.assign(chrom=d["chrom"].map(fix_chromosome),
                 pos=pd.to_numeric(d["start"], errors="coerce")).dropna(subset=["pos"])
    return {c: np.unique(g["pos"].to_numpy(dtype=np.int64))
            for c, g in d.groupby("chrom", sort=False)}


def universe(sample: str, species: str, mod: str = MOD) -> Pos:
    """Generic universe file through the evaluation-table loader (authoritative).

    Cached per sample: the files have 10^7 rows and are re-read by several
    panels, so each one is parsed once per run (~40 MB as int64 per sample).
    """
    key = (species, sample, mod)
    if key not in _UNI_CACHE:
        plain = UNIVERSE_ROOT / PLATFORM / species / f"{sample}__universe.tsv"
        gz = plain.with_suffix(".tsv.gz")
        path = plain if plain.exists() else gz
        if not path.exists():
            raise FileNotFoundError(path)
        _UNI_CACHE[key] = load_universe(path, mod, MIN_COV)
    return _UNI_CACHE[key]


def reference(species: str) -> Pos:
    """GLORI positions per chromosome (sorted), as used by the evaluation."""
    g = pd.read_csv(GLORI[species], sep="\t", header=None, usecols=[0, 1])
    pos: dict[str, list[int]] = {}
    for c, p in zip(g[0].astype(str).map(fix_chromosome), g[1].astype(int)):
        pos.setdefault(c, []).append(int(p))
    return {c: np.sort(np.asarray(v, dtype=np.int64)) for c, v in pos.items()}


def ref_in_universe(ref: Pos, uni: Pos) -> Pos:
    """Reference positions that lie in the universe (testable reference sites)."""
    out: Pos = {}
    for c, arr in ref.items():
        u = uni.get(c)
        if u is None or u.size == 0:
            continue
        i = np.searchsorted(u, arr)
        i = np.clip(i, 0, u.size - 1)
        keep = u[i] == arr
        if keep.any():
            out[c] = arr[keep]
    return out


def score_vs_ref(call: Pos, uni: Pos, ref: Pos, window: int = WINDOW
                 ) -> tuple[int, int]:
    """(calls inside the universe, calls within ``window`` bp of a reference site)."""
    n_calls = tp = 0
    for c, pos in call.items():
        u = uni.get(c)
        if u is None or u.size == 0:
            continue
        i = np.searchsorted(u, pos)
        i = np.clip(i, 0, u.size - 1)
        inside = u[i] == pos
        pos = pos[inside]
        if pos.size == 0:
            continue
        n_calls += int(pos.size)
        r = ref.get(c)
        if r is None or r.size == 0:
            continue
        j = np.searchsorted(r, pos)
        lo = r[np.clip(j - 1, 0, r.size - 1)]
        eq = r[np.clip(j, 0, r.size - 1)]
        hi = r[np.clip(j + 1, 0, r.size - 1)]
        dist = np.minimum(np.minimum(np.abs(pos - lo), np.abs(pos - eq)),
                          np.abs(pos - hi))
        tp += int((dist <= window).sum())
    return n_calls, tp


# --------------------------------------------------------------------------- #
# panel A
# --------------------------------------------------------------------------- #
def panel_a(out: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, ratio_rows = [], []
    ko_kd = pd.read_csv(TABLE_DIR / "ko_kd_metrics.tsv", sep="\t")

    for species in BLOCK_ORDER:
        spec = SPECIES_SPEC[species]
        wt_units, ctrl_units = units_of(spec["wt"]), units_of(spec["ctrl"])
        study_of = dict(zip(wt_units["sample"], wt_units["study"]))
        study_of.update(dict(zip(ctrl_units["sample"], ctrl_units["study"])))
        uni = {s: universe(s, species) for s in
               list(wt_units["sample"]) + list(ctrl_units["sample"])}
        sets: dict[tuple[str, str, str], Pos] = {}
        for _, u in wt_units.iterrows():
            for tool in TOOLS:
                sets[("wt", u["sample"], tool)] = calls(species, spec["wt"], tool, u["sample"])
        for _, u in ctrl_units.iterrows():
            for tool in TOOLS:
                sets[("ctrl", u["sample"], tool)] = calls(species, spec["ctrl"], tool, u["sample"])

        for tag, units, group, cond in (("wt", wt_units, spec["wt"], "WT"),
                                        ("ctrl", ctrl_units, spec["ctrl"], spec["ctrl_label"])):
            for _, u in units.iterrows():
                for tool in TOOLS:
                    s = sets[(tag, u["sample"], tool)]
                    rows.append(dict(species=species, group=group, condition=cond,
                                     tool=tool, sample=u["sample"],
                                     replicate_tag=u["replicate_tag"],
                                     unit=u["sequencing_unit"], study=u["study"],
                                     n_sites=n_pos(s),
                                     n_sites_in_universe=n_intersect(s, uni[u["sample"]])))

        if spec["pairing"] == "replicate_tag":
            pairs = [(w["sample"], c["sample"], str(w["replicate_tag"]))
                     for _, w in wt_units.iterrows()
                     for _, c in ctrl_units.iterrows()
                     if w["replicate_tag"] == c["replicate_tag"]]
        else:  # study-internal pairing, never across studies
            pairs = [(wr["sample"], cr["sample"], str(st))
                     for st in sorted(set(wt_units["study"]))
                     for _, wr in wt_units[wt_units.study == st].iterrows()
                     for _, cr in ctrl_units[ctrl_units.study == st].iterrows()]

        for ws, cs, tag in pairs:
            common = intersect(uni[ws], uni[cs])
            for tool in TOOLS:
                n_wt = n_intersect(sets[("wt", ws, tool)], common)
                n_ctrl = n_intersect(sets[("ctrl", cs, tool)], common)
                ratio_rows.append(dict(
                    species=species, tool=tool, wt_sample=ws, ctrl_sample=cs, pair=tag,
                    same_study=int(study_of.get(ws) == study_of.get(cs)),
                    n_common_universe=n_pos(common), n_wt_common=n_wt,
                    n_ctrl_common=n_ctrl,
                    ctrl_wt_ratio=(n_ctrl / n_wt) if n_wt else np.nan))

    counts = pd.DataFrame(rows)
    counts.to_csv(out / "fig2a_counts_by_unit.tsv", sep="\t", index=False)

    ratio = pd.DataFrame(ratio_rows)
    ev = ko_kd[["species", "tool", "sample_wt", "sample_ko", "same_study",
                "n_common_universe", "n_wt_common", "n_ko_common", "ko_wt_ratio"]]
    ratio = ratio.merge(
        ev.rename(columns={"sample_wt": "wt_sample", "sample_ko": "ctrl_sample",
                           "same_study": "eval_same_study",
                           "n_common_universe": "eval_n_common",
                           "n_wt_common": "eval_n_wt_common",
                           "n_ko_common": "eval_n_ctrl_common",
                           "ko_wt_ratio": "eval_ko_wt_ratio"}),
        on=["species", "tool", "wt_sample", "ctrl_sample"], how="left")
    ratio["eval_delta"] = ratio["ctrl_wt_ratio"] - ratio["eval_ko_wt_ratio"]
    ratio["eval_int_count_delta"] = (
        (ratio["n_wt_common"] - ratio["eval_n_wt_common"]).abs()
        + (ratio["n_ctrl_common"] - ratio["eval_n_ctrl_common"]).abs())
    ratio.to_csv(out / "fig2a_testable_ratio.tsv", sep="\t", index=False)

    summ = (ratio.groupby(["species", "tool"])
            .agg(n_pairs=("ctrl_wt_ratio", "size"),
                 ratio_mean=("ctrl_wt_ratio", "mean"),
                 ratio_sd=("ctrl_wt_ratio", "std"),
                 ratio_min=("ctrl_wt_ratio", "min"),
                 ratio_max=("ctrl_wt_ratio", "max"),
                 ratio_ci_lo=("ctrl_wt_ratio", lambda v: boot_ci(v)[0]),
                 ratio_ci_hi=("ctrl_wt_ratio", lambda v: boot_ci(v)[1]),
                 eval_max_abs_count_delta=("eval_int_count_delta", "max"),
                 eval_max_abs_ratio_delta=("eval_delta", lambda v: np.nanmax(np.abs(v))
                                           if v.notna().any() else np.nan))
            .reset_index())
    summ.to_csv(out / "fig2a_ratio_summary.tsv", sep="\t", index=False)

    n_cmp = int(ratio["eval_ko_wt_ratio"].notna().sum())
    worst = float(np.nanmax(ratio["eval_delta"])) if n_cmp else float("nan")
    worst_cnt = float(np.nanmax(ratio["eval_int_count_delta"])) if n_cmp else float("nan")
    log(f"A: {len(counts)} count rows, {len(ratio)} pair rows; {n_cmp} cross-checked "
        f"against ko_kd_metrics (max |ratio delta| = {worst:.2e}, "
        f"max |count delta| = {worst_cnt:.0f})")
    return counts, ratio


# --------------------------------------------------------------------------- #
# panel B
# --------------------------------------------------------------------------- #
def panel_b(out: Path) -> pd.DataFrame:
    rows = []
    for species in BLOCK_ORDER:
        spec = SPECIES_SPEC[species]
        for _, u in units_of(spec["wt"]).iterrows():
            avail, support = [], Counter()
            for t in TOOLS:
                s = calls(species, spec["wt"], t, u["sample"])
                if n_pos(s):
                    avail.append(t)
                    support.update(to_tuples(s))
            if not support:
                continue
            total = len(support)
            for k in range(1, len(avail) + 1):
                n_ge = sum(1 for v in support.values() if v >= k)
                rows.append(dict(species=species, sample=u["sample"],
                                 replicate_tag=u["replicate_tag"],
                                 unit=u["sequencing_unit"], study=u["study"], k=k,
                                 n_tools_available=len(avail), n_sites_total=total,
                                 n_sites_support_ge_k=n_ge,
                                 frac_support_ge_k=n_ge / total))
    d = pd.DataFrame(rows)
    d.to_csv(out / "fig2b_support_curves.tsv", sep="\t", index=False)
    summ = (d.groupby(["species", "k"])
            .agg(n_units=("frac_support_ge_k", "size"),
                 frac_mean=("frac_support_ge_k", "mean"),
                 frac_min=("frac_support_ge_k", "min"),
                 frac_max=("frac_support_ge_k", "max"),
                 n_sites_mean=("n_sites_support_ge_k", "mean"))
            .reset_index())
    summ.to_csv(out / "fig2b_support_summary.tsv", sep="\t", index=False)
    log(f"B: {d['species'].nunique()} species x {d['sample'].nunique()} units, "
        f"k up to {int(d['k'].max()) if len(d) else 0}")
    return d


# --------------------------------------------------------------------------- #
# gene annotation (cached) + overlap
# --------------------------------------------------------------------------- #
def gene_table(species: str, out: Path) -> pd.DataFrame:
    cache = out / "_refs" / f"{species}.genes.tsv.gz"
    gtf = GTF_PATH[species]
    if not gtf.exists():
        raise FileNotFoundError(gtf)
    if cache.exists() and cache.stat().st_mtime > gtf.stat().st_mtime:
        return pd.read_csv(cache, sep="\t")
    log(f"  parsing gene features from {gtf.name} (cached afterwards)")
    txt = subprocess.run(["awk", '-F\t', '$3=="gene"', str(gtf)], check=True,
                         capture_output=True, text=True).stdout.splitlines()
    recs = []
    for ln in txt:
        f = ln.split("\t")
        if len(f) < 9:
            continue
        attrs = {}
        for kv in f[8].strip().rstrip(";").split("; "):
            if " " in kv:
                k, v = kv.split(" ", 1)
                attrs[k] = v.strip().strip('";')
        recs.append((fix_chromosome(f[0]), int(f[3]) - 1, int(f[4]),
                     attrs.get("gene_id", ""), attrs.get("gene_name", ""),
                     attrs.get("gene_biotype", attrs.get("gene_type", ""))))
    d = pd.DataFrame(recs, columns=["chrom", "start", "end", "gene_id", "gene_name",
                                    "gene_biotype"]).drop_duplicates("gene_id")
    cache.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(cache, sep="\t", index=False)
    log(f"  {species}: {len(d):,} genes cached -> {cache.name}")
    return d


def gene_index(genes: pd.DataFrame) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    out: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for c, g in genes.groupby("chrom", sort=False):
        out[c] = (g["start"].to_numpy(), g["end"].to_numpy(), g["gene_id"].to_numpy())
    return out


def genes_of(pos: Pos,
             index: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]) -> set[str]:
    """Genes whose span ``[start, end)`` contains at least one queried position.

    Vectorised per chromosome: for a gene, the positions inside it are the
    contiguous slice ``[searchsorted(pos, start), searchsorted(pos, end))`` of
    the sorted position array, so the gene is a hit iff that slice is non-empty.
    """
    hit: set[str] = set()
    for c, p in pos.items():
        g = index.get(c)
        if g is None or p.size == 0:
            continue
        starts, ends, gids = g
        p = np.sort(p)
        lo = np.searchsorted(p, starts, side="left")
        hi = np.searchsorted(p, ends, side="left")
        hit |= set(gids[hi > lo])
    return hit


def panel_c(out: Path) -> pd.DataFrame:
    site_rows, contrib_rows = [], []
    for species in BLOCK_ORDER:
        spec = SPECIES_SPEC[species]
        units = units_of(spec["wt"])
        index = gene_index(gene_table(species, out))
        per_unit: dict[str, set] = {}
        per_unit_support: dict[str, Counter] = {}
        per_tool_sites: dict[str, list[set]] = {t: [] for t in TOOLS}
        for _, u in units.iterrows():
            sup: Counter = Counter()
            for t in TOOLS:
                s = to_tuples(calls(species, spec["wt"], t, u["sample"]))
                per_tool_sites[t].append(s)
                sup.update(s)
            per_unit_support[u["sample"]] = sup
            per_unit[u["sample"]] = {k for k, v in sup.items() if v >= SUPPORT_HIGH}
        thr = quorum(len(per_unit))
        cnt: Counter = Counter()
        for s in per_unit.values():
            cnt.update(s)
        consensus = {k for k, v in cnt.items() if v >= thr}

        # which tools keep their calls in the majority of the species' units
        for t in TOOLS:
            keep: Counter = Counter()
            for s in per_tool_sites[t]:
                keep.update(s)
            tool_cons = {k for k, v in keep.items() if v >= thr}
            contrib_rows.append(dict(
                species=species, tool=t, n_units_with_calls=len(per_tool_sites[t]),
                n_calls_majority=len(tool_cons),
                n_calls_in_high_conf=len(tool_cons & consensus)))

        common: Pos | None = None
        for s in units["sample"]:
            u = universe(s, species)
            common = u if common is None else intersect(common, u)
        bg_genes = genes_of(common or {}, index)
        fg_genes = genes_of(from_tuples(consensus), index)

        for _, u in units.iterrows():
            site_rows.append(dict(
                species=species, sample=u["sample"], replicate_tag=u["replicate_tag"],
                unit=u["sequencing_unit"], study=u["study"], quorum=thr,
                n_sites_total=len(per_unit_support[u["sample"]]),
                n_sites_support_ge5=len(per_unit[u["sample"]])))
        site_rows.append(dict(species=species, sample="<majority consensus>",
                              replicate_tag="", unit="", study="", quorum=thr,
                              n_sites_total=len(consensus),
                              n_sites_support_ge5=len(consensus)))
        pd.DataFrame([dict(species=species, chrom=c, pos=p,
                           n_units_supporting=int(cnt[(c, p)]), consensus_quorum=thr)
                      for c, p in sorted(consensus)]
                     ).to_csv(out / f"fig2c_high_conf_sites_{species}.tsv", sep="\t",
                              index=False)
        pd.DataFrame([dict(species=species, gene_id=g, role="foreground")
                      for g in sorted(fg_genes)]
                     + [dict(species=species, gene_id=g, role="background")
                        for g in sorted(bg_genes)]
                     ).to_csv(out / f"fig2c_genes_{species}.tsv", sep="\t", index=False)
        log(f"C: {species}: majority consensus (quorum {thr}/{len(per_unit)}) = "
            f"{len(consensus):,} high-confidence sites -> {len(fg_genes):,} foreground "
            f"genes; background {len(bg_genes):,} genes")

    sites_tab = pd.DataFrame(site_rows)
    sites_tab.to_csv(out / "fig2c_high_conf_by_unit.tsv", sep="\t", index=False)
    contrib = pd.DataFrame(contrib_rows)
    contrib.to_csv(out / "fig2c_tool_contribution.tsv", sep="\t", index=False)
    return contrib


# --------------------------------------------------------------------------- #
# panel D
# --------------------------------------------------------------------------- #
def panel_d(out: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    jac_rows, fp_rows = [], []
    for species in BLOCK_ORDER:
        spec = SPECIES_SPEC[species]
        wt_units, ctrl_units = units_of(spec["wt"]), units_of(spec["ctrl"])
        ref = reference(species)
        sets: dict[tuple[str, str, str], Pos] = {}
        for _, u in wt_units.iterrows():
            for t in TOOLS:
                sets[("wt", u["sample"], t)] = calls(species, spec["wt"], t, u["sample"])
        for _, u in ctrl_units.iterrows():
            for t in TOOLS:
                sets[("ctrl", u["sample"], t)] = calls(species, spec["ctrl"], t, u["sample"])

        for _, w in wt_units.iterrows():
            for _, c in ctrl_units.iterrows():
                for t in TOOLS:
                    a, b = sets[("wt", w["sample"], t)], sets[("ctrl", c["sample"], t)]
                    if not a and not b:
                        continue
                    jac_rows.append(dict(species=species, tool=t,
                                         wt_sample=w["sample"], ctrl_sample=c["sample"],
                                         wt_unit=w["sequencing_unit"],
                                         ctrl_unit=c["sequencing_unit"],
                                         same_study=int(w["study"] == c["study"]),
                                         n_wt=n_pos(a), n_ctrl=n_pos(b),
                                         n_shared=n_intersect(a, b), jaccard=jaccard(a, b)))

        for _, c in ctrl_units.iterrows():
            u = universe(c["sample"], species)
            ref_u = ref_in_universe(ref, u)
            for t in TOOLS:
                s = sets[("ctrl", c["sample"], t)]
                if not n_pos(s):
                    continue
                n_calls, tp = score_vs_ref(s, u, ref_u)
                fp = n_calls - tp
                n_u = n_pos(u)
                fp_rows.append(dict(
                    species=species, tool=t, ctrl_sample=c["sample"],
                    ctrl_unit=c["sequencing_unit"],
                    n_calls_total=n_pos(s), n_calls_in_universe=n_calls,
                    tp_glori_2bp=tp, fp_glori_unsupported=fp, n_universe=n_u,
                    fp_per_1e6_candidates=fp / n_u * 1e6,
                    fp_fraction_of_universe=fp / n_u,
                    # same convention as controls_ivt_fpr.tsv (every call counts)
                    fp_per_1e6_all_calls=n_pos(s) / n_u * 1e6))

    jac = pd.DataFrame(jac_rows)
    jac.to_csv(out / "fig2d_wt_ctrl_jaccard.tsv", sep="\t", index=False)
    summ = (jac.groupby(["species", "tool"])
            .agg(n_pairs=("jaccard", "size"),
                 jaccard_mean=("jaccard", "mean"), jaccard_sd=("jaccard", "std"),
                 jaccard_min=("jaccard", "min"), jaccard_max=("jaccard", "max"),
                 jaccard_ci_lo=("jaccard", lambda v: boot_ci(v)[0]),
                 jaccard_ci_hi=("jaccard", lambda v: boot_ci(v)[1]))
            .reset_index())
    summ.to_csv(out / "fig2d_jaccard_summary.tsv", sep="\t", index=False)

    fp = pd.DataFrame(fp_rows)
    fp.to_csv(out / "fig2d_ctrl_fp_by_unit.tsv", sep="\t", index=False)
    fp_summ = (fp.groupby(["species", "tool"])
               .agg(n_units=("fp_per_1e6_candidates", "size"),
                    fp_per_1e6_mean=("fp_per_1e6_candidates", "mean"),
                    fp_per_1e6_min=("fp_per_1e6_candidates", "min"),
                    fp_per_1e6_max=("fp_per_1e6_candidates", "max"),
                    tp_glori_2bp_sum=("tp_glori_2bp", "sum"))
               .reset_index())
    fp_summ.to_csv(out / "fig2d_fp_summary.tsv", sep="\t", index=False)

    # cross-check against the published-metric table for the samples it covers:
    # the *universe* is identical (same loader), while the call counts can differ
    # because that table counts raw callset rows and this one counts unique
    # positions (a site annotated to two overlapping transcripts appears twice).
    pub = pd.read_csv(TABLE_DIR / "controls_ivt_fpr.tsv", sep="\t")
    pub = pub[(pub.platform == PLATFORM) & (pub.mod_type == MOD)][
        ["sample", "tool", "n_calls", "n_calls_in_universe", "n_universe",
         "fp_per_1e6_candidates"]]
    chk = fp.merge(pub.rename(columns={"sample": "ctrl_sample",
                                       "n_calls": "eval_n_calls_raw_rows",
                                       "n_calls_in_universe": "eval_n_calls_in_universe_rows",
                                       "n_universe": "eval_n_universe",
                                       "fp_per_1e6_candidates": "eval_fp_per_1e6_all_rows"}),
                   on=["ctrl_sample", "tool"], how="left")
    chk["unique_over_row_ratio"] = chk["n_calls_total"] / chk["eval_n_calls_raw_rows"]
    chk.to_csv(out / "fig2d_ivt_eval_crosscheck.tsv", sep="\t", index=False)
    both = chk.dropna(subset=["eval_fp_per_1e6_all_rows"])
    if len(both):
        same_u = int((both["n_universe"] == both["eval_n_universe"]).sum())
        log(f"D: IVT cross-check vs controls_ivt_fpr.tsv on {len(both)} rows: "
            f"universe identical on {same_u}/{len(both)}; unique-position calls are "
            f"{np.median(both['unique_over_row_ratio']):.3f}x the raw row count "
            "(multi-transcript duplicate rows)")
    log(f"D: {len(jac)} unit pairs, {len(fp)} control-unit x tool FP rows")
    return jac, fp


# --------------------------------------------------------------------------- #
# panel E
# --------------------------------------------------------------------------- #
def panel_e(out: Path) -> pd.DataFrame:
    rows = []
    for species in BLOCK_ORDER:
        spec = SPECIES_SPEC[species]
        units = units_of(spec["wt"])
        for t in TOOLS:
            sets = [calls(species, spec["wt"], t, s) for s in units["sample"]]
            sets = [s for s in sets if n_pos(s)]
            if len(sets) < 2:
                continue
            pw = [jaccard(a, b) for a, b in combinations(sets, 2)]
            union: Pos = {}
            for s in sets:
                for c, v in s.items():
                    union[c] = np.union1d(union.get(c, np.zeros(0, np.int64)), v)
            inter: Pos | None = None
            for s in sets:
                inter = s if inter is None else intersect(inter, s)
            glob = n_pos(inter or {}) / n_pos(union) if n_pos(union) else np.nan
            for a, b in combinations(range(len(sets)), 2):
                rows.append(dict(species=species, tool=t, metric="pairwise",
                                 pair=f"{a + 1}-{b + 1}", value=jaccard(sets[a], sets[b])))
            rows.append(dict(species=species, tool=t, metric="global",
                             pair="intersection/union", value=glob))
            rows.append(dict(species=species, tool=t, metric="mean_pairwise",
                             pair=f"n_pairs={len(pw)}", value=float(np.mean(pw))))
    d = pd.DataFrame(rows)
    d.to_csv(out / "fig2e_replicate_consistency.tsv", sep="\t", index=False)
    log(f"E: {d['species'].nunique()} species x {d['tool'].nunique()} tools")
    return d


# --------------------------------------------------------------------------- #
# tool x panel inclusion matrix
# --------------------------------------------------------------------------- #
def tool_panel_matrix(out: Path, counts: pd.DataFrame, ratio: pd.DataFrame,
                      sup: pd.DataFrame, contrib: pd.DataFrame, jac: pd.DataFrame,
                      fp: pd.DataFrame, rep: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for species in BLOCK_ORDER:
        spec = SPECIES_SPEC[species]
        n_ctrl_units = len(units_of(spec["ctrl"]))
        for t in TOOLS:
            cnt = counts[(counts.species == species) & (counts.tool == t)]
            n_wt_units = int(cnt[cnt.condition == "WT"].shape[0])
            ct_called = cnt[(cnt.condition != "WT") & (cnt.n_sites > 0)]
            wt_called = cnt[(cnt.condition == "WT") & (cnt.n_sites > 0)]
            j = jac[(jac.species == species) & (jac.tool == t)]
            f = fp[(fp.species == species) & (fp.tool == t)]
            e = rep[(rep.species == species) & (rep.tool == t) & (rep.metric == "pairwise")]
            hi = contrib[(contrib.species == species) & (contrib.tool == t)]
            reasons = []
            if n_wt_units == 0:
                reasons.append("no WT unit")
            if wt_called.empty:
                reasons.append("no call in WT -> absent from A/B/C/E")
            if ct_called.empty:
                reasons.append(f"no call in {spec['ctrl_label']} -> no D pair")
            if j.empty:
                reasons.append("no WT x control pair with calls")
            if f.empty:
                reasons.append("control has no calls -> no control-side FP rate")
            if e.empty:
                reasons.append("< 2 WT units with calls -> no replicate consistency")
            n_hi = int(hi["n_calls_in_high_conf"].iloc[0]) if len(hi) else 0
            if n_hi == 0:
                reasons.append("no call inside the majority high-confidence set -> not in C")
            # which individual units contributed nothing (tool x sample coverage,
            # needed to justify inclusion/exclusion per panel -- R1-6)
            silent = cnt[cnt.n_sites == 0]
            units_without_calls = "; ".join(
                f"{r.condition}:{r.sample}" for r in silent.itertuples())
            n_silent = len(silent)
            if n_silent:
                reasons.append(f"{n_silent} unit(s) with zero calls")
            rows.append(dict(species=species, tool=t,
                             panel_A_counts=int(not wt_called.empty or not ct_called.empty),
                             panel_A_ratio=int(not ratio[(ratio.species == species)
                                                         & (ratio.tool == t)].empty),
                             panel_B_support=int(not sup[(sup.species == species)].empty),
                             panel_C_high_conf=int(n_hi > 0),
                             panel_D_jaccard=int(not j.empty),
                             panel_D_fp_rate=int(not f.empty),
                             panel_E_consistency=int(not e.empty),
                             n_wt_units=n_wt_units, n_ctrl_units=n_ctrl_units,
                             n_high_conf_sites=n_hi, n_units_without_calls=n_silent,
                             units_without_calls=units_without_calls,
                             inclusion_note="; ".join(reasons)))
    d = pd.DataFrame(rows)
    d.to_csv(out / "fig2_tool_panel_matrix.tsv", sep="\t", index=False)
    log(f"matrix: {len(d)} species x tool rows, "
        f"{int((d.inclusion_note != '').sum())} with an exclusion note")
    return d


LEGACY = [
    ("B", "fraction of sites detected by exactly ONE tool", "76.17%",
     "frac_single_tool"),
    ("B", "mean number of distinct sites per independent unit at k = 1",
     "277,819 (pooled across replicates)", "n_sites_mean_k1"),
    ("B", "fraction of sites supported by >= 4 tools", "< 1%", "frac_ge4"),
    ("A", "human control/WT ratio, ELIGOS2_diff (pooled files)", "0.062",
     "ratio_Human_ELIGOS2_diff"),
    ("D", "mean WT-vs-control Jaccard, Nanom6A (human)", "0.448", "jac_Human_Nanom6A"),
    ("D", "mean WT-vs-control Jaccard, MINES (human)", "0.424", "jac_Human_MINES"),
    ("D", "mean WT-vs-control Jaccard, CHEUI_m6A (human)", "0.015", "jac_Human_CHEUI_m6A"),
    ("D", "mean WT-vs-control Jaccard, Yanocomp (human)", "0.023", "jac_Human_yanocomp"),
    ("D", "mean WT-vs-control Jaccard, Nanocompore (human)", "0.029",
     "jac_Human_Nanocompore"),
]


def reconciliation(out: Path, sup: pd.DataFrame, ratio: pd.DataFrame,
                   jac: pd.DataFrame) -> pd.DataFrame:
    new: dict[str, float] = {}
    # k = 1 is 100 % by construction; "detected by a single tool" is 1 - P(>= 2)
    new["frac_single_tool"] = 1.0 - float(sup[sup.k == 2]["frac_support_ge_k"].mean())
    new["frac_ge4"] = float(sup[sup.k == 4]["frac_support_ge_k"].mean())
    new["n_sites_mean_k1"] = float(sup[sup.k == 1]["n_sites_support_ge_k"].mean())

    def pair_mean(df: pd.DataFrame, tool: str, col: str) -> float:
        d = df[(df.species == "Human") & (df.tool == tool)]
        return float(d[col].mean()) if len(d) else float("nan")

    new["ratio_Human_ELIGOS2_diff"] = pair_mean(ratio, "ELIGOS2_diff", "ctrl_wt_ratio")
    for tool, key in (("Nanom6A", "jac_Human_Nanom6A"), ("MINES", "jac_Human_MINES"),
                      ("CHEUI_m6A", "jac_Human_CHEUI_m6A"),
                      ("yanocomp", "jac_Human_yanocomp"),
                      ("Nanocompore", "jac_Human_Nanocompore")):
        new[key] = pair_mean(jac, tool, "jaccard")

    rows = []
    for panel, what, legacy, key in LEGACY:
        val = new.get(key, float("nan"))
        rows.append(dict(panel=panel, metric=what, legacy_published=legacy,
                         revision_value="" if not np.isfinite(val) else f"{val:.6g}",
                         note="legacy = pooled / single-replicate aggregate; revision = "
                              "per independent unit on the shared testable universe"))
    d = pd.DataFrame(rows)
    d.to_csv(out / "legacy_vs_revision_fig2.tsv", sep="\t", index=False)
    return d


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REV / "tables"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    log(f"platform {PLATFORM}, tools = {len(TOOLS)}, min coverage = {MIN_COV}, "
        f"support-high = {SUPPORT_HIGH}, window = {WINDOW}, seed = {SEED}")
    counts, ratio = panel_a(out)
    sup = panel_b(out)
    contrib = panel_c(out)
    jac, fp = panel_d(out)
    rep = panel_e(out)
    tool_panel_matrix(out, counts, ratio, sup, contrib, jac, fp, rep)
    reconciliation(out, sup, ratio, jac)

    (REV / "logs" / "01_fig2_panel_inputs.log").write_text("\n".join(log_lines) + "\n")
    log(f"done in {time.time() - t0:.1f}s -> {out}")


if __name__ == "__main__":
    main()
