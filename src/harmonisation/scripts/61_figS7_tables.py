#!/usr/bin/env python3
"""61 -- Figure S7 evidence tables (replicate-aware non-m6A ncRNA metagene).

Figure S7 (``$RNAMODBENCH_LOCAL/submission/02_AS_working_copy_and_revisions/sup/sup7.pdf``) shows the ncRNA
metagene of the six non-m6A tools in HeLa WT versus the unmodified IVT negative
control.  The published version pooled the three replicates (the legacy Guitar
inputs were an undocumented union) and carried no quantitative statement, which
is exactly what reviewer points R3-2/E6 (replicate structure) and R3-9 (why the
non-m6A tools over-call) ask to fix.  This script produces the numbers behind the
rebuilt figure:

* every profile is computed per independent sequencing unit (HeLa WT1-3,
  HeLa IVT rep1-3) plus the majority consensus (called in >= 2 of 3 units,
  ``common.consensus.quorum``);
* WT-versus-unmodified-IVT similarity is quantified by the Jensen-Shannon
  divergence of the mass profile (JSD in bits, [0, 1]), the L1 distance and the
  change of the ncRNA-body share, so "the distributions are similar" becomes a
  number instead of an eyeball call;
* the chance level is a **chromosome-stratified permutation null**: positions are
  redrawn from the sample's own candidate universe restricted to the ncRNA region
  model, with the same per-chromosome counts as the observed majority call set,
  R = 1000, seed = 20260919 (the R3-9 seed, so the analyses stay comparable).

Conventions (identical to R3-9 / Figure S6)
-------------------------------------------
* analysis layer = ``harmonisation/callsets`` (0-based BED, coordinate fixes);
* candidate universe = ``harmonisation/universe/<platform>/<species>/<canonical>__universe.tsv.gz``
  filtered to ``coverage >= 10`` and to bases that can carry the modification
  (Nm = A/C/G/T, Psi/m1Psi = T/A, m5C = C/G) -- the geometry helpers live in
  ``r39_build_evidence.py`` and are imported, never copied;
* annotation = Ensembl ``Homo_sapiens.GRCh38.112.chr.gtf`` **ncRNA** transcripts
  (everything that is not ``protein_coding``); GENCODE is banned project-wide.
  A non-coding transcript is drawn as one body plus 1 kb upstream/downstream
  windows (the ``1kb | ncRNA | 1kb`` axis of the published figure), so the
  UTR/CDS pieces of the few non-coding biotypes that carry a CDS annotation
  (e.g. ``nonsense_mediated_decay``) fold into the body, matching the Guitar
  ``pltTxType = "ncrna"`` semantics used to draw the row;
* a position is assigned to the **longest** transcript containing it, over both
  strands, by one implementation used for the observed calls and for the
  background pool alike (no dependence on the imputed call strand);
* no m6A-centred reference is used anywhere (R2-2); the unmodified IVT libraries
  are negative controls (E5), never a treatment arm;
* nothing is written outside ``$RNAMODBENCH_ROOT`` (no /tmp).

Outputs (``figures/figureS8/tables/``)
---------------------------------------------------------
s7_units.tsv               calls / in-universe / ncRNA-assigned counts per unit + majority
s7_density_profiles.tsv    200-bin density per tool x condition x profile x segment
s7_segment_shares.tsv      share of calls in each of the three metagene segments
s7_profile_distance.tsv    WT-IVT pairs (9 unit pairs + majority): JSD / L1 / shape JSD
s7_background_profiles.tsv candidate-universe profile of each tool x condition
s7_null_distribution.tsv   R = 1000 permutation draws per tool and statistic
s7_null_summary.tsv        observed vs null quantiles + empirical p
s7_anchor_check.tsv        every value vs the frozen R3-9 evidence tables
s7_geometry.tsv            run metadata (model, grids, pool sizes, timings)

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    src/harmonisation/scripts/61_figS7_tables.py
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
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

HERE = Path(__file__).resolve()
PROJECT = _RB
SITES_V2 = (_RB / "src/harmonisation")
sys.path.insert(0, str(SITES_V2))

from common.consensus import quorum                        # noqa: E402
from common.manifest import Inventory, setup_logger        # noqa: E402
from common.match import fix_chromosome                    # noqa: E402
from common.regionmodel import (CDS, F3, F5, NCR, UTR3,    # noqa: E402
                                UTR5, RegionIndex)

# --------------------------------------------------------------------------- #
# frozen inputs / outputs
# --------------------------------------------------------------------------- #
R39 = (_RB / "analysis/nonm6a_false_positives")
EV = (_RB / "analysis/nonm6a_false_positives/evidence")
OUT = (_RB / "figures/figureS8")
TAB = (_RB / "figures/figureS8/tables")
LOG = (_RB / "figures/figureS8/logs")

#: single implementation of the universe / geometry helpers (no drift, see S6)
_spec = importlib.util.spec_from_file_location(
    "r39_build_evidence", (_RB / "analysis/nonm6a_false_positives/analysis/r39_build_evidence.py"))
r39 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(r39)

MODEL_DIR = (_XB / "reference/regionmodels")
MIN_COV = 10
GRID = 200                    # bins per metagene segment (matches 22_*)
FINE = 800                    # histogram resolution before smoothing
BANDWIDTH = 0.05              # gaussian sigma in transcript units
SEED = 20260919               # same seed as the R3-9 permutation null
N_DRAW = 1000                 # permutation draws per tool
SPECIES = "Human"
#: every statistic is tested against a null computed at its own sample size,
#: otherwise a bigger call set would look "better than chance" for free
STAT_KEYS = ("jsd_wt_ivt_unit", "jsd_wt_ivt_majority", "jsd_wt_background",
             "jsd_ivt_background")

#: display names and the three classes resolved for R3-9 / Figure 7
CLASSES: list[tuple[str, list[str]]] = [
    ("FP-dominated", ["CHEUI_m5C", "NanoMUD_psi", "NanoMUD_m1psi"]),
    ("Intermediate", ["NanoNm"]),
    ("Sparse", ["NanoPsu", "NanoSPA_psU"]),
]
TOOL_ORDER: list[str] = [t for _, ts in CLASSES for t in ts]
CLASS_OF: dict[str, str] = {t: c for c, ts in CLASSES for t in ts}
DISPLAY: dict[str, str] = {
    "CHEUI_m5C": "CHEUI-m5C",
    "NanoMUD_psi": "NanoMUD-\u03a8",
    "NanoMUD_m1psi": "NanoMUD-m1\u03a8",
    "NanoNm": "NanoNm",
    "NanoPsu": "NanoPsu",
    "NanoSPA_psU": "NanoSPA-\u03a8",
}
MOD_OF: dict[str, str] = {tool: mod for mod, tool in r39.TOOLMOD}
UNITS: dict[str, list[str]] = {c: list(u) for c, u in r39.HELA_GROUPS.items()}
ALL_UNITS: list[str] = UNITS["WT"] + UNITS["IVT"]

#: metagene segments of the ncRNA axis, in plot order (regionmodel.NCR_ORDER)
SEGMENTS: list[tuple[int, str]] = [
    (F5, "five_prime_flank"), (NCR, "ncRNA_body"), (F3, "three_prime_flank")]
SEG_NAMES = [name for _, name in SEGMENTS]
#: exonic kinds that make up a non-coding transcript body (folded into NCR)
BODY_KINDS = [NCR, UTR5, CDS, UTR3]

#: genome-orientation base codes that can carry each modification
BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
MOD_BASE_CODES: dict[str, list[int]] = {
    "Nm": [0, 1, 2, 3], "m5C": [1, 2], "Psi": [0, 3], "m1Psi": [0, 3]}

anchors: list[dict] = []
_LOGGER = None


def _log():
    return _LOGGER


def record(quantity: str, value: float, ref: float, *, rtol: float = 1e-5,
           note: str = "") -> None:
    """Anchor one computed value against its frozen R3-9 counterpart."""
    ok = bool(np.isfinite(value) and np.isfinite(ref)
              and np.isclose(value, ref, rtol=rtol, atol=1e-12))
    anchors.append({"quantity": quantity, "s7_value": value, "r39_value": ref,
                    "rel_diff": (abs(value - ref) / abs(ref)) if ref else np.nan,
                    "status": "OK" if ok else "MISMATCH", "note": note})


# --------------------------------------------------------------------------- #
# density / distance helpers (density() matches 22_metagene_merge_density.py)
# --------------------------------------------------------------------------- #
def density(x: np.ndarray) -> np.ndarray:
    """Smooth density of ``x`` on [0, 1] (sum = ``GRID``), reflect at the edges."""
    if x.size == 0:
        return np.zeros(GRID)
    hist, _ = np.histogram(x, bins=FINE, range=(0.0, 1.0))
    sm = gaussian_filter1d(hist.astype(float), BANDWIDTH * FINE, mode="reflect")
    sm = sm.reshape(GRID, FINE // GRID).sum(axis=1)
    return sm / sm.sum() * GRID if sm.sum() > 0 else np.zeros(GRID)


def jsd(p: np.ndarray, q: np.ndarray) -> float:
    """Jensen-Shannon divergence in bits (inputs normalised to sum 1)."""
    p = np.asarray(p, float)
    q = np.asarray(q, float)
    if p.sum() <= 0 or q.sum() <= 0:
        return float("nan")
    p, q = p / p.sum(), q / q.sum()
    m = 0.5 * (p + q)
    tiny = np.finfo(float).tiny

    def _kl(a: np.ndarray, b: np.ndarray) -> float:
        keep = a > 0
        return float(np.sum(a[keep] * np.log2(a[keep] / np.maximum(b[keep], tiny))))

    return 0.5 * _kl(p, m) + 0.5 * _kl(q, m)


def l1_mass(p: np.ndarray, q: np.ndarray) -> float:
    """L1 (2 x total variation) distance of two mass profiles."""
    if p.sum() <= 0 or q.sum() <= 0:
        return float("nan")
    return float(np.abs(p / p.sum() - q / q.sum()).sum())


def profile(kinds: np.ndarray, nz: np.ndarray) -> dict:
    """Mass / shape profile of one placed site set over the three segments.

    ``mass``   sums to 1 over the whole axis (``share_seg * shape_seg``)
    ``shape``  each segment normalised to 1 (within-segment distribution)
    ``shares`` fraction of the sites in each of the three segments
    """
    out = {"n": int(kinds.size), "mass": np.zeros(GRID * len(SEGMENTS)),
           "shape": np.zeros(GRID * len(SEGMENTS)), "shares": {}}
    body = np.isin(kinds, BODY_KINDS)
    for i, (knd, name) in enumerate(SEGMENTS):
        vals = nz[body] if name == "ncRNA_body" else nz[kinds == knd]
        share = vals.size / out["n"] if out["n"] else 0.0
        dens = density(vals)                      # sums to GRID, or all zeros
        out["shares"][name] = share
        sl = slice(i * GRID, (i + 1) * GRID)
        if dens.sum() > 0:
            out["mass"][sl] = dens / GRID * share
            out["shape"][sl] = dens / GRID
    return out


# --------------------------------------------------------------------------- #
# genomic placement: ncRNA body / 1 kb flanks, longest transcript, both strands
# --------------------------------------------------------------------------- #
def chrom_table(model: RegionIndex, chrom: str) -> dict | None:
    """Merge both strands of one chromosome into one start-sorted table."""
    parts = [(s, model.tables[(chrom, s)]) for s in ("+", "-")
             if (chrom, s) in model.tables]
    if not parts:
        return None
    order = np.argsort(np.concatenate([t["start"] for _, t in parts]),
                       kind="stable")
    out: dict[str, np.ndarray] = {}
    for key in ("start", "end", "kind", "tx_start", "seg_len", "tx_len"):
        out[key] = np.concatenate([t[key] for _, t in parts])[order]
    out["plus"] = np.concatenate(
        [np.full(t["start"].size, s == "+") for s, t in parts])[order]
    return out


def place_positions(tab: dict, pos: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Place sorted positions: longest-transcript segment + transcript coordinate.

    Mirrors ``regionmodel.assign_sites(tx_rule="longest")`` but resolves the
    transcript over both strands (candidate-universe positions carry no strand)
    and folds every exonic piece of a non-coding transcript into the body
    segment, so the axis is ``1kb | ncRNA | 1kb``.
    """
    n = pos.size
    kind = np.full(n, -1, np.int8)
    nz = np.zeros(n, np.float32)
    best = np.zeros(n, np.int64)
    starts, ends = tab["start"], tab["end"]
    lo_all = np.searchsorted(pos, starts, side="left")
    hi_all = np.searchsorted(pos, ends, side="left")
    kinds, txs = tab["kind"], tab["tx_start"]
    lens, txlen, plus = tab["seg_len"], tab["tx_len"], tab["plus"]
    for j in range(starts.size):
        lo, hi = int(lo_all[j]), int(hi_all[j])
        if hi <= lo:
            continue
        sl = slice(lo, hi)
        better = txlen[j] > best[sl]
        if not better.any():
            continue
        idx = np.flatnonzero(better) + lo
        best[idx] = txlen[j]
        k = int(kinds[j])
        kind[idx] = NCR if k in BODY_KINDS else k
        p = pos[idx]
        local = p - starts[j] if plus[j] else ends[j] - 1 - p
        nz[idx] = (txs[j] + local + 0.5) / lens[j]
    return kind, nz


def load_sample(sample: str, model: RegionIndex, log) -> dict:
    """Candidate universe (coverage >= 10) of one library with its placement."""
    t0 = time.time()
    df = pd.read_csv(r39.universe_path(sample), sep="\t",
                     usecols=["chrom", "pos", "base", "coverage"],
                     dtype={"chrom": str, "base": str})
    df = df[df["coverage"] >= MIN_COV]
    df["chrom"] = [fix_chromosome(c) for c in df["chrom"]]
    df["pos"] = df["pos"].astype(np.int64)
    df["base"] = df["base"].map(BASE_CODE).fillna(4).astype(np.int8)
    out: dict[str, dict[str, np.ndarray]] = {}
    for chrom, sub in df.groupby("chrom", sort=False):
        sub = sub.sort_values("pos", kind="stable")
        pos = sub["pos"].to_numpy(np.int64)
        keep = np.concatenate([[True], pos[1:] != pos[:-1]])
        pos = pos[keep]
        base = sub["base"].to_numpy(np.int8)[keep]
        tab = chrom_table(model, chrom)
        if tab is None:
            kind = np.full(pos.size, -1, np.int8)
            nz = np.zeros(pos.size, np.float32)
        else:
            kind, nz = place_positions(tab, pos)
        out[chrom] = {"pos": pos, "base": base, "kind": kind, "nz": nz}
    del df
    n_cov = int(sum(v["pos"].size for v in out.values()))
    n_pool = int(sum((v["kind"] >= 0).sum() for v in out.values()))
    log.info("    %-14s cov>=%2d %10d | ncRNA pool %10d (%.1fs)",
             sample, MIN_COV, n_cov, n_pool, time.time() - t0)
    return out


def mod_pool(data: dict, mod: str) -> dict[str, dict[str, np.ndarray]]:
    """Subset of one library's universe that can carry ``mod`` inside ncRNA."""
    codes = MOD_BASE_CODES[mod]
    out: dict[str, dict[str, np.ndarray]] = {}
    for chrom, d in data.items():
        keep = np.isin(d["base"], codes) & (d["kind"] >= 0)
        if keep.any():
            out[chrom] = {"pos": d["pos"][keep], "kind": d["kind"][keep],
                          "nz": d["nz"][keep]}
    return out


def lookup(data: dict, chrom: str, pos: np.ndarray, codes: list[int]
           ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(in_universe, kind, nz) of ``pos`` in one library's universe.

    ``in_universe`` marks positions present with a base that can carry the
    modification (the R3-9 / Figure S6 definition); ``kind`` is -1 outside the
    ncRNA model (or outside the universe altogether).
    """
    n = pos.size
    inu = np.zeros(n, bool)
    kind = np.full(n, -1, np.int8)
    nz = np.zeros(n, np.float32)
    d = data.get(chrom)
    if d is None or d["pos"].size == 0 or n == 0:
        return inu, kind, nz
    j = np.searchsorted(d["pos"], pos)
    jj = np.clip(j, 0, d["pos"].size - 1)
    hit = d["pos"][jj] == pos
    hit &= np.isin(d["base"][jj], codes)
    inu[hit] = True
    kind[hit] = d["kind"][jj[hit]]
    nz[hit] = d["nz"][jj[hit]]
    return inu, kind, nz


def union_pool(pools: list[dict[str, dict[str, np.ndarray]]]
               ) -> dict[str, dict[str, np.ndarray]]:
    """Union of the units' pools of one condition (placement is genomic)."""
    out: dict[str, dict[str, np.ndarray]] = {}
    keys = set().union(*[set(p) for p in pools]) if pools else set()
    for chrom in keys:
        pos, kind, nz = [], [], []
        for p in pools:
            d = p.get(chrom)
            if d is not None:
                pos.append(d["pos"])
                kind.append(d["kind"])
                nz.append(d["nz"])
        allpos, allkind, allnz = (np.concatenate(pos), np.concatenate(kind),
                                  np.concatenate(nz))
        order = np.argsort(allpos, kind="stable")
        allpos, allkind, allnz = allpos[order], allkind[order], allnz[order]
        keep = np.concatenate([[True], allpos[1:] != allpos[:-1]])
        out[chrom] = {"pos": allpos[keep], "kind": allkind[keep],
                      "nz": allnz[keep]}
    return out


def draw_profile(pool: dict[str, dict[str, np.ndarray]],
                 counts: dict[str, int], rng: np.random.Generator) -> dict:
    """One chromosome-stratified draw from a condition's candidate pool."""
    kinds, nzs = [], []
    for chrom, n in counts.items():
        d = pool.get(chrom)
        if d is None or n <= 0 or d["pos"].size == 0:
            continue
        idx = rng.integers(0, d["pos"].size, size=int(n))
        kinds.append(d["kind"][idx])
        nzs.append(d["nz"][idx])
    kinds = np.concatenate(kinds) if kinds else np.zeros(0, np.int8)
    nzs = np.concatenate(nzs) if nzs else np.zeros(0, np.float32)
    return profile(kinds, nzs)


# --------------------------------------------------------------------------- #
# per-tool bookkeeping
# --------------------------------------------------------------------------- #
def collect_calls(res, tool: str, samples: dict) -> dict:
    """Raw calls, in-universe counts and ncRNA-placed sites, per unit."""
    mod = MOD_OF[tool]
    codes = MOD_BASE_CODES[mod]
    out: dict[str, dict[str, dict]] = {}
    for cond, units in UNITS.items():
        out[cond] = {}
        for sample in units:
            calls = res.calls(sample, mod, tool)
            n_in = n_on = 0
            kinds, nzs = [], []
            in_uni: dict[str, np.ndarray] = {}
            placed: dict[str, np.ndarray] = {}
            for chrom, pos in calls.items():
                if pos.size == 0:
                    continue
                inu, k, z = lookup(samples[sample], chrom, pos, codes)
                if inu.any():
                    in_uni[chrom] = pos[inu]
                    n_in += int(inu.sum())
                on = inu & (k >= 0)
                if on.any():
                    kinds.append(k[on])
                    nzs.append(z[on])
                    placed[chrom] = pos[on]
                    n_on += int(on.sum())
            out[cond][sample] = {
                "mod": mod, "n_raw": r39.n_positions(calls), "n_in": n_in,
                "n_on": n_on, "in_uni": in_uni, "placed": placed,
                "kind": np.concatenate(kinds) if kinds else np.zeros(0, np.int8),
                "nz": np.concatenate(nzs) if nzs else np.zeros(0, np.float32)}
    return out


def union_counts(per_chrom: dict[str, list[np.ndarray]]) -> dict[str, int]:
    """Union size per chromosome of the per-unit position lists."""
    out: dict[str, int] = {}
    for chrom, arrays in per_chrom.items():
        out[chrom] = int(np.unique(np.concatenate(arrays)).size)
    return out


def majority_sites(per_chrom: dict[str, list[np.ndarray]], samples: dict,
                   units: list[str], codes: list[int]
                   ) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    """Majority-consensus (>= quorum) ncRNA sites plus their per-chromosome count."""
    k = len(units)
    kinds, nzs = [], []
    counts: dict[str, int] = {}
    for chrom, arrays in per_chrom.items():
        allp = np.concatenate(arrays)
        uniq, c = np.unique(allp, return_counts=True)
        keep = c >= quorum(k)
        counts[chrom] = int(keep.sum())
        if not keep.any():
            continue
        mpos = uniq[keep]
        kind = np.full(mpos.size, -1, np.int8)
        nz = np.zeros(mpos.size, np.float32)
        todo = np.ones(mpos.size, bool)
        for s in units:                        # placement is genomic: any unit does
            if not todo.any():
                break
            inu, kk, zz = lookup(samples[s], chrom, mpos[todo], codes)
            got = inu & (kk >= 0)
            if got.any():
                idx = np.flatnonzero(todo)[got]
                kind[idx] = kk[got]
                nz[idx] = zz[got]
                todo[idx] = False
        ok = kind >= 0
        if ok.any():
            kinds.append(kind[ok])
            nzs.append(nz[ok])
    kinds = np.concatenate(kinds) if kinds else np.zeros(0, np.int8)
    nzs = np.concatenate(nzs) if nzs else np.zeros(0, np.float32)
    return kinds, nzs, counts


# --------------------------------------------------------------------------- #
def main() -> None:
    for d in (TAB, LOG):
        d.mkdir(parents=True, exist_ok=True)
    log = setup_logger("61_figS7_tables", log_dir=LOG)
    global _LOGGER
    _LOGGER = log
    inv = Inventory("61_figS7_tables")
    t0 = time.time()

    models = sorted(p for p in MODEL_DIR.glob(f"{SPECIES}.*.ncrna.regionmodel.pkl")
                    if "gencode" not in p.name)
    if len(models) != 1:
        raise SystemExit("expected exactly one Ensembl human ncRNA model, "
                         f"found {[p.name for p in models]}")
    log.info("Figure S7 ncRNA evidence | model=%s | seed=%d R=%d",
             models[0].name, SEED, N_DRAW)
    model = RegionIndex.load(models[0])

    log.info("[1/6] candidate universe + ncRNA placement (per library)")
    samples = {s: load_sample(s, model, log) for s in ALL_UNITS}

    log.info("[2/6] universe-size anchors vs the frozen R3-9 tables")
    per_rep = pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/hela_wt_ivt_per_replicate.tsv"), sep="\t")
    summ = pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/hela_wt_ivt_summary.tsv"), sep="\t")
    for mod in ("Nm", "m5C", "Psi", "m1Psi"):
        codes = MOD_BASE_CODES[mod]
        for sample in ALL_UNITS:
            mine = int(sum(np.isin(samples[sample][c]["base"], codes).sum()
                           for c in samples[sample]))
            row = per_rep[(per_rep.mod_type == mod) & (per_rep["sample"] == sample)]
            if row.empty:
                continue
            record(f"universe_n/{mod}/{sample}", float(mine),
                   float(row["universe_n"].iloc[0]),
                   note="coverage>=10 + modification-compatible base")

    log.info("[3/6] per tool: counts, profiles, distances, permutation null")
    res = r39.Resources(min_cov=MIN_COV)
    pool_cache: dict[tuple[str, str], dict] = {}   # (mod, condition) -> union pool
    unit_rows: list[dict] = []
    share_rows: list[dict] = []
    dens_rows: list[dict] = []
    dist_rows: list[dict] = []
    bg_rows: list[dict] = []
    null_rows: list[dict] = []
    null_summary: list[dict] = []
    geometry: list[dict] = []

    for tool in TOOL_ORDER:
        mod = MOD_OF[tool]
        codes = MOD_BASE_CODES[mod]
        calls = collect_calls(res, tool, samples)
        pools = {}
        for cond, units in UNITS.items():
            if (mod, cond) not in pool_cache:
                pool_cache[(mod, cond)] = union_pool(
                    [mod_pool(samples[s], mod) for s in units])
            pools[cond] = pool_cache[(mod, cond)]
        pool_size = {c: int(sum(v["pos"].size for v in pools[c].values()))
                     for c in pools}
        geometry.append({"tool": tool, "display": DISPLAY[tool], "mod_type": mod,
                         "pool_positions_WT": pool_size["WT"],
                         "pool_positions_IVT": pool_size["IVT"]})

        # ---- counts and profiles per unit, then the majority consensus ------
        prof: dict[str, dict[str, dict]] = {"WT": {}, "IVT": {}}
        obs_counts: dict[str, dict[str, int]] = {}
        for cond, units in UNITS.items():
            for sample in units:
                rec = calls[cond][sample]
                p = profile(rec["kind"], rec["nz"])
                prof[cond][sample] = p
                unit_rows.append(dict(
                    tool=tool, display=DISPLAY[tool], mod_type=mod,
                    tool_class=CLASS_OF[tool], condition=cond, profile=sample,
                    unit_kind="unit", n_units_used=1,
                    n_calls_raw=rec["n_raw"], n_calls_in_universe=rec["n_in"],
                    n_calls_on_ncrna=rec["n_on"],
                    share_five_prime_flank=p["shares"]["five_prime_flank"],
                    share_ncRNA_body=p["shares"]["ncRNA_body"],
                    share_three_prime_flank=p["shares"]["three_prime_flank"]))
                record(f"per-replicate in-universe calls/{tool}/{sample}",
                       float(rec["n_in"]),
                       float(per_rep[(per_rep.tool == tool) &
                                     (per_rep["sample"] == sample)]
                             ["n_calls_in_universe"].iloc[0]))
            # union of the in-universe calls (R3-9 summary anchor)
            for c, units_ in UNITS.items():
                uni = union_counts({chrom: [calls[c][s]["in_uni"][chrom]
                                            for s in units_
                                            if chrom in calls[c][s]["in_uni"]]
                                    for chrom in set().union(*[
                                        set(calls[c][s]["in_uni"]) for s in units_])})
                record(f"union in-universe calls/{tool}/{c}",
                       float(sum(uni.values())),
                       float(summ[summ.tool == tool]
                             [f"{c.lower()}_union_in_universe"].iloc[0]))
            # majority-consensus sites + per-chromosome permutation budget
            per_chrom = {chrom: [calls[cond][s]["placed"][chrom]
                                 for s in units if chrom in
                                 calls[cond][s]["placed"]]
                         for chrom in set().union(*[
                             set(calls[cond][s]["placed"]) for s in units])}
            maj_kind, maj_nz, maj_counts = majority_sites(per_chrom, samples,
                                                          units, codes)
            obs_counts[cond] = maj_counts
            p = profile(maj_kind, maj_nz)
            prof[cond]["majority"] = p
            unit_rows.append(dict(
                tool=tool, display=DISPLAY[tool], mod_type=mod,
                tool_class=CLASS_OF[tool], condition=cond, profile="majority",
                unit_kind="majority", n_units_used=len(units),
                n_calls_raw=np.nan, n_calls_in_universe=int(maj_kind.size),
                n_calls_on_ncrna=int(maj_kind.size),
                share_five_prime_flank=p["shares"]["five_prime_flank"],
                share_ncRNA_body=p["shares"]["ncRNA_body"],
                share_three_prime_flank=p["shares"]["three_prime_flank"]))
            log.info("    %-14s %-4s %-9s calls=%6d in-uni=%6d on-ncRNA=%6d "
                     "(majority %d)", tool, mod, cond,
                     sum(r["n_raw"] for r in [calls[cond][s] for s in units]),
                     sum(calls[cond][s]["n_in"] for s in units),
                     sum(calls[cond][s]["n_on"] for s in units),
                     int(maj_kind.size))

        # ---- WT-IVT distances: 9 unit pairs + the majority pair -------------
        for ws in UNITS["WT"]:
            for vs in UNITS["IVT"]:
                a, b = prof["WT"][ws], prof["IVT"][vs]
                dist_rows.append(dict(
                    tool=tool, display=DISPLAY[tool], mod_type=mod,
                    tool_class=CLASS_OF[tool], comparison="unit_pair",
                    profile_a=ws, profile_b=vs, n_a=a["n"], n_b=b["n"],
                    jsd_mass=jsd(a["mass"], b["mass"]),
                    l1_mass=l1_mass(a["mass"], b["mass"]),
                    jsd_shape=jsd(a["shape"], b["shape"]),
                    d_share_ncRNA_body=(a["shares"]["ncRNA_body"]
                                        - b["shares"]["ncRNA_body"])))
        a, b = prof["WT"]["majority"], prof["IVT"]["majority"]
        obs = dict(tool=tool, display=DISPLAY[tool], mod_type=mod,
                   tool_class=CLASS_OF[tool], comparison="majority",
                   profile_a="WT_majority", profile_b="IVT_majority",
                   n_a=a["n"], n_b=b["n"],
                   jsd_mass=jsd(a["mass"], b["mass"]),
                   l1_mass=l1_mass(a["mass"], b["mass"]),
                   jsd_shape=jsd(a["shape"], b["shape"]),
                   d_share_ncRNA_body=(a["shares"]["ncRNA_body"]
                                       - b["shares"]["ncRNA_body"]))
        dist_rows.append(obs)

        # ---- segment shares / density long tables ---------------------------
        for cond in ("WT", "IVT"):
            for name, p in prof[cond].items():
                share_rows.append(dict(
                    tool=tool, display=DISPLAY[tool], condition=cond,
                    profile=name, n_sites=p["n"],
                    share_five_prime_flank=p["shares"]["five_prime_flank"],
                    share_ncRNA_body=p["shares"]["ncRNA_body"],
                    share_three_prime_flank=p["shares"]["three_prime_flank"]))
                for i, seg in enumerate(SEG_NAMES):
                    for b_ in range(GRID):
                        dens_rows.append(dict(
                            tool=tool, display=DISPLAY[tool], condition=cond,
                            profile=name, segment=seg, bin=b_,
                            x_local=(b_ + 0.5) / GRID,
                            density=float(p["mass"][i * GRID + b_] * GRID),
                            mass=float(p["mass"][i * GRID + b_]),
                            shape=float(p["shape"][i * GRID + b_])))

        # ---- candidate-universe background profile per condition ------------
        bg = {}
        for cond in ("WT", "IVT"):
            kinds = (np.concatenate([v["kind"] for v in pools[cond].values()])
                     if pools[cond] else np.zeros(0, np.int8))
            nzs = (np.concatenate([v["nz"] for v in pools[cond].values()])
                   if pools[cond] else np.zeros(0, np.float32))
            bg[cond] = profile(kinds, nzs)
            for i, seg in enumerate(SEG_NAMES):
                for b_ in range(GRID):
                    bg_rows.append(dict(
                        tool=tool, display=DISPLAY[tool], mod_type=mod,
                        condition=cond, segment=seg, bin=b_,
                        x_local=(b_ + 0.5) / GRID,
                        mass=float(bg[cond]["mass"][i * GRID + b_]),
                        shape=float(bg[cond]["shape"][i * GRID + b_])))

        # ---- permutation null ------------------------------------------------
        # the majority statistic is tested at the majority sample size and the
        # unit statistic at the mean unit size, so neither can look better than
        # chance merely because it contains more sites
        unit_counts: dict[str, dict[str, int]] = {}
        for cond, units in UNITS.items():
            per_chrom: dict[str, list[int]] = {}
            for s in units:
                for chrom, pos in calls[cond][s]["placed"].items():
                    per_chrom.setdefault(chrom, []).append(int(pos.size))
            unit_counts[cond] = {
                chrom: int(round(float(np.mean(v + [0] * (len(units) - len(v))))))
                for chrom, v in per_chrom.items()}
        rng = np.random.default_rng(SEED + TOOL_ORDER.index(tool))
        draws: dict[str, list[float]] = {s: [] for s in STAT_KEYS}
        for _ in range(N_DRAW):
            d_maj = {cond: draw_profile(pools[cond], obs_counts[cond], rng)
                     for cond in ("WT", "IVT")}
            d_unit = {cond: draw_profile(pools[cond], unit_counts[cond], rng)
                      for cond in ("WT", "IVT")}
            draws["jsd_wt_ivt_majority"].append(jsd(d_maj["WT"]["mass"],
                                                   d_maj["IVT"]["mass"]))
            draws["jsd_wt_ivt_unit"].append(jsd(d_unit["WT"]["mass"],
                                                d_unit["IVT"]["mass"]))
            draws["jsd_wt_background"].append(jsd(d_maj["WT"]["mass"],
                                                  bg["WT"]["mass"]))
            draws["jsd_ivt_background"].append(jsd(d_maj["IVT"]["mass"],
                                                   bg["IVT"]["mass"]))
        obs_stat = {
            "jsd_wt_ivt_unit": float(np.mean(
                [r["jsd_mass"] for r in dist_rows
                 if r["tool"] == tool and r["comparison"] == "unit_pair"])),
            "jsd_wt_ivt_majority": obs["jsd_mass"],
            "jsd_wt_background":
                jsd(prof["WT"]["majority"]["mass"], bg["WT"]["mass"]),
            "jsd_ivt_background":
                jsd(prof["IVT"]["majority"]["mass"], bg["IVT"]["mass"])}
        for stat in STAT_KEYS:
            values = np.asarray(draws[stat], float)
            for d_i, v in enumerate(values):
                null_rows.append(dict(tool=tool, display=DISPLAY[tool],
                                      statistic=stat, draw=d_i, value=float(v)))
            o = obs_stat[stat]
            null_summary.append(dict(
                tool=tool, display=DISPLAY[tool], mod_type=mod,
                tool_class=CLASS_OF[tool], statistic=stat, observed=o,
                null_median=float(np.median(values)),
                null_p2_5=float(np.percentile(values, 2.5)),
                null_p97_5=float(np.percentile(values, 97.5)),
                null_max=float(values.max()),
                empirical_p=(1 + int((values >= o).sum())) / (N_DRAW + 1),
                n_draws=N_DRAW, seed=SEED + TOOL_ORDER.index(tool),
                pool_WT=pool_size["WT"], pool_IVT=pool_size["IVT"],
                counts_basis="per-chromosome counts of the observed majority "
                             "call set"))
        idx_u = len(null_summary) - len(STAT_KEYS)
        log.info("    %-14s (%-6s) pool WT=%7d IVT=%7d | unit obs=%.4f "
                 "null=%.4f p=%.3f | majority obs=%.4f null=%.4f p=%.3f",
                 tool, mod, pool_size["WT"], pool_size["IVT"],
                 obs_stat["jsd_wt_ivt_unit"],
                 float(np.median(draws["jsd_wt_ivt_unit"])),
                 null_summary[idx_u]["empirical_p"], obs["jsd_mass"],
                 float(np.median(draws["jsd_wt_ivt_majority"])),
                 null_summary[idx_u + 1]["empirical_p"])

    log.info("[4/6] writing tables")
    for df, name in (
            (pd.DataFrame(unit_rows), "s7_units.tsv"),
            (pd.DataFrame(dens_rows), "s7_density_profiles.tsv"),
            (pd.DataFrame(share_rows), "s7_segment_shares.tsv"),
            (pd.DataFrame(dist_rows), "s7_profile_distance.tsv"),
            (pd.DataFrame(bg_rows), "s7_background_profiles.tsv"),
            (pd.DataFrame(null_rows), "s7_null_distribution.tsv"),
            (pd.DataFrame(null_summary), "s7_null_summary.tsv"),
            (pd.DataFrame(geometry), "s7_geometry.tsv")):
        df.to_csv(TAB / name, sep="\t", index=False, float_format="%.6g")
        inv.record(TAB / name, n_rows=len(df))
        log.info("     %-30s %8d rows", name, len(df))

    log.info("[5/6] anchor check")
    adf = pd.DataFrame(anchors)
    adf.to_csv((_RB / "figures/figureS8/tables/s7_anchor_check.tsv"), sep="\t", index=False,
               float_format="%.6g")
    inv.record((_RB / "figures/figureS8/tables/s7_anchor_check.tsv"))
    bad = adf[adf.status != "OK"]
    log.info("     %d/%d anchors OK", len(adf) - len(bad), len(adf))
    for r in bad.itertuples():
        log.warning("     MISMATCH %s: s7=%.10g r39=%.10g", r.quantity,
                    r.s7_value, r.r39_value)

    log.info("[6/6] done in %.1f s; tables -> %s", time.time() - t0, TAB)
    inv.flush()
    if len(bad):
        raise SystemExit(f"{len(bad)} anchor mismatches; see "
                         f"{TAB / 's7_anchor_check.tsv'}")


if __name__ == "__main__":
    main()
