#!/usr/bin/env python
"""Replicate-aware GUITAR metagene, stage 1: annotate every m6A callset site
with its transcript region and normalised position.

Why this exists
---------------
The published metagene panels (Fig. 3C/D and the RNA004 equivalent) were drawn
from per-tool BED files that the legacy assembly produced from a single
replicate and with ``strand <- "+"`` forced, so neither replicate agreement nor
strand could be read off the figure.  ``harmonisation/callsets`` now holds one file
per replicate, so the same plot can be rebuilt with an explicit merge rule.

This script writes ONE row per unique site of a (group x tool) callset, already
merged over the group's *independent* units, with

    n_units / units   how many units called it, and which ones
    tx_class          mrna (protein-coding) or ncrna transcript model used
    strand_mode       aware  = site strand honoured (primary)
                      legacy = strand forced to "+" (published behaviour)
    kind              five_prime_flank / five_prime_UTR / CDS / three_prime_UTR
                      / three_prime_flank / ncRNA_body
    norm_pos          position within that segment, 0-1 in transcript direction

Merge rules, densities and figures are stage 2 (``22_metagene_merge_density.py``).

Outputs -> $RNAMODBENCH_LOCAL/guitar_metagene/tables/metagene_sites.tsv.gz
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
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from common.config import CALLSET_ROOT, SITES_ROOT                    # noqa: E402
from common.regionmodel import MODEL_DIR, RegionIndex, assign_sites    # noqa: E402

OUT = (_XB / "harmonisation/guitar_metagene")
TAB = OUT / "tables"
LOG = OUT / "logs"
REGISTRY = (_RB / "metadata/sample_registry.csv")

SPECIES = ("Arabidopsis", "Mouse", "Human")
MOD = "m6A"
ATTRS = ["is_drach", "coverage", "dist_glori", "mod_ratio"]


# --------------------------------------------------------------------------- #
# replicate structure
# --------------------------------------------------------------------------- #
def unit_table() -> pd.DataFrame:
    """One row per callable unit, plus the basis on which units may be merged.

    ``nested_subset`` rows are dropped: they are a depth-matched slice of another
    run, not an extra replicate.
    """
    reg = pd.read_csv(REGISTRY, sep="\t", dtype=str)
    reg = reg[reg.independence_class != "nested_subset"].copy()
    reg = reg.drop_duplicates(["platform", "dataset_group", "sequencing_unit"])

    def basis(row) -> str:
        if row.independence_class == "cross_study":
            return "studies"
        if row.species == "E.coli":
            return "sequencing_runs"
        return "replicates"

    reg["merge_basis"] = reg.apply(basis, axis=1)
    return reg[["platform", "species", "dataset_group", "sample",
                "replicate_tag", "sequencing_unit", "merge_basis"]]


# --------------------------------------------------------------------------- #
# callsets
# --------------------------------------------------------------------------- #
def load_unit_sites(path: Path, tag: str) -> pd.DataFrame | None:
    if not path.exists():
        return None
    d = pd.read_csv(path, sep="\t", low_memory=False)
    if d.empty:
        return None
    raw = pd.to_numeric(d["pos_raw"], errors="coerce")
    if "pos_center" in d.columns:
        pc = pd.to_numeric(d["pos_center"], errors="coerce")
        # the centre correction is only defined with a known strand
        pc = pc.where(d["strand"].astype(str).isin(["+", "-", "-1"]))
        raw = pc.fillna(raw)
    d["pos"] = raw
    d = d.dropna(subset=["pos"])
    d["pos"] = d["pos"].astype("int64")
    d["strand"] = np.where(d["strand"].astype(str).isin(["-", "-1"]), "-", "+")
    d["chrom"] = d["chrom"].astype(str)
    if "is_drach_center" in d.columns:
        d["is_drach"] = d["is_drach_center"].astype(str).isin(["True", "TRUE", "1"])
    else:
        d["is_drach"] = np.nan
    for src, dst in (("coverage", "coverage"), ("dist_glori", "dist_glori"),
                     ("mod_ratio", "mod_ratio")):
        d[dst] = pd.to_numeric(d.get(src), errors="coerce") if src in d.columns else np.nan
    cols = ["chrom", "pos", "strand", "is_drach", "coverage", "dist_glori",
            "mod_ratio"]
    out = d[cols].dropna(subset=["pos"]).copy()
    out["unit"] = tag
    return out


def merge_units(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Collapse per-unit callsets to one row per unique site."""
    a = pd.concat(frames, ignore_index=True)
    a["is_drach"] = a["is_drach"].astype(float)
    agg = a.groupby(["chrom", "pos", "strand"], sort=False).agg(
        n_units=("unit", "nunique"),
        units=("unit", lambda s: "|".join(sorted(set(s)))),
        is_drach=("is_drach", "max"),
        coverage=("coverage", "median"),
        dist_glori=("dist_glori", "min"),
        mod_ratio=("mod_ratio", "median"),
    ).reset_index()
    return agg


# --------------------------------------------------------------------------- #
# region annotation
# --------------------------------------------------------------------------- #
_cache: dict[tuple[str, str], RegionIndex] = {}


def get_model(species: str, cls: str) -> RegionIndex:
    key = (species, cls)
    if key not in _cache:
        hits = sorted(MODEL_DIR.glob(f"{species}.*.{cls}.regionmodel.pkl"))
        if not hits:
            raise SystemExit(f"no region model for {key}; run 20_build_region_model.py")
        _cache[key] = RegionIndex.load(hits[0])
    return _cache[key]


def annotate(sites: pd.DataFrame, idx: RegionIndex, *,
             strand_mode: str) -> tuple[np.ndarray, np.ndarray]:
    """Return (kind, norm_pos) aligned to ``sites``; kind = -1 when unassigned."""
    kind = np.full(len(sites), -1, dtype=np.int8)
    norm = np.full(len(sites), np.nan)
    for strand, rows in sites.groupby("strand", sort=False):
        use = strand if strand_mode == "aware" else "+"
        for chrom, sub in rows.groupby("chrom", sort=False):
            si, ki, nz, _tl = assign_sites(
                idx, chrom, use, sub["pos"].to_numpy(np.int64))
            if len(si) == 0:
                continue
            pos = sub.index.to_numpy()[si]
            kind[pos] = ki
            norm[pos] = nz
    return kind, norm


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--platform", default="RNA002")
    ap.add_argument("--strand-mode", default="both",
                    choices=["aware", "legacy", "both"])
    ap.add_argument("--groups", default=None, help="comma-separated subset")
    ap.add_argument("--tools", default=None, help="comma-separated subset")
    args = ap.parse_args()

    TAB.mkdir(parents=True, exist_ok=True)
    LOG.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    reg = unit_table()
    modes = ["aware", "legacy"] if args.strand_mode == "both" else [args.strand_mode]
    root = CALLSET_ROOT / args.platform
    groups = [g for g in reg.loc[reg.platform == args.platform, "dataset_group"]
              .unique()
              if any((root / s / g / MOD).is_dir() for s in SPECIES)]
    if args.groups:
        want = set(args.groups.split(","))
        groups = [g for g in groups if g in want]

    chunks = []
    for group in groups:
        species = next(s for s in SPECIES if (root / s / group / MOD).is_dir())
        units = reg[(reg.platform == args.platform)
                    & (reg.dataset_group == group)].to_dict("records")
        gdir = root / species / group / MOD
        tools = sorted(p.name for p in gdir.iterdir() if p.is_dir())
        if args.tools:
            want = set(args.tools.split(","))
            tools = [t for t in tools if t in want]
        for tool in tools:
            frames = [f for f in (load_unit_sites(gdir / tool / f"{u['sample']}.tsv",
                                                  u["replicate_tag"])
                                  for u in units) if f is not None and len(f)]
            if not frames:
                continue
            sites = merge_units(frames)
            n_units = len(units)
            for mode in modes:
                for tx_class in ("mrna", "ncrna"):
                    k, nz = annotate(sites, get_model(species, tx_class),
                                     strand_mode=mode)
                    ok = k >= 0
                    if not ok.any():
                        continue
                    out = sites.loc[ok].copy()
                    out["kind"] = k[ok]
                    out["norm_pos"] = np.round(nz[ok], 5)
                    for c in ("platform", "species", "group", "tool",
                              "tx_class", "strand_mode", "merge_basis",
                              "n_units_group"):
                        out[c] = {"platform": args.platform, "species": species,
                                  "group": group, "tool": tool,
                                  "tx_class": tx_class, "strand_mode": mode,
                                  "merge_basis": units[0]["merge_basis"],
                                  "n_units_group": n_units}[c]
                    chunks.append(out)
            print(f"{group:<18} {tool:<16} units={n_units} "
                  f"sites={len(sites):>7} "
                  f"maj={int((sites.n_units >= n_units // 2 + 1).sum()):>7} "
                  f"all={int((sites.n_units == n_units).sum()):>7}", flush=True)

    res = pd.concat(chunks, ignore_index=True)
    res = res[["platform", "species", "group", "tool", "merge_basis",
               "n_units_group", "tx_class", "strand_mode", "chrom", "pos",
               "strand", "units", "n_units", "kind", "norm_pos", "is_drach",
               "coverage", "dist_glori", "mod_ratio"]]
    res = res.drop_duplicates()
    out = TAB / "metagene_sites.tsv.gz"
    res.to_csv(out, sep="\t", index=False, compression="gzip")
    print(f"wrote {len(res):,} annotated site rows -> {out} "
          f"({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
