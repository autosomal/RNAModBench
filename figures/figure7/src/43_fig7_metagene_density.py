#!/usr/bin/env python3
"""43 -- Fig.7 panel D: replicate-aware non-m6A metagene densities (HeLa WT vs IVT).

Reuses ``common.regionmodel`` (Ensembl region model) and the density function of
``22_metagene_merge_density.py``.  For every non-m6A tool, each HeLa replicate's
clean callset is annotated to the transcript metagene coordinate; per-unit and
per-majority 5-segment densities are then computed.  The long table produced here
is stitched into panel D by ``47_fig7_rebuild.py``.

Data layer (read-only): ``harmonisation/callsets/RNA002/Human/{HeLa_WT,HeLa_IVT}/
{Psi,m1Psi,Nm,m5C}/<tool>/<sample>.tsv``  -- the cleaned, coordinate-correct
site layer (per the harmonisation README; 0-based BED, end = start + 1).

Annotation: Human Ensembl region model (GRCh38p14 / ensembl112), protein-coding
mRNA only -- the same model the m6A metagene pipeline uses, never GENCODE.

Output: figures/figure7/tables/fig7_metagene_density.tsv
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
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

HERE = _RB / "src/harmonisation"
sys.path.insert(0, str(HERE))

from common.config import SITES_ROOT                              # noqa: E402
from common.match import fix_chromosome                           # noqa: E402
from common.regionmodel import (KINDS, RegionIndex, SEGMENT_ORDER,  # noqa: E402
                                assign_sites)
from common.consensus import quorum                               # noqa: E402

SITES_CLEAN = (_RB / "data/callsets")
OUT_DIR = (_RB / "figures/figure7")
TAB = OUT_DIR / "tables"
LOG = OUT_DIR / "logs"

PLATFORM = "RNA002"
SPECIES = "Human"
MOD_DIRS = ["Psi", "m1Psi", "Nm", "m5C"]
CONDITIONS = {
    "HeLa_WT": ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"],
    "HeLa_IVT": ["HeLa_IVT_rep1", "HeLa_IVT_rep2", "HeLa_IVT_rep3"],
}

# density grid (must match 22_metagene_merge_density.py)
GRID = 200
FINE = 800
BANDWIDTH = 0.05


# --------------------------------------------------------------------------- #
def density(x: np.ndarray) -> np.ndarray:
    """Smooth density of ``x`` on [0, 1] (integral 1), reflect at the edges."""
    if x.size == 0:
        return np.zeros(GRID)
    hist, _ = np.histogram(x, bins=FINE, range=(0.0, 1.0))
    sm = gaussian_filter1d(hist.astype(float), BANDWIDTH * FINE, mode="reflect")
    sm = sm.reshape(GRID, FINE // GRID).sum(axis=1)
    total = sm.sum()
    return sm / total * GRID if total > 0 else np.zeros(GRID)


def load_replicate(path: Path, tag: str) -> pd.DataFrame | None:
    if not path.exists():
        return None
    d = pd.read_csv(path, sep="\t", low_memory=False)
    if d.empty:
        return None
    d["pos"] = pd.to_numeric(d["start"], errors="coerce")
    d = d.dropna(subset=["pos"])
    if d.empty:
        return None
    d["pos"] = d["pos"].astype("int64")
    d["strand"] = np.where(d["strand"].astype(str).isin(["-", "-1"]), "-", "+")
    d["chrom"] = d["chrom"].astype(str).map(fix_chromosome)
    out = d[["chrom", "pos", "strand"]].copy()
    out["unit"] = tag
    return out


def merge_replicates(frames: list[pd.DataFrame]) -> pd.DataFrame:
    a = pd.concat(frames, ignore_index=True)
    agg = a.groupby(["chrom", "pos", "strand"], sort=False).agg(
        n_units=("unit", "nunique"),
        units=("unit", lambda s: "|".join(sorted(set(s)))),
    ).reset_index()
    return agg


def masks_for(blk: pd.DataFrame) -> dict[str, np.ndarray]:
    n = blk["n_units"].to_numpy()
    k = int(n.size and blk["n_units_group"].iloc[0]) if "n_units_group" in blk else len(
        {t for row in blk["units"] for t in str(row).split("|")})
    out = {"union": n >= 1}
    if k >= 2:
        out["majority"] = n >= quorum(k)
        out["intersection"] = n >= k
    for t in sorted({t for row in blk["units"] for t in str(row).split("|")}):
        out[f"unit:{t}"] = blk["units"].str.contains(t, regex=False).to_numpy()
    return out


# --------------------------------------------------------------------------- #
def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    LOG.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    model = RegionIndex.load(sorted(
        ((_XB / "reference/regionmodels")).glob(
            f"{SPECIES}.*.mrna.regionmodel.pkl"))[0])

    # discover tools per mod per condition by scanning the clean tree
    tool_map: dict[str, str] = {}        # tool -> mod_type
    rows = []
    for cond, samples in CONDITIONS.items():
        gdir = SITES_CLEAN / PLATFORM / SPECIES / cond
        if not gdir.is_dir():
            print(f"missing {gdir}", flush=True)
            continue
        for mod in MOD_DIRS:
            mdir = gdir / mod
            if not mdir.is_dir():
                continue
            for tool_dir in sorted(p.name for p in mdir.iterdir() if p.is_dir()):
                tool_map.setdefault(tool_dir, mod)
                frames = [f for f in (load_replicate(
                    mdir / tool_dir / f"{s}.tsv", s) for s in samples)
                    if f is not None and len(f)]
                if not frames:
                    print(f"{cond:<9} {tool_dir:<16} ({mod}) -- no sites",
                          flush=True)
                    continue
                sites = merge_replicates(frames)
                k = len(frames)
                sites["n_units_group"] = k
                # annotate
                kind = np.full(len(sites), -1, dtype=np.int8)
                norm = np.full(len(sites), np.nan)
                for strand, sub in sites.groupby("strand", sort=False):
                    for chrom, ss in sub.groupby("chrom", sort=False):
                        si, ki, nz, _tl = assign_sites(
                            model, chrom, strand, ss["pos"].to_numpy(np.int64))
                        if len(si) == 0:
                            continue
                        pos = ss.index.to_numpy()[si]
                        kind[pos] = ki
                        norm[pos] = nz
                sites["kind"] = kind
                sites["norm_pos"] = norm
                ok = sites["kind"] >= 0
                if not ok.any():
                    print(f"{cond:<9} {tool_dir:<16} ({mod}) -- 0 assigned",
                          flush=True)
                    continue
                sites = sites.loc[ok].reset_index(drop=True)
                masks = masks_for(sites)
                for merge, mask in masks.items():
                    sub = sites.loc[mask]
                    tot = len(sub)
                    if tot == 0:
                        continue
                    for knd in SEGMENT_ORDER:
                        name = KINDS[knd]
                        vals = sub.loc[sub.kind == knd, "norm_pos"].to_numpy()
                        if vals.size == 0:
                            continue
                        dval = "|".join(
                            f"{v:.6g}" for v in density(vals) * (vals.size / tot))
                        rows.append({
                            "platform": PLATFORM, "species": SPECIES,
                            "mod_type": mod, "tool": tool_dir,
                            "condition": cond, "merge": merge, "kind": name,
                            "n_sites": int(vals.size),
                            "density": dval})
                maj = int((sites["n_units"] >= quorum(k)).sum())
                print(f"{cond:<9} {tool_dir:<16} ({mod}) units={k} "
                      f"sites={len(sites):>6} majority={maj:>6}", flush=True)

    res = pd.DataFrame(rows)
    out = TAB / "fig7_metagene_density.tsv"
    res.to_csv(out, sep="\t", index=False)
    print(f"wrote {len(res):,} density rows -> {out} "
          f"({time.time() - t0:.0f}s); tools={sorted(set(res.tool))}")


if __name__ == "__main__":
    main()
