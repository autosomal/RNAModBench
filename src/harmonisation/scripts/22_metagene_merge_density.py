#!/usr/bin/env python
"""Replicate-aware GUITAR metagene, stage 2: merge rules, densities, region shares.

Reads ``metagene_sites.tsv.gz`` (one row per unique site of a group x tool
callset, already merged over the group's independent units) and writes

  metagene_density.tsv.gz   per-unit / union / majority / intersection curves
  region_proportions.tsv    share of sites per transcript region, with the
                            across-unit spread for every merge rule
  merge_site_counts.tsv     how many sites each merge rule keeps, plus the
                            pairwise Jaccard between units (R3-2)

``majority`` = present in strictly more than half of the units (n // 2 + 1); it
is the primary curve of the revision figures because it keeps a site only when
the replicate structure supports it, while ``union`` reproduces the legacy
"pile all replicates together" behaviour and ``intersection`` is the strictest
view.  For a two-unit group the two rules coincide by definition.
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
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from common.config import SITES_ROOT                       # noqa: E402
from common.consensus import quorum                        # noqa: E402
from common.regionmodel import KINDS, NCR_ORDER, SEGMENT_ORDER  # noqa: E402
OUT = (_XB / "harmonisation/guitar_metagene")
TAB = OUT / "tables"

GRID = 200
FINE = 800
BANDWIDTH = 0.05
#: tool label for the pooled-across-tools curve that the published metagene shows
POOLED = "ALL (pooled)"
KEYS = ["platform", "species", "group", "tool", "merge_basis", "n_units_group",
        "tx_class", "strand_mode"]


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


def merge_masks(df: pd.DataFrame) -> dict[str, np.ndarray]:
    n = df["n_units"].to_numpy()
    k = int(df["n_units_group"].iloc[0])
    out = {"union": n >= 1}
    if k >= 2:
        out["majority"] = n >= quorum(k)
        out["intersection"] = n >= k
    tags = sorted({t for row in df["units"] for t in str(row).split("|")})
    for t in tags:
        out[f"unit:{t}"] = df["units"].str.contains(t, regex=False).to_numpy()
    return out


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", default=str(TAB / "metagene_sites.tsv.gz"))
    args = ap.parse_args()
    TAB.mkdir(parents=True, exist_ok=True)

    sites = pd.read_csv(args.sites, sep="\t", dtype={"units": str})
    sites["kind_name"] = sites["kind"].map(KINDS)
    print(f"{len(sites):,} annotated site rows", flush=True)

    dens, props, counts = [], [], []
    jac: dict[tuple, list[float]] = {}
    for key, blk in sites.groupby(KEYS, sort=False):
        kd = dict(zip(KEYS, key))
        k = int(blk["n_units_group"].iloc[0])
        order = NCR_ORDER if kd["tx_class"] == "ncrna" else SEGMENT_ORDER
        names = [KINDS[i] for i in order]
        masks = merge_masks(blk)
        unit_masks = {m: v for m, v in masks.items() if m.startswith("unit:")}

        for merge, mask in masks.items():
            sub = blk.loc[mask]
            tot = len(sub)
            if tot == 0:
                continue
            counts.append({**kd, "merge": merge, "n_sites": tot,
                           "sites_per_unit": round(tot / max(1, k), 1)})
            for name in names:
                vals = sub.loc[sub.kind_name == name, "norm_pos"].to_numpy()
                if vals.size == 0:
                    continue
                dens.append({**kd, "merge": merge, "kind": name,
                             "n_sites": int(vals.size),
                             "density": "|".join(
                                 f"{v:.6g}"
                                 for v in density(vals) * (vals.size / tot))})
                props.append({**kd, "merge": merge, "kind": name,
                              "n_sites": int(vals.size),
                              "share": round(vals.size / tot, 5)})
        # spread of the region shares across the individual units
        if len(unit_masks) >= 2:
            shares = {m: _shares(blk.loc[mask], names) for m, mask in unit_masks.items()}
            for name in names:
                v = [s[name] for s in shares.values() if not np.isnan(s[name])]
                if len(v) < 2:
                    continue
                props.append({**kd, "merge": "unit_mean", "kind": name,
                              "n_sites": int(round(np.mean([
                                  int((blk.loc[msk, "kind_name"] == name).sum())
                                  for msk in unit_masks.values()]))),
                              "share": round(float(np.mean(v)), 5),
                              "share_sd": round(float(np.std(v, ddof=1)), 5),
                              "n_units_observed": len(v)})
            sets = {m: set(blk.index[mask]) for m, mask in unit_masks.items()}
            jac[key] = [len(sets[a] & sets[b]) / max(1, len(sets[a] | sets[b]))
                        for a, b in itertools.combinations(sets, 2)]

    # pooled view: all tools of a library concatenated into one curve, as in
    # the published metagene (a directory of per-tool BEDs per condition)
    pooled_key = [k for k in KEYS if k != "tool"]
    for key, blk in sites.groupby(pooled_key, sort=False):
        kd = dict(zip(pooled_key, key))
        kd["tool"] = POOLED
        k = int(blk["n_units_group"].iloc[0])
        order = NCR_ORDER if kd["tx_class"] == "ncrna" else SEGMENT_ORDER
        names = [KINDS[i] for i in order]
        for merge, mask in merge_masks(blk).items():
            sub = blk.loc[mask]
            tot = len(sub)
            if tot == 0:
                continue
            counts.append({**kd, "merge": merge, "n_sites": tot,
                           "sites_per_unit": round(tot / max(1, k), 1),
                           "n_tools": int(sub.tool.nunique())})
            for name in names:
                vals = sub.loc[sub.kind_name == name, "norm_pos"].to_numpy()
                if vals.size == 0:
                    continue
                dens.append({**kd, "merge": merge, "kind": name,
                             "n_sites": int(vals.size),
                             "density": "|".join(
                                 f"{v:.6g}"
                                 for v in density(vals) * (vals.size / tot))})
                props.append({**kd, "merge": merge, "kind": name,
                              "n_sites": int(vals.size),
                              "share": round(vals.size / tot, 5)})

    for row in counts:
        j = jac.get(tuple(row[c] for c in KEYS))
        if j:
            row["mean_pairwise_jaccard"] = round(float(np.mean(j)), 4)
            row["min_pairwise_jaccard"] = round(float(np.min(j)), 4)

    pd.DataFrame(dens).to_csv(TAB / "metagene_density.tsv.gz", sep="\t",
                              index=False, compression="gzip")
    p = pd.DataFrame(props)
    for c in ("share_sd", "n_units_observed"):
        if c not in p.columns:
            p[c] = np.nan
    p.to_csv(TAB / "region_proportions.tsv", sep="\t", index=False)
    c = pd.DataFrame(counts)
    for col in ("mean_pairwise_jaccard", "min_pairwise_jaccard"):
        if col not in c.columns:
            c[col] = np.nan
    c.to_csv(TAB / "merge_site_counts.tsv", sep="\t", index=False)
    print(f"wrote {len(dens):,} density rows, {len(p):,} proportion rows, "
          f"{len(c):,} merge-count rows")


def _shares(sub: pd.DataFrame, names: list[str]) -> dict[str, float]:
    tot = len(sub)
    if tot == 0:
        return {n: np.nan for n in names}
    vc = sub["kind_name"].value_counts()
    return {n: float(vc.get(n, 0)) / tot for n in names}


if __name__ == "__main__":
    main()
