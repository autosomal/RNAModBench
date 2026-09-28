"""Score call sets on ONE shared measurable universe.

``evaluation/tables/reproducibility.tsv`` reports per-replicate rows on each
replicate's own universe and merged rows on the *common* universe, so recall
from one cannot be compared with recall from the other.  The revision figures
that compare "a single replicate" with "the consensus of the replicates"
therefore recompute every set on the same universe here:

    universe      exonic, reference-base-compatible positions with coverage
                  >= MIN_COV in EVERY independent unit of the group
    call set      per-replicate, or the majority consensus
                  (``common.consensus.quorum`` units)
    scoring       nearest GLORI site within ``window`` bp

Used by ``24_replicate_structure.py`` and
``27_window_combination.py``.
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
import numpy as np
import pandas as pd

from .config import CALLSET_ROOT, GLORI, SITES_ROOT, UNIVERSE_ROOT
from .consensus import quorum
from .match import fix_chromosome

MIN_COV = 10


def units(platform: str, group: str) -> pd.DataFrame:
    """Independent units of a group (nested subsets removed, runs deduplicated)."""
    reg = pd.read_csv((_RB / "metadata/sample_registry.csv"),
                      sep="\t", dtype=str)
    reg = reg[(reg.platform == platform) & (reg.dataset_group == group)
              & (reg.independence_class != "nested_subset")]
    return reg.drop_duplicates("sequencing_unit")


def site_keys(path):
    d = pd.read_csv(path, sep="\t",
                    usecols=lambda c: c in ("chrom", "pos_center", "pos_raw", "strand"))
    if d.empty:
        return pd.DataFrame(columns=["chrom", "pos"])
    pos = pd.to_numeric(d["pos_raw"], errors="coerce")
    if "pos_center" in d.columns:
        # ``pos_center`` is the centre-corrected coordinate and is only defined
        # when the transcript strand is known -- for a strand-less row it would
        # silently move the site by up to +/- 5 bp (Nanom6A wrote '*' for every
        # row).  ``32_impute_strand`` fills the strand; whatever stays unknown
        # keeps ``pos_raw``.
        centred = pd.to_numeric(d["pos_center"], errors="coerce")
        if "strand" in d.columns:
            centred = centred.where(d["strand"].astype(str).isin(["+", "-"]))
        pos = centred.fillna(pos)
    out = pd.DataFrame({"chrom": d["chrom"].astype(str).map(fix_chromosome),
                        "pos": pos}).dropna()
    out["pos"] = out["pos"].astype("int64")
    return out.drop_duplicates()


def replicate_sets(platform: str, species: str, group: str, tool: str,
                   mod: str = "m6A") -> dict[str, set]:
    """{replicate_tag: site set} for one tool of one group."""
    root = CALLSET_ROOT / platform / species / group / mod / tool
    out = {}
    for r in units(platform, group).to_dict("records"):
        f = root / f"{r['sample']}.tsv"
        if f.exists():
            out[r["replicate_tag"]] = set(map(tuple, site_keys(f).to_numpy()))
    return out


def consensus_sets(platform: str, species: str, group: str,
                   mod: str = "m6A") -> dict[str, dict]:
    """{tool: {"sites": majority consensus set, "n_units": units with data}}."""
    root = CALLSET_ROOT / platform / species / group / mod
    out = {}
    for tool in sorted(p.name for p in root.iterdir() if p.is_dir()):
        per = replicate_sets(platform, species, group, tool, mod)
        if not per:
            continue
        count: dict[tuple, int] = {}
        for s in per.values():
            for k in s:
                count[k] = count.get(k, 0) + 1
        thr = quorum(len(per))
        out[tool] = {"sites": {k for k, v in count.items() if v >= thr},
                    "n_units": len(per)}
    return out


def common_universe(platform: str, species: str, group: str,
                    mod: str = "m6A") -> set:
    """Positions measurable in every independent unit of the group."""
    common = None
    for r in units(platform, group).to_dict("records"):
        hits = sorted((UNIVERSE_ROOT / platform / species)
                      .glob(f"{r['sample']}__{mod}.tsv"))
        if not hits:
            continue
        u = pd.read_csv(hits[0], sep="\t", usecols=["chrom", "pos", "coverage"])
        u = u[u.coverage >= MIN_COV]
        keys = set(zip(u.chrom.astype(str).map(fix_chromosome),
                       u.pos.astype("int64")))
        common = keys if common is None else (common & keys)
    return common or set()


def reference(species: str) -> dict[str, np.ndarray]:
    """GLORI positions per chromosome, sorted, for one species."""
    g = pd.read_csv(GLORI[species], sep="\t", header=None, usecols=[0, 1])
    pos: dict[str, list[int]] = {}
    for c, p in zip(g[0].astype(str).map(fix_chromosome), g[1].astype(int)):
        pos.setdefault(c, []).append(p)
    return {c: np.sort(np.asarray(v)) for c, v in pos.items()}


def score(sites: set, universe: set, ref: dict[str, np.ndarray],
          ref_in_universe: set, window: int = 2) -> dict:
    """TP/FP and universe-restricted precision/recall of one call set."""
    in_u = {s for s in sites if s in universe}
    tp = 0
    for chrom, p in in_u:
        arr = ref.get(chrom)
        if arr is None or arr.size == 0:
            continue
        i = int(np.searchsorted(arr, p))
        cand = {max(0, i - 1), min(i, len(arr) - 1), min(i + 1, len(arr) - 1)}
        tp += min(abs(int(arr[j]) - p) for j in cand) <= window
    n = len(ref_in_universe)
    return dict(n_calls=len(in_u), tp=tp,
                precision=tp / len(in_u) if in_u else np.nan,
                n_reference=n, recall=tp / n if n else np.nan)


def reference_positions_in_universe(ref: dict[str, np.ndarray],
                                    universe: set) -> set:
    """Reference sites that lie in the shared measurable universe."""
    out = set()
    uni_by_chrom: dict[str, set] = {}
    for c, p in universe:
        uni_by_chrom.setdefault(c, set()).add(int(p))
    for c, arr in ref.items():
        u = uni_by_chrom.get(c)
        if not u:
            continue
        out |= {(c, int(p)) for p in arr if int(p) in u}
    return out
