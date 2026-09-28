#!/usr/bin/env python
"""61 -- Figure S5 panel E evidence: the quality of the sites behind the trade-off.

Reviewer R3-4 asks what the *practical advantage* of the proposed tool
combinations is while the absolute recall stays low, and reviewer R3-m2 asks how
union-based coverage differs from intersection-based precision.  Both claims are
about *which sites* a strategy keeps, so this script characterises those sites on
the analysis-ready ``harmonisation/callsets`` layer.

Site sets, per independent unit, restricted to the shared measurable universe
(``consensus_eval.common_universe``: exonic, reference-base compatible,
coverage >= 10 in every unit of the group) with the 2-bp GLORI window:

    single           calls of the selected k = 1 tool (``fig6_combination_selected``)
    union            union of the selected k = 5 tools
    marginal         union - single   -- what the four added tools contribute
    intersection k=2 sites called by both tools of the selected pair
    intersection k=5 sites called by all five selected tools
    reference        GLORI sites inside the universe (the recall denominator)

Per set: n, GLORI-anchored fraction (2 bp), ``coverage`` (Nanom6A BAM read
support) median / quartiles / fraction >= 20, and DRACH fraction.  For the two
tools shared by all four selected combinations (m6Anet: probability_modified,
MINES: mod_ratio) the native score of the sites that tool called is summarised
separately -- native scores are *not* comparable across tools and are therefore
never pooled.

Outputs -> ``figures/figure6/tables/``
    figS5_site_quality.tsv         per unit x set metrics
    figS5_site_quality_scores.tsv  per unit x set x tool native-score stats
    figS5_coverage_hist.tsv        per unit x set log2 coverage counts (CDF panel)
    figS5_tool_quality.tsv         per unit x **tool** (all 13 configurations):
                                   n, GLORI-anchored fraction, recall inside the
                                   universe, coverage quartiles, DRACH fraction and
                                   the tool's own native-score stats -- the
                                   single-tool background against which the
                                   combinations are read (2026-09-21: both the
                                   rebuilt main Figure 6 and S5 need the per-tool
                                   plane, so every tool is loaded now, not only the
                                   members of the selected combinations)

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    figures/figureS6/src/61_figS6_sitequality.py
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
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

HERE = _RB / "src/harmonisation"
sys.path.insert(0, str(HERE))

from common import consensus_eval as ev                      # noqa: E402
from common.config import CALLSET_ROOT, SITES_ROOT, UNIVERSE_ROOT  # noqa: E402
from common.io_utils import write_table                      # noqa: E402
from common.manifest import Inventory, setup_logger          # noqa: E402
from common.match import fix_chromosome                      # noqa: E402

CLEAN_ROOT = (_RB / "data/callsets")
MIN_COV_GOOD = 20            #: "well covered" threshold used for the quality panel
SCORE_TOOLS = ("m6Anet", "MINES")
SET_ORDER = ("single (k=1)", "union (k=5)", "union marginal",
             "intersection (k=2)", "intersection (k=5)", "GLORI reference")
TABLE_COLUMNS = ["species", "group", "unit", "study", "set", "quality_source",
                 "n_sites", "n_glori", "frac_glori",
                 "n_cov", "coverage_median", "coverage_q1", "coverage_q3",
                 "frac_cov_ge20",
                 "n_drach", "frac_drach"]
SCORE_COLUMNS = ["species", "group", "unit", "set", "tool", "score_type",
                 "n_scored", "score_median", "score_q1", "score_q3"]
#: per-tool table (all 13 configurations): the single-tool background
TOOL_TABLE_COLUMNS = ["species", "group", "unit", "study", "tool",
                      "n_sites", "n_glori", "frac_glori", "recall_in_universe",
                      "n_cov", "coverage_median", "coverage_q1", "coverage_q3",
                      "frac_cov_ge20", "n_drach", "frac_drach",
                      "score_type", "n_scored", "score_median", "score_q1",
                      "score_q3"]
#: log2 coverage bins (0, 1, 2, 4 ... 262144 reads) for the coverage-CDF panel
COV_BINS = np.concatenate(([0.0], 2.0 ** np.arange(19)))
HIST_COLUMNS = ["species", "group", "unit", "set", "bin_lo", "bin_hi", "n_sites"]


def _load_fig6():
    """Import ``40_fig6_combination`` so both scripts share one scoring rule."""
    path = Path(__file__).resolve().parent / "40_fig6_combination.py"
    spec = importlib.util.spec_from_file_location("fig6_combination", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


FIG6 = _load_fig6()
OUT = FIG6.OUT
TAB, LOG = OUT / "tables", OUT / "logs"


# --------------------------------------------------------------------------- #
# per-callset loading: callset keys + callsets quality, row-aligned
# --------------------------------------------------------------------------- #
def _callset_keys(callsets: Path, clean: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read one callset twice (keys, quality) and verify they are row-aligned.

    ``callsets`` is exported from the callsets in the same row order, so the
    row-wise check below is both the cheapest and the strictest way to make sure
    the quality columns belong to the same calls.
    """
    keys = pd.read_csv(callsets, sep="\t", dtype=str, keep_default_na=False,
                       usecols=lambda c: c in ("chrom", "pos_raw", "pos_center",
                                               "strand"))
    qual = pd.read_csv(clean, sep="\t", dtype=str, keep_default_na=False,
                       usecols=lambda c: c in ("chrom", "start", "coverage",
                                               "drach", "score", "score_type"))
    if len(keys) != len(qual):
        raise SystemExit(f"row count differs: {callsets.name} "
                         f"({len(keys)}) vs {clean.name} ({len(qual)})")
    if not (keys["chrom"].to_numpy() == qual["chrom"].to_numpy()).all() or \
       not (pd.to_numeric(keys["pos_raw"]).to_numpy()
            == pd.to_numeric(qual["start"]).to_numpy()).all():
        raise SystemExit(f"row order differs: {callsets} vs {clean}")
    return keys, qual


def load_calls_with_quality(species: str, group: str, tool: str, sample: str,
                            universe: set, logger) -> dict:
    """{(chrom, pos): (coverage, drach, score)} inside the universe for one unit.

    ``pos`` is the evaluation coordinate (``pos_center`` where the strand is
    known, else ``pos_raw``) exactly as :func:`consensus_eval.site_keys` uses it;
    quality is taken from the row with the highest native score when a coordinate
    is duplicated by overlapping transcripts.
    """
    callsets = CALLSET_ROOT / FIG6.PLATFORM / species / group / FIG6.MOD / tool / f"{sample}.tsv"
    clean = CLEAN_ROOT / FIG6.PLATFORM / species / group / FIG6.MOD / tool / f"{sample}.tsv"
    if not callsets.exists():
        return {}
    if not clean.exists():
        raise SystemExit(f"callsets export missing for {callsets}")
    keys, qual = _callset_keys(callsets, clean)
    if keys.empty:
        return {}

    pos_raw = pd.to_numeric(keys["pos_raw"], errors="coerce")
    pos = pos_raw
    if "pos_center" in keys.columns:                     # same rule as site_keys()
        centred = pd.to_numeric(keys["pos_center"], errors="coerce")
        centred = centred.where(keys["strand"].astype(str).isin(["+", "-"]))
        pos = centred.fillna(pos_raw)
    score = pd.to_numeric(qual["score"], errors="coerce")
    drach = qual["drach"].astype(str).str.lower() if "drach" in qual.columns else None
    cov = pd.to_numeric(qual["coverage"], errors="coerce")
    stype = qual["score_type"].astype(str) if "score_type" in qual.columns else None

    df = pd.DataFrame({"chrom": keys["chrom"].astype(str).map(fix_chromosome),
                       "pos": pos, "cov": cov, "score": score,
                       "drach": (drach.eq("true") if drach is not None else pd.NA),
                       "stype": (stype if stype is not None else "")})
    df = df.dropna(subset=["pos"])
    df = (df.sort_values("score", na_position="first")
            .drop_duplicates(["chrom", "pos"], keep="last"))
    out: dict = {}
    # bracket access: ``df.cov`` would hit ``DataFrame.cov`` (the covariance method)
    for c, p, cov_v, sc, dr, st in zip(df["chrom"], df["pos"], df["cov"],
                                       df["score"], df["drach"], df["stype"]):
        key = (c, int(p))
        if key in universe:
            out[key] = (None if pd.isna(cov_v) else int(cov_v),
                        None if pd.isna(dr) else bool(dr),
                        None if pd.isna(sc) else float(sc), str(st))
    return out


def load_reference_quality(species: str, sample: str, ref_u: set,
                           logger) -> dict:
    """Coverage / DRACH of the GLORI-in-universe sites from the universe layer.

    Reading is chunked and filters on the reference positions while streaming:
    the per-sample universe files hold 3.5-5.4 M rows, the reference subset only
    ~1e4-1e5 of them.
    """
    path = UNIVERSE_ROOT / FIG6.PLATFORM / species / f"{sample}__{FIG6.MOD}.tsv"
    if not path.exists():
        raise SystemExit(f"universe file missing: {path}")
    by_chrom: dict[str, set] = defaultdict(set)
    for c, p in ref_u:
        by_chrom[c].add(p)
    arr_by_chrom = {c: np.sort(np.fromiter(v, dtype=np.int64, count=len(v)))
                    for c, v in by_chrom.items()}
    out: dict = {}
    for chunk in pd.read_csv(path, sep="\t", chunksize=1_000_000,
                             usecols=["chrom", "pos", "is_drach", "coverage"],
                             dtype={"chrom": str, "pos": np.int64,
                                    "coverage": np.int64}):
        # normalise only the (few) distinct labels of each chunk, so the join key
        # is the same ``fix_chromosome`` space the reference uses
        chunk["chrom"] = chunk["chrom"].map({c: fix_chromosome(c)
                                             for c in chunk["chrom"].unique()})
        for c, sub in chunk.groupby("chrom", sort=False):
            arr = arr_by_chrom.get(c)
            if arr is None or arr.size == 0:
                continue
            pos = sub["pos"].to_numpy(np.int64)
            idx = np.searchsorted(arr, pos)
            np.clip(idx, 0, arr.size - 1, out=idx)
            hit = arr[idx] == pos
            if not hit.any():
                continue
            dr = sub["is_drach"].astype(str).str.lower().eq("true").to_numpy()
            for p, d, cov in zip(pos[hit], dr[hit],
                                 sub["coverage"].to_numpy()[hit]):
                out[(c, int(p))] = (int(cov), bool(d), None, "")
    logger.info("[%s/%s] GLORI reference quality loaded: %d / %d sites",
                species, sample, len(out), len(ref_u))
    return out


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #
def set_metrics(keys: set, quality: dict, ref: dict, ref_n: int) -> dict:
    """Counts, GLORI anchoring and quality summaries of one site set."""
    n = len(keys)
    tp = 0
    by_chrom: dict[str, list] = defaultdict(list)
    for c, p in keys:
        by_chrom[c].append(p)
    for c, plist in by_chrom.items():
        arr = ref.get(c)
        if arr is None or arr.size == 0:
            continue
        p = np.asarray(plist, dtype=np.int64)
        tp += int(FIG6._hit_mask(p, arr, FIG6.WINDOW).sum())
    cov = np.asarray([quality[k][0] for k in keys
                      if quality.get(k) and quality[k][0] is not None],
                     dtype=float)
    dr = np.asarray([quality[k][1] for k in keys
                     if quality.get(k) and quality[k][1] is not None], dtype=bool)
    q1, med, q3 = (np.percentile(cov, [25, 50, 75]) if cov.size else (np.nan,) * 3)
    return {"n_sites": n,
            "n_glori": tp,
            "frac_glori": (tp / n) if n else np.nan,
            "n_cov": int(cov.size),
            "coverage_median": float(med),
            "coverage_q1": float(q1),
            "coverage_q3": float(q3),
            "frac_cov_ge20": float((cov >= MIN_COV_GOOD).mean()) if cov.size else np.nan,
            "n_drach": int(dr.size),
            "frac_drach": float(dr.mean()) if dr.size else np.nan}


def score_stats(values: list[float]) -> dict:
    if not values:
        return {"n_scored": 0, "score_median": np.nan,
                "score_q1": np.nan, "score_q3": np.nan}
    q1, med, q3 = np.percentile(np.asarray(values, dtype=float), [25, 50, 75])
    return {"n_scored": len(values), "score_median": float(med),
            "score_q1": float(q1), "score_q3": float(q3)}


# --------------------------------------------------------------------------- #
def main() -> None:
    logger = setup_logger("61_figS6_sitequality", log_dir=LOG)
    inv = Inventory("61_figS6_sitequality")
    sel = pd.read_csv(TAB / "fig6_combination_selected.tsv", sep="\t")
    rows: list[dict] = []
    score_rows: list[dict] = []
    hist_rows: list[dict] = []
    tool_rows: list[dict] = []

    for sp, (group, style) in FIG6.GROUPS.items():
        universe = ev.common_universe(FIG6.PLATFORM, sp, group)
        ref = ev.reference(sp)
        ref_u = ev.reference_positions_in_universe(ref, universe)
        ref_n = len(ref_u)
        units = ev.units(FIG6.PLATFORM, group).to_dict("records")
        logger.info("[%s/%s] universe=%d, GLORI-in-universe=%d, units=%d",
                    sp, group, len(universe), ref_n, len(units))

        for u in units:
            sample, tag, study = u["sample"], u["replicate_tag"], u["study"]
            gname = group if style == "replicates" else tag
            sub = sel[(sel.species == sp) & (sel.group == gname)]
            if len(sub) != 5:
                raise SystemExit(f"selection table incomplete for {sp}/{gname}")
            tools1 = tuple(sub[sub.k == 1].iloc[0]["combination"].split("+"))
            tools2 = tuple(sub[sub.k == 2].iloc[0]["combination"].split("+"))
            tools5 = tuple(sub[sub.k == 5].iloc[0]["combination"].split("+"))

            # every configuration is loaded, not only the members of the selected
            # combinations: the per-tool quality table needs all 13.  Tools that
            # are not needed for the set metrics are summarised and dropped right
            # away so the per-unit peak stays at the sets' own tool count.
            need = set(tools5) | set(tools2) | {tools1[0]}
            per_tool: dict[str, dict] = {}
            for t in FIG6.TOOLS:
                q = load_calls_with_quality(sp, group, t, sample, universe,
                                            logger)
                if not q:
                    logger.warning("[%s/%s] no callset for %s", sp, sample, t)
                    continue
                m = set_metrics(set(q), q, ref, ref_n)
                vals = [q[k][2] for k in q if q[k][2] is not None]
                tool_rows.append({
                    "species": sp, "group": gname, "unit": sample,
                    "study": study, "tool": t,
                    "recall_in_universe": (m["n_glori"] / ref_n) if ref_n else np.nan,
                    "score_type": next((q[k][3] for k in q if q[k][3]), ""),
                    **m, **score_stats(vals)})
                if t in need:
                    per_tool[t] = q
                else:
                    del q, m, vals
            calls5 = [set(per_tool[t]) for t in tools5]
            keys_single = set(per_tool[tools1[0]]) if tools1[0] in per_tool else set()
            keys_union = set().union(*calls5) if calls5 else set()
            keys_isect = set.intersection(*calls5) if calls5 else set()
            keys_isect2 = (set.intersection(*[set(per_tool[t]) for t in tools2])
                           if tools2 else set())
            keys_marg = keys_union - keys_single

            # merged quality of the union: coverage / DRACH are properties of the
            # position (Nanom6A BAM support, local 5-mer), so tools must agree --
            # disagreement means the merge is buggy and is reported, not hidden.
            merged: dict = {}
            conflicts = 0
            for t in tools5:
                for k, q in per_tool[t].items():
                    if k in merged:
                        if (q[0] != merged[k][0]) or (q[1] != merged[k][1]):
                            conflicts += 1
                    else:
                        merged[k] = q
            if conflicts:
                logger.warning("[%s/%s] %d coverage/DRACH conflicts across tools",
                               sp, sample, conflicts)

            refq = load_reference_quality(sp, sample, ref_u, logger)
            sets = {"single (k=1)": (keys_single, per_tool.get(tools1[0], {})),
                    "union (k=5)": (keys_union, merged),
                    "union marginal": (keys_marg, merged),
                    "intersection (k=2)": (keys_isect2, merged),
                    "intersection (k=5)": (keys_isect, merged),
                    "GLORI reference": (ref_u, refq)}
            for name in SET_ORDER:
                keys, qual = sets[name]
                m = set_metrics(keys, qual, ref, ref_n)
                rows.append({"species": sp, "group": gname, "unit": sample,
                             "study": study, "set": name,
                             "quality_source": ("universe"
                                                if name == "GLORI reference"
                                                else "callsets"),
                             **m})
                # coverage histogram (log2 bins) -> coverage-CDF panel
                covs = [qual[k][0] for k in keys
                        if qual.get(k) and qual[k][0] is not None]
                counts, _ = np.histogram(np.asarray(covs, dtype=float),
                                         bins=COV_BINS)
                for lo, hi, cnt in zip(COV_BINS[:-1], COV_BINS[1:], counts):
                    if cnt:
                        hist_rows.append({"species": sp, "group": gname,
                                          "unit": sample, "set": name,
                                          "bin_lo": float(lo), "bin_hi": float(hi),
                                          "n_sites": int(cnt)})
            for name in set(SET_ORDER) - {"GLORI reference"}:
                keys = sets[name][0]
                for t in SCORE_TOOLS:
                    q = per_tool.get(t)
                    if not q:
                        continue
                    vals = [q[k][2] for k in keys if k in q and q[k][2] is not None]
                    if not vals:          # the tool did not call any site of this set
                        continue
                    stype = next((q[k][3] for k in keys if k in q and q[k][3]), "")
                    score_rows.append({"species": sp, "group": gname, "unit": sample,
                                       "set": name, "tool": t,
                                       "score_type": stype, **score_stats(vals)})
            logger.info("[%s/%s] single=%d union=%d marginal=%d isect2=%d isect5=%d",
                        sp, sample, len(keys_single), len(keys_union),
                        len(keys_marg), len(keys_isect2), len(keys_isect))
            del per_tool, merged

    quality = pd.DataFrame(rows)[TABLE_COLUMNS]
    scores = pd.DataFrame(score_rows)[SCORE_COLUMNS]
    hist = pd.DataFrame(hist_rows)[HIST_COLUMNS]
    toolq = pd.DataFrame(tool_rows)[TOOL_TABLE_COLUMNS]
    write_table(quality, TAB / "figS5_site_quality.tsv")
    write_table(scores, TAB / "figS5_site_quality_scores.tsv")
    write_table(hist, TAB / "figS5_coverage_hist.tsv")
    write_table(toolq, TAB / "figS5_tool_quality.tsv")
    inv.flush()
    logger.info("tables -> %s (%d rows) / %s (%d rows) / %s (%d rows) / "
                "%s (%d rows)",
                TAB / "figS5_site_quality.tsv", len(quality),
                TAB / "figS5_site_quality_scores.tsv", len(scores),
                TAB / "figS5_coverage_hist.tsv", len(hist),
                TAB / "figS5_tool_quality.tsv", len(toolq))


if __name__ == "__main__":
    main()
