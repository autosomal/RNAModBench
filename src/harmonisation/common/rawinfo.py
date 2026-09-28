"""Per-site raw-output enrichment for the clean site collection.

The working callsets were extracted from the ``converted_callsets`` layer,
whose files carry only 7 columns (Chr/Start/End/Status/Prob/Strand/mod_ratio).
The *raw* tool outputs under ``$RNAMODBENCH_LOCAL/raw/result/<tool>/<sample>/`` hold
much more per-site information -- coverage, stoichiometry, probability,
k-mer/motif, p-values, odds ratios, read counts (the author's "5th column and
beyond").

Two join paths
--------------
**Genomic tools** (Dorado/modkit, ELIGOS2, yanocomp, EpiNano_Error) publish
genomic coordinates themselves -- the callset ``start`` is joined *directly*
against the raw position, trying a small set of candidate offsets / score
columns and keeping only the combination that reproduces the callset score on
>= ``MIN_MATCH`` of the matched rows.

**Transcript-space tools** (CHEUI, DENA, MINES, m6Anet, DRUMMER, xPore,
Nanocompore, NanoMUD) are bridged by the ``*_liftover.txt`` R2Dtool junction
table that sits next to the converted file (or, when the callset's
``source_file`` points at the raw output, one directory up in
``converted_callsets/<tool>/<sample>/``): it holds both the genomic block
(cols 1-6) and the transcript block (cols 7-9) of the very same sites, so

    callset.pos_raw  ==  bridge.genomic_start            (exact, 0-based)
    bridge.tx_start  ==  raw.position + offset           (per-callset constant)

The offset is **inferred and verified per callset**: the verification scores
are the bridge's own Prob / mod_ratio columns against the callset score
(rounding-tolerant, R2Dtool writes 2 decimals), falling back to raw-vs-bridge
score agreement.  A tool/callset whose scores do not verify is skipped --
enrichment may be missing, never wrong.
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
from pathlib import Path

import numpy as np
import pandas as pd

from .config import RESULT_RNA002
from .match import fix_chromosome

#: R2Dtool converted layer -- the ``*_liftover.txt`` junction tables live here
CONVERTED_ROOT = Path(str(_XB / "raw/converted_callsets"))

_SUMMARY: pd.DataFrame | None = None


def _summary_lookup() -> dict:
    """(sample, tool, mod) -> absolute source path from 01's own manifest."""
    global _SUMMARY
    if _SUMMARY is None:
        p = Path(str(_XB / "metadata/callsets_summary.csv"))
        _SUMMARY = {}
        if p.exists():
            df = _read(p)
            for _, r in df.iterrows():
                if r.get("source_file"):
                    _SUMMARY[(str(r["sample"]), str(r["tool"]), str(r["mod_type"]))] = \
                        str(Path(str(_RB)) / r["source_file"])
        return _SUMMARY
    return _SUMMARY


def _resolve_source(source_file: str, tool: str, sample: str, mod: str) -> str:
    """Absolute source path: the callset value, else 01's manifest."""
    p = Path(source_file)
    if p.is_absolute() and p.exists():
        return str(p)
    hit = _summary_lookup().get((sample, tool, mod))
    if hit and Path(hit).exists():
        return hit
    return hit or source_file


_RTOL, _ATOL = 1e-4, 1e-6
#: 2-dp rounding tolerance -- R2Dtool writes probabilities with 2 decimals
_ROUND = 0.00505
#: minimum share of score-matched rows that must agree on one offset
MIN_MATCH = 0.95
#: minimum number of score-matched rows needed to infer an offset at all
MIN_SCORED = 10

EXTRA_COLUMNS = [
    "coverage", "counts_mod", "counts_canonical", "counts_other_mod",
    "counts_delete", "counts_fail", "counts_diff", "counts_nocall",
    "counts_total", "kmer", "kmer7", "probability", "prob_unmodified",
    "predicted_rate", "calibrated_rate", "stoichiometry", "pvalue",
    "pvalue_ttest", "adj_pvalue", "odds_ratio", "logit_lor", "logodds",
    "esb_test", "esb_ctrl", "sum_err_ko", "sum_err_wt", "lm_z_score",
    "stat_11", "stat_12", "stat_13", "stat_14", "stat_15", "stat_16",
    "stat_17", "stat_18", "stat_19",
]

BRIDGE_SCORE_COLS = (10, 12)  # Prob (col 11) and mod_ratio/frac_diff (col 13)

#: tools whose raw outputs are too large to enrich interactively
#: (ELIGOS2 ``*_baseExt0.txt`` ≈ 410 MB/sample x 47 callsets,
#:  EpiNano ``*prediction.csv`` ≈ 615 MB x 4 files x 22 callsets).
#: Enrichment is skipped for them -- add back here once a cheaper join exists.
SLOW_TOOLS = {"ELIGOS2_diff", "ELIGOS2_solo", "EpiNano_Error"}


def _read(path: Path, sep: str = "\t", header="infer") -> pd.DataFrame:
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, header=header)


def _num(s) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def _close(a, b) -> np.ndarray:
    a, b = _num(a), _num(b)
    return np.isclose(a, b, rtol=_RTOL, atol=_ATOL) | (a - b).abs().le(_ROUND)


def _keys(s) -> pd.Series:
    """Normalise chromosome/transcript labels on *both* sides of a join --
    callsets spell ``chr1`` where R2Dtool bridges spell ``1``."""
    return s.astype(str).map(fix_chromosome)


def _bridge_path(source_file: str, tool: str, sample: str, mod: str) -> Path | None:
    """The R2Dtool junction table that belongs to a converted callset source."""
    p = Path(source_file)
    if p.exists():
        for cand in (p.name.replace("_remove_chr.txt", "_liftover.txt"),
                     p.name.replace(".txt", "_liftover.txt")):
            q = p.with_name(cand)
            if q.exists():
                return q
    # the callset may point at the *raw* output -- the junction table sits in
    # the converted_callsets layer instead
    d = CONVERTED_ROOT / tool / sample
    if d.is_dir():
        hits = sorted(d.glob("*_liftover.txt"))
        if mod:
            filt = [h for h in hits if mod in h.name]
            if filt:
                hits = filt
            elif any(t in h.name for t in ("m6A", "m5C", "m1psi", "psi", "inosine")
                     for h in hits):
                return None  # bridge exists but belongs to a different mod
        if hits:
            return hits[0]
    return None


# --------------------------------------------------------------------------- #
# transcript-space raw readers -> (df, key_col, pos_col, score_cols, extras)
# --------------------------------------------------------------------------- #
def _raw_cheui(sample, mod):
    df = _read(RESULT_RNA002 / "CHEUI" / sample / f"site_level_{mod}_predictions.txt")
    return df, "contig", "position", ("stoichiometry",), {
        "coverage": "coverage", "kmer": "site", "probability": "probability",
        "stoichiometry": "stoichiometry"}


def _raw_dena(sample, mod):
    df = _read(RESULT_RNA002 / "DENA" / sample / f"{sample}.tsv", header=None)
    df.columns = ["tid", "position", "kmer", "counts_mod", "counts_total",
                  "ratio"][:df.shape[1]]
    return df, "tid", "position", ("ratio",), {
        "kmer": "kmer", "counts_mod": "counts_mod", "coverage": "counts_total"}


def _raw_mines(sample, mod):
    df = _read(RESULT_RNA002 / "MINES" / sample / f"{sample}.bed", header=None)
    df.columns = ["tid", "start", "end", "kmer", "name", "strand", "frac",
                  "coverage", "stat"][:df.shape[1]]
    return df, "tid", "start", ("frac",), {
        "kmer": "kmer", "probability": "frac", "coverage": "coverage"}


def _raw_m6anet(sample, mod):
    base = RESULT_RNA002 / "m6Anet"
    for cand in (base / sample / "data.site_proba.csv",
                 base / f"{sample}_transcripts_inference" / "data.site_proba.csv"):
        if cand.exists():
            df = _read(cand, sep=",")
            return df, "transcript_id", "transcript_position", \
                ("mod_ratio", "probability_modified"), {
                "coverage": "n_reads", "probability": "probability_modified",
                "kmer": "kmer"}
    return None, None, None, None, {}


def _raw_drummer(sample, mod):
    base = RESULT_RNA002 / "DRUMMER"
    for cand in (base / sample / "summary.txt",
                 base / f"{sample}_result1" / "summary.txt",
                 base / f"{sample}_result2" / "summary.txt"):
        if cand.exists():
            df = _read(cand)
            return df, "transcript_id", "transcript_pos", ("frac_diff",), {
                "coverage": "depth_treat", "kmer": "eleven_bp_motif",
                "pvalue": "OR_padj", "odds_ratio": "odds_ratio"}
    return None, None, None, None, {}


def _raw_xpore(sample, mod):
    df = _read(RESULT_RNA002 / "xPore" / sample / "diffmod.table", sep=",")
    return df, "id", "position", ("diff_mod_rate_ko_vs_wt",), {
        "kmer": "kmer", "pvalue": "pval_ko_vs_wt"}


def _raw_nanocompore(sample, mod):
    df = _read(RESULT_RNA002 / "Nanocompore" / sample / "outnanocompore_results.tsv")
    return df, "ref_id", "pos", (), {
        "kmer": "ref_kmer", "pvalue": "GMM_logit_pvalue", "logit_lor": "Logit_LOR"}


def _raw_nanomud(sample, mod):
    name = "psi_out.csv" if mod == "psi" else "m1psi_out.csv"
    df = _read(RESULT_RNA002 / "NanoMUD" / sample / name, sep=",")
    return df, "chrom", "chrom_pos", ("predicted_rate", "calibrated_rate"), {
        "kmer": "motif", "coverage": "coverage",
        "predicted_rate": "predicted_rate", "calibrated_rate": "calibrated_rate"}


#: tool -> (reader, bridge score column index to verify against, mod suffix)
SPECS = {
    "CHEUI_m6A": (_raw_cheui, 12, "m6A"),
    "CHEUI_m5C": (_raw_cheui, 12, "m5C"),
    "DENA": (_raw_dena, 12, "m6A"),
    "MINES": (_raw_mines, 12, "m6A"),
    "m6Anet": (_raw_m6anet, 12, "m6A"),
    "DRUMMER": (_raw_drummer, 12, "m6A"),
    "xPore": (_raw_xpore, 12, "m6A"),
    "Nanocompore": (_raw_nanocompore, None, "m6A"),
    "NanoMUD_psi": (_raw_nanomud, 12, "psi"),
    "NanoMUD_m1psi": (_raw_nanomud, 12, "m1psi"),
}


# --------------------------------------------------------------------------- #
# genomic direct readers -> (path resolver, key_col, pos_col, offsets,
#                            score candidates [(col_idx_or_name, factor)],
#                            extras {out: src})
# --------------------------------------------------------------------------- #
def _src_dorado(source_file, sample, mod):
    p = Path(_resolve_source(source_file, "Dorado", sample, mod))
    return [p] if p.exists() else []


def _src_eligos2(source_file, sample, mod):
    for root in ("ELIGOS2_solo", "ELIGOS2_diff"):
        d = RESULT_RNA002 / root / sample
        if d.is_dir():
            hits = sorted(d.glob("*_baseExt0.txt"))
            if hits:
                return hits
    return []


def _src_yanocomp(source_file, sample, mod):
    d = RESULT_RNA002 / "yanocomp" / sample
    hits = sorted(d.glob("*_output.bed")) if d.is_dir() else []
    if not hits:
        d = (_XB / "raw/converted_callsets/yanocomp") / sample
        hits = sorted(d.glob("*_output.bed")) if d.is_dir() else []
    return hits


def _src_epinano(source_file, sample, mod):
    d = RESULT_RNA002 / "EpiNano_DiffErr" / sample
    return sorted(d.glob("*prediction.csv")) if d.is_dir() else []


DORADO_EXTRAS = {"coverage": "valid_coverage", "counts_mod": "count_modified",
                 "counts_canonical": "count_canonical",
                 "counts_other_mod": "count_other_mode",
                 "counts_delete": "count_delete", "counts_fail": "count_fail",
                 "counts_diff": "count_diff", "counts_nocall": "count_nocall"}

GENOMIC_SPECS = {
    "Dorado_": (_src_dorado, "chrom", "chromStart", (0,), None, DORADO_EXTRAS, True),
    "ELIGOS2_diff": (_src_eligos2, "chrom", "start_loc", (-1, 0), None,
                     {"kmer": "kmer5", "kmer7": "kmer7", "odds_ratio": "oddR",
                      "pvalue": "pval", "adj_pvalue": "adjPval",
                      "counts_total": "total_reads", "esb_test": "ESB_test",
                      "esb_ctrl": "ESB_ctrl"}, False),
    "ELIGOS2_solo": (_src_eligos2, "chrom", "start_loc", (-1, 0), None,
                     {"kmer": "kmer5", "kmer7": "kmer7", "odds_ratio": "oddR",
                      "pvalue": "pval", "adj_pvalue": "adjPval",
                      "counts_total": "total_reads", "esb_test": "ESB_test",
                      "esb_ctrl": "ESB_ctrl"}, False),
    "yanocomp": (_src_yanocomp, 0, 1, (-2, -1, 0, 1, 2), None,
                 {"kmer": "_kmer", "counts_total": 4, "logodds": 6,
                  "pvalue": 7, "pvalue_ttest": 8, "odds_ratio": 9,
                  "stat_11": 10, "stat_12": 11, "stat_13": 12, "stat_14": 13,
                  "stat_15": 14, "stat_16": 15, "stat_17": 16, "stat_18": 17,
                  "stat_19": 18}, False),
    "EpiNano_Error": (_src_epinano, 0, 1, (-1, 0), None,
                      {"sum_err_ko": "ko_sum_err", "sum_err_wt": "wt_sum_err",
                       "lm_z_score": "lm_residuals_z_score"}, True),
}


def _genomic_key(tool: str):
    for prefix, spec in GENOMIC_SPECS.items():
        if tool.startswith(prefix.rstrip("_")) if prefix.endswith("_") else tool == prefix:
            return prefix, spec
    return None, None


def _join_genomic(callset: pd.DataFrame, tool: str, sample: str, source_file: str):
    """Direct genomic join for tools that publish genome coordinates.

    The (file, offset, score-column, factor) combination is selected by score
    agreement on a *sample* of the callset rows, then a single full merge
    attaches the extras -- trying every candidate on the full frame was the
    reason the first enrichment run needed >1 h.
    """
    prefix, (srcfn, key_col, pos_col, offsets, _scores, extras, headered) = _genomic_key(tool)
    extra = pd.DataFrame(index=callset.index, columns=EXTRA_COLUMNS, dtype=object)
    files = srcfn(source_file, sample, "m6A")
    if not files:
        return extra, "no_raw"
    cs = callset[["chrom", "start", "score"]].copy()
    cs["start"] = _num(cs["start"]).astype("Int64")
    sample_cs = cs.dropna().sample(n=min(4000, len(cs)), random_state=0) \
        if len(cs) > 4000 else cs.dropna()

    def _prep(fp: Path):
        raw = _read(fp, sep="," if fp.suffix == ".csv" else "\t",
                    header="infer" if headered else None)
        if raw is None or raw.empty:
            return None
        if not headered and tool.startswith("yanocomp"):
            raw["_kmer"] = raw[3].astype(str).str.split(":").str[-1]
        k = raw.columns[key_col] if key_col not in raw.columns else key_col
        p = raw.columns[pos_col] if pos_col not in raw.columns else pos_col
        keep = {c for c in extras.values() if c in raw.columns}
        cand = (["percent_modified"] if tool.startswith("Dorado")
                else [c for c in raw.columns
                      if c not in (k, p) and c not in keep and
                      pd.to_numeric(raw[c], errors="coerce").notna().mean() > 0.8])
        rr = pd.DataFrame({"rchrom": _keys(raw[k]),
                           "rpos": _num(raw[p]).astype("Int64")})
        for c in keep:
            rr[c] = raw[c].to_numpy()
        for c in cand:
            rr[c] = _num(rr[c]) if c in keep else _num(raw[c])
        return rr, cand

    best = None  # (agree, n, file, off, col, fac)
    for fp in files:
        prepped = _prep(fp)
        if prepped is None:
            continue
        rr, cand = prepped
        for off in offsets:
            r = rr.copy()
            r["rpos"] = r["rpos"] + off
            j = sample_cs.merge(r, left_on=["chrom", "start"],
                                right_on=["rchrom", "rpos"], how="inner")
            if j.empty or len(j) < MIN_SCORED:
                continue
            for sc in cand:
                for fac in (1.0, 100.0, 0.01):
                    same = _close(j["score"], j[sc] * fac)
                    n = int(same.sum())
                    if n < MIN_SCORED:
                        continue
                    agree = n / len(j)
                    if best is None or (agree, n) > (best[0], best[1]):
                        best = (agree, n, fp, off, sc, fac)
            if best and best[0] >= MIN_MATCH:
                break
        if best and best[0] >= MIN_MATCH:
            break
    if best is None or best[0] < MIN_MATCH:
        return extra, f"unverified(match={best[0]:.2f},n={best[1]})" if best else "unverified"
    agree, n, fp, off, sc, fac = best
    rr, _ = _prep(fp)
    rr = rr.copy()
    rr["rpos"] = rr["rpos"] + off
    j = cs.merge(rr, left_on=["chrom", "start"],
                 right_on=["rchrom", "rpos"], how="inner")
    j = j[_close(j["score"], j[sc] * fac)]
    if j.empty:
        return extra, f"unverified(match={agree:.2f},n={n})"
    j = j.drop_duplicates(subset=["chrom", "start"])
    j = j.set_index(j["chrom"].astype(str) + ":" + j["start"].astype(str))
    ckey = (callset["chrom"].astype(str) + ":" +
            _num(callset["start"]).astype(str))
    rev = {v: k for k, v in extras.items()}
    for out_col, src_col in rev.items():
        if src_col in j.columns:
            extra[out_col] = j[src_col].reindex(ckey).to_numpy()
    return extra, (f"genomic-joined(offset={off}, match={agree:.2f}, n={n}, "
                   f"hit={int(ckey.isin(j.index).sum())}/{len(callset)})")


# --------------------------------------------------------------------------- #
# bridge-based inference
# --------------------------------------------------------------------------- #
def _infer_offset(bridge: pd.DataFrame, raw: pd.DataFrame, key: str, pos: str,
                  score_cols: tuple, bridge_col: int | None):
    """(offset, match_rate, n_scored) -- or (None, 0, n)."""
    b = bridge.rename(columns={bridge.columns[6]: "bkey",
                               bridge.columns[7]: "bpos"})
    b["bkey"] = _keys(b["bkey"])
    r = raw.rename(columns={key: "bkey", pos: "rpos"})
    r["bkey"] = _keys(r["bkey"])
    m = b.merge(r, on="bkey", how="inner")
    if m.empty:
        return None, 0.0, 0
    deltas = (_num(m["bpos"]) - _num(m["rpos"])).astype("Int64")
    scored = None
    if score_cols:
        for bc in BRIDGE_SCORE_COLS if bridge_col is None else (bridge_col,):
            if bc >= len(bridge.columns):
                continue
            bc_name = bridge.columns[bc]
            if bc_name not in m.columns:
                continue
            for sc in score_cols:
                if sc not in m.columns:
                    continue
                same = _close(m[bc_name], m[sc])
                if scored is None or same.sum() > scored.sum():
                    scored = pd.Series(same, index=m.index)
    if scored is None or scored.sum() < MIN_SCORED:
        return None, 0.0, int(scored.sum()) if scored is not None else 0
    counts = deltas[scored].value_counts()
    best = int(counts.index[0])
    rate = float(counts.iloc[0] / max(scored.sum(), 1))
    return (best, rate, int(scored.sum())) if rate >= MIN_MATCH else (None, rate, int(scored.sum()))


def _bridge_callset_offset(bridge: pd.DataFrame, callset: pd.DataFrame,
                           raw: pd.DataFrame, key: str, pos: str):
    """Offset inferred by verifying the bridge's own score columns against the
    callset score at identical genomic positions (works even when the raw
    output has no score at all, e.g. Nanocompore)."""
    b = bridge.rename(columns={bridge.columns[0]: "gchrom",
                               bridge.columns[1]: "gpos",
                               bridge.columns[6]: "bkey",
                               bridge.columns[7]: "bpos"}).copy()
    keep = ["gchrom", "gpos", "bkey", "bpos"] + \
        [bridge.columns[c] for c in BRIDGE_SCORE_COLS if c < len(bridge.columns)]
    b = b[keep]  # the genomic block's own name/score/strand would collide with
                 # the callset's columns during the merge
    b["gchrom"] = _keys(b["gchrom"])
    b["bkey"] = _keys(b["bkey"])
    b["gpos"] = _num(b["gpos"])
    j = callset.merge(b, left_on=["chrom", "start"],
                      right_on=["gchrom", "gpos"], how="inner")
    if j.empty:
        return None, 0.0, 0
    scored = None
    for bc in BRIDGE_SCORE_COLS:
        if bc >= len(bridge.columns):
            continue
        bc_name = bridge.columns[bc]
        if bc_name in j.columns:
            same = _close(j["score"], j[bc_name])
            if scored is None or same.sum() > scored.sum():
                scored = pd.Series(same.to_numpy(), index=j.index)
    if scored is None or scored.sum() < MIN_SCORED:
        return None, 0.0, int(scored.sum()) if scored is not None else 0
    js = j[scored.to_numpy()]
    # map the verified genomic keys to (bkey, bpos) and diff against the raw
    r = raw.rename(columns={key: "bkey", pos: "rpos"})
    m = js.merge(r[["bkey", "rpos"]].drop_duplicates("bkey"), on="bkey",
                 how="inner")
    if m.empty:
        return None, 0.0, int(scored.sum())
    deltas = (_num(m["bpos"]) - _num(m["rpos"])).astype("Int64")
    counts = deltas.value_counts()
    best = int(counts.index[0])
    rate = float(counts.iloc[0] / len(m))
    return (best, rate, len(m)) if rate >= MIN_MATCH else (None, rate, len(m))


def enrich(callset: pd.DataFrame, tool: str, sample: str,
           source_file: str) -> tuple[pd.DataFrame, str]:
    """Attach raw per-site columns to one callset frame (frame, note)."""
    extra = pd.DataFrame(index=callset.index, columns=EXTRA_COLUMNS, dtype=object)
    callset = callset.copy()
    callset["chrom"] = _keys(callset["chrom"])
    if tool in SLOW_TOOLS:
        return extra, "skipped_slow_raw"
    if _genomic_key(tool)[0] is not None:
        try:
            return _join_genomic(callset, tool, sample, source_file)
        except Exception as exc:  # noqa: BLE001 - enrichment must never break export
            return extra, f"raw_error({exc})"
    spec = SPECS.get(tool)
    if spec is None:
        return extra, ""
    reader, bridge_col, mod = spec
    source_file = _resolve_source(source_file, tool, sample, mod)
    bridge_p = _bridge_path(source_file, tool, sample, mod)
    if bridge_p is None:
        return extra, "no_bridge"
    try:
        raw, key, pos, score_cols, extras = reader(sample, mod)
    except Exception as exc:  # noqa: BLE001
        return extra, f"raw_error({exc})"
    if raw is None or raw.empty:
        return extra, "no_raw"
    bridge = _read(bridge_p)
    if bridge.shape[1] < 11:
        return extra, "bad_bridge"

    offset, rate, n = _bridge_callset_offset(bridge, callset, raw, key, pos)
    if offset is None:
        offset, rate, n = _infer_offset(bridge, raw, key, pos, score_cols, bridge_col)
    if offset is None:
        return extra, f"unverified(match={rate:.2f},n={n})"

    b = bridge.rename(columns={bridge.columns[0]: "gchrom",
                               bridge.columns[1]: "gpos",
                               bridge.columns[6]: "bkey",
                               bridge.columns[7]: "bpos"})
    b["gchrom"] = _keys(b["gchrom"])
    b["bkey"] = _keys(b["bkey"])
    b["gpos"] = _num(b["gpos"]).astype("Int64")
    b["bpos"] = _num(b["bpos"]).astype("Int64")
    b["rpos"] = b["bpos"] - offset
    r = raw.rename(columns={key: "bkey", pos: "rpos"})
    r["bkey"] = _keys(r["bkey"])
    r["rpos"] = _num(r["rpos"]).astype("Int64")
    b = b.merge(r.drop(columns=["rpos"]), on="bkey", how="left")
    b = b.drop_duplicates(subset=["gchrom", "gpos"])
    b["gkey"] = b["gchrom"].astype(str) + ":" + b["gpos"].astype(str)
    b = b.set_index("gkey")
    ckey = (callset["chrom"].astype(str) + ":" + callset["start"].astype(str))
    for out_col, src_col in extras.items():
        if src_col in b.columns:
            extra[out_col] = b[src_col].reindex(ckey).to_numpy()
    hit = int(ckey.isin(b.index).sum())
    return extra, f"raw-joined(offset={offset}, match={rate:.2f}, n={n}, hit={hit}/{len(callset)})"
