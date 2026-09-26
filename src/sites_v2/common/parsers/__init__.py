"""Per-tool parsers: raw tool output -> canonical callset rows.

Canonical row schema (all parsers return these columns):

``chrom pos_raw end strand status score score_type mod_ratio frac_diff
src_coverage mod_label note``

* ``pos_raw`` is the 0-based coordinate *as reported by the tool* after the
  documented base conversion only.  Offsets are never silently corrected; the
  annotation step adds ``pos_center``/``dist_center``/``dist_drach_a`` and the
  offset audit reports systematic shifts per tool.
* ``end`` is always ``pos_raw + 1`` (the differr zero-width bug is fixed here).
* ``score``/``score_type`` carry the tool's primary confidence value in its own
  units (probability, p-value, FDR, delta error ...), never re-scaled.
* ``mod_ratio`` is only set when the tool reports an actual modification ratio.

Thresholds that the legacy conversion notebooks applied upstream are applied
here for raw formats (and recorded in the callset manifest), so that the
rebuilt callsets are comparable with the published ones.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..config import MODKIT_CODE
from ..io_utils import read_table
from ..match import fix_chromosome

CANONICAL_COLUMNS = [
    "chrom", "chrom_src", "pos_raw", "end", "strand", "status", "score",
    "score_type", "mod_ratio", "frac_diff", "src_coverage", "mod_label", "note",
]

#: recognized score columns, in priority order.
_SCORE_CANDIDATES = ["Prob", "Pvalue", "Padj", "FDR", "delta_sum_err", "score", "pval"]
_MOD_RATIO_CANDIDATES = ["mod_ratio", "Mod_rate", "m6a_ratio", "ratio"]
_FRAC_DIFF_CANDIDATES = ["FracDiff", "diff_mod_rate", "frac_diff"]

_STRAND_MAP = {"+": "+", "-": "-", "1": "+", "-1": "-", "*": "*", ".": "*",
               "": "*", "nan": "*"}


def _find_col(df: pd.DataFrame, names: list[str]) -> str | None:
    lower = {c.strip().lower(): c for c in df.columns}
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    return None


def _empty() -> pd.DataFrame:
    return pd.DataFrame(columns=CANONICAL_COLUMNS)


def _col_or_default(df: pd.DataFrame, name: str, default):
    """Column as a Series, or a constant-filled Series when the column is absent."""
    if name in df.columns:
        return df[name]
    return pd.Series([default] * len(df), index=df.index)


def _finalize(df: pd.DataFrame, *, note: str = "") -> pd.DataFrame:
    """Normalise dtypes / chromosome names and enforce the canonical schema."""
    if df.empty:
        return _empty()
    work = pd.DataFrame(index=df.index)
    work["chrom"] = _col_or_default(df, "chrom", "").astype(str)
    work["pos_raw"] = pd.to_numeric(_col_or_default(df, "pos_raw", np.nan), errors="coerce")
    work = work[work["pos_raw"].notna()]
    if work.empty:
        return _empty()
    out = pd.DataFrame(index=work.index)
    out["chrom_src"] = work["chrom"].str.strip()
    out["chrom"] = [fix_chromosome(c) for c in out["chrom_src"]]
    out["pos_raw"] = work["pos_raw"].astype(np.int64)
    out["end"] = out["pos_raw"].astype(np.int64) + 1
    strand = _col_or_default(df, "strand", "*").astype(str).str.strip()
    out["strand"] = [(_STRAND_MAP.get(s, s)) for s in strand.reindex(work.index)]
    out["status"] = _col_or_default(df, "status", "Mod").astype(str).reindex(work.index)
    out["score"] = pd.to_numeric(_col_or_default(df, "score", np.nan),
                                 errors="coerce").reindex(work.index)
    out["score_type"] = _col_or_default(df, "score_type", "").astype(str).reindex(work.index)
    out["mod_ratio"] = pd.to_numeric(_col_or_default(df, "mod_ratio", np.nan),
                                     errors="coerce").reindex(work.index)
    out["frac_diff"] = pd.to_numeric(_col_or_default(df, "frac_diff", np.nan),
                                     errors="coerce").reindex(work.index)
    out["src_coverage"] = pd.to_numeric(_col_or_default(df, "src_coverage", np.nan),
                                        errors="coerce").reindex(work.index)
    out["mod_label"] = _col_or_default(df, "mod_label", "").astype(str).reindex(work.index)
    out["note"] = note
    return out[CANONICAL_COLUMNS].reset_index(drop=True)


# --------------------------------------------------------------------------- #
# standard 6-7 column callsets
# --------------------------------------------------------------------------- #
def parse_standard(path: Path, **kw) -> pd.DataFrame:
    """``Chr Start End Status <score> Strand [mod_ratio]`` files.

    ``End`` may equal ``Start`` (differr); it is always rewritten to
    ``pos_raw + 1``.
    """
    df = read_table(path)
    if df.empty:
        return _empty()
    c_chr = _find_col(df, ["Chr", "chrom", "chromosome"])
    c_start = _find_col(df, ["Start", "start", "chromStart"])
    if c_chr is None or c_start is None:
        return _empty()
    c_end = _find_col(df, ["End", "end", "chromEnd"])
    c_strand = _find_col(df, ["Strand", "strand"])
    c_status = _find_col(df, ["Status", "status"])
    c_score = _find_col(df, _SCORE_CANDIDATES)
    c_ratio = _find_col(df, _MOD_RATIO_CANDIDATES)
    c_frac = _find_col(df, _FRAC_DIFF_CANDIDATES)

    out = pd.DataFrame(index=df.index)
    out["chrom"] = df[c_chr]
    out["pos_raw"] = pd.to_numeric(df[c_start], errors="coerce")
    if c_end is not None:
        end = pd.to_numeric(df[c_end], errors="coerce")
        # repair the zero-width differr convention
        bad = end == out["pos_raw"]
        if bad.any():
            out.loc[bad, "pos_raw"] = end[bad] - 1
    out["strand"] = df[c_strand] if c_strand else "*"
    out["status"] = df[c_status] if c_status else "Mod"
    if c_score is not None:
        out["score"] = pd.to_numeric(df[c_score], errors="coerce")
        out["score_type"] = c_score
    if c_ratio is not None:
        out["mod_ratio"] = pd.to_numeric(df[c_ratio], errors="coerce")
        if c_score is None:
            out["score"] = out["mod_ratio"]
            out["score_type"] = c_ratio
    if c_frac is not None:
        out["frac_diff"] = pd.to_numeric(df[c_frac], errors="coerce")
    return _finalize(out)


# --------------------------------------------------------------------------- #
# liftover junction tables (``*_liftover.txt``)
# --------------------------------------------------------------------------- #
def parse_liftover(path: Path, **kw) -> pd.DataFrame:
    """R2Dtool junction table: ``chromosome start end name score strand Chr Start
    End Status <score> Strand [extra]``.

    Column 1-6 are the **genomic** block (0-based, verified against the legacy
    ``_remove_chr.txt`` files: DRUMMER chr4:8572415, m6Anet chr1:633838), columns
    7-13 repeat the transcript coordinate and the tool score.

    IMPORTANT: the genomic block is lowercase (``start``) while the transcript
    block is capitalised (``Start``); a case-insensitive column lookup silently
    picks the *transcript* coordinate, so the columns are addressed by position.
    """
    df = read_table(path)
    if df.empty or df.shape[1] < 6:
        return _empty()
    out = pd.DataFrame(index=df.index)
    out["chrom"] = df.iloc[:, 0]
    out["pos_raw"] = pd.to_numeric(df.iloc[:, 1], errors="coerce")
    out["strand"] = df.iloc[:, 5]
    out["status"] = df["Status"] if "Status" in df.columns else "Mod"
    if df.shape[1] >= 11:
        out["score"] = pd.to_numeric(df.iloc[:, 10], errors="coerce")
        out["score_type"] = str(df.columns[10])
    if df.shape[1] >= 12:
        out["mod_ratio"] = pd.NA
    if df.shape[1] >= 13:
        extra_col = str(df.columns[12])
        vals = pd.to_numeric(df.iloc[:, 12], errors="coerce")
        if extra_col in ("mod_ratio", "Mod_rate", "m6a_ratio", "modRatio"):
            out["mod_ratio"] = vals
        else:
            out["frac_diff"] = vals
    return _finalize(out, note="liftover junction table (genomic block = cols 1-6)")


# --------------------------------------------------------------------------- #
# Nanom6A raw ``ratio.0.5.tsv``
# --------------------------------------------------------------------------- #
def parse_nanom6a_ratio_tsv(path: Path, *, support: int = 20, min_ratio: float = 0.1,
                            **kw) -> pd.DataFrame:
    """Header-less ``<gene>|<chrom>  <start>|<mod>|<total>|<ratio> ...`` file.

    Applies the documented notebook thresholds ``total >= 20 & ratio > 0.1``.
    """
    rows = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            head = fields[0].split("|")
            if len(head) < 2:
                continue
            chrom = head[1]
            for cell in fields[1:]:
                parts = cell.split("|")
                if len(parts) < 4:
                    continue
                try:
                    start, mod, total, ratio = (int(float(parts[0])), int(float(parts[1])),
                                                int(float(parts[2])), float(parts[3]))
                except ValueError:
                    continue
                if total < support or ratio <= min_ratio:
                    continue
                rows.append({"chrom": chrom, "pos_raw": start - 1,
                             "mod_ratio": ratio, "src_coverage": total})
    if not rows:
        return _empty()
    df = pd.DataFrame(rows)
    df["score"] = df["mod_ratio"]
    df["score_type"] = "mod_ratio"
    return _finalize(df, note=f"filter: total>={support} & ratio>{min_ratio}")


# --------------------------------------------------------------------------- #
# ELIGOS2 ``*_combine.txt`` (genomic, no liftover needed)
# --------------------------------------------------------------------------- #
def parse_eligos2(path: Path, *, adjp_max: float = 1e-4, oddr_min: float = 1.2,
                  min_reads: int = 20, **kw) -> pd.DataFrame:
    """Legacy ELIGOS2 filter: ref=='A' & total_reads>20 & adjPval<1e-4 & oddR>1.2."""
    df = read_table(path)
    if df.empty or "start_loc" not in df.columns:
        return _empty()
    adjp = pd.to_numeric(df["adjPval"], errors="coerce")
    oddr = pd.to_numeric(df["oddR"], errors="coerce")
    reads = pd.to_numeric(df["total_reads"], errors="coerce")
    keep = (df["ref"].astype(str) == "A") & (reads > min_reads) & \
           (adjp < adjp_max) & (oddr > oddr_min)
    sub = df[keep]
    if sub.empty:
        return _empty()
    out = pd.DataFrame(index=sub.index)
    out["chrom"] = sub["chrom"]
    out["pos_raw"] = pd.to_numeric(sub["start_loc"], errors="coerce")
    out["strand"] = sub["strand"]
    out["status"] = "Mod"
    out["score"] = adjp[keep]
    out["score_type"] = "adjPval"
    return _finalize(out, note=f"filter: ref==A & reads>{min_reads} & adjPval<{adjp_max} & oddR>{oddr_min}")


# --------------------------------------------------------------------------- #
# yanocomp ``*_output.bed`` (genomic 5-mer spans)
# --------------------------------------------------------------------------- #
def parse_yanocomp_bed(path: Path, **kw) -> pd.DataFrame:
    """19-column yanocomp bed; keep 5-mers whose middle base is A.

    Legacy convention: ``Start = bed_start + 2``, ``End = bed_start + 3`` (the
    central base of the 5-mer), score = column 9 (FDR).
    """
    df = read_table(path, header=None)
    if df.empty or df.shape[1] < 9:
        return _empty()
    kmer = df[3].astype(str).str.split(":").str[-1]
    keep = kmer.str.len().eq(5) & kmer.str[2].eq("A")
    sub = df[keep]
    if sub.empty:
        return _empty()
    start = pd.to_numeric(sub[1], errors="coerce")
    out = pd.DataFrame(index=sub.index)
    out["chrom"] = sub[0]
    out["pos_raw"] = start + 2
    out["strand"] = sub[5]
    out["status"] = "Mod"
    out["score"] = pd.to_numeric(sub[8], errors="coerce")
    out["score_type"] = "FDR"
    return _finalize(out, note="output.bed fallback (middle base of the 5-mer)")


# --------------------------------------------------------------------------- #
# mAFiA ``out_dir/mAFiA.sites.bed``
# --------------------------------------------------------------------------- #
def parse_mafia_bed(path: Path, **kw) -> pd.DataFrame:
    df = read_table(path)
    if df.empty:
        return _empty()
    c_chr = _find_col(df, ["chrom", "Chr"])
    c_start = _find_col(df, ["chromStart", "Start", "start"])
    c_strand = _find_col(df, ["strand", "Strand"])
    c_cov = _find_col(df, ["coverage", "depth"])
    c_ratio = _find_col(df, ["modRatio", "mod_ratio"])
    c_five = _find_col(df, ["ref5mer", "5mer"])
    out = pd.DataFrame(index=df.index)
    out["chrom"] = df[c_chr]
    out["pos_raw"] = pd.to_numeric(df[c_start], errors="coerce")
    out["strand"] = df[c_strand] if c_strand else "*"
    out["status"] = "Mod"
    if c_ratio is not None:
        out["mod_ratio"] = pd.to_numeric(df[c_ratio], errors="coerce")
        out["score"] = out["mod_ratio"]
        out["score_type"] = "modRatio(count)"
    if c_cov is not None:
        out["src_coverage"] = pd.to_numeric(df[c_cov], errors="coerce")
    if c_five is not None:
        out["note"] = ""
    res = _finalize(out, note="mAFiA bed; modRatio is a COUNT")
    if c_five is not None:
        res["note"] = res["note"] + "|ref5mer=" + df[c_five].astype(str).values[:len(res)]
    return res


# --------------------------------------------------------------------------- #
# RNA004 dorado / modkit pileup bed
# --------------------------------------------------------------------------- #
def parse_dorado_pileup(path: Path, *, mod_code: str | None = None,
                        min_cov: int | None = None, min_pct: float | None = None,
                        **kw) -> pd.DataFrame:
    """modkit pileup bed; rows are split by the ``name`` modification code.

    ``mod_code`` keeps only the rows carrying that modkit code.  It is set by the
    registry when one pileup file mixes several modifications (the Curlcake
    ``*_all_pileup.bed`` carries ``17802``/``a``/``m``), so that each callset
    holds exactly one modification -- previously every row of such a file went
    into a single "m6A" callset.

    Call filter (2026-09-18)
    ------------------------
    ``modkit pileup`` emits one row per *covered position* x *modification code*,
    **including positions where nothing was called** (``percent_modified = 0``).
    The HeLa RNA004 family is already filtered at source (every row has
    ``valid_coverage >= 20`` and ``percent_modified > 0``), but the Curlcake family
    (``RNA004_result/dorado_model/*_pileup.bed``) is not: there the no-call rows are
    the overwhelming majority, so the "callset" collapsed onto the construct's base
    composition (m6A/m5C/Psi all looked ~58 % expected-base, identical to a random
    position set).  ``min_cov``/``min_pct`` are declared per source in
    ``config.RNA004_SOURCES`` and keep rows with ``valid_coverage >= min_cov`` and
    ``percent_modified > min_pct`` (``None`` = that filter is off, so every source
    that does not declare them keeps its previous behaviour).  ``score`` /
    ``mod_ratio`` stay ``percent_modified/100`` so the downstream threshold scan
    (5/10/20/50 %) is unaffected.
    """
    df = read_table(path)
    if df.empty:
        return _empty()
    c_chr = _find_col(df, ["chrom", "#chrom", "chromosome"])
    c_start = _find_col(df, ["chromStart", "start"])
    c_end = _find_col(df, ["chromEnd", "end"])
    c_name = _find_col(df, ["name", "mod"])
    c_score = _find_col(df, ["score"])
    c_strand = _find_col(df, ["strand"])
    c_cov = _find_col(df, ["valid_coverage", "coverage"])
    c_pct = _find_col(df, ["percent_modified", "pct_modified"])
    if c_chr is None or c_start is None:
        return _empty()
    n_in = len(df)
    filter_note = ""
    if c_cov is not None and min_cov is not None:
        keep_cov = pd.to_numeric(df[c_cov], errors="coerce").fillna(0) >= min_cov
        n_before = len(df)
        df = df[keep_cov]
        filter_note += f"|cov>={min_cov} dropped {n_before - len(df)}/{n_in}"
    if c_pct is not None and min_pct is not None and not df.empty:
        keep_pct = pd.to_numeric(df[c_pct], errors="coerce").fillna(0) > min_pct
        n_before = len(df)
        df = df[keep_pct]
        filter_note += f"|pct>{min_pct} dropped {n_before - len(df)}/{n_in}"
    if df.empty:
        # every row was filtered out: a real 0-detection result (header-only callset)
        return _empty()
    name = df[c_name].astype(str) if c_name is not None else pd.Series("NA", index=df.index)
    code = name.str.split("#").str[0]
    if mod_code is not None:
        keep = code == mod_code
        df, code = df[keep], code[keep]
        if df.empty:
            # the file does not carry this modification after all
            return _empty()
    out = pd.DataFrame(index=df.index)
    out["chrom"] = df[c_chr]
    out["pos_raw"] = pd.to_numeric(df[c_start], errors="coerce")
    out["strand"] = df[c_strand] if c_strand else "*"
    out["status"] = "Mod"
    if c_pct is not None:
        out["score"] = pd.to_numeric(df[c_pct], errors="coerce") / 100.0
        out["score_type"] = "percent_modified/100"
        out["mod_ratio"] = out["score"]
    elif c_score is not None:
        out["score"] = pd.to_numeric(df[c_score], errors="coerce")
        out["score_type"] = "score"
    if c_cov is not None:
        out["src_coverage"] = pd.to_numeric(df[c_cov], errors="coerce")
    out["mod_label"] = code.map(lambda c: MODKIT_CODE.get(c, c))
    return _finalize(out, note=f"source={path.name}{filter_note}")


# --------------------------------------------------------------------------- #
# m6Anet site_proba csv (RNA004 Curlcake; RNA002 raw)
# --------------------------------------------------------------------------- #
def parse_m6anet_csv(path: Path, *, platform: str = "RNA004", prob: float = 0.5,
                     min_reads: int = 20, min_ratio: float = 0.1, **kw) -> pd.DataFrame:
    df = read_table(path, sep=",")
    if df.empty:
        return _empty()
    c_chr = _find_col(df, ["transcript_id", "transcript_ID", "chrom"])
    c_pos = _find_col(df, ["transcript_position", "position"])
    c_reads = _find_col(df, ["n_reads", "coverage"])
    c_prob = _find_col(df, ["probability_modified", "prob_modified"])
    c_ratio = _find_col(df, ["mod_ratio"])
    out = pd.DataFrame(index=df.index)
    out["chrom"] = df[c_chr].astype(str).str.replace(r"_[FR]$", "", regex=True)
    out["pos_raw"] = pd.to_numeric(df[c_pos], errors="coerce") - 1  # 1-based input
    out["strand"] = "*"
    out["status"] = "Mod"
    out["score"] = pd.to_numeric(df[c_prob], errors="coerce")
    out["score_type"] = c_prob
    if c_ratio is not None:
        out["mod_ratio"] = pd.to_numeric(df[c_ratio], errors="coerce")
    if c_reads is not None:
        out["src_coverage"] = pd.to_numeric(df[c_reads], errors="coerce")

    mask = out["score"] >= prob
    reads = out["src_coverage"]
    if platform == "RNA002" and c_ratio is not None:
        mask &= out["mod_ratio"].fillna(0) > min_ratio
    elif reads.notna().any():
        mask &= reads.fillna(0) >= min_reads
    out = out[mask.fillna(False)]
    note = f"filter: prob>={prob}"
    note += f" & mod_ratio>{min_ratio}" if platform == "RNA002" else f" & n_reads>={min_reads}"
    return _finalize(out, note=note)


# --------------------------------------------------------------------------- #
# NanoSPA / NanoPsu csv (Curlcake RNA004)
# --------------------------------------------------------------------------- #
def parse_nanospa_csv(path: Path, *, prob: float = 0.5, **kw) -> pd.DataFrame:
    """Header-less ``transcript[_F/_R], position, base, coverage, prob`` csv."""
    df = read_table(path, sep=",", header=None)
    if df.empty:
        return _empty()
    if df.shape[1] < 5:
        return _empty()
    df = df.iloc[:, :5]
    df.columns = ["tid", "position", "base", "coverage", "prob"]
    out = pd.DataFrame(index=df.index)
    tid = df["tid"].astype(str)
    out["chrom"] = tid.str.replace(r"_[FR]$", "", regex=True)
    out["strand"] = np.where(tid.str.endswith("_F"), "+",
                             np.where(tid.str.endswith("_R"), "-", "*"))
    out["pos_raw"] = pd.to_numeric(df["position"], errors="coerce") - 1
    out["status"] = "Mod"
    out["score"] = pd.to_numeric(df["prob"], errors="coerce")
    out["score_type"] = "prob"
    out["src_coverage"] = pd.to_numeric(df["coverage"], errors="coerce")
    out = out[out["score"] >= prob]
    return _finalize(out, note=f"filter: prob>={prob}")


# --------------------------------------------------------------------------- #
# dispatch
# --------------------------------------------------------------------------- #
PARSERS = {
    "standard": parse_standard,
    "liftover": parse_liftover,
    "eligos2": parse_eligos2,
    "yanocomp_bed": parse_yanocomp_bed,
    "eligos2": parse_eligos2,
    "yanocomp_bed": parse_yanocomp_bed,
    "nanom6a_ratio_tsv": parse_nanom6a_ratio_tsv,
    "mafia_bed": parse_mafia_bed,
    "dorado_pileup": parse_dorado_pileup,
    "m6anet_csv": parse_m6anet_csv,
    "nanospa_csv": parse_nanospa_csv,
}


def parse_callset(path: Path, parser: str, **kw) -> pd.DataFrame:
    if parser not in PARSERS:
        raise KeyError(f"unknown parser '{parser}'")
    return PARSERS[parser](Path(path), **kw)
