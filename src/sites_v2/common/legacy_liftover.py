"""Transcript-space raw  ->  legacy-style liftover input  ->  genomic callsets.

Several tools (CHEUI, DENA, MINES, m6Anet, DRUMMER, Tombo_com) emit
**transcript coordinates**; the legacy pipeline converted them with R2Dtool
(``r2d liftover -H -g <GTF> -i <input>``).  The converters below reproduce the
legacy preprocessing *exactly* (verified row-for-row against the existing
``*_liftover.txt`` files), so that replicates which the old pipeline never
converted can be filled in with the same conventions:

================  ==========================================================
tool              legacy input recipe (from code_user/python_postprocessing)
================  ==========================================================
CHEUI_m6A/m5C     ``site_level_<mod>_predictions.txt``; Prob>0.999 &
                  mod_ratio(stoichiometry)>0.1; **Start=pos+4, End=pos+5**
                  (legacy manual shift; recorded by the offset audit)
DENA              ``<sample>.tsv``; coverage(col3+col4)>20 & ratio>0.1;
                  Start=pos-1, End=pos
MINES             ``<sample>.bed``; coverage>=20 & mod_ratio>0.1;
                  Start/End 0-based as-is, strand from the bed
m6Anet            ``data.site_proba.csv``; prob>0.5 & mod_ratio>0.1;
                  Start=End=transcript_position+1
DRUMMER           ``summary.txt``; reference_base=='A' & depth_ctrl>20 &
                  |frac_diff|>0.1; Start=transcript_pos-1
Tombo_com         ``<sample>.txt``  (already transcript-space, pass through)
================  ==========================================================

The produced file has the legacy column layout ``Chr Start End Status <score>
Strand [mod_ratio]`` because R2Dtool reads columns 2-3 as coordinates and
column 6 as strand.
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

from .config import GENOMES, canonical_mod_type
from .io_utils import FASTA, normalize_chrom_for_fasta, read_table
from .match import fix_chromosome


def _out(chrom, start, end, status, score, strand, score_name: str = "Score",
         mod_ratio=None, extra_name: str = "mod_ratio") -> pd.DataFrame:
    """Legacy input layout: Chr Start End Status <score_name> Strand [extra]."""
    n = len(np.asarray(chrom, dtype=object))
    df = pd.DataFrame({
        "Chr": np.asarray(chrom, dtype=object).astype(str),
        "Start": np.asarray(start, dtype="int64"),
        "End": np.asarray(end, dtype="int64"),
        "Status": np.asarray(status, dtype=object).astype(str),
        score_name: np.asarray(score, dtype=float),
        "Strand": np.asarray(strand, dtype=object).astype(str),
    })
    assert len(df) == n
    if mod_ratio is not None:
        df[extra_name] = np.asarray(mod_ratio, dtype=float)
    return df[df["Chr"] != ""].reset_index(drop=True)


# --------------------------------------------------------------------------- #
# per-tool builders
# --------------------------------------------------------------------------- #
def build_cheui(sample_dir: Path, mod_type: str) -> pd.DataFrame | None:
    mod = "m6A" if mod_type == "m6A" else "m5C"
    f = sample_dir / f"site_level_{mod}_predictions.txt"
    if not f.exists():
        return None
    df = read_table(f)
    if df.empty:
        return None
    prob = pd.to_numeric(df["probability"], errors="coerce")
    stoich = pd.to_numeric(df["stoichiometry"], errors="coerce")
    keep = (prob > 0.999) & (stoich > 0.1)
    sub = df[keep]
    if sub.empty:
        return None
    pos = pd.to_numeric(sub["position"], errors="coerce").to_numpy(dtype=np.int64)
    # legacy convention: Start = position + 4, End = position + 5 (audited later)
    return _out(np.asarray(sub["contig"]), pos + 4, pos + 5, "Mod",
                prob[keep].to_numpy(), "*", score_name="Prob",
                mod_ratio=stoich[keep].to_numpy())


def build_dena(sample_dir: Path, mod_type: str) -> pd.DataFrame | None:
    hits = sorted(sample_dir.glob("*.tsv"))
    if not hits:
        return None
    df = read_table(hits[0], header=None)
    if df.empty or df.shape[1] < 6:
        return None
    cov = (pd.to_numeric(df[3], errors="coerce") + pd.to_numeric(df[4], errors="coerce"))
    ratio = pd.to_numeric(df[5], errors="coerce")
    keep = (cov > 20) & (ratio > 0.1)
    sub = df[keep]
    if sub.empty:
        return None
    pos = pd.to_numeric(sub[1], errors="coerce").to_numpy(dtype=np.int64)
    return _out(np.asarray(sub[0]), pos - 1, pos, "mod", ratio[keep].to_numpy(),
                "*", score_name="m6a_ratio")


def build_mines(sample_dir: Path, mod_type: str) -> pd.DataFrame | None:
    hits = sorted(sample_dir.glob("*.bed"))
    if not hits:
        return None
    df = read_table(hits[0], header=None)
    if df.shape[1] < 8 or df.empty:
        return None
    ratio = pd.to_numeric(df[6], errors="coerce")
    cov = pd.to_numeric(df[7], errors="coerce")
    keep = (cov >= 20) & (ratio > 0.1)
    sub = df[keep]
    if sub.empty:
        return None
    # legacy file ground truth: the MINES bed is 1-based, Start/End = columns 2-3 - 1
    return _out(np.asarray(sub[0]),
                pd.to_numeric(sub[1]).to_numpy(dtype=np.int64) - 1,
                pd.to_numeric(sub[2]).to_numpy(dtype=np.int64) - 1, "Mod",
                ratio[keep].to_numpy(), np.asarray(sub[5]), score_name="mod_ratio")


def build_m6anet(sample_dir: Path, mod_type: str) -> pd.DataFrame | None:
    f = sample_dir / "data.site_proba.csv"
    if not f.exists():
        return None
    df = read_table(f, sep=",")
    if df.empty:
        return None
    prob = pd.to_numeric(df["probability_modified"], errors="coerce")
    ratio = pd.to_numeric(df["mod_ratio"], errors="coerce")
    keep = (prob > 0.5) & (ratio > 0.1)
    sub = df[keep]
    if sub.empty:
        return None
    # legacy file ground truth: Start = transcript_position, End = position + 1
    pos = pd.to_numeric(sub["transcript_position"], errors="coerce").to_numpy(dtype=np.int64)
    return _out(np.asarray(sub["transcript_id"]), pos, pos + 1, "Mod",
                prob[keep].to_numpy(), "*", score_name="Prob",
                mod_ratio=ratio[keep].to_numpy())


def build_drummer(sample_dir: Path, mod_type: str) -> pd.DataFrame | None:
    f = sample_dir / "summary.txt"
    if not f.exists():
        return None
    df = read_table(f)
    if df.empty or "reference_base" not in df.columns:
        return None
    frac = pd.to_numeric(df["frac_diff"], errors="coerce")
    depth = pd.to_numeric(df["depth_ctrl"], errors="coerce")
    padj = pd.to_numeric(df["OR_padj"], errors="coerce")
    keep = (df["reference_base"] == "A") & (depth > 20) & (frac.abs() > 0.1)
    sub = df[keep]
    if sub.empty:
        return None
    pos = pd.to_numeric(sub["transcript_pos"], errors="coerce").to_numpy(dtype=np.int64)
    status = np.where(padj[keep].to_numpy() < 0.05, "Mod", "Unmod")
    return _out(np.asarray(sub["transcript_id"]), pos - 1, pos, status,
                padj[keep].to_numpy(), "*", score_name="Pvalue",
                mod_ratio=None if "frac_diff" not in sub else frac[keep].to_numpy(),
                extra_name="FracDiff")


def build_nanocompore(sample_dir: Path, mod_type: str,
                      species: str | None = None) -> pd.DataFrame | None:
    """``outnanocompore_results.tsv`` -> legacy input (two coordinate conventions).

    Legacy recipe (``code/python_postprocessing/*/Nanocompore.ipynb``)::

        ref_kmer[2] == 'A'
        Status = 'Mod' if GMM_logit_pvalue < 0.05 and |Logit_LOR| > 0.5   ('NC' -> NaN)
        keep only Status == 'Mod'

    * transcript reference (Arabidopsis / Mouse / Human / E. coli):
      ``Chr = ref_id``, ``Start = pos + 2``, ``End = pos + 3``  (legacy ``pos += 3``
      then ``Start -= 1``), then R2Dtool liftover and the chromosome whitelist of
      :data:`POST_LIFTOVER_KEEP`.
    * Curlcake: the reference was ``cc.fasta`` itself, so the notebook used the
      ``chr`` / ``genomicPos`` columns directly (``genomicPos += 2``, zero-width
      interval) with **no liftover**.
    """
    f = sample_dir / "outnanocompore_results.tsv"
    if not f.exists():
        return None
    df = read_table(f)
    if df.empty or "ref_kmer" not in df.columns:
        return None
    kmer = df["ref_kmer"].astype(str)
    df = df[kmer.str.len().ge(3) & kmer.str[2].eq("A")]
    if df.empty:
        return None
    pval = pd.to_numeric(df["GMM_logit_pvalue"], errors="coerce")
    lor = pd.to_numeric(df["Logit_LOR"].replace("NC", np.nan), errors="coerce")
    keep = pval.notna() & lor.notna() & (pval < 0.05) & (lor.abs() > 0.5)
    sub = df[keep]
    if sub.empty:
        return None
    construct = (species == "Curlcake") or (
        "genomicPos" in sub.columns and str(sub["ref_id"].iloc[0]).lower().startswith("curlcake"))
    if construct:
        pos = pd.to_numeric(sub["genomicPos"], errors="coerce").to_numpy(dtype=np.int64) + 2
        chrom = np.asarray(sub["chr"])
        return _out(chrom, pos, pos, "Mod", pval[keep].to_numpy(),
                    np.asarray(sub["strand"]), score_name="Pvalue")
    pos = pd.to_numeric(sub["pos"], errors="coerce").to_numpy(dtype=np.int64) + 3
    return _out(np.asarray(sub["ref_id"]), pos - 1, pos, "Mod",
                pval[keep].to_numpy(), np.asarray(sub["strand"]), score_name="Pvalue")


def build_xpore(sample_dir: Path, mod_type: str) -> pd.DataFrame | None:
    """``diffmod.table`` (transcript space) -> FDR-adjusted legacy input."""
    f = sample_dir / "diffmod.table"
    if not f.exists():
        return None
    df = read_table(f, sep=",")
    if df.empty or "kmer" not in df.columns:
        return None
    kmer = df["kmer"].astype(str)
    keep_kmer = kmer.str.len().eq(5) & kmer.str[2].eq("A")
    df = df[keep_kmer]
    if df.empty:
        return None
    from statsmodels.stats.multitest import multipletests
    cols = list(df.columns)
    p_col = [c for c in cols if c.startswith("pval")]
    d_col = [c for c in cols if c.startswith("diff_mod_rate")]
    if not p_col or not d_col:
        return None
    pvals = pd.to_numeric(df[p_col[0]], errors="coerce").to_numpy()
    ok = ~np.isnan(pvals)
    fdr = np.full(pvals.shape, np.nan)
    fdr[ok] = multipletests(pvals[ok], method="fdr_bh")[1]
    diff = pd.to_numeric(df[d_col[0]], errors="coerce").to_numpy()
    keep = (fdr < 0.05) & (np.abs(diff) > 0.1)
    if not keep.any():
        return None
    pos = pd.to_numeric(df[cols[1]], errors="coerce").to_numpy(dtype=np.int64)[keep]
    chrom = np.asarray(df[cols[0]])[keep]
    return _out(chrom, pos, pos + 1, "Mod", fdr[keep], "*",
                score_name="FDR", mod_ratio=diff[keep], extra_name="diff_mod_rate")


# --------------------------------------------------------------------------- #
# genomic tools (no liftover): NanoSPA / EpiNano_Error
#
# These emit *genomic* coordinates directly, so they skip r2d/exon liftover; the
# legacy notebook chain is reproduced verbatim (filter -> chromosome whitelist ->
# 5-mer centre filter for NanoSPA).  Both were validated row-for-row against the
# archived rep3 files (NanoSPA_m6A Arabidopsis_WT_rep3 = 1455/1455;
# EpiNano_Error Arabidopsis_WT_result3 = 5833/5833).
# --------------------------------------------------------------------------- #
#: species -> chromosome whitelist in ``fix_chromosome``-normalised form.  The
#: legacy notebooks hard-coded the whitelist in each tool's own naming -- bare
#: ``{1..22,X}`` for Arabidopsis, ``{chr1..chr22,chrX}`` for Mouse/Human,
#: ``{1..22,Chromosome}`` for E. coli -- which all collapse onto the same
#: normalised set (scaffolds / organelle contigs like ``chrUn``/``Mt`` drop out).
_EUK_WHITELIST_NORM: set[str] = {f"chr{i}" for i in range(1, 23)} | {"chrx"}
CHROMOSOME_WHITELIST: dict[str, set[str]] = {
    "Arabidopsis": _EUK_WHITELIST_NORM,
    "Mouse": _EUK_WHITELIST_NORM,
    "Human": _EUK_WHITELIST_NORM,
    "E.coli": _EUK_WHITELIST_NORM | {"chromosome"},
}

#: modification -> (plus-strand centre base, minus-strand centre base) kept by the
#: legacy ``*_5mer.txt`` filter (genome orientation, extract_5mer.py centre).
MOD_CENTRE_BASES: dict[str, tuple[str, str]] = {
    "m6A": ("A", "T"), "inosine": ("A", "T"),
    "m5C": ("C", "G"), "Psi": ("T", "A"), "m1Psi": ("T", "A"),
}


def _empty_legacy(score_name: str = "Score") -> pd.DataFrame:
    """Empty legacy-layout frame (raw present, 0 rows pass the tool threshold)."""
    return _out([], [], [], [], [], [], score_name=score_name)


def _whitelist(df: pd.DataFrame, species: str | None) -> pd.DataFrame:
    """Keep only the chromosomes the legacy ``*_remove_chr`` step kept.

    The membership test runs on ``fix_chromosome``-normalised names so the bare
    (Arabidopsis ``1``) and ``chr``-prefixed (Mouse ``chr1``) tool namings both
    work; the original ``Chr`` value is preserved so the output still matches the
    legacy ``*_remove_chr.txt`` / ``*_5mer.txt`` byte-for-byte.
    """
    wl = CHROMOSOME_WHITELIST.get(species) if species else None
    if not wl or df.empty:
        return df
    norm = df["Chr"].astype(str).map(fix_chromosome)
    return df[norm.isin(wl)].reset_index(drop=True)


_FA_NAME_CACHE: dict[tuple[str, str], str | None] = {}


def _genome_center_base(genome: Path, chrom: str, pos0: int) -> str:
    """Single reference base at 0-based ``pos0`` (extract_5mer.py centre).

    Uses the shared :data:`~common.io_utils.FASTA` cache (whole-chromosome string,
    loaded once per process -- the same mechanism ``03_annotate``/``06_eval`` use),
    so a batch fill over many Mouse samples pays the cold CephFS genome read once.
    The chromosome-name normalisation is cached per (genome, chrom).
    """
    ck = (str(genome), str(chrom))
    if ck not in _FA_NAME_CACHE:
        _FA_NAME_CACHE[ck] = normalize_chrom_for_fasta(FASTA.keys(genome), chrom)
    fa_name = _FA_NAME_CACHE[ck]
    if fa_name is None:
        return "N"
    seq = FASTA.chrom(genome, fa_name)
    if seq is None or pos0 < 0 or pos0 >= len(seq):
        return "N"
    return seq[pos0].upper()


def build_nanospa(sample_dir: Path, mod_type: str,
                  species: str | None = None) -> pd.DataFrame | None:
    """``alignment/prediction_<mod>.csv`` -> genomic callset (NanoSPA notebook).

    Recipe (``code_user/python_postprocessing/<species>/NanoSPA_m6A.ipynb`` +
    ``code/next_postprocessing/extract_5mer.py``)::

        header-less csv: chr_strand, Position(1-based), Base, Coverage, Prob
        keep Coverage > 20 & Prob > 0.5  (Status = Mod)
        Chr = chr_strand.split('_')[0]; Strand = '+' if suffix=='F' else '-'
        Start = Position - 1; End = Position
        chromosome whitelist; then the genome 5-mer CENTRE base must be
        A (on '+') or T (on '-') for m6A  -- the legacy ``*_5mer.txt`` filter.

    Returns ``None`` only when the raw csv is absent (vs an empty frame for a
    present-but-0-site csv).  Already genomic -> no liftover.
    """
    mod = canonical_mod_type(mod_type)
    suffix = "m6A" if mod == "m6A" else "psU"
    f = sample_dir / "alignment" / f"prediction_{suffix}.csv"
    if not f.exists():
        return None
    df = read_table(f, sep=",", header=None)
    if df.empty or df.shape[1] < 5:
        return _empty_legacy("Prob")
    df = df.iloc[:, :5]
    chr_strand = df[0].astype(str)
    cov = pd.to_numeric(df[3], errors="coerce")
    prob = pd.to_numeric(df[4], errors="coerce")
    keep = (cov > 20) & (prob > 0.5)
    if not bool(keep.any()):
        return _empty_legacy("Prob")
    cs = chr_strand[keep]
    chrom = cs.str.split("_").str[0].to_numpy()
    strand = np.where(cs.str.split("_").str[1].to_numpy() == "F", "+", "-")
    pos = pd.to_numeric(df[1][keep], errors="coerce").to_numpy(dtype="float")
    prob_v = prob[keep].to_numpy(dtype="float")
    ok = np.isfinite(pos)
    chrom, strand, pos, prob_v = chrom[ok], strand[ok], pos[ok], prob_v[ok]
    start = pos.astype(np.int64) - 1  # 0-based
    out = _out(chrom, start, start + 1, "Mod", prob_v, strand, score_name="Prob")
    out = _whitelist(out, species)
    # 5-mer centre filter (legacy extract_5mer + *_5mer.txt)
    genome = GENOMES.get(species) if species else None
    plus_b, minus_b = MOD_CENTRE_BASES.get(mod, ("A", "T"))
    if genome is not None and not out.empty:
        centres = np.array([_genome_center_base(genome, c, int(s))
                            for c, s in zip(out["Chr"].to_numpy(),
                                            out["Start"].to_numpy())])
        strands = out["Strand"].to_numpy()
        m = np.array([(st == "+" and b == plus_b) or (st == "-" and b == minus_b)
                      for st, b in zip(strands, centres)])
        out = out[m].reset_index(drop=True)
    return out


def _epinano_error_csvs(sample_dir: Path) -> list[Path]:
    """Aggregated ``<dir>.csv`` if present, else the fwd+rev delta-sum_err preds.

    The aggregated csv and the ``*_fwd/rev.delta-sum_err.prediction.csv`` pair have
    identical columns (``chr_pos,ko_feature,wt_feature,delta_sum_err,z_scores,
    z_score_prediction``); some rep1/rep2 dirs only ever got the per-direction
    prediction files, so they are concatenated to rebuild the aggregated input.
    """
    agg = sample_dir / f"{sample_dir.name}.csv"
    if agg.exists():
        return [agg]
    return (sorted(sample_dir.glob("*_fwd.delta-sum_err.prediction.csv"))
            + sorted(sample_dir.glob("*_rev.delta-sum_err.prediction.csv")))


def dedupe_epinano_sites(df: pd.DataFrame | None,
                         subset: tuple[str, ...] = ("Chr", "Start", "Strand")
                         ) -> pd.DataFrame | None:
    """Collapse repeated site rows of the EpiNano DiffErr chain (2026-09-18).

    The legacy ``*_fwd/rev.delta-sum_err.prediction.csv`` files were **appended to
    by repeated R runs**: the same ``chr_pos`` key appears up to 5x, with identical
    values for the key (verified: 0 keys with differing scores in the Arabidopsis
    KD rep1 and both mouse callsets).  The legacy notebook keeps every copy, so the
    published EpiNano call counts are inflated ~2-3.3x:

    * Arabidopsis KD rep1: 31 981 rows -> **9 561 sites** (rep3 13 996 -> 4 213);
      WT rep3: 5 833 -> 1 728;
    * mouse KO: 30 478 -> 21 420; mouse WT (Mettl3): 24 458 -> 17 186.

    HeLa / E.coli / Curlcake chains are clean (1 copy per site).  A callset is a
    *set* of sites, so the copies are collapsed here -- for both the rebuilt
    (``result/``) and the converted-legacy (``01_extract``) paths, because the
    mouse chains have no per-site tables left to re-run from.
    """
    if df is None or len(df) == 0:
        return df
    cols = [c for c in subset if c in df.columns]
    if not cols:
        return df
    return df.drop_duplicates(subset=cols, keep="first").reset_index(drop=True)


def build_epinano_error(sample_dir: Path, mod_type: str,
                        species: str | None = None) -> pd.DataFrame | None:
    """EpiNano DiffErr ``<sample>.csv`` -> genomic callset (EpiNano_Error notebook).

    Recipe (``code_user/python_postprocessing/<species>/EpiNano_Error.ipynb``)::

        chr_pos = "<chrom> <pos> <base> <strand>" (space separated)
        keep pos>0 & strand in {+,-} & z_score_prediction in {mod,unm} & delta_sum_err notna
        condition: z_score_prediction=='mod' & delta_sum_err>0.1 &
                   ((strand=='+' & base=='A') | (strand=='-' & base=='T'))
        Start = pos (NO -1 shift, unlike DRUMMER/NanoSPA); End = pos+1; whitelist.

    Returns ``None`` only when no raw csv is present.  Already genomic -> no liftover.
    """
    csvs = _epinano_error_csvs(sample_dir)
    if not csvs:
        return None
    frames = []
    for c in csvs:
        d = read_table(c, sep=",")
        if not d.empty and "chr_pos" in d.columns:
            frames.append(d)
    if not frames:
        return _empty_legacy("delta_sum_err")
    df = pd.concat(frames, ignore_index=True)
    split = df["chr_pos"].astype(str).str.split(r"\s+", expand=True)
    if split.shape[1] < 4:
        return _empty_legacy("delta_sum_err")
    chrom = split[0]
    posn = pd.to_numeric(split[1], errors="coerce")
    base = split[2].astype(str)
    strand = split[3].astype(str)
    dse = pd.to_numeric(df["delta_sum_err"], errors="coerce")
    zpred = df["z_score_prediction"].astype(str)
    keep = (posn.notna() & (posn > 0) & strand.isin(["+", "-"])
            & zpred.isin(["mod", "unm"]) & dse.notna()
            & (zpred == "mod") & (dse > 0.1)
            & (((strand == "+") & (base == "A")) | ((strand == "-") & (base == "T"))))
    if not bool(keep.any()):
        return _empty_legacy("delta_sum_err")
    pos_v = posn[keep].to_numpy(dtype=np.int64)   # Start = Pos (no -1 shift)
    out = _out(chrom[keep].to_numpy(), pos_v, pos_v + 1, zpred[keep].to_numpy(),
               dse[keep].to_numpy(dtype="float"), strand[keep].to_numpy(),
               score_name="delta_sum_err")
    #: collapse the multi-appended prediction rows (see dedupe_epinano_sites)
    return dedupe_epinano_sites(_whitelist(out, species))


#: tools whose builder output is already in *genomic* coordinates for every
#: species (skip the transcript->genome liftover entirely).  NOTE: only *m6A*
#: tools are raw-filled; non-m6A tools (NanoSPA_psU / NanoPsu / NanoMUD / NanoNm /
#: EpiNano_SVM / differr) are deliberately NOT registered as builders -- their
#: per-species notebooks use tool-specific thresholds that are not ported, their
#: in-scope (HeLa / Curlcake) callsets already come from the archived converted
#: files via ``01_extract``, and their Arabidopsis/Mouse/E.coli callsets are
#: deleted by ``11``.  Raw-filling them would inject wrongly-thresholded sites.
GENOMIC_TOOLS: set[str] = {"NanoSPA_m6A", "EpiNano_Error"}

#: for a genomic (no-liftover) builder, the archived final file 02b validates
#: against (NanoSPA keeps the 5-mer-filtered ``*_5mer.txt``; EpiNano the
#: chromosome-whitelisted ``*_remove_chr.txt``; Nanocompore the construct file).
GENOMIC_LEGACY_GLOB: dict[str, str] = {
    "Nanocompore": "*_Nanocompore.txt",
    "NanoSPA_m6A": "*_NanoSPA_m6A_5mer.txt",
    "EpiNano_Error": "*_EpiNano_Error_remove_chr.txt",
}


def _first_existing(sample_dir: Path, *rels: str) -> Path | None:
    for r in rels:
        p = sample_dir / r
        if p.exists():
            return p
    return None


#: tool -> probe returning the raw input path when it is present (regardless of
#: how many sites survive the tool threshold).  ``01b`` uses this to tell a real
#: "raw present, Total_Detected = 0" result apart from "raw never produced".
RAW_PROBE: dict[str, object] = {
    "DRUMMER": lambda d, m: _first_existing(d, "summary.txt"),
    "CHEUI_m6A": lambda d, m: _first_existing(d, "site_level_m6A_predictions.txt"),
    "CHEUI_m5C": lambda d, m: _first_existing(d, "site_level_m5C_predictions.txt"),
    "DENA": lambda d, m: next(iter(sorted(d.glob("*.tsv"))), None),
    "MINES": lambda d, m: next(iter(sorted(d.glob("*.bed"))), None),
    "m6Anet": lambda d, m: _first_existing(d, "data.site_proba.csv"),
    "Nanocompore": lambda d, m: _first_existing(d, "outnanocompore_results.tsv"),
    "xPore": lambda d, m: _first_existing(d, "diffmod.table"),
    "NanoSPA_m6A": lambda d, m: _first_existing(d, "alignment/prediction_m6A.csv"),
    "EpiNano_Error": lambda d, m: (next(iter(_epinano_error_csvs(d)), None)),
}


BUILDERS = {
    "CHEUI_m6A": build_cheui,
    "CHEUI_m5C": build_cheui,
    "DENA": build_dena,
    "MINES": build_mines,
    "m6Anet": build_m6anet,
    "DRUMMER": build_drummer,
    "Nanocompore": build_nanocompore,
    "xPore": build_xpore,
    "NanoSPA_m6A": build_nanospa,
    "EpiNano_Error": build_epinano_error,
}

#: builders that need to know the species (different coordinate conventions,
#: chromosome whitelists, or a genome 5-mer filter).
SPECIES_AWARE: set[str] = {"Nanocompore", "NanoSPA_m6A", "EpiNano_Error"}

#: (tool, species) pairs whose legacy callset is already in *genome* space and
#: must therefore NOT be passed through ``r2d liftover``.
#:
#: Curlcake has no transcriptome reference -- every tool ran directly on
#: ``cc.fasta``, so each tool's raw output is already in construct coordinates
#: and ``build_<tool>`` reproduces the user's notebook convention verbatim
#: (e.g. CHEUI keeps the unconditional ``Start = position + 4`` of the HeLa
#: notebook; m6Anet/DENA/MINES converted files exist for all five runs).
NO_LIFTOVER: set[tuple[str, str]] = {
    ("Nanocompore", "Curlcake"),
    ("CHEUI_m6A", "Curlcake"),
    ("CHEUI_m5C", "Curlcake"),
}

#: Chromosome whitelist applied by the legacy ``*_remove_chr`` step, per tool and
#: species.  The notebooks hard-coded ``{1..22, X}`` (E. coli: ``Chromosome``), so
#: scaffolds / organelle contigs were silently dropped - reproduce it to stay
#: identical to the published callsets.
POST_LIFTOVER_KEEP: dict[str, dict[str, set[str]]] = {
    "Nanocompore": {
        "Arabidopsis": {str(i) for i in range(1, 23)} | {"X"},
        "Mouse": {str(i) for i in range(1, 23)} | {"X"},
        "Human": {str(i) for i in range(1, 23)} | {"X"},
        "E.coli": {str(i) for i in range(1, 23)} | {"Chromosome"},
    },
}

#: Raw directories whose name cannot be canonicalised to a sample, mapped to the
#: sample they belong to (mostly the Curlcake Nanocompore comparisons, where the
#: directory is named after the run index rather than the sample).
EXPLICIT_DIR_SAMPLE: dict[str, dict[str, str]] = {
    "Nanocompore": {
        # ``nanocompore sampcomp`` writes its output into ONE directory for the
        # E.coli IVT ctrl-vs-treat pair (label1=E_IVT_neg1, label2=E_IVT_neg2);
        # for every other pair the legacy output directory is named after
        # ``label2`` (the treated sample), so this comparison belongs to
        # E_IVT_neg2.
        # The four Curlcake comparisons used to need an entry here
        # (``nanocompore_result1..4``); their directories are now named
        # ``<test>_vs_<ctrl>`` and resolve through the alias index instead.
        "E_IVT_neg": "E_IVT_neg2",
    },
    # Group-level comparison directories of the E.coli IVT (negative control)
    # pair.  Each tool ran ctrl=E_IVT_neg1 / test=E_IVT_neg2 in a single output
    # directory; the legacy copy map only knew the group (``output/E.coli_IVT``),
    # so the sample is recovered from the tool's own "test side", exactly as the
    # legacy pipeline named the WT output directories:
    #   EpiNano_DiffErr  ``-w``      -> E_ss_rd_RNA*_result
    #   xPore            ``wt``      -> E_ss_rd_RNA*_result
    #   yanocomp         ``-t``      -> E_ss_rd_RNA*_result_output.bed
    #   DRUMMER          ``-t``      -> kept on E_IVT_neg1 (script: -t neg1 -c neg2)
    "EpiNano_Error": {"E_IVT_neg_result": "E_IVT_neg2"},
    "xPore": {"E_IVT_neg_result": "E_IVT_neg2"},
    "yanocomp": {"E_IVT_neg_result": "E_IVT_neg2"},
    "DRUMMER": {"E_IVT_neg": "E_IVT_neg1"},
}

#: Explicit **comparison-pair** directories (``<A>_vs_<B>``) and the sample each
#: tool's output is filed under (user decision 2026-09-16: rename only, ownership
#: unchanged, so every evaluation number stays identical).
#:
#: ``registry.canonical_from_dir`` deliberately refuses to attribute a ``_vs_``
#: directory by name (the owning side differs per tool), so every pair directory
#: MUST be declared here -- otherwise the tool's callset for that sample silently
#: becomes "missing".  These are *declarations*, not guesses: ``sample_of_dir``
#: returns them with ``ambiguous=False`` (they must not appear in the
#: ``00_build_registry`` ambiguous-directory warning list).
#:
#: * ``E_IVT_neg1_vs_E_IVT_neg2`` -- the single E.coli IVT null-vs-null
#:   comparison (both halves of one library, SRR27228854).  Owners follow the
#:   detection scripts: DRUMMER ``-t neg1`` and ELIGOS2_diff (first side) ->
#:   ``E_IVT_neg1``; EpiNano_DiffErr ``-w``, xPore yml ``wt``, yanocomp ``-t``
#:   and Nanocompore ``label2`` -> ``E_IVT_neg2``.
#: * ``mESCs_Mettl3_KO_vs_mES_KO`` -- the cross-study mouse KO-vs-KO comparison.
#:   DRUMMER ``-t mESCs_Mettl3_KO``, ELIGOS2_diff (first side), xPore yml
#:   ``wt: mES_KO`` and yanocomp (legacy ``mES_KO`` output name) -> ``mES_KO``;
#:   EpiNano_DiffErr (``-k`` output name) and Nanocompore ``label1`` ->
#:   ``mESCs_Mettl3_KO``.  Keeping the legacy ownership is what keeps the
#:   published KO-side numbers reproducible.
PAIR_DIR_SAMPLE: dict[str, dict[str, str]] = {
    "ELIGOS2_diff": {"E_IVT_neg1_vs_E_IVT_neg2": "E_IVT_neg1",
                     "mESCs_Mettl3_KO_vs_mES_KO": "mES_KO"},
    "DRUMMER": {"E_IVT_neg1_vs_E_IVT_neg2": "E_IVT_neg1",
                "mESCs_Mettl3_KO_vs_mES_KO": "mES_KO"},
    "EpiNano_Error": {"E_IVT_neg1_vs_E_IVT_neg2": "E_IVT_neg2",
                      "mESCs_Mettl3_KO_vs_mES_KO": "mESCs_Mettl3_KO"},
    "xPore": {"E_IVT_neg1_vs_E_IVT_neg2": "E_IVT_neg2",
              "mESCs_Mettl3_KO_vs_mES_KO": "mES_KO"},
    "yanocomp": {"E_IVT_neg1_vs_E_IVT_neg2": "E_IVT_neg2",
                 "mESCs_Mettl3_KO_vs_mES_KO": "mES_KO"},
    "Nanocompore": {"E_IVT_neg1_vs_E_IVT_neg2": "E_IVT_neg2",
                    "mESCs_Mettl3_KO_vs_mES_KO": "mESCs_Mettl3_KO"},
}

#: raw results that the legacy ``output/`` assembly NEVER used.  Recorded so the
#: completeness audit can report them explicitly instead of leaving a "silent
#: gap": DRUMMER ran both Curlcake *IVT-vs-IVT* comparisons and ELIGOS2_diff
#: wrote two empty (1-byte) comparison files, but the legacy copy map only took
#: DRUMMER result3/4 (m6A) and ``output/Curlcake_IVT/{DRUMMER,ELIGOS2_diff}/``
#: stayed empty.  Deliberately not filled, matching the published assembly.
#:
#: Keys are canonical sample names (``config.SAMPLES_BY_NAME``) and the paths are
#: relative to :data:`config.PROJECT`, pointing at the **current** post-tidy
#: location under ``02_raw_results/result/``.  The ``(sample, tool)`` pairing
#: follows the legacy file-stem index: ``Curlcake_IVT_result1_*`` -> ``rep1``,
#: ``Curlcake_IVT_result2_*`` -> ``rep2_partial`` (the same library index the
#: registry uses for these tools), while the ELIGOS2_diff files name their first
#: side (the treatment library) in the file name itself.
LEGACY_UNUSED_RAW: dict[tuple[str, str], str] = {
    ("Curlcake_IVT_rep1", "DRUMMER"):
        "02_raw_results/result/DRUMMER/Curlcake_IVT_rep2_partial_vs_IVT_rep1/summary.txt",
    ("Curlcake_IVT_rep2_partial", "DRUMMER"):
        "02_raw_results/result/DRUMMER/Curlcake_IVT_rep3_vs_IVT_rep1/summary.txt",
    ("Curlcake_IVT_rep2_partial", "ELIGOS2_diff"):
        "02_raw_results/result/ELIGOS2_diff/Curlcake_IVT_rep2_partial_vs_IVT_rep1/"
        "RNA081120181part_vs_RNAAB089716_on_Curlcake_baseExt0.txt",
    ("Curlcake_IVT_rep3", "ELIGOS2_diff"):
        "02_raw_results/result/ELIGOS2_diff/Curlcake_IVT_rep3_vs_IVT_rep1/"
        "RNA081120181_vs_RNAAB089716_on_Curlcake_baseExt0.txt",
}


def needs_liftover(tool: str, species: str) -> bool:
    if tool in GENOMIC_TOOLS:
        return False
    return (tool, species) not in NO_LIFTOVER


def build(tool: str, sample_dir: Path, mod_type: str,
          species: str | None = None) -> pd.DataFrame | None:
    """Dispatch to the tool builder (species-aware when required)."""
    fn = BUILDERS.get(tool)
    if fn is None:
        return None
    if tool in SPECIES_AWARE:
        return fn(sample_dir, mod_type, species)
    return fn(sample_dir, mod_type)


def post_liftover_filter(tool: str, species: str, df: pd.DataFrame) -> pd.DataFrame:
    """Apply the legacy chromosome whitelist to a liftover output (or callset).

    Works on both the raw R2Dtool frame (column ``chromosome``) and the parsed
    canonical frame (column ``chrom``, ``chr`` prefix normalised away).
    """
    keep = POST_LIFTOVER_KEEP.get(tool, {}).get(species)
    if not keep or df.empty:
        return df
    col = next((c for c in ("chromosome", "chrom", "Chr") if c in df.columns), None)
    if col is None:
        return df
    # case/prefix-insensitive: normalise BOTH sides to lower-case with an
    # optional 'chr' prefix removed.  E.coli's contig is literally named
    # 'Chromosome' -- stripping '^chr' unconditionally mangles it to
    # 'omosome' and silently dropped every rebuilt E.coli row (and would
    # drop X-chromosome sites for human/mouse).
    def _canon(name: str) -> str:
        s = str(name).lower()
        return s[3:] if s.startswith("chr") else s

    keep_ci = {_canon(k) for k in keep}
    norm = df[col].astype(str).str.lower().map(_canon)
    return df[norm.isin(keep_ci)].reset_index(drop=True)

#: species -> GTF used by the legacy R2Dtool invocations
SPECIES_GTF = {
    "Arabidopsis": Path(str(_XB / "reference/arabidopsis/Arabidopsis_thaliana.TAIR10.61.gtf")),
    "Mouse": Path(str(_XB / "reference/GRCm39/ensembl/Mus_musculus.GRCm39.114.gtf")),
    "Human": Path(str(_XB / "reference/GRCh38p14/ensembl112/Homo_sapiens.GRCh38.112.chr.gtf")),
    "E.coli": Path(str(_XB / "reference/K_12/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.61.gtf")),
}

R2D_BIN = Path(str(_XB / "source_code/nanopore/R2Dtool/target/release/r2d"))
