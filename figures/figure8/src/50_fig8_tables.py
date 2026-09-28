#!/usr/bin/env python3
"""50 -- Figure 8 / Figure S8 revision tables (reviewer R3-8, editor E8).

Everything the revised Figure 8 (main) and Figure S8 (supplementary) plot is
recomputed here from the clean callset layer (``harmonisation/callsets``), the
frozen evaluation tables (``harmonisation/evaluation/tables``), the per-sample
candidate universes (``harmonisation/universe``) and the raw RNA004 ORCA output.

Site-level semantics follow the 06 confusion tables: one row per
``(chrom, start)`` position, keeping the row with the highest score (callsets
keeps overlapping-transcript duplicates of the same physical position).

FPR definitions (they are repeated verbatim in the figure legends)
------------------------------------------------------------------
* ``n_fp``                 -- calls on a genuinely unmodified control; every
                              call there is a false positive by construction;
* ``fp_per_10kb``           -- ``n_fp / mappable_bp x 1e4`` with the mappable
                              length anchored on ``controls_ivt_fpr.tsv``
                              (Curlcake 10,135 bp; HeLa 157,105,601 bp);
* ``fp_per_1e6_candidates`` -- ``n_fp / n_candidates x 1e6`` where the candidate
                              set is the sample's own universe (coverage >= 10,
                              reference base compatible with the modification
                              type) -- the same universe the paper's TP/FP
                              definitions use.

Outputs -> ``figures/figure8/tables``
--------------------------------------------------------
* fig8_counts.tsv                 detected sites per sample x tool (WT / IVT)
* fig8_fpr_curlcake_scan.tsv      Dorado threshold sweep, unmodified Curlcake
* fig8_fpr_hela_ivt.tsv          HeLa IVT control per tool (delivered cutoff)
* fig8_ppv_glori.tsv             PPV vs. GLORI (2 bp), RNA004 HeLa WT
* fig8_tradeoff.tsv              PPV vs. FP-per-10 kb operating points
* fig8_jaccard_wt.tsv / _ivt.tsv Dorado m6A concordance (moves to Fig. S8)
* figS9_chem_compare.tsv         Curlcake FPR, RNA002 vs RNA004, per unit
* figS9_orca_counts.tsv          ORCA + tool counts per modification type
* figS9_modratio_agreement.tsv   modification ratio vs GLORI association
* fig8_candidate_denominators.tsv candidate-site denominators actually used
* fig8_quoted_numbers.tsv  numbers quoted in the reply, with sources

Usage
-----
conda run -n benchmark-revision --no-capture-output \
    python $RNAMODBENCH_ROOT/figures/figure8/src/50_fig8_tables.py
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
import gzip
import math
import re
import sys
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(str(_RB))
SITES = (_RB / "data")
sys.path.insert(0, str((_RB / "src/harmonisation")))

from common.match import fix_chromosome  # noqa: E402  (project chromosome normaliser)
CLEAN = (_RB / "data/callsets")
EVAL = (_RB / "data/evaluation/tables")
UNIVERSE = (_XB / "harmonisation/universe")
GLORI_HELA = (_XB / "third_party/NGS/GLORI/Hela_GLORI.bed")
ORCA_DIR = (_XB / "raw/result_RNA004/HeLa/raw_calls/RNA004_result/ORCA_filtered")
ORCA_COUNTS = (_XB / "raw/result_RNA004/HeLa/figures/barplot/tool_counts_orca.csv")
OUT = (_RB / "figures/figure8")
TABLES = (_RB / "figures/figure8/tables")

DORADO_THRESHOLDS = (5, 10, 20, 50, 90)  # percent_modified
PLOT_THRESHOLDS = (5, 10, 20, 50)

#: RNA004 samples that feed the figure (other mod types included for S8A/S8C)
RNA004_SAMPLES = {
    "HeLa_RNA004_WT": ("RNA004", "Human", "RNA004_HeLa_WT"),
    "HeLa_RNA004_IVT": ("RNA004", "Human", "RNA004_HeLa_IVT"),
    "Curlcake_RNA004_IVT": ("RNA004", "Curlcake", "RNA004_Curlcake_IVT"),
}
#: RNA002 controls used by the chemistry comparison (S8B)
RNA002_SAMPLES = {
    "Curlcake_IVT_rep1": ("RNA002", "Curlcake", "Curlcake_IVT"),
    "Curlcake_IVT_rep2_partial": ("RNA002", "Curlcake", "Curlcake_IVT"),
    "Curlcake_IVT_rep3": ("RNA002", "Curlcake", "Curlcake_IVT"),
    "HeLa_IVT_rep1": ("RNA002", "Human", "HeLa_IVT"),
    "HeLa_IVT_rep2": ("RNA002", "Human", "HeLa_IVT"),
    "HeLa_IVT_rep3": ("RNA002", "Human", "HeLa_IVT"),
}

#: Dorado model-set families (parsed from the tool directory name)
FAMILY_LABEL = {
    "m6A_DRACH": "m6A DRACH",
    "m6A": "m6A general",
    "pseU_m6A": "pseU + m6A",
    "inosine_m6A": "inosine + m6A",
    "all": "all-modification",
    "m5C": "m5C",
    "pseU": "pseU",
    "inosine": "inosine",
}
#: base sets compatible with each modification type (mirrors common.evaluation)
MOD_BASES = {"m6A": ("A", "T"), "inosine": ("A", "T"), "Psi": ("A", "T"),
             "m1Psi": ("A", "T"), "m5C": ("C", "G")}


def log(msg: str) -> None:
    print(f"[50_fig8_tables] {msg}", flush=True)


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def read_tsv(path: Path, **kw) -> pd.DataFrame:
    if not Path(path).exists() or Path(path).stat().st_size == 0:
        return pd.DataFrame()
    kw.setdefault("sep", "\t")
    kw.setdefault("dtype", str)
    kw.setdefault("keep_default_na", False)
    return pd.read_csv(path, **kw)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, float_format="%.6g")
    log(f"wrote {path.relative_to(PROJECT)}  ({len(df)} rows)")


def parse_dorado(tool: str) -> dict:
    """``Dorado_hac@v5.0.0_m6A_DRACH@v1`` -> mode/version/family."""
    m = re.match(r"Dorado_(?P<mode>hac|sup)@(?P<ver>v[\d.]+)_(?P<rest>.+)", tool)
    if not m:
        return {"mode": "", "version": "", "family": ""}
    rest = re.sub(r"@v\d+$", "", m["rest"]).replace("_otherMod", "")
    return {"mode": m["mode"], "version": m["ver"], "family": rest}


def list_callsets(sample: str) -> list[tuple[str, str, str, Path]]:
    """[(mod_type, tool, path)] for one sample of the clean layer."""
    platform, species, group = (
        RNA004_SAMPLES.get(sample) or RNA002_SAMPLES[sample]
    )
    base = CLEAN / platform / species / group
    out = []
    if not base.exists():
        return out
    for mod_dir in sorted(p for p in base.iterdir() if p.is_dir()):
        for tool_dir in sorted(p for p in mod_dir.iterdir() if p.is_dir()):
            f = tool_dir / f"{sample}.tsv"
            if f.exists():
                out.append((mod_dir.name, tool_dir.name, f))
    return out


def load_positions(path: Path, *, keep_score: bool = False) -> pd.DataFrame:
    """Site-level callset: one row per (chrom, start), highest score kept."""
    df = read_tsv(path, usecols=lambda c: c in ("chrom", "start", "end", "strand",
                                                "score", "score_type", "coverage",
                                                "mod_ratio", "ref_base"))
    if df.empty:
        return df
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["start"])
    if "score" in df.columns:
        df["score_num"] = pd.to_numeric(df["score"], errors="coerce")
    else:
        df["score_num"] = np.nan
    if keep_score:
        df = df.sort_values("score_num", ascending=False, na_position="last")
    df = df.drop_duplicates(["chrom", "start"])
    return df.reset_index(drop=True)


# --------------------------------------------------------------------------- #
# denominators
# --------------------------------------------------------------------------- #
def curlcake_universe(uni_path: Path) -> dict:
    """Candidate-site denominators of one Curlcake sample universe."""
    df = read_tsv(uni_path,
                  usecols=lambda c: c in ("chrom", "pos", "base", "coverage"))
    df["coverage"] = pd.to_numeric(df["coverage"], errors="coerce")
    df = df[df["coverage"] >= 10]
    out = {"mappable_bp": None}
    for mod, bases in MOD_BASES.items():
        out[f"n_cand_{mod}"] = int(df["base"].isin(bases).sum())
    out["n_cand_Nm"] = int(len(df))
    return out


def hela_universe(sample: str, cache: Path) -> dict:
    """Candidate-site denominators of a HeLa sample universe (one awk pass)."""
    if cache.exists():
        df = read_tsv(cache)
        row = df[df["sample"] == sample]
        if not row.empty:
            return {k: float(row.iloc[0][k]) for k in
                    ("at_cov10", "cg_cov10", "n_all_cov10")}
    gz = (_XB / "harmonisation/universe/RNA004/Human") / f"{sample}__universe.tsv.gz"
    cmd = (f"zcat {gz} | awk -F'\\t' 'NR>1 && $4>=10 "
           "{b=$3; if(b==\"A\"||b==\"T\") at++; else if(b==\"C\"||b==\"G\") cg++; "
           "tot++} END{print at+0, cg+0, tot+0}'")
    at, cg, tot = (int(v) for v in subprocess.check_output(
        cmd, shell=True, text=True).split())
    log(f"universe {sample}: A/T cov>=10 = {at:,}, C/G = {cg:,}, all = {tot:,}")
    return {"at_cov10": at, "cg_cov10": cg, "n_all_cov10": tot}


def candidates_for(den: dict, mod_type: str) -> int:
    if mod_type in ("m6A", "inosine", "Psi", "m1Psi"):
        return int(den["at_cov10"])
    if mod_type == "m5C":
        return int(den["cg_cov10"])
    return int(den["n_all_cov10"])


# --------------------------------------------------------------------------- #
# 1. counts per sample x tool (Fig. 8A / S8A)
# --------------------------------------------------------------------------- #
def counts_table() -> pd.DataFrame:
    rows = []
    for sample in list(RNA004_SAMPLES) + list(RNA002_SAMPLES):
        for mod_type, tool, path in list_callsets(sample):
            df_raw = read_tsv(path, usecols=lambda c: c in ("chrom", "start", "score", "coverage"))
            n_rows = len(df_raw)
            df = load_positions(path, keep_score=True)
            fam = parse_dorado(tool) if tool.startswith("Dorado") else {}
            rows.append({
                "platform": (RNA004_SAMPLES.get(sample) or RNA002_SAMPLES[sample])[0],
                "sample": sample,
                "mod_type": mod_type,
                "tool": tool,
                "family": fam.get("family", ""),
                "mode": fam.get("mode", ""),
                "version": fam.get("version", ""),
                "n_rows": n_rows,
                "n_sites": len(df),
                "median_coverage": (float(pd.to_numeric(df["coverage"], errors="coerce")
                                          .median()) if "coverage" in df.columns and len(df) else np.nan),
            })
    out = pd.DataFrame(rows).sort_values(["platform", "sample", "mod_type", "tool"])
    write_tsv(out, (_RB / "figures/figure8/tables/fig8_counts.tsv"))
    return out


# --------------------------------------------------------------------------- #
# 2. Curlcake IVT FPR, Dorado threshold sweep (Fig. 8B / S8D)
# --------------------------------------------------------------------------- #
def fpr_curlcake_scan(region_bp: float) -> pd.DataFrame:
    sample = "Curlcake_RNA004_IVT"
    den = curlcake_universe((_XB / "harmonisation/universe/RNA004/Curlcake/Curlcake_RNA004_IVT__universe.tsv"))
    den["mappable_bp"] = region_bp
    rows = []
    for mod_type, tool, path in list_callsets(sample):
        df = load_positions(path, keep_score=True)
        if df.empty:
            continue
        n_cand = den.get(f"n_cand_{mod_type}", den["n_cand_Nm"])
        fam = parse_dorado(tool) if tool.startswith("Dorado") else {}
        score_pct = df["score_num"] * 100 if "score_num" in df.columns else pd.Series(dtype=float)
        thresholds = (list(DORADO_THRESHOLDS) if tool.startswith("Dorado") else [None])
        for thr in thresholds:
            n_fp = int((score_pct >= thr).sum()) if thr is not None else int(len(df))
            rows.append({
                "sample": sample, "mod_type": mod_type, "tool": tool,
                "family": fam.get("family", ""), "mode": fam.get("mode", ""),
                "version": fam.get("version", ""),
                "threshold_pct": thr,
                "n_fp": n_fp,
                "n_tested": int(len(df)),
                "fp_per_10kb": 1e4 * n_fp / den["mappable_bp"],
                "fp_per_1e6_candidates": 1e6 * n_fp / n_cand if n_cand else np.nan,
                "n_candidates": n_cand,
                "mappable_bp": den["mappable_bp"],
            })
    out = pd.DataFrame(rows)
    write_tsv(out, (_RB / "figures/figure8/tables/fig8_fpr_curlcake_scan.tsv"))
    return out


# --------------------------------------------------------------------------- #
# 3. HeLa IVT / WT FPR control (Fig. 8B right / 8D / S8C)
# --------------------------------------------------------------------------- #
def fpr_hela(region_bp: float, den_wt: dict, den_ivt: dict) -> pd.DataFrame:
    rows = []
    for sample, den in (("HeLa_RNA004_WT", den_wt), ("HeLa_RNA004_IVT", den_ivt)):
        for mod_type, tool, path in list_callsets(sample):
            df = load_positions(path, keep_score=True)
            if df.empty:
                continue
            fam = parse_dorado(tool) if tool.startswith("Dorado") else {}
            n_cand = candidates_for(den, mod_type)
            rows.append({
                "sample": sample, "mod_type": mod_type, "tool": tool,
                "family": fam.get("family", ""), "mode": fam.get("mode", ""),
                "version": fam.get("version", ""),
                "n_fp": int(len(df)),
                "min_score_pct": float(df["score_num"].min() * 100) if df["score_num"].notna().any() else np.nan,
                "fp_per_10kb": 1e4 * len(df) / region_bp,
                "fp_per_1e6_candidates": 1e6 * len(df) / n_cand if n_cand else np.nan,
                "n_candidates": n_cand,
                "mappable_bp": region_bp,
            })
    out = pd.DataFrame(rows)
    write_tsv(out, (_RB / "figures/figure8/tables/fig8_fpr_hela_ivt.tsv"))
    return out


# --------------------------------------------------------------------------- #
# 4. PPV vs. GLORI (2 bp), RNA004 HeLa WT
# --------------------------------------------------------------------------- #
def ppv_glori() -> pd.DataFrame:
    df = read_tsv((_RB / "data/evaluation/tables/m6a_glori_confusion.tsv"))
    df = df[(df["platform"] == "RNA004") & (df["sample"] == "HeLa_RNA004_WT")
            & (df["window"].astype(int) == 2)].copy()
    for c in ("tp", "fp", "fn", "tn", "universe", "n_calls_in_universe",
              "precision", "recall", "mcc", "precision_ci_lo", "precision_ci_hi"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.rename(columns={"precision": "ppv_glori_w2", "recall": "glori_coverage_w2",
                            "n_calls_in_universe": "n_calls_in_universe"})
    fam = [parse_dorado(t) if t.startswith("Dorado") else {} for t in df["tool"]]
    df["family"] = [f.get("family", "") for f in fam]
    df["mode"] = [f.get("mode", "") for f in fam]
    keep = ["platform", "species", "dataset_group", "sample", "tool", "family",
            "mode", "tp", "fp", "fn", "tn", "universe", "n_calls_in_universe",
            "ppv_glori_w2", "glori_coverage_w2", "mcc",
            "precision_ci_lo", "precision_ci_hi", "n_calls_total"]
    out = df[[c for c in keep if c in df.columns]].sort_values(
        "ppv_glori_w2", ascending=False)
    write_tsv(out, (_RB / "figures/figure8/tables/fig8_ppv_glori.tsv"))
    return out


# --------------------------------------------------------------------------- #
# 5. PPV vs. FPR trade-off (Fig. 8D)
# --------------------------------------------------------------------------- #
def tradeoff(ppv: pd.DataFrame, hela: pd.DataFrame) -> pd.DataFrame:
    ivt = hela[hela["sample"] == "HeLa_RNA004_IVT"].copy()
    merged = ppv.merge(ivt[["tool", "n_fp", "fp_per_10kb", "fp_per_1e6_candidates"]],
                       on="tool", how="inner")
    merged["operating_point"] = np.where(merged["tool"].str.startswith("Dorado"),
                                         "delivered calls, percent_modified >= 90",
                                         "tool default threshold")
    out = merged.sort_values("fp_per_10kb")
    write_tsv(out, (_RB / "figures/figure8/tables/fig8_tradeoff.tsv"))
    return out


# --------------------------------------------------------------------------- #
# 6. Jaccard concordance of the Dorado m6A models (moves to S8E)
# --------------------------------------------------------------------------- #
def jaccard_tables() -> None:
    for sample, stem in (("HeLa_RNA004_WT", "fig8_jaccard_wt.tsv"),
                         ("HeLa_RNA004_IVT", "fig8_jaccard_ivt.tsv")):
        sets: dict[str, set] = {}
        for mod_type, tool, path in list_callsets(sample):
            if mod_type != "m6A" or not tool.startswith("Dorado"):
                continue
            df = load_positions(path)
            sets[tool] = set(zip(df["chrom"], df["start"].astype("Int64")))
        names = sorted(sets)
        mat = pd.DataFrame(index=names, columns=names, dtype=float)
        for a in names:
            for b in names:
                inter = len(sets[a] & sets[b])
                union = len(sets[a] | sets[b])
                mat.loc[a, b] = inter / union if union else np.nan
        mat = mat.reset_index().rename(columns={"index": "tool"})
        write_tsv(mat, TABLES / stem)


# --------------------------------------------------------------------------- #
# 7. Curlcake FPR, RNA002 vs RNA004, per sequencing unit (S8B)
# --------------------------------------------------------------------------- #
def chem_compare(region: dict) -> pd.DataFrame:
    rows = []
    for sample in list(RNA002_SAMPLES) + ["Curlcake_RNA004_IVT"]:
        platform = "RNA004" if sample.startswith(("HeLa_RNA004", "Curlcake_RNA004")) else "RNA002"
        species = "Curlcake" if sample.startswith("Curlcake") else "Human"
        region_bp = float(region.get(species, np.nan))
        for mod_type, tool, path in list_callsets(sample):
            if mod_type not in ("m6A", "Psi"):
                continue
            df = load_positions(path)
            n_fp = int(len(df))
            rows.append({
                "platform": platform, "sample": sample, "species": species,
                "mod_type": mod_type, "tool": tool, "n_fp": n_fp,
                "fp_per_10kb": 1e4 * n_fp / region_bp,
            })
    out = pd.DataFrame(rows)
    write_tsv(out, (_RB / "figures/figure8/tables/figS9_chem_compare.tsv"))
    return out


# --------------------------------------------------------------------------- #
# 8. ORCA counts per modification type (S8A)
# --------------------------------------------------------------------------- #
def orca_counts() -> pd.DataFrame:
    rows = []
    if ORCA_DIR.exists():
        for mod_dir in sorted(p for p in ORCA_DIR.iterdir() if p.is_dir()):
            for f in sorted(mod_dir.glob("*.txt")):
                df = read_tsv(f, usecols=lambda c: c in ("Chr", "Start", "End", "Status"))
                n_rows = len(df)
                if df.empty:
                    n_sites = 0
                else:
                    df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
                    n_sites = int(df.dropna(subset=["Start"])
                                  .drop_duplicates(["Chr", "Start"]).shape[0])
                sample = ("WT" if "WT" in f.name else "IVT")
                rows.append({"source": "ORCA_filtered", "sample": sample,
                             "mod_type": mod_dir.name, "tool": f"ORCA_{mod_dir.name}",
                             "n_rows": n_rows, "n_sites": n_sites})
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.pivot_table(index=["source", "mod_type", "tool"], columns="sample",
                              values="n_sites", aggfunc="sum").reset_index()
        for c in ("WT", "IVT"):
            if c not in out.columns:
                out[c] = 0
        out = out[["source", "mod_type", "tool", "WT", "IVT"]]
    write_tsv(out, (_RB / "figures/figure8/tables/figS9_orca_counts.tsv"))
    return out


# --------------------------------------------------------------------------- #
# 9. modification ratio vs. GLORI (S8F)
# --------------------------------------------------------------------------- #
def modratio_agreement() -> pd.DataFrame:
    glori = read_tsv(GLORI_HELA, header=None,
                     usecols=[0, 1, 3], names=["chrom", "start", "glori_ratio"])
    glori["start"] = pd.to_numeric(glori["start"], errors="coerce")
    glori["glori_ratio"] = pd.to_numeric(glori["glori_ratio"], errors="coerce")
    # the GLORI BED carries bare contig names ("1", "x"); normalise both sides
    glori["chrom"] = [fix_chromosome(c) for c in glori["chrom"]]
    glori = glori.dropna(subset=["start"]).drop_duplicates(["chrom", "start"])
    rows = []
    for mod_type, tool, path in list_callsets("HeLa_RNA004_WT"):
        if mod_type != "m6A":
            continue
        df = read_tsv(path, usecols=lambda c: c in ("chrom", "start", "score", "mod_ratio"))
        if df.empty:
            continue
        df["start"] = pd.to_numeric(df["start"], errors="coerce")
        ratio = pd.to_numeric(df["mod_ratio"] if "mod_ratio" in df.columns else
                              df["score"], errors="coerce")
        df = df.assign(ratio=ratio).dropna(subset=["start", "ratio"])
        df["chrom"] = [fix_chromosome(c) for c in df["chrom"]]
        j = df.merge(glori, on=["chrom", "start"], how="inner").drop_duplicates(["chrom", "start"])
        if len(j) < 10:
            continue
        r = float(np.corrcoef(j["ratio"], j["glori_ratio"])[0, 1])
        rho = float(pd.Series(j["ratio"]).corr(pd.Series(j["glori_ratio"]), method="spearman"))
        ccc = 2 * np.cov(j["ratio"], j["glori_ratio"])[0, 1] / (
            j["ratio"].var() + j["glori_ratio"].var()
            + (j["ratio"].mean() - j["glori_ratio"].mean()) ** 2)
        rows.append({"tool": tool, "n_matched": int(len(j)),
                     "pearson_r": r, "spearman_rho": rho, "ccc": float(ccc),
                     "ratio_column": "mod_ratio" if "mod_ratio" in df.columns else "score",
                     "glori_sites_in_bed": int(len(glori))})
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values("pearson_r", ascending=False)
    write_tsv(out, (_RB / "figures/figure8/tables/figS9_modratio_agreement.tsv"))
    return out


# --------------------------------------------------------------------------- #
# 10. response-letter anchors
# --------------------------------------------------------------------------- #
def anchors(ppv: pd.DataFrame, cc: pd.DataFrame, hela: pd.DataFrame) -> pd.DataFrame:
    rows = []

    def add(name, value, unit, source, detail):
        rows.append({"anchor": name, "value": value, "unit": unit,
                     "source_table": source, "detail": detail})

    for _, r in ppv.iterrows():
        add(f"PPV vs GLORI (2 bp) | HeLa WT | {r['tool']}", round(r["ppv_glori_w2"], 6),
            "fraction", "fig8_ppv_glori.tsv",
            f"tp={int(r['tp'])} fp={int(r['fp'])} n_calls_in_universe={int(r['n_calls_in_universe'])}")
    dr = cc[(cc["family"] == "m6A_DRACH") & (cc["threshold_pct"].isin([10, 50]))]
    for _, r in dr.iterrows():
        add(f"Curlcake IVT FP/10kb | {r['tool']} @ {int(r['threshold_pct'])}%",
            round(r["fp_per_10kb"], 4), "per 10 kb", "fig8_fpr_curlcake_scan.tsv",
            f"n_fp={int(r['n_fp'])}, per 1e6 candidates={r['fp_per_1e6_candidates']:.0f}")
    gen = cc[(cc["family"] == "m6A") & (cc["threshold_pct"] == 50)]
    for _, r in gen.iterrows():
        add(f"Curlcake IVT FP/10kb | {r['tool']} @ 50%", round(r["fp_per_10kb"], 4),
            "per 10 kb", "fig8_fpr_curlcake_scan.tsv",
            f"n_fp={int(r['n_fp'])}, per 1e6 candidates={r['fp_per_1e6_candidates']:.0f}")
    ivt = hela[hela["sample"] == "HeLa_RNA004_IVT"]
    for _, r in ivt[ivt["mod_type"] == "m6A"].iterrows():
        add(f"HeLa IVT FP/10kb | {r['tool']} (delivered cutoff)",
            round(r["fp_per_10kb"], 5), "per 10 kb", "fig8_fpr_hela_ivt.tsv",
            f"n_fp={int(r['n_fp'])}, min_score_pct={r['min_score_pct']:.0f}")
    out = pd.DataFrame(rows)
    write_tsv(out, (_RB / "figures/figure8/tables/fig8_quoted_numbers.tsv"))
    return out


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-heavy", action="store_true",
                    help="skip the HeLa universe awks (use cache / frozen fallback)")
    args = ap.parse_args()
    TABLES.mkdir(parents=True, exist_ok=True)

    # ---- region lengths anchored on the frozen control table ---------------
    ctrl = read_tsv((_RB / "data/evaluation/tables/controls_ivt_fpr.tsv"))
    ctrl["region_bp"] = pd.to_numeric(ctrl["region_bp"], errors="coerce")
    region = (ctrl.dropna(subset=["region_bp"])
              .groupby("species")["region_bp"].max().to_dict())
    cc_bp = float(region.get("Curlcake", 10135))
    hela_bp = float(region.get("Human", 157105601))
    log(f"mappable bp: Curlcake={cc_bp:.0f}  HeLa={hela_bp:.0f}")

    cache = (_RB / "figures/figure8/tables/fig8_candidate_denominators.tsv")
    if args.skip_heavy and cache.exists():
        log("universe denominators: cache only")
        den_ivt = {"at_cov10": 15817652, "cg_cov10": 14155581, "n_all_cov10": np.nan}
        den_wt = {"at_cov10": 13563513, "cg_cov10": 12217411, "n_all_cov10": np.nan}
    else:
        den_ivt = hela_universe("HeLa_RNA004_IVT", cache)
        den_wt = hela_universe("HeLa_RNA004_WT", cache)
        write_tsv(pd.DataFrame([
            {"sample": "HeLa_RNA004_IVT", **den_ivt},
            {"sample": "HeLa_RNA004_WT", **den_wt},
        ]), cache)

    counts_table()
    cc = fpr_curlcake_scan(cc_bp)
    hela = fpr_hela(hela_bp, den_wt, den_ivt)
    ppv = ppv_glori()
    tradeoff(ppv, hela)
    jaccard_tables()
    chem_compare(region)
    orca_counts()
    modratio_agreement()
    anchors(ppv, cc, hela)
    log("done")


if __name__ == "__main__":
    main()
