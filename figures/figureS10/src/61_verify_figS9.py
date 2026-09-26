#!/usr/bin/env python3
"""61 -- acceptance checks for the revised Supplementary Figure S9.

Checks
------
1. the page is exactly the page of the replaced ``sup9.pdf`` (1152 x 864 pt)
   and the 300 dpi PNG sibling exists;
2. every embedded font is Arial (no DejaVu fallback, nothing unembedded);
3. house rules in the plotting script: no grid element, no font below 7 pt, the
   six-model white list, the inosine models never drawn, THREE bold block
   letters A/B/C and no figure number (A = the six Guitar density panels as one
   block), and the bottom-left half on the purple/green version pair (no grey);
4. the figure inputs (``figS9_panel_inputs.tsv``, written by the R script: one
   row per drawn curve, block A + grid slot) match
   the source layer: BED line counts equal the de-duplicated ``callsets``
   call counts, and those equal the frozen 2026-09-20 anchors;
5. the numbers in ``figS9_region_shares.tsv`` are internally consistent
   (shares sum to 1, body counts add up) and the assigned-site counts equal the
   anchors that were frozen before the figure was drawn;
6. ``figS9_wt_ivt_contrast.tsv`` follows the stated +/-5 pp verdict rule and its
   numbers are reproduced verbatim in ``figS9_key_numbers.md``;
7. the drawn text actually present in the PDF: six model titles, the key labels,
   the region labels -- and no excluded model.

Usage
-----
conda run -n benchmark-revision --no-capture-output \
    python $RNAMODBENCH_ROOT/figures/figureS10/src/61_verify_figS9.py
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
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(str(_RB))
OUT = (_RB / "figures/figureS10")
TABLES, FIGS, SCRIPTS, LOGS = (_RB / "figures/figureS10/tables"), (_RB / "figures/figureS10/figures"), (_RB / "figures/figureS10/src"), (_RB / "figures/figureS10/logs")
CLEAN = (_RB / "data/callsets")
R_SCRIPT = (_RB / "src/harmonisation/scripts/23e_figS9_guitar.R")
PDF = FIGS / "FigureS9_rev.pdf"
PAGE = (1152.0, 864.0)                      # = the replaced sup9.pdf
EXPECTED_ASSIGNED = {                       # frozen 2026-09-20 (pre-figure)
    ("Dorado_hac@v5.0.0_pseU@v1_otherMod", "WT"): 94,
    ("Dorado_hac@v5.0.0_pseU@v1_otherMod", "IVT"): 56,
    ("Dorado_hac@v5.1.0_pseU@v1_otherMod", "WT"): 62,
    ("Dorado_hac@v5.1.0_pseU@v1_otherMod", "IVT"): 27,
    ("Dorado_sup@v5.0.0_pseU@v1_otherMod", "WT"): 221,
    ("Dorado_sup@v5.0.0_pseU@v1_otherMod", "IVT"): 175,
    ("Dorado_sup@v5.1.0_pseU@v1_otherMod", "WT"): 201,
    ("Dorado_sup@v5.1.0_pseU@v1_otherMod", "IVT"): 96,
    ("Dorado_hac@v5.1.0_m5C@v1_otherMod", "WT"): 396,
    ("Dorado_hac@v5.1.0_m5C@v1_otherMod", "IVT"): 136,
    ("Dorado_sup@v5.1.0_m5C@v1_otherMod", "WT"): 118,
    ("Dorado_sup@v5.1.0_m5C@v1_otherMod", "IVT"): 52,
}
MODELS = sorted({k[0] for k in EXPECTED_ASSIGNED})
SAMPLES = {"WT": ("RNA004_HeLa_WT", "HeLa_RNA004_WT"),
           "IVT": ("RNA004_HeLa_IVT", "HeLa_RNA004_IVT")}

FAILS: list[str] = []
CHECKS = [0]


def check(name: str, ok: bool, detail: str = "") -> None:
    CHECKS[0] += 1
    print(f"[{'OK  ' if ok else 'FAIL'}] {name}{'  -- ' + detail if detail else ''}")
    if not ok:
        FAILS.append(name)


def close(a: float, b: float, tol: float = 1e-6) -> bool:
    return abs(a - b) <= tol * max(abs(b), 1e-9)


def pdf_page_size() -> tuple[float, float]:
    txt = subprocess.check_output(["pdfinfo", str(PDF)], text=True)
    m = re.search(r"Page size:\s+([\d.]+) x ([\d.]+)", txt)
    return float(m.group(1)), float(m.group(2))


def pdf_fonts() -> list[tuple[str, str]]:
    txt = subprocess.check_output(["pdffonts", str(PDF)], text=True)
    out = []
    for line in txt.splitlines()[2:]:
        parts = line.split()
        if len(parts) < 6:
            continue
        emb = next((p for p in reversed(parts) if p in ("yes", "no")), "?")
        out.append((parts[0], emb))
    return out


def pdf_text() -> str:
    return subprocess.check_output(["pdftotext", "-layout", str(PDF), "-"], text=True)


def word_boxes(pdf: Path) -> list[tuple[float, float, float, float, str]]:
    """Word bounding boxes (xMin, yMin, xMax, yMax, text) from ``pdftotext -bbox``."""
    out = subprocess.run(["pdftotext", "-bbox", str(pdf), "-"],
                         capture_output=True, text=True, check=True).stdout
    return [(float(a), float(b), float(c), float(d), w) for a, b, c, d, w in
            re.findall(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" '
                       r'yMax="([\d.]+)">([^<]*)</word>', out)]


def clean_rows(tool: str, condition: str) -> int:
    group, sample = SAMPLES[condition]
    mod = "Psi" if "pseU" in tool else "m5C"
    path = (_RB / "data/callsets/RNA004/Human") / group / mod / tool / f"{sample}.tsv"
    d = pd.read_csv(path, sep="\t", usecols=["chrom", "start"], low_memory=False)
    return int(d.drop_duplicates(["chrom", "start"]).shape[0])


def main() -> None:
    LOGS.mkdir(parents=True, exist_ok=True)

    # 1 + 2 --------------------------------------------------------------- #
    size = pdf_page_size()
    check(f"page size {PAGE[0]:g} x {PAGE[1]:g} pt (replaced sup9.pdf)",
          close(size[0], PAGE[0], 1e-4) and close(size[1], PAGE[1], 1e-4),
          f"got {size}")
    check("PNG sibling exists", PDF.with_suffix(".png").exists())
    fonts = pdf_fonts()
    bad = [f for f in fonts if f[1] != "yes" or "Arial" not in f[0]]
    check("all fonts embedded Arial", bool(fonts) and not bad, f"offending: {bad}")

    # 3 ------------------------------------------------------------------- #
    src = R_SCRIPT.read_text()
    check("R script: panel grid blanked, no grid element",
          "panel.grid = element_blank()" in src
          and "panel.grid.major = element_line" not in src
          and "panel.grid.minor = element_line" not in src)
    fs_block = re.findall(r"FS <- list\((.*?)\)\n", src, re.S)[0]
    fs = [float(x) for x in re.findall(r"=\s*([\d.]+)", fs_block)]
    check("R script: smallest explicit font >= 7 pt", min(fs) >= 7.0, f"FS={fs}")
    check("R script: three bold block letters A/B/C, no figure number",
          'plot_annotation(tag_levels = "A")' in src
          and 'BLOCK_LETTERS <- c("A", "B", "C")' in src
          and 'face = "bold"' in src and "FIGURE_NUM" not in src
          and "MODELS$letter" not in src and "m$letter" not in src)
    check("R script: bottom-left half (block B) on the purple/green version pair",
          'B_V500 <- "#6A3D9A"' in src and 'B_V510 <- "#1B9E77"' in src
          and "grey20" not in src and "grey55" not in src)
    check("R script: white list = the six drawn models, no inosine panel",
          all(m in src for m in MODELS) and "inosine" in src
          and "EXCLUDED <- c(" in src)
    check("R script: callsets source, Ensembl annotation, >=90 % point stated",
          "callsets" in src and "Ensembl" in src and ">= 90 %" in src
          and "gencode" not in src.lower())

    # 4 ------------------------------------------------------------------- #
    plan = pd.read_csv(TABLES / "figS9_panel_inputs.tsv", sep="\t")
    check("panel inputs: 6 models x 2 conditions = 12 curves",
          len(plan) == 12 and plan["model"].nunique() == 6
          and set(plan["condition"]) == {"WT", "IVT"}, f"{len(plan)} rows")
    slots = [f"R{r}C{c}" for r in (1, 2) for c in (1, 2, 3)]
    check("panel inputs: every density panel sits in block A at its own slot",
          set(plan["block"]) == {"A"} and set(plan["slot"]) == set(slots)
          and plan.groupby("label")["slot"].nunique().eq(1).all()
          and plan.groupby("slot")["label"].nunique().eq(1).all(),
          f"slots {sorted(set(plan['slot']))}")
    ok_bed = True
    for _, r in plan.iterrows():
        n = len(Path(r["path"]).read_text().splitlines())
        ok_bed &= n == int(r["n_sites"])
    check("panel inputs: BED line counts equal the recorded n_sites", ok_bed)
    ok_src = all(int(r["n_sites"]) == clean_rows(r["model"], r["condition"])
                 for _, r in plan.iterrows())
    check("panel inputs: BED counts equal de-duplicated callsets counts", ok_src)

    # 5 ------------------------------------------------------------------- #
    shares = pd.read_csv(TABLES / "figS9_region_shares.tsv", sep="\t")
    segs = ["five_prime_UTR", "CDS", "three_prime_UTR"]
    share_sum = shares[[f"share_{s}" for s in segs]].sum(axis=1)
    check("region shares: 5'UTR + CDS + 3'UTR = 1 for every row",
          bool(((share_sum - 1).abs() < 1e-9).all()), f"max dev {abs(share_sum-1).max():.2e}")
    body = shares[[f"n_{s}" for s in segs]].sum(axis=1)
    check("region shares: segment counts add up to n_body",
          bool((body == shares["n_body"]).all()))
    ok_anchor = all(int(r.n_assigned) == EXPECTED_ASSIGNED[(r.model, r.condition)]
                    for _, r in shares.iterrows())
    check("region shares: assigned-site counts equal the frozen anchors",
          ok_anchor, f"{len(shares)} rows")
    check("region shares: delivered threshold >= 0.90 on every input",
          bool((shares["score_min"] >= 0.90).all()),
          f"min {shares['score_min'].min():.2f}")
    check("region shares: one sequencing unit per condition, stated",
          bool((shares["n_units"] == 1).all()))

    # 5b -- site-level bootstrap intervals (R1-8 / R1-9 / E7) -------------- #
    ci_ok = True
    for s in segs:
        lo = shares[f"share_{s}_ci_lo"].to_numpy(float)
        hi = shares[f"share_{s}_ci_hi"].to_numpy(float)
        point = shares[f"share_{s}"].to_numpy(float)
        # all-zero segments (e.g. a 5'UTR share of exactly 0) legitimately give a
        # degenerate interval, hence >= rather than > on the width
        ci_ok &= bool(np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))
                      and np.all(lo <= point + 1e-12)
                      and np.all(point <= hi + 1e-12) and np.all(hi >= lo))
    check("region shares: every segment share carries a 95 % CI bracketing it",
          ci_ok, f"{len(shares)} rows x {len(segs)} segments")
    t3 = shares["share_three_prime_UTR"].to_numpy(float)
    w3 = (shares["share_three_prime_UTR_ci_hi"]
          - shares["share_three_prime_UTR_ci_lo"]).to_numpy(float)
    check("region shares: every 3'UTR interval is non-degenerate", bool(np.all(w3 > 0)),
          f"narrowest {100*w3.min():.1f} pp")

    # independent recomputation of the very first drawn curve, same seed: the
    # interval must reproduce bit-for-bit (deterministic bootstrap contract)
    spec = importlib.util.spec_from_file_location(
        "figs9_tables", SCRIPTS / "60_figS9_tables.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mp = [p for p in sorted(((_XB / "reference/regionmodels")).glob("Human.*.mrna.regionmodel.pkl"))
          if "gencode" not in p.name.lower()][0]
    region = mod.RegionIndex.load(mp)
    tool0, cond0 = mod.MODELS[0][1], "WT"
    sites0 = mod.load_sites(mod.clean_path(mod.MODELS[0][0], tool0, cond0))
    kind0 = mod.assign_kinds(region, sites0)
    inside0 = np.flatnonzero(kind0 >= 0)
    body0 = kind0[inside0][np.isin(kind0[inside0], mod.BODY_KINDS)]
    ci0 = mod.bootstrap_ci({mod.KINDS[k]: (body0 == k) for k in mod.BODY_KINDS},
                           np.random.default_rng(mod.RNG_SEED))
    row0 = shares[(shares.model == tool0) & (shares.condition == cond0)].iloc[0]
    rec = max(abs(ci0[s][0] - row0[f"share_{s}_ci_lo"]) +
              abs(ci0[s][1] - row0[f"share_{s}_ci_hi"]) for s in segs)
    check("bootstrap interval reproduces bit-for-bit with the frozen seed",
          rec <= 1e-12, f"{tool0}/{cond0} max |diff| = {rec:.2e}")

    # 6 ------------------------------------------------------------------- #
    con = pd.read_csv(TABLES / "figS9_wt_ivt_contrast.tsv", sep="\t")
    rule_ok = all(
        (abs(r.three_prime_UTR_delta_pp) <= 5.0) ==
        (r.verdict == "reproduced on unmodified IVT")
        for _, r in con.iterrows())
    check("contrast: verdict follows the stated +/-5 pp rule", rule_ok)
    delta_ok = all(
        close(r.three_prime_UTR_delta_pp,
              100 * (r.share_three_prime_UTR_WT - r.share_three_prime_UTR_IVT), 1e-6)
        for _, r in con.iterrows())
    check("contrast: deltas equal the share differences", delta_ok)
    kn = (TABLES / "figS9_key_numbers.md").read_text()
    quoted = all(f"{100*r.share_three_prime_UTR_WT:.1f}%" in kn
                 and f"{100*r.share_three_prime_UTR_IVT:.1f}%" in kn
                 and f"{r.three_prime_UTR_delta_pp:+.1f}" in kn
                 for _, r in con.iterrows())
    check("key numbers: every WT/IVT 3'UTR share and delta is quoted", quoted)
    check("key numbers: states the single-unit structure (R3-2/E6)",
          "n = 1" in kn and "no biological replication" in kn)

    # 6b -- the difference intervals and the intervals quoted in the legend -- #
    dcols = ("delta_pp_ci_lo", "delta_pp_ci_hi", "delta_pp_ci_width")
    marg = {(r.model, r.condition): r for _, r in shares.iterrows()}
    ci_ok = all(c in con.columns for c in dcols)
    cover_zero = 0
    for _, r in con.iterrows():
        w, v = marg[(r.model, "WT")], marg[(r.model, "IVT")]
        ci_ok &= bool(r.delta_pp_ci_lo - 1e-9 <= r.three_prime_UTR_delta_pp
                      <= r.delta_pp_ci_hi + 1e-9)
        ci_ok &= bool(close(r.delta_pp_ci_width,
                            r.delta_pp_ci_hi - r.delta_pp_ci_lo, 1e-9))
        # a difference of two independent bootstraps can never be wider than the
        # extremes of the two marginal intervals
        ci_ok &= bool(r.delta_pp_ci_lo >= 100 * (w.share_three_prime_UTR_ci_lo
                                                 - v.share_three_prime_UTR_ci_hi) - 1e-9
                      and r.delta_pp_ci_hi <= 100 * (w.share_three_prime_UTR_ci_hi
                                                     - v.share_three_prime_UTR_ci_lo) + 1e-9)
        cover_zero += int(r.delta_pp_ci_lo <= 0 <= r.delta_pp_ci_hi)
    check("contrast: difference CIs bracket the point estimate and the marginals",
          ci_ok, f"{len(con)} models")
    check("contrast: the difference interval covers zero in all six models",
          cover_zero == 6, f"{cover_zero}/6")

    legend = (FIGS / "FigS9_legends.md").read_text()
    quoted_ci = all(f"[{100*r.share_three_prime_UTR_ci_lo:.1f}, "
                    f"{100*r.share_three_prime_UTR_ci_hi:.1f}]" in legend
                    for _, r in shares.iterrows())
    quoted_ci &= all(f"[{r.delta_pp_ci_lo:+.1f}, {r.delta_pp_ci_hi:+.1f}]" in legend
                     for _, r in con.iterrows())
    check("legend: every quoted CI resolves to the frozen tables", quoted_ci,
          "12 marginal + 6 difference intervals")

    # 6c -- reported values table: strict FPR at the delivered HeLa point -- #
    fpr = pd.read_csv(TABLES / "figS9_ivt_fpr.tsv", sep="\t")
    src = pd.read_csv((_RB / "data/evaluation/tables/controls_ivt_fpr.tsv"), sep="\t")
    src = src[src["sample"] == "HeLa_RNA004_IVT"]
    rel_ok, denom_ok = True, True
    for _, r in fpr.iterrows():
        # NB: bracket access -- ``r.mod`` is Series.mod, not the column
        s = src[(src["tool"] == r["model"]) & (src["mod_type"] == r["mod"])]
        rel_ok &= bool(len(s) == 1)
        if len(s) == 1:
            rel_ok &= bool(close(float(r["fp_per_10kb"]),
                                 float(s["fp_per_10kb"].iloc[0]), 1e-5)
                           and close(float(r["fp_per_1e6_candidates"]),
                                     float(s["fp_per_1e6_candidates"].iloc[0]), 1e-5))
        denom_ok &= bool(int(r["n_universe"])
                         == (15817652 if r["mod"] == "Psi" else 14155581)
                         and int(r["region_bp"]) == 157105601)
    check("reported values: six models, one row each, in the frozen drawing order",
          len(fpr) == 6 and set(fpr["block"]) == {"A"}
          and list(fpr["slot"]) == [f"R{r}C{c}" for r in (1, 2) for c in (1, 2, 3)],
          f"{len(fpr)} rows")
    check("reported values: every FPR equals the frozen evaluation table", rel_ok)
    check("reported values: denominators frozen candidate sets / mappable length",
          denom_ok)
    check("reported values: IVT call counts equal de-duplicated callsets counts",
          all(int(r["n_calls"]) == clean_rows(r["model"], "IVT")
              for _, r in fpr.iterrows()))

    # 6d -- block B (bottom left): Curlcake threshold scan (R3-8 / E8) ------- #
    scan = pd.read_csv(TABLES / "figS9_curlcake_scan.tsv", sep="\t")
    frac = len(scan) / scan["label"].nunique()
    check("block B: six models x ten thresholds",
          scan["label"].nunique() == 6 and abs(frac - 10) < 1e-9, f"{len(scan)} rows")
    mono = all(scan[scan["label"] == lab].sort_values("threshold_pct")["n_calls"]
               .diff().dropna().le(0).all() for lab in scan["label"].unique())
    check("block B: the threshold profile is non-increasing for every model", mono)
    def_ok = bool(np.allclose(scan["fp_per_10kb"], 1e4 * scan["n_calls"] / 10135)
                  and np.allclose(scan["fp_per_1e6_candidates"],
                                  1e6 * scan["n_calls"] / scan["n_candidates"])
                  and set(scan["n_candidates"].unique()) <= {4928, 3963})
    check("block B: densities follow the 10,135-bp / candidate definitions", def_ok)
    anchors_ok = True
    for lab, (raw, n5, n50, n90) in {
            "hac@v5.0.0_pseU": (2422, 89, 2, 0), "hac@v5.1.0_pseU": (2400, 107, 5, 0),
            "sup@v5.0.0_pseU": (2406, 70, 3, 0), "sup@v5.1.0_pseU": (2368, 111, 7, 0),
            "hac@v5.1.0_m5C": (2488, 246, 4, 0), "sup@v5.1.0_m5C": (2435, 106, 2, 0),
    }.items():
        sub = scan[scan["label"] == lab]
        c = dict(zip(sub["threshold_pct"].astype(int), sub["n_calls"].astype(int)))
        anchors_ok &= (c[5], c[50], c[90]) == (n5, n50, n90)
    check("block B: the frozen 5/50/90 % call anchors hold", anchors_ok)
    # independent recomputation straight from the call-set layer
    spec60 = importlib.util.spec_from_file_location(
        "figs9_tables_b", SCRIPTS / "60_figS9_tables.py")
    mod60 = importlib.util.module_from_spec(spec60)
    spec60.loader.exec_module(mod60)
    probe_model, probe_mod = "Dorado_sup@v5.1.0_all_m5C", "m5C"
    probe_path = (PROJECT / "data/callsets/RNA004/"
                  f"Curlcake/RNA004_Curlcake_IVT/{probe_mod}/{probe_model}/"
                  "Curlcake_RNA004_IVT.tsv")
    s = pd.read_csv(probe_path, sep="\t", usecols=["score"]).iloc[:, 0].to_numpy(float) * 100
    sub = scan[scan["model"] == probe_model].sort_values("threshold_pct")
    rec_ok = all(int(r["n_calls"]) == int((s >= r["threshold_pct"]).sum())
                 for _, r in sub.iterrows())
    check("block B: recomputed call counts match the table", rec_ok,
          f"{probe_model}")

    # 6e -- block C (bottom right): score validity in HeLa (R3-9) ----------- #
    validity = pd.read_csv(TABLES / "figS9_score_validity.tsv", sep="\t")
    vsum = pd.read_csv(TABLES / "figS9_score_validity_summary.tsv", sep="\t")
    check("block C: long table covers both conditions of all six models",
          validity["label"].nunique() == 6 and set(validity["condition"]) == {"WT", "IVT"}
          and len(validity) > 3000, f"{len(validity)} calls")
    n_ok = all(len(validity[(validity["label"] == r["label"])
                            & (validity["condition"] == "WT")]) == int(r["n_WT"])
               and len(validity[(validity["label"] == r["label"])
                                & (validity["condition"] == "IVT")]) == int(r["n_IVT"])
               for _, r in vsum.iterrows())
    check("block C: per-condition call counts agree with the summary", n_ok)
    base = shares.set_index(["model", "condition"])["n_calls_site"]
    call_ok = all(int(r["n_WT"]) == int(base.loc[(r["model"], "WT")])
                  and int(r["n_IVT"]) == int(base.loc[(r["model"], "IVT")])
                  for _, r in vsum.iterrows())
    check("block C: call counts equal the de-duplicated region-share counts",
          call_ok)
    exp_auc = {"hac@v5.0.0_pseU": 0.4767, "hac@v5.1.0_pseU": 0.5799,
               "sup@v5.0.0_pseU": 0.2394, "sup@v5.1.0_pseU": 0.6074,
               "hac@v5.1.0_m5C": 0.5486, "sup@v5.1.0_m5C": 0.5036}
    auc_ok = all(abs(float(r["auc_WT_vs_IVT"]) - exp_auc[r["label"]]) < 5e-4
                 for _, r in vsum.iterrows())
    check("block C: AUC values equal the frozen anchors", auc_ok)
    rule_ok = all((r["verdict"].startswith("inverted") if r["auc_WT_vs_IVT"] < 0.4
                   else r["verdict"] == "no discrimination" if r["auc_WT_vs_IVT"] <= 0.6
                   else r["verdict"] == "wild type scores higher")
                  for _, r in vsum.iterrows())
    check("block C: the AUC verdict rule is applied consistently", rule_ok)
    d_wt = validity[(validity["label"] == "hac@v5.0.0_pseU")
                    & (validity["condition"] == "WT")]["score"].to_numpy()
    d_ivt = validity[(validity["label"] == "hac@v5.0.0_pseU")
                     & (validity["condition"] == "IVT")]["score"].to_numpy()
    from scipy.stats import mannwhitneyu
    u = mannwhitneyu(d_wt, d_ivt, alternative="two-sided")
    auc_re = u.statistic / (len(d_wt) * len(d_ivt))
    rec = vsum[vsum["label"] == "hac@v5.0.0_pseU"]["auc_WT_vs_IVT"].iloc[0]
    check("block C: an independent AUC recomputation reproduces the table",
          abs(auc_re - float(rec)) < 5e-5, f"recomputed {auc_re:.4f} vs table {rec:.4f}")

    # 6f -- drawing script and the text actually on the page ----------------- #
    src_r = R_SCRIPT.read_text()
    check("blocks B/C: the bottom row reads the two tables, no in-figure numbers",
          "figS9_curlcake_scan.tsv" in src_r and "figS9_score_validity.tsv" in src_r
          and "stat_ecdf" in src_r and "facet_wrap" in src_r
          and "scale_y_log10" in src_r and "geom_text" not in src_r)
    check("blocks A/B/C: design AA/AA/BC, each block wrapped as one element",
          "PAGE_DESIGN  <- \"AA\\nAA\\nBC\"" in src_r
          and "PAGE_HEIGHTS <- c(1, 1, 0.75)" in src_r
          and "blk_a <- wrap_elements(wrap_plots(panels" in src_r
          and "blk_b <- wrap_elements(draw_block_b())" in src_r
          and "blk_c <- wrap_elements(draw_block_c())" in src_r
          and "wrap_plots(list(blk_a, blk_b, blk_c)" in src_r)
    txt_g = pdf_text()
    check("blocks B/C: half titles, axes and the six model labels are drawn",
          all(s in txt_g for s in ["pseU models", "m5C models",
                                   "modified-read threshold (%)",
                                   "reported modified fraction", "cumulative fraction",
                                   "FP per 10 kb"])
          and all(l in txt_g for l in vsum["label"]))
    fpr_strings = {f"{v:.5g}" for v in fpr["fp_per_10kb"]} | \
                  {f"{v:.3f}" for v in fpr["fp_per_1e6_candidates"]} | \
                  {f"{v:.4f}" for v in vsum["auc_WT_vs_IVT"]}
    check("blocks B/C: no FPR or AUC value is printed inside the figure",
          not any(s in txt_g for s in fpr_strings), f"{len(fpr_strings)} values checked")
    quoted = all(f"{v:.5g}" in legend for v in fpr["fp_per_10kb"]) and \
             all(f"{v:.4f}".rstrip("0")[:6] in legend or f"{v:.3f}" in legend
                 for v in vsum["auc_WT_vs_IVT"])
    check("legend: the reported FPR and AUC values resolve to the frozen tables", quoted)

    # 7 ------------------------------------------------------------------- #
    txt = pdf_text()
    labels = ["hac@v5.0.0_pseU", "hac@v5.1.0_pseU", "sup@v5.0.0_pseU",
              "sup@v5.1.0_pseU", "hac@v5.1.0_m5C", "sup@v5.1.0_m5C"]
    check("PDF text: the six panel titles are drawn",
          all(l in txt for l in labels))
    check("PDF text: key and axes labels are drawn",
          all(s in txt for s in ["Density", "WT", "unmodified IVT", "5'UTR",
                                 "CDS", "3'UTR", "1kb"]))
    check("PDF text: no excluded (inosine) panel is drawn", "inosine" not in txt)
    boxes_l = [b for b in word_boxes(FIGS / "FigureS9_rev.pdf")
               if len(b[4]) == 1 and b[4] in "ABCDEFGHIJ"]
    check("PDF text: the three block letters A, B, C exactly once each",
          sorted(b[4] for b in boxes_l) == ["A", "B", "C"],
          f"letters: {[b[4] for b in boxes_l]}")
    check("PDF text: no figure number is drawn on the page", "S9" not in txt)
    lb = {b[4]: b for b in boxes_l}
    check("PDF text: A/B/C sit at the top-left of their own block",
          lb["A"][0] < PAGE[0] * 0.1 and lb["A"][1] < PAGE[1] * 0.1        # page top-left
          and lb["B"][0] < PAGE[0] * 0.1 and lb["B"][1] > PAGE[1] * 0.6    # bottom, left
          and PAGE[0] * 0.45 < lb["C"][0] < PAGE[0] * 0.6                  # bottom, right
          and abs(lb["B"][1] - lb["C"][1]) < 6.0,
          "A={:.0f},{:.0f} B={:.0f},{:.0f} C={:.0f},{:.0f}".format(
              lb["A"][0], lb["A"][1], lb["B"][0], lb["B"][1],
              lb["C"][0], lb["C"][1]))

    boxes = word_boxes(FIGS / "FigureS9_rev.pdf")
    overlaps = near = 0
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i], boxes[j]
            vy = min(a[3], b[3]) - max(a[1], b[1])
            vx = min(a[2], b[2]) - max(a[0], b[0])
            hmin = min(a[3] - a[1], b[3] - b[1])
            if vy > 0.5 * hmin:
                if vx > 0.3:
                    overlaps += 1
                elif 0 < -vx < 2.0:
                    near += 1
    check("typography: no overlapping or near-colliding text on the page",
          overlaps == 0 and near == 0, f"{len(boxes)} words")

    print()
    if FAILS:
        print(f"{len(FAILS)} FAILED of {CHECKS[0]} checks: {FAILS}")
        sys.exit(1)
    print(f"all {CHECKS[0]} checks passed")


if __name__ == "__main__":
    LOGS.mkdir(parents=True, exist_ok=True)
    main()
