#!/usr/bin/env python3
"""55 -- Figure S7 verification (rebuilt 2026-09-21, R3-9-centred contract).

Every quantity in the rebuilt Figure S7 must be traceable to
``figures/figureS7/tables/`` and reconcile with the frozen
R3-9 evidence; the page must obey the house rules.  The checks encode the
defects the author rejected on 2026-09-21, so they cannot come back:

1. tables: the six-tool / two-condition evidence set, the per-pair Jaccard table
   (90 pairs), the k-of-n table, the per-unit score table (1,800 rows) and the
   Curlcake table (five tools, 20 rows); the pooled score table must be *gone*
   (retired) because pooling replicates is not allowed;
2. anchors: ``s7_anchor_check.tsv`` has no MISMATCH, and the counts / ratios /
   Curlcake densities are re-checked here directly against the frozen R3-9
   evidence at 1e-5 relative tolerance;
3. the numbers the manuscript quotes (883 / 888 WT unions, 0.77 / 0.78 union
   ratios, CHEUI 47,747 / 51,171) are still supported, and Table S5's NanoNm row is
   anchored to the depth diagnosis of the R3-9b package -- 2026-09-30 (author decision): the
   plain mean-of-counts ratio is retired from the figure and from Table S5, so this
   script now also requires that no printed ratio exceeds one;
4. figure contract (``logs/54_layout.json``): 170 x 240 mm canvas, >= 7 pt text,
   exactly the panel letters A-C, **no tick labelled "0" on a log axis**, no
   scientific-notation ticks, every tool facet drawn with >= 3 lines per
   condition (no pooled curve), 90 pair points, **no axes title anywhere** and
   **no legend covering a data point**;
4b. 2026-10-05 (author decision): the figure is A/B/C.  Panel D was retired after five passes in two days -- a log axis of per-unit call
   counts with the two union lines, one axis of depth-matched pair marks plus a pooled
   bootstrap diamond, all six tools as Figure-3B boxes, NanoNm alone with the published
   per-site set as a band, and NanoNm alone without it.  Every version either repeated
   panel C, duplicated Figure 7A, or could not be read against a control (the published
   set is one pooled run with no unmodified control).  The row and its 3.06 in are gone,
   the three remaining rows took the height, and C carries the 45 degree tool names
   because it is the bottom row now.  The tables behind the retired cross-check are
   still verified here because Table S5's note quotes them;
5. PDF: one A4 page, ArialMT/Arial-BoldMT only and embedded, no in-figure value
   annotations, numeric tokens limited to the axis ticks;
6. residue scan: no legacy tool names, no "T/WT", no "treatment", no retired
   figure title, **no panel/facet title and no class-header label** (the house
   rule of 2026-09-19 deletes every title; the classes live in the caption) --
   in the figure text objects and in ``FigS7_legends.md``;
7. the legend file declares the no-m6A-reference basis, the unmodified-IVT
   negative-control role and the CHEUI-m5C synthetic-control gap.

Usage: conda run -n benchmark-revision --no-capture-output python \
    figures/figureS7/src/55_verify_figS7.py
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
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PROJECT = _RB
OUT = (_RB / "figures/figureS7")
TAB = (_RB / "figures/figureS7/tables")
FIG = (_RB / "figures/figureS7/figures")
LOG = (_RB / "figures/figureS7/logs")
EV = (_RB / "analysis/nonm6a_false_positives/evidence")

#: R3-9b package, so the anchor below reads the same evidence file they do
MATCHED_EV = ((_RB / "analysis/R3-9b_nonm6a_replicate_depth_20260928/evidence/band_matched_rate_ratios.tsv"))

#: row per depth-matched WT/IVT replicate pair, computed with the pooled-band weights
#: the bootstrap above resamples -- not from nearest_depth_pairs.tsv
PAIRS_EV = MATCHED_EV.with_name("band_matched_pair_ratios.tsv")
#: 2026-10-06 (author decision, second task): the two cells the delivered Table S5 prints beside the
#: ratio are the depth-standardised rate of each condition, so the gate reads them from the
#: same diagnosis (one row per unit and condition) and from the arm's own run
RATE_UNITS_EV = MATCHED_EV.with_name("band_matched_unit_rates.tsv")
#: the ratio range panel D used to draw.  Kept as the band the cross-check
#: evidence has to sit inside: the ratios are no longer drawn on the page (2026-10-04,
#: user moved D to the tool's own modification ratio), but the caption still quotes
#: the pooled intervals, so the range they are checked against stays recorded here.
RATIO_RANGE = (0.40, 1.45)

#: and Table S5's note.  Panel D drew it as a band for one pass and no longer does, so
#: these tables are verified here only because that note quotes them; the drawing itself
#: is required to carry no trace of them (see the figure-only residue below).
PUBLISHED_SCORE = ((_RB / "analysis/figS11_revision/tables/s11_published_hela_score.tsv"))
SRC_MM6 = ((_RB / "analysis/figS11_revision/sources/mmc6_fingerprint.tsv"))
XLSX = (_XB / "third_party/NanoNm/mmc6.xlsx")
#: the mean of three pair ratios and the bootstrap mean are different functionals;
#: the build refuses to draw past this gap and the gate repeats the bound
PAIR_SUMMARY_TOL = 0.12
PDF = (_RB / "figures/figureS7/figures/FigureS7_rev.pdf")
SCRIPT = HERE.parent / "54_figS7_figure.py"

TOOLS = ["CHEUI_m5C", "NanoMUD_psi", "NanoMUD_m1psi", "NanoNm", "NanoPsu",
         "NanoSPA_psU"]

#: the SI print width, which lifts the font/page ratio from 1.6 % to 2.1 %
PAGE_PT = (481.89, 680.315)
PAGE_TOL = 1.0
#: class labels of Fig. 7A -- they belong in the caption, never in the figure
CLASSES = ("FP-dominated", "Intermediate", "Sparse")
WT_UNITS = ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"]
IVT_UNITS = ["HeLa_IVT_rep1", "HeLa_IVT_rep2", "HeLa_IVT_rep3"]

failures: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"[{'ok' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""),
          flush=True)
    if not ok:
        failures.append(name)


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def close(a, b, rtol: float = 1e-5) -> bool:
    try:
        return bool(np.isclose(float(a), float(b), rtol=rtol, atol=1e-12))
    except (TypeError, ValueError):
        return False


def main() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)

    # ---- 1. evidence tables ------------------------------------------------- #
    need = {
        "s7_counts_per_replicate.tsv": 36,
        "s7_counts_summary.tsv": 12,
        "s7_ratio_ci.tsv": 6,
        "s7_jaccard_within_cross.tsv": 6,
        "s7_jaccard_pairs.tsv": 90,
        "s7_replicate_support.tsv": 36,
        "s7_score_density_per_unit.tsv": 1800,
        "s7_score_separation.tsv": 6,
        "s7_score_location_per_unit.tsv": 36,
        "s7_curlcake_per_construct.tsv": 20,
        "TableS5_reconciliation.tsv": 6,
        
        #: three samples for NanoNm (HeLa_WT2, HeLa_WT3, the published pooled set), so the
        #: cut that produced those cells is part of the evidence this gate owns: two cuts x
        #: (three unit universes + the union), the same three units and the union of this
        #: study's own call set, and the arm's three samples plus their union
        "s7_nanonm_published_wt.tsv": 16,
        #: the seven coverage strata of the published-arm ratio Table S5 prints for
        #: NanoNm (one row per stratum, both arms and the weights)
        "s7_nanonm_published_matched.tsv": 7,
        
        #: conditions restricted to the sites every unit of both conditions covers,
        #: priced for all six tools (6 tools x 2 bases)
        "s7_shared_universe.tsv": 12,
    }
    for name, rows in need.items():
        path = TAB / name
        if not path.exists():
            check(f"table {name} exists", False, "missing")
            continue
        df = pd.read_csv(path, sep="\t")
        check(f"table {name} has {rows} rows", len(df) == rows, f"{len(df)} rows")
    check("pooled score table retired (not regenerated)",
          not ((_RB / "figures/figureS7/tables/s7_score_density.tsv")).exists()
          and ((_RB / "figures/figureS7/tables/s7_score_density_RETIRED_pooled.tsv")).exists(),
          "s7_score_density.tsv must not be rewritten")
    src = SCRIPT.read_text(encoding="utf-8")
    check("figure reads the anchored S7 tables (location, not densities)",
          all(n in src for n in ("s7_score_location_per_unit.tsv",
                                 "s7_jaccard_pairs.tsv",
                                 "s7_curlcake_per_construct.tsv"))
          and "s7_score_density.tsv" not in src,
          "panel C must read the per-unit location table, never a density table")
    check("figure keeps the page gate and margin letters",
          "assert_page_clean" in src and "margin_letter" in src)
    check("the drawing code sets no axes title at all",
          "set_title" not in src and "_row_title" not in src)
    check("the retired fake-zero log tick is gone from the source",
          not any(('"0"' in ln and "10$^{" in ln)
                  for ln in src.splitlines() if "set_yticklabels" in ln))

    anchors = pd.read_csv((_RB / "figures/figureS7/tables/s7_anchor_check.tsv"), sep="\t")
    bad = anchors[anchors.status != "OK"]
    check("anchor check has no mismatch", len(bad) == 0,
          f"{len(anchors) - len(bad)}/{len(anchors)} OK"
          + (f"; first bad: {bad.iloc[0].quantity}" if len(bad) else ""))

    # ---- 2. anchors re-checked against the frozen R3-9 evidence ------------- #
    frozen = pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/hela_wt_ivt_summary.tsv"), sep="\t").set_index("tool")
    counts = pd.read_csv((_RB / "figures/figureS7/tables/s7_counts_summary.tsv"), sep="\t")
    ratios = pd.read_csv((_RB / "figures/figureS7/tables/s7_ratio_ci.tsv"), sep="\t")
    matched = pd.read_csv(MATCHED_EV, sep="\t")
    pairs = pd.read_csv((_RB / "figures/figureS7/tables/s7_jaccard_pairs.tsv"), sep="\t")
    support = pd.read_csv((_RB / "figures/figureS7/tables/s7_replicate_support.tsv"), sep="\t")
    scor = pd.read_csv((_RB / "figures/figureS7/tables/s7_score_density_per_unit.tsv"), sep="\t")
    sep = pd.read_csv((_RB / "figures/figureS7/tables/s7_score_separation.tsv"), sep="\t")
    loc = pd.read_csv((_RB / "figures/figureS7/tables/s7_score_location_per_unit.tsv"), sep="\t")
    cc = pd.read_csv((_RB / "figures/figureS7/tables/s7_curlcake_per_construct.tsv"), sep="\t")
    disc = pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/score_distributions.tsv"), sep="\t")
    disc = disc[disc["sample"] == "__discrimination__"].set_index("tool")
    frozen_cc = pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/curlcake_ivt_fp_per_construct.tsv"), sep="\t")

    ok_union = True
    for tool in TOOLS:
        row = counts[(counts.tool == tool) & (counts.condition == "WT")].iloc[0]
        if not close(row.union_raw, frozen.loc[tool, "wt_union_raw"]):
            ok_union = False
        if not close(row.union_in_universe,
                     frozen.loc[tool, "wt_union_in_universe"]):
            ok_union = False
    check("WT unions (raw and in-universe) match the frozen R3-9 summary", ok_union)

    ok_ratio = True
    for tool in TOOLS:
        r = ratios[ratios.tool == tool].iloc[0]
        c = counts[(counts.tool == tool) & (counts.condition == "IVT")].iloc[0]
        if not close(r.ratio_union_raw, frozen.loc[tool, "ratio_twt_union"]):
            ok_ratio = False
        if not close(r.ratio_union_raw,
                     float(c.union_raw)
                     / float(counts[(counts.tool == tool)
                                    & (counts.condition == "WT")].iloc[0].union_raw)):
            ok_ratio = False
    check("reference raw-union ratios match the frozen R3-9 summary", ok_ratio)

    ok_pairs = True
    for tool in TOOLS:
        d = pairs[pairs.tool == tool]
        if len(d) != 15 or set(d.comparison) != {"within_WT", "within_IVT", "cross"}:
            ok_pairs = False
        if not close(d[d.comparison == "within_WT"].jaccard_uni.mean(),
                     frozen.loc[tool, "wt_mean_pairwise_jaccard_uni"]):
            ok_pairs = False
        if not close(d[d.comparison == "within_IVT"].jaccard_uni.mean(),
                     frozen.loc[tool, "ivt_mean_pairwise_jaccard_uni"]):
            ok_pairs = False
    check("per-pair Jaccard table aggregates to the frozen pairwise means", ok_pairs)

    ok_k = True
    for (tool, cond), grp in support.groupby(["tool", "condition"]):
        if int(grp.n_sites.sum()) != int(grp.union_n.iloc[0]):
            ok_k = False
        for r in grp.itertuples():
            if not close(r.n_sites, frozen.loc[tool, f"{cond.lower()}_k{r.k}"]):
                ok_k = False
    check("k-of-n support matches the frozen replicate support", ok_k)

    ok_scor = True
    for tool in TOOLS:
        d = scor[scor.tool == tool]
        for cond, units in (("WT", WT_UNITS), ("IVT", IVT_UNITS)):
            for sample in units:
                g = d[(d.condition == cond) & (d["sample"] == sample)]
                if g.empty or int(g["count"].sum()) != int(g.n_sites.iloc[0]):
                    ok_scor = False
    check("per-unit score densities are complete and never pooled", ok_scor)

    ok_cc = False
    if set(cc.tool) == set(TOOLS) - {"CHEUI_m5C"}:
        merged = cc[cc.row_type == "construct"].merge(
            frozen_cc[["tool", "construct", "n_calls_in_universe",
                       "fp_per_1e6_candidates"]],
            on=["tool", "construct"], how="left")
        ok_cc = bool(np.allclose(merged.n_calls_in_universe_x,
                                 merged.n_calls_in_universe_y, rtol=1e-5)
                     and np.allclose(merged.fp_per_1e6_candidates_x,
                                     merged.fp_per_1e6_candidates_y, rtol=1e-5))
    check("Curlcake table (five tools, no CHEUI-m5C) matches the frozen values",
          ok_cc)

    # ---- 3. the numbers quoted by the manuscript --------------------------- #
    def get(tool: str, col: str) -> float:
        return float(counts[(counts.tool == tool)
                            & (counts.condition == "WT")].iloc[0][col])

    # 2026-09-28 (union drop) / 2026-10-06 (author decision, second task): the manuscript quotes the
    # per-unit means and SDs; from 2026-10-06 Table S5 prints the depth-standardised rate
    # cells instead, so the raw unions below are the deposited evidence's, not a column.
    quoted = {
        "NanoPsu WT mean": (get("NanoPsu", "count_mean_in_universe"), 201.0),
        "NanoPsu WT SD": (get("NanoPsu", "count_sd_in_universe"), 77.0),
        "NanoSPA_psU WT mean": (get("NanoSPA_psU", "count_mean_in_universe"), 204.0),
        "NanoSPA_psU WT SD": (get("NanoSPA_psU", "count_sd_in_universe"), 81.0),
        "CHEUI_m5C WT mean": (get("CHEUI_m5C", "count_mean_in_universe"), 15845.0),
        "CHEUI_m5C WT SD": (get("CHEUI_m5C", "count_sd_in_universe"), 3827.0),
        "CHEUI_m5C IVT mean": (float(counts[(counts.tool == "CHEUI_m5C")
                                            & (counts.condition == "IVT")]
                                     .iloc[0].count_mean_in_universe), 17129.0),
        "CHEUI_m5C IVT SD": (float(counts[(counts.tool == "CHEUI_m5C")
                                          & (counts.condition == "IVT")]
                                   .iloc[0].count_sd_in_universe), 7910.0),
        
        #: from Table S5, so the anchor is the coverage-matched ratio the D strip and
        #: the table now print; the retired value (CHEUI_m5C 1.081) stays in
        #: `s7_ratio_ci.tsv` for the deposited evidence and is deliberately not
        #: re-asserted here.
        "CHEUI_m5C coverage-matched ratio":
            (round(float(matched[matched.tool == "CHEUI_m5C"]
                         .iloc[0].band_matched_rate_ratio), 3), 0.785),
    }
    ok_quoted = True
    for label, (val, expect) in quoted.items():
        if abs(val - expect) > 1.0:
            ok_quoted = False
            print(f"       {label}: {val} != {expect}")
    check("manuscript-quoted per-unit means/SDs (201/204/15,845/17,129) hold",
          ok_quoted)
    check("the in-universe WT unions behind the panel are 488 / 497",
          close(get("NanoPsu", "union_in_universe"), 488.0)
          and close(get("NanoSPA_psU", "union_in_universe"), 497.0))
    
    #: printed for the unmodified control, in the strip or in Table S5
    check("no coverage-matched ratio exceeds one",
          bool((matched["band_matched_rate_ratio"] <= 1.0).all()),
          str(matched[["tool", "band_matched_rate_ratio"]].round(3)
              .to_dict("records")))
    
    #: NanoNm -- two of this study's units merged with the published per-site set.  Five
    #: things have to hold at once, and each is a defect the author would otherwise have to
    #: catch by eye: the recorded union of the arm is the union of its samples the evidence
    #: holds, the published sample is within a factor of three of this study's own
    #: re-analysis (so the note may call it a reproduction), the published set really is a
    #: member of the arm the mean and SD are taken over, the third unit of this study is
    #: *out* of that arm (the author asked for the smallest-count sample to be dropped, so a
    #: rebuild may not silently put it back), and the unmodified-IVT side stays this study's
    #: own numbers, because the published set has no control and nothing about it may leak
    #: into the control columns.
    nano_pub = pd.read_csv((_RB / "figures/figureS7/tables/s7_nanonm_published_wt.tsv"), sep="\t")
    pub_cell = nano_pub[(nano_pub.source.str.startswith("published"))
                        & (nano_pub.cut == "matched_ratio_support")
                        & (nano_pub.unit == "union_of_the_three_universes")]
    own_cell = nano_pub[nano_pub.unit == "union_of_the_three_units"]
    arm_cell = nano_pub[(nano_pub.source == "merged wild-type arm")
                        & (nano_pub.unit == "union_of_the_three_samples")]
    arm_samples = nano_pub[(nano_pub.source == "merged wild-type arm")
                           & (nano_pub.unit != "union_of_the_three_samples")]
    anchors = pd.read_csv((_RB / "figures/figureS7/logs/76_anchors.tsv"), sep="\t").set_index("key")["value"]
    recon = pd.read_csv((_RB / "figures/figureS7/tables/TableS5_reconciliation.tsv"), sep="\t")
    recon_nano = recon[recon.tool == "NanoNm"]
    #: 2026-10-06: the table prints the arm's *rate* cells now, so the union this check reads
    #: is the arm's own union as the evidence holds it, not a printed cell any more
    check("Table S5: the NanoNm arm's union (no longer printed) is the three-sample union",
          len(pub_cell) == 1 and len(own_cell) == 1 and len(arm_cell) == 1
          and int(recon_nano.iloc[0].t5_wt_union) == int(arm_cell.iloc[0].sites_in_universe)
          and int(recon_nano.iloc[0].s7_wt_union_raw) == 4181,
          f"recorded {int(recon_nano.iloc[0].t5_wt_union)}, arm "
          f"{int(arm_cell.iloc[0].sites_in_universe)}, our own re-analysis "
          f"{int(own_cell.iloc[0].sites_in_universe)}")
    check("Table S5: the arm is HeLa_WT2 + HeLa_WT3 + the published set",
          len(arm_samples) == 3
          and set(arm_samples.sites_in_universe.astype(int))
          == {3433, 1648, int(pub_cell.iloc[0].sites_in_universe)}
          and str(anchors["merged_wt_units_excluded"]) == "HeLa_WT1"
          and round(float(arm_samples.sites_in_universe.mean()), 1)
          == float(anchors["merged_wt_mean_per_unit"])
          and round(float(arm_samples.sites_in_universe.std(ddof=1)), 1)
          == float(anchors["merged_wt_sd_per_unit"]),
          f"samples {sorted(arm_samples.sites_in_universe.astype(int).tolist())}, "
          f"excluded {anchors['merged_wt_units_excluded']}")
    check("Table S5: the published NanoNm set reproduces our own wild-type count",
          1.0 / 3.0 <= float(pub_cell.iloc[0].sites_in_universe)
          / float(own_cell.iloc[0].sites_in_universe) <= 3.0,
          f"{float(pub_cell.iloc[0].sites_in_universe):.0f} against "
          f"{float(own_cell.iloc[0].sites_in_universe):.0f} sites "
          f"({float(pub_cell.iloc[0].rate_per_1e6):.0f} against "
          f"{float(own_cell.iloc[0].rate_per_1e6):.0f} per million candidate sites)")
    check("Table S5: NanoNm's control side stays this study's own",
          int(recon_nano.iloc[0].t5_ivt_union) == int(recon_nano.iloc[0].s7_ivt_union_raw)
          == 6399 and close(recon_nano.iloc[0].s7_ratio_raw, 6399 / 4181)
          and "has no unmodified control" in recon_nano.iloc[0].note,
          f"IVT union {int(recon_nano.iloc[0].t5_ivt_union)}, count ratio "
          f"{float(recon_nano.iloc[0].s7_ratio_raw):.4f}")
    
    #: the four samples of the arm.  The gate recomputes the weighted estimator from the
    #: per-stratum evidence, requires the reconciliation to record that basis, requires the
    #: note to carry the within-study value beside it (so the row cannot be read as if all
    #: six rows came from one arm), and requires the arm to stay under one -- a value above
    #: one would put the row in the same direction as the retired mean-of-counts artefact.
    band = pd.read_csv((_RB / "figures/figureS7/tables/s7_nanonm_published_matched.tsv"), sep="\t")
    recomputed = float((band.pooled_weight * band.ivt_this_study_rate_per_1e6).sum()
                       / (band.pooled_weight * band.wt_merged_arm_rate_per_1e6).sum())
    check("Table S5: the printed NanoNm ratio is the arm's weighted estimator",
          len(band) == 7
          and abs(float(recon_nano.iloc[0].matched_ratio_printed) - recomputed) < 5e-4
          and abs(recomputed - float(anchors["merged_ratio_matched"])) < 5e-4
          and str(recon_nano.iloc[0].matched_ratio_basis)
          == "this study's HeLa_WT2 and HeLa_WT3 units plus the published pooled set"
          and 0.0 < recomputed < 1.0,
          f"printed {float(recon_nano.iloc[0].matched_ratio_printed):.4f}, recomputed "
          f"{recomputed:.4f}, basis {recon_nano.iloc[0].matched_ratio_basis}")
    #: 2026-10-06 (author decision, option C): the delivered Table S5 prints a counting interval for
    #: every row, so this row's has to be the counting one: named as such in the anchors,
    #: bracketing the printed point estimate, and narrower than the unit-level interval that
    #: stays on the record for comparison.  A swap of the two would print a second
    #: convention beside the five rows and would not be visible from the cell alone.
    check("Table S5: that ratio's interval is the counting (site-level) one",
          float(anchors["merged_ratio_matched_ci_lo"]) < recomputed
          < float(anchors["merged_ratio_matched_ci_hi"])
          and str(anchors["merged_ratio_matched_ci_convention"]).startswith("site-level")
          and float(anchors["merged_ratio_matched_ci_hi"])
          < float(anchors["merged_ratio_matched_ci_hi_unit_level"])
          and "counting" in recon_nano.iloc[0].note,
          f"interval {float(anchors['merged_ratio_matched_ci_lo']):.3f}-"
          f"{float(anchors['merged_ratio_matched_ci_hi']):.3f}, unit-level "
          f"{float(anchors['merged_ratio_matched_ci_hi_unit_level']):.3f}")
    check("Table S5: the same row still carries the within-study ratio in its note",
          "0.962 (0.931-0.994)" in recon_nano.iloc[0].note,
          "note must name the arm the other five rows use")
    
    #: page any more -- panel D was retired, so Figure S8 carries A-C only -- but Table
    #: S5's note and the manuscript point at Figure S11, which is built from these same
    #: tables, so they are still checked here: 6 tools x 3 depth-matched pairs, every
    #: value finite, and the pooled intervals finite as well.
    pair_ratios = pd.read_csv(PAIRS_EV, sep="\t")
    
    #: parity, the same presentation requirement as "no coverage-matched ratio above one"
    #: for the other five.  It is the counting convention that makes it hold: the two
    #: unit-level intervals kept in the anchors are wider, and the all-resampled one crosses
    #: one, so a rebuild that silently reverted the convention would fail here.
    check("Table S5: the NanoNm interval stays below parity",
          float(anchors["merged_ratio_matched_ci_hi"]) < 1.0
          < float(anchors["merged_ratio_matched_ci_hi_all_resampled"]),
          f"printed upper bound {float(anchors['merged_ratio_matched_ci_hi']):.3f}, "
          f"all-resampled {float(anchors['merged_ratio_matched_ci_hi_all_resampled']):.3f}")
    
    #: The gate owns the shape, the reproduction of the delivered column from the depth
    #: diagnosis, and the equality of the shared column with that run's anchors.  It does
    #: *not* require the shared ratios to stay below one: on that site space NanoNm's does not
    #: (1.036), which is why the numbers are deposited evidence and only an aggregate sentence
    #: is printable.
    shared = pd.read_csv((_RB / "figures/figureS7/tables/s7_shared_universe.tsv"), sep="\t")
    a77 = pd.read_csv((_RB / "figures/figureS7/logs/77_anchors.tsv"), sep="\t").set_index("key")["value"]
    own = shared[shared.basis == "own_candidate_sets"].set_index("tool")
    sh = shared[shared.basis == "shared_universe"].set_index("tool")
    check("shared universe: 12 rows, every tool on both bases",
          sorted(own.index.tolist()) == sorted(TOOLS) == sorted(sh.index.tolist()),
          f"{sorted(own.index.tolist())}")
    check("shared universe: the delivered column reproduces the depth diagnosis",
          all(abs(float(own.loc[t].matched_ratio)
                  - float(matched[matched.tool == t].band_matched_rate_ratio.iloc[0])) < 5e-4
              for t in TOOLS),
          str(own.matched_ratio.to_dict()))
    check("shared universe: the shared column equals the run's anchors",
          all(abs(float(sh.loc[t].matched_ratio)
                  - float(a77[f"{t}_{sh.loc[t].mod_type}_shared_match"])) < 5e-4
              for t in TOOLS),
          str(sh.matched_ratio.to_dict()))
    check("shared universe: every tool moves toward parity on the shared site space",
          all(float(sh.loc[t].matched_ratio) >= float(own.loc[t].matched_ratio) - 0.03
              for t in TOOLS)
          and float(sh.loc["NanoNm"].matched_ratio) > 1.0,
          "NanoNm is the only tool the shared site space pushes above one")
    check("cross-check evidence: 18 depth-matched pair ratios, all finite",
          sorted(pair_ratios.tool.unique().tolist()) == sorted(TOOLS)
          and len(pair_ratios) == 18
          and bool(np.isfinite(pair_ratios["matched_rate_ratio_ivt_over_wt"]).all()),
          str(pair_ratios.groupby("tool").size().to_dict()))
    check("cross-check evidence: every pooled interval is finite",
          bool(np.isfinite(matched["boot_matched_ci_lo"]).all()
               and np.isfinite(matched["boot_matched_ci_hi"]).all()
               and (matched["boot_matched_ci_hi"]
                    >= matched["boot_matched_ci_lo"]).all()),
          f"intervals {matched['boot_matched_ci_lo'].min():.3f}-"
          f"{matched['boot_matched_ci_hi'].max():.3f}")
    #: 2026-10-06 (author decision, option C): Table S5 prints the counting interval of each of the six
    #: ratios, so the depth diagnosis has to carry one for every tool, it has to bracket the
    #: printed point estimate, and -- the requirement the author set for this table -- every one
    #: of them has to stay below parity, which is what the change of convention buys.
    check("cross-check evidence: the counting interval of every tool brackets its ratio",
          bool(np.isfinite(matched["site_boot_ci_lo"]).all()
               and np.isfinite(matched["site_boot_ci_hi"]).all()
               and (matched["site_boot_ci_lo"] <= matched["band_matched_rate_ratio"]).all()
               and (matched["band_matched_rate_ratio"] <= matched["site_boot_ci_hi"]).all()),
          f"intervals {matched['site_boot_ci_lo'].min():.3f}-"
          f"{matched['site_boot_ci_hi'].max():.3f}")
    check("cross-check evidence: every printed interval stays below parity",
          bool((matched["site_boot_ci_hi"] < 1.0).all()),
          str(matched[["tool", "site_boot_ci_hi"]].round(3).to_dict("records")))
    #: 2026-10-06 (author decision, second task): the two count columns beside the ratio are the
    #: depth-standardised rate of each condition now, and the contract the author set for them
    #: is that the wild type reads *above* the unmodified control in every row and that the
    #: ratio column is their quotient.  The gate recomputes all three from the evidence --
    #: per-unit rates for the five study rows, the 76 anchors for NanoNm's arm -- and requires
    #: the reconciliation to have recorded the same two cells, so a rebuild cannot print a
    #: rate pair that no evidence file holds.
    rates = pd.read_csv(RATE_UNITS_EV, sep="\t")
    t5_anchors = pd.read_csv((_RB / "figures/figureS7/logs/76_anchors.tsv"), sep="\t").set_index("key")["value"]
    rate_ok, rate_msg = True, []
    for _, row in recon.iterrows():
        if row.tool == "NanoNm":
            wt_r = float(t5_anchors["merged_wt_rate_mean_per_1e6"])
            ivt_r = float(t5_anchors["ivt_rate_mean_per_1e6"])
        else:
            g = rates[rates.tool == row.tool]
            wt_r = float(g[g.condition == "WT"].band_matched_rate_per_1e6.mean())
            ivt_r = float(g[g.condition == "IVT"].band_matched_rate_per_1e6.mean())
        ok = (wt_r > ivt_r
              and abs(round(wt_r, 1) - float(row.printed_wt_rate)) < 0.051
              and abs(round(ivt_r, 1) - float(row.printed_ivt_rate)) < 0.051
              and abs(ivt_r / wt_r - float(row.matched_ratio_printed)) < 1e-3)
        rate_ok = rate_ok and ok
        rate_msg.append(f"{row.tool} {wt_r:.0f}/{ivt_r:.0f} -> {ivt_r / wt_r:.3f}"
                        f"{'' if ok else ' MISMATCH'}")
    check("Table S5: the printed rate cells put the wild type above the control and are the"
          " ratio's two sides", bool(rate_ok), "; ".join(rate_msg))
    #: the published distribution Figure S11 draws for NanoNm, and the source it was
    #: parsed from: the table must exist, its quantiles must be ordered, and the
    #: workbook's md5 must still be the one recorded when it was parsed.  Table S5's
    #: note quotes this cross-check, which is why the S7 gate carries it.
    pub = PUBLISHED_SCORE
    if not pub.exists():
        check("published NanoNm distribution table present", False, str(pub))
    else:
        prow = pd.read_csv(pub, sep="\t")
        need_q = ["n", "min", "q05", "q1", "median", "q3", "q95", "max"]
        ok_q = False
        detail = str(list(prow.columns))
        if len(prow) == 1 and all(c in prow.columns for c in need_q):
            v = {c: float(prow[c].iloc[0]) for c in need_q}
            ok_q = (v["min"] <= v["q05"] <= v["q1"] <= v["median"] <= v["q3"]
                    <= v["q95"] <= v["max"] <= 1.0
                    and str(prow["label"].iloc[0]) == "HeLa"
                    and str(prow["source"].iloc[0]) == "published")
            detail = (f"n={v['n']:.0f} q05={v['q05']:.3f} q1={v['q1']:.3f} "
                      f"med={v['median']:.3f} q3={v['q3']:.3f} q95={v['q95']:.3f}")
        check("published NanoNm distribution: one ordered HeLa row inside [0, 1]",
              ok_q, detail)
        fp_rows = pd.read_csv(SRC_MM6, sep="\t", comment=None) \
            if SRC_MM6.exists() else None
        if fp_rows is None or "md5" not in fp_rows.columns:
            check("published table still carries its source fingerprint", False,
                  str(SRC_MM6))
        else:
            src = fp_rows[fp_rows.role == "source"].iloc[0]
            check("published table still carries its source fingerprint",
                  src.md5 == md5(XLSX), f"{src.md5} vs {md5(XLSX)}")

    # ---- 4. figure contract ------------------------------------------------- #
    rep_path = (_RB / "figures/figureS7/logs/54_layout.json")
    if not rep_path.exists():
        check("layout report exists", False, str(rep_path.name))
    else:
        rep = json.loads(rep_path.read_text())
        w, h = rep["page_pt"]
        check("page is the compact 170 x 240 mm canvas",
              abs(w - PAGE_PT[0]) < PAGE_TOL and abs(h - PAGE_PT[1]) < PAGE_TOL,
              f"{w} x {h} pt (expected {PAGE_PT})")
        check("smallest text >= 7 pt", rep["min_font_pt"] >= 7.0,
              f"{rep['min_font_pt']} pt")
        panels = {p["panel"]: p for p in rep["panels"]}
        
        #: Panel D went through five passes and every version either repeated panel C,
        #: duplicated Figure 7A, or could not be read against a control (a published
        #: set with no unmodified control of its own).  The gate therefore requires the
        #: row to be gone, requires C to have taken the 45 degree tool names with it
        #: (it is the bottom row now), and keeps the retired-element scan for the whole
        #: page so none of D's five designs can come back through the layout report.
        check("panels A-C are all present and the retired D row is gone",
              set(panels) == {"A", "B", "C"} and "D" not in panels,
              f"{sorted(panels)}")
        check("C carries the tool names (it is the bottom row now)",
              panels["C"].get("kind") == "score_location_strip"
              and float(panels["C"].get("xlabel_band_in", 0.0)) > 0.8,
              f"xlabel band {panels['C'].get('xlabel_band_in')} in")
        check("the retired D designs are absent from every panel",
              not ({"counts_ylim", "pair_ylim", "ratio_ylim", "pair_marks", "ylim_d",
                    "pooled_marks", "axes_count", "n_boxes", "units_per_box",
                    "published_band", "published_n", "published_bands", "unit_medians_d",
                    "n_our_units"}
                   & set().union(*(set(p) for p in panels.values()))),
              str({k: sorted(set(p) - {"panel"}) for k, p in panels.items()}))
        check("C: score-location strip, six tools, 36 unit marks",
              panels["C"].get("kind") == "score_location_strip"
              and panels["C"].get("n_marks") == 36
              and panels["C"].get("marks_per_tool") == 6,
              str({k: v for k, v in panels["C"].items()
                   if k in ("kind", "n_marks", "marks_per_tool")}))
        check("B: 90 pair points, none at the floor",
              panels["B"]["n_pair_points"] == 90
              and panels["B"]["n_zero_pairs"] == 0)
        check("A: five tools with constructs, no CHEUI-m5C column",
              panels["A"]["n_positive_points"] == 12
              and panels["A"]["n_zero_points"] == 3
              and panels["A"].get("n_tools") == 5
              and panels["A"].get("gap_marker") is False,
              str({k: v for k, v in panels["A"].items()
                   if k in ("n_tools", "gap_marker")}))
        logs = rep["log_axes"]
        check("log axes carry no tick labelled 0",
              all("0" not in a["ytick_labels"] for a in logs),
              str([a["ytick_labels"] for a in logs]))
        check("no scientific-notation tick labels",
              not any(re.search(r"\de[+-]?\d", t)
                      for a in logs for t in a["ytick_labels"]))
        check("no axes carries a title", rep["n_titles"] == 0,
              f"n_titles={rep['n_titles']}")
        check("no legend covers a data point",
              not rep["legend_data_violations"],
              str(rep["legend_data_violations"])[:120])
        check("no legend overflows its own panel",
              not rep.get("legend_out_of_panel"),
              str(rep.get("legend_out_of_panel"))[:120])
        check("C: every unit's median lies inside its interquartile range",
              all(float(r.q25) <= float(r.median) <= float(r.q75)
                  for r in loc.itertuples()),
              f"{len(loc)} unit rows")
        saturated = {t for t in loc.tool
                     if float(loc[loc.tool == t]["frac_ge_0p9"].min()) >= 1.0}
        check("C: the four saturated tools are flagged and pinned at the ceiling",
              saturated == {"NanoMUD_psi", "NanoMUD_m1psi", "NanoPsu",
                            "NanoSPA_psU"}
              and set(panels["C"].get("saturated_tools", [])) == saturated
              and all(abs(float(loc[loc.tool == t]["median"].min()) - 1.0) < 1e-9
                      for t in ("NanoMUD_psi", "NanoMUD_m1psi")),
              f"saturated={sorted(saturated)}")
        check("separation table: pair range ordered, pooled column = Fig. 7F",
              all(float(r.pair_auc_min) <= float(r.pair_auc_mean)
                  <= float(r.pair_auc_max) for r in sep.itertuples())
              # the table stores %.6g, so compare at the 1e-5 anchor tolerance
              and all(abs(float(r.auc_pooled)
                          - float(disc.loc[r.tool, "auc"])) <= 1e-5
                      for r in sep.itertuples()),
              str({r.tool: (round(float(r.pair_auc_min), 4),
                            round(float(r.pair_auc_mean), 4),
                            round(float(r.pair_auc_max), 4))
                   for r in sep.itertuples()}))

    # ---- 5. PDF / PNG ------------------------------------------------------- #
    try:
        info = subprocess.run(["pdfinfo", str(PDF)], capture_output=True,
                              text=True, check=True).stdout
        pages = re.search(r"Pages:\s+(\d+)", info)
        size = re.search(r"Page size:\s+([\d.]+) x ([\d.]+)", info)
        check("PDF is one page of the compact canvas",
              pages is not None and int(pages.group(1)) == 1
              and size is not None
              and abs(float(size.group(1)) - PAGE_PT[0]) < PAGE_TOL
              and abs(float(size.group(2)) - PAGE_PT[1]) < PAGE_TOL,
              size.group(0) if size else "no page size")
    except FileNotFoundError:
        check("pdfinfo available", False, "pdfinfo not on PATH")

    fonts = subprocess.run(["pdffonts", str(PDF)], capture_output=True,
                           text=True).stdout.splitlines()[2:]
    rows = [ln.split() for ln in fonts if ln.strip()]
    names = {r[0].split("+")[-1] for r in rows}          # drop the subset prefix
    check("fonts are Arial only and embedded",
          bool(rows) and all(n in {"ArialMT", "Arial-BoldMT"} for n in names)
          and all(r[-5] == "yes" for r in rows), f"{sorted(names)}")

    txt = subprocess.run(["pdftotext", "-layout", str(PDF), "-"],
                         capture_output=True, text=True, check=True).stdout
    # The panel letters sit in the left page margin (x < 12 pt).  They cannot be
    # found in the text flow: pdftotext splits the 45-degree tool labels into
    # single glyphs, so "NanoMUD-m1Psi" contributes a stray "D" and "CHEUI-m5C"
    # a stray "C".  The bounding boxes separate them unambiguously.
    xml = subprocess.run(["pdftotext", "-bbox", str(PDF), "-"],
                         capture_output=True, text=True, check=True).stdout
    words = re.findall(r'<word xMin="([\d.]+)"[^>]*>([^<]*)</word>', xml)
    margin = sorted(t for x, t in words if float(x) < 12.0)
    letters = {L: margin.count(L) for L in "ABC"}
    # 2026-09-28: the rotated axis title of the counts panel ("Calls in the
    # candidate-site set") sits in the same left band as the panel letters in
    # the delivered page, so the band is checked as "the letters exactly once plus
    # at most the known axis-title words", never as a count of three.
    
    # the whitelist again and the retired letter D must not appear at all.
    AXIS_TITLE_WORDS = ("Calls", "in", "the")
    stray = [w for w in margin if w not in "ABC"]
    check("panel letters A-C appear once each in the page margin",
          all(n == 1 for n in letters.values())
          and not [w for w in stray if w not in AXIS_TITLE_WORDS],
          f"{letters} (margin words: {margin})")
    check("no retired or further panel letter is left over",
          "D" not in margin and "E" not in margin)
    titles = ["Unmodified Curlcake controls",
              "Overlap within and between conditions",
              "Score distribution per sequencing unit",
              "HeLa calls and the unmodified-control ratio"]
    left = [t for t in titles + list(CLASSES) if t in txt]
    check("no title and no class label survives in the figure", not left,
          str(left))

    numeric = re.findall(r"\d+(?:[.,]\d+)?", txt)
    suspicious = [n for n in numeric
                  if ("," in n) or (len(n.split(".")[-1]) > 2 if "." in n
                                    else len(n) >= 4)]
    check("no value annotations inside the panels", not suspicious,
          str(suspicious[:6]))
    check("figure text names both conditions",
          "Within WT" in txt and "Within IVT" in txt)

    png = (_RB / "figures/figureS7/figures/FigureS7_rev.png")
    if png.exists():
        from PIL import Image
        with Image.open(png) as im:
            px, dpi = im.size, im.info.get("dpi", (0, 0))
        check("PNG is 300 dpi at the compact page size",
              abs(px[0] - round(PAGE_PT[0] / 72 * 300)) <= 4
              and int(round(dpi[0])) == 300,
              f"{px[0]} px, dpi={dpi}")
    else:
        check("PNG written", False, png.name)
    prev = (_RB / "figures/figureS7/figures/FigureS7_print_preview.png")
    if prev.exists():
        from PIL import Image
        with Image.open(prev) as im:
            check("print preview at 169 mm / 300 dpi",
                  abs(im.size[0] - 1996) <= 4, f"{im.size[0]} px")
    else:
        check("print preview written", False, prev.name)

    # ---- 6. naming / wording residue --------------------------------------- #
    banned = {
        "NanoSPA_Psu": re.compile(r"NanoSPA_Psu"),
        "NanoMUD_m1psi": re.compile(r"NanoMUD_m1psi"),
        "NanoMUD_psi": re.compile(r"NanoMUD_psi"),
        "T/WT": re.compile(r"(?<![A-Za-z])T/WT"),
        "treatment": re.compile(r"treatment", re.I),
        "Stoichiometric analysis": re.compile(r"Stoichiometric analysis"),
        
        #: page or in the legend file -- the strip and Table S5 print the matched ratio
        "mean-of-counts": re.compile(r"mean-of-counts", re.I),
        
        
        #: 2026-10-04 (author decision): the published per-site set is no longer drawn
        #: on the page; the cross-check that quotes it is Table S5's note.  A key or axis
        #: label naming it would be a leftover of the retired band.
        "published (figure only)": re.compile(r"published"),
    }
    #: wording that belongs in the legend file (which must point the reader at the
    #: published-set evidence) but never in the drawing itself
    figure_only = {"published (figure only)"}
    residues = []
    for label, pat in banned.items():
        if pat.search(txt):
            residues.append(f"figure: {label}")
    legend_md = (_RB / "figures/figureS7/figures/FigS7_legends.md")
    body = legend_md.read_text(encoding="utf-8") if legend_md.exists() else ""
    for label, pat in banned.items():
        if label in figure_only:
            continue
        if pat.search(body):
            residues.append(f"legends.md: {label}")
    check("no legacy naming / retired wording", not residues, str(residues))

    # ---- 7. legend file ----------------------------------------------------- #
    # collapse the line wrapping: a phrase must be found whatever the reflow
    # (the CHEUI-m5C "never run" clause was split across two lines on 2026-09-21
    # and the substring check failed although the sentence was there)
    low = re.sub(r"\s+", " ", body.lower())
    check("legend file declares the no-m6A-reference basis",
          "no m6a-centred reference" in low or "no m6a reference" in low)
    check("legend file frames the IVT libraries as unmodified controls",
          "unmodified ivt" in low and "negative control" in low)
    check("legend file states the CHEUI-m5C synthetic-control gap",
          "cheui" in low and ("not run" in low or "never run" in low
                              or "no construct run" in low))
    check("legend file documents panels A-C and drops the retired D",
          sum(1 for p in "abc" if f"({p})" in low) == 3 and "(d)" not in low)

    print(f"\n{'ALL CHECKS PASSED' if not failures else 'FAILURES: ' + str(failures)}",
          flush=True)
    with LOG.joinpath("55_verify.log").open("w") as fh:
        fh.write("failures: " + (", ".join(failures) if failures else "none") + "\n")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
