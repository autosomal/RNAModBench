#!/usr/bin/env python3
"""55 -- Figure S6 verification (rebuilt 2026-09-21, R3-9-centred contract).

Every quantity in the rebuilt Figure S6 must be traceable to
``figures/figureS7/tables/`` and reconcile with the frozen
R3-9 evidence; the page must obey the house rules.  The checks encode the
defects the user rejected on 2026-09-21, so they cannot come back:

1. tables: the six-tool / two-condition evidence set, the per-pair Jaccard table
   (90 pairs), the k-of-n table, the per-unit score table (1,800 rows) and the
   Curlcake table (five tools, 20 rows); the pooled score table must be *gone*
   (retired) because pooling replicates is not allowed;
2. anchors: ``s7_anchor_check.tsv`` has no MISMATCH, and the counts / ratios /
   Curlcake densities are re-checked here directly against the frozen R3-9
   evidence at 1e-5 relative tolerance;
3. the numbers the manuscript quotes (883 / 888 WT unions, 0.77 / 0.78 union
   ratios, CHEUI 47,747 / 51,171) are still supported;
4. figure contract (``logs/54_layout.json``): 170 x 240 mm canvas, >= 7 pt text,
   exactly the panel letters A-D, **no tick labelled "0" on a log axis**, no
   scientific-notation ticks, every tool facet drawn with >= 3 lines per
   condition (no pooled curve), 90 pair points, **no axes title anywhere** and
   **no legend covering a data point**;
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
PDF = (_RB / "figures/figureS7/figures/FigureS6_rev.pdf")
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
    }
    for name, rows in need.items():
        path = TAB / name
        if not path.exists():
            check(f"table {name} exists", False, "missing")
            continue
        df = pd.read_csv(path, sep="\t")
        check(f"table {name} has {rows} rows", len(df) == rows, f"{len(df)} rows")
    check("pooled score table retired (not regenerated)",
          not ((_RB / "figures/figureS7/tables/s6_score_density.tsv")).exists()
          and ((_RB / "figures/figureS7/tables/s7_score_density_RETIRED_pooled.tsv")).exists(),
          "s6_score_density.tsv must not be rewritten")
    src = SCRIPT.read_text(encoding="utf-8")
    check("figure reads the anchored S6 tables (location, not densities)",
          all(n in src for n in ("s7_score_location_per_unit.tsv",
                                 "s7_jaccard_pairs.tsv",
                                 "s7_curlcake_per_construct.tsv"))
          and "s6_score_density.tsv" not in src,
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

    # 2026-09-28 (union drop): the manuscript now quotes the per-unit means and
    # SDs; the raw unions are reference columns of Table S5 only.
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
        "CHEUI_m5C mean-of-counts ratio": (round(float(ratios[ratios.tool == "CHEUI_m5C"]
                                                        .iloc[0].ratio_mean_counts), 3), 1.081),
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
        check("panels A-D are all present",
              set(panels) == {"A", "B", "C", "D"}, f"{sorted(panels)}")
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
    letters = {L: margin.count(L) for L in "ABCD"}
    # 2026-09-28: the rotated axis title of the counts panel ("Calls in the
    # candidate-site set") sits in the same left band as the panel letters in
    # the delivered page, so the band is checked as "the four letters exactly
    # once plus at most the known axis-title words", never as a count of four.
    stray = [w for w in margin if w not in "ABCD"]
    check("panel letters A-D appear once each in the page margin",
          all(n == 1 for n in letters.values())
          and not [w for w in stray if w not in ("Calls", "in", "the")],
          f"{letters} (margin words: {margin})")
    check("no fifth panel letter is left over", "E" not in margin)
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

    png = (_RB / "figures/figureS7/figures/FigureS6_rev.png")
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
    prev = (_RB / "figures/figureS7/figures/FigureS6_print_preview.png")
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
    }
    residues = []
    for label, pat in banned.items():
        if pat.search(txt):
            residues.append(f"figure: {label}")
    legend_md = (_RB / "figures/figureS7/figures/FigS7_legends.md")
    body = legend_md.read_text(encoding="utf-8") if legend_md.exists() else ""
    for label, pat in banned.items():
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
    check("legend file documents panels A-D",
          sum(1 for p in "abcd" if f"({p})" in low) == 4)

    print(f"\n{'ALL CHECKS PASSED' if not failures else 'FAILURES: ' + str(failures)}",
          flush=True)
    with LOG.joinpath("55_verify.log").open("w") as fh:
        fh.write("failures: " + (", ".join(failures) if failures else "none") + "\n")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
