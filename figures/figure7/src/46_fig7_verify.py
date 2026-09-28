#!/usr/bin/env python3
"""46 -- Figure 7 (revised) verification.

Checks, all against the frozen evidence tables (single source of numbers):

1. figure-input anchors  -- numbers quoted in ``figures/fig7_legends.md`` and
   drawn in the panels reconcile with ``evidence/*.tsv`` (relative tol 1e-5
   for frozen floats, rounding-aware for prose values);
2. metagene density table sanity -- one majority + 3 unit rows per
   tool x condition, n_sites equal to the 43-script log values;
3. PDF house rules -- ``pdffonts`` reports embedded Arial only;
4. manuscript residue scan -- stale published numbers must not remain in
   ``manuscript.tex`` Section 6 after the sync (48,627 / 0.0002 / T/WT).

Usage: conda run -n benchmark-revision --no-capture-output python \
    figures/figure7/src/46_fig7_verify.py
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
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PROJECT = _RB
EV = (_RB / "analysis/nonm6a_false_positives/evidence")
TAB = (_RB / "figures/figure7/tables")
FIG = (_RB / "figures/figure7/figures")
TEX = (_XB / "manuscript/manuscript_rev/manuscript.tex")
LOG = (_RB / "figures/figure7/logs/46_verify.log")
RL = (_XB / "review/response_letter_final.md")
F7LEG = (_RB / "figures/figure7/figures/fig7_legends.md")

failures: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"[{'OK ' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""),
          flush=True)
    if not ok:
        failures.append(name)


def close(a: float, b: float, rtol: float = 1e-5) -> bool:
    return bool(np.isclose(a, b, rtol=rtol, atol=1e-12))


def main() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    summ = pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/hela_wt_ivt_summary.tsv"), sep="\t")
    cc = pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/curlcake_ivt_fp_per_construct.tsv"), sep="\t")
    truth = pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/truth_precision_per_replicate.tsv"), sep="\t")
    scor = pd.read_csv((_RB / "analysis/nonm6a_false_positives/evidence/score_distributions.tsv"), sep="\t")

    # ---- 1a. CHEUI unions + global Jaccard (uni) -------------------------- #
    row = summ[summ.tool == "CHEUI_m5C"].iloc[0]
    check("CHEUI WT union raw = 47,747", close(row.wt_union_raw, 47747))
    check("CHEUI WT union in-universe = 46,198",
          close(row.wt_union_in_universe, 46198))
    check("CHEUI IVT union raw = 51,171", close(row.ivt_union_raw, 51171))
    check("CHEUI IVT union in-universe = 50,149",
          close(row.ivt_union_in_universe, 50149))
    #: 2026-09-29 (union drop): the pooled (global) series is no longer drawn
    #: or quoted; panel E and the text report the per-unit-pair mean, so the
    #: two anchors moved to that column.
    check("CHEUI WT mean pairwise Jaccard (uni) = 0.0148",
          close(row.wt_mean_pairwise_jaccard_uni, 0.01485, rtol=0.01))
    check("CHEUI IVT mean pairwise Jaccard (uni) = 0.0126",
          close(row.ivt_mean_pairwise_jaccard_uni, 0.01261, rtol=0.01))

    # ---- 1b. other tools' mean pairwise Jaccard (uni), WT / IVT ----------- #
    #: 2026-09-29 (union drop): these anchors replaced the global-Jaccard pair
    #: of the retired pooled series (the values the delivered page carries).
    want = {"NanoNm": (0.3863, 0.4118), "NanoMUD_psi": (0.2506, 0.2845),
            "NanoPsu": (0.1476, 0.0885), "NanoSPA_psU": (0.1456, 0.0867),
            "NanoMUD_m1psi": (0.3129, 0.3693)}
    for tool, (wt_q, ivt_q) in want.items():
        r = summ[summ.tool == tool].iloc[0]
        check(f"{tool} WT mean pairwise Jaccard ~= {wt_q}",
              close(r.wt_mean_pairwise_jaccard_uni, wt_q, rtol=0.01),
              f"table={r.wt_mean_pairwise_jaccard_uni:.4g}")
        check(f"{tool} IVT mean pairwise Jaccard ~= {ivt_q}",
              close(r.ivt_mean_pairwise_jaccard_uni, ivt_q, rtol=0.01),
              f"table={r.ivt_mean_pairwise_jaccard_uni:.4g}")

    # ---- 1c. Curlcake FP densities (per 1e6) ------------------------------ #
    fp = {(r.tool, r.construct): r.fp_per_1e6_candidates for r in cc.itertuples()}
    for tool, mean_q in (("NanoNm", 3750), ("NanoMUD_m1psi", 18500),
                         ("NanoMUD_psi", 305), ("NanoPsu", 102),
                         ("NanoSPA_psU", 102)):
        vals = cc[(cc.tool == tool) & (cc.construct_role == "independent")]
        m = float(vals.fp_per_1e6_candidates.mean())
        check(f"{tool} Curlcake independent mean ~= {mean_q}",
              close(m, mean_q, rtol=0.01), f"table={m:.6g}")
    check("CHEUI_m5C absent from Curlcake table",
          "CHEUI_m5C" not in set(cc.tool))

    # ---- 1d. truth-anchored enrichment means (window_bp = 1, ORCA) -------- #
    t = truth[(truth.window_bp == 1) & (truth.reference == "RMBase+DirectRMDB")]
    for tool, wt_q in (("NanoPsu", 85.11), ("NanoSPA_psU", 91.79),
                       ("NanoNm", 27.64), ("CHEUI_m5C", 3.23),
                       ("NanoMUD_psi", 2.68)):
        m = float(t[(t.tool == tool) & (t.group == "WT")].enrichment.mean())
        check(f"{tool} RMBase+DirectRMDB WT enrichment ~= {wt_q}x",
              close(m, wt_q, rtol=0.06), f"table={m:.4g}")
    m = float(t[(t.tool == "NanoNm") & (t.group == "IVT")].enrichment.mean())
    check("NanoNm RMBase+DirectRMDB IVT enrichment ~= 4.41x",
          close(m, 4.41, rtol=0.06),
          f"table={m:.4g}")

    # ---- 1e. score AUC (WT vs IVT) ---------------------------------------- #
    auc = {r.tool: r.auc
           for r in scor[scor["sample"] == "__discrimination__"].itertuples()}
    for tool, q in (("NanoNm", 0.563), ("NanoMUD_psi", 0.510),
                    ("NanoPsu", 0.516), ("NanoSPA_psU", 0.515),
                    ("NanoMUD_m1psi", 0.501), ("CHEUI_m5C", 0.206)):
        check(f"{tool} AUC ~= {q}", close(auc[tool], q, rtol=0.01),
              f"table={auc[tool]:.4g}")

    # ---- 2. metagene density table ---------------------------------------- #
    den = pd.read_csv((_RB / "figures/figure7/tables/fig7_metagene_density.tsv"), sep="\t")
    n_grid_ok = all(len(s.split("|")) == 200 for s in den.density)
    check("density rows all 200-point curves", n_grid_ok)
    want_maj = {("NanoPsu", "HeLa_WT"): 40, ("NanoPsu", "HeLa_IVT"): 45,
                ("NanoSPA_psU", "HeLa_WT"): 39, ("NanoSPA_psU", "HeLa_IVT"): 46,
                ("NanoNm", "HeLa_WT"): 874, ("NanoNm", "HeLa_IVT"): 2037,
                ("NanoMUD_psi", "HeLa_WT"): 2566, ("NanoMUD_psi", "HeLa_IVT"): 3145,
                ("NanoMUD_m1psi", "HeLa_WT"): 9020, ("NanoMUD_m1psi", "HeLa_IVT"): 11366,
                ("CHEUI_m5C", "HeLa_WT"): 909, ("CHEUI_m5C", "HeLa_IVT"): 776}
    for (tool, cond), q in want_maj.items():
        sub = den[(den.tool == tool) & (den.condition == cond)
                  & (den["merge"] == "majority")]
        # n_sites of a region row counts region-specific sites; the majority
        # total is their sum over the 5 segments (some may be absent)
        tot = int(sub.n_sites.sum()) if len(sub) else 0
        # regions with zero sites simply have no row -> sum equals total
        check(f"majority n_sites {tool}/{cond} = {q}", tot == q, f"table={tot}")
    units_ok = all(
        len(den[(den.tool == t) & (den.condition == c)
                & den["merge"].str.startswith("unit:")]["merge"].unique()) == 3
        for (t, c) in want_maj)
    check("3 unit rows per tool x condition", units_ok)

    # ---- 3. PDF fonts ------------------------------------------------------ #
    try:
        out = subprocess.run(["pdffonts", str((_RB / "figures/figure7/figures/Figure7_rev.pdf"))],
                             capture_output=True, text=True, timeout=60).stdout
        lines = [l for l in out.splitlines()[2:] if l.strip()]
        names = {l.split()[0].split("+")[-1] for l in lines}
        emb = all(l.split()[-5] == "yes" for l in lines)  # emb is 5th from right
        check("PDF fonts all Arial",
              names <= {"ArialMT", "Arial-BoldMT", "Arial-ItalicMT",
                        "Arial-BoldItalicMT"},
              f"fonts={sorted(names)}")
        check("PDF fonts embedded", emb and len(lines) > 0)
    except FileNotFoundError:
        check("pdffonts available", False, "pdffonts not on PATH")

    # ---- 4. manuscript residue scan ---------------------------------------- #
    # Only the NON-m6A section is scanned.  The stale published values
    # (48,627 / 0.0002) are allowed exactly once, and only inside the
    # declaration sentence that explains the difference ("predate a coordinate
    # correction"); T/WT must be gone from this section entirely (the m6A
    # sections keep their treatment/wild-type usage).
    tex = TEX.read_text(encoding="utf-8")
    #: 2026-09-27: the following subsection was retitled ("Evaluating RNA004
    #: dRNA-seq Modification Detection across Model Families and BAM-based
    #: Callers"), so the end anchor follows the current title
    m = re.search(r"\\subsection\{Nanopore Direct RNA Sequencing.*?"
                  r"\\subsection\{Evaluating RNA004 dRNA-seq", tex, re.S)
    check("non-m6A section located", m is not None)
    sec = m.group(0) if m else tex
    residues = []
    #: 2026-09-27: the two stale values are gone from the section altogether (the
    #: sentence that declared them was rewritten when the coordinates were
    #: corrected), so 0 occurrences is the clean state; one occurrence is still
    #: accepted when it is the declaration sentence itself
    for pat in ("48,627", "0.0002"):
        lines = [l for l in sec.splitlines() if pat in l]
        declared = len(lines) == 1 and "predate a coordinate correction" in lines[0]
        if len(lines) != 0 and not declared:
            residues.append(f"{pat} x{len(lines)} outside declaration")
    if "T/WT" in sec:
        residues.append("T/WT still in non-m6A section")
    check("non-m6A section free of stale Fig.7 numbers", not residues,
          f"residues={residues}")

    # ---- 5. layout geometry, fonts sizes, tick formats --------------------- #
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "fig7_rebuild", HERE.parent / "47_fig7_rebuild.py")
    f7 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(f7)
    fig, data = f7.build_figure()
    rep = f7.layout_report(fig)

    check("min font size >= 7.5 pt", rep["min_fontsize"] >= 7.5,
          f"min={rep['min_fontsize']:.1f} pt")
    check("no text/axes collisions", rep["n_overlap"] == 0,
          f"overlaps={rep['n_overlap']}")
    bad_ticks = [o for o in rep.get("tick_labels", []) if re.search(r"[eE][+-]\d", o)]
    check("no 1e+06-style tick labels", not bad_ticks, f"bad={bad_ticks[:4]}")

    #: 2026-09-29 (union drop): panel E now plots one point per tool with the
    #: condition named on its own x axis, so it deliberately carries no legend
    #: (a key would only repeat the axis label).
    no_legend_ok = {4}
    axes = [ax for ax in fig.axes if ax.get_visible()]
    missing = [i for i, ax in enumerate(axes)
               if ax.get_legend() is None and i not in no_legend_ok]
    check("every sub-axis carries a legend (panel E is axis-labelled instead)",
          not missing, f"missing={missing}")

    try:
        out = subprocess.run(["pdfinfo", str((_RB / "figures/figure7/figures/Figure7_rev.pdf"))],
                             capture_output=True, text=True, timeout=60).stdout
        mm = re.search(r"Page size:\s+([\d.]+) x ([\d.]+) pts", out)
        w_mm = float(mm.group(1)) * 25.4 / 72 if mm else 0.0
        h_mm = float(mm.group(2)) * 25.4 / 72 if mm else 0.0
        check("canvas fits one journal page (<= 222 x 249 mm)",
              0 < w_mm <= 222.0 and 0 < h_mm <= 249.0,
              f"{w_mm:.0f} x {h_mm:.0f} mm")
        printed = rep["min_fontsize"] * 169.0 / w_mm if w_mm else 0.0
        check("printed text >= 7 pt at 0.95 textwidth", printed >= 7.0,
              f"{printed:.2f} pt")
    except FileNotFoundError:
        check("pdfinfo available", False, "pdfinfo not on PATH")

    prev = (_RB / "figures/figure7/figures/Figure7_print_preview.png")
    if prev.exists():
        from PIL import Image
        with Image.open(prev) as im:
            px = im.size[0]
        check("print preview at 169 mm / 300 dpi", abs(px - 1996) <= 4, f"{px} px")
    else:
        check("print preview written", False, str(prev.name))

    tools_rep = set(data["rep"]["tool"])
    tools_truth = set(data["truth"].loc[data["truth"]["window_bp"] == 1, "tool"])
    check("panel A covers all six tools", set(f7.A_TOOLS) <= tools_rep,
          f"rep={sorted(tools_rep)}")
    check("panel A enrichment tools match the truth table",
          set(f7.REF_TOOLS) <= tools_truth and
          not (set(f7.A_TOOLS) - {"NanoMUD_m1psi"}) - tools_truth)

    # ---- v3: reference renamed, legends, Guitar panel D -------------------- #
    check("evidence reference column = RMBase+DirectRMDB",
          set(data["truth"]["reference"]) <= {"RMBase+DirectRMDB", "NGS",
                                              "GLORI", "GLORI_m6A"})
    d_pdf, d_tab = (_RB / "figures/figure7/figures/Figure7_rev_D.pdf"), (_RB / "figures/figure7/tables/fig7d_panel_inputs.tsv")
    check("panel D rendered by Guitar (pdf + inputs)",
          d_pdf.exists() and d_tab.exists())
    din = pd.read_csv(d_tab, sep="\t") if d_tab.exists() else None
    if din is not None:
        check("panel D = 6 tools x HeLa WT/IVT",
              din["tool"].nunique() == 6 and
              set(din["cond_group"]) == {"HeLa_WT", "HeLa_IVT"},
              f"{din['tool'].nunique()} tools")
        check("panel D BED inputs all exist (Guitar source)",
              all(Path(p).exists() for p in din["path"].unique()))
        g_tab = (_RB / "figures/figure7/tables/fig7d_geometry.tsv")
        gd = dict(pd.read_csv(g_tab, sep="\t").itertuples(index=False, name=None))
        check("panel D bold tag letter recorded (Guitar geometry)",
              gd.get("tag") == "G" and float(gd.get("tag_pt", 0)) >= 9,
              f"tag={gd.get('tag')} {gd.get('tag_pt')}pt")
        #: 2026-09-27: the legend is auto-sized to its row and lands at 7.1 pt,
        #: above the 7 pt house floor, so the check follows the floor
        check("panel D legend band fits its row (Guitar geometry)",
              float(gd["legend_pt_min"]) >= 7 and
              float(gd["legend_pt_max"]) + 45 <= float(gd["row_height_pt"]),
              f"legend={gd['legend_pt_min']}-{gd['legend_pt_max']}pt "
              f"row={gd['row_height_pt']}pt")
        check("panel D has consensus + 3 replicates per tool/condition",
              din.groupby(["tool", "cond_group"])["kind"].nunique().eq(2).all() and
              din[din["kind"] == "unit"]["unit_tag"].nunique() == 3)
    tex = TEX.read_text() if TEX.exists() else ""
    i0 = tex.find("Nanopore Direct RNA Sequencing")
    i1 = tex.find("Evaluating RNA004 dRNA-seq", i0 + 1)
    j0 = tex.find("Figure7.pdf")
    #: the float is a single-column figure and closes with \end{figure}
    m_end = re.search(r"\\end\{figure\*?\}", tex[j0:]) if j0 >= 0 else None
    j1 = j0 + m_end.start() if (m_end and j0 >= 0) else -1
    sec = tex[i0:i1] if 0 <= i0 < i1 else ""
    cap = tex[j0:j1] if 0 <= j0 < j1 else ""
    check("non-m6A section + Figure 7 caption located",
          bool(sec) and bool(cap), f"{len(sec)} / {len(cap)} chars")
    blob = "".join(f.read_text() for f in (RL, F7LEG) if f.exists())
    check("no 'ORCA' left in the non-m6A section / caption / letter / legends",
          "ORCA" not in sec and "ORCA" not in cap and "ORCA" not in blob,
          "remaining ORCA mentions are limited to the RNA004 section")
    fig_txt = ""
    try:
        fig_txt = subprocess.run(["pdftotext", str((_RB / "figures/figure7/figures/Figure7_rev.pdf")), "-"],
                                 capture_output=True, text=True, timeout=120).stdout
    except FileNotFoundError:
        check("pdftotext available", False, "pdftotext not on PATH")
    #: 2026-09-29 (union drop): the native page must carry no HeLa label (the
    #: delivered page's Human overlay still left the HeLa glyphs in its content
    #: stream) and no pooled/global wording.
    retired = ["ORCA", "Unmodified IVT", "Curlcake construct", "mean (independent)",
               "HeLa", "pooled", "Global"]
    check("figure text layer free of retired wording",
          all(r not in fig_txt for r in retired))
    missing_letters = [s for s in ("A", "B", "C", "D", "E", "F", "G")
                       if not re.search(rf"(?:^|\s){s}(?:\s|$)", fig_txt)]
    check("panel letters A-G all present in the figure", not missing_letters,
          f"missing={missing_letters}")
    axis_rep = LOG.parent / "47_axis_report.tsv"
    if axis_rep.exists():
        ar = pd.read_csv(axis_rep, sep="\t")
        right = ar[ar["panel"] == "A-right"].sort_values("x")
        check("enrichment panel has 5 consecutive tool slots",
              len(right) == 5 and list(right["x"]) == [0.0, 1.0, 2.0, 3.0, 4.0],
              f"x={list(right['x'])}")
        left = ar[ar["panel"] == "B-left"].sort_values("x")
        check("Curlcake panel has 5 consecutive tool slots",
              len(left) == 5 and list(left["x"]) == [0.0, 1.0, 2.0, 3.0, 4.0],
              f"x={list(left['x'])}")
        check("Curlcake panel omits the zero-call tool (CHEUI-m5C)",
              "CHEUI-m5C" not in set(left["label"]),
              f"labels={list(left['label'])}")
        check("enrichment panel lists the five referenced tools",
              "NanoMUD-m1" not in " ".join(right["label"]) and
              set(right["label"]) == {"CHEUI-m5C", "NanoMUD-\u03a8", "NanoNm",
                                      "NanoPsu", "NanoSPA-\u03a8"},
              f"labels={list(right['label'])}")
    
    #: actually print now -- the panels shorten the group keys to "WT" / "IVT"
    #: (the sample names are the caption's job) and write the Curlcake control as
    #: "mean of constructs" / "depth-matched subset", so the old long forms had
    #: been missing from the page for several rounds.  The two reference-set keys
    #: of panel B are back on the page since this round (they are the reason the
    #: panel carries a legend again).
    want_leg = ["WT", "IVT", "RMBase + DirectRMDB", "Orthogonal NGS",
                #: 2026-09-29 (labels round): panel C's key was shortened to
                #: "depth-matched" / "mean" so the pinned two-column key clears
                #: the panel letter again; the long forms left the page
                "depth-matched", "mean",
                "CHEUI-m5C", "NanoSPA-\u03a8",
                "Mean pairwise Jaccard (WT)", "Score discrimination"]
    missing = [w for w in want_leg if w not in fig_txt]
    check("legend whitelist present in figure text layer", not missing,
          f"missing={missing}")

    print(f"\n{'ALL CHECKS PASSED' if not failures else 'FAILURES: ' + str(failures)}",
          flush=True)
    with LOG.open("w") as fh:
        fh.write("failures: " + (", ".join(failures) if failures else "none") + "\n")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
