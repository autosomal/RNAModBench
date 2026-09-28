#!/usr/bin/env python3
"""Assemble Supplementary Tables S1-S9 for the revision (2026-09-21).

Every cell is read from a frozen table of the delivered figures or from the
registered analysis outputs -- nothing is recomputed here; the few curated text
tables (S1, S6-S9) are transcribed from the delivered supplementary legend /
response-letter wording and are marked as curated in `README.md`.

Outputs (this directory)
------------------------
tables/TableS1_tools.tsv ... TableS9_evaluation_boundary.tsv   machine-readable
tables/_index.tsv                                              what each table is

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/tables/build_supp_tables.py
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
import csv
import re
from pathlib import Path

import pandas as pd

BENCH = Path(str(_RB))
RA = (_RB / "analysis")
OUT = (_RB / "tables/tables")
OUT.mkdir(parents=True, exist_ok=True)

MAXCELL = 150  # characters kept per cell in the printed per-tool table


#: 2026-09-26 (user, de-AI pass): two curated cells of the deposited evidence CSV carry an
#: em dash.  The printed table uses plain punctuation instead, so the delivered SI carries
#: none; the evidence CSV is deliberately left byte-identical, because it is deposited with
#: the source code.  The assertion in `write_table` fails if any other dash ever appears.
DASH_FIX = (
    ("none — output has no p-value/q-value column",
     "none: output has no p-value/q-value column"),
    ("0.2 (setup.py) — editable install of the vendored tree",
     "0.2 (setup.py), editable install of the vendored tree"),
    
    #: curated cells reach the printed page; the source tables stay byte-identical
    #: (they are deposited) and the SI prints plain punctuation.
    ("not applicable — comparative stats path (no learned checkpoint)",
     "not applicable: comparative stats path (no learned checkpoint)"),
    ("none — GMM fitted per position at runtime",
     "none: GMM fitted per position at runtime"),
    ("none — statistical tool (no model file, no --model option)",
     "none: statistical tool (no model file, no --model option)"),
)


#: 2026-09-27: the tool is Yanocomp, and two of the upstream tables (TableS2, TableS3) carry
#: it lower-cased.  The printed table is normalised; the evidence files stay as they are.
CASE_FIX = (("yanocomp", "Yanocomp"),)



#: not everything that happens to sit in the working tree.  The manuscript
#: benchmarks 15 dRNA-seq tools -- 13 m6A configurations, because ELIGOS2 ran in
#: two modes -- plus the RNA004 built-in Dorado models, and
#: `code/harmonisation/common/config.py` names the callsets that are **not part of the
#: manuscript** and deletes them in `11_scope_split.py`: differr, EpiNano_SVM,
#: Tombo_com, CHEUI-diff and mAFiA.  Rows for those tools had nevertheless been
#: printed in S6 and S7; they are filtered out here, and the asserts after each
#: table fail if one ever comes back.
PAPER_TOOLS = frozenset({
    # 13 m6A configurations (Fig. 1-6, Fig. 8, Fig. S8)
    "CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
    "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
    "NanoSPA_m6A", "xPore", "Yanocomp",
    # non-m6A tools of Fig. 7
    "CHEUI_m5C", "NanoMUD_psi", "NanoMUD_m1psi", "NanoNm", "NanoPsu",
    "NanoSPA_psU",
    # RNA004 built-in modification models (Fig. 8, Fig. S8)
    "Dorado",
})

#: how a row of an inventory table names one of the tools above
SCOPE_ALIAS = {
    "CHEUI": "CHEUI_m6A", "ELIGOS2": "ELIGOS2_diff", "EpiNano": "EpiNano_Error",
    "NanoSPA": "NanoSPA_m6A", "NanoMUD": "NanoMUD_psi", "yanocomp": "Yanocomp",
    "Dorado (RNA004 modification models)": "Dorado",
}


def in_scope(name: object) -> bool:
    """True when a row of an inventory table names a tool of this study."""
    key = str(name).strip()
    return SCOPE_ALIAS.get(key, key) in PAPER_TOOLS


def plain(text: str) -> str:
    """Strip the curated em dashes from a cell, a title or a note.

    Also fixes the casing of the tool name Yanocomp, which a few source tables lower-case.
    """
    for old, new in DASH_FIX:
        text = text.replace(old, new)
    for old, new in CASE_FIX:
        text = text.replace(old, new)
    return text



#: Two cells of the CHEUI_m6A row still described the differential mode
#: (CHEUI-diff), which is declared out of scope by `harmonisation/common/config.py` and
#: whose row the scope filter removes, and one of them pointed the reader to it
#: ("For CHEUI-diff see below.") though nothing follows.  The cells now describe
#: the single-sample mode alone.
MODE_FIX = (
    ("model probability (default); CHEUI-diff upper_cutoff/lower_cutoff",
     "model probability (default)"),
    ("None for CHEUI-solo (model1/model2 emit probability only, no p-value). "
     "For CHEUI-diff see below.",
     "None: model1/model2 emit probability only, no p-value"),
)


def cell(text: object) -> str:
    """A printed cell, **complete**.

    2026-09-27 (user): the per-tool table clipped every long cell with an
    ellipsis (34 of them), which hid exactly the details the table exists for --
    the dependency list, the model checkpoints, the thresholds that were run.
    Cells are now printed in full: whitespace is collapsed, the curated cell
    separator ``|`` that stands for a line break becomes ``/``, the em-dash and
    casing fixes are applied, and nothing is cut.  The reader who wants the raw
    cell plus its provenance opens the deposited machine-readable table.
    """
    s = re.sub(r"\s+", " ", str(text)).strip()
    s = s.replace("|", "/").replace("\t", " ")
    for old, new in MODE_FIX:
        s = s.replace(old, new)
    return plain(s)


def write_table(name: str, title: str, columns: list[str], rows: list[list[object]],
                note: str = "") -> None:
    path = OUT / f"{name}.tsv"
    title, note = plain(title), plain(note)
    with path.open("w") as fh:
        fh.write(f"# {title}\n")
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            line = plain("\t".join(str(c) for c in r))
            assert "—" not in line, f"{name}: em dash in a rendered cell: {line[:90]}"
            assert "yanocomp" not in line, f"{name}: lower-cased tool name: {line[:90]}"
            fh.write(line + "\n")
        if note:
            assert "—" not in note, f"{name}: em dash in the note: {note[:90]}"
            fh.write(f"# note: {note}\n")
    print(f"{name}: {len(rows)} rows -> {path}")


# ---------------------------------------------------------------- S1 tools
TOOL_ROWS = [
    # tool, classification strategy, classification method, algorithms, target, modes evaluated here
    ("m6Anet", "Signal (current)", "Deep learning", "Multiple-instance-learning neural network",
     "m6A", "de novo (RNA002 and RNA004-optimised models)"),
    ("MINES", "Signal (current)", "Machine learning", "Random forest",
     "m6A", "de novo"),
    ("DENA", "Signal (current)", "Deep learning", "Bidirectional LSTM",
     "m6A", "de novo"),
    ("Nanom6A", "Signal (current)", "Machine learning", "XGBoost",
     "m6A", "de novo"),
    ("Yanocomp", "Signal (current)", "Machine learning + statistical test",
     "Gaussian mixture model, G-test", "m6A", "de novo"),
    ("xPore", "Signal (current)", "Clustering",
     "Multi-sample two-Gaussian mixture model", "m6A", "comparative (WT vs control)"),
    ("DRUMMER", "Basecalling errors", "Statistical test",
     "Comparative profiling of basecall error rates, G-test", "m6A",
     "comparative (WT vs control)"),
    ("ELIGOS2", "Basecalling errors", "Statistical test", "Fisher's exact test",
     "m6A, m5C", "differential (ELIGOS2_diff) and single-sample (ELIGOS2_solo)"),
    ("EpiNano", "Basecalling errors", "Machine learning", "Support vector machine (error features)",
     "m6A", "error-based mode (EpiNano_Error)"),
    ("Nanocompore", "Signal (current)", "Machine learning + statistical test",
     "Two-GMM followed by logit test", "m6A", "comparative (WT vs control)"),
    ("CHEUI", "Signal (current)", "Deep learning + statistical test",
     "Convolutional neural network, two-tailed Mann-Whitney U-test", "m6A, m5C",
     "CHEUI_m6A, CHEUI_m5C"),
    ("NanoSPA", "Signal (current)", "Machine learning", "Feedforward neural network",
     "m6A, Psi", "NanoSPA_m6A, NanoSPA_Psi"),
    ("NanoMUD", "Signal (current)", "Machine learning + deep learning",
     "Bidirectional LSTM, regression and motif-specific models", "Psi, m1Psi",
     "NanoMUD_Psi, NanoMUD_m1Psi"),
    ("NanoPsu", "Basecalling errors", "Machine learning + statistical test",
     "U-to-C mismatch error analysis, EXT", "Psi", "de novo"),
    ("NanoNm", "Signal (current)", "Machine learning", "XGBoost", "Nm", "de novo"),
]
write_table(
    "TableS1_tools",
    "Table S1: Modification detection tools benchmarked in this study.",
    ["Tool", "Classification strategy", "Classification method", "Specific algorithms",
     "Targeted modification", "Modes evaluated in this study"],
    TOOL_ROWS,
    "RNA002 datasets for all tools; Dorado built-in modification models (RNA004) and the "
    "RNA004-optimised m6Anet model are reported separately in Fig. 8 / Fig. S8. Tool names "
    "follow the standardized nomenclature of the manuscript; versions and command lines are "
    "in Table S6.",
)

# ------------------------------------------------- S2 HeLa per-replicate counts
summ = pd.read_csv((_RB / "figures/figure1/inputs/tables/tool_counts_group_summary.tsv"), sep="\t")
rows = []
for cond in ("HeLa_WT", "HeLa_IVT"):
    sub = summ[summ.dataset_group == cond].sort_values("sites_mean", ascending=False)
    for r in sub.itertuples():
        each = [int(float(x)) for x in str(r.sites_each).split("|")]
        rows.append([cond.replace("HeLa_", "HeLa "), r.tool, *each,
                     round(r.sites_mean, 0), round(r.sites_sd, 0) if r.sites_sd == r.sites_sd else "n.a.",
                     round(r.sites_sd / r.sites_mean, 3) if r.sites_sd == r.sites_sd else "n.a."])
write_table(
    "TableS2_hela_replicates",
    "Table S2: m6A sites detected in the three independent HeLa sequencing units per tool.",
    ["Condition", "Tool", "Unit 1", "Unit 2", "Unit 3", "Mean", "SD", "CV"],
    rows,
    "One row per tool and condition; each unit is an independent sequencing unit (HeLa_WT1-3, "
    "HeLa_IVT_rep1-3). The values are those plotted in Fig. 1B and are deposited in the "
    "machine-readable summary table that accompanies the source code.",
)

# ------------------------------------------- S3 T/WT on the shared testable space
tw = pd.read_csv((_RB / "figures/figure2/tables/fig2a_testable_ratio.tsv"), sep="\t")
rows = []
for (sp, tool), sub in tw.groupby(["species", "tool"], sort=True):
    for r in sub.sort_values("pair").itertuples():
        rows.append([sp, tool, r.pair, int(r.n_common_universe), int(r.n_wt_common),
                     int(r.n_ctrl_common), round(r.ctrl_wt_ratio, 3)])
    rows.append([sp, tool, "mean", "", "", "", round(sub.ctrl_wt_ratio.mean(), 3)])
write_table(
    "TableS3_shared_universe_ratio",
    "Table S3: Wild-type versus modification-deficient calls on the shared testable universe.",
    ["Species", "Tool", "Pair", "Shared universe (positions)", "WT calls",
     "Deficient-condition calls", "Deficient/WT ratio"],
    rows,
    "The shared testable universe of a matched pair is the set of exonic positions with "
    "coverage >= 10 and a modification-compatible reference base in BOTH samples; the ratio is "
    "the deficient-condition call count over the wild-type call count on that universe "
    "(Fig. 2A). Mouse pairs are one per independent study and are never pooled. Deficient "
    "samples are fip37 KD (Arabidopsis), Mettl3 KO (mouse) and the unmodified IVT library "
    "(HeLa); KO/KD are partial negatives, so the ratio is a relative enrichment change, not a "
    "specificity measure. The per-pair values are deposited in the machine-readable table that "
    "accompanies the source code.",
)

# ------------------------------------------------------------- S4 top-5 5-mers
t5 = pd.read_csv((_RB / "figures/figureS2/analysis/figS2_top5_full.tsv"), sep="\t")
rows = []
for r in t5.sort_values(["Species", "Tool", "Rank"]).itertuples():
    rows.append([r.Species, r.Tool_display, int(r.Rank), r.Kmer,
                 int(r.Count_pooled), round(r.RelFreq_repmix, 4),
                 "yes" if str(r.Is_RRACH).lower() in ("true", "1") else "no"])
write_table(
    "TableS4_top5_motifs",
    "Table S4: Five most frequently detected 5-mer motifs per tool and species.",
    ["Species", "Tool", "Rank", "5-mer", "Count", "Relative frequency", "RRACH"],
    rows,
    "Motifs are counted on the strand-normalized 5-mers centred on each called A of the "
    "cleaned call set, with equal weight per independent unit; the relative frequency is the "
    "replicate-mixed frequency within the tool x species 5-mer pool (Fig. S2B). Mouse uses the "
    "mES_WT study only; E. coli has no DENA or MINES run. The full per-tool table is deposited "
    "with the source code.",
)

# ---------------------------------------------------------- S5 non-m6A detection
# 2026-09-28 (union drop): the per-unit means are the primary columns; the raw
# unions of the three replicates stay as labelled reference columns, and the
# coverage-matched ratio (calls compared inside matched coverage strata) is read
# from the depth-diagnosis evidence of the R3-9b package.
cnt = pd.read_csv((_RB / "figures/figureS7/tables/s7_counts_summary.tsv"), sep="\t")
ci = pd.read_csv((_RB / "figures/figureS7/tables/s7_ratio_ci.tsv"), sep="\t")
mat = pd.read_csv((_RB / "analysis/R3-9b_nonm6a_replicate_depth_20260928/evidence/band_matched_rate_ratios.tsv"),
                  sep="\t")
order = ["CHEUI_m5C", "NanoMUD_psi", "NanoMUD_m1psi", "NanoNm", "NanoSPA_psU", "NanoPsu"]
disp = {"CHEUI_m5C": "CHEUI-m5C", "NanoMUD_psi": "NanoMUD-Psi", "NanoMUD_m1psi": "NanoMUD-m1Psi",
        "NanoNm": "NanoNm", "NanoSPA_psU": "NanoSPA-Psi", "NanoPsu": "NanoPsu"}
rows = []
for tool in order:
    wt = cnt[(cnt.tool == tool) & (cnt.condition == "WT")].iloc[0]
    ivt = cnt[(cnt.tool == tool) & (cnt.condition == "IVT")].iloc[0]
    c = ci[ci.tool == tool].iloc[0]
    m = mat[mat.tool == tool].iloc[0]
    rows.append([disp[tool], ivt.mod_type,
                 f"{wt.count_mean_in_universe:.0f} ({wt.count_sd_in_universe:.0f})",
                 f"{ivt.count_mean_in_universe:.0f} ({ivt.count_sd_in_universe:.0f})",
                 f"{c.ratio_mean_counts:.3f} ({c.ci_lo:.3f}-{c.ci_hi:.3f})",
                 f"{m.band_matched_rate_ratio:.3f} ({m.boot_matched_ci_lo:.3f}-{m.boot_matched_ci_hi:.3f})",
                 int(wt.union_raw), int(ivt.union_raw)])
write_table(
    "TableS5_non_m6a",
    "Table S5: Non-m6A calls in HeLa and the unmodified-IVT-to-WT ratio.",
    ["Tool", "Modification", "WT calls per unit, mean (SD)",
     "Unmodified IVT calls per unit, mean (SD)", "Mean-of-counts ratio (95% CI)",
     "Coverage-matched IVT/WT ratio (95% CI)", "WT union, reference only",
     "Unmodified IVT union, reference only"],
    rows,
    "Counts are means (SD) over the three independent HeLa sequencing units of each "
    "condition, inside the candidate-site set, and are the primary measure. The "
    "mean-of-counts ratio is the unmodified-IVT-to-WT ratio of those unit means with its "
    "unit-level percentile bootstrap 95% CI (B = 1000, seeds 20260920/20260921). The "
    "coverage-matched ratio compares calls from candidate sites within matched coverage "
    "strata, which removes the difference in coverage composition between the two "
    "libraries. The raw unions of the three replicates are listed for reference only and "
    "must not be read as replicate-level quantities. WT = wild type, unmodified IVT = "
    "negative control (not a treatment arm). CHEUI-m5C uses the corrected call set "
    "(coordinate fix and reference-base filter, 2026-09-18): the union of 48,627 WT calls "
    "quoted in the original submission predates that correction. Values for the other "
    "tools are unchanged with respect to the published table. The per-unit counts and "
    "ratio intervals are deposited in the machine-readable tables that accompany the "
    "source code, reconciled row by row against the frozen evidence set behind this table.",
)

# ---------------------------------------------- S6 per-tool implementation table
na4 = pd.read_csv((_RB / "analysis/supp_table_inputs/per_tool_implementation.csv"))
rows = []
dropped_s6: list[str] = []
for r in na4.itertuples():
    if not in_scope(r.tool_display or r.tool):
        dropped_s6.append(str(r.tool_display or r.tool))
        continue
    rows.append([
        cell(r.tool_display or r.tool),
        cell(r.software_or_version),
        
        #: table reports the model/checkpoint and the filtering parameters of every
        #: tool, and the reviewer asked for the exact configuration -- so the two
        #: columns are printed, not left to the deposited table
        cell(r.model_checkpoint),
        cell(r.required_input),
        cell(r.min_read_coverage),
        cell(r.filtering_parameters),
        cell(r.calling_threshold),
        cell(r.multiple_testing_correction),
        cell(r.default_vs_optimised),
        cell(r.coordinate_harmonisation),
    ])
write_table(
    "TableS6_per_tool_implementation",
    "Table S6: Per-tool implementation details (software, model, input, coverage floor, "
    "filtering, thresholds, multiple-testing correction, coordinate harmonisation).",
    ["Tool", "Software / version", "Model / checkpoint", "Required input",
     "Minimum read coverage", "Filtering parameters", "Calling threshold",
     "Multiple-testing correction", "Default or optimised",
     "Coordinate harmonisation"],
    rows,
    "One row per tool and configuration used in this study (13 m6A configurations, the non-m6A "
    "tools of the modification panels, and the RNA004 built-in modification models); callsets "
    "that exist in the working tree but are outside the study are not part of the manuscript and "
    "are not listed. Generated from the artifacts of this benchmark and from the tool "
    "documentation; every cell is printed in full and every curated cell is marked in the "
    "evidence column of the machine-readable table deposited with the source code.",
)

#: gate -- the table carries the study's tools and nothing else
_printed_s6 = [r[0] for r in rows]
assert sorted(dropped_s6) == ["Tombo", "differr", "mAFiA"], dropped_s6

#: a reader who opens the machine-readable copy meets the tools the paper never
#: used.  The filtered rows are written as a CSV beside the printed TSV.
with ((_RB / "tables/tables/TableS6_per_tool_implementation.csv")).open("w", encoding="utf-8",
                                                         newline="") as _fh:
    _w = csv.writer(_fh, quoting=csv.QUOTE_MINIMAL)
    _w.writerow(["Tool", "Software / version", "Model / checkpoint", "Required input",
                 "Minimum read coverage", "Filtering parameters", "Calling threshold",
                 "Multiple-testing correction", "Default or optimised",
                 "Coordinate harmonisation"])
    _w.writerows(rows)
for _name in ("CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
              "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
              "NanoSPA", "xPore", "Yanocomp", "NanoMUD", "NanoNm", "NanoPsu",
              "Dorado (RNA004 modification models)"):
    assert _name in _printed_s6, f"S6 lost {_name}"
assert len(_printed_s6) == 17, _printed_s6
for _bad in ("Tombo", "differr", "mAFiA", "CHEUI-diff", "Epinano_SVM"):
    assert _bad not in _printed_s6, f"S6 still lists {_bad}"

#: and no cell may point the reader below to a row that the scope filter removed
_blob_s6 = "\n".join(c for r in rows for c in r)
for _gone in ("CHEUI-diff", "CHEUI-solo", "see below"):
    assert _gone not in _blob_s6, f"S6 cell still says {_gone!r}"

# --------------------------------------------- S7 chemistry / inclusion matrix
# 2026-09-26: the "no / no" rows used to read "not run - BAM-based, applicable", which
# contradicted the two chemistry columns.  Each row now carries one of four explicit
# failure or exclusion classes: (i) raw signal that cannot be produced from RNA004 data,
# (ii) an RNA002-trained model with no RNA004 checkpoint evaluated here, (iii) a software
# component that is not maintained for the new chemistry, (iv) a configuration that was
# run but whose output is not used in this study.
SIGNAL_MODEL_CLASS = ("input format: raw signal (a nanopolish/f5c eventalign or Tombo "
                      "resquiggle) cannot be produced from RNA004 data")
TRAINED_MODEL_CLASS = ("basecalling model: released model trained on RNA002 features; "
                       "no RNA004 checkpoint evaluated in this study")
MAINTENANCE_CLASS = ("software maintenance: the resquiggle component is not maintained "
                     "for the new chemistry")
OUT_OF_SCOPE_CLASS = "outside the scope of this study: output generated but not used"


def _s7_class(status: str, reason: str, used_002: bool) -> tuple[str, str]:
    """(status, failure/exclusion class) for one compatibility-matrix row.

    2026-09-28 (user): the ``used_002`` fallback used to answer "basecalling
    model" whenever none of the patterns matched, so MINES and Yanocomp -- whose
    recorded reason is "Tombo resquiggle not available for RNA004" -- came out as
    model-limited while their own Input class column says raw signal, and while
    the Results section lists them among the six tools that need a resquiggle.
    The resquiggle reason now maps to the input-format class, and the fallback
    stays only for a genuinely unmatched reason.
    """
    s, r = str(status).lower(), str(reason).lower()
    if s == "included":
        return "included in the RNA004 benchmark", ""
    if "trained on rna002" in r:
        return "not applied to RNA004", TRAINED_MODEL_CLASS
    if "resquiggle unsupported" in r or "not maintained" in r:
        return "not applied to RNA004", MAINTENANCE_CLASS
    if "resquiggle" in r or "signal model" in r or "same as above" in r:
        return "not applied to RNA004", SIGNAL_MODEL_CLASS
    if used_002:
        # benchmarked on RNA002; the BAM-based caller itself would run, but its
        # released model is RNA002-trained and no RNA004 checkpoint was evaluated
        return "not applied to RNA004", TRAINED_MODEL_CLASS
    return "not included in this benchmark", OUT_OF_SCOPE_CLASS


#: The compatibility matrix establishes "run on a chemistry" from the call-directory name;
#: EpiNano's RNA002 export lives in `$RNAMODBENCH_LOCAL/raw/result/EpiNano_DiffErr`, so the automatic
#: check missed it (2026-09-26).  EpiNano is one of the 15 benchmarked tools (Table S1).
RNA002_OVERRIDE = {"EpiNano_Error": True}

na5 = pd.read_csv((_RB / "analysis/supp_table_inputs/chemistry_compatibility.csv"))
rows = []
dropped_s7: list[str] = []
for r in na5.itertuples():
    if not in_scope(r.tool_display or r.tool):
        dropped_s7.append(str(r.tool_display or r.tool))
        continue
    ran_002 = (str(r.ran_on_RNA002).lower() == "true"
               or RNA002_OVERRIDE.get(str(r.tool_display), False))
    used_002 = "yes" if ran_002 else "no"
    used_004 = "yes" if str(r.ran_on_RNA004).lower() == "true" else "no"
    status, cls = _s7_class(r.status_RNA004, r.reason_if_excluded, ran_002)
    rows.append([cell(r.tool_display or r.tool), cell(r.input_class),
                 used_002, used_004, status, cell(cls)])

#: study's tool be applied to the new chemistry", so a tool the study never used
#: has no row here; the response letter promises the matrix for the 15 tools.
write_table(
    "TableS7_chemistry_matrix",
    "Table S7: Applicability of each tool to the RNA002 and RNA004 chemistries, with the "
    "failure or exclusion class where a tool was not applied or not included.",
    ["Tool", "Input class", "Used on RNA002", "Used on RNA004", "Status",
     "Failure or exclusion class"],
    rows,
    "One row per tool and configuration used in this study. Presence on a chemistry is read from "
    "the file system (call directories per tool and sample), not from memory. Input classes: "
    "raw-signal (FAST5 plus a nanopolish/Tombo resquiggle) and basecalled (BAM/FASTQ) callers; "
    "RNA004 native output is POD5 and Dorado BAM with 9-mer models, so RNA002 signal models cannot "
    "be produced for it. Failure and exclusion classes: input format (raw signal not producible "
    "from RNA004 data), basecalling model (RNA002-trained model, no RNA004 checkpoint evaluated), "
    "software maintenance (resquiggle not maintained for the new chemistry), and outside the scope "
    "of this study (output generated but not used).",
)

#: gate -- the table carries the study's tools and nothing else
_printed_s7 = [r[0] for r in rows]
assert sorted(dropped_s7) == ["CHEUI-diff", "Epinano_SVM", "Tombo",
                             "Tombo (comparative)", "differr", "mAFiA"], dropped_s7
#: the deposited input of this table, same rows as the printed one (see S6)
with ((_RB / "tables/tables/TableS7_chemistry_matrix.csv")).open("w", encoding="utf-8",
                                                 newline="") as _fh:
    _w = csv.writer(_fh, quoting=csv.QUOTE_MINIMAL)
    _w.writerow(["Tool", "Input class", "Used on RNA002", "Used on RNA004",
                 "Status", "Failure or exclusion class"])
    _w.writerows(rows)
assert len(_printed_s7) == 16, _printed_s7
for _bad in ("Tombo", "differr", "mAFiA", "CHEUI-diff", "Epinano_SVM",
             "Tombo (comparative)"):
    assert _bad not in _printed_s7, f"S7 still lists {_bad}"

#: story -- a raw-signal tool cannot be excluded for a basecalling model, and a
#: tool that was not applied to RNA004 must carry a class at all.
for _tool, _inp, _u002, _u004, _status, _cls in rows:
    if _inp.startswith("raw-signal"):
        assert _cls != TRAINED_MODEL_CLASS, f"{_tool}: raw-signal tool with a model class"
    else:
        assert _cls != SIGNAL_MODEL_CLASS, f"{_tool}: basecalled tool with an input-format class"
    if _u004 == "no":
        assert _cls, f"{_tool}: excluded from RNA004 without a failure/exclusion class"

# ------------------------------------------- S8 GLORI-to-nanopore matching table
write_table(
    "TableS8_glori_matching",
    "Table S8: Biological source of each GLORI reference relative to the matched nanopore "
    "dataset.",
    ["Species", "Nanopore libraries", "GLORI reference", "Biological relationship",
     "GLORI sites"],
    [
        ["Human", "HeLa WT (SRP393373, RNA002) and HeLa WT/IVT (ERP164259, RNA004)",
         "HeLa, GSE210563 (Liu et al. 2023)",
         "Same cell line, independent laboratory and independent RNA extract",
         "112,451"],
        ["Arabidopsis", "Col-0 WT (SRP329449) and fip37 KD (SRP363914), RNA002",
         "Col-0 seedlings, GSE246632 (Xie et al. 2025)",
         "Same species and accession from independent studies and independent RNA extracts",
         "80,624"],
        ["Mouse", "mESC Mettl3 KO (SRP166020) and mESC WT (SRP357195), RNA002",
         "mESC, GSE246632",
         "Same cell type (mESC) from independent studies and independent RNA extracts",
         "41,961"],
    ],
    "Counts are the replicate-overlapping GLORI sites with modification ratio > 0.1 in both "
     "replicates (Fig. S3A). GSE210563 also contains HEK293T GLORI, and its mouse samples are "
     "MEF; neither was used. Nanopore and GLORI libraries were never derived from the same RNA "
     "extract, so a GLORI-positive site is not assumed to be modified in the corresponding "
     "nanopore sample and the PPV against GLORI is a precision-like measure conditioned on this "
     "site-level partial reference, not sensitivity.",
)

# ------------------------------------------- S9 non-m6A evaluation boundary table
write_table(
    "TableS9_evaluation_boundary",
    "Table S9: Positive and negative references available for each modification type, and the "
    "metrics that can therefore be reported.",
    ["Modification", "Positive reference", "Negative control", "Metrics reported",
     "Limitation stated in the manuscript"],
    [
        ["m6A", "GLORI (site-level partial reference; Table S8)",
         "Unmodified Curlcake (SRP174366, ERP162788) and HeLa IVT",
         "PPV against GLORI, recall, F1, MCC, false-positive burden, replicate consistency",
         "GLORI is an independently generated, DRACH-centred and coverage-limited reference; "
         "absolute values are reference-bounded."],
        ["Psi", "None (no orthogonal Psi map is used for scoring)",
         "Unmodified Curlcake and HeLa IVT",
         "False-positive burden, replicate consistency, enrichment over a chromosome-stratified "
         "permutation null against RMBase + DirectRMDB",
         "Sensitivity and recall are not reported; the m6A-centred GLORI map carries no Psi "
         "labels."],
        ["m1Psi", "None", "Unmodified Curlcake and HeLa IVT",
         "False-positive burden, replicate consistency",
         "No external reference covers m1Psi (RMBase + DirectRMDB covers Psi, m5C and Nm only)."],
        ["m5C", "None (orthogonal UBS-seq m5C used only as an external overlap set)",
         "Unmodified Curlcake and HeLa IVT",
         "False-positive burden, replicate consistency, enrichment over the permutation null",
         "The synthetic constructs are unmodified or fully m6A-substituted, so no m5C-positive "
         "control is available."],
        ["Nm", "None (orthogonal Nm-Mut-seq used only as an external overlap set)",
         "Unmodified Curlcake and HeLa IVT",
         "False-positive burden, replicate consistency, enrichment over the permutation null",
         "Same as m5C: no synthetic Nm-positive control, and GLORI provides no Nm labels."],
    ],
    "The Curlcake libraries used here are the unmodified IVT set and the fully m6A-substituted "
     "set (SRP174366 / GSE124309) plus the unmodified RNA004 set (ERP162788); later Curlcake "
     "collections that spike in Psi, m5C or m1Psi were not among the deposited accessions and "
     "were not used.",
)

index = [
    
    #: neither the reviewer-point tags of the revision (E4 / R1-5, R3-10 / R1-6,
    #: ...) nor the pre-deposit working paths (superseded_output/, fig*_revision/).
    #: It also stopped at S9 while the delivered SI has twelve tables, so every
    #: row now names the released file a reader can open.  The gate below fails
    #: if a tag or a working path ever comes back.
    ["S1", "Modification detection tools benchmarked in this study",
     "tables/tables/TableS1_tools.tsv"],
    ["S2", "HeLa m6A calls per independent sequencing unit",
     "tables/tables/TableS2_hela_replicates.tsv"],
    ["S3", "Wild-type vs deficient calls on the shared testable universe",
     "tables/tables/TableS3_shared_universe_ratio.tsv"],
    ["S4", "Five most frequent 5-mer motifs per tool and species",
     "tables/tables/TableS4_top5_motifs.tsv"],
    ["S5", "Non-m6A calls and IVT/WT ratios",
     "tables/tables/TableS5_non_m6a.tsv"],
    ["S6", "Per-tool implementation details",
     "tables/tables/TableS6_per_tool_implementation.tsv; "
     "analysis/supp_table_inputs/per_tool_implementation.csv"],
    ["S7", "Chemistry applicability matrix",
     "tables/tables/TableS7_chemistry_matrix.tsv; "
     "analysis/supp_table_inputs/chemistry_compatibility.csv"],
    ["S8", "GLORI reference vs matched nanopore dataset",
     "tables/tables/TableS8_glori_matching.tsv"],
    ["S9", "Evaluation boundary per modification type",
     "tables/tables/TableS9_evaluation_boundary.tsv"],
    ["S10", "Datasets and replicate structure",
     "tables/tables/TableS10_datasets.tsv"],
    ["S11", "Tool nomenclature used throughout the manuscript",
     "tables/tables/TableS11_nomenclature.tsv"],
    ["S12", "Sensitivity of the site-level ranking to coverage and to the "
     "reference composition",
     "tables/tables/TableS12_coverage_rank_stability.tsv; "
     "analysis/coverage_rank_stability/"],
]
_index_blob = "\n".join("\t".join(r) for r in index)
for _tag in ("R1-", "R2-", "R3-", "E1 /", "E4 /", "superseded_output", "figures/figure2",
             "figures/figureS2", "figures/figureS7", "figure1b_replicates", "curated"):
    assert _tag not in _index_blob, f"internal trace in the table index: {_tag}"
assert len(index) == 12, len(index)
with ((_RB / "tables/tables/_index.tsv")).open("w") as fh:
    fh.write("# Supplementary table index\n")
    fh.write("table\tcontent\tsource\n")
    for r in index:
        fh.write("\t".join(r) + "\n")
print("_index.tsv written")
