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


#: 2026-09-26 (author decision, de-AI pass): two curated cells of the deposited evidence CSV carry an
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



#: manuscript -- underscores rather than hyphen variants, and ELIGOS2 rather than
#: ELIGOS -- and the manuscript writes pseudouridine as the Greek Psi, so the
#: printed tables do too (the tables are rendered with an embedded Arial that
#: carries the glyph; the earlier Helvetica could not).  The same pass folds the
#: modification labels, which are what is left once the tool names have been
#: rewritten: "m1Psi" first, then the bare "Psi".  The upstream evidence tables
#: are deposited byte-identical and keep whatever they carry.
NAME_FIX = (
    ("ELIGOS_diff", "ELIGOS2_diff"),
    ("ELIGOS_solo", "ELIGOS2_solo"),
    ("CHEUI-m5C", "CHEUI_m5C"),
    ("NanoMUD-Psi", "NanoMUD_Psi"),
    ("NanoMUD-m1Psi", "NanoMUD_m1Psi"),
    ("NanoSPA-Psi", "NanoSPA_Psi"),
    ("m1Psi", "m1\u03a8"),
    ("Psi", "\u03a8"),
)



#: not everything that happens to sit in the working tree.  The manuscript
#: benchmarks 15 dRNA-seq tools -- 13 m6A configurations, because ELIGOS2 ran in
#: two modes -- plus the RNA004 built-in Dorado models, and
#: `code/harmonisation/common/config.py` names the callsets that are **not part of the
#: manuscript** and deletes them in `11_scope_split.py`: differr, EpiNano_SVM,
#: Tombo_com, CHEUI-diff and mAFiA.  Rows for those tools had nevertheless been
#: printed in S6 and S7; they are filtered out here, and the asserts after each
#: table fail if one ever comes back.
PAPER_TOOLS = frozenset({
    # 13 m6A configurations (Fig. 1-6, Fig. 8, Fig. S9)
    "CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
    "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
    "NanoSPA_m6A", "xPore", "Yanocomp",
    # non-m6A tools of Fig. 7
    "CHEUI_m5C", "NanoMUD_psi", "NanoMUD_m1psi", "NanoNm", "NanoPsu",
    "NanoSPA_psU",
    # RNA004 built-in modification models (Fig. 8, Fig. S9)
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
    for old, new in NAME_FIX:
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

    2026-09-27 (author decision): the per-tool table clipped every long cell with an
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
    
    #: itself) and Fig. S10 (the quantitative accuracy panels); Fig. S9 is the
    #: non-coding metagene of the non-m6A tools, and the built-in models for the
    #: remaining modifications are in Fig. S5, so the old pointer to Fig. S9 was
    #: wrong on both counts.
    "RNA002 datasets for all tools; Dorado built-in modification models (RNA004) and the "
    "RNA004-optimised m6Anet model are reported separately in Fig. 8 and Fig. S10, with the "
    "built-in models for the other modifications in Fig. S5. Tool names "
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
#: 2026-09-29 (R1-4 round): the coverage-matched recomputation of the same ratio
#: (figures/figure2/src/01b_fig2_matched_ratio.py) joins Table S3 as a second
#: ratio column; the per-pair values come from the pair table, the summary and
#: its bootstrap interval from the species x tool table.
mt = pd.read_csv((_RB / "figures/figure2/tables/fig2a_matched_ratio.tsv"), sep="\t")
ms = pd.read_csv((_RB / "figures/figure2/tables/fig2a_matched_summary.tsv"), sep="\t")
mk = mt.set_index(["species", "tool", "pair"])["matched_ratio"]
sk = ms.set_index(["species", "tool"])
rows = []
for (sp, tool), sub in tw.groupby(["species", "tool"], sort=True):
    for r in sub.sort_values("pair").itertuples():
        m = float(mk.loc[(sp, tool, str(r.pair))]) if (sp, tool, str(r.pair)) in mk.index \
            else float("nan")
        rows.append([sp, tool, r.pair, int(r.n_common_universe), int(r.n_wt_common),
                     int(r.n_ctrl_common), round(r.ctrl_wt_ratio, 3), round(m, 3)])
    s = sk.loc[(sp, tool)]
    ci = "n.a." if s.matched_ci_lo != s.matched_ci_lo \
        else f"{s.matched_ci_lo:.3f}-{s.matched_ci_hi:.3f}"
    rows.append([sp, tool, "mean", "", "", "", round(sub.ctrl_wt_ratio.mean(), 3),
                 f"{s.matched_mean:.3f} ({ci})"])
write_table(
    "TableS3_shared_universe_ratio",
    "Table S3: Wild-type versus modification-deficient calls on the shared testable universe.",
    ["Species", "Tool", "Pair", "Shared universe (positions)", "WT calls",
     "Deficient-condition calls", "Deficient/WT ratio",
     "Coverage-matched ratio (95% CI)"],
    rows,
    "The shared testable universe of a pair is the set of exonic positions with coverage >= "
"10 and a modification-compatible reference base in both samples; the ratio is the "
"deficient-condition count over the wild-type count there (Fig. 2A). Deficient samples "
"are fip37 KD (Arabidopsis), Mettl3 KO (mouse) and the unmodified HeLa IVT library; KO/KD "
"are partial negatives, so the ratio is an enrichment change, not a specificity measure. "
"The coverage-matched ratio column recomputes the counts inside matched coverage strata "
"(10/20/30/50/100/200/400 reads); mouse has one pair per study and no interval.",
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
# 2026-09-28 (union drop) / 2026-10-06 (author decision, second task): the per-unit counts and the
# raw unions have left this table altogether; its two count columns are the
# depth-standardised rate of each condition and the coverage-matched ratio (calls
# compared inside matched coverage strata) is their quotient, read from the
# depth-diagnosis evidence of the R3-9b package.

# value exceeds one for four of the six tools (NanoNm 1.71, NanoMUD-m1Psi 1.16,
# NanoMUD-Psi 1.16, CHEUI-m5C 1.08), which reads as the unmodified control calling
# more sites than the wild type, while the R3-9b depth diagnosis shows the excess is
# the coverage composition of the two libraries and not the modification state: inside
# matched coverage strata all six tools fall to 0.71-0.96.  The matched ratio is now
# the leading quantitative column and the mean-of-counts ratio is no longer printed;
# `band_matched_rate_ratios.tsv` keeps `count_ratio_ivt_over_wt` for the record, and
# `s7_ratio_ci.tsv` keeps the retired interval for the deposited evidence.
cnt = pd.read_csv((_RB / "figures/figureS7/tables/s7_counts_summary.tsv"), sep="\t")
mat = pd.read_csv((_RB / "analysis/nonm6a_replicate_depth/evidence/band_matched_rate_ratios.tsv"),
                  sep="\t")
order = ["CHEUI_m5C", "NanoMUD_psi", "NanoMUD_m1psi", "NanoNm", "NanoSPA_psU", "NanoPsu"]
disp = {"CHEUI_m5C": "CHEUI-m5C", "NanoMUD_psi": "NanoMUD-Psi", "NanoMUD_m1psi": "NanoMUD-m1Psi",
        "NanoNm": "NanoNm", "NanoSPA_psU": "NanoSPA-Psi", "NanoPsu": "NanoPsu"}

#: stoichiometry published with the model, and as of the second pass that set is *merged
#: into* the wild-type arm rather than standing in for it: the arm carries three samples --
#: this study's HeLa_WT2 and HeLa_WT3 (HeLa_WT1, the smallest call count of the four
#: candidates, is left out so the arm holds three samples as every other row does) plus the
#: published pooled run -- so its mean (SD), its union and its coverage-matched ratio all
#: keep the estimator the other five rows use.  Every
#: number is read from what `76_nanonm_published_wt.py` writes and the anchors of that run
#: are cross-checked here, so the table cannot print a value the evidence does not hold.
#: The unmodified-IVT columns stay this study's own: the published set has no unmodified
#: control and no replicates.
NANONM_PUB = (_RB / "figures/figureS7/tables/s7_nanonm_published_wt.tsv")
NANONM_MATCHED = (_RB / "figures/figureS7/tables/s7_nanonm_published_matched.tsv")
NANONM_ANCHORS = (_RB / "figures/figureS7/logs/76_anchors.tsv")
for _need in (NANONM_PUB, NANONM_MATCHED, NANONM_ANCHORS):
    if not _need.exists():
        raise SystemExit(f"missing {_need}; run src/harmonisation/scripts/"
                         f"76_nanonm_published_wt.py before building the tables")
npub = pd.read_csv(NANONM_PUB, sep="\t")
_per_unit = pd.read_csv((_RB / "figures/figureS7/tables/s7_counts_per_replicate.tsv"), sep="\t")
_anch = pd.read_csv(NANONM_ANCHORS, sep="\t").set_index("key")["value"]


def _anchor(key: str) -> float:
    """One number of the 76 run, so this table prints what that run computed."""
    if key not in _anch.index:
        raise SystemExit(f"the 76 anchors carry no {key!r}; rebuild before the table")
    return float(_anch[key])


#: the arm Table S5 prints: the three samples the wild-type cells are computed over
NANONM_ARM = npub[(npub.source == "merged wild-type arm")
                  & (npub.cut == "matched_ratio_support")]
_arm_union = NANONM_ARM[NANONM_ARM.unit == "union_of_the_three_samples"]
if len(_arm_union) != 1:
    raise SystemExit(f"expected one merged-arm union row, found {len(_arm_union)}")
NANONM_WT_UNION = int(_arm_union.iloc[0].sites_in_universe)
_arm_samples = NANONM_ARM[NANONM_ARM.unit != "union_of_the_three_samples"]
NANONM_WT_MEAN = float(_arm_samples.sites_in_universe.mean())
NANONM_WT_SD = float(_arm_samples.sites_in_universe.std(ddof=1))
if (round(NANONM_WT_MEAN, 1) != _anchor("merged_wt_mean_per_unit")
        or round(NANONM_WT_SD, 1) != _anchor("merged_wt_sd_per_unit")):
    raise SystemExit("the merged-arm samples of the evidence table do not give the mean "
                     "(SD) the 76 run anchored; the arm printed here is not the one cut")
#: the ratio of that arm, from the per-stratum rates the evidence table carries
_matched = pd.read_csv(NANONM_MATCHED, sep="\t")
NANONM_MATCHED_ARM = float(
    (_matched.pooled_weight * _matched.ivt_this_study_rate_per_1e6).sum()
    / (_matched.pooled_weight * _matched.wt_merged_arm_rate_per_1e6).sum())
if abs(NANONM_MATCHED_ARM - _anchor("merged_ratio_matched")) > 5e-4:
    raise SystemExit(f"the band table gives {NANONM_MATCHED_ARM:.4f}, the 76 run "
                     f"anchored {_anchor('merged_ratio_matched'):.4f}")
NANONM_CI_LO = _anchor("merged_ratio_matched_ci_lo")
NANONM_CI_HI = _anchor("merged_ratio_matched_ci_hi")
#: 2026-10-06 (author decision, option C): those two anchors now hold the *counting* interval.  If an
#: older evidence file is on disk they would hold the unit-level one, and the table would
#: silently print a second convention beside the five rows, so the convention named in the
#: anchors is checked first.
if not str(_anch["merged_ratio_matched_ci_convention"]).startswith("site-level"):
    raise SystemExit(f"the 76 anchors were written by an older run: the printed interval's "
                     f"convention is {_anch['merged_ratio_matched_ci_convention']!r}")
#: the single-source arms stay on the record: the note names the published one, and the
#: within-study one is the value the other five rows and Figure 7A report
_pub_sel = npub[(npub.source.str.startswith("published"))
                & (npub.cut == "matched_ratio_support")
                & (npub.unit == "union_of_the_three_universes")]
if len(_pub_sel) != 1:
    raise SystemExit(f"expected one published NanoNm union row, found {len(_pub_sel)}")
NANONM_PUB_CELL = int(_pub_sel.iloc[0].sites_in_universe)
NANONM_PUB_RATE = float(_pub_sel.iloc[0].rate_per_1e6)
NANONM_PUB_MATCHED = float(
    (_matched.pooled_weight * _matched.ivt_this_study_rate_per_1e6).sum()
    / (_matched.pooled_weight * _matched.wt_published_sample_rate_per_1e6).sum())
_own_row = cnt[(cnt.tool == "NanoNm") & (cnt.condition == "WT")].iloc[0]
NANONM_OWN_UNION = int(_own_row.union_in_universe)
_own_union_row = npub[npub.unit == "union_of_the_three_units"].iloc[0]
NANONM_OWN_RATE = float(_own_union_row.rate_per_1e6)
NANONM_OWN_UNITS = _per_unit[(  # per-unit rates, from the replicate-level table
    _per_unit.tool == "NanoNm") & (_per_unit.condition == "WT")]
NANONM_OWN_RATE_MIN = float((NANONM_OWN_UNITS.n_calls_in_universe
                             / NANONM_OWN_UNITS.universe_n * 1e6).min())
NANONM_OWN_RATE_MAX = float((NANONM_OWN_UNITS.n_calls_in_universe
                             / NANONM_OWN_UNITS.universe_n * 1e6).max())
#: the within-study arm and its interval, so the note can name what the merged value is
#: made of and what the other five rows are computed with
_nano_own = mat[mat.tool == "NanoNm"].iloc[0]
_nano_own_ratio = float(_nano_own.band_matched_rate_ratio)
#: 2026-10-06 (author decision, option C): every interval this table prints is the counting one, so
#: the within-study value the note names carries the counting interval of that arm too
_nano_own_lo = _anchor("matched_ci_this_study_site_lo")
_nano_own_hi = _anchor("matched_ci_this_study_site_hi")

#: *depth-standardised* call count of each condition -- the seven coverage strata combined
#: with the pooled site weights of the six HeLa units and the result read on the row's
#: wild-type candidate-site count (see base_sites) -- as a mean (SD) over the condition's
#: samples, and the ratio is the quotient of its two columns.  The raw per-unit counts
#: and the two union columns are gone: they showed the unmodified side *above* the wild type
#: for four of the six tools (the unmodified libraries are deeper), which a reader takes for
#: a contradiction of the ratio.  The raw numbers stay in the manuscript's prose and in the
#: deposited evidence; the (SD) of the standardised count is the between-library spread the
#: counting interval deliberately leaves out, and it is what Figure 7A draws per unit.
RATE_UNITS_EV = ((_RB / "analysis/nonm6a_replicate_depth/evidence/band_matched_unit_rates.tsv"))
if not RATE_UNITS_EV.exists():
    raise SystemExit(f"missing {RATE_UNITS_EV}; run the R3-9b depth diagnosis before the "
                     f"tables")
_rate_units = pd.read_csv(RATE_UNITS_EV, sep="\t")


def rate_cells(tool: str) -> tuple[float, float, float, float]:
    """(WT mean, WT SD, IVT mean, IVT SD) of one tool's per-unit standardised rate."""
    g = _rate_units[_rate_units.tool == tool]
    wt = g[g.condition == "WT"].band_matched_rate_per_1e6.to_numpy(dtype=float)
    ivt = g[g.condition == "IVT"].band_matched_rate_per_1e6.to_numpy(dtype=float)
    if wt.size != 3 or ivt.size != 3:
        raise SystemExit(f"{tool}: expected three wild-type and three IVT units in "
                         f"{RATE_UNITS_EV.name}, found {wt.size} and {ivt.size}")
    return (float(wt.mean()), float(wt.std(ddof=1)),
            float(ivt.mean()), float(ivt.std(ddof=1)))




#: A rate is per candidate site, so one candidate-site count turns it into counts of the
#: estimator the ratio is built from.  The row's own wild-type units supply that count, and
#: the *same* count serves both columns, which is what keeps the printed ratio the quotient
#: of the two printed cells.  NanoNm is the one row whose wild-type arm is not simply the
#: three units (HeLa_WT1 is left out of it, as the note says), so its base is the mean of
#: the two units the arm keeps.
def base_sites(tool: str) -> tuple[float, tuple[str, ...]]:
    """(mean candidate sites per wild-type unit of that row, the units it averages)."""
    g = _per_unit[(_per_unit.tool == tool) & (_per_unit.condition == "WT")]
    if tool == "NanoNm":
        g = g[g["sample"] != "HeLa_WT1"]
    units = tuple(sorted(g["sample"].astype(str)))
    if not units:
        raise SystemExit(f"{tool}: no wild-type unit in the per-replicate counts table")
    if g.universe_n.isna().any() or (g.universe_n <= 0).any():
        raise SystemExit(f"{tool}: a wild-type unit has no candidate-site count")
    return float(g.universe_n.mean()), units


rows = []
for tool in order:
    ivt = cnt[(cnt.tool == tool) & (cnt.condition == "IVT")].iloc[0]
    m = mat[mat.tool == tool].iloc[0]
    wt_rate, wt_rate_sd, ivt_rate, ivt_rate_sd = rate_cells(tool)
    ratio = float(m.band_matched_rate_ratio)
    tol = 1e-6
    if tool == "NanoNm":
        #: the wild-type side is the arm of `76_nanonm_published_wt.py` (two of this study's
        #: units plus the published sample), so its two cells come from that run's anchors --
        #: the row cannot print a rate the arm's own run did not compute.  The unmodified
        #: side is unchanged: this study's three IVT units, whose mean the anchors hold too.
        if abs(ivt_rate - _anchor("ivt_rate_mean_per_1e6")) > 0.05:
            raise SystemExit(f"the per-unit rates give an unmodified-IVT mean of "
                             f"{ivt_rate:.1f}, the 76 run anchored "
                             f"{_anchor('ivt_rate_mean_per_1e6'):.1f}")
        wt_rate = _anchor("merged_wt_rate_mean_per_1e6")
        wt_rate_sd = _anchor("merged_wt_rate_sd_per_1e6")
        ratio = NANONM_MATCHED_ARM
        tol = 5e-4
    #: the two printed cells: that row's standardised rate on the row's own candidate-site
    #: count (see base_sites).  Both columns carry the same base on purpose -- the printed
    #: ratio is the quotient of the rates, so a base the two do not share
    #: would make the three cells disagree in the last digit.
    base, base_units = base_sites(tool)
    wt_sites, wt_sites_sd = wt_rate / 1e6 * base, wt_rate_sd / 1e6 * base
    ivt_sites, ivt_sites_sd = ivt_rate / 1e6 * base, ivt_rate_sd / 1e6 * base
    if not wt_sites > ivt_sites:
        raise SystemExit(
            f"{tool}: the wild-type cell {wt_sites:.0f} sites is not above the "
            f"unmodified-IVT cell {ivt_sites:.0f}; the two count columns would read the "
            f"other way round, which is the reading this table was changed to remove")
    if abs(round(ivt_sites) / round(wt_sites) - ratio) > 5e-3:
        raise SystemExit(
            f"{tool}: the printed cells give {round(ivt_sites)}/{round(wt_sites)} = "
            f"{round(ivt_sites) / round(wt_sites):.4f} against the printed ratio "
            f"{ratio:.4f}; rounding to whole sites moved the two apart")
    if not wt_rate > ivt_rate:
        raise SystemExit(
            f"{tool}: the printed wild-type rate {wt_rate:.1f} per million is not above the "
            f"unmodified-IVT rate {ivt_rate:.1f}; the two rate columns would read the other "
            f"way round, which is the reading this table was changed to remove")
    if abs(ivt_rate / wt_rate - ratio) > tol:
        raise SystemExit(
            f"{tool}: the rate columns give a quotient of {ivt_rate / wt_rate:.6f}, the "
            f"printed ratio is {ratio:.6f}; they are not the same estimator")
    lo, hi = (NANONM_CI_LO, NANONM_CI_HI) if tool == "NanoNm" else \
        (float(m.site_boot_ci_lo), float(m.site_boot_ci_hi))
    if not lo <= ratio <= hi:
        raise SystemExit(
            f"{tool}: the printed ratio {ratio:.4f} falls outside the counting interval "
            f"{lo:.4f}-{hi:.4f}; the evidence file is stale or the two are different "
            f"estimators")
    matched_cell = f"{ratio:.3f} ({lo:.3f}-{hi:.3f})"
    print(f"  S5 {disp[tool]:14s} base {base:>12,.0f} sites over {'+'.join(base_units)}; "
          f"cells {wt_sites:,.0f} / {ivt_sites:,.0f}")
    
    #: divides; the call headers drop to "calls per unit" and keep the Depth-standardised
    #: prefix, and the ratio's own header is unchanged.
    rows.append([disp[tool], ivt.mod_type,
                 matched_cell,
                 f"{wt_sites:.0f} ({wt_sites_sd:.0f})",
                 f"{ivt_sites:.0f} ({ivt_sites_sd:.0f})"])
write_table(
    "TableS5_non_m6a",
    "Table S5: Depth-standardised non-m6A call sites in HeLa and the unmodified-IVT-to-WT "
    "ratio.",
    ["Tool", "Modification",
     "Coverage-matched IVT/WT ratio (95% CI)",
     "Depth-standardised WT calls per unit, mean (SD)",
     "Depth-standardised unmodified-IVT calls per unit, mean (SD)"],
    rows,
    
#: bootstrap (B = 1000, one seed per tool, recorded in `s7_ratio_ci.tsv`) are no longer
#: printed.  The interval above is `rd2_matched_bands.py::site_boot_ratio` (BOOT = 20_000,
#: SEED = 20260928), the *counting* interval of the two arms; the same estimator read with
#: whole units resampled is `boot_matched_ci_*` in the same evidence file and is what the
#: R3-9b figure draws, so the two must not be mixed up when the text is edited.  The "must
#: not be read as reproducibility" construction stays retired (2026-09-28).
#: 2026-10-06 (author decision, second task): the raw per-unit counts and the two union columns left the
#: table -- they put the unmodified side above the wild type for four of six tools -- and the
#: two call columns are now the standardised rate the ratio is built from, so the printed
#: ratio is the quotient of the two printed cells.  The note says so, and says why
#: the raw counts read the other way; the raw numbers stay in the manuscript's prose, in the
#: Figure S8 legend and in the deposited evidence.
f"Calls are the sites a unit calls: the seven coverage strata combined by the pooled site "
f"weights, read on the tool's wild-type candidate-site count per unit; the ratio is their "
f"quotient. Raw counts would mislead: the unmodified libraries hold 1.8-fold the share of "
f"sites above coverage 100, putting the raw control above the wild type for four tools; "
f"standardisation removes it. WT = wild type, unmodified IVT = negative control. CHEUI_m5C "
f"uses the corrected call set. Intervals are counting (Poisson) intervals: they carry "
f"counting error inside these libraries, not the spread between them, which the (SD) "
f"reports. NanoNm's arm is our HeLa_WT2 and HeLa_WT3 units plus the model's published "
f"per-site Nm set, which has no control; alone, our units give "
f"{_nano_own_ratio:.3f} ({_nano_own_lo:.3f}-{_nano_own_hi:.3f}), the published set "
f"{NANONM_PUB_MATCHED:.3f}.",
)

# ---------------------------------------------- S6 per-tool implementation table

#: not take, which contradicted both the two chemistry columns of Table S7 and the
#: Experimental Section of the manuscript ("the BAM-based callers").  Evidence: the
#: tools' own argument parsers and READMEs (`$RNAMODBENCH_LOCAL/source_code/benchmark/DRUMMER`,
#: `.../eligos2-v2.1.0/_option_parsers.py`, `README.md:144-151`), the curated
#: inventory (`tools/inventory/curated/tool_inventory_curated.csv`,
#: rows DRUMMER / ELIGOS2_diff / ELIGOS2_solo / m6Anet) and the run scripts of this
#: study (`$RNAMODBENCH_LOCAL/tool_scripts/detection/*/m6Anet.sh` for RNA002,
#: `$RNAMODBENCH_LOCAL/nanopore/user_code/HeLa_WT_004{m6anet,f5c}.sh` for RNA004).
#: DRUMMER and ELIGOS2 take basecalled, aligned BAMs; only m6Anet reads eventalign,
#: which the RNA002 libraries produced with nanopolish and the RNA004 libraries with
#: f5c on the BLOW5 signal.  The upstream NA4 table is deposited unchanged.
REQUIRED_INPUT_FIX = {
    "DRUMMER": "coordinate-sorted, indexed BAM files for both conditions (one or more per "
               "condition) aligned to the same reference FASTA, plus that reference FASTA; "
               "isoform mode also needs a transcript-ID list matching the FASTA headers; "
               "no signal files and no annotation",
    "ELIGOS2_diff": "sorted, indexed BAM files for both conditions, a BED6/BED12 region file "
                    "and a faidx-indexed genome FASTA; the tool refuses to run when the BAM, "
                    "the regions or the reference is missing",
    "ELIGOS2_solo": "sorted, indexed dRNA BAM, a BED6/BED12 region file and a faidx-indexed "
                    "genome FASTA; an optional cDNA BCF serves as the control, otherwise the "
                    "built-in k-mer error model is used",
    "m6Anet": "eventalign table built from a transcriptome-aligned BAM: nanopolish eventalign "
              "on the RNA002 libraries and f5c eventalign on the RNA004 BLOW5 signal; no "
              "control sample, no motif file and no annotation are required",
}


#: because Table S11 used to carry a separate "not recorded" version for each.  They are
#: not: the same source tree (`$RNAMODBENCH_LOCAL/source_code/benchmark/CHEUI/scripts`, commit
#: c1c9ab6 above) serves both modes, and only the bundled checkpoint differs.  The run
#: commands of this study record which file belongs to which mode --
#: `tools/inventory/tables/command_lines.csv` calls
#: `CHEUI_predict_model1.py ... -m .../CHEUI_m6A_model1.h5 -o .../read_level_m6A_predictions.txt`
#: and `CHEUI_predict_model2.py ... -m .../CHEUI_m6A_model2.h5 -o .../site_level_m6A_predictions.txt`,
#: with the same two files per mode for m5C (`CHEUI_m5C_model1.h5` / `CHEUI_m5C_model2.h5`).
#: The cell is extended, never rewritten, so the upstream NA4 wording stays visible; the
#: one local m5C fix recorded in `$RNAMODBENCH_LOCAL/code/code/cheui_m5c_refix/README.md` is a coordinate
#: repair of the conversion chain (HeLa_WT2 lost the +4 centre offset), not a second build.
MODEL_FIX = {
    "CHEUI_m6A": (", the software version above serving both modes: the read-level and site-level "
                  "weights are CHEUI_m6A_model1.h5 and CHEUI_m6A_model2.h5 in m6A mode, and "
                  "CHEUI_m5C_model1.h5 and CHEUI_m5C_model2.h5 in m5C mode"),
}

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
        #: columns are printed, not left to the deposited table.
        
        #: says which mode each pair belongs to (see MODEL_FIX above).
        cell(r.model_checkpoint) + MODEL_FIX.get(str(r.tool_display or r.tool), ""),
        cell(REQUIRED_INPUT_FIX.get(str(r.tool_display or r.tool), r.required_input)),
        cell(r.min_read_coverage),
        cell(r.filtering_parameters),
        cell(r.calling_threshold),
        cell(r.multiple_testing_correction),
        cell(r.default_vs_optimised),
        cell(r.coordinate_harmonisation),
    ])
write_table(
    "TableS6_per_tool_implementation",
    
    
    #: The manuscript sentence that enumerates the same fields stays as it is: that is the
    #: manuscript's own prose, not the table title.
    "Table S6: Per-tool implementation details.",
    ["Tool", "Software / version", "Model / checkpoint", "Required input",
     "Minimum read coverage", "Filtering parameters", "Calling threshold",
     "Multiple-testing correction", "Default or optimised",
     "Coordinate harmonisation"],
    rows,
    "One row per tool and configuration used in this study (13 m6A configurations, the "
"non-m6A tools of the modification panels, and the RNA004 built-in modification models); "
"configurations outside this study are not listed. The calling thresholds of representative "
"tools are swept in Figure S11.",
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

    2026-09-28 (author decision): the ``used_002`` fallback used to answer "basecalling
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


#: tool, so the three BAM-based callers read "raw-signal (Tombo resquiggle)" while the
#: same rows said they ran on RNA004 -- and this table's own note says a resquiggle
#: cannot be produced from RNA004 data.  The input class is what the tool consumes,
#: and the RNA004 runs fed these callers the Dorado basecalls aligned with minimap2
#: (`$RNAMODBENCH_LOCAL/nanopore/user_code/HeLa_WT_004{DRUMMER.sh,ELIGOS.sh,ELIGOS_diff.sh}`),
#: which is the same class the manuscript names ("the BAM-based callers", Experimental
#: Section, "RNA Modification Detection Tools").  m6Anet keeps the eventalign class;
#: the earlier wording called its RNA004 checkpoint a "BAM model", and the checkpoint
#: is the RNA004-optimised one rather than a model over BAM features.
INPUT_CLASS_OVERRIDE = {
    "DRUMMER": "basecalled (BAM pair)",
    "ELIGOS2_diff": "basecalled (BAM pair)",
    "ELIGOS2_solo": "basecalled (BAM)",
    "m6Anet": "raw-signal (eventalign)",
}

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
    rows.append([cell(r.tool_display or r.tool),
                 cell(INPUT_CLASS_OVERRIDE.get(str(r.tool_display or r.tool), r.input_class)),
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
    "One row per tool and configuration used in this study. Input classes: raw-signal (FAST5 "
"plus a nanopolish/Tombo resquiggle) and basecalled (BAM/FASTQ) callers; RNA004 native "
"output is POD5 and Dorado BAM with 9-mer models, so RNA002 signal models cannot be "
"produced for it, and the three BAM-based callers were re-run from Dorado basecalls. "
"Failure and exclusion classes: input format (no raw signal), basecalling model (no "
"RNA004 checkpoint), software maintenance, and outside the scope of this study.",
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
