#!/usr/bin/env python3
"""Create / update the manual curation template.

The template is the **only** file in this directory that is meant to be edited
by hand.  It is generated once and then preserved: re-running this script adds
rows for newly discovered (tool, field) pairs but never overwrites what a human
has already filled in.

Columns
-------
tool_canonical, category, modification, role, field, value, evidence, status,
how_to_obtain

``status`` is one of
* ``CONFIRMED``    -- taken from the manuscript or a tool's own output, safe to cite
* ``NEEDS_CHECK``  -- best current knowledge, please verify before submission
* ``TODO``         -- unknown; the cell must be filled (or explicitly justified
                      as "not applicable") before the response letter is sent

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/tools/inventory/scripts/make_curated_template.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ti_common as tc  # noqa: E402

TEMPLATE = tc.CURATED_DIR / "tool_inventory_curated.csv"
COLS = ["tool_canonical", "category", "modification", "role", "field", "value",
        "evidence", "status", "how_to_obtain"]

#: facts already established by the manuscript or by the tools themselves.
#: (tool, field) -> (value, evidence, status)
KNOWN: dict[tuple[str, str], tuple[str, str, str]] = {
    ("guppy", "software_or_version"):
        ("6.3.9", "manuscript_benchmark/manuscript.tex:174", "CONFIRMED"),
    ("guppy", "required_input"):
        ("POD5/FAST5 raw signal", "manuscript_benchmark/manuscript.tex:174",
         "CONFIRMED"),
    ("Dorado", "software_or_version"):
        ("0.9.1", "result_RNA004/documents/IMPORTANT_INFO.md (dorado-0.9.1 path); " "a 2.0.0 install also exists under $RNAMODBENCH_LOCAL/tool", "NEEDS_CHECK"),
    ("Dorado", "required_input"):
        ("POD5 -> basecalling with --modified-bases-models",
         "result_RNA004/documents/IMPORTANT_INFO.md", "CONFIRMED"),
    ("Dorado", "min_read_coverage"):
        ("10 (applied in this study for the Curlcake FPR analysis, NA-9)",
         "code/revision/common/config.py:DORADO_MIN_COVERAGE", "CONFIRMED"),
    ("Dorado", "calling_threshold"):
        ("percent_modified >= 5 / 10 / 20 / 50 (threshold scan, NA-9)",
         "code/revision/common/config.py:DORADO_PCT_THRESHOLDS", "CONFIRMED"),
    ("Dorado", "multiple_testing_correction"):
        ("not applicable (model probability, no per-site test)",
         "tool documentation", "NEEDS_CHECK"),
    ("Dorado", "default_vs_optimised"):
        ("default model checkpoints; the modification models are explicitly "
         "selected per run", "result_RNA004/documents/IMPORTANT_INFO.md",
         "CONFIRMED"),
    ("Dorado", "coordinate_harmonisation"):
        ("genomic (modkit pileup against GRCh38)", "IMPORTANT_INFO.md; NA-9",
         "CONFIRMED"),
    ("Nanom6A", "calling_threshold"):
        ("0.5", "manuscript (reviewer R1-5 quotes this value)", "CONFIRMED"),
    ("Nanom6A", "required_input"):
        ("basecalled BAM + reference genome (minimap2, -x map-ont / splice)",
         "code_user/detection/*/Nanom6A.sh", "NEEDS_CHECK"),
    ("Nanom6A", "min_read_coverage"):
        ("20 (f5c-mode re-run: predict_sites --support 20); "
         "the original RNA002 run did not record this value",
         "result/f5c/*/logs/pipeline.log", "NEEDS_CHECK"),
    ("DENA", "calling_threshold"):
        ("0.1", "manuscript (reviewer R1-5 quotes this value)", "CONFIRMED"),
    ("DENA", "required_input"):
        ("nanopolish eventalign output; needs a WT + modification-deficient pair",
         "code_user/detection/Hela/IVT/DENA.sh", "NEEDS_CHECK"),
    ("m6Anet", "calling_threshold"):
        ("probability >= 0.5 (tool default)", "m6Anet documentation", "NEEDS_CHECK"),
    ("m6Anet", "min_read_coverage"):
        ("20 (tool default)", "m6Anet documentation", "NEEDS_CHECK"),
    ("Nanocompore", "min_read_coverage"):
        ("5", "result/Nanocompore/*/out_sampcomp.log:min_coverage", "CONFIRMED"),
    ("Nanocompore", "software_or_version"):
        ("1.0.4", "result/Nanocompore/*/out_sampcomp.log:package_version",
         "CONFIRMED"),
    ("Nanocompore", "calling_threshold"):
        ("GMM + KS test p-value (comparison_methods)",
         "result/Nanocompore/*/out_sampcomp.log", "CONFIRMED"),
    ("xPore", "min_read_coverage"):
        ("20", "result/xPore/*.yml:readcount_min", "CONFIRMED"),
    ("xPore", "filtering_parameters"):
        ("readcount_max = 2000000", "result/xPore/*.yml", "CONFIRMED"),
    ("xPore", "multiple_testing_correction"):
        ("FDR on diff_mod_rate (tool default)", "xPore documentation",
         "NEEDS_CHECK"),
    ("DRUMMER", "multiple_testing_correction"):
        ("Benjamini-Hochberg FDR (tool default)", "DRUMMER documentation",
         "NEEDS_CHECK"),
    ("ELIGOS2_diff", "multiple_testing_correction"):
        ("FDR (tool default)", "ELIGOS2 documentation", "NEEDS_CHECK"),
    ("ELIGOS2_solo", "multiple_testing_correction"):
        ("FDR (tool default)", "ELIGOS2 documentation", "NEEDS_CHECK"),
    ("Tombo", "required_input"):
        ("raw FAST5 (resquiggle); used here for coverage and as input to "
         "DRUMMER / MINES / ELIGOS2 / Yanocomp",
         "code_user/detection/*/Tombo.sh", "NEEDS_CHECK"),
    ("Tombo", "multiple_testing_correction"):
        ("FDR (Tombo default)", "Tombo documentation", "NEEDS_CHECK"),
}


#: required input per tool, as evidenced by the command lines that were run.
#: Kept deliberately short and marked NEEDS_CHECK: it states *which data class*
#: a tool consumes, never a numeric value.
INPUT_CLASS: dict[str, str] = {
    "Nanom6A": "basecalled BAM (minimap2, genomic) + reference genome",
    "m6Anet": "nanopolish/f5c eventalign + resquiggle (transcript)",
    "CHEUI_m6A": "nanopolish/f5c eventalign (transcript)",
    "CHEUI_m5C": "nanopolish/f5c eventalign (transcript)",
    "CHEUI-diff": "eventalign output of two conditions",
    "DENA": "eventalign output; requires a WT + modification-deficient pair",
    "DRUMMER": "Tombo-resquiggled FAST5; requires WT vs KO/KD",
    "ELIGOS2_diff": "Tombo-resquiggled FAST5 + matched control",
    "ELIGOS2_solo": "Tombo-resquiggled FAST5 (de novo, no control)",
    "EpiNano_Error": "basecalled FASTQ/BAM + reference (mismatch/error features)",
    "EpiNano_SVM": "basecalled FASTQ/BAM + reference (k-mer + error features)",
    "differr": "basecalled BAM of two conditions + reference",
    "MINES": "Tombo-resquiggled FAST5 (RRACH contexts)",
    "Nanocompore": "eventalign_collapse output of two conditions",
    "NanoMUD": "basecalled BAM (Psi / m1Psi models)",
    "NanoNm": "basecalled BAM + FAST5-derived features",
    "NanoPsu": "basecalled FASTQ/BAM (k-mer features)",
    "NanoSPA": "basecalled FASTQ/BAM (9-mer features)",
    "Tombo": "raw FAST5 (resquiggle)",
    "Tombo_com": "raw FAST5 of two conditions (comparative)",
    "xPore": "f5c/nanopolish eventalign dataprep of two conditions",
    "Yanocomp": "Tombo-resquiggled FAST5 of two conditions",
    "mAFiA": "sorted BAM + reference",
    "f5c": "BLOW5/SLOW5 signal + BAM + reference (eventalign / resquiggle)",
    "nanopolish": "FAST5 + BAM + transcriptome (eventalign)",
    "minimap2": "FASTQ + reference (genome or transcriptome)",
    "samtools": "SAM/BAM",
    "slow5tools": "FAST5 <-> BLOW5/SLOW5 conversion",
    "modkit": "modified-base BAM from Dorado/Guppy",
    "guppy": "raw FAST5/POD5 signal",
}

#: coordinate space in which the tool calls sites, and how it was harmonised
COORD: dict[str, str] = {
    "Nanom6A": "genomic (minimap2 genomic BAM)",
    "m6Anet": "transcript -> genomic via BED12/gtf2bed12 transcript alignment",
    "CHEUI_m6A": "transcript -> genomic via BED12/gtf2bed12 transcript alignment",
    "CHEUI_m5C": "transcript -> genomic via BED12/gtf2bed12 transcript alignment",
    "CHEUI-diff": "transcript -> genomic via BED12/gtf2bed12 transcript alignment",
    "DENA": "transcript -> genomic via BED12/gtf2bed12 transcript alignment",
    "DRUMMER": "transcript -> genomic via BED12/gtf2bed12 transcript alignment",
    "ELIGOS2_diff": "transcript -> genomic via BED12/gtf2bed12 transcript alignment",
    "ELIGOS2_solo": "transcript -> genomic via BED12/gtf2bed12 transcript alignment",
    "EpiNano_Error": "genomic (minimap2 alignments)",
    "EpiNano_SVM": "genomic (minimap2 alignments)",
    "differr": "genomic (minimap2 alignments)",
    "MINES": "transcript -> genomic via BED12/gtf2bed12 transcript alignment",
    "Nanocompore": "transcript (eventalign against the transcriptome)",
    "NanoMUD": "genomic (minimap2 alignments)",
    "NanoNm": "genomic (minimap2 alignments)",
    "NanoPsu": "transcript -> genomic via BED12/gtf2bed12 transcript alignment",
    "NanoSPA": "transcript -> genomic via BED12/gtf2bed12 transcript alignment",
    "Tombo": "transcript (Tombo wiggle) — NOT usable for genomic stratification",
    "Tombo_com": "transcript (Tombo wiggle)",
    "xPore": "transcript -> genomic via BED12/gtf2bed12 transcript alignment",
    "Yanocomp": "transcript -> genomic via BED12/gtf2bed12 transcript alignment",
    "mAFiA": "genomic (minimap2 alignments)",
}

#: tools for which the benchmark did not change any calling parameter
DEFAULT_PARAMS: dict[str, str] = {
    "m6Anet": "default (RNA004 run uses the RNA004-optimised model)",
    "EpiNano_Error": "default", "EpiNano_SVM": "default", "differr": "default",
    "NanoMUD": "default", "NanoNm": "default", "NanoPsu": "default",
    "NanoSPA": "default", "mAFiA": "default", "Tombo": "default",
    "Tombo_com": "default", "MINES": "default", "DRUMMER": "default",
    "ELIGOS2_diff": "default", "ELIGOS2_solo": "default", "Yanocomp": "default",
    "Nanocompore": "default (min_coverage 5 = tool default)",
    "CHEUI_m6A": "default", "CHEUI_m5C": "default", "CHEUI-diff": "default",
}

EVIDENCE_CMD = "code_user/detection/<species>/<sample>/<tool>.sh (see TI2)"


def expand_known() -> None:
    """Add the per-tool input / coordinate / default statements to KNOWN."""
    for tool, value in INPUT_CLASS.items():
        KNOWN.setdefault((tool, "required_input"),
                         (value, EVIDENCE_CMD, "NEEDS_CHECK"))
    for tool, value in COORD.items():
        KNOWN.setdefault((tool, "coordinate_harmonisation"),
                         (value, "post-processing scripts in code_user/*postprocessing",
                          "NEEDS_CHECK"))
    for tool, value in DEFAULT_PARAMS.items():
        KNOWN.setdefault((tool, "default_vs_optimised"),
                         (value, EVIDENCE_CMD, "NEEDS_CHECK"))


def tools_to_curate() -> list[str]:
    """Canonical tools: the catalogue plus everything seen in the raw scan."""
    tools = list(tc.TOOL_META)
    cmd_csv = tc.RAW_DIR / "commands_raw.csv"
    if cmd_csv.exists():
        df = pd.read_csv(cmd_csv)
        for t in df["tool_canonical"].dropna().unique():
            t = str(t)
            if t in tc.TOOL_META and t not in tools:
                tools.append(t)
    return sorted(set(tools))


def build_rows() -> pd.DataFrame:
    rows = []
    for tool in tools_to_curate():
        cat, mod, role = tc.TOOL_META.get(tool, ("unclassified", "n/a", "n/a"))
        for field in tc.FIELDS:
            value, evidence, status = KNOWN.get(
                (tool, field), ("", tc.HOW_TO_OBTAIN.get(field, ""), "TODO"))
            rows.append({
                "tool_canonical": tool, "category": cat, "modification": mod,
                "role": role, "field": field, "value": value,
                "evidence": evidence, "status": status,
                "how_to_obtain": tc.HOW_TO_OBTAIN.get(field, ""),
            })
    return pd.DataFrame(rows, columns=COLS)


def main() -> None:
    logger = tc.setup_logger("make_curated_template")
    tc.CURATED_DIR.mkdir(parents=True, exist_ok=True)
    expand_known()
    fresh = build_rows()

    if TEMPLATE.exists():
        old = pd.read_csv(TEMPLATE).fillna("")
        have = set(zip(old["tool_canonical"], old["field"]))
        add = fresh[~fresh.set_index(["tool_canonical", "field"]).index.isin(have)]
        merged = pd.concat([old, add], ignore_index=True) if len(add) else old
        merged = merged[COLS]
        merged.to_csv(TEMPLATE, index=False)
        logger.info("template exists: kept %d rows, appended %d new rows",
                    len(old), len(add))
    else:
        fresh.to_csv(TEMPLATE, index=False)
        logger.info("template created: %d rows", len(fresh))

    df = pd.read_csv(TEMPLATE)
    logger.info("status counts: %s", df["status"].value_counts().to_dict())
    logger.info("template -> %s", TEMPLATE)


if __name__ == "__main__":
    main()
