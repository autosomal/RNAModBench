#!/usr/bin/env python3
"""Write Table S12 (sensitivity of the site-level tool ranking to coverage).

Source (frozen, read-only): the 2026-09-25 coverage analysis
`04_revision_analysis/R2-1_coverage_rank_stability/tables/`

  * `r21b_threshold_stability.tsv` -- ranking recomputed with a different coverage
    floor (cov >= 5/20/50/100 reads) against the primary floor (cov >= 10);
  * `r21_rank_stability.tsv`       -- ranking recomputed *inside* a coverage
    stratum (5-9, 10-19, 20-49, >= 50 reads per site) against the pooled ranking;
  * `R2-1_reference_and_transcript_strata/tables/r21c_stoich_rank_stability.tsv`
    -- ranking recomputed inside a reference modification-ratio stratum
    (0.1-0.3, 0.3-0.6, > 0.6) against the full high-confidence reference.

Both files report the same statistic: the Spearman correlation between the
per-unit mean rank of the 13 m6A tool configurations and the primary ranking,
with the top-ranked tool of the subset and its BH-FDR.  Only the MCC family is
tabulated here, because MCC is the metric that orders the rows of Fig. 5A.

Usage: conda run -n benchmark-revision --no-capture-output python build_S12.py
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

import pandas as pd

HERE = Path(__file__).resolve().parent
RA = Path(str(_RB / "analysis"))
SRC = (_RB / "analysis/R2-1_coverage_rank_stability/tables")
OUT = HERE / "tables" / "TableS12_coverage_rank_stability.tsv"

ORDER = ["Arabidopsis", "Mouse", "Human"]
UNITS = {"Arabidopsis": "3", "Mouse": "1", "Human": "3"}
FLOORS = ["cov>=5", "cov>=20", "cov>=50", "cov>=100"]
FLOOR_LABEL = {"cov>=5": "coverage floor >= 5 reads",
               "cov>=20": "coverage floor >= 20 reads",
               "cov>=50": "coverage floor >= 50 reads",
               "cov>=100": "coverage floor >= 100 reads"}
STRATA = ["5-9", "10-19", "20-49", ">=50"]
STRATUM_LABEL = {"5-9": "stratum 5-9 reads",
                 "10-19": "stratum 10-19 reads",
                 "20-49": "stratum 20-49 reads",
                 ">=50": "stratum >= 50 reads"}


REF_LEVELS = ["ratio0.1-0.3", "ratio0.3-0.6", "ratio>0.6"]
REF_LABEL = {"ratio0.1-0.3": "reference ratio 0.1-0.3",
             "ratio0.3-0.6": "reference ratio 0.3-0.6",
             "ratio>0.6": "reference ratio > 0.6"}


def fmt_rho(v: float) -> str:
    return f"{v:.3f}"


def fmt_fdr(v: float) -> str:
    if v == 0:
        return "< 1e-16"
    return f"{v:.2g}"


def main() -> None:
    thr = pd.read_csv((_RB / "analysis/R2-1_coverage_rank_stability/tables/r21b_threshold_stability.tsv"), sep="\t")
    strat = pd.read_csv((_RB / "analysis/R2-1_coverage_rank_stability/tables/r21_rank_stability.tsv"), sep="\t")
    ref = pd.read_csv((_RB / "analysis/R2-1_reference_and_transcript_strata/tables/r21c_stoich_rank_stability.tsv"), sep="\t")
    thr = thr[thr["metric"] == "mcc"].set_index(["species", "stratum"])
    strat = strat[strat["metric"] == "mcc"].set_index(["species", "stratum"])
    ref = ref[ref["metric"] == "mcc"].set_index(["species", "stratum"])

    rows: list[list[str]] = []
    for species in ORDER:
        for level in FLOORS:
            r = thr.loc[(species, level)]
            rows.append([species, UNITS[species], "coverage floor", FLOOR_LABEL[level],
                         str(int(r["n_tools"])), fmt_rho(float(r["spearman_rho"])),
                         fmt_fdr(float(r["p_value_fdr_bh"])),
                         "yes" if bool(r["top1_same"]) else "no",
                         str(r["top1_tool_stratum"])])
        for level in STRATA:
            r = strat.loc[(species, level)]
            rows.append([species, UNITS[species], "coverage stratum", STRATUM_LABEL[level],
                         str(int(r["n_tools"])), fmt_rho(float(r["spearman_rho"])),
                         fmt_fdr(float(r["p_value_fdr_bh"])),
                         "yes" if bool(r["top1_same"]) else "no",
                         str(r["top1_tool_stratum"])])
        for level in REF_LEVELS:
            r = ref.loc[(species, level)]
            rows.append([species, UNITS[species], "reference ratio stratum", REF_LABEL[level],
                         str(int(r["n_tools"])), fmt_rho(float(r["spearman_rho"])),
                         fmt_fdr(float(r["p_value_fdr_bh"])),
                         "yes" if bool(r["top1_same"]) else "no",
                         str(r["top1_tool_stratum"])])

    # anchors: the two floors flanking the primary one, and m6Anet's stability
    def val(species: str, axis: str, level: str, col: str) -> str:
        for row in rows:
            if row[0] == species and row[2] == axis and row[3].endswith(level):
                return row[{"spearman_rho": 5, "top1_same": 7, "top1_tool": 8}[col]]
        raise KeyError((species, axis, level, col))

    assert val("Arabidopsis", "coverage floor", "20 reads", "spearman_rho") == "1.000"
    assert val("Mouse", "coverage floor", "20 reads", "top1_tool") == "m6Anet"
    assert val("Human", "coverage floor", "5 reads", "spearman_rho") == "1.000"
    assert val("Arabidopsis", "coverage stratum", ">= 50 reads", "spearman_rho") == "0.876"
    assert val("Human", "reference ratio stratum", "> 0.6", "spearman_rho") == "0.904"
    assert val("Mouse", "reference ratio stratum", "0.1-0.3", "top1_tool") == "m6Anet"
    assert val("Arabidopsis", "reference ratio stratum", "0.1-0.3", "top1_tool") == "MINES"
    n_ref = [r for r in rows if r[2] == "reference ratio stratum"]
    assert len(n_ref) == 9, len(n_ref)

    with OUT.open("w") as fh:
        fh.write("# Table S12: Sensitivity of the site-level tool ranking to transcript coverage "
                 "and to the reference's own composition.\n")
        fh.write("Species\tUnits\tSensitivity axis\tLevel\tn tools\tSpearman rho vs primary\t"
                 "BH-FDR\tTop-ranked tool identical\tTop-ranked tool\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
        fh.write(
            "# note: Ranking statistic MCC (the metric that orders the rows of Fig. 5A), computed "
            "per independent sequencing unit inside each sample's own candidate-site set and "
            "averaged across units, exactly as in the main analysis; each row restricts the "
            "evaluation to the listed level and compares the resulting per-unit mean rank with the "
            "primary ranking (13 m6A tool configurations, candidate-site set with coverage >= 10 "
            "reads, full high-confidence reference). Spearman rho is that rank correlation and "
            "BH-FDR its correction within the axis. Arabidopsis and HeLa contribute three units, "
            "the mouse study one, so no interval is reported for mouse. The three coverage floors "
            "that leave the evaluation set comparable (>= 5, >= 10 and >= 20 reads) give the same "
            "ordering (rho = 0.99-1.00, top-ranked tool unchanged in every species); inside the "
            "informative coverage strata (20-49 and >= 50 reads) the ranking is reproduced at rho "
            "= 0.80-0.98. The sparse strata (5-9 and 10-19 reads) carry almost no information - "
            "the mean recall of every tool is <= 0.9 % there - and are reported as sensitivity "
            "checks only. Stratifying the reference by modification ratio changes the ordering "
            "more strongly (rho = 0.48-0.90, with the top-ranked tool changing in seven of the "
            "nine rows), which is why the primary metrics are reported against the full "
            "high-confidence reference.\n")

    print(f"TableS12: {len(rows)} rows -> {OUT}")


if __name__ == "__main__":
    main()
