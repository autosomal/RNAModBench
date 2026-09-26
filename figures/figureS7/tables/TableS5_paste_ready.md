# Table S5 — suggested replacement values (paste-ready)

The submitted `05_submission_work/02_AS_working_copy_and_revisions/Supplementary_Table.pdf` has **no
editable source inside the repository** (verified 2026-09-20 by a repository-wide
search), so this file gives the values to paste into the external layout.

## What changes

Five of the six tools reproduce the published numbers exactly. Only **CHEUI_m5C**
changes, and only through the 2026-09-18 coordinate correction + centre-base hard
filter of the call set (WT union 48,627 → 47,747; ratio 1.0523 → 1.0717). If the
table is refreshed, use the current call-set values below and add a footnote:
*"CHEUI_m5C values use the corrected call set (coordinate fix and reference-base
filter, 2026-09-18); the originally published union of 48,627 WT calls predates
that correction. All other tools are unchanged."*

## Recommended table (raw unions of the three replicates)

| Tool | IVT union | WT union | IVT/WT ratio |
|---|---|---|---|
| NanoNm | 6,399 | 4,181 | 1.5305 |
| NanoMUD-Ψ | 14,794 | 13,711 | 1.0790 |
| NanoMUD-m1Ψ | 43,063 | 40,656 | 1.0592 |
| CHEUI-m5C | 51,171 | **47,747** | **1.0717** |
| NanoSPA-Ψ | 694 | 888 | 0.7815 |
| NanoPsu | 678 | 883 | 0.7678 |

Column semantics (unchanged): union over the three independent HeLa replicates of
the condition; the ratio is unmodified-IVT over WT. The manuscript's 883 / 888 WT
calls and the 0.77 / 0.78 ratios are these values.

## Optional additions

* the **unit-level bootstrap 95 % CI** of the mean-of-counts ratio and the spread
  of the nine WT × IVT replicate pairs — `s6_ratio_ci.tsv` (columns
  `ratio_mean_counts`, `ci_lo`, `ci_hi`, `ratio_pair_min`, `ratio_pair_max`);
* a note that the unmodified IVT libraries are **negative controls, not a
  treatment arm** (the retired "T/WT" label should not reappear).

Source of every number: `tables/s6_counts_summary.tsv`, `tables/s6_ratio_ci.tsv`,
`tables/TableS5_reconciliation.tsv` (all reconciled with the frozen R3-9 evidence
tables; `tables/s6_anchor_check.tsv`, 120/120 OK).
