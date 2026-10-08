# Table S5 — suggested replacement values (paste-ready)

The submitted `$RNAMODBENCH_LOCAL/submission/Supplementary_Table.pdf` has **no
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

**NanoNm** changes for a different reason (2026-10-04, twice extended 2026-10-05): its
wild-type arm **merges** the per-site Nm stoichiometry published with the model with two
of this study's three units, so the row carries three samples -- as every other row does.
The arm is HeLa_WT2 and HeLa_WT3 (3,433 and 1,648 calls) plus the published set; HeLa_WT1
(1,630 calls, the smallest of the four candidates) is left out on the author's instruction
of 2026-10-05. The published set is cut at this study's own thresholds (ratio 0.1 or more
and at least 20 modified reads per site, inside this study's candidate-site set; it is
keyed by `pos - 1`, the sheet being 1-based while every layer of this study is 0-based
BED). It is a single pooled run with no unmodified control and no replicates, so:

* the wild-type **count cell** is the arm's depth-standardised call count, mean (SD) over its
  three samples, read on this row's wild-type candidate-site count (**9,644,510** sites, the
  mean of HeLa_WT2 and HeLa_WT3): **4,230 (1,520) sites** (the rate behind it is 438.6
  (157.6) per million candidate sites; the published sample contributes its 3,948 sites
  inside this study's candidate set against HeLa_WT2's 3,433 and HeLa_WT3's 1,648);
* the unmodified-IVT **count cell** is this study's own three IVT units under the same
  standardisation and the same base: **3,159 (809) sites** (327.5 (83.8) per million), the
  published set having no control;
* the coverage-matched ratio is the quotient of those two cells over the arm's three
  samples (0.747; 2026-10-06: the printed interval is 0.727-0.767, the counting interval of
  `site_boot_ratio`, which is conditional on the libraries and needs no replicate spread
  from the pooled sample; the unit-level intervals stay in the anchors at 0.609-0.966 with
  the pooled sample held fixed and 0.515-1.111 with everything resampled). The footnote
  names the other two arms: 0.962 (0.931-0.994) with this study's three units, the value the
  other five rows print, and 0.529 with the published set alone;
* the arm's **raw** numbers are no longer printed: its union of 6,030 sites and its per-unit
  counts (3,010 with SD 1,207) stay in `s7_nanonm_published_wt.tsv` and `76_anchors.tsv`
  (this study's own three-unit re-analysis alone gives 4,061 sites, the published set alone
  3,948, i.e. 290 against 298 per million candidate sites, 4,061 against 3,948 sites on the
  same base).

2026-10-06 (author decision, second task): the two call columns left the raw per-unit counts
behind -- they carried the depth-standardised rate of each condition (the seven coverage strata
combined by the pooled site weights, mean (SD) over the units) -- and the ratio is their
quotient, so the wild type reads above the unmodified control in every row. 2026-10-07 (author decision):
those two cells are counts of sites again, but not raw counts -- the same standardised rate
read on that row's wild-type candidate-site count (`build_supp_tables.py::base_sites`), which
is what keeps the printed ratio the quotient of the two cells and what Figure 7A draws. The raw per-unit counts and unions showed the unmodified side above the wild
type for four of the six tools (the unmodified libraries carry 1.8-fold the share of
candidate sites above coverage 100), which reads as a contradiction of the ratio; those
numbers are now quoted in the manuscript prose and deposited with the call sets. 2026-10-08
(author decision): the ratio column moved ahead of the two call columns -- the table now reads Tool,
Modification, ratio, wild-type calls, unmodified-IVT calls -- and the two call headers dropped
to "calls per unit" (they keep the Depth-standardised prefix); no value, interval or row moved.

The published set is a reproducibility reference, not a positive reference.
Numbers: `s7_nanonm_published_wt.tsv` (counts), `s7_nanonm_published_matched.tsv`
(per-stratum rates and weights) and `logs/76_anchors.tsv`, all cut by
`76_nanonm_published_wt.py`.

## Retired: the raw-union recommendation (2026-09-20 to 2026-10-05)

2026-10-06: the delivered table no longer prints these columns -- see "What changes" above for the two count cells it prints now. The union numbers below are kept for provenance.

| Tool | IVT union | WT union | IVT/WT ratio |
|---|---|---|---|
| NanoNm | 6,399 | 4,181 *(ours alone; the delivered table prints 6,030, the three-sample arm)* | 1.5305 |
| NanoMUD-Ψ | 14,794 | 13,711 | 1.0790 |
| NanoMUD-m1Ψ | 43,063 | 40,656 | 1.0592 |
| CHEUI-m5C | 51,171 | **47,747** | **1.0717** |
| NanoSPA-Ψ | 694 | 888 | 0.7815 |
| NanoPsu | 678 | 883 | 0.7678 |

Column semantics: union over the three independent HeLa replicates of
the condition; the ratio is unmodified-IVT over WT. Neither is printed any more. The manuscript's 883 / 888 WT
calls and the 0.77 / 0.78 ratios are these values.

## Optional additions

* the **unit-level bootstrap 95 % CI** of the mean-of-counts ratio and the spread
  of the nine WT × IVT replicate pairs — `s7_ratio_ci.tsv` (columns
  `ratio_mean_counts`, `ci_lo`, `ci_hi`, `ratio_pair_min`, `ratio_pair_max`);
* a note that the unmodified IVT libraries are **negative controls, not a
  treatment arm** (the retired "T/WT" label should not reappear).

Source of every number: `tables/s7_counts_summary.tsv`, `tables/s7_ratio_ci.tsv`,
`tables/TableS5_reconciliation.tsv` (all reconciled with the frozen R3-9 evidence
tables; `tables/s7_anchor_check.tsv`, 120/120 OK).
