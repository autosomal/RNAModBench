# Figure S7 -- legend (rebuilt 2026-09-21, R3-9-centred contract)

> Paste-ready replacement for the Figure S7 entry of the SI legend compilation
> (`$RNAMODBENCH_LOCAL/submission/Supplementary_Figures.pdf`). That compilation
> and `Supplementary_Table.pdf` have no editable source inside the repository, so
> both need an external re-layout pass (the S4 legend and the missing Figure S5 entry
> are outstanding there as well). The figure page itself is
> `figures/FigureS7_rev.pdf`; `$RNAMODBENCH_LOCAL/manuscript/sup/sup6.pdf` is the
> staged copy for the submission tree.

**Figure S7. Non-m6A modification detection against unmodified controls.** Every
quantity is computed per independent sequencing unit (HeLa WT rep1-3, *n* = 3;
HeLa unmodified IVT rep1-3, *n* = 3) inside each sample's candidate-site
universe, and **no m6A-centred reference is used anywhere in this figure** (R2-2):
an m6A-centred reference cannot validate tools that predict other modification
types, so these six tools are characterised by the unmodified synthetic
controls, their own negative controls and replicate-to-replicate consistency
only. Only the unmodified IVT libraries serve as the control arm; they are never
treated as a perturbed or treated sample (E5). **Panels B, C and D show all six
non-m6A tools of the study in the same order**, grouped as in Fig. 7A
(false-positive-dominated: CHEUI-m5C, NanoMUD-Ψ, NanoMUD-m1Ψ; intermediate:
NanoNm; specific but sparse: NanoPsu, NanoSPA-Ψ). **Panel A shows only the five
tools that were actually run on the synthetic constructs**: CHEUI-m5C was never
run on them, so its column is absent rather than drawn blank (R1-6). Values
behind every panel, with the frozen-table reconciliation, are in `../tables/`.
Following the journal-wide house style, the figure carries no titles of any
kind: the four panels are identified here and by their bold letters (A, B, C, D
from top to bottom), and the tools by their axis labels. B, C and D share one
six-column tool grid, with the tool names drawn at 45 degrees under each of
them; A has its own five-column axis, so A's columns do not sit under B's. The
page is 170 x 240 mm, i.e. the printed width of this journal's supplementary
pages.

**(A)** Unmodified Curlcake controls: calls inside each construct's candidate
universe, per construct, for the five tools that were run on the synthetic
constructs. Filled dots are the two independent constructs
(*Curlcake_IVT_rep1*, *rep3*), the open dot is the depth-matched subset
(*Curlcake_IVT_rep2_partial*, excluded from the mean), and the bar is the mean of
the independent constructs. The five tools that were run on the constructs span
three orders of magnitude: NanoMUD-m1Ψ reports 89 and 93 calls (18,100 and 18,899
per 10⁶ candidate sites), NanoNm 36 and 39 (3,598 and 3,903 per 10⁶), whereas
NanoPsu and NanoSPA-Ψ report 0 and 1 call (0 and 203 per 10⁶) and NanoMUD-Ψ 0 and
3 (0 and 610 per 10⁶). The synthetic constructs carry ~5,000 candidate sites each
(10,135 bp), so the call counts and the per-10⁶ densities are proportional;
exact counts, per-10⁶ and per-10-kb densities are tabulated. CHEUI-m5C was never
run on the synthetic constructs and therefore has no column in this panel at
all (R1-6); its negative-control behaviour is documented by the HeLa
unmodified-IVT libraries in panels B-D instead.

**(B)** Overlap at the level of individual replicate pairs (Jaccard index, log
axis). Circles are the mean pairwise overlap of the three pairs within HeLa WT
(blue) and within unmodified IVT (orange); the grey circles are the nine
WT × IVT pairs, each computed on the shared candidate universe of that pair;
bars give the mean of each group. For the false-positive-dominated tools the
between-condition overlap sits at the level of the within-replicate overlap
(CHEUI-m5C: 0.0125-0.0148, i.e. its calls carry no condition-dependent signal),
whereas the two Ψ tools remain clearly below their own replicate overlap
(NanoPsu and NanoSPA-Ψ: 0.063-0.088 between conditions versus 0.125-0.283 within
replicates). The set-level global Jaccard quoted in the main text (CHEUI-m5C
5.2 × 10⁻⁴ in WT, 6.8 × 10⁻⁴ in IVT) is tabulated
(`s7_jaccard_within_cross.tsv`, `s7_jaccard_pairs.tsv`).

**(C)** Where each tool places its calls inside its own 0-1 score range -- the
mechanism behind the failures R3-9 asks about. Every mark is one independent
sequencing unit: the point is that unit's median score and the whisker spans its
interquartile range (q25-q75); filled marks are HeLa WT and open marks the
unmodified-IVT libraries, three units each (the colours and the filled/open
convention of panels B and D). Each tool is drawn on its own native score -- the
modification ratio for CHEUI-m5C and NanoNm, the reported probability for the
other four -- so the strip is *not* a separation statistic: the pooled per-tool
AUC that summarises score behaviour is the quantity already drawn in Fig. 7F,
and the per-unit-pair AUCs are tabulated. Four signal-level reasons are visible
in this one strip.
(i) *No dynamic range:* NanoMUD-Ψ and NanoMUD-m1Ψ report their maximum
probability for the typical call (median 1.000, interquartile range 0.997-1.000
and 1.000-1.000, and 100 % of their calls at ≥0.9) in *both* conditions, so
their score cannot separate anything.
(ii) *A narrow band:* NanoPsu and NanoSPA-Ψ report 0.965 (0.955-0.980) for
essentially every call in both conditions.
(iii) *Overlap:* NanoNm's modification ratio sits low and the two conditions
overlap (WT unit medians 0.204-0.224 versus unmodified IVT 0.180-0.188).
(iv) *Inversion:* CHEUI-m5C scores the **unmodified** control at or above the
wild type (WT unit medians 0.316-0.444 versus unmodified IVT 0.414-0.842), and
its spread across units is the widest of the six tools -- its nine unit-pair
AUCs run from 0.00004 to 0.552 (mean 0.166) against 0.50-0.56 for every other
tool.
The reproducibility half of the same critique is quantitative as well: 97.2 % of
the 46,198 in-universe sites of the pooled CHEUI-m5C WT set are supported by a
single unit (unmodified IVT 97.6 %; the other five tools 56-88 %,
`../tables/s7_replicate_support.tsv`; the pooled set serves to describe this
structure only and is not a replicate-level quantity). Numbers behind this panel: medians and
quartiles in `s7_score_location_per_unit.tsv` (drawn here), per-unit-pair AUCs
and KS statistics in `s7_score_separation.tsv`, per-unit score densities
(not plotted) in `s7_score_density_per_unit.tsv`.

**(D)** HeLa calls per independent unit (dots, log axis: blue filled = WT, orange
open = unmodified IVT) with the union of the three units of each condition
(dashed lines, each drawn on its own condition's half of the tool slot, shown for
reference only), and,
below, the ratio of calls on
the unmodified IVT libraries to calls on WT (linear axis, reference line at one --
one means as many calls on unmodified RNA as on WT). Circles with whiskers are
the mean-of-counts ratio with its unit-level percentile bootstrap 95 % confidence
interval (B = 1000, seed = 20260920); the coverage-matched ratio, computed inside
matched coverage strata, is tabulated in Table S5 and is the quantity the text
reports. The mark of the ratio strip
is defined here rather than by an in-figure key: the compressed strip leaves no
empty corner for a legend. NanoPsu and NanoSPA-Ψ report 201 ± 77 and
204 ± 81 calls per unit in WT (coverage-matched IVT/WT ratios 0.71 and 0.72,
95 % CI 0.64-0.79 in both), whereas CHEUI-m5C reports 15,845 ± 3,827 calls per
unit in WT and 17,129 ± 7,910 on the unmodified control (coverage-matched ratio
0.79, 95 % CI 0.49-1.14; the raw unions of the three replicates, 47,747 and
51,171, are listed in Table S5 for reference only). CHEUI-m5C's
union departs from the published Table S5 value (48,627) because that value
predates the coordinate correction and the reference-base filter of the call set;
five of the six tools reproduce Table S5 exactly
(`TableS5_reconciliation.tsv`). The confidence intervals resample three units per
condition and are therefore indicative, not precise.

---

## Definitions and provenance

* independent units: `HeLa_WT1..3` and `HeLa_IVT_rep1..3`
  (`harmonisation/manifest/sample_registry.csv`); Curlcake independent constructs =
  `Curlcake_IVT_rep1` and `Curlcake_IVT_rep3`, with `Curlcake_IVT_rep2_partial`
  reported separately because it is a depth-matched subset of `rep3`.
* call sets: `harmonisation/callsets` (BED, 0-based; co-duplicated rows collapse to
  one site, keeping the best score).
* candidate universe: `harmonisation/universe/RNA002/...` filtered to coverage >= 10
  and to reference bases that can carry the modification (Nm = A/C/G/T,
  Ψ/m1Ψ = T/A, m5C = C/G).
* tables: `s7_curlcake_per_construct.tsv`, `s7_jaccard_pairs.tsv`,
  `s7_jaccard_within_cross.tsv`, `s7_replicate_support.tsv`,
  `s7_score_location_per_unit.tsv` (panel C), `s7_score_separation.tsv`,
  `s7_score_density_per_unit.tsv`, `s7_counts_per_replicate.tsv`,
  `s7_counts_summary.tsv`, `s7_ratio_ci.tsv`, `s7_anchor_check.tsv`
  (552/552 values reconciled with the frozen R3-9 tables at 1e-5),
  `TableS5_reconciliation.tsv`.
* the earlier pooled density table (`s7_score_density.tsv`) was retired on
  2026-09-21: pooling the replicates of a condition into one curve contradicts
  the per-unit rule of this study. It is kept as
  `s7_score_density_RETIRED_pooled.tsv` for provenance only.

## Caveats to keep with the numbers

1. The unit-level bootstrap resamples only three independent units per condition;
   the intervals are indicative, not precise, and the nine WT × IVT pairs are
   shown individually in panel B for the same reason.
2. No external reference is used in this figure. The truth-anchored enrichment
   against the database compilation (RMBase + DirectRMDB), the orthogonal NGS
   reference and the CHEUI authors' own *E. coli* data set are shown in Fig. 7
   and its evidence tables and are deliberately not duplicated here.
3. The Curlcake constructs are short synthetic RNAs; their false-positive density
   is a per-candidate-site rate, not a genome-wide estimate.
4. In-universe WT unions behind panel A/D: NanoPsu 488 and NanoSPA-Ψ 497; the
   883 / 888 quoted in the main text are the raw unions of three replicates.
