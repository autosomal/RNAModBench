# Figure S9 (revised) — legend, provenance and audit

File: `figures/FigureS9_rev.pdf` (+ 300 dpi `FigureS9_rev.png`), page
**1152 x 864 pt** (identical to the page of the file it replaces), Arial
embedded (cairo_pdf, `pdffonts`-verifiable), no gridlines, minimum drawn font
12 pt (axis text 12.5-15 pt, panel titles 19 pt, three bold block letters 24 pt,
no figure number). This is the replacement for `$RNAMODBENCH_LOCAL/submission/02_AS_working_copy_and_revisions/sup/sup9.pdf`;
every replacement kept a `.bak_<timestamp>` copy of the previous file
(`.bak_20260921_0157` = original Illustrator version, `.bak_20260921_0215` =
six-panel interim version, `.bak_20260921_0445_preABC` = the single **S9**
figure-number version that preceded the three block letters).

## Legend (paste-ready)

> **Figure S9. mRNA region distribution and false-positive behaviour of the
> Dorado other-modification models on the RNA004 HeLa libraries.** Metagene
> density of the calls reported by six Dorado built-in modification models —
> four pseudouridine models (hac@v5.0.0_pseU, hac@v5.1.0_pseU,
> sup@v5.0.0_pseU, sup@v5.1.0_pseU) and two 5-methylcytosine models
> (hac@v5.1.0_m5C, sup@v5.1.0_m5C), drawn as the six density panels of **block A**
> (read left to right, top to bottom; the page letters the three blocks, not the
> individual panels; **B** = the bottom-left block and **C** = the bottom-right
> block) — in the wild-type HeLa library
> (blue, WT) and in the unmodified in vitro transcribed library (orange,
> unmodified IVT; each panel carries its own key). Model names are written
> `<mode>@<basecaller version>_<modification>` (for example hac@v5.1.0_pseU); the
> modification-model version is `@v1` in every case, so the six full Dorado model
> ids are hac@v5.0.0_pseU@v1, hac@v5.1.0_pseU@v1, sup@v5.0.0_pseU@v1,
> sup@v5.1.0_pseU@v1, hac@v5.1.0_m5C@v1 and sup@v5.1.0_m5C@v1. Every call is taken from the
> analysis-ready call-set layer (one call per genomic position, the sample's own
> strand) at the delivered operating point of the RNA004 Dorado models
> (modified-read fraction >= 90 %), and is placed on the **Ensembl** GRCh38p14
> release 112 protein-coding transcript model (1 kb upstream flank, 5'UTR, CDS,
> 3'UTR, 1 kb downstream flank, transcript direction); the density is computed
> with the Bioconductor *Guitar* package (site intervals padded by 1 bp, three
> equidistant points per interval, Gaussian smoothing), y-axes are scaled per
> panel. **No m6A-centred or other positive reference is used in this figure;
> the unmodified IVT library is the negative control.** In every panel the
> unmodified control reproduces the 3'UTR-dominated profile of the wild type
> (fraction of the transcript-body calls that fall in the 3'UTR: 53.7-91.8 % on
> IVT versus 57.6-84.7 % in WT), so the transcript-position preference of these
> models is not modification-specific; the largest WT-minus-IVT difference is
> +14.0 pp (sup@v5.1.0_m5C, last of the six), in three of the six models it is within +/-5 pp,
> and in **all six models the 95 % site-level bootstrap interval of the
> difference spans zero** (percentile bootstrap, B = 1000, seed = 20260920,
> resampling the transcript-body calls of one library — call sampling, not
> replicate uncertainty). WT and unmodified IVT are **one sequencing unit each
> (*n* = 1)**; no biological replication exists for this library pair, so no
> replicate spread is shown. The two inosine models of the same family are not
> drawn here because they report too few calls (15 and 32 in WT, 0 and 6 in
> IVT); their counts are reported in Figure S8A.
>
> **Blocks B (bottom left) and C (bottom right).** Two complementary
> characterisations of the same six models, side by
> side. **B** — false-positive density on the unmodified Curlcake control
> (4 synthetic constructs, 10,135 bp — every call there is a false positive) as
> a function of the modified-read threshold (5-90 %), one curve per model on a
> log y-axis; HAC calls are solid and sup calls dashed, v5.0.0 purple circles
> and v5.1.0 green triangles, and hollow symbols mark zero calls (drawn at the floor of the log
> axis, 0.5 per 10 kb; one call = 0.987 per 10 kb). The density falls steeply
> with the threshold — 69.07-242.72 per 10 kb at 5 %, 1.97-6.91 at 50 %
> (between the DRACH-specific and the general m6A models of Fig. 8B) and zero
> for every model at 90 %. **C** — cumulative distribution of each model's
> own reported modified fraction in the wild-type (blue) and unmodified HeLa
> IVT (orange) libraries, one facet per model; the delivered call sets already
> satisfy >= 90 % modified read, so the axis spans 0.90-1.00. The two
> distributions overlap almost completely (Mann-Whitney AUC 0.5036-0.6074,
> except sup@v5.0.0_pseU, where the unmodified control scores *higher*: AUC
> 0.2394 with 63.4 % of its calls saturated at >= 0.99), i.e. the models' own
> confidence does not separate modified from unmodified RNA.

## Data provenance and rules (internal record / response letter)

1. **Inputs.** `data/callsets/RNA004/Human/`
   `RNA004_HeLa_{WT,IVT}/{Psi,m5C}/<model>/<sample>.tsv` — the analysis-ready
   layer written by `harmonisation/scripts/34_export_callsets.py`. The metagene
   panels are drawn from the Guitar BEDs exported from that layer by
   `harmonisation/scripts/21b_export_guitar_bed.py` (own strand preserved, intervals
   padded 1 bp on each side, chromosome labels in the Ensembl GTF spelling,
   `bed/RNA004/majority/<group>/<mod>/<model>.bed`). With a single sequencing
   unit per condition the "majority" merge is that unit itself — the labels
   never claim a consensus of replicates.
2. **Block B (bottom left)** reads `../tables/figS10_curlcake_scan.tsv`, computed from
   `callsets/RNA004/Curlcake/RNA004_Curlcake_IVT/{Psi,m5C}/<model>/`
   (the synthetic control is **not** pre-filtered, so the threshold can be
   swept). The Curlcake run names the same model families differently from
   HeLa (v5.0.0: `Dorado_hac@v5.0.0_pseU_m6A_Psi` / `..._sup@v5.0.0_pseU_m6A_Psi`;
   v5.1.0: `Dorado_{hac,sup}@v5.1.0_all_Psi` and `..._all_m5C`; there is no plain
   `pseU@v1` model), so the curves are labelled with the HeLa model names. Raw call numbers and both denominators are pinned to the frozen
   evaluation table at 1e-5 relative tolerance.
3. **Block C (bottom right)** reads `../tables/figS10_score_validity.tsv` (one row per
   call) and `..._summary.tsv` (one row per model, written by
   `scripts/60_figS10_tables.py`); the AUC is the Mann-Whitney common-language
   effect size between the wild-type and unmodified-IVT calls, with the verdict
   rule < 0.4 inverted, 0.4-0.6 no discrimination, > 0.6 wild type higher.
4. **What was wrong with the published sup9.pdf.** It was drawn by
   `result_RNA004/scripts/R/dorado_other_mods_guitar_wt_ivt.R`, which forced
   every BED line to strand `"+"` (minus-strand sites read on the plus-strand
   coordinate), used one pooled curve per condition without any call-set
   provenance, and carried no numbers — which is why its legend could only say
   "the distribution of m5C and pseU modifications". The revised panels are
   strand-aware, traceable to per-unit call sets and quantified
   (`tables/figS10_region_shares.tsv`, `tables/figS10_wt_ivt_contrast.tsv`,
   `tables/figS10_curlcake_scan.tsv`, `tables/figS10_score_validity_summary.tsv`).
5. **Density kernel.** The Bioconductor *Guitar* package itself
   (`samplePoints -> normalize -> .generateDensity_CI`, `enableCI = FALSE`,
   `CI_ResamplingTime = 20`; identical to `= 1000` while CI is off — probe in
   `figures/figureS1/logs/figS1_probe.log`). Component rectangles/labels mirror
   `Guitar:::.RNAPlotStructure`; both flanks are labelled `1kb`.
6. **Definitions used in the numbers.** `3'UTR share` = share of the calls
   placed in the 5'UTR + CDS + 3'UTR body (the 1 kb flanks are tabulated
   separately); `WT - IVT difference` is in percentage points; the
   "reproduced on unmodified IVT" verdict is the pre-registered rule
   |delta 3'UTR share| <= 5 pp. `coverage >= 10` replicates of every share are
   in the tables (they move no share by more than 4 pp). Block B:
   1 call = 0.987 per 10 kb; block C: score = the model's own reported
   modified fraction.
7. **Reviewer mapping.** R2-2/E1 (no m6A-centred reference for non-m6A tools;
   unmodified IVT is the control), R3-9 (non-m6A calls are background-dominated:
   the distribution is reproduced on the control and the models' own scores do
   not separate the conditions), R3-8/E8 (quantified, threshold-resolved
   false-positive metric and moderated wording instead of a qualitative
   "distribution of modifications"), R3-2/E6 (sequencing-unit structure stated,
   single unit admitted), R1-6 (inclusion/exclusion: why only six models, why
   inosine is absent), E10/R3-m4 (legend covers every panel and both axes).

## Numbers behind the figure

Block A — six density panels, top two rows (WT and unmodified IVT). The `slot`
column gives each panel's position inside block A (R1C1 ... R2C3, left to right,
top to bottom); no per-panel letter is printed anywhere:

| slot | model | WT calls (raw / on mRNA) | IVT calls (raw / on mRNA) | 3'UTR share WT [95 % CI] | 3'UTR share IVT [95 % CI] | delta (pp) | delta 95 % CI (pp) |
|---|---|---|---|---|---|---|---|
| R1C1 | hac@v5.0.0_pseU | 279 / 94 | 138 / 56 | 84.7 % [77.6, 92.9] | 91.8 % [83.7, 98.0] | -7.1 | [-18.0, +4.6] |
| R1C2 | hac@v5.1.0_pseU | 235 / 62 | 70 / 27 | 84.3 % [74.5, 92.2] | 90.9 % [77.3, 100.0] | -6.6 | [-21.6, +8.4] |
| R1C3 | sup@v5.0.0_pseU | 549 / 221 | 953 / 175 | 83.2 % [78.1, 87.8] | 85.0 % [79.1, 90.2] | -1.8 | [-9.5, +6.3] |
| R2C1 | sup@v5.1.0_pseU | 610 / 201 | 233 / 96 | 81.2 % [75.6, 86.4] | 77.0 % [69.0, 85.1] | +4.2 | [-5.5, +14.5] |
| R2C2 | hac@v5.1.0_m5C | 725 / 396 | 261 / 136 | 57.6 % [52.4, 62.2] | 54.5 % [45.5, 63.6] | +3.0 | [-7.0, +13.3] |
| R2C3 | sup@v5.1.0_m5C | 218 / 118 | 97 / 52 | 67.6 % [57.8, 75.5] | 53.7 % [36.6, 68.3] | +14.0 | [-2.6, +31.6] |

Intervals: site-level percentile bootstrap, B = 1000, seed = 20260920 (unit =
one called site within the single library of that condition).

Block B (bottom left) — threshold scan on the unmodified Curlcake control (FP = calls on
the unmodified construct; 10,135 bp):

| slot | Curlcake model | calls @5 % | @50 % | @90 % | FP per 10 kb @5 % | @50 % |
|---|---|---|---|---|---|---|
| R1C1 | Dorado_hac@v5.0.0_pseU_m6A_Psi | 89 | 2 | 0 | 87.81 | 1.97 |
| R1C2 | Dorado_hac@v5.1.0_all_Psi | 107 | 5 | 0 | 105.57 | 4.93 |
| R1C3 | Dorado_sup@v5.0.0_pseU_m6A_Psi | 70 | 3 | 0 | 69.07 | 2.96 |
| R2C1 | Dorado_sup@v5.1.0_all_Psi | 111 | 7 | 0 | 109.52 | 6.91 |
| R2C2 | Dorado_hac@v5.1.0_all_m5C | 246 | 4 | 0 | 242.72 | 3.95 |
| R2C3 | Dorado_sup@v5.1.0_all_m5C | 106 | 2 | 0 | 104.59 | 1.97 |

Source: `../tables/figS10_curlcake_scan.tsv` (all ten thresholds per model).

Block C (bottom right) — score validity in HeLa (does the model's own score separate the
conditions?):

| slot | model | WT / IVT calls | median score WT / IVT | >= 0.99 WT / IVT | median coverage WT / IVT | AUC | verdict |
|---|---|---|---|---|---|---|---|
| R1C1 | hac@v5.0.0_pseU | 279 / 138 | 0.929 / 0.930 | 6.5 % / 14.5 % | 135 / 69 | 0.4767 | no discrimination |
| R1C2 | hac@v5.1.0_pseU | 235 / 70 | 0.933 / 0.925 | 5.5 % / 1.4 % | 174 / 48 | 0.5799 | no discrimination |
| R1C3 | sup@v5.0.0_pseU | 549 / 953 | 0.949 / 1.000 | 14.8 % / 63.4 % | 145 / 6 | 0.2394 | inverted |
| R2C1 | sup@v5.1.0_pseU | 610 / 233 | 0.935 / 0.921 | 5.6 % / 2.6 % | 112 / 50 | 0.6074 | wild type higher |
| R2C2 | hac@v5.1.0_m5C | 725 / 261 | 0.930 / 0.926 | 5.9 % / 5.4 % | 65 / 60 | 0.5486 | no discrimination |
| R2C3 | sup@v5.1.0_m5C | 218 / 97 | 0.932 / 0.933 | 5.1 % / 5.2 % | 79 / 102 | 0.5036 | no discrimination |

Source: `../tables/figS10_score_validity_summary.tsv` (p-values and full summary
in the table).

Reported values — strict false-positive rate on the unmodified HeLa IVT control
at the delivered operating point (source data of the text and reply letter;
reported values, not drawn):

| slot | model | IVT calls (in universe) | FP per 10 kb | FP per 10^6 candidates |
|---|---|---|---|---|
| R1C1 | hac@v5.0.0_pseU | 138 (97) | 0.0087839 | 8.724 |
| R1C2 | hac@v5.1.0_pseU | 70 (54) | 0.0044556 | 4.425 |
| R1C3 | sup@v5.0.0_pseU | 953 (243) | 0.06066 | 60.249 |
| R2C1 | sup@v5.1.0_pseU | 233 (187) | 0.014831 | 14.730 |
| R2C2 | hac@v5.1.0_m5C | 261 (216) | 0.016613 | 18.438 |
| R2C3 | sup@v5.1.0_m5C | 97 (87) | 0.0061742 | 6.852 |

Source: `../tables/figS10_ivt_fpr.tsv`, cross-checked row by row against the
frozen evaluation table `harmonisation/evaluation/tables/controls_ivt_fpr.tsv`
(the same table that feeds Fig. 8B/D) at 1e-5 relative tolerance.

All numbers above: `../tables/figS10_region_shares.tsv`,
`../tables/figS10_wt_ivt_contrast.tsv`, `../tables/figS10_key_numbers.md`.

## Audit

`../logs/61_verify_figS10.log` — **55/55 checks pass**: page size equals the
replaced file, Arial-only
embedded fonts, no grid element in the drawing script, minimum font >= 7 pt, the
six-model white list with the inosine models excluded and absent from the PDF
text, BED line counts and de-duplicated `callsets` counts equal the frozen
anchors, every density panel in block A at its own grid slot, region shares sum
to 1, the +/-5 pp verdict rule holds, every number
quoted here resolves to the frozen tables, the uncertainty layer reproduces
bit-for-bit from the frozen seed, the Curlcake threshold profiles are monotone
and hit their frozen call anchors, the score-validity AUCs reproduce
independently, and no FPR or AUC value is printed inside the figure, no two text elements on the page overlap or come closer than 2 pt, the page carries the three block letters A, B, C exactly once each, at the top-left of their own block (A on the six Guitar panels, B on the threshold scan, C on the score validity), with no figure number, and the bottom-left half uses the purple/green version pair.
