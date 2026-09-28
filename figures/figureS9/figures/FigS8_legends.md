# Figure S8 (rebuilt) - legend, provenance and audit

**Delivered as Figure S9.**  This directory is the analysis-numbered Figure S8
workflow; the manuscript and the delivered Supplementary Figures number this page
**Figure S9** (md5 of `$RNAMODBENCH_LOCAL/submission/new_submission/sup/FigureS9_rev.pdf` equals
`figures/FigureS8_rev.pdf`).  The caption printed in the SI is the `Figure S9.`
block of `tables/sup_figure_captions.md`.

**Figure S9. Detection counts, explicit false-positive control, window-resolved
accuracy and calibration of the RNA004 calling.**
The figure uses the RNA004 chemistry only.  Every call made on an unmodified
control is a false positive by construction, so the control panels report the
false-positive rate with its denominator and the matching specificity, and the
statistical panels report association, agreement, effect size and rank stability
rather than claims of quantitative accuracy.

**Layout (2026-09-27, user).**  Nine cells on one A4 **landscape** sheet
(842.4 x 595.44 pt): three columns of 3.740 in and three rows of 3.400 / 2.195 /
2.195 in (margins 0.10 in, gutters 0.14 in), so A | B | C sit on the first row,
D | E | H on the second and F | G | I on the third.  Every panel draws its data in
a **closed box** (top and right spines visible) that **fills its own cell**: the
box runs from the common label column to the cell's right margin and from its
bottom offset to 0.075 in below the cell top.  **All nine cells share one label
column** (1.638 in, set by A's longest model name at 7 pt and by the 0.23 in of
clearance F's longest tool name needs for the 24 pt panel letter), so all nine
boxes start and end on the same two vertical lines; their heights are uniform per
row (2.620 in on row 1, 1.400 in on rows 2 and 3) and their bottom offsets are
0.414 in (row 1, which carries no key) and 0.720 in (rows 2 and 3, whose strip
carries the x title and the panel's key).  A cell-filling box is a wide rectangle:
a square box in a 3.740 x 2.195 cell would leave ~2.4 in of white on its right.
The
sheet is composed 1:1 by `71_figs8_page.py`
(`common.panelpage.compose_page` asserts that every piece stays inside the page
and that no two pieces overlap), so the point sizes in the panel PDFs are the
printed sizes.  Arial only, no grid, no panel titles (the sample names live in
this legend), minimum drawn type 7 pt.  **A key prints under the panel it belongs
to**: the species key of the two stability cells prints under the box of each of
them (H and I) and the key of the two effect-size facets prints under G's box.  The
thirteen-entry key of the window sweep is the one departure -- thirteen entries do
not fit in a single cell's strip -- so each facet prints one half of it in the
strip under its own box (D the first six entries, E the remaining seven) and the
two halves read as one band under the pair.

**Definitions.** *False positives per 10 kb* = reported sites / mappable control
sequence (Curlcake 10,135 bp).  *False positives per 10^6 candidate adenosines* =
reported sites / the 4,928 Curlcake candidate adenosines (annotated exons,
coverage >= 10x, modification-compatible base); *specificity* = 1 - calls /
candidates.  *Exact-nucleotide fraction* = share of the window-matched calls whose
single-nucleotide position matches a GLORI site.  *Effect size* = OLS slope of the
tool modification ratio on the GLORI ratio at shared sites (1.0 = proportional,
> 1 = systematic over-estimation).  *Spearman rho vs. primary ranking* (H, I) = the
rank correlation of the 13 configuration ordering recomputed inside one coverage
level or one reference-ratio stratum against the primary ordering of the
manuscript (candidate-site set with coverage >= 10 reads, full high-confidence
GLORI reference), computed per independent sequencing unit and averaged across
units; the numbers are those of the delivered Table S12.

**A.** Detected sites in the HeLa RNA004 libraries as grouped horizontal bars on a
log axis: for every entry the wild-type count (blue) and the unmodified-IVT count
(orange) are drawn as one pair, an entry with no call carries an open left
triangle on the axis floor instead of a zero-length bar (a measured zero is the
open circle).  Rows hold the eight Dorado m6A models, the five other m6A tools and
the ten models for other modification types (Dorado pseU x 4, m5C x 2, inosine x 2,
NanoPsu, NanoSPA_psU), sorted by wild-type count inside each block.  m6Anet
reported the most wild-type sites (30,327); in the unmodified IVT library the
DRACH-specific models fell to 2-10 sites against 13,028-46,963 for the non-DRACH
m6A models.

**B.** The eight ORCA channels on the same HeLa RNA004 library, one bar pair per
channel (wild type in blue, unmodified IVT in orange): m5C 14,784/20,580, m6A
8,389/10,190, Nm 3,389/1,873, pseU 533/568, m1A 507/758, inosine 252/128, m7G
53/55 and m6Am 39/58.  ORCA reports several modification types from one run, so its
channels are a facet of their own rather than rows of A.

**C.** Explicit false-positive control on the unmodified RNA004 Curlcake library,
at the 50 % modified-read operating point of the Dorado models: false positives
per 10^6 candidate adenosines (log axis) with the matching specificity on the top
axis.  The DRACH-specific models produced 0-4 false positives (0.0-3.9 per 10 kb;
two models with no call at all), the non-DRACH m6A models 10-13 (9.9-12.8 per
10 kb), the inosine + m6A models 5-9, and of the two m6A tools run on this library
NanoSPA_m6A 17 (16.8 per 10 kb) and m6Anet 5 (4.9 per 10 kb); the full range is
99.66-100 % specificity.  The rate is not zero and it is strongly
threshold-dependent (main-text Fig. 8B), so no statement of perfect
false-positive control is made anywhere in the figure.

**D.** Matching window versus detection quality on the HeLa RNA004 library, for
every one of the thirteen entries measured there: left cell, positive predictive
value against GLORI; right cell (E), the exact-nucleotide fraction; the dashed
guide marks the 2-bp working point of the primary analysis.  Every entry has a
colour of its own -- no two curves share a hue -- and every entry is named in the
thirteen-entry key printed in the strip under the two sweep cells.  Reading
0 -> 50 bp: the
best DRACH model 82.08 -> 85.73 % PPV while its exact fraction falls 1.000 ->
0.958; m6Anet 47.63 -> 52.68 % with 1.000 -> 0.904; ELIGOS2_solo 0.43 -> 14.29 %
with 1.000 -> 0.030.  A wide window buys sensitivity at a measured cost in
single-nucleotide localisation, which is the resolution-sensitivity trade-off the
manuscript reports.

**F.** Modification-ratio agreement per tool: only the nine tools whose score is a
ratio are drawn (ELIGOS2_solo, ELIGOS2_diff and NanoSPA_m6A have no ratio), one
lollipop per tool with the dot marking the Pearson r against GLORI and the whisker
its Fisher-z 95 % confidence interval.  m6Anet r = 0.653 (0.643-0.662, n =
14,559); the DRACH-specific models r = 0.37-0.45.

**G.** The agreement and calibration of the same nine tools: the least-squares
calibration slope of the tool ratio on the GLORI ratio (filled dot) and Lin's
concordance correlation coefficient (open dot), each with its bootstrap 95 %
confidence interval, the two estimates dodged within the row.  m6Anet slope 0.612
(0.599-0.623) and CCC 0.645 (0.634-0.655); the DRACH-specific models slope
0.20-0.24 and CCC 0.15-0.19.  Every tool sits well below proportional tracking
(a slope of 1), i.e. the tools rank sites consistently but compress the dynamic
range of the ratio.  The key of the two facets prints in one row under G's own
box (F and G draw the same two handles; the pair shares the row).

**H.** Rank stability of the tool ordering against transcript coverage (added
2026-09-27 for reviewer point R2-1).  One row per coverage level -- the comparable
floors (>= 5 and >= 20 reads) and the informative strata (20-49 and >= 50 reads) --
one marker per species (Arabidopsis open circle, mouse square, HeLa triangle) at
its Spearman rho against the primary ranking; the dashed guide marks rho = 1.0.
The ordering is reproduced at rho = 0.80-1.00 and 0.995-1.000 at the comparable
floors, and the top-ranked tool stays m6Anet in every species at those floors.  The
sparse strata (5-9 and 10-19 reads) are not drawn: the mean recall of every tool
is <= 0.9 % there, so they carry no ranking information (Table S12).  The three
species markers this cell draws are keyed in one row under its own box.

**I.** The same statistic against the reference stratified by its own modification
ratio (added 2026-09-27 for reviewer point R2-1): one row per stratum (ratio
0.1-0.3, 0.3-0.6, > 0.6), the same three species markers, and the same species
key under its own box.  The ordering moves more strongly here (rho = 0.48-0.90)
and the
top-ranked tool is no longer always m6Anet -- it differs from the primary one in
6 of the 9 species-by-stratum cells -- which is why the primary metrics are
reported against the full high-confidence reference.

**Sources.** Frozen tables in `../tables/source/` (`s8in_counts.tsv`,
`figS8_orca_counts.tsv`, `s8in_fpr_curlcake_scan.tsv`, `figS8_modratio_agreement.tsv`,
`figS8F_pairs.tsv`) plus `../tables/S8F_glori_agreement.tsv` for the Lin's CCC
series (the same values quoted in the manuscript and the reply letter),
`data/evaluation/tables/m6a_localization_curve.tsv` for
the window sweep, and `tables/tables/
TableS12_coverage_rank_stability.tsv` for the two stability cells (read, never
recomputed); rendered by `70_figs8_panels.py`, composed 1:1 by
`71_figs8_page.py`, accepted by `72_verify_figs8.py`; anchors in
`../tables/S8_anchors.tsv`.

**Differences to the originally submitted S8.** The figure is RNA004 only: the
Curlcake block of RNA002-only tools, the grey RNA002 control points, the
depth-matched nested subset, the per-unit means, the chemistry-comparison arrows
and the 13-tool colour stripe are gone.  The Jaccard matrices (R3-8 asked for
defined false-positive metrics instead of overlap) and the meta-gene band
(main-text Fig. 8D) are not part of the figure.  Since 2026-09-27 the page is a
three-by-three grid of nine lettered cells; the two stability cells (H coverage,
I reference ratio) are new and answer the reviewer's question about how the
coverage and the reference's own composition move the ranking.
