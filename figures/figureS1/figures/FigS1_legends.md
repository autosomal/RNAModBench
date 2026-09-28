# Figure S1 (revised) — legend and caveats

File: `figures/FigureS1_rev.pdf` (+ 300 dpi `FigureS1_rev.png`), page A4
(595.276 x 841.89 pt), Arial embedded (cairo_pdf, `pdffonts`-verifiable), no
gridlines, minimum font 7 pt. This file replaces the published
`$RNAMODBENCH_LOCAL/submission/sup/sup1.pdf` **as a new file; the original is
untouched.**

## Suggested caption (English)

> **Figure S1. Metagene distribution of m6A sites across mRNA regions, drawn
> replicate-aware.** For each of the ten tool configurations of the published
> figure (CHEUI_m6A, DRUMMER, ELIGOS_diff, ELIGOS_solo, EpiNano_Error, MINES,
> Nanocompore, NanoSPA_m6A, xPore, Yanocomp), the metagene shows the normalized
> position of detected sites along mRNA transcripts (1 kb upstream flank,
> 5'UTR, CDS, 3'UTR, 1 kb downstream flank; transcript direction), computed
> with the Bioconductor *Guitar* package on the majority consensus of each
> library's independent sequencing units (thick line, shaded area), with one
> thin line per biological replicate (Arabidopsis, Human; *n* = 3) or per
> independent study (Mouse; *n* = 2, never pooled). Consensus curves built
> from fewer than 100 sites are omitted and only the replicate/study lines are
> shown (gated list: `tables/figS1_gated_consensus.tsv`) — for the
> differential tools (ELIGOS2_diff, DRUMMER) the IVT side is nearly empty by
> construction, because significant sites are almost exclusively "up in WT".
> Wild type in blue, the modification-deficient condition (fip37 KD / Mettl3
> KO / IVT) in orange. Density is the Gaussian-smoothed site occupancy in
> arbitrary units; y-axes are scaled per panel. Unlike the published panel,
> which was drawn from single-replicate (Arabidopsis, Mouse) or
> undocumented-union (Human) inputs, every curve here is traceable to explicit
> per-unit call sets (`tables/figS1_panel_inputs.tsv`). The published tool
> label "ELIGOS_diff" / "ELIGOS_solo" / "Yanocomp" corresponds to
> ELIGOS2_diff / ELIGOS2_solo / yanocomp in the revised pipeline.

## Data provenance and rules (for the response letter / internal record)

1. **Inputs.** `$RNAMODBENCH_LOCAL/guitar_metagene/bed/RNA002/`
   written by `harmonisation/scripts/21b_export_guitar_bed.py` from the per-unit
   call sets under `harmonisation/callsets/`:
   - `majority/<Condition>/m6A/<Tool>.bed` — sites supported by *strictly more
     than half* of the group's independent units (`common.consensus.quorum`);
   - `rep_<tag>/<Condition>/m6A/<Tool>.bed` — one curve per independent unit.
   Intervals are padded 1 bp on each side exactly as the original
   `Guitar_*.r` scripts did; chromosome labels match the Ensembl GTF handed to
   Guitar (TAIR10.61 / GRCm39.114 / GRCh38.112 — Ensembl-only, per house rule).
2. **Density kernel.** The Bioconductor *Guitar* package itself
   (`samplePoints -> normalize -> .generateDensity_CI`, `enableCI = FALSE`,
   `CI_ResamplingTime = 20`; verified identical to `= 1000` to the last digit,
   `logs/figS1_probe.log`). Component rectangles/labels reproduce
   `Guitar:::.RNAPlotStructure` with all five region labels on one line —
   `1kb / 5'UTR / CDS / 3'UTR / 1kb` (both flanks are labelled `1kb`) — and
   the flank labels clamped inside `[0.05, 0.95]` of the panel so they are
   never clipped at the panel edge; label size pinned to 8 pt.
3. **Mouse caveat (house rule: two studies are never averaged).** The
   consensus of `Mouse_WT`/`Mouse_KO` is a *cross-study concordance filter*
   (a site must be called in both studies). Where a tool was only run on one
   study, the consensus equals that single study; the missing study's thin
   curve is correspondingly absent:
   - Mouse_KO single-study tools: DRUMMER, ELIGOS2_diff, EpiNano_Error,
     Nanocompore, xPore, yanocomp (study A = mESCs_Mettl3_KO only).
   - Mouse_WT: all ten tools have both studies since the 2026-09-20 data
     refresh (Nanocompore study A = mESCs_Mettl3_WT was added that day; its
     consensus dropped from 3,023 sites = study B alone to 214 sites = the
     two-study intersection, see `tables/figS1_nanocompore_refresh.tsv`).
   Study A = mESCs_Mettl3_WT / mESCs_Mettl3_KO, study B = mES_WT / mES_KO
   (`sample_registry.csv`).
4. **Low-input curves, kept and flagged** (`figS1_panel_inputs.tsv`, column
   `note`; never silently dropped):
   - Arabidopsis KD DRUMMER consensus: 6 sites;
   - Arabidopsis KD ELIGOS2_diff replicate 3: 8 sites;
   - **Arabidopsis KD ELIGOS2_diff has no majority consensus at all** (no site
     supported by >= 2 of 3 replicates; its KD panel therefore shows the three
     replicate curves without the orange consensus), and
   - Human IVT DRUMMER replicate 2 (8 sites) crashed inside
     `Guitar::samplePoints` ("subscript out of bounds"); the curve is skipped
     and the failure recorded.
5. **Panel counts are BED rows**, i.e. all sites of the call set before
   transcript assignment; the share that maps onto an annotated mRNA (what
   Guitar actually plots) is smaller and recorded per group x tool in
   `guitar_metagene_replicates/tables/merge_site_counts.tsv` and
   `metagene_density.tsv.gz` (38-47% of the m6A union rows).
6. **Geometry self-check.** `tables/figS1_geometry.tsv`: 3 blocks x 10 panels,
   5 columns; fonts 7/7/7.5/8/11 pt (min 7 pt); `CI_ResamplingTime = 20`;
   consensus rule string.

## What changed vs. the published sup1.pdf

| aspect | published sup1.pdf | FigureS1_rev |
|---|---|---|
| replicate structure | none (Arabidopsis rep3 only; Mouse mES_WT only; HeLa undocumented 3-replicate union) | majority consensus over independent units + one curve per unit |
| strand | forced `+` | site strand honoured upstream (`21b`, `strand_mode = aware`) |
| drawing | Illustrator-traced Guitar panels, 3.9-5 pt text, non-searchable | drawn by R/ggplot2 on Guitar densities, >= 7 pt Arial, embedded and searchable |
| naming | ELIGOS_diff / ELIGOS_solo / Yanocomp | ELIGOS2_diff / ELIGOS2_solo / yanocomp (labels keep the published spelling, mapping stated here) |
| tools | 10 | same 10 (per reviewer-facing comparability) |
| provenance | none | per-curve table `figS1_panel_inputs.tsv` |

## Reproduce

```bash
conda run -n guitar_asm --no-capture-output Rscript \
  figures/figureS1/src/23d_figS1_guitar.R          # all species
# flags: --species Mouse   --recompute   --merge majority   --min-sites 10   --rt 20
#        --mode panels|assemble|all        --only-panel Mouse_Nanocompore
```
Density is cached per species in `tables/figS1_density_<Species>.rds`; re-styling
does not recompute sampling.

## Two-stage rendering (2026-09-20, natural-size stitch, not A4)

1. **Standalone panels** — each species x tool subplot is drawn ONCE at
   3.4 x 2.4 in with 8-11 pt type and **no in-panel legend**, and saved to
   `figures/panels/FigureS1_rev_<Species>_<Tool>.{pdf,png}` (30 files).
   `--only-panel Mouse_Nanocompore` re-renders a single panel.
2. **Page assembly** — `--mode assemble` re-runs the same `draw_panel()` and
   the code stitches the 30 panels by the fixed pattern — one bold
   species-title row (top, horizontally centred) + a 2 x 5 tool grid per
   species, three species blocks — into a single page at the sum of the panel
   sizes (**~17.6 x 17.2 in; the author waived the A4 constraint for this
   supplement figure**). Because the panels are placed 1:1 (no scaling),
   curves cannot be stretched; the legend (thick solid = majority consensus,
   thin dashed = individual replicate/study; blue = WT, orange = fip37 KD /
   Mettl3 KO / IVT) is collected once per species block.

A rendering audit (`23d_check_clipping.py` -> `tables/figS1_clipping_check.tsv`)
asserts that no text or vector element of any panel or of the page extends
beyond the page box (numeric 0-1 x ticks are removed on purpose — the region
labels are the x axis).

## Consensus gating (< 100 sites)

A thick filled consensus curve is drawn **only when the majority consensus of
that condition keeps >= 100 sites** (`--gate 100`); below the threshold the
panel shows the per-replicate/study thin lines alone, because a smoothed curve
over a handful of sites would look like signal. Gated conditions
(`tables/figS1_gated_consensus.tsv`):

| species | tool | condition | consensus sites |
|---|---|---|---|
| Arabidopsis | DRUMMER | fip37 KD | 6 |
| Mouse | NanoSPA_m6A | Mettl3 KO | 18 |
| Human | DRUMMER | IVT | 11 |
| Human | ELIGOS2_diff | IVT | 24 |
| Human | Nanocompore | IVT | 55 |
| Human | NanoSPA_m6A | IVT | 94 |

Arabidopsis KD ELIGOS2_diff has no majority consensus at all (no site supported
by >= 2 of 3 replicates), so its KD side also shows replicate lines only.

## Why the differential tools' IVT side is almost empty (not a data error)

ELIGOS2_diff (and DRUMMER) compare two conditions and return sites that differ
between them; the benchmark attributes each side's calls by tool semantics
("ELIGOS2 take their `-t` side"). For an unmodified IVT library almost nothing
is "up in IVT", so the IVT-side call sets are tiny while the full comparison
tables are large (Human: combine tables 96k-166k rows vs IVT-side result files
28/83/73 rows). The published figure drew these panels from the same sparse
data — as an undocumented 3-replicate **union** (WT 2,530 / IVT 157 rows) —
while this revision uses the replicate-aware **majority** (WT 816 / IVT 24,
i.e. sites supported by >= 2 of 3 replicates). The thinner curves are the
honest replicate-aware view, not lost data; per-curve provenance is in
`tables/figS1_panel_inputs.tsv`.

BED inputs were refreshed the same day with `21b_export_guitar_bed.py`
(`--groups <6 S1 groups> --mods m6A`) after Nanocompore data were added;
before/after row counts for all 278 exported BEDs are in
`tables/figS1_nanocompore_refresh.tsv` (only the two Mouse_WT/Nanocompore
files changed).
