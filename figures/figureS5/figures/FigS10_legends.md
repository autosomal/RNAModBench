# Figure S10 (new, revision) — legends

**Figure S10.** Independent characterisation of the purified-site definition
(R3-7). Every WT/deficient pair is split into three site groups inside the
pair's common testable universe: **def-only** (called only in the deficient
sample), **shared** (called in both samples) and **purified** (called in WT
only). Deficient samples: *fip37* knockdown (Arabidopsis), Mettl3 knockout
(mouse, two studies) and the matched IVT negative control (HeLa). All data come
from the per-replicate `harmonisation` call sets, the per-sample universes
(annotated exons, coverage ≥ 10, reference-base compatible) and the unified
GLORI reference (2-bp matching window). Dots are the individual pairs
(Arabidopsis and HeLa n = 3; mouse one pair per study, filled = study A
SRP166020, open = study B SRP357195, **never averaged**) and the black tick is
their mean; error bars are the SD across pairs for A and B (the mean IQR of the
per-pair score distributions for C). The y ticks are labelled on every column.
Fonts use the **enlarged scale of this figure** (tick 12 pt, axis 14 pt, species
title 17 pt, panel letter 24 pt, key 10.5 pt at the final print size): S10 is
printed on the same A4 page as Figure S4, but its panels are taller, so a scale
of its own keeps the type proportional to the drawing (Figure S4 keeps its own,
smaller scale). The page uses the same layout system as Figure S4: **A4 portrait**
(8.27 × 11.70 in), its shared key drawn on its own line inside the panel-A canvas
(no page-level legend stripe), panel letters in the page margin, one canvas per
panel (2.13 × 2.84 in sub-panels, panel C split into two halves of three tools),
assembled with `common/panelpage.py`, every canvas gated by
`common/pagelayout.py` at render time (no text overlap, no text entering another
panel, no text on or within 3 pt of an axis line).

- **(A)** Overlap with the GLORI reference (2 bp). Purified sites are strongly
  enriched relative to both control groups: 0.348 (purified) vs. 0.158 (shared)
  and 0.060 (def-only) as means over the 104 tool–pair combinations — a ~6-fold
  enrichment of GLORI-supported positions over the deficient-only calls, and
  ~2-fold over the calls that persist in the deficient sample. Per species
  (purified / shared / def-only): Arabidopsis 0.437 / 0.180 / 0.051, mouse
  0.278 / 0.219 / 0.076, HeLa 0.306 / 0.105 / 0.059.
- **(B)** DRACH-motif fraction of the same groups, against the DRACH fraction
  of each sample's own universe (dashed line, per sample; ≈ 0.067 Arabidopsis,
  0.077 mouse, 0.071 HeLa — i.e. the universes are motif-agnostic, ~7-fold below
  the site groups). Purified sites move furthest away from that background:
  0.516 (purified) vs. 0.483 (shared) and 0.440 (def-only) as means over the
  104 tool–pair combinations. Per species: Arabidopsis 0.476 / 0.452 / 0.416,
  HeLa 0.524 / 0.437 / 0.416; mouse 0.564 / 0.625 / 0.534 is the only column in
  which shared sites are the most motif-enriched group.
- **(C)** Stoichiometry-semantic score (mod-ratio: CHEUI_m6A, DENA, MINES,
  Nanom6A; probability: m6Anet, NanoSPA_m6A) of purified vs. shared sites. The
  direction is **tool-dependent** and therefore reported descriptively: purified
  minus shared is +0.004 (CHEUI), +0.038 (Nanom6A), −0.047 (m6Anet), −0.053
  (DENA), −0.125 (MINES) and −0.008 (NanoSPA_m6A). Purified sites therefore do
  **not** carry a systematically higher modification ratio or probability; the
  independent support for the group comes from the GLORI overlap (A) and the
  DRACH context (B), not from this panel. No claim is made for tools whose
  native score is a *P*-value-like statistic.

**Relation to Figure S4.** These three panels were the "E–G" panels of the
first draft of the revised Figure S4. They answer a different reviewer point
(R3-7, circularity of the purified definition) from the panels of Figure S4
(counts after purification, window sweep, PPV before/after purification), and
Figure S4 has to keep the row layout of the submitted figure, so they were split
into this separate, appended figure (S10). Panel letters restart at A.

## Caveats (also stated in the manuscript)

1. **Circularity (R3-7):** the groups are defined by presence/absence in the
   deficient call set. This figure characterises the resulting sites; it does
   not by itself prove their modification status, and GLORI is only a partial
   (site-level) reference.
2. **Group sizes differ by orders of magnitude** (purified ≫ def-only ≫ shared
   for most tools), so panel C compares score distributions of very different
   sizes; it is descriptive and no statistical test is attached.
3. **HeLa IVT pairing:** the IVT library is an unmodified copy of the same
   transcriptome, so "purified" for HeLa means "not recalled from an unmodified
   copy" — a stronger criterion than the KO/KD case, not a weaker one.
4. **Mouse independence:** the two mouse studies are independent datasets and
   are never collapsed into n = 2.
