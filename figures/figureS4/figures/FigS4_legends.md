# Figure S4 (revision) — legends

**Figure S4.** Purified-site counts and the window sensitivity of the m6A
benchmark under the revised definition. The page keeps the row layout of the
submitted Figure S4 (four lettered panels, A–D): every row splits into the three
species columns (Arabidopsis | mouse | human) and the 13 tool configurations
share one palette, explained in a **centred colour-key stripe between panel A
and the B/C canvas** (two rows of seven entries, the middle of the page). The
keys that belong to a single panel sit on their own line inside it: the species /
mouse-study / mean key at the top of panel A and the line-type key (which also
names the two mouse study styles and the two panel-D entries) at the top of the
B/C canvas. The window sweep is split into **two stacked facets per species
column** — **B** = PPV vs. GLORI and **C** = Exact-nucleotide fraction — which
share one x axis and each carry a panel letter of their own in the page margin;
the window ticks and the axis name are printed once, under C. The page is A4
**portrait** (8.27 × 11.70 in, the page of the submitted figure): every species
column carries its own y tick labels, the panel letters sit in the page margin
(left of the y axis labels), and the species sub-panels are close to square —
2.16 × 1.83 in for each of B and C and 2.16 × 1.90 in for D. The 13 tool names
on the x axis of panel A are tilted by 45 degrees (ha=right), exactly as in the
submitted figure, and the layout gate measures the printed ink of the rotated
names (a 1.5 pt gap between two neighbours at 8 pt). The four pieces (A, the key
stripe, the B/C canvas, D) are drawn one at a time on their own canvases at
their final print size and then assembled onto the page
(`common/panelpage.py`, pypdf translation, vector preserved); every canvas has
to pass
`common/pagelayout.py`, which re-measures all text boxes at render time and
aborts on any text-text overlap, any text entering another panel's drawing area,
**any text sitting on an axis line or tick mark and any text closer than 3 pt to
a frame line** (the sub-row titles of panel A are indented by 0.07 in so that
they clear the top of their left spine, and the column gap is 0.50 in so that the
neighbouring column's tick labels clear the previous frame). All panels use the
clean per-replicate call sets
(`harmonisation/callsets`, RNA002 chemistry), the per-sample candidate universes
(annotated exons, coverage ≥ 10, reference-base compatible) and the unified
GLORI reference (Arabidopsis 80,624 / mouse liftover 41,961 / HeLa 112,451
sites). A reported site counts as a true positive when it lies within the
matching window *w* of a GLORI site of that sample's own universe, so the
quantity labelled **PPV vs. GLORI** is the fraction of a tool's reported sites
that are GLORI-concordant; GLORI enters as an orthogonal **site-level partial**
reference rather than as a complete ground truth. Every statistic is computed
per independent sequencing unit; the two
mouse datasets (SRP166020 = study A, SRP357195 = study B) are independent
cross-study datasets and are **never averaged** — they are drawn as separate
markers (filled = study A, open = study B) and separate mean lines (solid =
study A, dashed = study B).

- **(A)** Calls per tool inside the common testable universe of each
  WT/deficient pair (Arabidopsis WT × *fip37* KD, n = 3 pairs; mouse one pair
  per study; HeLa WT × IVT, n = 3 pairs), shown as three stacked sub-rows:
  WT-common counts (log scale), purified counts (log scale) and the purified/WT
  ratio. "Purified" = called in the WT sample and absent from the matched
  deficient call set, both restricted to the pair's common universe; the
  deficient samples are the *fip37* knockdown, the Mettl3 knockout (two
  studies) and the matched HeLa IVT negative control (an unmodified copy of the
  same transcriptome background rather than a knockdown of the same RNA). Dots
  are the individual pairs of a column (three for Arabidopsis and HeLa) and the
  black tick is their mean. The mouse column deliberately gets **no** mean tick:
  the two cross-study datasets are drawn as two dots joined by a thin line, so
  nothing is ever averaged. Most tools keep nearly all of their WT-common calls
  (CHEUI, ELIGOS2_diff and DRUMMER > 98 %), whereas MINES keeps only 36 %,
  Nanom6A 49 % and DENA 57 % (main text).
- **(B)** The window dependence of the benchmark, upper facet: the PPV vs. GLORI
  as a function of the matching window *w* (0–50 bp), per tool, in the 13 tool
  colours used by every other panel. Two layers per tool: thin lines are the
  individual sequencing units and the thick line is their group mean. For mouse
  the two cross-study datasets each contribute one sequencing unit, so one thick
  line is drawn per study (study A solid, study B dashed) and the two are never
  merged. The apparent gain of a wide window is dominated by localisation
  tolerance rather than detection (compare Fig. 5C, D; facet **C** below shows
  the localisation side of the same trade-off); *w* = 2 bp is the primary
  definition used throughout the revised manuscript, at which the tool-averaged
  PPV is 41.1 % (Arabidopsis), 27.7 % (HeLa) and 26.3 % pooled over the two mouse
  studies (24.9 % study A, drawn solid; the main-text mouse values are those of
  study B, drawn dashed). Expanding the window from 0 to ±50 bp raises the
  tool-averaged PPV from 28.9 % to 51.4 % (Arabidopsis), 24.3 % to 36.7 % (HeLa)
  and 22.1 % to 37.3 % (mouse, two studies pooled; the main-text mouse endpoints,
  23.5 % → 37.7 %, are again those of study B).
- **(C)** The same window sweep as panel B, with the same axes, the same colours
  and the same two layers per tool, but plotting the **Exact-nucleotide
  fraction**: among the calls that match a GLORI site at *w*, the fraction that
  falls on the exact reference nucleotide. Its decline is what makes the
  wide-window PPV an effect of tolerance rather than of detection (R1-3 / E2 /
  R3-3) — 1.00 at 0 bp → ≈0.52 at ±50 bp for all three species, the equivalent
  main-text panel being Fig. 5D. Panels B and C share one x axis, so the window
  ticks and the axis name are printed once, under C, and the trade-off is read
  by going down each species column: PPV up, exact fraction down. For mouse, as
  in B, one line per study and facet (study A solid, study B dashed), never
  merged.
- **(D)** PPV (2 bp) of the full WT-common call set versus the purified subset
  (this panel was lettered C before the window sweep was split into B and C on
  2026-09-21); each thin line connects the same tool in the same pair
  (colour = tool), the connectors are ordered by their shift (purified − WT)
  inside the species column so they do not cross, and the black points with
  error bars are the mean ± SD over all tool–pair combinations of a species. An
  alternative box-plot rendering of the same quantity (two boxes per species
  with the raw tool–pair points and a black mean diamond) is kept as a variant
  (`figures/variants/FigureS4_rev_Dbox.pdf`); the dumbbell version is the
  delivered one because it keeps the pairing. **Descriptive only**: conditioning on absence from a deficient
  call set is circular for accuracy assessment (absence can reflect coverage or
  expression rather than the absence of modification, R3-7), so this panel is
  not used as evidence of improved accuracy: the subset is selected on each
  tool's own calls rather than on the reference, so its higher PPV is expected
  rather than a measure of higher accuracy. Averaged over the 104 tool–pair
  combinations, PPV rises only modestly, 0.326 → 0.348, and 80 of the 104
  combinations shift upwards.

**Panel identity.** The main text cites Fig. S4A (counts after purification)
and Fig. S4D (precision on the purified subset); the width of the matching window
is examined in Fig. S4B (PPV) and Fig. S4C (exact-nucleotide fraction), the two
lettered facets of the window sweep. Before 2026-09-21 the window sweep was a
single panel B and the purified comparison was panel C; the re-lettering moved
the purified comparison to **D**, and the manuscript, this legend and the
response letter were updated together. The *independent*
validation of the purified sites requested by R3-7 (GLORI overlap, DRACH
context, stoichiometry-semantic scores of the three site groups) is **not** part
of this figure any more — it is the separate **Figure S10**
(`figures/figureS5/`).

## Caveats (also stated in the manuscript)

1. **Circularity (R3-7):** the purified definition conditions on absence from a
   deficient call set. Panels A and D are descriptive and are not used to claim
   improved accuracy.
2. **HeLa IVT pairing:** IVT libraries are true negative controls (unmodified
   RNA), so "purified" for HeLa means WT-called sites not recalled from an
   unmodified copy of the same transcriptome — a stronger, not a weaker,
   criterion than the KO/KD case.
3. **Mouse independence:** SRP166020 and SRP357195 are separate studies, never
   collapsed into n = 2; they are drawn separately in every panel.
