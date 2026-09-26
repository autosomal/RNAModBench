# REPLICATION_SPEC — Figure 8 / Figure S8 restyled to the published figure language

Single construction reference for `scripts/51_fig8_figure.py` and
`scripts/52_figS8_figure.py`.  Extracted 2026-09-21 from the published panel
sources and the approved revision implementations (file:line noted); the user
rejected the 2026-09-20 draft as "too ugly / not the published look", so every
value below is taken from the published code, not invented.

## Sources

| element | source |
|---|---|
| A dumbbell | `the original submission's figure code/Figure8/dorado_model_lollipop_plot.py:34-207` |
| B Jaccard heatmaps (S8E) | `…/result_RNA004/scripts/analysis/dorado_m6a_jaccard_analysis.py:120-148` |
| C dot plot + % callouts | `the original submission's figure code/Figure8/rna004_m6a_ngs_analysis.py:346-407` (dot version) and `:239-279` (bar version) |
| D guitar (published) | `the original submission's figure code/Figure8/dorado_m6a_guitar_wt_ivt.R:58-278` |
| D metagene (approved revision) | `src/harmonisation/scripts/47_fig7_rebuild.py:353-404` (`metagene_frame`, `panel_d`) + `common/figstyle.py:53-87` |
| S8 mod-ratio scatter (published) | `the original submission's figure code/Figure8/rna004_m6a_ngs_analysis.py` (mod-ratio vs GLORI section) |
| value-callout permissions | user decision 2026-09-21 (key numbers, >= 12 pt) |
| house rules | `common/figstyle.py`, `figure_text_audit/README.md`, user rules 2026-09-18/19 |

## Global constants

* Fonts: Arial everywhere; `pdf.fonttype=42`, `ps.fonttype=42`,
  `mathtext.fontset="custom"` + Arial (never DejaVu).
* Size scale on the final page (this is what makes the figure look "printed"):
  ticks 8.5, row labels 8.8 (A, C), axis labels 10.0, panel letters 15.0 bold,
  legends 8.5, **key value callouts 12.0 bold** (user decision).
* Page sizes (must equal the replaced files): Figure 8 = 595.276 x 740.0 pt;
  Figure S8 = 595.276 x 770.419 pt.  Save without `bbox_inches="tight"` so the
  pt sizes above are the printed sizes.
* Axis furniture: `top`/`right` spines off, left/bottom lw 1.1, ticks out,
  `FixedLocator` + `NullLocator` on log axes (no minors), no `1e+06` style
  labels, legends `frameon=False`.

### Palette

| meaning | value | note |
|---|---|---|
| family m6A DRACH | `#E83C1E` | published `MOD_COLORS` |
| family m6A (non-DRACH) | `#F5B264` | published |
| family m5C | `#2E86AB` | published |
| family pseU | `#A23B72` | published |
| family inosine (+ m6A) | `#F18F01` | published |
| family fallback | `#333333` | published |
| WT | `#4B81B8` | house convention (published guitar used the reverse, see deviations) |
| IVT | `#E8A76B` | house convention |
| C dot | `#1E888B` | published dot plot |
| metagene fill WT | `0.86` (grey) | approved Fig. 7 implementation |
| metagene fill IVT | `#F6DCC3` | approved Fig. 7 implementation |
| Jaccard heatmap | `YlOrRd`, vmin 0, vmax 1 | published seaborn heatmap |
| trend / connector grey | `0.72`–`0.75`; dumbbell rod `k` alpha 0.4 lw 2.0 | published |

## Panel specifications

### Figure 8A — Dorado model counts (published dumbbell)

* Rows: every Dorado model of the HeLa RNA004 WT/IVT pair, sorted by WT count
  descending (published sort); label = `mode@version_model` with `@v1`/`_otherMod`
  stripped and the inosine-m6A channels marked `#m6A` / `#inosine`.
* Colour = model family (5-colour palette), marker = condition
  (WT circle, IVT triangle), rod = black alpha 0.4 lw 2.0 between the two calls.
* Log x axis `Counts (log scale)` with `10^0 … 10^5` labels (published used
  `Counts (Log)`), xlim `floor … 2.2 x max`.
* Models with **zero IVT calls** sit at the axis floor as an open triangle and
  carry a bold `0` callout (12 pt) — replaces the published "invisible on log"
  behaviour and satisfies "0 must not be drawn as a fallen line".
* Legend: 2 sample handles (WT/IVT) + 5 family patches, upper-left corner (the
  only free region once the dumbbells are drawn).
* Deviations: no x grid (house rule forbids any grid; the published panel had a
  faint dashed x grid), counts are site-level deduplicated.

### Figure 8B — strict FPR, family x threshold grouped bars (**new panel**)

The reviewed claim needs an explicit rate; the published figure had no such
panel, so this is specified from scratch in the same visual language.

* Nested 2 x 1 axes inside panel B (shared x semantics, separate y ranges):
  * top: unmodified Curlcake IVT, x = Dorado `percent_modified` threshold
    {5, 10, 20, 50} %, 3 family bars per threshold
    (`#E83C1E` / `#F5B264` / `#F18F01`), bar = sup, open circle = hac;
  * bottom: unmodified HeLa IVT (real library), x = family
    (`DRACH`, `non-DRACH`, `inosine + m6A`), bar = family mean at the delivered
    `>= 90 %` cutoff, black whisker = model range (min–max).
* y = false positives per 10 kb, log scale; ranges 3e-5–2e3 (Curlcake) and
  6e-5–6 (HeLa); y ticks 1e-4 / 1e-2 / 1 / 100 / 10000 (and 1e-4 … 1 in the
  HeLa block); no `1e+0x` style labels.
* Zero bars are drawn at the floor and annotated; the four **key range callouts**
  use the 12 pt bold rule: `0-4`, `10-13` (50 %), `10-38`, `193-368` (10 %),
  `2-10` (HeLa DRACH), `13,028-46,963` (HeLa non-DRACH).
* Legend: 3 family patches + `hac (bar: sup)` marker entry, lower-left of the
  Curlcake block.
* Data: `tables/fig8_fpr_curlcake_scan.tsv`, `tables/fig8_fpr_hela_ivt.tsv`
  (`MOD_SLOT`: Curlcake `pseU_m6A` represents the non-DRACH m6A family, HeLa uses
  `m6A@v1`).

### Figure 8C — PPV vs. GLORI (published dot plot)

* Horizontal grey guide line `0 → value` (`color 0.72`, lw 1.6), dot
  `#1E888B`, s = 95, black edge lw 1.0 (published `s=200` on a wider canvas;
  scaled to this page).
* Right-hand `%.1f%%` callout per row at 12 pt bold (published used 14 pt on a
  10-inch canvas = same printed size).
* x axis `PPV vs. GLORI (2 bp)`, 0–1.22 with `0% … 100%` ticks, y rows = tools
  sorted ascending, row labels 8.8 pt, `y` ticks length 0.

### Figure 8D — filled metagene of the three Dorado m6A models (published guitar look)

* Three stacked axes (one per model: hac 5.0 m6A / hac 5.1 inosine + m6A /
  sup 5.0 m6A), each with the five-segment axis
  `1kb | 5'UTR | CDS | 3'UTR | 1kb` (equal widths, segment names only under the
  bottom axis), dotted black separators `ls=(0,(1,1.8))` lw 0.8.
* Per model: IVT fill `#F6DCC3` + IVT curve `#E8A76B` lw 1.5 (drawn first),
  WT fill `0.86` alpha 0.85 + WT curve `#4B81B8` lw 1.5 (drawn on top).
* Transcript schematic below each axis: baseline `0.35` lw 1.0 at
  `y=-0.055`, grey body `0.80` lw 6 at `y=-0.014`, black flank bars lw 3.2
  (published `48_fig7d_guitar.R:182-213` / `47_fig7_rebuild.py:353-373`).
* y = `Density` (middle axis only), per-axis autoscale; shared 2-entry legend
  (WT / IVT patches) under the bottom axis.
* Data: `tables/figS8G_metagene_density.tsv` (Ensembl GRCh38p14 release 112
  mRNA region model; never GENCODE).

### Figure S8 panels (same language)

| panel | spec |
|---|---|
| A counts (WT vs IVT, incl. ORCA) | paired bars per tool/modification, WT `#4B81B8` / IVT `#E8A76B`, log x, `int` tick labels, legend lower-right |
| B Curlcake chemistry | per-unit points joined by a thin line, RNA002 = `#8C8C8C`, RNA004 = `#2E86AB`, open markers for 0 calls, y = FP per 10 kb (log) |
| C HeLa IVT control | horizontal bars coloured by family (3-colour rule as B), log x, tick labels `0.0001 … 10` |
| D threshold detail | curves with `o`/`s` markers, sup solid / hac dashed, twin right axis for per-10^6 candidates (published secondary-axis style) |
| E Jaccard WT/IVT | two `imshow` matrices, `YlOrRd`, vmin 0 vmax 1, colorbar label `Jaccard`, no in-cell numbers (12 pt rule), short row/col labels |
| F mod-ratio vs GLORI | scatter + least-squares line per tool class, printed-class colours, r in the legend, x = GLORI ratio, y = tool ratio |
| G PPV–FPR trade-off | scatter, marker by family, log x (FP per 10 kb on HeLa IVT), y = PPV, legend with class counts |

## Deliberate deviations from the published graphics (approved or mandated)

1. **WT/IVT colours** — the published guitar R used WT = orange `#F5B264` and
   IVT = steel blue `#3778A0`, i.e. the reverse of every other panel and of the
   whole revision series; kept house semantics (WT blue / IVT orange) so the
   figure is internally consistent.
2. **No grids anywhere** — the published A panel had a faint dashed x grid;
   user rule 2026-09-18 forbids any grid.
3. **Site-level counts** — published counts were raw call rows; the revision
   counts distinct genomic positions (same semantics as the confusion tables).
   The published A-panel ordering therefore can differ slightly.
4. **Value callouts** restored only for the key numbers (>= 12 pt); all other
   numbers stay in the legends/tables.
5. **Panel B is new** (no published counterpart); it replaces the Jaccard
   matrices that the reviewer rejected as false-positive evidence.  Jaccard
   survives in S8E as model concordance only.
6. **HeLa IVT operating point** — the Dorado HeLa calls arrive pre-filtered at
   `percent_modified >= 90 %` and `valid_coverage >= 20`, so the HeLa block is
   drawn at that delivered cutoff and labelled as such (a full threshold sweep
   is only possible on Curlcake).

## Acceptance checklist (extends `53_verify_fig8_s8.py`)

* page size exact (595.276 x 740.0 / 595.276 x 770.419 pt, 1e-4 tolerance);
* only embedded Arial fonts (`pdffonts`);
* no `grid` call in the plotting scripts;
* minimum explicit font >= 7 pt, **value callouts >= 12 pt**;
* every plotted number reproduces the frozen tables at 1e-5 relative tolerance;
* panel letters A–D (main) / A–G (S8) match the legends file, no retired wording
  (`hit rate`, `near-perfect`).
