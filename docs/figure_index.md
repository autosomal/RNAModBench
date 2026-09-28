# Figure index

Every published figure, the script that draws it, and the tables it reads. Paths
are relative to the repository root. **Run order** is the order within a figure:
producers first, then panels, then the page assembly, then the layout gate.

All of a figure's scripts - producer, panels, page assembly, gate - live together in
that figure's own `figures/<figure>/src/`, so a directory is self-contained: its code,
then the frozen tables it reads. What is not there is the callset pipeline, which is
the numbered stages in `src/harmonisation/scripts/` writing `data/callsets/` and
`data/evaluation/tables/`; `docs/pipeline.md` covers those.

Environments: `py` = `envs/analysis/env-benchmark-revision.yml`,
`guitar` = `env-guitar_asm.yml` (R + Bioconductor Guitar), `motif` =
`env-motif_analysis.yml`, `enrich` = `env-enrich_r.yml`.

## Main figures

| figure | panels | scripts, in run order | env | reads |
|---|---|---|---|---|
| **Figure 1** | A–D overview, tool call counts, Curlcake density, DRACH counts | `figures/figure1/src/28_fig_tool_counts_replicates.py` → `figures/figure1/src/35_fig1C_curlcake_density.py` → `figures/figure1/src/36_fig1D_rrach_counts.py` → `figures/figure1/src/67_fig1_panels.py` → `figures/figure1/src/72_fig1_A_artwork.py` → `figures/figure1/src/69_fig1_page.py`; gate `figures/figure1/src/70_verify_fig1_page.py` | py | `figures/figure1/inputs/`, `figures/figure1/tables/` |
| **Figure 2** | A–F | `figures/figure2/src/01_fig2_panel_inputs.py` → `figures/figure2/src/02_fig2_go_enrichment.R` → `figures/figure2/src/03_fig2_figure.py`; gate `figures/figure2/src/04_verify_fig2.py` | py, enrich | `figures/figure2/tables/`, `data/evaluation/tables/` |
| **Figure 3** | A MDS, B modification-ratio agreement, C–D metagene | `figures/figure3/src/56_fig3a_tool_similarity_mds.py` → `figures/figure3/src/57_fig3b_modratio_wt_treatment.py` → `figures/figure3/src/58_fig3cd_guitar.R` → `figures/figure3/src/61_fig3_legend_band.py` → `figures/figure3/src/59_fig3_assembled.py`; gate `figures/figure3/src/60_verify_fig3.py` | py, guitar | `figures/figure3/tables/`, `figures/figure3/bed/` |
| **Figure 4** | A–D motif / KL | `figures/figure4/src/fig4_kl.py` → `figures/figure4/src/fig4_figures.py` | motif | `figures/figure4/tables/`, `figures/figure4/analysis/fig4_kmer_counts.tsv.gz` |
| **Figure 5** | A metric ranks, B ratio–hit-rate, C–D window sweep | `figures/figure5/src/35_fig5a_metric_ranks.py` → `figures/figure5/src/36_fig5b_modratio_hitrate.py` → `figures/figure5/src/37_fig5cd_window_sweep.py` → `figures/figure5/src/39_fig5_assembled.py` | py | `figures/figure5/tables/`, `data/evaluation/tables/m6a_localization_curve.tsv` |
| **Figure 6** | A–F tool combinations | `figures/figure6/src/40_fig6_combination.py` → `figures/figure6/src/65_fig6_combination_page.py`; gate `figures/figure6/src/66_verify_fig6_page.py` | py | `figures/figure6/tables/` |
| **Figure 7** | A–G known-site recovery, permutation null, metagene | `figures/figure7/src/47_fig7_rebuild.py` → `figures/figure7/src/48_fig7d_guitar.R` → `figures/figure7/src/43_fig7_metagene_density.py` → page is LaTeX: `pdflatex figures/figure7/src/assemble_fig7.tex`; gate `figures/figure7/src/46_fig7_verify.py` | py, guitar, TeX | `figures/figure7/tables/`, `analysis/nonm6a_false_positives/evidence/` |
| **Figure 8** | A–E RNA004 counts, FPR, agreement, trade-off, metagene | `figures/figure8/src/50_fig8_tables.py` → `figures/figure8/src/57_fig8a_counts.py` → `figures/figure8/src/59_fig8b_fpr.py` → `figures/figure8/src/58_fig8c_ppv.py` → `figures/figure8/src/67_fig8d_tradeoff.py` → `figures/figure8/src/69_fig8e_guitar.R` → `figures/figure8/src/68_assemble_a4.py`; gate `figures/figure8/src/53_verify_fig8.py` | py, guitar | `figures/figure8/tables/`, `data/callsets/RNA004/` |

## Supplementary figures

| figure | panels | scripts, in run order | env | reads |
|---|---|---|---|---|
| **Figure S1** | per-species metagene traces | `src/harmonisation/scripts/21b_export_guitar_bed.py` → `figures/figureS1/src/23d_figS1_guitar.R`; QC `src/harmonisation/scripts/23d_check_clipping.py` | guitar | Ensembl GTF via `$RNAMODBENCH_LOCAL` (regenerate `figures/figure3/bed/`) |
| **Figure S2** | A–E motif composition, bias control | `figures/figureS2/src/01_figS2_motif_inputs.py` → `figures/figureS2/src/08_figS2v4_stats.py` → `figures/figureS2/src/09_figS2v4_panels.py` → `figures/figureS2/src/10_figS2v4_assemble.py` | motif | `figures/figureS2/analysis/figS2_kl_full.tsv`, `figS2_pwm_per_rep.npz`, `analysis/motif_bias_control/` |
| **Figure S3** | A overlap, B PPV vs GLORI, C ratio agreement | `figures/figureS3/src/make_s3a_venn.py` → `figures/figureS3/src/make_s3b_ppv.py` → `figures/figureS3/src/make_s3c_reuse_other.py` → `figures/figureS3/src/make_figs3_page.py`; gate `figures/figureS3/src/verify_s3.py` | py | `figures/figureS3/tables/`, `analysis/mod_ratio_replicates/tables/` |
| **Figure S4** | A counts, B–C window sweeps, D purified vs WT | `figures/figureS4/src/41_figS4_tables.py` → `figures/figureS4/src/42_figS4_figure.py`; gate `figures/figureS4/src/49_verify_figS4_S10.py` | py | `figures/figureS4/tables/`, `data/evaluation/tables/` |
| **Figure S5** | A–C purified-site definition | `figures/figureS5/src/48_figS10_validation.py`; gate `figures/figureS4/src/49_verify_figS4_S10.py` | py | `figures/figureS4/tables/figS4_validation_groups.tsv` |
| **Figure S6** | A–E combination enumerations beyond the five selected tools | `figures/figureS6/src/61_figS5_sitequality.py` → `figures/figureS6/src/62_figS5_figure.py`; gate `figures/figureS6/src/63_verify_figS5.py` | py | `figures/figure6/tables/` (shared with Figure 6) |
| **Figure S7** | A–D non-m6A per independent unit | `figures/figureS7/src/53_figS6_tables.py` → `figures/figureS7/src/54_figS6_figure.py`; gate `figures/figureS7/src/55_verify_figS6.py` | py | `figures/figureS7/tables/`, `analysis/nonm6a_false_positives/` |
| **Figure S8** | A–F ncRNA metagene and background | `figures/figureS8/src/61_figS7_tables.py` → `figures/figureS8/src/62_figS7_guitar.R` → `figures/figureS8/src/63_figS7_page.R`; gate `figures/figureS8/src/64_verify_figS7.py` | py, guitar | `figures/figureS8/tables/`, `analysis/nonm6a_false_positives/` |
| **Figure S9** | A–F RNA004 model calling | `figures/figureS9/src/70_figs8_panels.py` → `figures/figureS9/src/71_figs8_page.py`; gate `figures/figureS9/src/72_verify_figs8.py` | py | `figures/figureS9/tables/source/` (mapping in its `SOURCE.md`) |
| **Figure S10** | Dorado other-modification models | `figures/figureS10/src/60_figS9_tables.py` → `figures/figureS10/src/23e_figS9_guitar.R`; gate `figures/figureS10/src/61_verify_figS9.py` | py, guitar | `data/callsets/RNA004/`, Ensembl annotation |

The page assembly of the combined `Supplementary_Figures.pdf` and of
`Supplementary_Tables.pdf` is in `tables/`: `tables/build_supp_tables.py` (S1-S9),
`tables/build_S10.py`, `tables/build_S11.py`, `tables/build_S12.py`, then `tables/order_supp_tables.py`,
which puts every table in the order the manuscript uses, and
`tables/render_supp_tables.py` (`tables/merge_supplementary_figures.py` assembles the figure
album).

## Naming note

Folders are named by the **published** figure number. Internal working folders
`figS5`–`figS10` were renumbered when the supplementary material was ordered for
publication, so inside a few file names and comments the older numbers still
appear. The mapping is:

| published | was internally |
|---|---|
| Figure S5 | `figures/figureS5` |
| Figure S6 | `figures/figureS6` |
| Figure S7 | `figures/figureS7` |
| Figure S8 | `figures/figureS8` |
| Figure S9 | `figures/figureS9` |
| Figure S10 | `figures/figureS10` |

## Two panels are not code-generated

* **Figure 1A** is a hand-drawn schematic, in the original submission as well.
  `figures/figure1/src/72_fig1_A_artwork.py`
  only places the artwork into the page, and the artwork file itself is not
  redistributed. The
  code-only alternative is `figures/figure1/src/71_fig1_A_redraw.py` (a pure
  re-draw, `Figure1_A_redrawn.pdf`); `figures/figure1/src/69_fig1_page.py` documents the switch.
* **The graphical abstract** is likewise hand-made and is not part of this deposit.

## Figure layout conventions

Figures are printed at 1:1 size as vector PDF (plus a 300 dpi PNG for the
publisher, not deposited here). Font sizes, page budgets and the collision checks
are enforced by the shared library: `src/harmonisation/common/figstyle.py` (typography),
`src/harmonisation/common/pagelayout.py` (A4 portrait/landscape budget, tick-label and grid checks),
`src/harmonisation/common/panelpage.py` (panel composition into a page). Each figure's `verify_*.py` gate
calls into them and exits non-zero on a violation.
