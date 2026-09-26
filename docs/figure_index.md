# Figure index

Every delivered figure, the script that draws it, and the tables it reads. Paths
are relative to the repository root. **Run order** is the order within a figure:
producers first, then panels, then the page assembly, then the layout gate.

Environments: `py` = `envs/analysis/env-benchmark-revision.yml`,
`guitar` = `env-guitar_asm.yml` (R + Bioconductor Guitar), `motif` =
`env-motif_analysis.yml`, `enrich` = `env-enrich_r.yml`.

## Main figures

| figure | panels | scripts, in run order | env | reads |
|---|---|---|---|---|
| **Figure 1** | A–D overview, tool call counts, Curlcake density, DRACH counts | `src/sites_v2/scripts/28_fig_tool_counts_replicates.py` → `35_fig1C_curlcake_density.py` → `36_fig1D_rrach_counts.py` → `67_fig1_panels.py` → `72_fig1_A_from_ai.py` → `69_fig1_page.py`; gate `70_verify_fig1_page.py` | py | `figures/figure1/inputs/`, `figures/figure1/tables/` |
| **Figure 2** | A–F | `figures/figure2/src/01_fig2_panel_inputs.py` → `02_fig2_go_enrichment.R` → `03_fig2_figure.py`; gate `04_verify_fig2.py` | py, enrich | `figures/figure2/tables/`, `data/evaluation/tables/` |
| **Figure 3** | A MDS, B modification-ratio agreement, C–D metagene | `src/sites_v2/scripts/56_fig3a_tool_similarity_mds.py` → `57_fig3b_modratio_wt_treatment.py` → `58_fig3cd_guitar.R` → `61_fig3_legend_band.py` → `59_fig3_assembled.py`; gate `60_verify_fig3.py` | py, guitar | `figures/figure3/tables/`, `figures/figure3/bed/` |
| **Figure 4** | A–D motif / KL | `figures/figure4/src/fig4_kl.py` → `fig4_figures.py` | motif | `figures/figure4/tables/`, `figures/figure4/analysis/fig4_kmer_counts.tsv.gz` |
| **Figure 5** | A metric ranks, B ratio–hit-rate, C–D window sweep | `src/sites_v2/scripts/35_fig5a_metric_ranks.py` → `36_fig5b_modratio_hitrate.py` → `37_fig5cd_window_sweep.py` → `39_fig5_assembled.py` | py | `figures/figure5/tables/`, `data/evaluation/tables/m6a_localization_curve.tsv` |
| **Figure 6** | A–F tool combinations | `src/sites_v2/scripts/40_fig6_combination.py` → `65_fig6_combination_page.py`; gate `66_verify_fig6_page.py` | py | `figures/figure6/tables/` |
| **Figure 7** | A–G known-site recovery, permutation null, metagene | `src/sites_v2/scripts/47_fig7_rebuild.py` → `48_fig7d_guitar.R` → `43_fig7_metagene_density.py` → page is LaTeX: `pdflatex figures/figure7/src/assemble_fig7.tex`; gate `46_fig7_verify.py` | py, guitar, TeX | `figures/figure7/tables/`, `analysis/nonm6a_false_positives/evidence/` |
| **Figure 8** | A–E RNA004 counts, FPR, agreement, trade-off, metagene | `figures/figure8/src/50_fig8_tables.py` → `57_fig8a_counts.py` → `59_fig8b_fpr.py` → `58_fig8c_ppv.py` → `67_fig8d_tradeoff.py` → `69_fig8e_guitar.R` → `68_assemble_a4.py`; gate `53_verify_fig8.py` | py, guitar | `figures/figure8/tables/`, `data/sites_clean/RNA004/` |

## Supplementary figures

| figure | panels | scripts, in run order | env | reads |
|---|---|---|---|---|
| **Figure S1** | per-species metagene traces | `src/sites_v2/scripts/21b_export_guitar_bed.py` → `23d_figS1_guitar.R`; QC `23d_check_clipping.py` | guitar | Ensembl GTF via `$RNAMODBENCH_LOCAL` (regenerate `figures/figure3/bed/`) |
| **Figure S2** | A–E motif composition, bias control | `figures/figureS2/src/01_figS2_motif_inputs.py` → `08_figS2v4_stats.py` → `09_figS2v4_panels.py` → `10_figS2v4_assemble.py` | motif | `figures/figureS2/analysis/figS2_kl_full.tsv`, `figS2_pwm_per_rep.npz`, `analysis/motif_bias_control/` |
| **Figure S3** | A overlap, B PPV vs GLORI, C ratio agreement | `figures/figureS3/src/make_s3a_venn.py` → `make_s3b_ppv.py` → `make_s3c_reuse_other.py` → `make_figs3_page.py`; gate `verify_s3.py` | py | `figures/figureS3/tables/`, `analysis/mod_ratio_replicates/tables/` |
| **Figure S4** | A counts, B–C window sweeps, D purified vs WT | `src/sites_v2/scripts/41_figS4_tables.py` → `42_figS4_figure.py`; gate `49_verify_figS4_S10.py` | py | `figures/figureS4/tables/`, `data/evaluation/tables/` |
| **Figure S5** | A–C purified-site definition | `src/sites_v2/scripts/48_figS10_validation.py`; gate `49_verify_figS4_S10.py` | py | `figures/figureS4/tables/figS4_validation_groups.tsv` |
| **Figure S6** | A–E combination enumerations beyond the five selected tools | `src/sites_v2/scripts/61_figS5_sitequality.py` → `62_figS5_figure.py`; gate `63_verify_figS5.py` | py | `figures/figure6/tables/` (shared with Figure 6) |
| **Figure S7** | A–D non-m6A per independent unit | `src/sites_v2/scripts/53_figS6_tables.py` → `54_figS6_figure.py`; gate `55_verify_figS6.py` | py | `figures/figureS7/tables/`, `analysis/nonm6a_false_positives/` |
| **Figure S8** | A–F ncRNA metagene and background | `src/sites_v2/scripts/61_figS7_tables.py` → `62_figS7_guitar.R` → `63_figS7_page.R`; gate `64_verify_figS7.py` | py, guitar | `figures/figureS8/tables/`, `analysis/nonm6a_false_positives/` |
| **Figure S9** | A–F RNA004 model calling | `figures/figureS9/src/70_figs8_panels.py` → `71_figs8_page.py`; gate `72_verify_figs8.py` | py | `figures/figureS9/tables/source/` (mapping in its `SOURCE.md`) |
| **Figure S10** | Dorado other-modification models | `figures/figureS10/src/60_figS9_tables.py` → `src/sites_v2/scripts/23e_figS9_guitar.R`; gate `61_verify_figS9.py` | py, guitar | `data/sites_clean/RNA004/`, Ensembl annotation |

The page assembly of the combined `Supplementary_Figures.pdf` and of
`Supplementary_Tables.pdf` is in `tables/` (`merge_supplementary_figures.py`,
`render_supp_tables.py`, `build_supp_tables.py`).

## Naming note

Folders are named by the **delivered** figure number. Internal working folders
`figS5`–`figS10` were renumbered when the supplementary material was ordered for
publication, so inside a few file names and comments the older numbers still
appear. The mapping is:

| delivered | was internally |
|---|---|
| Figure S5 | `figS10_revision` |
| Figure S6 | `figS5_revision` |
| Figure S7 | `figS6_revision` |
| Figure S8 | `figS7_revision` |
| Figure S9 | `figS8_revision` |
| Figure S10 | `figS9_revision` |

## Two panels are not code-generated

* **Figure 1A** is a hand-drawn schematic, in the original submission as well.
  `72_fig1_A_from_ai.py`
  only places the artwork into the page, and the artwork file itself is not
  redistributed. The
  code-only alternative is `src/sites_v2/scripts/71_fig1_A_redraw.py` (a pure
  re-draw, `Figure1_A_redrawn.pdf`); `69_fig1_page.py` documents the switch.
* **The graphical abstract** is likewise hand-made and is not part of this deposit.

## Figure layout conventions

Figures are printed at 1:1 size as vector PDF (plus a 300 dpi PNG for the
publisher, not deposited here). Font sizes, page budgets and the collision checks
are enforced by the shared library: `src/sites_v2/common/figstyle.py` (typography),
`pagelayout.py` (A4 portrait/landscape budget, tick-label and grid checks),
`panelpage.py` (panel composition into a page). Each figure's `verify_*.py` gate
calls into them and exits non-zero on a violation.
