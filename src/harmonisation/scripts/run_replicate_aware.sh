#!/usr/bin/env bash
# Replicate-aware revision analyses (GUITAR metagene + reviewer figures).
#
#   bash scripts/run_replicate_aware.sh [RNA002|RNA004]
#
# Prerequisite: the harmonisation pipeline itself (scripts/run_all.sh) must have
# finished, because steps 24/26/27 read its evaluation tables and 28 reads
# callsets/ + manifest/.
set -euo pipefail
cd "$(dirname "$0")/.."
PLATFORM="${1:-RNA002}"

python scripts/20_build_region_model.py
python scripts/21_metagene_replicates.py --platform "$PLATFORM"
python scripts/22_metagene_merge_density.py
python scripts/23_guitar_metagene.py
python scripts/24_replicate_structure.py
python scripts/25_motif_metagene_bias.py
python scripts/26_negative_controls.py
python scripts/27_window_combination.py
python "$RB/figures/figure1/src/28_fig_tool_counts_replicates.py"
python "$RB/figures/figure5/src/35_fig5a_metric_ranks.py"       # Figure 5A (per-unit metric ranks)
python "$RB/figures/figure5/src/36_fig5b_modratio_hitrate.py"   # Figure 5B (PPV vs. mod ratio)
python "$RB/figures/figure5/src/37_fig5cd_window_sweep.py"      # Figure 5C/5D (window sweep)
python "$RB/figures/figure5/src/38_legacy_fig5_compare.py"      # Figure 5 legacy vs revision values
python "$RB/figures/figure5/src/39_fig5_assembled.py"           # Figure 5 assembled (4 rows, page width)
python "$RB/figures/figureS7/src/53_figS6_tables.py"             # Figure S6 evidence tables (callsets)
python "$RB/figures/figureS7/src/54_figS6_figure.py"             # Figure S6 page (A4, panels A-D)
python "$RB/figures/figureS7/src/55_verify_figS6.py"             # Figure S6 verification (must pass)
python "$RB/figures/figureS8/src/61_figS7_tables.py"             # Figure S7 evidence tables (callsets)
# Figure S7 rows A-F and the page are drawn with Guitar/ggplot under the
# guitar_asm environment (see figures/figureS8/README.md):
#   conda run -n guitar_asm --no-capture-output Rscript "$RB/figures/figureS8/src/62_figS7_guitar.R"
#   conda run -n guitar_asm --no-capture-output Rscript "$RB/figures/figureS8/src/63_figS7_page.R"
python "$RB/figures/figureS8/src/64_verify_figS7.py"             # Figure S7 verification (must pass)
python "$RB/figures/figure6/src/40_fig6_combination.py"         # Figure 6 / S5 evidence tables (incl. 2026-09-21 additions)
python "$RB/figures/figureS6/src/61_figS5_sitequality.py"        # per-set and per-tool site quality (callsets)
python "$RB/figures/figure6/src/65_fig6_combination_page.py"    # main Figure 6, printed 1:1 (6.66 in wide)
python "$RB/figures/figureS6/src/62_figS5_figure.py"             # Figure S5 assembled (A-J, A4 1:1)
python "$RB/figures/figure6/src/66_verify_fig6_page.py"         # Figure 6 verification (scale + >= 7 pt, must pass)
python "$RB/figures/figureS6/src/63_verify_figS5.py"             # Figure S5 verification (must pass)
