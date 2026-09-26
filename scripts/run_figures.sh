#!/usr/bin/env bash
# Re-render the benchmark figures from the deposited tables.
#
#   bash scripts/run_figures.sh              # every figure that needs no private input
#   bash scripts/run_figures.sh figure4 S3   # only these
#   GUITAR=1 bash scripts/run_figures.sh     # also the R/Guitar metagene panels
#
# Figures drawn in R need the Guitar BED inputs (or, for Figure S1, the Ensembl
# GTFs behind $RNAMODBENCH_LOCAL); they are skipped unless GUITAR=1 is set, and the
# reason is printed.  Order within a figure is producer -> panels -> page -> gate.
set -euo pipefail

if [ -z "${RNAMODBENCH_ROOT:-}" ]; then
  RNAMODBENCH_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
  while [ "$RNAMODBENCH_ROOT" != "/" ] && [ ! -f "$RNAMODBENCH_ROOT/RNAMOD_BENCH_ROOT" ]; do
    RNAMODBENCH_ROOT="$(dirname -- "$RNAMODBENCH_ROOT")"
  done
fi
RB="$RNAMODBENCH_ROOT"
XB="${RNAMODBENCH_LOCAL:-$RB/_local}"
export RNAMODBENCH_ROOT RB XB

PY="${PY:-python}"
RENVR="${PYENV:-benchmark-revision}"
GUITAR="${GUITAR:-}"
WANT=("$@")
FAILED=()
# conda is preferred when it is on PATH; otherwise run with an already activated
# environment (export CONDA= to force that path).
CONDA="${CONDA-$(command -v conda 2>/dev/null || true)}"

with_env() {  # with_env <conda-env> <command...>
  local env="$1"; shift
  if [ -n "$CONDA" ]; then
    "$CONDA" run -n "$env" --no-capture-output "$@"
  else
    "$@"
  fi
}

wanted() {  # wanted <label> -- false when an explicit figure list was given and this is not in it
  local label="$1" w
  [ ${#WANT[@]} -eq 0 ] && return 0
  for w in "${WANT[@]}"; do
    [ "$w" = "all" ] && return 0
    case "$label" in *"$w"*) return 0 ;; esac
  done
  return 1
}

# Every step is attempted in dependency order; a failure is reported and the run
# continues, and at the end the failures are grouped by the reason each one actually
# printed.  Two groups are expected on a public checkout, where the intermediate
# layer under $RNAMODBENCH_LOCAL (raw tool output, candidate universes, region
# models, GLORI BEDs) is absent:
#   * producers, which recompute the frozen tables from that layer;
#   * gates that cross-check a caption or number in the manuscript or the
#     peer-review correspondence, which are not deposited.
# A producer that dies part-way can leave a rewritten table behind; the tables are
# committed, so restore them with:  git checkout -- figures/<figure>/tables
skipped_ran=(); failed=()

classify() {  # classify <script> <captured output>
  local script="$1" msg="$2" name
  name=$(basename "$script")
  if printf '%s' "$msg" | grep -qE '_local|harmonisation/callsets|universe/|third_party|regionmodel|guitar_metagene'; then
    failed+=("$name  <- reads the intermediate layer (\$RNAMODBENCH_LOCAL)")
  elif printf '%s' "$msg" | grep -qE 'review/|manuscript/|submission/|response_letter'; then
    failed+=("$name  <- cross-checks a non-deposited manuscript/review file")
  elif printf '%s' "$msg" | grep -qE 'KeyError|IndexError|AssertionError|EmptyDataError|no files to open|not found under'; then
    failed+=("$name  <- empty input from the intermediate layer")
  else
    failed+=("$name")
  fi
}

run_env() {  # run_env <conda-env> <label> <script> [args...]
  local env="$1" label="$2" script="$3"; shift 3
  wanted "$label" || return 0
  if [ ! -f "$script" ]; then
    echo "-- $label: $script not found (skipped)"; return 0
  fi
  echo "== $label: $(basename "$script") (env: $env)"
  local msg
  if msg=$(with_env "$env" "$PY" "$script" "$@" 2>&1); then
    return 0
  fi
  printf '%s\n' "$msg" | tail -2 | sed 's/^/    | /'
  classify "$script" "$msg"
  return 0
}

run() {  # run <label> <script> [args...] -- default analysis environment
  run_env "$RENVR" "$@"
}

runr() {  # R/Guitar panels, opt-in
  local label="$1" script="$2" env="${3:-guitar_asm}"
  wanted "$label" || return 0
  if [ -z "$GUITAR" ]; then
    echo "-- $label: skipped (needs Guitar; re-run with GUITAR=1)"; return 0
  fi
  echo "== $label: $(basename "$script") (R)"
  if ! with_env "$env" Rscript "$script"; then
    echo "   .. failed"
    FAILED+=("$label -- $(basename "$script")")
  fi
}

cd "$RB"

echo "### Figure 1"
run "fig1 counts"   src/harmonisation/scripts/28_fig_tool_counts_replicates.py
run "fig1 panelC"   src/harmonisation/scripts/35_fig1C_curlcake_density.py
run "fig1 panelD"   src/harmonisation/scripts/36_fig1D_rrach_counts.py
run "fig1 panels"   src/harmonisation/scripts/67_fig1_panels.py
run "fig1 page"     src/harmonisation/scripts/69_fig1_page.py
run "fig1 gate"     src/harmonisation/scripts/70_verify_fig1_page.py

echo "### Figure 2"
run "fig2 inputs"   figures/figure2/src/01_fig2_panel_inputs.py
runr "fig2 GO:BP"   figures/figure2/src/02_fig2_go_enrichment.R enrich_r
run "fig2 page"     figures/figure2/src/03_fig2_figure.py
run "fig2 gate"     figures/figure2/src/04_verify_fig2.py

echo "### Figure 3"
run "fig3 A"        src/harmonisation/scripts/56_fig3a_tool_similarity_mds.py
run "fig3 B"        src/harmonisation/scripts/57_fig3b_modratio_wt_treatment.py
runr "fig3 C-D"     src/harmonisation/scripts/58_fig3cd_guitar.R
run "fig3 legend"   src/harmonisation/scripts/61_fig3_legend_band.py
run "fig3 page"     src/harmonisation/scripts/59_fig3_assembled.py
run "fig3 gate"     src/harmonisation/scripts/60_verify_fig3.py

echo "### Figure 4"
run_env motif_analysis "fig4 kl"   figures/figure4/src/fig4_kl.py
run_env motif_analysis "fig4 page" figures/figure4/src/fig4_figures.py

echo "### Figure 5"
run "fig5 A"        src/harmonisation/scripts/35_fig5a_metric_ranks.py
run "fig5 B"        src/harmonisation/scripts/36_fig5b_modratio_hitrate.py
run "fig5 C-D"      src/harmonisation/scripts/37_fig5cd_window_sweep.py
run "fig5 page"     src/harmonisation/scripts/39_fig5_assembled.py

echo "### Figure 6"
run "fig6 tables"   src/harmonisation/scripts/40_fig6_combination.py
run "fig6 page"     src/harmonisation/scripts/65_fig6_combination_page.py
run "fig6 gate"     src/harmonisation/scripts/66_verify_fig6_page.py

echo "### Figure 7"
run "fig7 A-C"      src/harmonisation/scripts/47_fig7_rebuild.py
runr "fig7 D"       src/harmonisation/scripts/48_fig7d_guitar.R
run "fig7 density"  src/harmonisation/scripts/43_fig7_metagene_density.py
if [ -n "$GUITAR" ] && command -v pdflatex >/dev/null; then
  echo "== fig7 page (LaTeX)"
  (cd figures/figure7/src && pdflatex -interaction=nonstopmode assemble_fig7.tex >/dev/null)
else
  echo "-- fig7 page: skipped (needs the R panels above and pdflatex)"
fi
run "fig7 gate"     src/harmonisation/scripts/46_fig7_verify.py

echo "### Figure 8"
run "fig8 tables"   figures/figure8/src/50_fig8_tables.py
run "fig8 A"        figures/figure8/src/57_fig8a_counts.py
run "fig8 B"        figures/figure8/src/59_fig8b_fpr.py
run "fig8 C"        figures/figure8/src/58_fig8c_ppv.py
run "fig8 D"        figures/figure8/src/67_fig8d_tradeoff.py
runr "fig8 E"       figures/figure8/src/69_fig8e_guitar.R
run "fig8 page"     figures/figure8/src/68_assemble_a4.py
run "fig8 gate"     figures/figure8/src/53_verify_fig8.py

echo "### Supplementary figures"
runr "figure S1"    src/harmonisation/scripts/23d_figS1_guitar.R
run_env motif_analysis "figure S2" figures/figureS2/src/01_figS2_motif_inputs.py
run_env motif_analysis "figure S2" figures/figureS2/src/08_figS2v4_stats.py
run_env motif_analysis "figure S2" figures/figureS2/src/09_figS2v4_panels.py
run_env motif_analysis "figure S2" figures/figureS2/src/10_figS2v4_assemble.py
run "figure S3"     figures/figureS3/src/make_s3a_venn.py
run "figure S3"     figures/figureS3/src/make_s3b_ppv.py
run "figure S3"     figures/figureS3/src/make_s3c_reuse_other.py
run "figure S3"     figures/figureS3/src/make_figs3_page.py
run "figure S3 gate" figures/figureS3/src/verify_s3.py
run "figure S4"     src/harmonisation/scripts/41_figS4_tables.py
run "figure S4"     src/harmonisation/scripts/42_figS4_figure.py
run "figure S5"     src/harmonisation/scripts/48_figS10_validation.py
run "figure S5/S4 gate" src/harmonisation/scripts/49_verify_figS4_S10.py
run "figure S6"     src/harmonisation/scripts/61_figS5_sitequality.py
run "figure S6"     src/harmonisation/scripts/62_figS5_figure.py
run "figure S6 gate" src/harmonisation/scripts/63_verify_figS5.py
run "figure S7"     src/harmonisation/scripts/53_figS6_tables.py
run "figure S7"     src/harmonisation/scripts/54_figS6_figure.py
run "figure S7 gate" src/harmonisation/scripts/55_verify_figS6.py
run "figure S8"     src/harmonisation/scripts/61_figS7_tables.py
runr "figure S8"    src/harmonisation/scripts/62_figS7_guitar.R
runr "figure S8"    src/harmonisation/scripts/63_figS7_page.R
run "figure S8 gate" src/harmonisation/scripts/64_verify_figS7.py
run "figure S9"     figures/figureS9/src/70_figs8_panels.py
run "figure S9"     figures/figureS9/src/71_figs8_page.py
run "figure S9 gate" figures/figureS9/src/72_verify_figs8.py
run "figure S10"    figures/figureS10/src/60_figS9_tables.py
runr "figure S10"   src/harmonisation/scripts/23e_figS9_guitar.R
run "figure S10 gate" figures/figureS10/src/61_verify_figS9.py

echo
if [ ${#failed[@]} -gt 0 ]; then
  echo "failed (${#failed[@]}):"; printf '  %s\n' "${failed[@]}"
fi
echo "done.  Panels land in figures/<figure>/figures/; compare against figures/<figure>/delivered/."
