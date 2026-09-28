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
#: copy of the committed tables, taken before the first step and re-installed after
#: any step that fails, so one aborted producer cannot blank a renderer's input
SNAP="${TMPDIR:-/tmp}/rnamodbench_tables_$$"

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
    guard_tables
    return 0
  fi
  printf '%s\n' "$msg" | tail -2 | sed 's/^/    | /'
  classify "$script" "$msg"
  restore_tables
  return 0
}

snapshot_tables() {  # keep a copy of every frozen table, before any producer runs
  rm -rf "$SNAP"; mkdir -p "$SNAP"
  (cd "$RB" && tar -cf "$SNAP/tables.tar" \
      figures/*/tables figures/*/inputs figures/*/analysis data/evaluation/tables analysis 2>/dev/null) || true
  # the size of every table that holds data at the start.  A step that leaves one
  # with only its header line has produced nothing, whatever its exit status said;
  # sizes are compared so that no table has to be re-read after every step.
  : > "$SNAP/size.txt"
  while IFS= read -r f; do
    case "$f" in *.tsv|*.csv) ;; *) continue ;; esac
    [ -f "$RB/$f" ] || continue
    n=$(stat -c %s "$RB/$f")
    [ "$n" -gt 200 ] && printf '%s\t%s\n' "$f" "$n" >> "$SNAP/size.txt"
  done < <(tar -tf "$SNAP/tables.tar" 2>/dev/null) || true
}

restore_tables() {  # undo a producer that died part-way through rewriting a table
  if [ -f "$SNAP/tables.tar" ]; then
    tar -C "$RB" -xf "$SNAP/tables.tar" 2>/dev/null || true
    RESTORED=1
  fi
}

guard_tables() {  # run after every step: a frozen input must not end up emptied
  [ -s "$SNAP/size.txt" ] || return 0
  local f n now lost=0
  while IFS=$'\t' read -r f n; do
    [ -f "$RB/$f" ] || continue
    now=$(stat -c %s "$RB/$f")
    [ "$now" -lt "$((n / 10))" ] && lost=$((lost + 1))
  done < "$SNAP/size.txt"
  if [ "$lost" -gt 0 ]; then
    echo "   .. $lost table(s) left almost empty; re-installed from the snapshot"
    restore_tables
  fi
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
    restore_tables
  fi
}

cd "$RB"
snapshot_tables
trap 'rm -rf "$SNAP"' EXIT

echo "### Figure 1"
run "fig1 counts"   figures/figure1/src/28_fig_tool_counts_replicates.py
run "fig1 panelC"   figures/figure1/src/35_fig1C_curlcake_density.py
run "fig1 panelD"   figures/figure1/src/36_fig1D_rrach_counts.py
run "fig1 panels"   figures/figure1/src/67_fig1_panels.py
# panel A: the published figure uses the supplied Illustrator strip (72, which reads
# $RNAMODBENCH_LOCAL/figures_original/Fig1a.pdf); the redraw below is what a checkout
# without that artwork composes from, and it is the variant meeting the 7 pt floor.
run "fig1 panelA"   figures/figure1/src/71_fig1_A_redraw.py
run "fig1 page"     figures/figure1/src/69_fig1_page.py
run "fig1 gate"     figures/figure1/src/70_verify_fig1_page.py

echo "### Figure 2"
run "fig2 inputs"   figures/figure2/src/01_fig2_panel_inputs.py
runr "fig2 GO:BP"   figures/figure2/src/02_fig2_go_enrichment.R enrich_r
run "fig2 page"     figures/figure2/src/03_fig2_figure.py
run "fig2 gate"     figures/figure2/src/04_verify_fig2.py

echo "### Figure 3"
run "fig3 A"        figures/figure3/src/56_fig3a_tool_similarity_mds.py
run "fig3 B"        figures/figure3/src/57_fig3b_modratio_wt_treatment.py
runr "fig3 C-D"     figures/figure3/src/58_fig3cd_guitar.R
run "fig3 legend"   figures/figure3/src/61_fig3_legend_band.py
run "fig3 page"     figures/figure3/src/59_fig3_assembled.py
run "fig3 gate"     figures/figure3/src/60_verify_fig3.py

echo "### Figure 4"
run_env motif_analysis "fig4 kl"   figures/figure4/src/fig4_kl.py
run_env motif_analysis "fig4 page" figures/figure4/src/fig4_figures.py

echo "### Figure 5"
run "fig5 A"        figures/figure5/src/35_fig5a_metric_ranks.py
run "fig5 B"        figures/figure5/src/36_fig5b_modratio_hitrate.py
run "fig5 C-D"      figures/figure5/src/37_fig5cd_window_sweep.py
run "fig5 page"     figures/figure5/src/39_fig5_assembled.py

echo "### Figure 6"
run "fig6 tables"   figures/figure6/src/40_fig6_combination.py
run "fig6 page"     figures/figure6/src/65_fig6_combination_page.py
run "fig6 gate"     figures/figure6/src/66_verify_fig6_page.py

echo "### Figure 7"
run "fig7 A-C"      figures/figure7/src/47_fig7_rebuild.py
runr "fig7 D"       figures/figure7/src/48_fig7d_guitar.R
run "fig7 density"  figures/figure7/src/43_fig7_metagene_density.py
if [ -n "$GUITAR" ] && command -v pdflatex >/dev/null; then
  echo "== fig7 page (LaTeX)"
  (cd figures/figure7/src && pdflatex -interaction=nonstopmode assemble_fig7.tex >/dev/null)
else
  echo "-- fig7 page: skipped (needs the R panels above and pdflatex)"
fi
run "fig7 gate"     figures/figure7/src/46_fig7_verify.py

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
runr "figure S1"    figures/figureS1/src/23d_figS1_guitar.R
run_env motif_analysis "figure S2" figures/figureS2/src/01_figS2_motif_inputs.py
run_env motif_analysis "figure S2" figures/figureS2/src/08_figS2v4_stats.py
run_env motif_analysis "figure S2" figures/figureS2/src/09_figS2v4_panels.py
run_env motif_analysis "figure S2" figures/figureS2/src/10_figS2v4_assemble.py
run "figure S3"     figures/figureS3/src/make_s3a_venn.py
run "figure S3"     figures/figureS3/src/make_s3b_ppv.py
run "figure S3"     figures/figureS3/src/make_s3c_reuse_other.py
run "figure S3"     figures/figureS3/src/make_figs3_page.py
run "figure S3 gate" figures/figureS3/src/verify_s3.py
run "figure S4"     figures/figureS4/src/41_figS4_tables.py
run "figure S4"     figures/figureS4/src/42_figS4_figure.py
run "figure S5"     figures/figureS5/src/48_figS5_validation.py
run "figure S5/S4 gate" figures/figureS4/src/49_verify_figS4_S10.py
run "figure S6"     figures/figureS6/src/61_figS6_sitequality.py
run "figure S6"     figures/figureS6/src/62_figS6_figure.py
run "figure S6 gate" figures/figureS6/src/63_verify_figS6.py
run "figure S7"     figures/figureS7/src/53_figS7_tables.py
run "figure S7"     figures/figureS7/src/54_figS7_figure.py
run "figure S7 gate" figures/figureS7/src/55_verify_figS7.py
run "figure S8"     figures/figureS8/src/61_figS8_tables.py
runr "figure S8"    figures/figureS8/src/62_figS8_guitar.R
runr "figure S8"    figures/figureS8/src/63_figS8_page.R
run "figure S8 gate" figures/figureS8/src/64_verify_figS8.py
run "figure S9"     figures/figureS9/src/70_figs9_panels.py
run "figure S9"     figures/figureS9/src/71_figs9_page.py
run "figure S9 gate" figures/figureS9/src/72_verify_figs9.py
run "figure S10"    figures/figureS10/src/60_figS10_tables.py
runr "figure S10"   figures/figureS10/src/23e_figS10_guitar.R
run "figure S10 gate" figures/figureS10/src/61_verify_figS10.py

echo
if [ ${#failed[@]} -gt 0 ]; then
  echo "failed (${#failed[@]}):"; printf '  %s\n' "${failed[@]}"
  if [ -n "${RESTORED:-}" ]; then
    echo "the committed tables were re-installed wherever a step had altered them"
    echo "(a failure, or a table left with only its header line), so no step here drew on"
    echo "input that the deposit does not ship; the panels and pages it did build stay in"
    echo "figures/<figure>/figures/."
  fi
fi
echo "done.  Panels and pages land in figures/<figure>/figures/ as vector PDF;"
echo "compare them against the figures in the manuscript and its supplementary info."
