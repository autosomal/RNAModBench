#!/bin/bash
# --- RNAModBench path bootstrap (added when this file was deposited) ----------
if [ -z "${RNAMODBENCH_ROOT:-}" ]; then
  RNAMODBENCH_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)"
  while [ "$RNAMODBENCH_ROOT" != "/" ] && [ ! -f "$RNAMODBENCH_ROOT/RNAMOD_BENCH_ROOT" ]; do
    RNAMODBENCH_ROOT="$(dirname -- "$RNAMODBENCH_ROOT")"
  done
fi
RB="$RNAMODBENCH_ROOT"
XB="${RNAMODBENCH_LOCAL:-$RB/_local}"
export RNAMODBENCH_ROOT RB XB
# --------------------------------------------------------------------------- #
# Human universe Ensembl-112 refresh chain (2026-09-19, Fig6 rebuild plan stage 1-2)
# Stage A: 05b generic universe export for 8 Human samples (06-09 read this family)
# Stage B: backup eval tables -> 06 07 08 09 10 export_figure_ready -> diff report
# Stage C: dependent figures 24 26 27(--recompute) 35
# Strictly serial; each step logs its own rc into chain_summary.txt.
set -u
ROOT=$RNAMODBENCH_ROOT
SV2=$ROOT/01_code/code/sites_v2
LOG=$ROOT/04_revision_analysis/fig6_revision/logs
EVAL=$ROOT/04_revision_analysis/sites_v2/evaluation
SUM=$LOG/chain_summary.txt
PY="conda run -n benchmark-revision --no-capture-output python"
HUMAN=(HeLa_WT1 HeLa_WT2 HeLa_WT3 HeLa_IVT_rep1 HeLa_IVT_rep2 HeLa_IVT_rep3 HeLa_RNA004_WT HeLa_RNA004_IVT)
cd "$SV2" || exit 1
: > "$SUM"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$SUM"; }

log "Stage A: 05b generic universe export (8 Human samples, out=UNIVERSE_ROOT)"
ARGS=(); for s in "${HUMAN[@]}"; do ARGS+=(--sample "$s"); done
$PY scripts/05b_export_generic_universe.py "${ARGS[@]}" > "$LOG/05b_human_generic.log" 2>&1
log "05b rc=$?"

log "Stage B: backup eval tables + refresh 06-09/10/export"
TS=$(date +%Y%m%d_%H%M%S)
BK=$EVAL/tables/_bak_human_ensembl_$TS
mkdir -p "$BK"
cp -p "$EVAL"/tables/*.tsv "$BK"/ 2>/dev/null
echo "backup=$BK" > "$LOG/eval_refresh_backup.txt"
log "backup -> $BK"
for s in 06_eval_m6a_glori 07_eval_controls 08_eval_nonm6a 09_eval_rna004 10_qc_reconcile; do
  $PY scripts/$s.py > "$LOG/${s}_refresh.log" 2>&1
  log "$s rc=$?"
done
$PY scripts/export_figure_ready.py > "$LOG/export_figure_ready_refresh.log" 2>&1
log "export_figure_ready rc=$?"

log "diff report old->new"
$PY "$ROOT/01_code/reorganization/diff_evaluation_tables.py" --old "$BK" --new "$EVAL/tables" --show 8 \
  > "$LOG/eval_tables_diff_20260919.txt" 2>&1
log "diff rc=$?"

log "Stage C: dependent figures 24 26 27(--recompute) 35"
for s in 24_fig_replicate_structure 26_fig_negative_controls; do
  $PY scripts/$s.py > "$LOG/${s}_refresh.log" 2>&1
  log "$s rc=$?"
done
$PY scripts/27_fig_window_combination.py --recompute > "$LOG/27_fig_window_combination_refresh.log" 2>&1
log "27 rc=$?"
$PY scripts/35_fig5a_metric_ranks.py > "$LOG/35_fig5a_metric_ranks_refresh.log" 2>&1
log "35 rc=$?"

log "verify paths"
$PY "$ROOT/01_code/reorganization/verify_sites_v2_paths.py" > "$LOG/verify_paths_after_refresh.log" 2>&1
log "verify rc=$?"
log "CHAIN DONE"
