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
# Completion chain v2 (2026-09-20): wait for the concurrent pipeline to quiesce,
# then refresh only stale eval tables (mtime vs newest callset), diff v2,
# affected figures and the final Fig6/S5 run (40).
set -u
ROOT=$RNAMODBENCH_ROOT
SV2=$ROOT/01_code/code/sites_v2
LOG=$ROOT/04_revision_analysis/fig6_revision/logs
EVAL=$ROOT/04_revision_analysis/sites_v2/evaluation
SUM=$LOG/completion2_summary.txt
PY="conda run -n benchmark-revision --no-capture-output python"
cd "$SV2" || exit 1
: > "$SUM"
log() { echo "[$(date '+%F %T')] $*" >> "$SUM"; }
log "start (v2: quiescence wait + mtime-smart refresh)"

# --- guard: 3 consecutive clean minutes (no pipeline process, no fresh writes) ---
# "--now" skips the wait (used when the concurrent pipeline is known to be done).
clean=0; waited=0
if [ "${1:-}" != "--now" ]; then
while [ "$clean" -lt 3 ] && [ "$waited" -lt 120 ]; do
  sleep 60; waited=$((waited+1))
  busy=0
  pgrep -u "$USER" -f "conda run -n benchmark-revision" >/dev/null 2>&1 && busy=1
  if [ "$busy" -eq 0 ]; then
    f=$(find "$ROOT/04_revision_analysis/sites_v2/callsets" "$ROOT/04_revision_analysis/sites_v2/evaluation" -type f -newermt "-6 minutes" 2>/dev/null | head -1)
    [ -n "$f" ] && busy=1
  fi
  if [ "$busy" -eq 0 ]; then clean=$((clean+1)); else clean=0; fi
done
if [ "$waited" -ge 120 ]; then log "GAVE UP: pipeline still busy after 2h"; exit 1; fi
log "quiesced after ${waited} min"
fi

# --- stage R1: mtime-smart eval refresh (only stale scripts) ---
NEWCS=$(find "$ROOT/04_revision_analysis/sites_v2/callsets" -name "*.tsv" -printf '%T@ %p\n' 2>/dev/null | sort -rn | head -1 | cut -d' ' -f1)
log "newest callset epoch: ${NEWCS:-none}"
fresh() { [ -z "$NEWCS" ] && return 0; m=$(stat -c '%Y' "$1" 2>/dev/null || echo 0); awk -v a="$m" -v b="$NEWCS" 'BEGIN{exit !(a>b)}'; }
declare -A OUT=(
  [06_eval_m6a_glori]="$EVAL/tables/m6a_glori_confusion.tsv"
  [07_eval_controls]="$EVAL/tables/controls_ivt_fpr.tsv"
  [08_eval_nonm6a]="$EVAL/tables/hela_nonm6a.tsv"
  [09_eval_rna004]="$EVAL/tables/rna004_tool_eval.tsv"
  [10_qc_reconcile]="$EVAL/tables/legacy_published_check.tsv"
)
for s in 06_eval_m6a_glori 07_eval_controls 08_eval_nonm6a 09_eval_rna004 10_qc_reconcile; do
  if fresh "${OUT[$s]}"; then
    log "$s SKIP (output already newer than callsets)"
  else
    $PY scripts/$s.py > "$LOG/${s}_refresh3.log" 2>&1 < /dev/null
    log "$s rc=$? (was stale)"
  fi
done
$PY scripts/export_figure_ready.py > "$LOG/export_refresh3.log" 2>&1 < /dev/null
log "export rc=$?"

# --- stage R2: diff v2 ---
BK=$(ls -d "$EVAL"/tables/_bak_human_ensembl_* 2>/dev/null | head -1)
$PY "$ROOT/01_code/reorganization/diff_evaluation_tables.py" --old "$BK" --new "$EVAL/tables" --show 8 \
  > "$LOG/eval_tables_diff_20260920.txt" 2>&1 < /dev/null
log "diff rc=$? (1 = differences exist, expected)"

# --- stage R3: affected figures + final Fig6/S5 ---
for s in 24_fig_replicate_structure 26_fig_negative_controls; do
  $PY scripts/$s.py > "$LOG/${s}_refresh3.log" 2>&1 < /dev/null
  log "$s rc=$?"
done
$PY scripts/27_fig_window_combination.py --recompute > "$LOG/27_refresh3.log" 2>&1 < /dev/null
log "27 rc=$?"
$PY scripts/35_fig5a_metric_ranks.py > "$LOG/35_refresh3.log" 2>&1 < /dev/null
log "35 rc=$?"
$PY scripts/40_fig6_combination.py > "$LOG/40_fig6_final.log" 2>&1 < /dev/null
log "40 rc=$?"
touch "$LOG/completion_done.txt"
log "COMPLETION DONE"
