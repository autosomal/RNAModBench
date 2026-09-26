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
# Completion watcher for refresh-eval-tables (2026-09-19).
# Waits for the CONCURRENT pipeline (other session) to quiesce, then:
#   1. re-runs any of 06-09/10/export whose tables are older than the newest callset
#   2. re-runs the diff report vs the 07:00 backup baseline
#   3. re-runs 40_fig6_combination.py against the final callsets
set -u
ROOT=$RNAMODBENCH_ROOT
SV2=$ROOT/01_code/code/sites_v2
LOG=$ROOT/04_revision_analysis/fig6_revision/logs
EVAL=$ROOT/04_revision_analysis/sites_v2/evaluation
SUM=$LOG/completion_summary.txt
PY="conda run -n benchmark-revision --no-capture-output python"
cd "$SV2" || exit 1
: > "$SUM"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$SUM"; }

log "watching for concurrent pipeline to quiesce"
quiet=0
while [ "$quiet" -lt 5 ]; do
  sleep 60
  # NOTE: must not match this script's own cmdline (its path contains "sites_v2"),
  # hence we look for the conda-run wrappers the pipeline actually uses.
  if pgrep -u "$USER" -f "conda run -n benchmark-revision" >/dev/null 2>&1; then
    quiet=0
    continue
  fi
  newest=$(find "$ROOT/04_revision_analysis/sites_v2" -type f -newermt "-5 minutes" 2>/dev/null | head -1)
  if [ -z "$newest" ]; then
    quiet=$((quiet+1))
  else
    quiet=0
  fi
done
log "quiesced (5 min no related process, no writes) -- proceeding"

log "stage R1: refresh stale eval tables (06-09/10/export, only if older than newest callset)"
NEWCS=$(find "$ROOT/04_revision_analysis/sites_v2/callsets" -name "*.tsv" -newermt "2026-09-19 00:00" -printf '%T@ %p\n' 2>/dev/null | sort -rn | head -1 | cut -d' ' -f1)
log "newest callset epoch: ${NEWCS:-none}"
fresh() { # fresh <table> -> rc0 if table newer than newest callset
  [ -z "$NEWCS" ] && return 0
  m=$(stat -c '%Y' "$1" 2>/dev/null || echo 0)
  awk -v a="$m" -v b="$NEWCS" 'BEGIN{exit !(a>b)}'
}
declare -A STEPS=(
  [06_eval_m6a_glori]="$EVAL/tables/m6a_glori_confusion.tsv"
  [07_eval_controls]="$EVAL/tables/controls_ivt_fpr.tsv"
  [08_eval_nonm6a]="$EVAL/tables/hela_nonm6a.tsv"
  [09_eval_rna004]="$EVAL/tables/rna004_tool_eval.tsv"
)
for s in 06_eval_m6a_glori 07_eval_controls 08_eval_nonm6a 09_eval_rna004; do
  t=${STEPS[$s]}
  if fresh "$t"; then
    log "$s SKIP (table $t already newer than callsets)"
  else
    $PY scripts/$s.py > "$LOG/${s}_refresh2.log" 2>&1
    log "$s rc=$? (was stale)"
  fi
done
$PY scripts/export_figure_ready.py > "$LOG/export_figure_ready_refresh2.log" 2>&1
log "export_figure_ready rc=$?"

log "stage R2: diff report vs 07:00 backup baseline"
BK=$(ls -d "$EVAL"/tables/_bak_human_ensembl_* 2>/dev/null | head -1)
$PY "$ROOT/01_code/reorganization/diff_evaluation_tables.py" --old "$BK" --new "$EVAL/tables" --show 8 \
  > "$LOG/eval_tables_diff_20260919_v2.txt" 2>&1
log "diff rc=$? (1 = differences exist, expected)"

log "stage R3: Fig6/S5 final rerun against settled callsets"
$PY scripts/40_fig6_combination.py > "$LOG/40_fig6_final.log" 2>&1
log "40 rc=$?"

log "COMPLETION DONE"
