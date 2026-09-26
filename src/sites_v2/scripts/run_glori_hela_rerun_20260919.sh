#!/usr/bin/env bash
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
# GLORI HeLa recompute chain v2 — 2026-09-19
# Trigger: Hela_GLORI.bed replaced 06:31 (113,485 -> 112,451).
# v2 = wait-for-concurrency + resume chain: another session was observed
# running the same recomputation (03 done at 07:05, 06 in flight, 04 universe
# rebuild).  This chain WAITS until no related process is alive, then runs
# every step whose output is NOT newer than the new bed.  Already-refreshed
# outputs (produced by the other session) are skipped; the rest are recomputed.
# Serial, fail-fast, singleton.
set -u

BENCH=$RNAMODBENCH_ROOT
S=$BENCH/01_code/code/sites_v2/scripts
R=$BENCH/01_code/code/revision
BED=$BENCH/07_third_party/NGS/GLORI/Hela_GLORI.bed
TS=$(date +%Y%m%d_%H%M%S)
LOG=$BENCH/04_revision_analysis/_archive/glorihela_rerun_logs_$TS
mkdir -p "$LOG"
PIDFILE=$LOG/chain.pid
CHAINLOG=$LOG/chain.log
# nohup shells do not load the conda shell function -> use the env python directly
RUN="$XB/miniconda3/envs/benchmark-revision/bin/python"
MAX_WAIT_SEC=7200   # give the other session up to 2 h

log() { echo "[$(date '+%F %T')] $*" >> "$CHAINLOG"; }

# ---- singleton guard --------------------------------------------------------
if [ -f "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
    echo "another chain is running (pid $(cat "$PIDFILE"))" >&2; exit 1
fi
echo $$ > "$PIDFILE"

# ---- phase 1: wait until no related process is alive ------------------------
log "waiting for other session processes (max ${MAX_WAIT_SEC}s)"
elapsed=0
while [ "$elapsed" -lt "$MAX_WAIT_SEC" ]; do
    if pgrep -f "conda run -n benchmark-revision.*(sites_v2/scripts|revision/NA)" >/dev/null 2>&1 \
       || pgrep -f "sites_v2/scripts/[0-3][0-9]_" >/dev/null 2>&1; then
        sleep 60; elapsed=$((elapsed+60))
    else
        break
    fi
done
if pgrep -f "sites_v2/scripts/" >/dev/null 2>&1; then
    log "WARN: related processes still alive after ${MAX_WAIT_SEC}s; proceeding anyway"
fi
log "phase 1 done after ${elapsed}s"

# ---- helpers ----------------------------------------------------------------
step() {
    local name="$1"; shift
    log "START $name"
    "$RUN" "$@" > "$LOG/${name}.log" 2>&1
    local rc=$?
    if [ "$rc" -ne 0 ]; then
        log "FAIL $name rc=$rc (see $LOG/${name}.log)"; exit "$rc"
    fi
    log "DONE  $name"
}

# run only when outfile is missing or older than the new bed
step_if() {
    local name="$1" out="$2"; shift 2
    if [ -f "$out" ] && [ "$out" -nt "$BED" ]; then
        log "SKIP  $name ($out already newer than bed)"
        return 0
    fi
    step "$name" "$@"
}

EV=$BENCH/04_revision_analysis/sites_v2/evaluation
MR=$BENCH/04_revision_analysis/mod_ratio_replicates
RO=$BENCH/04_revision_analysis/revision_output
F5=$BENCH/04_revision_analysis/fig5_revision
FIG=$BENCH/03_figures/revision_output/figures

# ---- phase 2: resume chain (dependency order) --------------------------------
# 03 was already run on the HeLa samples by the other session (07:05);
# Ath/Mouse callsets never depended on the changed bed.
step 05_offset_audit.py       "$S/05_offset_audit.py"
step 34_export_sites_clean.py "$S/34_export_sites_clean.py"
step_if 06_eval_m6a_glori.py  "$EV/tables/m6a_glori_confusion.tsv"      "$S/06_eval_m6a_glori.py"
step_if 07_eval_controls.py   "$EV/tables/controls_ivt_fpr.tsv"         "$S/07_eval_controls.py"
step_if 08_eval_nonm6a.py     "$EV/tables/hela_nonm6a.tsv"              "$S/08_eval_nonm6a.py"
step_if 09_eval_rna004.py     "$EV/tables/rna004_tool_eval.tsv"         "$S/09_eval_rna004.py"
step 10_qc_reconcile.py       "$S/10_qc_reconcile.py"
step 11_scope_split.py        "$S/11_scope_split.py"
step_if 15_mod_ratio.py       "$MR/tables/mod_ratio_summary_by_group.tsv" "$S/15_mod_ratio_replicate_agreement.py"
step_if 16_mod_ratio_fig.py   "$MR/figures/mod_ratio_regression_S3C_style.pdf" "$S/16_mod_ratio_regression_fig.py"
step_if NA1_window_sweep.py   "$RO/tables/NA1_window_sweep.csv"         "$R/NA1_window_sweep.py"
step_if NA2_glori_stratified.py "$RO/tables/NA2_stratified_recall.csv"  "$R/NA2_glori_stratified.py"
step_if NA2_stratified_figure.py "$FIG/NA2_stratified_sensitivity.pdf"  "$R/NA2_stratified_figure.py"
step_if 35_fig5a.py           "$F5/tables/fig5a_mean_rank_by_group_tool.tsv" "$S/35_fig5a_metric_ranks.py"
step_if 36_fig5b.py           "$F5/tables/fig5b_bin_summary.tsv"        "$S/36_fig5b_modratio_hitrate.py"
step_if 37_fig5cd.py          "$F5/tables/fig5cd_window_sweep_summary.tsv" "$S/37_fig5cd_window_sweep.py"
step_if 38_legacy_compare.py  "$F5/tables/fig5_legacy_vs_revision.tsv"  "$S/38_legacy_fig5_compare.py"
step_if 39_fig5_assembled.py  "$F5/figures/Figure5_rev.pdf"             "$S/39_fig5_assembled.py"
step 24_fig_replicate.py      "$S/24_fig_replicate_structure.py"
step 25_fig_motif_bias.py     "$S/25_fig_motif_metagene_bias.py"
step 27_fig_window_comb.py    "$S/27_fig_window_combination.py" --recompute

log "ALL DONE"
rm -f "$PIDFILE"
