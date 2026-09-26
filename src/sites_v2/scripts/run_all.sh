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
# sites_v2 -- full rebuild, in order.
#
#   conda activate benchmark-revision
#   bash $RNAMODBENCH_ROOT/src/sites_v2/scripts/run_all.sh
#
# Steps (each is idempotent and can be re-run alone):
#   00 registry  -> manifest/sample_registry.csv, sample_tool_registry.csv, pending.csv
#   01 extract   -> callsets/<platform>/<species>/<group>/<mod_type>/<tool>/<sample>.tsv
#   01b fill     -> per-replicate callsets rebuilt from raw via the code_user
#                   conversion recipes (common/legacy_liftover.BUILDERS); writes a
#                   header-only callset when a raw file yields 0 sites (Total_Detected=0)
#   02b validate -> manifest/liftover_validation.csv (rebuilt == legacy, row-for-row)
#   32 impute    -> fills the transcript strand for callsets that report none
#                   (Nanom6A writes '*' for every row) from the tool's own
#                   read-alignment BED (config.READ_STRAND_BED) and, failing
#                   that, from the species exon annotation.  MUST run after 02b
#                   (which compares 01's output against the legacy layer) and
#                   before 03 (which is strand-aware).
#   03 annotate  -> adds 5mer / DRACH / centre / coverage / GLORI columns in place
#                   plus ``center_status`` (ok / off_base / unknown_strand /
#                   no_expectation) -- the strict centre-base verdict.
#   01 note      -> RNA004 Dorado sources declare mod_type="auto": the
#                   modification is read from the pileup's modkit `name` code
#                   (registry.peek_modkit_codes), a file with several codes yields
#                   one callset per code, and the m6A channel keeps its historical
#                   tool label while the other channels get _Psi/_m5C/_inosine
#                   suffixes (see sites_v2/README.md §4a)
#   11 scope     -> DELETE everything outside the manuscript's tool scope (non-m6A
#                   on Arabidopsis/Mouse/E.coli AND tools the paper never used:
#                   differr / EpiNano_SVM / Tombo_com / CHEUI-diff / mAFiA /
#                   CHEUI on Curlcake) + prune every empty dir;
#                   manifest/{nonm6a_scope,nonm6a_deleted,out_of_scope_tools_deleted}.csv
#                   (moved 2026-09-16: runs BEFORE 05/04/06-09, see ORDER MATTERS)
#   33 filter    -> HARD centre-base filter (2026-09-18): deletes every call whose
#                   reference base is not the modification's base in transcript
#                   orientation (Nm exempt).  Writes
#                   evaluation/tables/{center_base_strict_audit.tsv,
#                   center_filter_removed.csv} + manifest/center_base_filter.csv.
#                   This is the only step that removes rows from a callset.
#   05 audit     -> evaluation/tables/offset_audit.tsv (+ per-row offset_flag)
#   04 universe  -> universe/<platform>/<species>/<sample>__universe.tsv[.gz]
#                   (only needs the coverage BAM + annotation, never a callset:
#                   re-run it only when the BAMs or the annotation change)
#   06-09 eval   -> evaluation/tables/*.tsv
#   12 fig7      -> evaluation/tables/nonm6a_fig7_summary.tsv (HeLa + Curlcake)
#   13 audit     -> manifest/completeness_audit.csv: every (sample,tool) is
#                   filled / ok_zero / raw_missing / out_of_scope_deleted (+ flags
#                   any raw_present_unfilled anomaly)
#   14 coverage  -> manifest/legacy_coverage_audit.csv: independent reconciliation
#                   against the 7 legacy cp.sh copy maps + the assembled output/
#                   tree (any legacy-used file without a rebuilt callset)
#   29 anchor    -> evaluation/tables/anchor_audit.tsv (+ nanom6a_anchor_audit.tsv):
#                   report-only check that every callset sits on a
#                   modification-compatible base; flags sources whose coordinates
#                   are off-axis (added 2026-09-18 after the f5c +7-bp anchor bug,
#                   see code/f5c_mode/README.md).  Never filters anything.
#   30 pileup    -> evaluation/tables/pileup_call_filter_audit.tsv: report every
#                   declared Dorado pileup source (rows / no-call rows / rows below
#                   the declared coverage floor / retained share) and assert that no
#                   Dorado callset carries a no-call row; the Curlcake RNA004 family
#                   was extracted unfiltered until 2026-09-18 (m6A/m5C/Psi all sat
#                   on the expected base only ~58 % of the time).  ``--no-fail`` in
#                   the chain (report first, products stay complete); run it WITHOUT
#                   the flag for the acceptance assertion.
#   31 persite   -> evaluation/tables/persite_reference_audit.tsv: for every
#                   EpiNano ``*per.site.csv`` check that the ``base`` column follows
#                   the genome its BAM was aligned to, and say whether a low
#                   agreement is a coordinate offset or a wrong reference (the
#                   fip37 tables were built against a human FASTA; 2026-09-18).
#   10 qc        -> evaluation/qc_report.md + legacy reconciliation
#
# ORDER MATTERS: 11 (scope delete + prune) runs AFTER 01b/03 so no empty non-m6A
# directory can survive a clean rebuild, and -- since 2026-09-16 -- BEFORE
# 05/04/06-09, so the audit and every evaluation table are built from the scope-
# cleaned callset tree and can never contain a callset that is deleted later.
# (The old order ran 11 after 09, which left rows for deleted callsets in
# offset_audit.tsv / offset_histogram.tsv / legacy_reconciliation.tsv.)  01 does
# not know about scope, so it re-creates the out-of-scope files on every clean
# rebuild; 11 is idempotent (--dry-run / --undo) and removes them again.
# 13/14/10 reconcile last.
set -euo pipefail

CODE=$RNAMODBENCH_ROOT/src/sites_v2
LOGS=$RNAMODBENCH_ROOT/data/evaluation/logs
mkdir -p "$LOGS"

run() {
  echo "=== $1 $(date '+%F %T') ==="
  python "$CODE/scripts/$1" "${@:2}"
}

run 00_build_registry.py
run 01_extract_callsets.py
run 01b_liftover_missing.py
run 02b_validate_liftover.py
run 32_impute_strand.py
run 03_annotate_callsets.py
run 11_scope_split.py
run 33_center_base_filter.py --apply
run 05_offset_audit.py
# 04 only rebuilds universes from the coverage BAMs; skip it unless those change.
# run 04_build_universe.py
run 06_eval_m6a_glori.py
run 07_eval_controls.py
run 08_eval_nonm6a.py
run 09_eval_rna004.py
run 12_nonm6a_fig7.py
run 13_completeness_audit.py
run 14_legacy_coverage_audit.py
run 29_anchor_audit.py --fail-on-bad
run 30_pileup_call_filter_audit.py
run 31_persite_reference_audit.py
run 10_qc_reconcile.py
run export_figure_ready.py

echo "sites_v2 rebuild finished $(date '+%F %T')"
