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
# Rebuild the whole tool inventory (R1-5 / R1-10 evidence).
#
#   conda activate benchmark-revision
#   bash $RNAMODBENCH_ROOT/tools/inventory/scripts/run_all.sh
#
# Every output stays inside $RNAMODBENCH_ROOT/tools/inventory
# The curated template is only created/updated, never overwritten.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "== 1/5 collect exact command lines (R1-10) =="
python "$HERE/collect_commands.py"

echo "== 2/5 collect parameters / models from artefacts (R1-5) =="
python "$HERE/collect_params.py"

echo "== 3/5 probe installed versions (conda / binary / static) =="
python "$HERE/collect_versions.py"

echo "== 4/5 create or update the manual curation template =="
python "$HERE/make_curated_template.py"

echo "== 5/5 build the supplementary tables =="
python "$HERE/build_tables.py"

echo "Done -> $RNAMODBENCH_ROOT/tools/inventory/tables"
