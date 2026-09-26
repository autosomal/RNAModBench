#!/usr/bin/env bash
# Redraw the GUITAR metagene panels with the R Guitar package, on the
# replicate-merged BEDs written by 21b_export_guitar_bed.py.
#
#   conda activate guitar_asm
#   bash run_guitar_R.sh [Arabidopsis Mouse Human ...]
#
# One species at a time would serialise the slow part (building the transcript
# model), so the species run in parallel; inside a species the calls are
# sequential because each reuses the same cached TxDb.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT=$RNAMODBENCH_LOCAL/sites_v2/guitar_metagene/figures
LOGS=$RNAMODBENCH_LOCAL/sites_v2/guitar_metagene/logs
mkdir -p "$OUT" "$LOGS"
if [ "$#" -gt 0 ]; then set -- "$@"; else set -- Arabidopsis Mouse Human; fi

run() { Rscript "$HERE/23b_guitar_metagene.R" "$@"; }

for species in "$@"; do
  (
    set -e
    # panel C: libraries pooled over tools, majority consensus
    run --species "$species" --mode condition --merge majority --txtype mrna
    run --species "$species" --mode condition --merge majority --txtype ncrna
    # merge-rule sensitivity: union / majority / intersection in one axis
    run --species "$species" --mode condition \
        --merge union,majority,intersection --txtype mrna
    # panel D: one curve per tool x library
    run --species "$species" --mode pertool --merge majority --txtype mrna
    # replicate spread before any consensus
    run --species "$species" --mode replicates --merge majority --txtype mrna
  ) > "$LOGS/guitar_R_${species}.log" 2>&1 &
done
wait
echo "GUITAR R figures done $(date '+%F %T')"
