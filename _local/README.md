# `_local/` -- private inputs, not redistributed

Code in this repository resolves two roots:

| variable | environment | holds |
|---|---|---|
| `_RB` / `$RNAMODBENCH_ROOT` | repository root | everything deposited here |
| `_XB` / `$RNAMODBENCH_LOCAL` | `<root>/_local` | inputs we could not publish |

The second column contains, on the machine that ran the analysis:

- `raw/` -- the per-tool result trees (`result/`, `result_RNA004/`,
  `converted_callsets/`) that stages `00`-`34` parse, and the BAM/FASTQ
  files the coverage columns were read from.  The accessions of those reads
  are deposited: `metadata/samples.csv`, `metadata/runs.csv`.
- `harmonisation/callsets/`, `harmonisation/callsets_extended/`,
  `harmonisation/universe/` -- intermediate extraction and candidate-universe
  output (the universe layer alone is ~15 GB); regenerate with the pipeline
  scripts in `src/harmonisation/`.
- `reference/` -- Ensembl/TAIR/GRC annotations, GTF-derived Guitar BED
  inputs, the pickled region models, and the GLORI BED files behind the
  overlap panels of Figure S3.
- `third_party/` -- other groups' processed data used for cross-checks.
- `superseded_output/` -- tables written by the earlier form of this
  analysis, which stage `10` still cross-checks against.
- `review/`, `manuscript/`, `submission/` -- peer-review and manuscript
  files that a few verification gates also cross-check.

What each layer unblocks:

| layer | needed by |
|---|---|
| (nothing) | the page and panel renderers, which read `figures/*/tables/`, `figures/*/analysis/` and `data/` |
| `harmonisation/callsets/` | the table-recomputing stages (`28_`, `36_`, `40_`, `41_`, `50_`, `53_figS6`, `61_figS5`, `60_figS9`, `fig4_kl`, S2 inputs) and `src/harmonisation/scripts/run_all.sh` |
| `reference/` | the Guitar metagene panels (Figures 3, 7, 8, S1, S8, S10) and the GLORI overlap panels of Figure S3 |
| `review/`, `manuscript/` | the caption and number cross-checks inside some `verify_*` gates |

How much of the repository runs with this directory empty?  Every panel and page renderer that needs no third-party reference layer builds its figure from the frozen tables committed here -- that is, all of them except the Guitar metagene panels (Figures 3, 7, 8, S1, S8, S10) and Figure S3's GLORI overlap panels, which read data other groups published. The steps that do need this directory are the producers (they recompute those tables from raw tool output and candidate universes) and a few `verify_*` gates that also re-read the manuscript; the driver reports each failure with its reason instead of aborting. Restore the layer here to run those too.
