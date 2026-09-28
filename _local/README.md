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
- `figures_original/` -- the artwork the published figures were assembled
  from; `72_fig1_A_artwork.py` reads Figure 1's panel A here, and without it
  the page composes from the code-only redraw (`71_fig1_A_redraw.py`).
- `superseded_output/` -- tables written by the earlier form of this
  analysis, which stage `10` still cross-checks against.
- `review/`, `manuscript/`, `submission/` -- peer-review and manuscript
  files that a few verification gates also cross-check.

What each layer unblocks:

| layer | needed by |
|---|---|
| (nothing) | the page and panel renderers, which read `figures/*/tables/`, `figures/*/analysis/` and `data/` |
| `harmonisation/callsets/` | the table-recomputing stages (`28_`, `36_`, `40_`, `41_`, `50_`, `53_figS6`, `61_figS5`, `60_figS9`, `fig4_kl`, S2 inputs) and `src/harmonisation/scripts/run_all.sh` |
| `reference/` | the Guitar metagene panels, and with them the pages of Figures 3, 7, 8, S1, S8 and S10, plus the GLORI overlap panels of Figure S3 |
| `figures_original/` | Figure 1's panel A artwork; without it the page composes from the code-only redraw (`71_fig1_A_redraw.py`) |
| `review/`, `manuscript/` | the caption and number cross-checks inside some `verify_*` gates |

How much of the repository runs with this directory empty?  The panel and page renderers that need no third-party reference layer build their figure from the frozen tables committed here -- measured on such a checkout, that is Figures 1, 2, 4, 5 and 6 and Figures S2, S4, S5, S6, S7 and S9, with Figure 1 composing panel A from the code-only redraw instead of the supplied artwork and Figure 7 stopping at its A-C rows. What does not build is everything downstream of this directory: the Guitar metagene panels, and with them the pages of Figures 3, 8, S1, S8 and S10; Figure S3's GLORI overlap panels; the producers, which recompute those tables from raw tool output and candidate universes; and the `verify_*` gates that also re-read the manuscript. The driver reports each failure with its reason instead of aborting, and re-installs any table a step left empty. Restore the layer here to run those too.
