# Re-running the analysis

## Two roots

Every script resolves two directories and nothing else:

| variable | environment override | default | holds |
|---|---|---|---|
| `_RB` (Python) · `.RB` (R) · `$RB` (shell) | `$RNAMODBENCH_ROOT` | the directory containing `RNAMOD_BENCH_ROOT`, found by walking up from the script | this repository |
| `_XB` · `.XB` · `$XB` | `$RNAMODBENCH_LOCAL` | `$RNAMODBENCH_ROOT/_local` | inputs that are not redistributed (`_local/README.md`) |

A bootstrap block, added when the code was deposited, computes both. No script in
this repository contains an absolute path to a particular machine: set
`$RNAMODBENCH_ROOT` if you move the checkout somewhere else.

## Environments

```bash
conda env create -f envs/analysis/env-benchmark-revision.yml   # python + matplotlib
conda env create -f envs/analysis/env-guitar_asm.yml           # R + Guitar (metagene panels)
conda env create -f envs/analysis/env-motif_analysis.yml       # Figure 4, Figure S2
conda env create -f envs/analysis/env-enrich_r.yml             # Figure 2 GO:BP
```

`envs/` holds the rest: `envs/*.yaml` are the per-tool environments the Snakemake
pipeline creates when it runs the detection tools, and `envs/as_run/*.yaml` are the
specifications the reported runs actually used. Neither is needed to re-render a
figure. A TeX Live installation is needed for the two LaTeX assemblies (Figure 7
page, combined supplementary figures); `pdflatex` is invoked from `PATH`.

## Re-render the figures

Each figure directory carries both kinds of step, and they differ in what they
need:

* **renderers** build panels and pages from the frozen tables shipped beside them
  (`figures/<figure>/tables`, `.../analysis`) and from `data/` and `analysis/`; they
  run on a public checkout;
* **producers** recompute those tables from the intermediate layer (raw tool output,
  candidate universes, region models, the GLORI BEDs) and need
  `$RNAMODBENCH_LOCAL`.

```bash
# one page renderer, straight from the deposited tables
conda run -n benchmark-revision --no-capture-output \
    python figures/figure5/src/39_fig5_assembled.py
# every figure, in dependency order
bash scripts/run_figures.sh
```

`scripts/run_figures.sh` attempts every step in the order given by
[`figure_index.md`](figure_index.md) - producers, panel renderers, page assembly, then
the `verify_*` layout gate. It does not abort on a failing step, and it groups the
failures afterwards by the reason each step actually printed. On a checkout without
`$RNAMODBENCH_LOCAL` two groups are expected: the producers, which recompute the
frozen tables from raw tool output, candidate universes, region models or the GLORI
BEDs; and the gates that additionally re-read a caption or number from the
manuscript or the peer-review correspondence.

What this actually yields on a checkout with no `$RNAMODBENCH_LOCAL`, measured by
running the driver on one: the page renderers of Figures 1, 2, 4, 5 and 6 and of
Figures S2, S4, S5, S6, S7 and S9 each produced their PDF from the committed tables,
with Figure 1 composing panel A from the code-only redraw; Figure 7 assembled its
A-C rows and stopped before the metagene row. What does not build is everything
downstream of a step that needs the private layer: the Guitar metagene panels and
with them the pages of Figures 3, 7 (row D), 8, S1, S8 and S10, Figure S3's GLORI
overlap panels, the producers, and the `verify_*` gates that also re-read the
manuscript. Those are the same steps the table above attributes to
`$RNAMODBENCH_LOCAL`; `GUITAR=1` runs the R panels once the reference layer is there.

Two notes on the gates. They measure the *rendered* file - page budget, tick-label
collisions, the smallest effective font - so their verdict depends on the local
matplotlib and font-metric versions, and a re-rendered page is equivalent to the
published figure rather than byte-identical to it; the figures as accepted are the
ones in the manuscript and its supplementary information, which are not deposited
here. And because several steps read and write the same figure directory, run them
in the documented order: a producer that dies part-way can leave a rewritten table
behind, which is why the driver snapshots the committed tables first and re-installs
them after any failing step. Without the driver,
`git checkout -- figures/<figure>/tables` restores them by hand.

## Rebuild the callsets (needs the private inputs)

`src/harmonisation/` is the pipeline that produced `data/callsets/`. It reads the
per-tool result trees and the reference layer under `$RNAMODBENCH_LOCAL`:

```bash
bash src/harmonisation/scripts/run_all.sh                 # stages 00 -> 10 (extraction, evaluation, QC)
python src/harmonisation/scripts/34_export_callsets.py # writes data/callsets/
bash src/harmonisation/scripts/run_replicate_aware.sh RNA002   # replicate-aware revision analyses
```

Stage order matters and is documented in the header of `src/harmonisation/scripts/run_all.sh`: extraction
(`00`–`03`) → scope split (`11`) → centre-base filter (`33`) → audits (`05`) →
evaluation (`06`–`09`, `12`) → completeness and reconciliation (`13`, `14`,
`29`–`31`, `10`) → figure-ready export. `src/harmonisation/scripts/04_build_universe.py` is skipped by
default because the candidate universes only change when the coverage BAMs change.

Stages `20`/`21` (region model, Guitar BED) and therefore the Guitar panels need
the Ensembl GTFs named in `metadata/annotation_summary.csv`; the pickled region
models and the exported BED files are not redistributed.

## Verifying a checkout

```bash
conda run -n benchmark-revision --no-capture-output python scripts/verify_deposit.py
```

It re-hashes the deposit, checks that no personal absolute path is present, that
every source file parses, and that the row counts of the 406 deposited callsets
agree with `metadata/callsets_index.tsv`.

## Conventions the code assumes

* Match windows and coverage floors are the constants in
  `src/harmonisation/common/config.py` (`w = 2` bp primary, `c = 10` reads); metrics
  are precision against GLORI at that window.
* Replicates are never merged before evaluation: cross-replicate agreement is
  reported as pairwise Jaccard / k-of-n / within-group mean ± SD, and a union over
  replicates is never treated as a sample. The two mouse studies are never pooled.
* Randomised quantities carry explicit seeds (bootstrap, permutation, jitter);
  they are named at the top of the script that uses them.
* `common/pagelayout.py` is a *gate*, not a helper: figure scripts call it and fail
  when text collides with other text, with an axis, or with the page edge, or when
  a font falls below the size floor for the target page.
