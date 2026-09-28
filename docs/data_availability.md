# Data availability

This is the list of what the repository contains. It is organised around the three
things a benchmark has to expose to be checkable:
**the processed callsets**, **the exact commands and configurations that produced
them**, and **the code that draws the figures**.

---

## 1. Processed modification callsets - `data/callsets/`

Every modification call retained for the analysis, one TSV per
tool × sample × modification, in a single harmonised convention:

```
data/callsets/<platform>/<species>/<dataset_group>/<modification>/<tool>/<sample>.tsv
platform      RNA002 | RNA004
species       Human | Mouse | Arabidopsis | E.coli | Curlcake   (Curlcake = synthetic constructs)
dataset_group HeLa_WT, HeLa_IVT, Arabidopsis_WT, Arabidopsis_KD, Mouse_WT, Mouse_KO,
              E.coli_WT, E.coli_IVT, Curlcake_IVT, Curlcake_m6A, RNA004_*, …
modification  m6A | m5C | Psi | m1Psi | Nm | inosine
```

* 406 files, 4,451,964 site calls (49 tool configurations; the count is
  reproduced per file by `scripts/verify_deposit.py`).
* Coordinates are BED: 0-based, half-open, `end = start + 1`, one base per call.
* Chromosome names are each species' native contig names (`chr1`… for human,
  mouse and *Arabidopsis*; `chromosome` for *E. coli*; `Curlcake1…4` for the
  synthetic constructs).
* A callset with zero retained calls is still present, as a header-only file, so
  that the file set maps one-to-one onto the tool × sample grid.
* Column semantics, including what each tool's `score` actually means, are
  defined in [`sites_columns.md`](sites_columns.md).

The per-file provenance record - which source file each callset was extracted
from, the parser that read it, its SHA-256 fingerprint and its row counts in and
out - is [`../metadata/callsets_index.tsv`](../metadata/callsets_index.tsv).
Row counts were re-verified against the deposited files by
`scripts/verify_deposit.py`.

## 2. Commands, configurations, versions - `tools/`

The workflow that *runs* the tools is also here: `Snakefile` (59 rules),
`config/config.yaml`, `envs/*.yaml` and the `scripts/postprocess_*.py` normalisers -
with [third_party_scripts.md](third_party_scripts.md) covering the tool-owned scripts
the rules call but which cannot be redistributed. The tables below record what was
actually run against it.

| file | content |
|---|---|
| `tools/inventory/tables/per_tool_implementation.csv` / `.md` | for each of the 30 tools and configurations this study used: software/version, model/checkpoint, required input, minimum read coverage, filtering parameters, probability/p-value threshold, multiple-testing correction, default-vs-optimised, and how coordinates were converted and harmonised - each field with its evidence pointer |
| `tools/inventory/tables/command_lines.csv` | 3,045 exact command lines, keyed by tool, dataset, conda environment and the script+line number they were taken from |
| `tools/inventory/tables/software_versions.csv` | 3,479 version records: tool binaries, conda environments and package pins, with the method each version was probed by |
| `tools/inventory/tables/model_checkpoints.csv` | 66 model/checkpoint files used, with the tool and dataset they served |
| `tools/inventory/tables/coverage_report.md` | field-by-field, tool-by-tool completeness of the implementation ledger |
| `tools/configs/` | the configuration files tools were actually driven from (28 xPore YAMLs). xPore is the only tool of this study driven from a file-based configuration; every other tool was driven by command-line flags, which is why they appear in `command_lines.csv` and not here |
| `../envs/` | `envs/*.yaml` are the per-tool environments the Snakemake pipeline creates; `envs/as_run/*.yaml` are the 6 specifications the reported runs actually used (cited by `software_versions.csv`); `envs/analysis/env-*.yml` lock the five environments the benchmark code itself runs in |

How these tables were assembled - and where their evidence is weakest - is in
[`tool_inventory_notes.md`](tool_inventory_notes.md).

## 3. Figure and table code - `figures/`, `tables/`, `src/`

Every main figure (1-8) and supplementary figure (S1-S10) has a directory under
`figures/` holding **all** of its code - the stage that computes its inputs, the panel
scripts, the page assembly and the layout gate - in `src/`, next to the frozen tables it
reads (`tables/`, plus `analysis/` for inputs shared between figures). One figure, one
directory: nothing that draws a published figure lives anywhere else. The panel and page
PDFs are what that code produces, so they are not deposited here:
`bash scripts/run_figures.sh` writes them into `figures/<figure>/figures/`.
`docs/figure_index.md` maps figure - scripts - inputs - outputs, including which
script version is the current one. Supplementary-table builders (Table S1-S12) and
the SI assembly live in `tables/`.

The callset pipeline itself - the stages that turn raw tool output into
`data/callsets/` and the metrics in `data/evaluation/tables/` - is
`src/harmonisation/` (`scripts/run_all.sh`, stages `00`→`34`), and the shared library
the figures import (`config`, `figstyle`, `pagelayout`, `panelpage`, `match`,
`evaluation`, …) is `src/harmonisation/common/`.

## 4. Sample, dataset and run metadata - `metadata/`

* `samples.csv` - one row per analysed sequencing sample: species, dataset group,
  condition, chemistry (RNA002/RNA004), replicate tag, **sequencing unit** and
  **independence class**, SRA/ENA study, BioProject, BioSample and run accessions,
  and which BAM the coverage column was taken from.
* `datasets.csv` - the 12 source studies with their publication and DOI.
* `runs.csv` - 157 run-level records fetched from SRA/ENA metadata.
* `glori_reference.csv` - the GLORI reference sets, their GEO accessions, site
  counts and the selection criterion used to call a site high confidence.
* `replicate_structure.csv` - per tool × dataset group: how many
  replicates, how many studies, whether the group is cross-study and what counts
  as independent.
* `callsets_index.tsv`, `completeness_audit.csv`, `file_inventory.csv` -
  provenance and completeness of the deposited callsets.
* QC ledgers: strand imputation, centre-base filtering, 5-mer centre checks,
  chromosome checks, mouse liftover validation, non-m6A scope.

Two columns in these ledgers refer to the *intermediate* extraction layer rather
than to deposited files - `completeness_audit.csv:callset` and
`file_inventory.csv:source` name working-tree paths under
`$RNAMODBENCH_LOCAL`. They are kept because they are the audit trail from a raw
tool output to a deposited callset; `metadata/callsets_index.tsv` is the
deposited-file-side view of the same relationship.

## 5. Frozen metrics behind the figures - `data/evaluation/tables/`

24 tables: precision/recall/F1/MCC against GLORI at each match window with
bootstrap intervals (`m6a_glori_confusion.tsv`), localisation accuracy versus
window (`m6a_localization_curve.tsv`), replicate agreement
(`reproducibility.tsv`), negative-control false-positive rates
(`controls_ivt_fpr.tsv`), synthetic-truth scoring (`curlcake_truth.tsv`),
knock-down metrics, RNA004 evaluations, non-m6A panels, and the offset/anchor/
pileup QC audits that justify the filtering. [`pipeline.md`](pipeline.md)
describes which stage writes which table.

## 6. Supporting analyses behind individual figures - `analysis/`

Five directories, each named for what it shows rather than for the review point it
answered:

| directory | read by |
|---|---|
| `analysis/nonm6a_false_positives/` | Figure 7 and Figures S7-S8: non-m6A calls against unmodified and knock-out controls |
| `analysis/motif_bias_control/` | Figure 4 and Figure S2: motif preference with the algorithmic-bias control |
| `analysis/mod_ratio_replicates/` | Figure 3 and Figure S3: modification-ratio agreement between replicates |
| `analysis/coverage_rank_stability/`, `analysis/reference_transcript_strata/` | Table S12: the ranking recomputed inside coverage and reference-ratio strata |
| `analysis/supp_table_inputs/` | Tables S6 and S7: per-tool implementation and chemistry-applicability ledgers |

These are frozen inputs, read by the renderers and by `tables/`; their own producers
belong to the `producers` group in [`reproducing.md`](reproducing.md).

## 7. What is deliberately not in this repository

| not deposited | why | where it comes from instead |
|---|---|---|
| raw FASTQ / pod5 / fast5, basecalled BAM/CRAM | size and the upstream data-use terms of each accession | SRA/ENA accessions in `metadata/samples.csv` |
| reference genomes, transcripts and GTF annotations | third-party redistribution | Ensembl / UCSC / TAIR / GRC releases named in `metadata/annotation_summary.csv` |
| tool sources, weights and containers | each tool keeps its own licence | version + checkpoint records in `tools/inventory/tables/software_versions.csv` and `model_checkpoints.csv` |
| GLORI peak files | generated by other groups | GEO accessions in `metadata/glori_reference.csv` |
| the intermediate candidate-universe layer | ~15 GB of derivable data | regenerate with `src/harmonisation/scripts/04_build_universe.py` |
| per-tool raw output directories | intermediate copies of what `callsets` already expresses | regenerate the extraction with stages `00`–`34` |
| the rendered figure and panel PDFs | they are the output of the deposited code, not an input to it | `bash scripts/run_figures.sh`; the published versions are the figures in the manuscript and its SI |

Code in this repository resolves those inputs against a single directory,
`_local/` (or `$RNAMODBENCH_LOCAL`); see [`../_local/README.md`](../_local/README.md).
Most panel and page renderers do not read it - they rebuild from the tables deposited
here - but the producers, the Guitar metagene panels, Figure 1's Illustrator-derived
strip and a few `verify_*` gates do, and `docs/reproducing.md` says which is which.

## 8. Integrity

`scripts/verify_deposit.py` re-reads the deposit and checks that: no personal
absolute paths remain, no private working-tree folder name and no non-English
working note survives in a text source, every Python/R/shell file parses, the row
counts of each deposited callset match `metadata/callsets_index.tsv`, the figure
tree holds 8 main and 10 supplementary directories each named in
`docs/figure_index.md`, and every script path that index names exists. It writes
`metadata/deposited_files.sha256`.
