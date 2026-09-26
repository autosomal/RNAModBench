# Data availability

This file is the authoritative list of what the repository contains. It is
organised around the three things a benchmark has to expose to be checkable:
**the processed callsets**, **the exact commands and configurations that produced
them**, and **the code that draws the figures**.

---

## 1. Processed modification callsets — `data/sites_clean/`

Every modification call retained for the analysis, one TSV per
tool × sample × modification, in a single harmonised convention:

```
data/sites_clean/<platform>/<species>/<dataset_group>/<modification>/<tool>/<sample>.tsv
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

The per-file provenance record — which source file each callset was extracted
from, the parser that read it, its SHA-256 fingerprint and its row counts in and
out — is [`../metadata/callsets_summary.csv`](../metadata/callsets_summary.csv).
Row counts were re-verified against the deposited files by
`scripts/verify_deposit.py`.

## 2. Commands, configurations, versions — `tools/`

The workflow that *runs* the tools is also here: `Snakefile` (59 rules),
`config/config.yaml`, `envs/*.yaml` and the `scripts/postprocess_*.py` normalisers —
with [third_party_scripts.md](third_party_scripts.md) covering the tool-owned scripts
the rules call but which cannot be redistributed. The tables below record what was
actually run against it.

| file | content |
|---|---|
| `tools/inventory/tables/TI1_per_tool_implementation.csv` / `.md` | for each of 42 tool configurations: software/version, model/checkpoint, required input, minimum read coverage, filtering parameters, probability/p-value threshold, multiple-testing correction, default-vs-optimised, and how coordinates were converted and harmonised — each field with its evidence pointer |
| `tools/inventory/tables/TI2_command_lines.csv` | 3,608 exact command lines, keyed by tool, dataset, conda environment and the script+line number they were taken from |
| `tools/inventory/tables/TI3_software_versions.csv` | 4,272 version records: tool binaries, conda environments and package pins, with the method each version was probed by |
| `tools/inventory/tables/TI4_model_checkpoints.csv` | 99 model/checkpoint files used, with the tool and dataset they served |
| `tools/inventory/tables/TI8_coverage_report.md` | field-by-field, tool-by-tool completeness of TI1 |
| `tools/configs/` | the configuration files tools were actually driven from (22 xPore YAMLs, CHEUI-diff/curlcake YAMLs). Only these two tools used file-based configuration; every other tool was driven by command-line flags, which is why they appear in TI2 and not here |
| `../envs/` | `envs/*.yaml` are the per-tool environments the Snakemake pipeline creates; `envs/as_run/*.yaml` are the 11 specifications the reported runs actually used (cited by `TI3`); `envs/analysis/env-*.yml` lock the five environments the benchmark code itself runs in |

How these tables were assembled — and where their evidence is weakest — is in
[`tool_inventory_notes.md`](tool_inventory_notes.md).

## 3. Figure and table code — `figures/`, `tables/`, `src/`

Every main figure (1–8) and supplementary figure (S1–S10) has a directory under
`figures/` holding its renderer code (`src/`), the tables it reads and writes
(`tables/`), its panels (`figures/`) and the delivered vector figure
(`delivered/`). `docs/figure_index.md` maps figure → scripts → inputs → outputs,
including which script version is the current one. Supplementary-table builders
(Table S1–S12) and the SI assembly live in `tables/`.

The callset pipeline itself — the stages that turn raw tool output into
`data/sites_clean/` — is `src/sites_v2/` (`scripts/run_all.sh`, stages
`00`→`34`), and the shared library the figures import (`config`, `figstyle`,
`pagelayout`, `panelpage`, `match`, `evaluation`, …) is `src/sites_v2/common/`.

## 4. Sample, dataset and run metadata — `metadata/`

* `samples.csv` — one row per analysed sequencing sample: species, dataset group,
  condition, chemistry (RNA002/RNA004), replicate tag, **sequencing unit** and
  **independence class**, SRA/ENA study, BioProject, BioSample and run accessions,
  and which BAM the coverage column was taken from.
* `datasets.csv` — the 12 source studies with their publication and DOI.
* `runs.csv` — 157 run-level records fetched from SRA/ENA metadata.
* `glori_reference.csv` — the GLORI reference sets, their GEO accessions, site
  counts and the selection criterion used to call a site high confidence.
* `sites_v2_replicate_structure.csv` — per tool × dataset group: how many
  replicates, how many studies, whether the group is cross-study and what counts
  as independent.
* `callsets_summary.csv`, `completeness_audit.csv`, `file_inventory.csv` —
  provenance and completeness of the deposited callsets.
* QC ledgers: strand imputation, centre-base filtering, 5-mer centre checks,
  chromosome checks, mouse liftover validation, non-m6A scope.

Two columns in these ledgers refer to the *intermediate* extraction layer rather
than to deposited files — `completeness_audit.csv:callset` and
`file_inventory.csv:source` name working-tree paths under
`$RNAMODBENCH_LOCAL`. They are kept because they are the audit trail from a raw
tool output to a deposited callset; `metadata/callsets_index.tsv` is the
deposited-file-side view of the same relationship.

## 5. Frozen metrics behind the figures — `data/evaluation/tables/`

24 tables: precision/recall/F1/MCC against GLORI at each match window with
bootstrap intervals (`m6a_glori_confusion.tsv`), localisation accuracy versus
window (`m6a_localization_curve.tsv`), replicate agreement
(`reproducibility.tsv`), negative-control false-positive rates
(`controls_ivt_fpr.tsv`), synthetic-truth scoring (`curlcake_truth.tsv`),
knock-down metrics, RNA004 evaluations, non-m6A panels, and the offset/anchor/
pileup QC audits that justify the filtering. [`pipeline.md`](pipeline.md)
describes which stage writes which table.

## 6. What is deliberately not in this repository

| not deposited | why | where it comes from instead |
|---|---|---|
| raw FASTQ / pod5 / fast5, basecalled BAM/CRAM | size and the upstream data-use terms of each accession | SRA/ENA accessions in `metadata/samples.csv` |
| reference genomes, transcripts and GTF annotations | third-party redistribution | Ensembl / UCSC / TAIR / GRC releases named in `metadata/annotation_summary.csv` |
| tool sources, weights and containers | each tool keeps its own licence | version + checkpoint records in `tools/inventory/tables/TI3`, `TI4` |
| GLORI peak files | generated by other groups | GEO accessions in `metadata/glori_reference.csv` |
| the intermediate candidate-universe layer | ~15 GB of derivable data | regenerate with `src/sites_v2/scripts/04_build_universe.py` |
| per-tool raw output directories | intermediate copies of what `sites_clean` already expresses | regenerate the extraction with stages `00`–`34` |

Code in this repository resolves those inputs against a single directory,
`_local/` (or `$RNAMODBENCH_LOCAL`); see [`../_local/README.md`](../_local/README.md).
Nothing in the figure code reads `_local/` — the figures rebuild from the tables
deposited here.

## 7. Integrity

`scripts/verify_deposit.py` re-reads the deposit and checks that: no personal
absolute paths remain, every Python/R/shell file parses, the row counts of each
deposited callset match `metadata/callsets_summary.csv`, and every figure named in
`docs/figure_index.md` has both a script and a delivered file. It writes
`metadata/deposited_files.sha256`.
