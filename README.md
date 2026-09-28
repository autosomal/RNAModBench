# RNAModBench

Nanopore direct-RNA **RNA-modification detection benchmark**. One repository, two
roles:

* **a pipeline** - Snakemake rules, per-tool conda environments and post-processing
  scripts that run 15+ modification callers over any species and normalise their
  output into one BED-like format;
* **a deposit** - everything needed to verify the published benchmark without
  running anything: the processed modification callsets, the exact command lines /
  configurations / versions behind them, and the code that draws every figure and
  table.

Raw reads are not redistributed; they are public in SRA/ENA/GEO and every sample's
accession is recorded in [`metadata/samples.csv`](metadata/samples.csv).

## What is in here

| | count | where |
|---|---|---|
| per-site callsets (TSV) | 406 files / 4,451,964 site calls | [`data/callsets/`](data/callsets/README.md) |
| modifications covered | m6A, m5C, Ψ, m1Ψ, Nm, inosine | path component |
| tool configurations | 49 (Dorado model variants included) | `data/callsets/<platform>/<species>/<group>/<mod>/<tool>/` |
| frozen evaluation tables | 24 | `data/evaluation/tables/` |
| sample / study / run metadata | 28 samples, 12 studies, 157 runs | `metadata/` |
| exact command lines recorded | 3,045 | `tools/inventory/tables/command_lines.csv` |
| per-tool implementation fields | 42 configurations × 10 fields | `tools/inventory/tables/per_tool_implementation.md` |
| figure, table and callset code | 176 Python + 14 R + 5 shell sources | `src/harmonisation/`, `figures/*/src/`, `tables/`, `tools/inventory/scripts/` |
| rendered figure PDFs | produced by that code, not deposited | `bash scripts/run_figures.sh` |
| pipeline rules / tool envs / post-processors | 59 rules, 15 envs, 20 scripts | `Snakefile`, `envs/`, `scripts/` |

## Repository map

```
RNAModBench/
├── RNAMOD_BENCH_ROOT      repository-root marker (scripts locate the root through it)
├── Snakefile  config/  envs/  setup.py
│                          the detection pipeline: run the tools yourself
├── scripts/               pipeline post-processors, plus deposit helpers
│   ├── postprocess_*.py …           (tool output -> common format)
│   ├── run_figures.sh               re-render the deposited figures
│   └── verify_deposit.py            integrity + de-identification check
├── data/
│   ├── callsets/       processed callsets (the deposit's core)
│   └── evaluation/tables/ frozen metrics behind the tables and figures
├── metadata/              samples, studies, accessions, callset provenance and QC
├── tools/
│   ├── inventory/         per-tool implementation · command lines · versions · model checkpoints
│   └── configs/           the configuration files tools were driven from
├── src/harmonisation/     the callset, evaluation and QC stages (00-34) + shared library
├── figures/figure{1..8}, figures/figureS{1..10}
│                          one directory per figure: src (its whole code) · tables (its inputs)
├── tables/                supplementary-table builders (Table S1-S12) and the SI assembly
├── analysis/              evidence tables individual figures re-use
├── envs/                  per-tool envs (pipeline) · as_run/ (reported runs) · analysis/ (figure code)
├── docs/                  read these first
└── _local/                NOT redistributed; see _local/README.md
```

## Reading order

| document | answers |
|---|---|
| [docs/data_availability.md](docs/data_availability.md) | exactly what is deposited, what is not, and why |
| [docs/sites_columns.md](docs/sites_columns.md) | column-by-column definition of a callset, incl. what each tool's `score` means |
| [docs/figure_index.md](docs/figure_index.md) | every figure → the scripts that draw it → the tables they read |
| [docs/reproducing.md](docs/reproducing.md) | path roots, environments, re-rendering a figure, rebuilding the callsets |
| [docs/pipeline.md](docs/pipeline.md) | the `src/harmonisation` stages and which table each one writes |
| [docs/tool_inventory_notes.md](docs/tool_inventory_notes.md) | how commands/versions were collected, and where the evidence is weakest |
| [docs/third_party_scripts.md](docs/third_party_scripts.md) | tool-owned scripts the Snakefile expects but does not ship |

## A. Running the detection pipeline

```bash
python setup.py                          # checks snakemake/conda, creates the directory skeleton
conda env create -f envs/cheui.yaml      # …one environment per tool you want to run
$EDITOR config/config.yaml               # samples, tools, reference paths, thresholds
snakemake --use-conda --cores 40 --dry-run
snakemake --use-conda --cores 40
```

Put raw pod5/fast5 under `raw_data/<sample>/` and the reference FASTA/GTF under
`reference/`. Dorado, Modkit and the tool binaries are installed separately;
[docs/third_party_scripts.md](docs/third_party_scripts.md) lists the tool-owned scripts
the rules call and where they ship. Every command line actually used for the reported
runs is in `tools/inventory/tables/command_lines.csv`, so a rule can be compared
against what really ran.

## B. Using the deposit

Read a callset - no installation:

```bash
head -3 data/callsets/RNA002/Human/HeLa_WT/m6A/CHEUI_m6A/HeLa_WT1.tsv
```

Re-render a figure from the deposited tables (each script resolves the repository root
on its own, so any working directory works):

```bash
conda env create -f envs/analysis/env-benchmark-revision.yml
conda run -n benchmark-revision --no-capture-output python figures/figure4/src/fig4_figures.py
bash scripts/run_figures.sh figure4
```

Rebuild the callsets themselves from raw tool output with
`bash src/harmonisation/scripts/run_all.sh`, which needs the private inputs listed in
[`_local/README.md`](_local/README.md).

## Checking a checkout

```bash
python scripts/verify_deposit.py
```

It asserts that no machine-specific absolute path is present, that no folder name from
the private working tree and no non-English working note survives in a text source,
that every Python/R/shell source parses, that each of the 406 callsets is listed in
`metadata/callsets_index.tsv` with a matching row count, and that every script named
in the figure index exists - then writes `metadata/deposited_files.sha256`, which the
second CI step compares against the committed copy. The same check runs on every push
(`.github/workflows/deposit-integrity.yml`).

## Licence and reuse

The pipeline code, the derived callsets and the analysis code are released under the MIT
License ([LICENSE](LICENSE)). The underlying sequencing data remain subject to the terms
of their original accessions, listed in `metadata/datasets.csv` with each study's
publication and DOI; third-party tools keep their own licences and are not redistributed
here. Please cite the source publication of the benchmark when reusing the callsets or
the evaluation definitions.
