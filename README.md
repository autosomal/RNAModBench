# RNAModBench — A Comprehensive Benchmark Pipeline for Nanopore-Based RNA Modification Detection

[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Pipeline-2D72B8.svg?logo=data:image/svg%2Bxml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciPjwvc3ZnPg==)]()
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue.svg)]()
[![Snakemake 7+](https://img.shields.io/badge/Snakemake-7%2B-brightgreen.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)]()

---

## Abstract

RNAModBench is an end-to-end, species-agnostic pipeline for the detection of RNA modifications from direct RNA sequencing (DRS) data produced on Oxford Nanopore Technologies (ONT) platforms. The pipeline integrates **12 state-of-the-art modification-calling methods**, spanning both **signal-level** (current-intensity) and **alignment-level** (base-quality/mismatch) approaches, and harmonises their heterogeneous outputs into a single, standardised, BED-compatible tabular format. Downstream modules implement consensus-calling, inter-tool concordance analysis, motif enrichment, metagene profiling, and automated reporting. The pipeline is implemented in [Snakemake](https://snakemake.github.io/) and requires only a sample directory of FAST5 files and standard reference resources (genome FASTA, transcriptome FASTA, and a GTF gene model).

**Pipeline diagram:** see [docs/FLOWCHART.md](docs/FLOWCHART.md) for publication-quality Mermaid and Graphviz flowcharts.

---

## Table of Contents

1. [Features](#1-features)
2. [Installation](#2-installation)
3. [Reference Files](#3-reference-files)
4. [Input Data Structure](#4-input-data-structure)
5. [Running the Pipeline](#5-running-the-pipeline)
6. [Tool Overview & Data Requirements](#6-tool-overview--data-requirements)
7. [Output Format Specification](#7-output-format-specification)
8. [Interpreting Results](#8-interpreting-results)
9. [Computational Resource Estimates](#9-computational-resource-estimates)
10. [Quality Control](#10-quality-control)
11. [Troubleshooting](#11-troubleshooting)
12. [Documentation Roadmap](#12-documentation-roadmap)
13. [Citing the Integrated Tools](#13-citing-the-integrated-tools)
14. [Contributing](#14-contributing)
15. [License & Support](#15-license--support)

---

## 1. Features

- **Multi-tool integration** — 12 different RNA-modification detection methods, spanning signal-level and alignment-level approaches.
- **Single-sample and contrast modes** — tools accept both single-sample inputs and paired (treatment vs control) experimental designs.
- **Species-agnostic** — works with any genome/transcriptome provided compatible reference files are supplied.
- **Standardised output** — all 12 tools are harmonised into a single 7-column TSV format with consistent coordinate semantics (BED 0-based).
- **Consensus calling** — automated derivation of high-confidence modification sites supported by ≥3 tools.
- **Comprehensive reporting** — HTML report with embedded per-tool QC, Venn diagrams of tool overlap, k-mer motif analysis, and metagene profiles.
- **Reproducible execution** — Snakemake workflow with deterministic rule DAG; all environments defined in `envs/*.yaml`.
- **Modular design** — individual tools (or subsets) can be executed independently.

---

## 2. Installation

### 2.1 Prerequisites

| Software | Version | Notes |
|---|---|---|
| Python | ≥ 3.8 | Core interpreter for all tool wrappers |
| Snakemake | ≥ 7.0 | Workflow orchestration |
| Conda / Miniconda | latest | Environment management |
| ONT Guppy | ≥ 6.4 | Base-calling (GPU version recommended) |
| minimap2 | ≥ 2.24 | Transcriptome / genome alignment |
| samtools | ≥ 1.15 | BAM processing & indexing |
| nanopolish | ≥ 0.14 | Event-level signal alignment |
| R | ≥ 4.2 | For Guitar plots & depth plots |
| R2Dtool | latest | Transcriptome → genome coordinate liftover |

### 2.2 Conda environment setup

```bash
# Clone the repository
git clone https://github.com/autosomal/RNAModBench.git
cd RNAModBench

# Create the primary Snakemake environment
conda env create -f envs/cheui.yaml
conda env create -f envs/m6anet.yaml
conda env create -f envs/nanocompore.yaml
conda env create -f envs/eligos2.yaml
conda env create -f envs/epinano.yaml
conda env create -f envs/drummer.yaml
conda env create -f envs/dena.yaml
conda env create -f envs/mines.yaml
conda env create -f envs/xpore.yaml
conda env create -f envs/yanocomp.yaml
conda env create -f envs/nanospa.yaml
conda env create -f envs/r.yaml
```

Snakemake automatically activates the appropriate environment for each rule. **You do not need to manually activate any environment before launching the workflow.**

### 2.3 External binaries

The following tools are *not* distributed via conda and must be installed separately:

- **Guppy** — obtain from the [ONT community portal](https://community.nanoporetech.com/downloads). Ensure `guppy_basecaller` is on `PATH`.
- **R2Dtool** — used for coordinate liftover. Compile with:
  ```bash
  git clone https://github.com/bartongroup/R2Dtool.git
  cd R2Dtool && cargo build --release
  echo "export PATH=\$PATH:$(pwd)/target/release" >> ~/.bashrc
  ```

Detailed, step-by-step installation instructions are available in [docs/INSTALL.md](docs/INSTALL.md).

---

## 3. Reference Files

Four reference files are required. **Chromosome identifiers must be identical** across all files.

| File | Format | Purpose |
|---|---|---|
| `genome.fa` | Multi-FASTA | Genome reference for alignment-based tools; must be bgzip-decompressed (plain text). |
| `transcriptome.fa` | Multi-FASTA | Spliced transcript sequences; headers must match `transcript_id` values in the GTF. |
| `genes.gtf` | Ensembl / Gencode GTF | Gene model — required for Guitar plots, R2Dtool liftover, and feature-aware tools. |
| `genes.bed` | 12-column BED | Optional; derived from the GTF for ELIGOS2 region-based calling. |

> **Important:** `genome.fa` must be indexed before execution (`samtools faidx genome.fa`). The pipeline also expects a matching `genome.fa.fai`.

### 3.1 Example: Human (GRCh38 / Ensembl 112)

```bash
mkdir -p reference && cd reference
wget ftp://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa
samtools faidx genome.fa

wget ftp://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
gunzip Homo_sapiens.GRCh38.112.gtf.gz
mv Homo_sapiens.GRCh38.112.gtf genes.gtf

# Transcriptome — either download pre-built Ensembl cDNA or generate from the GTF:
gffread -w transcriptome.fa -g genome.fa genes.gtf
cd ..
```

For non-human species, substitute the FTP path with the appropriate species directory.

---

## 4. Input Data Structure

Place each sample's FAST5 files (multi- or single-read format) in a subdirectory under `data/`. The directory name *is* the sample identifier. Ensure it contains only letters, digits, and underscores.

```
RNAModBench/
├── config/
│   └── config.yaml
├── data/
│   ├── treatment_rep1/
│   │   └── fast5/
│   │       ├── batch0_0.fast5
│   │       ├── batch0_1.fast5
│   │       └── ...
│   ├── treatment_rep2/
│   │   └── fast5/
│   ├── control_rep1/
│   │   └── fast5/
│   └── control_rep2/
│       └── fast5/
├── reference/
│   ├── genome.fa
│   ├── genome.fa.fai
│   ├── transcriptome.fa
│   └── genes.gtf
└── Snakefile
```

Contrast tools (Nanocompore, xPore, yanocomp, Epinano_DiffErr, DRUMMER) require at least **one treatment and one control sample** and will pair them automatically (order in `config.yaml` matters).

---

## 5. Running the Pipeline

### 5.1 Configure

Edit `config/config.yaml` to specify:
- the list of samples (directory names under `data/`)
- the subset of tools to run
- paths to reference files
- per-tool thresholds

Refer to [docs/CONFIG_REFERENCE.md](docs/CONFIG_REFERENCE.md) for the complete parameter reference.

### 5.2 Dry-run first

Always validate the workflow before computing:

```bash
snakemake --dry-run --printshellcmds --cores 1
# Inspect DAG visually
snakemake --dag | dot -Tsvg > docs/figures/dag.svg
```

### 5.3 Full pipeline execution

```bash
# 40-core server, 256 GB RAM (typical for event-align + signal tools)
snakemake --use-conda --cores 40 --resources mem_mb=200000 --keep-going
```

### 5.4 Running subsets of tools

```bash
# Signal tools only (require event-align)
snakemake --use-conda --cores 40 \
    results/CHEUI/sample1/sample1_CHEUI_processed.txt \
    results/m6Anet/sample1/sample1_m6Anet_processed.txt

# Alignment tools only (faster, no FAST5 needed after base-calling)
snakemake --use-conda --cores 20 \
    results/ELIGOS2_solo/sample1/sample1_ELIGOS2_solo_processed.txt \
    results/Epinano/sample1/sample1_Epinano_processed.txt

# Contrast tools (pair treatment vs control)
snakemake --use-conda --cores 20 \
    results/Nanocompore/treatment_vs_control/treatment_vs_control_Nanocompore_processed.txt

# Generate only the summary / HTML report (skips re-calling tools)
snakemake --use-conda --cores 4 results/report/RNAModBench_report.html
```

See [docs/TUTORIAL.md](docs/TUTORIAL.md) for a step-by-step tutorial with a test dataset.

---

## 6. Tool Overview & Data Requirements

| Tool | Type | Coordinates | Requires `nanopolish eventalign` | Requires control sample | Typical runtime¹ |
|---|---|---|:---:|:---:|---|
| **CHEUI** | Signal, single-sample | Transcriptome | ✅ | ❌ | 30–60 min |
| **m6Anet** | Signal, single-sample | Transcriptome | ✅ | ❌ | 15–30 min |
| **DENA** | Signal, single-sample | Transcriptome | ✗² | ❌ | 20–40 min |
| **MINES** | Signal, single-sample | Transcriptome | ✗² | ❌ | 15–30 min |
| **Nanocompore** | Signal, contrast | Transcriptome | ✅ | ✅ | 20–40 min |
| **xPore** | Signal, contrast | Transcriptome | ✅ | ✅ | 15–30 min |
| **yanocomp** | Signal, contrast | Transcriptome | ✅ | ✅ | 10–30 min |
| **ELIGOS2** | Alignment, single-sample | Genome | ✗ | ❌ | 10–20 min |
| **Epinano** | Alignment, single-sample | Genome | ✗ | ❌ | 5–15 min |
| **Epinano_DiffErr** | Alignment, contrast | Genome | ✗ | ✅ | 10–20 min |
| **DRUMMER** | Alignment, contrast | Genome | ✗ | ✅ | 10–15 min |
| **NanoSPA** | Alignment, single-sample | Genome | ✗ | ❌ | 10–30 min |

> ¹ Per sample, ~1M reads, 40 cores.  
> ² DENA/MINES use **Tombo re-squiggle** instead of `nanopolish eventalign`.

Full technical details are provided in [docs/TOOLS_OVERVIEW.md](docs/TOOLS_OVERVIEW.md).

---

## 7. Output Format Specification

Every tool produces one output file per sample (or per contrast pair):

```
results/<tool>/<sample>/<sample>_<tool>_processed.txt
```

All files share a **standard 7-column tab-separated format** with a header row:

| Column | Header | Type | Description |
|---|---|---|---|
| 1 | `Chr` | `string` | Chromosome or transcript identifier. **Transcriptome-coordinate tools** populate this column with transcript IDs (e.g., `ENST00000367770`); **genome-coordinate tools** populate it with chromosome names (e.g., `chr1`). |
| 2 | `Start` | `int` | **BED 0-based** start coordinate of the modification site. |
| 3 | `End` | `int` | **BED 0-based, exclusive** end coordinate. Equals `Start + 1` for single-base modifications. |
| 4 | `Status` | `string` | Either `Mod` (site passes significance / probability filter) or `Unmod`. |
| 5 | `Prob` | `float` | Probability, score, or adjusted P-value. Semantics vary by tool — see §6. |
| 6 | `Strand` | `{+ , − , *}` | Genomic strand. `*` is used when strand is undetermined by a particular tool. |
| 7 | `mod_ratio` | `float` | Estimated stoichiometry (fraction of reads supporting the modification), in the range [0, 1]. |

A complete specification, including column semantics per tool, is provided in [docs/OUTPUT_FORMAT.md](docs/OUTPUT_FORMAT.md).

---

## 8. Interpreting Results

We recommend the following exploratory workflow after running RNAModBench:

1. **Per-tool modification counts** — inspect `results/summary/modification_summary.tsv`. Tools reporting orders-of-magnitude more (or fewer) sites than expected should prompt a review of filter thresholds.  
2. **Inter-tool concordance** — tools agreeing on a site substantially increase confidence. Filter on `Tool_Count ≥ 3` from `modification_summary.tsv` for a high-confidence call-set.  
3. **Metagene profile** — `results/summary/guitar_plots_mrna.png` shows the positional distribution of modifications along mature mRNA (5′UTR → CDS → 3′UTR). m⁶A is typically enriched near the stop codon.  
4. **Motif analysis** — run `scripts/extract_5mer.py` on any output TSV to obtain 5-mer context; enrichment of the canonical `RRACH` motif for m⁶A is a sanity check for signal quality.  
5. **Liftover** — for transcriptome-coordinate tools, run `scripts/r2d_liftover.py` to project coordinates to genome-space for compatibility with browser tracks and genome-wide analyses.

---

## 9. Computational Resource Estimates

Benchmarks are for **~1 million reads / ~10 GB of FAST5 per sample** on a 40-core Intel Xeon, 256 GB RAM, with a single NVIDIA A100 GPU available:

| Stage | CPU threads | RAM (GB) | Wall-clock time |
|---|:---:|:---:|---|
| Guppy basecalling (CPU) | 40 | 16 | 30–60 min |
| Guppy basecalling (GPU) | 20 | 8 GPU + 16 CPU | 5–10 min |
| minimap2 (transcriptome + genome) | 40 | 8 | <5 min |
| nanopolish eventalign | 40 | 64 | **2–4 h** |
| CHEUI inference | 40 | 32 | 30–60 min |
| m6Anet inference | 40 | 16 | 15–30 min |
| ELIGOS2 + Epinano + NanoSPA | 12 each | 8 | 20–40 min combined |
| Nanocompore / xPore / yanocomp (per contrast) | 20 | 16 | 20–40 min |
| Summary & report generation | 4 | 4 | <5 min |

**`nanopolish eventalign` is the compute bottleneck.** Consider:
- down-sampling very large samples;
- substituting with `tombo resquiggle` if signal tools DENA/MINES suffice;
- using GPU-accelerated alternatives where available.

---

## 10. Quality Control

RNAModBench generates several QC artefacts automatically:

- **Depth-of-coverage plots** (`results/qc/*_depth.png`) generated from BAM files, enabling inspection of coverage biases per gene/transcript.
- **Per-tool modification-call histograms** embedded in the HTML report.
- **k-mer composition tables** produced by `extract_5mer.py`, which reveal whether called sites cluster at canonical motifs (e.g., `DRACH` / `RRACH` for m⁶A).

---

## 11. Troubleshooting

| Error / symptom | Likely cause | Remedy |
|---|---|---|
| `MissingInputException` / Missing input files for rule | Sample directory name in `config.yaml` does not match filesystem, or reference files are absent | Double-check `samples:` list; confirm `reference/genome.fa`, `genes.gtf` exist and are bgzip-decompressed |
| `MemoryError` / `Killed` | nanopolish eventalign or CHEUI exceeding available RAM | Reduce `nanopolish.threads`; add `--resources mem_mb=32000` |
| `[E::fai_build_core] Truncated file` | `genome.fa` not fully gunzip'd, or corrupted FASTA | Re-download / re-index with `samtools faidx` |
| `OSError: File not found … R2Dtool/target/release/r2d` | `r2d` binary missing or not on `PATH` | Compile with `cargo build --release`; add to `PATH` |
| Per-tool `*_processed.txt` is empty | Filter thresholds too strict, or no modification-candidate sites | Relax `prob_threshold` / `pvalue_threshold` / `ratio_threshold` in `config.yaml`; check base-call quality |
| CUDA out-of-memory | GPU card lacks VRAM for CHEUI / deep-learning tools | Switch to CPU execution; reduce batch size |

Additional advice and a worked example are in [docs/TUTORIAL.md](docs/TUTORIAL.md) and [docs/FAQ.md](docs/FAQ.md).

---

## 12. Documentation Roadmap

| Document | Content |
|---|---|
| [docs/INSTALL.md](docs/INSTALL.md) | Step-by-step installation, conda environments, external binaries |
| [docs/TUTORIAL.md](docs/TUTORIAL.md) | End-to-end walkthrough with a real dataset |
| [docs/TOOLS_OVERVIEW.md](docs/TOOLS_OVERVIEW.md) | Per-tool input requirements, output semantics, run-time estimates |
| [docs/OUTPUT_FORMAT.md](docs/OUTPUT_FORMAT.md) | Full column-by-column specification + example records |
| [docs/CONFIG_REFERENCE.md](docs/CONFIG_REFERENCE.md) | `config/config.yaml` parameter reference |
| [docs/FLOWCHART.md](docs/FLOWCHART.md) | Publication-quality Mermaid + Graphviz pipeline diagrams |
| [docs/FAQ.md](docs/FAQ.md) | Frequently-asked questions |
| [`CHANGELOG.md`](CHANGELOG.md) | Version history |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | Bug reports, pull-request guidelines, code style |

---

## 13. Citing the Integrated Tools

If you use RNAModBench in published work, **please cite the individual tools on which calls depend** rather than RNAModBench itself:

| Tool | Reference |
|---|---|
| CHEUI | Liu et al. *Nature Biotechnology* (2023) |
| m6Anet | Hendra et al. *Nature Methods* (2022) |
| Nanocompore | Leger et al. *Nature Communications* (2021) |
| DENA | Li et al. *Genome Biology* (2022) |
| Epinano | Liu et al. *Nature Communications* (2019) |
| MINES | Begik et al. *Nucleic Acids Research* (2019) |
| xPore | Pratanwanich et al. *Nature Biotechnology* (2021) |
| yanocomp | Ma et al. *Genome Biology* (2022) |
| ELIGOS2 | Jenjaroenpun et al. *Genome Biology* (2021) |
| DRUMMER | Ortiz et al. *Nucleic Acids Research* (2021) |
| NanoSPA | Dong et al. *Cell Genomics* (2023) |
| R2Dtool | Barton group (GitHub, ongoing) |
| Guitar (R package) | Li et al. *Bioinformatics* (2021) |

For tool-specific DOIs, consult the original publications.

---

## 14. Contributing

We welcome bug reports, feature requests, and pull requests. Before contributing, please read [`CONTRIBUTING.md`](CONTRIBUTING.md), which outlines code style, commit conventions, and the procedure for adding a new detection tool to the pipeline.

---

## 15. License & Support

RNAModBench is released under the **MIT License**. See `LICENSE` for the full text.

For bug reports, questions, or feature requests, please open an issue on the [GitHub repository](https://github.com/autosomal/RNAModBench).
