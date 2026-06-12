# Tool-Level Overview — RNAModBench

This document summarises the **12 RNA-modification detection tools** integrated into RNAModBench, their input requirements, the coordinate systems they operate on, and their run-time profiles. It is intended both as a reference when planning experiments (e.g., which tools can be run given available data) and as a quick technical index for method comparison.

---

## Contents

1. [Complete tool matrix](#1-complete-tool-matrix)
2. [Signal-level vs alignment-level tools](#2-signal-level-vs-alignment-level-tools)
3. [Single-sample vs contrast (paired) tools](#3-single-sample-vs-contrast-paired-tools)
4. [Coordinate systems](#4-coordinate-systems)
5. [Typical run times and memory usage](#5-typical-run-times-and-memory-usage)
6. [Tool-specific reference and model files](#6-tool-specific-reference-and-model-files)
7. [How to select tools for your experiment](#7-how-to-select-tools-for-your-experiment)

---

## 1. Complete tool matrix

| Tool | Primary modification | Input type | Coordinates | Requires `nanopolish eventalign` | Requires a paired control | Conda environment |
|---|---|---|---|:---:|:---:|---|
| CHEUI | m⁶A | Signal-level (raw current) | Transcriptome | ✅ | ❌ | `envs/cheui.yaml` |
| m6Anet | m⁶A | Signal-level (per-site likelihood) | Transcriptome | ✅ | ❌ | `envs/m6anet.yaml` |
| DENA | m⁶A (motif-restricted) | Signal-level (Tombo re-squiggle) | Transcriptome | ❌¹ | ❌ | `envs/dena.yaml` |
| MINES | m⁶A (machine-learning model) | Signal-level (Tombo fraction-modified) | Transcriptome | ❌¹ | ❌ | `envs/mines.yaml` |
| Nanocompore | Any significant signal shift | Signal-level (per-sample k-mer model) | Transcriptome | ✅ | ✅ | `envs/nanocompore.yaml` |
| xPore | Differential m⁶A / m⁶Am | Signal-level (Gaussian mixture model) | Transcriptome | ✅ | ✅ | `envs/xpore.yaml` |
| yanocomp | Any significant signal shift | Signal-level (GMM) | Transcriptome | ✅ | ✅ | `envs/yanocomp.yaml` |
| ELIGOS2 | Any (statistics-based) | Alignment-level (base-mismatch ratio) | Genome | ❌ | ❌ | `envs/eligos2.yaml` |
| Epinano | m⁶A (SVM-based) | Alignment-level (quality + indel) | Genome | ❌ | ❌ | `envs/epinano.yaml` |
| Epinano_DiffErr | Differential m⁶A (treatment vs control) | Alignment-level (delta error) | Genome | ❌ | ✅ | `envs/epinano.yaml` |
| DRUMMER | Differential modification (coverage-based) | Alignment-level | Genome | ❌ | ✅ | `envs/drummer.yaml` |
| NanoSPA | m⁶A, Ψ (pseudo-uridine), multi-mod | Alignment-level (Bayesian) | Genome | ❌ | ❌ | `envs/nanospa.yaml` |

> ¹ DENA and MINES rely on **Tombo re-squiggle** rather than `nanopolish eventalign`. Tombo is invoked via a separate Snakemake rule and does not require a separate conda environment.

---

## 2. Signal-level vs alignment-level tools

### Signal-level tools

Signal-level tools model the distribution of raw ionic current at each k-mer in each read. They therefore require:
1. FAST5 files (raw signal) from the sequencer;
2. A **per-read event alignment** produced by `nanopolish eventalign` (or Tombo `resquiggle` for DENA/MINES);
3. A transcriptome FASTA to anchor the events to known transcripts.

Advantages:
- More sensitive for m⁶A and m⁵C-like modifications.
- Some tools (m6Anet, CHEUI) can produce per-site stoichiometry estimates.

Disadvantages:
- Slow — event alignment takes 2–4 h per million reads.
- Requires FAST5 (which is increasingly phased out by ONT in favour of POD5).

**Tools:** CHEUI, m6Anet, DENA, MINES, Nanocompore, xPore, yanocomp.

### Alignment-level tools

Alignment-level tools operate on base-called reads aligned to the *genome*. They detect modification-associated changes in base-call quality, per-site mismatch rate, or local coverage. They therefore require:
1. A base-called FASTQ from Guppy / Dorado;
2. A genome-aligned BAM;
3. A reference genome FASTA.

Advantages:
- Fast — minutes rather than hours.
- Do not require FAST5 / raw signal.

Disadvantages:
- Cannot distinguish certain modifications from true sequence variants.
- Prone to false positives at low-coverage sites.

**Tools:** ELIGOS2, Epinano, Epinano_DiffErr, DRUMMER, NanoSPA.

---

## 3. Single-sample vs contrast (paired) tools

### Single-sample tools

| Tool | Statistical approach |
|---|---|
| CHEUI | Probabilistic (HMM + deep-learning classifier) |
| m6Anet | Probabilistic (site-level likelihood) |
| DENA | LSTM classifier on k-mer window |
| MINES | Ensemble ML on Tombo fraction-modified signal |
| ELIGOS2 | Logistic regression + rate-ratio test |
| Epinano | Support-vector machine (SVM) on quality/mismatch/indel |
| NanoSPA | Bayesian hierarchical model |

### Contrast (paired) tools

These tools compare a treated and an untreated sample and report sites whose modification state differs significantly. They are **not** designed to call modifications in a single sample.

| Tool | Statistical approach |
|---|---|
| Nanocompore | GMM-based log-odds ratio + KS-test |
| xPore | Gaussian mixture model (D.M. probability) |
| yanocomp | GMM test |
| Epinano_DiffErr | Fisher's exact test on per-site error delta |
| DRUMMER | Coverage-based rate-ratio test |

> RNAModBench automatically pairs samples in the order they appear in `config.yaml` (`samples[0]` vs `samples[1]`, `samples[2]` vs `samples[3]`, …). For experiments with unbalanced numbers of treatments and controls, run the contrast tools manually or use a wrapper script.

---

## 4. Coordinate systems

RNAModBench distinguishes two coordinate systems and keeps calls from each system separate until the summary stage. Downstream `r2d_liftover` can project transcriptome-coordinate results to genome-space.

### Transcriptome-coordinate tools
`Chr` column contains a **transcript ID** (e.g., `ENST00000367770`); `Start`/`End` are 0-based offsets from the 5′-end of the transcript.

*Tools:* CHEUI · m6Anet · DENA · MINES · Nanocompore · xPore · yanocomp.

### Genome-coordinate tools
`Chr` column contains a **chromosome name** (e.g., `chr1`). `Start`/`End` are BED 0-based genomic coordinates.

*Tools:* ELIGOS2 · Epinano · Epinano_DiffErr · DRUMMER · NanoSPA.

---

## 5. Typical run times and memory usage

Benchmarked on a 40-core Intel Xeon server, 256 GB RAM, ~1M reads / ~10 GB of FAST5:

| Stage | Cores | RAM (GB) | Wall time |
|---|---|---|---|
| Guppy basecalling (CPU) | 40 | 16 | 30–60 min |
| Guppy basecalling (GPU) | 20 | 8 GPU + 16 CPU | 5–10 min |
| minimap2 (transcriptome + genome) | 40 | 8 | <5 min |
| nanopolish eventalign | 40 | **64** | **2–4 h** |
| CHEUI | 40 | 32 | 30–60 min |
| m6Anet | 40 | 16 | 15–30 min |
| DENA | 40 | 16 | 20–40 min |
| MINES | 20 | 8 | 15–30 min |
| ELIGOS2 + Epinano + NanoSPA | 12 each | 8 | 20–40 min combined |
| Nanocompore / xPore / yanocomp (per contrast) | 20 | 16 | 20–40 min |
| Summary and reporting | 4 | 4 | <5 min |

**Take-home messages:**
- `nanopolish eventalign` is the pipeline bottleneck for signal-level tools.
- Alignment-only tools are orders-of-magnitude faster and suitable for exploratory / large-scale analysis.
- GPU acceleration helps only Guppy and CHEUI; the remaining tools are CPU-bound.

---

## 6. Tool-specific reference and model files

Beyond the standard reference files (genome FASTA, transcriptome FASTA, genes GTF), the following tools require additional per-site models, training data, or binary executables:

| Tool | Additional resource | Configuration field |
|---|---|---|
| CHEUI | HDF5 model files (`model_1.h5`, `model_2.h5`), k-mer model CSV | `cheui.model1`, `cheui.model2`, `cheui.kmer_model` |
| Epinano | Trained SVM model (linear), per-feature error columns | `epinano.model`, `epinano.columns` |
| DENA | Trained LSTM model, `RRACH` k-mer mask | `dena.model`, `dena.motif` |
| MINES | Model weights & names file | `mines.kmer_models` |
| NanoSPA | Per-mod base-model files | Defaults bundled |
| R2Dtool | Compiled `r2d` binary | `utilities.r2d_tool` |

Each path is specified in `config/config.yaml` and described in full in [docs/CONFIG_REFERENCE.md](CONFIG_REFERENCE.md).

---

## 7. How to select tools for your experiment

As a rule of thumb, RNAModBench users tend to select tools according to the following decision tree:

1. **Do you have FAST5?**  
   ✅ Yes → run all signal-level tools (CHEUI, m6Anet, DENA, MINES) + all alignment-level tools. Take the intersection (`Tool_Count ≥ 3` in the summary file) as the high-confidence call-set.  
   ❌ No → run the alignment-level tools only (ELIGOS2, Epinano, NanoSPA).

2. **Do you have paired treatment / control samples?**  
   ✅ Yes → additionally run the contrast tools (Nanocompore, xPore, yanocomp, Epinano_DiffErr, DRUMMER) to obtain sites whose modification state changes between conditions.  
   ❌ No → skip contrast tools.

3. **Are you targeting m⁶A specifically?**  
   If so, CHEUI / m6Anet / DENA / Epinano are the most sensitive tools. If you are interested in Ψ or multi-mod calling, add NanoSPA (alignment-level) and/or Nanocompore (signal-level, any modification).

4. **Compute budget?**  
   If resources are severely limited, start with ELIGOS2 and Epinano (fast, alignment-level) to verify modification signal exists, then add signal-level tools incrementally.
