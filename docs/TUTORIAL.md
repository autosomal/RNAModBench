# RNAModBench — End-to-End Tutorial

This document walks through running the RNAModBench pipeline on a real
dataset. It assumes the dependencies described in
[INSTALL.md](INSTALL.md) have been installed.

## Contents

- [Prepare reference files](#prepare-reference-files)
- [Prepare raw data (FAST5)](#prepare-raw-data-fast5)
- [Edit config.yaml](#edit-configyaml)
- [Dry-run the workflow](#dry-run-the-workflow)
- [Run the full pipeline](#run-the-full-pipeline)
- [Run only a subset of tools](#run-only-a-subset-of-tools)
- [Inspect outputs](#inspect-outputs)

---

## Prepare reference files

We use human GRCh38 / Ensembl 112 as the default reference assembly.
Substitute the FTP paths below for your species of interest.

```bash
mkdir -p reference
cd reference

# 1) Genome FASTA (~3 GB)
wget -O genome.fa.gz \
  ftp://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip genome.fa.gz
samtools faidx genome.fa        # produces genome.fa.fai

# 2) Transcriptome FASTA
gffread -w transcriptome.fa -g genome.fa Homo_sapiens.GRCh38.112.gtf

# 3) Gene model (GTF)
wget -O genes.gtf.gz \
  ftp://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
gunzip genes.gtf.gz

cd ..
```

> **Chromosome-name consistency.** The `>` lines of `genome.fa` and
> column 1 of `genes.gtf` must match exactly. If one uses `chr1` and
> the other uses `1`, reformat with `sed` before running the pipeline.

---

## Prepare raw data (FAST5)

### Option A — Public example data (small, fast)

A quick way to validate the pipeline end-to-end is to download a small
public FAST5 dataset. A few megabases of coverage per sample are enough
to exercise every rule:

```bash
mkdir -p data/treatment_rep1
# Download a public FAST5 file (e.g. from the ENA under accession
# SRRNNNNNN or similar) into data/treatment_rep1/fast5/.
```

### Option B — Use your own data

Place your FAST5 files under sample-specific directories as follows:

```
data/
├── treatment_rep1/fast5/*.fast5
├── treatment_rep2/fast5/*.fast5
├── control_rep1/fast5/*.fast5
└── control_rep2/fast5/*.fast5
```

> Directory names must match the entries in `config.yaml`
> exactly (case-sensitive).

---

## Edit config.yaml

Open `config/config.yaml` and adjust at least the first two blocks:

```yaml
samples:
  - treatment_rep1
  - control_rep1

tools:
  - CHEUI
  - ELIGOS2
  - m6Anet

data_dir:      "data"
reference_dir: "reference"
results_dir:   "results"
```

Refer to [CONFIG_REFERENCE.md](CONFIG_REFERENCE.md) for the full
parameter catalogue. Leave per-tool thresholds at their defaults for
an initial run.

---

## Dry-run the workflow

Always dry-run before launching a full compute. This validates the
rule-DAG and confirms that every required input file is present:

```bash
snakemake --dry-run --printshellcmds --cores 1
```

If you only want to see which rules fire, add `-n` for the summary:

```bash
snakemake -n --cores 1
```

The summary at the bottom of the output lists the per-rule counts
that will execute.

---

## Run the full pipeline

```bash
# Recommended: constrain concurrency to avoid nanopolish eventalign
# consuming all available RAM.
snakemake --use-conda --cores 24 --resources mem_mb=64000

# If no GPU is available, reduce cheui.threads or disable CHEUI in
# the tools list.
```

---

## Run only a subset of tools

```bash
# Run only CHEUI + ELIGOS2 (by limiting the tools list in config.yaml)
# Or, alternatively, tell Snakemake to build only the target files
# for specific tools:
snakemake --use-conda --cores 24 \
  results/CHEUI/treatment_rep1/treatment_rep1_CHEUI_processed.txt \
  results/ELIGOS2_solo/treatment_rep1/treatment_rep1_ELIGOS2_solo_processed.txt
```

---

## Inspect outputs

A successful pipeline run produces a directory layout like this:

```
results/
├── CHEUI/treatment_rep1/treatment_rep1_CHEUI_processed.txt
├── ELIGOS2_solo/treatment_rep1/treatment_rep1_ELIGOS2_solo_processed.txt
├── m6Anet/treatment_rep1/treatment_rep1_m6Anet_processed.txt
│
├── summary/
│   ├── modification_summary.tsv     # per-tool site counts
│   ├── tool_comparison.tsv          # site-level concordance matrix
│   ├── tool_overlap.png             # pairwise Venn diagram (if ≥2 tools)
│   ├── modification_distribution.png
│   ├── guitar_plots_mrna.png        # mRNA metagene profile
│   └── liftover_summary.tsv         # (if liftover is enabled)
│
├── report/
│   └── RNAModBench_report.html      # Bootstrap-style HTML report
│
└── qc/
    └── depth_plots.png              # per-sample coverage
```

### Suggested first steps for interpretation

1. Open `results/report/RNAModBench_report.html` for a visual overview.
2. Check `modification_summary.tsv` for the number of sites per tool;
   suspiciously high or low counts warrant revisiting thresholds.
3. For genome-browser visualisation, convert the processed files to
   BED6/bedGraph using standard Unix tools.
4. For targeted deep-dive on specific genes or transcripts, invoke
   `scripts/r2d_liftover.py` once to project calls from transcriptome
   coordinates back to genome space.

> **Troubleshooting hint.** If any rule fails, look at the log
> directory that Snakemake writes next to the rule, and consult the
> troubleshooting section of the README.
