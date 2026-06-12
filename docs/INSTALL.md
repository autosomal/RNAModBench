# RNAModBench — Installation Guide

## Contents

- [System requirements](#system-requirements)
- [Clone the repository](#clone-the-repository)
- [Install Conda environments](#install-conda-environments)
- [Install external binaries](#install-external-binaries)
- [R packages (visualisation)](#r-packages-visualisation)
- [Verify the installation](#verify-the-installation)
- [Known issues and workarounds](#known-issues-and-workarounds)

---

## System requirements

| Component | Minimum | Recommended |
|---|---|---|
| Operating system | Linux (any mainstream distribution) | Ubuntu 22.04 / CentOS 7+ |
| CPU | 8 cores | 32 cores or more |
| Memory | 16 GB | 64 GB+ (for `nanopolish eventalign`) |
| Disk | 100 GB | 500 GB+ (for FAST5 and intermediate files) |
| Python | 3.8+ | 3.10 |
| Conda | Any | Miniconda3 |
| GPU (optional) | NVIDIA V100 / A100 | ≥ 16 GB VRAM to accelerate CHEUI |

macOS / Windows users should use WSL2 or Docker.

---

## Clone the repository

```bash
git clone https://github.com/autosomal/RNAModBench.git
cd RNAModBench
```

---

## Install Conda environments

> Every tool in RNAModBench ships with a dedicated conda environment
> (`envs/*.yaml`). Snakemake activates the appropriate environment
> automatically as each rule runs.

### 0. Install Miniconda (if not already installed)

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
export PATH=$HOME/miniconda3/bin:$PATH
conda init bash
source ~/.bashrc
```

### 1. Configure the Bioconda / Conda-forge channels

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible
```

### 2. Install the main pipeline environment (with Snakemake)

```bash
conda create -n rnamodbench -c bioconda -c conda-forge \
    snakemake=7.32.4 minimap2 samtools nanopolish bedtools tabix \
    python=3.10 pandas=1.5 numpy=1.24 matplotlib=3.7 seaborn=0.12 \
    biopython=1.81
conda activate rnamodbench

# Additional PyPI packages
pip install matplotlib-venn
```

### 3. Install per-tool environments (as required)

Snakemake will `conda env create -f envs/xxx.yaml` automatically on
first invocation, but you can pre-install them to avoid network issues:

```bash
cd envs
for env in cheui m6anet nanocompore epinano dena mines xpore yanocomp \
           r eligos2 nanospa; do
  conda env create -f ${env}.yaml -n ${env}
done
cd ..
```

> If you only want a subset of tools (e.g. CHEUI + ELIGOS2), then only
> `envs/cheui.yaml` and `envs/eligos2.yaml` need to be installed.

---

## Install external binaries

The following tools are **not** distributed via conda and must be
installed separately.

### Guppy (basecaller)

Obtain from the Oxford Nanopore Technologies community portal (requires
an ONT account). Ensure the version matches your flow-cell type.

```bash
# Extract and put on PATH
tar -xzf ont-guppy_X.XX.X_linux64.tar.gz
export PATH=$(pwd)/ont-guppy/bin:$PATH

# Verify
guppy_basecaller --version
```

> ⚠ GPU versions require CUDA 11.2 or newer.

### R2Dtool (optional — used for coordinate liftover)

```bash
git clone https://github.com/bartongroup/R2Dtool.git
cd R2Dtool
cargo build --release
# Put ./target/release/r2d on PATH, or specify the absolute path in config.yaml
export PATH=$(pwd)/target/release:$PATH
```

### Tombo (required by DENA and MINES)

```bash
# Tombo is usually installed in its own conda environment
conda create -n tombo python=3.8
conda activate tombo
pip install ont-tombo
tombo --version
```

---

## R packages (visualisation)

Used by `scripts/create_guitar_plots.R` and
`scripts/generate_depth_plots.R`:

```r
# Run inside an R session:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Guitar")          # metagene profiles
install.packages("ggplot2")             # general plotting
install.packages("dplyr")
install.packages("readr")
```

---

## Verify the installation

Run the following commands in sequence. Each should print a version
number rather than an error:

```bash
conda activate rnamodbench
python        --version        # 3.8+
snakemake     --version        # 7.x
minimap2      --version        # 2.x
samtools      --version        # 1.x
nanopolish    --version        # 0.14+
guppy_basecaller --version     # 6.x / 7.x (skip if not installed)
Rscript       --version        # 4.x

# Python packages
python -c "import pandas, numpy, matplotlib, seaborn; print('OK')"

# R packages
Rscript -e 'library(Guitar); library(ggplot2); cat("OK\n")'
```

If everything passes, move on to [TUTORIAL.md](TUTORIAL.md) to run the
pipeline on a test dataset.

---

## Known issues and workarounds

| Symptom | Cause | Solution |
|---|---|---|
| `CondaHTTPError` or slow downloads | Institutional / academic network restricts conda access | Configure a mirror via `conda config --add channels`, or install `mamba` as a faster solver |
| `OSError: File not found: resources/R2Dtool/target/release/r2d` | `r2d` was not compiled | See "R2Dtool" section above |
| `std::bad_alloc` from nanopolish | eventalign parallelism too high | Reduce `nanopolish.threads` in `config.yaml` to 10 |
| CUDA out-of-memory in CHEUI | Large batch size exceeds GPU VRAM | Switch to CPU execution; reduce `cheui.threads` |
| Garbled characters in file paths | Python / system locale issue | Keep working-directory names in plain ASCII; `export LC_ALL=C.UTF-8` |
