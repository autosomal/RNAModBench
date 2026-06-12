# RNAModBench 安装指南

## 目录

- [系统要求](#系统要求)
- [克隆仓库](#克隆仓库)
- [安装 Conda 环境](#安装-conda-环境)
- [安装外部二进制工具](#安装外部二进制工具)
- [R 包（可视化用）](#r-包可视化用)
- [验证安装](#验证安装)
- [已知问题与 workaround](#已知问题与-workaround)

---

## 系统要求

| 组件 | 最低要求 | 推荐 |
|---|---|---|
| 操作系统 | Linux（任何主流发行版） | Ubuntu 22.04 / CentOS 7+ |
| CPU | 8 核 | 32 核或更多 |
| 内存 | 16 GB | 64 GB+（给 nanopolish eventalign） |
| 磁盘 | 100 GB | 500 GB+（存放 fast5 + 中间文件） |
| Python | 3.8+ | 3.10 |
| Conda | Any | Miniconda3 |
| GPU (可选) | NVIDIA V100 / A100 | 16GB 显存以上加速 CHEUI |

macOS / Windows 用户请使用 WSL2 或 Docker。

---

## 克隆仓库

```bash
git clone https://github.com/autosomal/RNAModBench.git
cd RNAModBench
```

---

## 安装 Conda 环境

> RNAModBench 的每个工具都有独立的 conda 环境配置（`envs/*.yaml`）。
> Snakemake 会在需要时自动激活对应环境。

### 0. 安装 Miniconda（如尚未安装）

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
export PATH=$HOME/miniconda3/bin:$PATH
conda init bash
source ~/.bashrc
```

### 1. 配置 Bioconda / Conda-forge channels

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible
```

### 2. 安装流水线主环境（含 Snakemake）

```bash
conda create -n rnamodbench -c bioconda -c conda-forge \
    snakemake=7.32.4 minimap2 samtools nanopolish bedtools tabix \
    python=3.10 pandas=1.5 numpy=1.24 matplotlib=3.7 seaborn=0.12 \
    biopython=1.81
conda activate rnamodbench

# 用 pip 安装 PyPI 上的额外包
pip install matplotlib-venn
```

### 3. 安装各工具的独立环境（按需要）

Snakemake 会在执行时自动 `conda env create -f envs/xxx.yaml`，
但你也可以预先装好，避免网络问题：

```bash
cd envs
for env in cheui m6anet nanocompore epinano dena mines xpore yanocomp \
           r eligos2 nanospa; do
  conda env create -f ${env}.yaml -n ${env}
done
cd ..
```

> 如果你只想跑子集工具（例如只跑 CHEUI + ELIGOS2），就只需要
> `cheui.yaml` 和 `eligos2.yaml`。

---

## 安装外部二进制工具

以下工具不在 Conda 中（或 Conda 的版本过旧），需手动安装：

### Guppy（basecaller）

从 Oxford Nanopore Technologies 官网下载对应的 Guppy CPU / GPU 版本：

```bash
# 解压并添加到 PATH
tar -xzf ont-guppy_X.XX.X_linux64.tar.gz
export PATH=$(pwd)/ont-guppy/bin:$PATH

# 验证
guppy_basecaller --version
```

> ⚠️ GPU 版本需 CUDA 11.2+。

### R2Dtool（可选 —— 用于坐标转换）

```bash
git clone https://github.com/chrisam/r2d-tool.git
cd r2d-tool
cargo build --release
# 把 ./target/release/r2d 放到 PATH 或在 config.yaml 指定路径
export PATH=$(pwd)/target/release:$PATH
```

### Tombo（DENA / MINES 需要）

```bash
# Tombo 通常通过 ONT 官网下载，并安装在单独的 conda 环境
conda create -n tombo python=3.8
conda activate tombo
pip install ont-tombo
tombo --version
```

---

## R 包（可视化用）

对于 `scripts/create_guitar_plots.R` 和 `scripts/generate_depth_plots.R`：

```r
# 在 R 控制台执行：
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Guitar")          # Guitar plots
install.packages("ggplot2")             # 通用画图
install.packages("dplyr")
install.packages("readr")
```

---

## 验证安装

依次执行以下命令，每个都应输出版本号而非报错：

```bash
conda activate rnamodbench
python        --version        # 3.8+
snakemake     --version        # 7.x
minimap2      --version        # 2.x
samtools      --version        # 1.x
nanopolish    --version        # 0.14+
guppy_basecaller --version     # 6.x / 7.x（可选，若未安装则跳过）
Rscript       --version        # 4.x

# Python 包
python -c "import pandas, numpy, matplotlib, seaborn; print('OK')"

# R 包
Rscript -e 'library(Guitar); library(ggplot2); cat("OK\n")'
```

若全部通过则安装完成。接下来按照 [TUTORIAL.md](TUTORIAL.md) 跑最小工作示例。

---

## 已知问题与 workaround

| 问题 | 原因 | 解决方案 |
|---|---|---|
| `CondaHTTPError` 或网络慢 | 公司 / 学术网络访问 conda 受限 | 使用 `conda config --add channels` 配置镜像站；或改用 `mamba` 作为 solver |
| `OSError: File not found: resources/R2Dtool/target/release/r2d` | `r2d` 未编译 | 见上文 "R2Dtool" 小节 |
| nanopolish `std::bad_alloc` | eventalign 默认并行度过高 | 调低 `config.yaml` 的 `nanopolish.threads` 到 10 |
| CHEUI `CUDA out of memory` | 大批次 → GPU 显存不足 | 关闭 GPU 模式；或在 `cheui.threads` 调小 |
| 中文路径出现乱码 | Python/系统 locale 问题 | 保持工作目录全英文；`export LC_ALL=C.UTF-8` |
