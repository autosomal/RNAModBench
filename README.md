# RNAModBench - Comprehensive Nanopore RNA Modification Detection Pipeline

## Overview

RNAModBench is a comprehensive, cross-species pipeline for detecting RNA modifications from nanopore direct RNA sequencing data. This pipeline integrates multiple state-of-the-art tools and provides a standardized workflow for comparative analysis.

## Features

- **Multi-tool Integration**: Supports 11 different RNA modification detection tools
- **Cross-species Compatibility**: Works with any species with available reference genome
- **Standardized Output**: Converts all tool outputs to unified BED-like format
- **Comparative Analysis**: Generates tool overlap statistics and consensus calls
- **Quality Control**: Includes comprehensive QC metrics and filtering
- **Scalable**: Supports high-throughput analysis with parallel processing
- **Modular Design**: Run individual tools or complete pipeline

## Installation

### Prerequisites

- Python 3.8+
- Conda or Miniconda
- Snakemake
- Git

### Quick Setup

```bash
# Clone the repository
git clone https://github.com/autosomal/RNAModBench.git
cd RNAModBench

# Create conda environments
conda env create -f envs/cheui.yaml
conda env create -f envs/m6anet.yaml
conda env create -f envs/nanocompore.yaml
# ... create other environments as needed

# Install pipeline dependencies
conda install -c conda-forge -c bioconda snakemake minimap2 samtools nanopolish
```

### Reference Data Setup

参考文件的**格式要求直接决定工具能否正常运行**，请严格按照下表准备：

| 文件 | 要求 | 用途 |
|---|---|---|
| `reference/genome.fa` | 标准 FASTA，染色体 ID 需与 GTF 一致；**必须 bgzip 解压为纯文本**（不可留 `.gz`）；需用 `samtools faidx` 建 `.fai` 索引 | 基因组比对（minimap2）、ELIGOS2、Epinano、NanoSPA |
| `reference/transcriptome.fa` | 标准 FASTA，序列 ID 格式：`ENST00000…` 或 `gene_id`；**由 pipeline 从 genome.fa + genes.gtf 自动生成**，也可手动准备 | 转录组比对（minimap2）、CHEUI、m6Anet、DENA、Nanocompore |
| `reference/genes.gtf` | Ensembl/GENCODE GTF 格式；**必须包含 `transcript_id`、`gene_id`、`gene_type` 属性**；染色体 ID 需与 genome.fa 一致 | 坐标转换、Guitar 图、yanocomp、xPore |
| `reference/genes.bed` | 标准 12 列 BED 或 6 列 BED；**由 pipeline 从 GTF 自动生成** | ELIGOS2 区域定义 |

> **染色体命名一致性**（常见坑）：如果 genome.fa 的 chr 列是 `1`（无前缀），genes.gtf 里也必须是 `1`，不能混合 `chr1`。可使用 `sed 's/^chr//'` 或 `awk '{print "chr"$0}'` 统一。

```bash
# 准备参考文件目录
mkdir -p reference
cd reference

# ====== 以人类 (GRCh38 / Ensembl 112) 为例 ======
# 1. 下载基因组
wget -O genome.fa.gz \
  ftp://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip genome.fa.gz                              # 必须解压为纯文本

# 2. 下载基因注释 (GTF)
wget -O genes.gtf.gz \
  ftp://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
gunzip genes.gtf.gz

# 3. 为基因组建立索引（后续 minimap2/samtools 需要）
samtools faidx genome.fa                        # 生成 genome.fa.fai

# 4. (可选) 手动生成转录组 — 也可让 pipeline 自动做
# gffread -w transcriptome.fa -g genome.fa genes.gtf
# awk '$3=="transcript" {print $1"\t"$4-1"\t"$5"\t"$12"\t0\t"$7}' genes.gtf \
#   | sed 's/"//g;s/;//g' > genes.bed
#
# cd ..
```

**以其他物种为例**：将上面 URL 中的 `homo_sapiens` 替换为对应物种目录名（如 `mus_musculus`、`arabidopsis_thaliana`）。对于非 Ensembl 注释源（如 Gencode、NCBI），请确保 GTF 中存在 `transcript_id` 和 `gene_id` 属性。

### Sample Naming Convention

```
data/
├── <sample_name>/         # 每个样本一个目录，目录名 = config.yaml 中的 samples 列表元素
│   └── fast5/             # 仅存放 .fast5 文件（可以多层子目录，Guppy 会递归）
│       ├── batch0_0.fast5
│       ├── batch0_1.fast5
│       └── …
├── <sample_name>/         # 例如：wt_rep1、ko_rep1、m6a_control …
│   └── fast5/
└── <control_name>/        # 对照样本同样放在 data/ 下，与处理样本目录结构一致
    └── fast5/
```

**重要规则**：
- 样本名**只能使用字母、数字、下划线**（`A–Z a–z 0–9 _`），避免空格、短横线、斜杠。
- **对照样本也必须写入 `config.yaml` 的 `samples:` 列表**，它们会被 Nanocompore / xPore / yanocomp / Epinano_DiffErr / DRUMMER 这些**对比型工具**自动配对。
- 对比型工具会按列表顺序**两两配对**：`sample1_vs_control1`、`sample2_vs_control2`。建议将「处理样本」和「对照样本」交错排列或通过配置文件显式指定。

## Usage

### Basic Usage

```bash
# Edit configuration file
nano config/config.yaml

# Run complete pipeline
snakemake --use-conda --cores 40

# Run specific tools only
snakemake --use-conda --cores 40 results/CHEUI/sample1/sample1_CHEUI_processed.txt

# Dry run to check workflow
snakemake --dry-run
```

### Tool Input Requirements Matrix

11 个工具对原始数据和前置步骤的依赖并不相同。如果你**只需要跑部分工具**，可以参考下面这张表来判断哪些前置步骤必须完成：

| Tool | 需要 genome BAM | 需要 transcriptome BAM | 需要 nanopolish eventalign | 需要对照样本 | 主要修饰类型 | 坐标系统 |
|---|:---:|:---:|:---:|:---:|---|---|
| CHEUI | ✗ | ✓ | ✓ | ✗ | m6A | transcriptome |
| m6Anet | ✗ | ✓ | ✓ | ✗ | m6A | transcriptome |
| Nanocompore | ✗ | ✓ | ✓ | ✓ | Any | transcriptome |
| ELIGOS2 | ✓ | ✗ | ✗ | ✗ | Any | genome |
| DENA | ✗ | ✓ | ✗ (uses fast5+Tombo) | ✗ | m6A (RRACH) | transcriptome |
| Epinano | ✓ | ✗ | ✗ | ✗ | m6A | genome |
| Epinano_DiffErr | ✓ | ✗ | ✗ | ✓ | Any (diff) | genome |
| MINES | ✗ | ✓ | ✗ (uses Tombo resquiggle) | ✗ | m6A | transcriptome |
| xPore | ✗ | ✓ | ✓ | ✓ | m6A | transcriptome |
| yanocomp | ✗ | ✓ | ✓ | ✓ | Any | transcriptome |
| NanoSPA | ✓ | ✗ | ✗ | ✗ | m6A, Ψ | genome |
| DRUMMER | ✓ | ✗ | ✗ | ✓ | Any | genome |

**简记规则**：以 "transcriptome 坐标" 工作的工具（CHEUI / m6Anet / Nanocompore / DENA / MINES / xPore / yanocomp）都需要 `transcriptome BAM`；以 "genome 坐标" 工作的工具（ELIGOS2 / Epinano / NanoSPA / DRUMMER）都需要 `genome BAM`。需要对比型工具（Nanocompore / xPore / yanocomp / Epinano_DiffErr / DRUMMER）必须同时提供至少两个样本（处理 + 对照）。

### Running Individual Tools

```bash
# 只跑 CHEUI（需要 transcriptome BAM + eventalign）
snakemake --use-conda --cores 40 \
  results/CHEUI/sample1/sample1_CHEUI_processed.txt

# 只跑 ELIGOS2（需要 genome BAM，不需要 eventalign）
snakemake --use-conda --cores 40 \
  results/ELIGOS2_solo/sample1/sample1_ELIGOS2_solo_processed.txt

# 只跑 Nanocompore 对比分析（需要两个样本 + eventalign）
snakemake --use-conda --cores 40 \
  results/Nanocompore/sample1_vs_control1/sample1_vs_control1_Nanocompore_processed.txt

# 只跑汇总和报告（假设各工具结果已经生成）
snakemake --use-conda --cores 4 \
  results/summary/modification_summary.tsv \
  results/report/RNAModBench_report.html
```

### Configuration

Edit `config/config.yaml` to specify:

- Sample names
- Tools to run
- Reference paths
- Filtering thresholds
- Computational resources

完整参数解释请参考 [docs/CONFIG_REFERENCE.md](docs/CONFIG_REFERENCE.md)。

### Sample Configuration

```yaml
# Sample configuration
samples:
  - sample1
  - sample2
  - control1

# Tools to run
tools:
  - CHEUI
  - m6Anet
  - Nanocompore
  - ELIGOS2

# Reference paths
reference_dir: "reference"
data_dir: "data"
results_dir: "results"

# Filtering thresholds
cheui:
  prob_threshold: 0.999
  ratio_threshold: 0.1

m6anet:
  prob_threshold: 0.5
  ratio_threshold: 0.1
```

## Input Data Structure

```
data/
├── sample1/
│   └── fast5/
│       ├── read1.fast5
│       ├── read2.fast5
│       └── ...
├── sample2/
│   └── fast5/
└── control1/
        └── fast5/
```

## Output Structure

```
results/
├── basecalling/
├── alignment/
├── nanopolish/
├── CHEUI/
│   └── sample1/
│       └── sample1_CHEUI_processed.txt
├── m6Anet/
├── Nanocompore/
├── summary/
│   ├── modification_summary.tsv
│   └── tool_comparison.tsv
└── report/
    └── RNAModBench_report.html
```

## Output Format

所有工具的 postprocess 脚本会把各自的原生输出统一成**相同的 7 列 TSV 格式**。

| 列 # | 字段名 | 类型 | 含义 |
|---|---|---|---|
| 1 | `Chr` | string | 染色体 ID 或转录本 ID（取决于工具使用的是 genome 还是 transcriptome 坐标，见上表） |
| 2 | `Start` | int | **BED 0-based 起始位置**（所有工具 postprocess 最后都做了 `Start -= 1` 对齐） |
| 3 | `End` | int | BED 0-based 结束位置（通常 = Start + 1，即修饰位点覆盖单个碱基） |
| 4 | `Status` | string | `Mod` / `Unmod`。只有被过滤阈值判为 "已修饰" 的位点会被保留（Status=Mod） |
| 5 | `Prob` | float | 概率或评分（0–1）。不同工具的含义不同：CHEUI/m6Anet/xPore/NanoSPA 是概率；ELIGOS2/Nanocompore 是 p-value；Epinano 是 delta error；DENA/MINES 则填入 mod_ratio |
| 6 | `Strand` | char | `+` / `-` / `*`（部分工具无法确定 strand，用 `*`） |
| 7 | `mod_ratio` | float | 修饰化学计量比（0–1）。表示被检测为修饰的 reads 占总 reads 的比例 |

> **坐标系统重要提示**：CHEUI、m6Anet、Nanocompore、DENA、MINES、xPore、yanocomp 的 `Chr` 列是**转录本 ID**（如 `ENST00000367770`），位置是**相对于该转录本 5' 端的 0-based 偏移**。如果你需要基因组坐标，请使用 `r2d_liftover.py`（R2Dtool）做坐标转换，结果会写到 `results/summary/liftover_summary.tsv`。

输出文件命名规则：
- 单样本工具：`results/<Tool>/<sample>/<sample>_<Tool>_processed.txt`
- 对比型工具：`results/<Tool>/<sample>_vs_<control>/<sample>_vs_<control>_<Tool>_processed.txt`

完整的输出格式说明见 [docs/OUTPUT_FORMAT.md](docs/OUTPUT_FORMAT.md)。

## Interpreting Results

拿到结果后建议按以下顺序检查：

1. **每个工具报告的修饰位点数量**是否在合理范围（m6A 通常占 mRNA 的 0.1%–1% 位点）。如果某个工具报告的位点显著偏离 10x 以上，检查该工具的过滤阈值是否合理。
2. **工具一致性**：查看 `results/summary/modification_summary.tsv` 中的 `Tool_Count` 列。若某个位点被多个工具同时报告，可信度显著更高。可按 `Tool_Count >= 3` 过滤得到高可信度共识位点。
3. **染色体/转录本分布**：Guitar 图（`results/summary/guitar_plots_mrna.png`）可展示修饰位点在 mRNA 结构上的偏好（m6A 倾向富集于 3'UTR 和终止密码子附近）。
4. **Motif 验证**：对输出文件运行 `scripts/extract_5mer.py` 添加 k-mer 上下文，检查是否富集于 `RRACH`（m6A 的经典 motif）。

## Hardware Requirements Estimation

以 1 个样本、~1M reads / ~10G fast5 为例的估算：

| 步骤 | 内存 | 核心数 | 预估时间 |
|---|---|---|---|
| Guppy basecalling (CPU) | 16G | 40 | 30–60 分钟 |
| Guppy basecalling (GPU) | 8G GPU + 16G CPU | 20 | 5–10 分钟 |
| minimap2 alignment | 16G | 40 | 5–10 分钟 |
| nanopolish eventalign | 64G | 40 | 2–4 小时（**最慢步骤**） |
| CHEUI | 32G | 40 | 30–60 分钟 |
| m6Anet | 16G | 40 | 15–30 分钟 |
| ELIGOS2 | 8G | 12 | 30–60 分钟 |
| Epinano | 8G | 12 | 10–20 分钟 |
| Nanocompore | 16G | 20 | 1–2 小时 |
| 汇总/报告/可视化 | 4G | 1 | < 5 分钟 |

> **nanopolish eventalign 是整体瓶颈**。如果你的数据集很大（>5M reads），考虑分批处理或使用 GPU。CHEUI 的深度学习模型在 GPU 上可获得约 5–10 倍加速。

## Quality Control

The pipeline includes several QC steps:

- Basecalling quality filtering
- Read alignment quality control
- Coverage-based filtering
- Statistical significance testing
- Tool-specific quality metrics

## Performance

- **Parallel Processing**: Utilizes multiple cores for speed
- **Memory Efficient**: Optimized for large datasets
- **GPU Support**: Accelerated computation for deep learning tools
- **Modular Execution**: Run only required tools



## Troubleshooting

### Common Issues with Concrete Error Messages

| 常见报错文本 | 原因 | 解决方案 |
|---|---|---|
| `MissingInputException` / `Missing input files for rule` | Snakefile 找不到输入文件（fast5 路径错、样本名拼写不一致、参考文件缺失） | 检查 `config.yaml` 的 `samples:` 列表与 `data/` 目录名是否严格一致；确认 `reference/genome.fa`、`genes.gtf` 存在 |
| `MemoryError` / `Segmentation fault` / `Killed` | 内存不足，通常发生在 nanopolish eventalign 或 CHEUI 大数据集上 | 1) 减少 `config.yaml` 中各工具的 threads；2) 减小样本规模；3) 用 `--resources mem_mb=32000` 限制 Snakemake 总内存 |
| `[E::fai_build_core]` / `could not parse` | genome.fa 格式异常（未解压、空行、非标准 FASTA 头） | 重新 `gunzip`、`samtools faidx` 重新建索引 |
| `OSError: File not found: resources/R2Dtool/.../r2d` | `utilities.r2d_tool` 路径不正确或未编译 R2Dtool | 克隆 [chrisam/r2d-tool](https://github.com/chrisam/r2d-tool) 并 `cargo build --release`，更新 `config.yaml` |
| `KeyError: 'transcript_id'` (R/Guitar 相关) | genes.gtf 缺少 `transcript_id` 属性 | 使用 Ensembl/Gencode 标准 GTF；不要用 NCBI RefSeq 的 GFF 直接改名充当 GTF |
| `[C++ Error] std::bad_alloc` （nanopolish） | eventalign 同时处理过多 reads | 调低 `nanopolish.threads` 到 10 以下 |
| `CUDA error: out of memory` （CHEUI / 深度学习工具） | GPU 显存不足 | 减小 batch size（在对应工具脚本中）或改为 CPU 模式 |
| 某工具输出 `_processed.txt` 为空 | 过滤阈值太严格，没有位点通过 | 适当放松 `prob_threshold` / `pvalue_threshold` / `ratio_threshold` |
| `FileNotFoundError: data/sample1/fast5/` | fast5 目录不存在或层级不对 | 确保结构是 `data/<sample_name>/fast5/*.fast5` |

### Debug Mode

```bash
# 先 dry-run，看规则链路是否正确
snakemake --dry-run --use-conda -n -p
# 查看完整的 DAG（有向无环图）
snakemake --dag | dot -Tpng > workflow_dag.png
#  verbose 输出
snakemake --use-conda --cores 40 --verbose
# 单步调试某个 rule
snakemake --use-conda --cores 1 results/CHEUI/sample1/sample1_CHEUI_processed.txt --debug
# 失败后保留临时文件以便排查
snakemake --use-conda --cores 40 --keep-going --rerun-incomplete
```

更多常见问题请见 [docs/FAQ.md](docs/FAQ.md)。

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## License

This pipeline is released under the MIT License. See LICENSE file for details.

## Support

For issues and questions:

1. Check existing GitHub issues
2. Create a new issue with detailed description
3. Include error logs and configuration files

## Acknowledgments & Citations

RNAModBench 集成了以下研究团队开发的工具。如果你在工作中使用了对应的工具，请引用原始文献：

| Tool | Reference |
|---|---|
| **CHEUI** | Liu et al., "Quantifying RNA modifications at single-molecule resolution", *Nat. Biotechnol.*, 2023 |
| **ELIGOS2** | Jenjaroenpun et al., "Detection of internal RNA modifications using nanopore sequencing", *Genome Biol.*, 2021 |
| **m6Anet** | Hendra et al., "m6Anet detects m6A modifications from direct RNA-Seq data", *Nat. Methods*, 2022 |
| **Nanocompore** | Leger et al., "Nanocompore identifies context-dependent RNA modifications", *Nat. Commun.*, 2021 |
| **DENA** | Li et al., "Detection of m6A modifications in direct RNA-Seq", *Genome Biol.*, 2022 |
| **Epinano** | Liu et al., "Accurate detection of m6A RNA modifications in native RNA sequences", *Nat. Commun.*, 2019 |
| **MINES** | Begik et al., "Machine learning for detection of m6A modifications", *NAR*, 2019 |
| **xPore** | Pratanwanich et al., "Identification of differential RNA modifications from nanopore sequencing", *Nat. Biotechnol.*, 2021 |
| **yanocomp** | Ma et al., "Detection of differential RNA modifications using GMMs", *Genome Biol.*, 2022 |
| **NanoSPA** | Dong et al., "Simultaneous profiling of multiple RNA modifications", *Cell Genomics*, 2023 |
| **DRUMMER** | Ortiz et al., "DRUMMER detects m6A modifications from nanopore data", *NAR*, 2021 |
| **R2Dtool** | Chrisam et al., "Transcript-to-genome coordinate liftover", GitHub, 2024 |
| **Guitar** (R package) | Li et al., "Guitar: an R/Bioconductor package for gene annotation-guided visualization", *Bioinformatics*, 2021 |

> 注意：以上引用信息为通用参考，实际使用时请查阅各工具的原始发表论文以获取最新 DOI 和完整作者信息。

---

## Documentation Navigation

本项目附带以下文档，按需求查阅：

| 文件 | 用途 |
|---|---|
| [docs/INSTALL.md](docs/INSTALL.md) | 完整的安装指南（Conda 环境、R 包、外部二进制） |
| [docs/TUTORIAL.md](docs/TUTORIAL.md) | 最小工作示例：用公开数据集跑通完整流程 |
| [docs/TOOLS_OVERVIEW.md](docs/TOOLS_OVERVIEW.md) | 11 个工具的输入需求、坐标系统、运行时间对比表 |
| [docs/OUTPUT_FORMAT.md](docs/OUTPUT_FORMAT.md) | 标准 7 列 TSV 输出的完整解释 |
| [docs/CONFIG_REFERENCE.md](docs/CONFIG_REFERENCE.md) | `config/config.yaml` 每个参数的详细解释 |
| [docs/FAQ.md](docs/FAQ.md) | 常见问题解答 |
| [CONTRIBUTING.md](CONTRIBUTING.md) | 贡献指南 |
| [CHANGELOG.md](CHANGELOG.md) | 版本变更记录 |

For more information, please visit the [GitHub repository](https://github.com/autosomal/RNAModBench).