# 工具输入需求矩阵

下表总结了 RNAModBench 中集成的 12 个 RNA 修饰检测工具对原始数据和前置步骤的依赖。

## 目录

- [工具输入需求矩阵](#工具输入需求矩阵)
- [单样本工具 vs 对比型工具](#单样本工具-vs-对比型工具)
- [坐标系统一览](#坐标系统一览)
- [信号工具 vs 碱基比对工具](#信号工具-vs-碱基比对工具)
- [典型运行时间估算](#典型运行时间估算)
- [工具特定的参考资源](#工具特定的参考资源)

---

## 工具输入需求矩阵

| 工具 | 样本类型 | 需要 fast5 | 需要 minimap2 alignment | 需要 nanopolish eventalign | 需要 Tombo resquiggle | 染色体命名 |
|---|---|---|---|---|---|---|
| **CHEUI** | single | Yes | transcriptome only | Yes | No | transcript id |
| **m6Anet** | single | Yes | transcriptome only | Yes | No | transcript id |
| **DENA** | single | Yes | transcriptome only | No | Yes | transcript id |
| **MINES** | single | Yes | transcriptome only | No | Yes | transcript id |
| **ELIGOS2** | single | No | genome only | No | No | chromosome |
| **Epinano** | single | No | genome only | No | No | chromosome |
| **NanoSPA** | single | No | genome only | No | No | chromosome |
| **Nanocompore** | case vs control | Yes | transcriptome only | Yes | No | transcript id |
| **xPore** | case vs control | Yes | transcriptome only | Yes | No | transcript id |
| **yanocomp** | case vs control | Yes | transcriptome only | Yes | No | transcript id |
| **Epinano_DiffErr** | case vs control | No | genome only | No | No | chromosome |
| **DRUMMER** | case vs control | No | genome only | No | No | chromosome |

---

## 单样本工具 vs 对比型工具

**单样本工具**（仅处理一个样本，无需对照）：

- CHEUI, m6Anet, DENA, MINES, ELIGOS2, Epinano, NanoSPA
- 每个样本独立产生结果
- `samples` 列表中的每个样本都会执行

**对比型工具**（需要同时提供处理样本和对照样本）：

- Nanocompore, xPore, yanocomp, Epinano_DiffErr, DRUMMER
- 工具会按 `samples` 列表顺序两两配对，自动生成 `<处理>_vs_<对照>` 输出目录
- 在 `config.yaml` 中必须至少包含 2 个样本

> 提示：若只关注「是否存在修饰」而非「处理 vs 对照的差异」，可以跳过对比型工具以节省计算时间。

---

## 坐标系统一览

RNAModBench 中使用 **两种坐标系统并存**：

### transcriptome 坐标
- 使用 **转录本 ID**（如 `ENST00000367770`、`XM_00000000.1`）作为 `Chr` 列
- position 为 **1-based** 位置（从转录本 5' 端开始计数）
- postprocess 后统一为 **BED 0-based** [start, end)
- 工具：CHEUI、m6Anet、Nanocompore、DENA、MINES、xPore、yanocomp

### genome 坐标
- 使用 **染色体 ID**（如 `chr1` / `1` / `MT`）作为 `Chr` 列
- position 为 **1-based** 位置（从染色体 p 端开始计数）
- postprocess 后统一为 **BED 0-based** [start, end)
- 工具：ELIGOS2、Epinano、NanoSPA、DRUMMER

> **坐标转换**：若需要将 transcriptome 坐标结果投影回 genome，使用
> `scripts/r2d_liftover.py`（调用 [R2Dtool](https://github.com/chrisam/r2d-tool)）。

---

## 信号工具 vs 碱基比对工具

### 信号工具（signal-level tools）
- **处理 raw nanopore signal**：直接利用 FAST5 中的 pico-ampere 时间序列
- 通常能提供更高的敏感度，尤其对低化学计量的修饰
- **CHEUI、m6Anet、Nanocompore、DENA、xPore、yanocomp**

### 碱基比对工具（alignment-level tools）
- **处理 base-called FASTQ**：基于与参考序列的比对结果寻找错配 / Q-score 偏差
- 运行速度更快，不依赖 FAST5 文件
- **ELIGOS2、Epinano、NanoSPA、DRUMMER**

> **建议**：数据量大且有 FAST5 → 优先跑信号工具；只有 FASTQ → 用比对工具。

---

## 典型运行时间估算

以 1 个样本，~1M reads / ~10GB fast5 / 32 CPU 为例：

| 步骤 | 时间 | 瓶颈 |
|---|---|---|
| Guppy basecalling (CPU) | 30-60 min | I/O + CPU |
| Guppy basecalling (GPU) | 5-10 min | GPU |
| minimap2 + samtools | < 5 min | I/O |
| nanopolish eventalign | 2-4 hr | 全局瓶颈 |
| CHEUI model inference | 30-60 min | I/O + CPU/GPU |
| m6Anet inference | 15-30 min | CPU |
| ELIGOS2 | 10-20 min | CPU |
| Epinano | 5-15 min | CPU |
| DENA | 20-40 min | CPU |
| Nanocompore sampcomp | 20-40 min | CPU + memory |
| MINES | 15-30 min | CPU |
| 汇总与画图 (summary) | < 5 min | CPU |

---

## 工具特定的参考资源

| 工具 | 需要的附加文件 | 用途 |
|---|---|---|
| CHEUI | `model_kmer.csv` + `.h5` 模型文件 | 5-mer lookup + deep learning inference |
| Epinano | `rrach.q3.mis3.del3.linear.dump` | SVM 模型（线性核） |
| MINES | `.names` + 模型文件 | 5-mer 指定分类器 |
| DENA | 目录内的 LSTM 模型文件 | LSTM 推断 |
| yanocomp | `—`（无） | 仅使用 eventalign 的统计 |
| R2Dtool | 编译好的 `r2d` 二进制 | transcript→genome 坐标转换 |

参见 `config.yaml` 中的 `cheui.*`、`epinano.model`、`mines.kmer_models`、
`dena.model`、`utilities.r2d_tool` 字段配置具体路径。
