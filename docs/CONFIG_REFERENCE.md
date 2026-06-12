# config.yaml 参数详解

本文档按功能分组说明 `config/config.yaml` 中每个字段的用途、类型和默认值。
对于标注 **[NOT YET USED]** 的字段，Snakefile 尚未接入，修改它们不会影响运行。

## 目录

- [基本配置](#基本配置)
- [Guppy basecalling](#guppy-basecalling)
- [比对与信号处理](#比对与信号处理)
- [CHEUI](#cheui)
- [ELIGOS2](#eligos2)
- [m6Anet](#m6anet)
- [Nanocompore](#nanocompore)
- [DENA](#dena)
- [Epinano / Epinano_DiffErr](#epinano--epinano_differr)
- [DRUMMER](#drummer)
- [MINES](#mines)
- [xPore](#xpore)
- [yanocomp](#yanocomp)
- [NanoSPA](#nanospa)
- [辅助工具参数](#辅助工具参数)
- [参考文件路径](#参考文件路径)
- [预留配置项（未接入）](#预留配置项未接入)

---

## 基本配置

| 字段 | 类型 | 说明 |
|---|---|---|
| `samples` | list[str] | 样本目录名。长度 = 每个单样本工具执行次数；对比型工具会按列表两两配对 |
| `tools` | list[str] | 要运行的工具列表。有效值：CHEUI, ELIGOS2, m6Anet, Nanocompore, DENA, Epinano, Epinano_DiffErr, DRUMMER, MINES, xPore, yanocomp, NanoSPA |
| `data_dir` | string | 存放 `<sample>/fast5/` 的父目录，通常为 `"data"` |
| `reference_dir` | string | 参考文件目录，通常为 `"reference"` |
| `results_dir` | string | 所有输出文件的父目录，通常为 `"results"` |

---

## Guppy basecalling

| 字段 | 类型 | 说明 | 默认 |
|---|---|---|---|
| `guppy.config` | string | 配置文件名。RNA 常用 `rna_r9.4.1_70bps_hac.cfg`（高准确度） | `rna_r9.4.1_70bps_hac.cfg` |
| `guppy.num_callers` | int | 并行 basecalling worker 数 | `4` |
| `guppy.threads_per_caller` | int | 每个 worker 的线程数 | `20` |
| `guppy.total_threads` | int | 总线程上限（给 Snakemake `--threads` 用） | `80` |

---

## 比对与信号处理

| 字段 | 类型 | 说明 | 默认 |
|---|---|---|---|
| `alignment.threads` | int | minimap2 同时用于 transcriptome 和 genome 比对的线程 | `40` |
| `nanopolish.threads` | int | `nanopolish eventalign` 线程数。大数据集（>10G fast5）建议 ≤ 20 防止内存崩溃 | `40` |

---

## CHEUI

| 字段 | 类型 | 说明 | 默认 |
|---|---|---|---|
| `cheui.kmer_model` | string | 5-mer 模型 CSV 路径 | `resources/CHEUI/kmer_models/model_kmer.csv` |
| `cheui.model1` | string | 阶段 1（初筛）的 .h5 模型 | `resources/CHEUI/.../model1.h5` |
| `cheui.model2` | string | 阶段 2（精筛）的 .h5 模型 | `resources/CHEUI/.../model2.h5` |
| `cheui.threads` | int | infer 阶段的 CPU 线程 | `40` |
| `cheui.prob_threshold` | float | 阶段 2 概率下限（高于此值的位点输出到 processed.txt） | `0.999` |
| `cheui.ratio_threshold` | float | 化学计量比下限 | `0.1` |

> CHEUI 的过滤很激进（`0.999`）。如果你希望更敏感，可以放宽到 `0.9`。

---

## ELIGOS2

| 字段 | 类型 | 说明 | 默认 |
|---|---|---|---|
| `eligos2.threads` | int | 并行线程数（每个 region 一个 worker） | `12` |
| `eligos2.max_depth` | int | 每个位点最大读深度（过大则 downsample） | `2000000` |
| `eligos2.min_depth` | int | 最小覆盖度，低于此值跳过 | `5` |
| `eligos2.padj_threshold` | float | BH 校正后 p-value 上限 | `0.0001` |
| `eligos2.oddr_threshold` | float | odds ratio 下限（修饰效应量） | `1.2` |

---

## m6Anet

| 字段 | 类型 | 说明 | 默认 |
|---|---|---|---|
| `m6anet.dataprep_threads` | int | `m6anet dataprep` 阶段线程数 | `12` |
| `m6anet.inference_threads` | int | `m6anet infer` 阶段线程数 | `40` |
| `m6anet.readcount_max` | int | 单 transcript 的 read 数上限（内存保护） | `2000000` |
| `m6anet.prob_threshold` | float | 位点 m6A 概率下限 | `0.5` |
| `m6anet.ratio_threshold` | float | 化学计量比下限 | `0.1` |

---

## Nanocompore

| 字段 | 类型 | 说明 | 默认 |
|---|---|---|---|
| `nanocompore.threads` | int | sampcomp 并行线程 | `20` |
| `nanocompore.min_coverage` | int | 每个位点最小 reads 覆盖 | `5` |
| `nanocompore.min_ref_length` | int | 最短 transcript 长度 | `10` |
| `nanocompore.pvalue_threshold` | float | GMM p-value 上限 | `0.05` |
| `nanocompore.lor_threshold` | float | log odds ratio 绝对值下限 | `0.5` |

---

## DENA

| 字段 | 类型 | 说明 | 默认 |
|---|---|---|---|
| `dena.motif` | string | 目标 motif（IUPAC 简并碱基，如 `RRACH`） | `RRACH` |
| `dena.corr_grp` | string | Tombo 校正后信号组名 | `RawGenomeCorrected_000` |
| `dena.windows` | string | 目标位点上下游窗口（空格分隔，如 `"2 2"`） | `"2 2"` |
| `dena.processes` | int | LSTM 并行进程数 | `40` |
| `dena.model` | string | LSTM 模型目录 | `resources/DENA/DENA_LSTM_Model` |
| `dena.ratio_threshold` | float | 化学计量比下限 | `0.1` |
| `dena.coverage_threshold` | int | 最小覆盖度 | `20` |

---

## Epinano / Epinano_DiffErr

| 字段 | 类型 | 说明 | 默认 |
|---|---|---|---|
| `epinano.threads` | int | 并行线程 | `12` |
| `epinano.model` | string | SVM 模型文件（linear.dump） | `resources/Epinano/models/rrach.q3.mis3.del3.linear.dump` |
| `epinano.columns` | string | 作为特征的列（1-based，逗号分隔） | `"8,13,23"` |
| `epinano.delta_threshold` | float | delta error 差异下限（仅 DiffErr） | `0.1` |

---

## DRUMMER

| 字段 | 类型 | 说明 | 默认 |
|---|---|---|---|
| `drummer.threads` | int | 并行线程 | `12` |
| `drummer.max_depth` | int | 最大读深度 | `2000000` |
| `drummer.min_depth` | int | 最小读深度 | `5` |
| `drummer.pvalue_threshold` | float | OR 校正 p-value 上限 | `0.05` |
| `drummer.frac_diff_threshold` | float | 处理/对照比例差下限 | `0.1` |

---

## MINES

| 字段 | 类型 | 说明 | 默认 |
|---|---|---|---|
| `mines.kmer_models` | string | k-mer 模型 `.names` 文件路径 | `resources/MINES/Final_Models/names.txt` |
| `mines.coverage_threshold` | int | 最小覆盖度 | `20` |
| `mines.ratio_threshold` | float | 修饰比例下限 | `0.1` |

---

## xPore

| 字段 | 类型 | 说明 | 默认 |
|---|---|---|---|
| `xpore.dataprep_threads` | int | 预处理阶段线程 | `40` |
| `xpore.readcount_max` | int | 单 transcript read 数上限（内存保护） | `2000000` |
| `xpore.prob_threshold` | float | 差异修饰概率下限 | `0.5` |

---

## yanocomp

| 字段 | 类型 | 说明 | 默认 |
|---|---|---|---|
| `yanocomp.prep_threads` | int | eventalign → HDF5 预处理阶段线程 | `8` |
| `yanocomp.test_threads` | int | GMM 检验阶段线程 | `8` |
| `yanocomp.n_components` | int | 高斯混合成分数（越大越精细但更慢） | `50` |
| `yanocomp.fdr` | float | Benjamini-Hochberg FDR 控制 | `0.05` |
| `yanocomp.pvalue_threshold` | float | p-value 硬过滤 | `0.05` |

---

## NanoSPA

| 字段 | 类型 | 说明 | 默认 |
|---|---|---|---|
| `nanospa.prob_threshold` | float | 修饰概率下限 | `0.5` |

---

## 辅助工具参数

| 字段 | 类型 | 说明 | 默认 |
|---|---|---|---|
| `utilities.r2d_tool` | string | R2Dtool 可执行文件的绝对路径 | `resources/R2Dtool/target/release/r2d` |

---

## 参考文件路径

| 字段 | 类型 | 说明 |
|---|---|---|
| `reference_files.genome` | string | 基因组 FASTA（已解压为 `.fa`） |
| `reference_files.transcriptome` | string | 转录组 FASTA |
| `reference_files.genes_gtf` | string | Ensembl 风格 GTF（需含 transcript_id / gene_id） |
| `reference_files.genes_bed` | string | BED 格式基因注释（ELIGOS2 使用） |

---

## 预留配置项（未接入）

以下字段已写入 `config.yaml`，但 **Snakefile 尚未实现**。修改它们不会影响当前流水线。

### 自动生成参考文件
- `create_transcriptome`：是否自动由 genome + GTF 生成 transcriptome
- `convert_gtf_to_bed`：是否自动将 GTF 转换为 BED
- `create_gene_annotations`：是否自动构建基因注释

### QC
- `qc.min_read_length` / `min_read_quality` / `min_coverage` / `max_coverage`

### 输出控制
- `output.bed_format` / `liftover_to_genome` / `include_5mer_context`

### 计算资源
- `resources.max_memory` / `max_cpus` / `gpu_available` / `gpu_devices`

### 日志
- `logging.level` / `save_intermediates` / `cleanup_temp_files`

### 物种
- `species.name` / `codon_table` / `mitochondrial_genome` / `annotation_source`

> 如果你希望接入这些字段，可 fork 本仓库并修改 `Snakefile` 中对应的 `config["xxx"]`
> 读取语句。欢迎提交 PR。
