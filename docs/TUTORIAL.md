# RNAModBench 最小工作示例

本文档演示如何用一组**公开数据**跑通整个流水线。所有示例使用的是可自由下载的
nanopore 数据，假设你已经完成了 [INSTALL.md](INSTALL.md) 中的依赖安装。

## 目录

- [准备参考文件](#准备参考文件)
- [准备原始数据（fast5）](#准备原始数据fast5)
- [编辑 config.yaml](#编辑-configyaml)
- [Dry-run 检查](#dry-run-检查)
- [执行完整流水线](#执行完整流水线)
- [只跑部分工具（例如只跑 CHEUI + ELIGOS2）](#只跑部分工具例如只跑-cheui--eligos2)
- [检查输出](#检查输出)

---

## 准备参考文件

以人类 GRCh38 / Ensembl 112 为例：

```bash
mkdir -p reference
cd reference

# 1) 基因组 FASTA（~3 GB）
wget -O genome.fa.gz \
  ftp://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip genome.fa.gz
samtools faidx genome.fa        # 生成 .fai 索引

# 2) 转录组 FASTA
gffread -w transcriptome.fa -g genome.fa Homo_sapiens.GRCh38.112.gtf

# 3) GTF 注释（直接下载）
wget -O genes.gtf.gz \
  ftp://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
gunzip genes.gtf.gz

cd ..
```

> **染色体命名一致性**：如果后续发现 `genome.fa` 的 `>` 行写的是 `1  dna:...`
> 而 `genes.gtf` 第一列写的也是 `1`，则没问题；不要混合使用 `chr1` 和 `1`。

---

## 准备原始数据（fast5）

### 方式 A：公开示例数据（小型、快速）

建议先从小型数据集入手，验证流水线能跑通：

```bash
mkdir -p data/treatment_rep1
cd data/treatment_rep1

# 从 ENA / SRA 下载几个 fast5 作为玩具数据
# 例如：SRA accession SRXxxxxxxxx 的 fast5 子集
```

> 对于真正的测试，建议下载 1-2 个 SRA runs（~10G 量级）。

### 方式 B：使用你自己的数据

将你的 fast5 文件放在对应样本子目录下即可：

```
data/
├── treatment_rep1/fast5/*.fast5
├── treatment_rep2/fast5/*.fast5
├── control_rep1/fast5/*.fast5
└── control_rep2/fast5/*.fast5
```

> 目录名必须与 `config.yaml` 中 `samples` 列表的元素**完全一致**（区分大小写）。

---

## 编辑 config.yaml

```yaml
samples:
  - treatment_rep1
  - control_rep1

tools:
  - CHEUI
  - ELIGOS2
  - Epinano

data_dir:      "data"
reference_dir: "reference"
results_dir:   "results"

# 其余参数保持默认即可；也可按 docs/CONFIG_REFERENCE.md 微调过滤阈值
```

---

## Dry-run 检查

**第一步务必先 dry-run**，确认规则链正确、所有输入文件存在：

```bash
snakemake --dry-run --printshellcmds -j 4
```

若你只想看哪些 rule 会被执行、总共需要多少 jobs，加 `-n` 即可：

```bash
snakemake -n -j 4
```

预期输出末端会显示类似 `Job counts: count jobs 1 ...` 的摘要。

---

## 执行完整流水线

```bash
# 推荐：限制并发，避免 nanopolish eventalign 占用所有内存
snakemake --use-conda -j 24 --resources mem_mb=64000

# 如果没有 GPU，则需在 cheui.threads 调低
```

---

## 只跑部分工具（例如只跑 CHEUI + ELIGOS2）

```bash
# 方式 1：在 config.yaml 的 tools 列表里只保留要跑的工具
# 方式 2：直接指定目标文件（Snakemake 会自动推断依赖）
snakemake --use-conda -j 24 \
  results/CHEUI/treatment_rep1/treatment_rep1_CHEUI_processed.txt \
  results/ELIGOS2_solo/treatment_rep1/treatment_rep1_ELIGOS2_solo_processed.txt
```

---

## 检查输出

完整流水线成功运行后，你会看到：

```
results/
├── CHEUI/treatment_rep1/treatment_rep1_CHEUI_processed.txt
├── ELIGOS2_solo/treatment_rep1/treatment_rep1_ELIGOS2_solo_processed.txt
├── Epinano/treatment_rep1/treatment_rep1_Epinano_processed.txt
│
├── summary/
│   ├── modification_summary.tsv     # 每个工具的修饰位点数汇总
│   ├── tool_comparison.tsv          # 工具间一致性
│   ├── tool_overlap.png             # Venn 图（若 ≥2 工具跑通）
│   ├── modification_distribution.png# 染色体分布柱状图
│   ├── guitar_plots_mrna.png        # mRNA 上的分布
│   └── liftover_summary.tsv         # (若启用)
│
├── report/
│   └── RNAModBench_report.html      # Bootstrap 风格的 HTML 报告
│
└── qc/
    └── depth_plots.png              # 每个样本的覆盖度
```

## 下一步

- 打开 `results/report/RNAModBench_report.html` 查看可视化汇总
- 根据 `modification_summary.tsv` 检查每个工具的位点数量是否合理
- 如需做基因组浏览器可视化，可将 `*_processed.txt` 转换为 BED12
  （列 1-3 对应 Chr/Start/End）
- 若需要针对特定基因/转录本深挖：在 `script/r2d_liftover.py`
  上单独调用一次即可

> **排错提示**：若某一步失败，查看 README 的
> [常见报错文本 → 解决方案](https://github.com/autosomal/RNAModBench#troubleshooting)
> 表，或直接看 Snakemake 打印的 `.snakemake/log/` 目录下的日志文件。
