# 常见问题解答 (FAQ)

## 运行相关

**Q1. 我只想跑 3 个工具，为什么还要把所有工具的 conda 环境都装好？**

不需要。Snakemake 只会创建 `config.yaml` 中 `tools:` 列表列出的工具所需的
conda 环境。如果你只写 `tools: [CHEUI, ELIGOS2]`，那么其他工具的环境完全不需要装。

**Q2. `snakemake --dry-run` 报 `MissingInputException`，怎么办？**

检查：
1. `data/<sample_name>/fast5/` 目录下是否存在 `.fast5` 文件
2. 目录名是否与 `samples:` 列表中的元素完全一致（区分大小写）
3. `reference/genome.fa`、`genes.gtf` 是否存在（注意 `.fa` 和 `.fasta` 的区别）
4. 是否忘记 `gunzip`（genome.fa.gz 是不能直接用的）

**Q3. nanopolish eventalign 运行非常慢 / 内存炸了，怎么办？**

nanopolish 是整个流水线最慢的步骤。优化：
1. 把 `config.yaml` 的 `nanopolish.threads` 调低到 10–16
2. 限制 Snakemake 总资源：`snakemake --resources mem_mb=32000`
3. 若只是 m6A 分析，其实可以跑 **ELIGOS2 / Epinano**（不需要 eventalign）

## 结果相关

**Q4. 为什么某个工具报告 0 个修饰位点？**

可能是过滤阈值过于严格。例如 CHEUI 默认 `prob_threshold=0.999`。尝试：
1. 调低对应工具的 `*_threshold`
2. 检查 basecalling 的质量 —— Q7 以下数据可能工具无法识别修饰
3. 检查是否有足够多的 reads 覆盖目标基因/转录本

**Q5. 不同工具结果完全不一致，怎么办？**

这是正常现象 —— 每种工具对 "修饰" 的定义和信号模型不同。建议：
1. 使用 `results/summary/modification_summary.tsv` 查看被 **3+ 工具**
   共同报告的位点，这些是最可信的
2. 关注 motif 是否合理（m6A 的话应富集于 `RRACH`）
3. 检查数据是否真的包含修饰（对照样本是否为 IVT / 去修饰的）

**Q6. 为什么 transcriptome 坐标的工具结果和 genome 坐标的工具结果没法直接对齐？**

这是设计选择 —— CHEUI 等工具输出的是转录本坐标，需要通过 R2Dtool 投影回
genome。执行：

```bash
python scripts/r2d_liftover.py \
    --input results/CHEUI/sample1/sample1_CHEUI_processed.txt \
    --gtf reference/genes.gtf \
    --output results/CHEUI/sample1/sample1_CHEUI_genomic.tsv
```

## 安装相关

**Q7. conda 安装失败 / 下载太慢**

- 切换镜像站
- 使用 `mamba` 代替 conda 作为 solver（`conda install mamba -n base -c conda-forge`）
- 或者按 `docs/INSTALL.md` 中的 `pip install` 方式只装你实际要用的工具

**Q8. R 包 `Guitar` 安装失败**

```r
BiocManager::install("Guitar")  # 需要 Bioconductor >= 3.15
packageVersion("Guitar")         # 输出 >= 2.12.0 即正常
```

如果还报缺包，再装：`BiocManager::install(c("GenomicFeatures", "rtracklayer"))`。

## 数据相关

**Q9. 我的样本没有 fast5（只有 fastq），还能用 RNAModBench 吗？**

可以跑以下工具：**ELIGOS2、Epinano、NanoSPA、DRUMMER**。
这些工具只需要 basecalled reads。在 `tools:` 里把 signal-level tools 去掉即可。

**Q10. 我可以跑人类以外的物种吗？**

可以。只需替换 `reference/genome.fa` 和 `genes.gtf` 为对应物种的文件，
参考文件的 chromosome 命名必须前后一致。

## 其它

**Q11. 如何把 RNAModBench 接入集群调度系统（SLURM / SGE）？**

Snakemake 原生支持 `--profile`。示例 SLURM profile：

```
# ~/.config/snakemake/slurm/config.yaml
jobs: 50
cluster: "sbatch --nodes=1 --ntasks={threads} --mem={resources.mem_mb}M -t {resources.runtime} -J {rule}"
default-resources: [mem_mb=4000, runtime=120]
```

然后：`snakemake --profile slurm`

**Q12. 如何重新跑某一步，而不是清掉所有结果？**

```bash
snakemake --forceall --rerun-incomplete results/CHEUI/sample1/sample1_CHEUI_processed.txt
```

或者用 `--touch` 把已存在文件标记为最新：

```bash
snakemake --touch --cores 1
```
