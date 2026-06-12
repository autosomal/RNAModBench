# 标准输出格式说明

RNAModBench 的所有工具经过 `scripts/postprocess_*.py` 之后，都会被转换为相同的
**7 列 TSV**（制表符分隔文本文件）。这种一致性让汇总、工具对比和下游分析变得
简单可靠。

## 目录

- [标准 7 列定义](#标准-7-列定义)
- [字段的工具差异](#字段的工具差异)
- [输出文件命名约定](#输出文件命名约定)
- [工具对比文件格式](#工具对比文件格式)
- [坐标系统说明](#坐标系统说明)
- [下游分析示例](#下游分析示例)

---

## 标准 7 列定义

每个 `*_processed.txt` 文件：

| 列 # | 字段 | 类型 | 含义 |
|---|---|---|---|
| 1 | `Chr` | string | 染色体 ID 或转录本 ID（具体见"坐标系统"） |
| 2 | `Start` | int | **BED 0-based** 起始位置 |
| 3 | `End` | int | BED 0-based 结束位置（通常 = Start + 1） |
| 4 | `Status` | string | 固定为 `Mod`（只有被判定为修饰的位点保留） |
| 5 | `Prob` | float | 概率 / p-value / 得分，视工具而定（0–1） |
| 6 | `Strand` | char | `+` / `-` / `*`（工具无法确定链时为 `*`） |
| 7 | `mod_ratio` | float | 修饰化学计量比（被修饰 reads 占总 reads 的比例，0–1） |

> 所有文件第一行为表头：`Chr\tStart\tEnd\tStatus\tProb\tStrand\tmod_ratio`

---

## 字段的工具差异

`Prob` 列的实际含义随工具不同：

| 工具 | `Prob` 列含义 | `mod_ratio` 含义 |
|---|---|---|
| **CHEUI** | 二分类模型的 m6A 概率（0–1） | 化学计量比（0–1） |
| **m6Anet** | m6A 概率（0–1） | 化学计量比 |
| **ELIGOS2** | adj P-value（越小越显著） | N/A（填 1.0） |
| **Epinano** | 原始 delta sum error 或 SVM 决策得分 | N/A（填 1.0） |
| **DENA** | LSTM 输出概率 | 化学计量比 |
| **MINES** | 随机森林 / 贝叶斯得分 | 化学计量比 |
| **Nanocompore** | GMM p-value（越小越显著） | log odds ratio 效应量 |
| **xPore** | 差异修饰概率 | 化学计量比差异 |
| **yanocomp** | GMM p-value | 效应量 |
| **NanoSPA** | 概率 / 评分 | 化学计量比 |
| **DRUMMER** | OR 校正 p-value | 比例差 frac_diff |

> 提示：做跨工具比较时，按 `mod_ratio` 而非 `Prob` 排序更有意义。

---

## 输出文件命名约定

### 单样本工具

```
results/<Tool>/<sample>/<sample>_<Tool>_processed.txt
```

示例：

```
results/CHEUI/treatment_rep1/treatment_rep1_CHEUI_processed.txt
results/ELIGOS2_solo/treatment_rep1/treatment_rep1_ELIGOS2_solo_processed.txt
```

### 对比型工具

```
results/<Tool>/<case>_vs_<control>/<case>_vs_<control>_<Tool>_processed.txt
```

示例：

```
results/Nanocompore/treatment_rep1_vs_control_rep1/treatment_rep1_vs_control_rep1_Nanocompore_processed.txt
```

---

## 工具对比文件格式

`results/summary/` 目录下会额外生成：

| 文件 | 说明 |
|---|---|
| `modification_summary.tsv` | 每个工具报告的位点总数、平均化学计量比、最高概率位点 |
| `tool_comparison.tsv` | 位点级别上的工具一致性矩阵（每个位点被多少工具报告） |
| `liftover_summary.tsv` | 若启用坐标转换，transcript→genome 的汇总 |

---

## 坐标系统说明

RNAModBench 同时输出两套坐标：

### transcriptome 坐标（7 个工具）
- `Chr` 列是**转录本 ID**（如 `ENST00000367770`、`NM_000014.5`）
- `Start` / `End` 是相对转录本 5' 端的 0-based 偏移
- 工具：CHEUI、m6Anet、Nanocompore、DENA、MINES、xPore、yanocomp

### genome 坐标（4 个工具）
- `Chr` 列是**染色体名**（如 `chr1`、`1`、`MT`）
- `Start` / `End` 是标准 BED 0-based 坐标
- 工具：ELIGOS2、Epinano、NanoSPA、DRUMMER

### 坐标转换
如需把 transcriptome 坐标投影回 genome，执行：

```bash
python scripts/r2d_liftover.py \
    --input results/summary/modification_summary.tsv \
    --gtf  reference/genes.gtf \
    --tool r2d \
    --output results/summary/liftover_summary.tsv
```

---

## 下游分析示例

### 1. 过滤高置信位点

```python
import pandas as pd

df = pd.read_csv("results/CHEUI/sample1/sample1_CHEUI_processed.txt", sep="\t")
high_conf = df[(df["Prob"] > 0.9) & (df["mod_ratio"] > 0.2)]
high_conf.to_csv("sample1_cheui_highconf.tsv", sep="\t", index=False)
```

### 2. 求 3 个信号工具的交集

```python
cheui  = pd.read_csv("results/CHEUI/sample1/sample1_CHEUI_processed.txt", sep="\t")
m6anet = pd.read_csv("results/m6Anet/sample1/sample1_m6Anet_processed.txt", sep="\t")
nano   = pd.read_csv("results/Nanocompore/sample1_vs_ctl/sample1_vs_ctl_Nanocompore_processed.txt", sep="\t")

# 以 Chr + Start 做 key 求交集
cheui["key"] = cheui.Chr + ":" + cheui.Start.astype(str)
m6anet["key"] = m6anet.Chr + ":" + m6anet.Start.astype(str)
nano["key"] = nano.Chr + ":" + nano.Start.astype(str)

common = set(cheui.key) & set(m6anet.key) & set(nano.key)
print(f"三工具共有 {len(common)} 个修饰位点")
```

### 3. 把 TSV 转为 BED6 做基因组浏览器可视化

```bash
awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$4"|"$5"\t0\t"$6}' \
    results/CHEUI/sample1/sample1_CHEUI_processed.txt \
    > sample1_cheui.bed
bedSort sample1_cheui.bed sample1_cheui.sorted.bed
bgzip sample1_cheui.sorted.bed && tabix -p bed sample1_cheui.sorted.bed.gz
```

> 现在你可以把 `.gz` / `.tbi` 载入 IGV 进行可视化。

### 4. 统计 motif 富集

```bash
python scripts/extract_5mer.py \
    results/CHEUI/sample1/sample1_CHEUI_processed.txt \
    reference/transcriptome.fa \
    sample1_cheui_with_5mer.tsv
# 然后统计以 "A" 为中心的 5-mer 分布
```
