# 贡献指南

感谢你对 RNAModBench 的关注！无论你是要**报告 Bug**、**添加新工具**、
**改进文档** 还是**重构代码**，都欢迎提交 Issue 或 Pull Request。

## 目录

- [报告 Bug / 提功能需求](#报告-bug--提功能需求)
- [开发环境搭建](#开发环境搭建)
- [代码风格约定](#代码风格约定)
- [添加一个新工具](#添加一个新工具)
- [提交 PR 之前的自检清单](#提交-pr-之前的自检清单)
- [Commit message 约定](#commit-message-约定)

---

## 报告 Bug / 提功能需求

在 GitHub 上新建 Issue 时，请附上：

1. **错误文本** —— 直接粘贴（用 ` ``` ` 包裹）
2. **config.yaml 的样本部分** —— 避免泄露隐私，但工具列表和参考文件路径很重要
3. **Snakemake 版本** —— `snakemake --version`
4. **最小复现步骤** —— 例如 "我把 reference 目录改成 `ref` 就报错"

---

## 开发环境搭建

```bash
# 1. fork 本仓库 → clone 你自己的 fork
git clone git@github.com:<YOUR_USER>/RNAModBench.git
cd RNAModBench

# 2. 创建 conda 主环境
conda create -n rnamodbench-dev -c bioconda -c conda-forge \
    snakemake=7.32.4 python=3.10 pandas=1.5 numpy=1.24 \
    matplotlib=3.7 seaborn=0.12 biopython=1.81
conda activate rnamodbench-dev
pip install matplotlib-venn pre-commit black flake8

# 3. 安装 pre-commit hook
pre-commit install
```

---

## 代码风格约定

### Python（`scripts/*.py`）

- **格式化**：使用 [black](https://github.com/psf/black)（行宽 100）
  ```bash
  black --line-length 100 scripts/my_script.py
  ```
- **linting**：flake8，忽略 E501（行宽）、W503（运算符换行）
- **docstring**：每个公共函数/脚本顶部都需要 docstring。输入文件格式需明确说明列名、分隔符、是否含表头。
- **Type hints**：可选，推荐对核心函数添加。

### Snakefile

- 所有规则的 `threads`、`resources` 必须从 `config` 读取，不可硬编码
- 每个规则的 shell 命令前写一行注释说明它的意图
- 输出文件命名遵循：`results/<Tool>/<sample>/...`

### YAML

- 使用 2 空格缩进（不要用 Tab）
- 键名用下划线风格（`snake_case`）
- 每个新参数必须附带注释

---

## 添加一个新工具

以添加工具 `MyNewTool` 为例（假设它需要 transcriptome BAM + eventalign）：

### 1. 在 `config/config.yaml` 添加参数块

```yaml
mynewtool:
  model: "resources/MyNewTool/model.pkl"
  threads: 12
  prob_threshold: 0.5
```

并在 `tools:` 列表里加上 `MyNewTool`。

### 2. 在 `envs/mynewtool.yaml` 创建 conda 环境文件

```yaml
name: mynewtool
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.10
  - numpy=1.24
  - pip
  - pip:
    - mynewtool==1.0.0
```

### 3. 在 `Snakefile` 增加对应规则

```python
rule mynewtool_predict:
    input:
        bam="results/alignment/{sample}_transcriptome.bam",
        eventalign="results/eventalign/{sample}.eventalign.txt"
    output:
        "results/MyNewTool/{sample}/{sample}_MyNewTool_raw.txt"
    params:
        model=config["mynewtool"]["model"],
        threads=config["mynewtool"]["threads"]
    shell:
        "mynewtool predict --bam {input.bam} --eventalign {input.eventalign} "
        "--model {params.model} --threads {params.threads} --output {output}"
```

### 4. 在 `scripts/postprocess_mynewtool.py` 编写后处理脚本

**必须输出 7 列 TSV**（Chr、Start、End、Status、Prob、Strand、mod_ratio）。
参考 `scripts/postprocess_cheui.py` 的结构。

### 5. 在 `scripts/generate_summary.py` 中注册工具名

在 `read_tool_results` 的 `priority_list` 中添加：

```python
("MyNewTool", "MyNewTool"),
```

### 6. 更新文档

- `README.md` 的 "Features" 工具列表
- `docs/TOOLS_OVERVIEW.md` 的工具矩阵
- `docs/CONFIG_REFERENCE.md` 增加新的参数块

---

## 提交 PR 之前的自检清单

- [ ] `black --check scripts/` 全部通过
- [ ] `snakemake --dry-run --cores 1` 在最小数据集上能生成正确的 DAG
- [ ] 新增后处理脚本输出符合 **7 列 TSV** 标准格式
- [ ] 更新了 README / docs 中相应的说明
- [ ] 添加或更新了相关的 docstring

## Commit message 约定

使用 `type(scope): subject` 风格，type 可选：

- `feat` —— 新功能 / 新工具
- `fix` —— Bug 修复
- `docs` —— 文档改动
- `refactor` —— 代码重构（不含行为变更）
- `perf` —— 性能优化
- `test` —— 测试相关
- `chore` —— 构建 / CI / 依赖升级

示例：

```
feat(eligos2): add min_depth filter from config
fix(postprocess_cheui): handle zero-coverage transcripts
docs(readme): add hardware requirements table
```
