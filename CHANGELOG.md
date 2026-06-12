# Changelog

本文件记录 RNAModBench 的版本变更。格式遵循
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/)，版本号遵循
[Semantic Versioning](https://semver.org/spec/v2.0.0.html)。

## [Unreleased]

### Added

- `docs/` 目录：新增 `INSTALL.md`、`TUTORIAL.md`、`TOOLS_OVERVIEW.md`、
  `OUTPUT_FORMAT.md`、`CONFIG_REFERENCE.md`、`FAQ.md`
- 项目根：新增 `CONTRIBUTING.md`、`CHANGELOG.md`
- `README.md`：补充参考文件格式、样本命名、工具矩阵、坐标系统、
  硬件估算、结果解读、常见报错对照表、引用信息等章节

### Changed

- `config/config.yaml` 重写：每个关键参数添加注释，标注
  `[NOT YET USED]` 未接入的配置项
- `Snakefile` 顶部增加 config 接入状态表，便于开发者判断哪些字段生效
- 所有 `scripts/postprocess_*.py` 补充顶层 docstring，描述原始输入格式
  （列名、分隔符、是否含表头）
- `scripts/generate_summary.py` 的 `read_tool_results` 增加 DRUMMER 注册，
  并明确按 priority_list 匹配工具名
- `scripts/extract_5mer.py` 重构：合并重复函数定义、单次遍历 FASTA 构建
  序列字典、增加 5-mer 列名注释

### Fixed

- `scripts/extract_5mer.py`：`extract_5mer_context` 被定义 3 次的问题
- `scripts/generate_summary.py`：`read_tool_results` 被定义 2 次的问题
- `scripts/extract_5mer.py`：边界外 position 现在返回 `'NNNNN'` 而非空字符串

### Known Issues

- `create_transcriptome` / `qc.*` / `output.*` / `resources.*` / `logging.*`
  / `species.*` 在 config.yaml 中已声明，但 Snakefile 尚未接入
- `create_report.py` 的 quality metrics 为硬编码 "N/A"，尚未接入工具层的
  真实统计

---

## [1.0.0] — 2025-06-01

### Added

- 初版发布
- 集成 12 个 RNA 修饰检测工具：CHEUI、ELIGOS2、m6Anet、Nanocompore、DENA、
  Epinano、Epinano_DiffErr、DRUMMER、MINES、xPore、yanocomp、NanoSPA
- Snakemake 驱动的统一工作流
- 标准 7 列 TSV 输出格式
- 汇总模块：`scripts/generate_summary.py`
- HTML 报告：`scripts/create_report.py`
- Guitar 图：`scripts/create_guitar_plots.R`
- 深度分布图：`scripts/generate_depth_plots.R`
- 坐标转换：`scripts/r2d_liftover.py`（调用 R2Dtool）
- 5-mer 提取：`scripts/extract_5mer.py`
- 示例配置 `config/config.yaml`
- 最小化 `README.md`
