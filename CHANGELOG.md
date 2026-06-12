# Changelog

All notable changes to RNAModBench will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and version numbers follow [Semantic Versioning](https://semver.org/).

---

## [Unreleased] — documentation and i18n pass

### Added

- Comprehensive Mermaid and Graphviz pipeline diagrams in
  `docs/FLOWCHART.md` (main workflow, post-processing, coordinate
  system overview, publication-ready DOT file).
- Publication-quality English rewrite of all documentation files:
  `docs/INSTALL.md`, `docs/TUTORIAL.md`,
  `docs/OUTPUT_FORMAT.md`, `docs/CONFIG_REFERENCE.md`,
  `docs/FAQ.md`, `docs/TOOLS_OVERVIEW.md`,
  `CONTRIBUTING.md`, and `CHANGELOG.md`.
- Detailed docstrings and input-format descriptions in
  `scripts/postprocess_*.py`, `scripts/extract_5mer.py`,
  `scripts/Epinano_DiffErr.R`, `scripts/create_guitar_plots.R`,
  and `scripts/generate_depth_plots.R`.
- `config.yaml` — rich inline comments for every parameter, with
  explicit `[NOT YET USED]` markers for reserved keys.
- `Snakefile` — header comment listing which config blocks are
  wired into rules and which remain reserved for future releases.

### Changed

- `README.md` rewritten in full scientific English.
- All Chinese-language docstrings and comments in Python and R
  scripts translated to scientific English.
- Unified naming conventions across the documentation (BED 0-based
  coordinates, transcriptome vs. genome coordinate systems, IUPAC
  motif notation).

---

## [1.1.0] — 2025-01-15

### Added

- Full English rewrite of `README.md`, all files under `docs/`,
  and all in-file documentation for `config.yaml` and the
  `Snakefile`.
- Publication-quality flowchart diagrams (Mermaid + Graphviz DOT)
  in `docs/FLOWCHART.md` and `docs/figures/pipeline.dot`.
- `docs/OUTPUT_FORMAT.md` — formal column-by-column specification
  of the unified 7-column TSV (Chr, Start, End, Status, Prob,
  Strand, mod_ratio) produced by every tool.
- `docs/CONFIG_REFERENCE.md` — full per-parameter reference for
  `config.yaml`.
- `docs/FAQ.md` — troubleshooting and HPC integration guidance.
- `CONTRIBUTING.md` — English-language contribution guide covering
  bug reporting, tool-integration instructions, and style
  conventions.

### Changed

- `scripts/postprocess_*.py` — module docstrings rewritten to
  describe native tool input format (column names, delimiter,
  presence/absence of a header row) in scientific English.
- `scripts/extract_5mer.py` — `extract_5mer_context` consolidated
  into a single definition; boundary positions now return `NNNNN`
  instead of an empty string.
- `scripts/generate_summary.py` — `read_tool_results` consolidated
  into a single definition with explicit tool-name priority list.
- `config/config.yaml` — parameter comments added; reserved
  parameters marked with `[NOT YET USED]`.
- `Snakefile` — header table describing which config blocks are
  wired into rules; inline comments clarified.

### Fixed

- `scripts/extract_5mer.py` — removed duplicate definitions of
  `extract_5mer_context` (previously defined three times).
- `scripts/generate_summary.py` — removed duplicate definition of
  `read_tool_results`.
- `config/config.yaml` — `tombo.basecall_grp` and
  `tombo.processes` now properly referenced from the Snakefile.

### Known issues

- `qc.*`, `output.*`, `resources.*`, `logging.*`, and
  `species.*` in `config.yaml` are declared but not yet wired into
  the Snakefile.
- Quality metrics in `scripts/create_report.py` are hard-coded to
  `N/A`; per-tool real statistics are not yet piped into the
  report.

---

## [1.0.0] — 2024-12-01

### Added

- Initial public release.
- 12 integrated RNA-modification detection tools: CHEUI, ELIGOS2,
  m6Anet, Nanocompore, DENA, Epinano, Epinano_DiffErr, DRUMMER,
  MINES, NanoSPA, xPore, yanocomp.
- Unified Snakemake-driven workflow.
- Standard 7-column TSV output format.
- `scripts/generate_summary.py` for aggregation and
  cross-tool comparison.
- `scripts/create_report.py` — automated HTML report.
- `scripts/create_guitar_plots.R` — gene-body metagene profiles.
- `scripts/generate_depth_plots.R` — per-sample depth-of-coverage
  QC.
- `scripts/r2d_liftover.py` — transcriptome-to-genome coordinate
  liftover (wraps R2Dtool).
- `scripts/extract_5mer.py` — k-mer context extraction around
  called modification sites.
- Sample configuration file at `config/config.yaml`.
- Minimal `README.md` with project overview.
