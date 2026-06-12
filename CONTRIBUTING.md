# Contributing to RNAModBench

Thank you for your interest in RNAModBench! Whether you want to
**report a bug**, **add a new tool**, **improve documentation**, or
**refactor code**, feel free to open an issue or submit a pull
request.

## Contents

- [Reporting a bug or requesting a feature](#reporting-a-bug-or-requesting-a-feature)
- [Setting up the development environment](#setting-up-the-development-environment)
- [Coding style conventions](#coding-style-conventions)
- [Adding a new tool](#adding-a-new-tool)
- [Pre-submission checklist](#pre-submission-checklist)
- [Commit message conventions](#commit-message-conventions)

---

## Reporting a bug or requesting a feature

When opening a new issue on GitHub, include:

1. **Error text** — paste directly, wrapped in triple backticks.
2. **Relevant sample of `config.yaml`** — redact private paths, but
   keep the `tools:` list and the reference-file paths intact.
3. **Snakemake version** — from `snakemake --version`.
4. **Minimal reproduction steps** — e.g., "when I rename the
   reference directory to `ref`, rule X fails with …".

---

## Setting up the development environment

```bash
# 1. fork the repository → clone your own fork
git clone git@github.com:<YOUR_USER>/RNAModBench.git
cd RNAModBench

# 2. create the main conda environment
conda create -n rnamodbench-dev -c bioconda -c conda-forge \
    snakemake=7.32.4 python=3.10 pandas=1.5 numpy=1.24 \
    matplotlib=3.7 seaborn=0.12 biopython=1.81
conda activate rnamodbench-dev
pip install matplotlib-venn pre-commit black flake8

# 3. install the pre-commit hook
pre-commit install
```

---

## Coding style conventions

### Python (`scripts/*.py`)

- **Formatting** — use [black](https://github.com/psf/black) with a
  100-character line width.
  ```bash
  black --line-length 100 scripts/my_script.py
  ```
- **Linting** — use flake8, ignoring E501 (line length) and W503
  (line break before binary operator).
- **Docstrings** — every public function / script needs a docstring
  describing the input file format (column names, delimiter, whether
  a header row is present).
- **Type hints** — optional; recommended for core functions.

### Snakefile

- `threads` and `resources` must be read from `config` — do not
  hard-code them.
- Precede each rule's shell command with a one-line comment describing
  the intent of the rule.
- Output files must follow the naming convention
  `results/<Tool>/<sample>/...`.

### YAML

- Indent with 2 spaces (no tabs).
- Keys use `snake_case`.
- Every new parameter must carry an inline or block comment.

---

## Adding a new tool

The example below adds a tool called `MyNewTool` (assume it requires
a transcriptome BAM + event-align output).

### 1. Add the parameter block to `config/config.yaml`

```yaml
mynewtool:
  model: "resources/MyNewTool/model.pkl"
  threads: 12
  prob_threshold: 0.5
```

Also add `MyNewTool` to the `tools:` list.

### 2. Create the conda environment file at `envs/mynewtool.yaml`

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

### 3. Add the rule to the `Snakefile`

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
        "mynewtool predict --bam {input.bam} "
        "--eventalign {input.eventalign} "
        "--model {params.model} --threads {params.threads} "
        "--output {output}"
```

### 4. Write the post-processing script at `scripts/postprocess_mynewtool.py`

The script **must emit a 7-column TSV** (Chr, Start, End, Status, Prob,
Strand, mod_ratio). Use `scripts/postprocess_cheui.py` as a structural
reference.

### 5. Register the tool name in `scripts/generate_summary.py`

Append to the priority list in `read_tool_results`:

```python
("MyNewTool", "MyNewTool"),
```

### 6. Update documentation

- "Features" list in `README.md`
- Tool matrix in `docs/TOOLS_OVERVIEW.md`
- Parameter reference in `docs/CONFIG_REFERENCE.md`

---

## Pre-submission checklist

- [ ] `black --check scripts/` passes.
- [ ] `snakemake --dry-run --cores 1` generates a valid DAG on a
  minimal dataset.
- [ ] The new post-processing script outputs the standard **7-column
  TSV** format.
- [ ] Related sections of the README / docs have been updated.
- [ ] Docstrings for new functions / scripts have been added.

---

## Commit message conventions

Use a `type(scope): subject` style, where `type` is one of:

- `feat` — new functionality / new tool.
- `fix` — bug fix.
- `docs` — documentation change.
- `refactor` — code restructure that does not change behaviour.
- `perf` — performance optimisation.
- `test` — test-related work.
- `chore` — build / CI / dependency upgrades.

Example:

```
feat(eligos2): add min_depth filter from config
fix(postprocess_cheui): handle zero-coverage transcripts
docs(readme): add hardware requirements table
```
