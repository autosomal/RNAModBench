# RNAModBench — Standardised Output Format Specification

Every tool in RNAModBench is post-processed by a tool-specific script
(`scripts/postprocess_*.py`) into a single common **7-column
tab-separated** format. This uniform structure simplifies tool
comparison and downstream analysis.

## Contents

- [Standard 7-column layout](#standard-7-column-layout)
- [Per-tool semantics](#per-tool-semantics)
- [Output file naming convention](#output-file-naming-convention)
- [Tool-comparison files](#tool-comparison-files)
- [Coordinate systems](#coordinate-systems)
- [Downstream analysis examples](#downstream-analysis-examples)

---

## Standard 7-column layout

Every `*_processed.txt` file uses the same header row and column
ordering:

| Column # | Header     | Type   | Description |
|---|---|---|---|
| 1 | `Chr`      | string | Chromosome or transcript identifier (see "Coordinate systems" below) |
| 2 | `Start`    | int    | **BED 0-based** start coordinate of the modification site |
| 3 | `End`      | int    | **BED 0-based, exclusive** end coordinate (equals `Start + 1` for single-base modifications) |
| 4 | `Status`   | string | Always `Mod` for sites that pass the tool's significance filters |
| 5 | `Prob`     | float  | Score, probability, or (adjusted) p-value; semantics vary by tool — see the next table |
| 6 | `Strand`   | char   | `+` / `-` / `*` (strand-unknown or not reported) |
| 7 | `mod_ratio` | float | Estimated stoichiometry — the fraction of reads supporting modification at the site |

> The first line of every `*_processed.txt` file is the header row:
> `Chr\tStart\tEnd\tStatus\tProb\tStrand\tmod_ratio`.

---

## Per-tool semantics

The precise meaning of `Prob` and `mod_ratio` differs between tools.
Downstream analysts should interpret these columns in light of the
tool they originate from:

| Tool | `Prob` column | `mod_ratio` column |
|---|---|---|
| **CHEUI** | Stage-2 m6A probability (0–1) | Estimated stoichiometry |
| **m6Anet** | Site-level m6A probability (0–1) | Estimated stoichiometry |
| **ELIGOS2** | Adjusted p-value (smaller is more significant) | Not reported (set to 1.0) |
| **Epinano** | SVM decision score or delta sum-of-errors | Not reported (set to 1.0) |
| **DENA** | LSTM probability | Estimated stoichiometry |
| **MINES** | Random-forest / Bayesian score | Estimated stoichiometry |
| **Nanocompore** | GMM p-value (smaller is more significant) | Log odds ratio (effect size) |
| **xPore** | Differential modification probability | Stoichiometry difference |
| **yanocomp** | GMM p-value | Effect size |
| **NanoSPA** | Modification probability | Estimated stoichiometry |
| **DRUMMER** | Odds-ratio-adjusted p-value | Treatment-vs-control proportion difference |

> For cross-tool ranking, we recommend sorting by `mod_ratio`
> (the estimated stoichiometry or effect size) rather than `Prob`.

---

## Output file naming convention

### Single-sample tools

```
results/<Tool>/<sample>/<sample>_<Tool>_processed.txt
```

For example:

```
results/CHEUI/treatment_rep1/treatment_rep1_CHEUI_processed.txt
results/ELIGOS2_solo/treatment_rep1/treatment_rep1_ELIGOS2_solo_processed.txt
```

### Contrast (paired) tools

```
results/<Tool>/<case>_vs_<control>/<case>_vs_<control>_<Tool>_processed.txt
```

---

## Tool-comparison files

The `results/summary/` directory is produced by
`scripts/generate_summary.py` and contains:

| File | Description |
|---|---|
| `modification_summary.tsv` | Per-tool site counts, average stoichiometry, and the most significant site |
| `tool_comparison.tsv` | Site-level concordance matrix — how many tools independently detected each site |
| `liftover_summary.tsv` | (If enabled) Transcriptome-to-genome liftover summary |

---

## Coordinate systems

RNAModBench tools natively report calls in one of two coordinate
spaces. Column `Chr` records which space each call belongs to.

### Transcriptome coordinates (7 tools)

- `Chr` contains a **transcript ID** (e.g. `ENST00000367770`,
  `NM_000014.5`)
- `Start` / `End` are 0-based offsets from the 5′ end of the
  transcript
- Tools: CHEUI, m6Anet, Nanocompore, DENA, MINES, xPore, yanocomp

### Genome coordinates (4 tools)

- `Chr` contains a **chromosome name** (e.g. `chr1`, `1`, `MT`)
- `Start` / `End` are standard BED 0-based coordinates
- Tools: ELIGOS2, Epinano, NanoSPA, DRUMMER

### Coordinate liftover

To project transcriptome-coordinate calls into genome space, run:

```bash
python scripts/r2d_liftover.py \
    --input results/summary/modification_summary.tsv \
    --gtf  reference/genes.gtf \
    --tool r2d \
    --output results/summary/liftover_summary.tsv
```

---

## Downstream analysis examples

### 1. Filter high-confidence sites

```python
import pandas as pd

df = pd.read_csv("results/CHEUI/sample1/sample1_CHEUI_processed.txt",
                 sep="\t")
high_conf = df[(df["Prob"] > 0.9) & (df["mod_ratio"] > 0.2)]
high_conf.to_csv("sample1_cheui_highconf.tsv", sep="\t", index=False)
```

### 2. Find the intersection of three signal tools

```python
cheui  = pd.read_csv("results/CHEUI/sample1/sample1_CHEUI_processed.txt",
                     sep="\t")
m6anet = pd.read_csv("results/m6Anet/sample1/sample1_m6Anet_processed.txt",
                     sep="\t")
nano   = pd.read_csv("results/Nanocompore/sample1_vs_ctl/sample1_vs_ctl_Nanocompore_processed.txt",
                     sep="\t")

# Compute per-tool site sets keyed by (Chr, Start)
cheui["key"] = cheui.Chr + ":" + cheui.Start.astype(str)
m6anet["key"] = m6anet.Chr + ":" + m6anet.Start.astype(str)
nano["key"]   = nano.Chr + ":" + nano.Start.astype(str)

common = set(cheui.key) & set(m6anet.key) & set(nano.key)
print(f"Three-tool intersection: {len(common)} modification sites")
```

### 3. Convert to BED6 for genome-browser visualisation

```bash
awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$4"|"$5"\t0\t"$6}' \
    results/CHEUI/sample1/sample1_CHEUI_processed.txt \
    > sample1_cheui.bed
bedSort sample1_cheui.bed sample1_cheui.sorted.bed
bgzip sample1_cheui.sorted.bed && tabix -p bed sample1_cheui.sorted.bed.gz
```

### 4. Summarise k-mer context

```bash
python scripts/extract_5mer.py \
    results/CHEUI/sample1/sample1_CHEUI_processed.txt \
    reference/transcriptome.fa \
    sample1_cheui_with_5mer.tsv
# Then tally the distribution of 5-mers centred on the modified base.
```
