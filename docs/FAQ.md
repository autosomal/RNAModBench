# RNAModBench — Frequently Asked Questions

## Running the pipeline

**Q1. I only want to run three tools. Why do I still need the full conda environment installed?**

You do not. Snakemake only builds the conda environments required for
the tools listed in the `tools` list in `config.yaml`. If your list
contains only `CHEUI` and `ELIGOS2`, Snakemake creates only those
environments.

**Q2. `snakemake --dry-run` reports a `MissingInputException` for a sample file. What should I check?**

Common causes:

1. The `data/<sample_name>/fast5/` directory contains no `.fast5` files.
2. Directory names differ by case from the `samples:` list.
3. `reference/genome.fa` or `genes.gtf` are absent or misnamed.
4. A gunzipped reference FASTA is still compressed (`.fa.gz` is not
   accepted by minimap2).

**Q3. `nanopolish eventalign` runs very slowly / exhausts memory. What can I do?**

- Reduce `nanopolish.threads` in `config.yaml` to 10–16.
- Add a memory constraint to Snakemake: `--resources mem_mb=32000`.
- If the analysis is focused on m6A only, consider running the
  alignment-only tools (Epinano / ELIGOS2) instead, since they do not
  depend on event-align.

---

## Interpreting results

**Q4. A tool reports zero modification sites. Is that expected?**

Probably not. First check whether filter thresholds are too strict;
CHEUI ships with a default `prob_threshold` of 0.999, which is very
conservative. Relax it to 0.9 and re-run. If sites are still missing,
investigate whether the base-called FASTQ has sufficient quality and
whether the reference annotations (GTF) cover the regions of interest.

**Q5. The tools disagree substantially on site sets. Is this a bug?**

No — tools model different signals (raw ionic current, base-quality,
mismatch-rate, indel-rate, or alignment-distribution), so their calls
will not always coincide. RNAModBench is deliberately designed to
aggregate heterogeneous outputs. To find the most reliable sites,
filter `results/summary/tool_comparison.tsv` for rows with
`Tool_Count ≥ 3`.

**Q6. Why can't I directly compare transcriptome-coordinate tool results against genome-coordinate tool results?**

That is by design. Transcriptome-coordinate tools output transcript
IDs (e.g. `ENST00000367770`) and offsets; genome-coordinate tools
output chromosome names and positions. To merge them, run
`scripts/r2d_liftover.py` to project the transcriptome calls back into
genome space:

```bash
python scripts/r2d_liftover.py \
    --input results/summary/modification_summary.tsv \
    --gtf  reference/genes.gtf \
    --output results/summary/liftover_summary.tsv
```

---

## Installation and configuration

**Q7. Conda environment creation fails / is slow.**

- Try configuring a conda mirror with
  `conda config --add channels https://...`.
- Install `mamba` as a faster solver: `conda install -n base mamba`.

**Q8. The `Guitar` R package fails to install.**

```r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Guitar")
packageVersion("Guitar")  # expect ≥ 2.12.0
```

If dependencies are missing, install them explicitly:

```r
BiocManager::install(c("GenomicFeatures", "rtracklayer"))
```

---

## Data and species

**Q9. I only have base-called FASTQ (no FAST5). Can I still run anything?**

Yes — the alignment-level tools work without raw signal. Run
`ELIGOS2`, `Epinano`, `NanoSPA`, and `DRUMMER`. Omit signal-level
tools (CHEUI, m6Anet, Nanocompore, DENA, MINES, xPore, yanocomp).

**Q10. Can I run RNAModBench on non-human species?**

Yes. Replace `reference/genome.fa` and `genes.gtf` with your species'
reference files. The only hard requirement is matching chromosome
identifiers across the genome FASTA, transcriptome FASTA, and GTF.

---

## HPC / cluster integration

**Q11. How do I run RNAModBench on a cluster (SLURM / SGE)?**

Snakemake supports profiles natively. Create a profile YAML:

```yaml
# ~/.config/snakemake/slurm/config.yaml
jobs: 50
cluster: "sbatch --nodes=1 --ntasks={threads} --mem={resources.mem_mb}M -t {resources.runtime} -J {rule}"
default-resources: [mem_mb=4000, runtime=120]
```

Then launch with:

```bash
snakemake --profile slurm --use-conda
```

**Q12. How do I force re-run of a single rule without wiping results?**

```bash
snakemake --forceall --rerun-incomplete \
    results/CHEUI/sample1/sample1_CHEUI_processed.txt
```

Or mark existing output files as up-to-date without re-running them:

```bash
snakemake --touch --cores 1
```
