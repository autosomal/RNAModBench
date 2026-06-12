# RNAModBench — config.yaml Parameter Reference

This document provides a detailed reference for every key in
`config/config.yaml`. Keys marked **[NOT YET USED]** are declared in
the YAML but not wired into the current Snakefile — modifying them
has no effect on pipeline runs.

## Contents

- [Base configuration](#base-configuration)
- [Guppy basecalling](#guppy-basecalling)
- [Alignment and signal handling](#alignment-and-signal-handling)
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
- [Utility-tool parameters](#utility-tool-parameters)
- [Reference-file paths](#reference-file-paths)
- [Reserved parameters](#reserved-parameters)

---

## Base configuration

| Key | Type | Description |
|---|---|---|
| `samples` | `list[str]` | Sample directory names. The length determines how many times single-sample tools run; contrast tools pair entries in order. |
| `tools` | `list[str]` | Tools to run. Valid values: CHEUI, ELIGOS2, m6Anet, Nanocompore, DENA, Epinano, Epinano_DiffErr, DRUMMER, MINES, xPore, yanocomp, NanoSPA. |
| `data_dir` | `string` | Parent directory holding `<sample>/fast5/` — normally `"data"`. |
| `reference_dir` | `string` | Reference directory — normally `"reference"`. |
| `results_dir` | `string` | Output directory — normally `"results"`. |

---

## Guppy basecalling

| Key | Type | Description | Default |
|---|---|---|---|
| `guppy.config` | `string` | Guppy configuration file; RNA kits typically use `rna_r9.4.1_70bps_hac.cfg`. | `rna_r9.4.1_70bps_hac.cfg` |
| `guppy.num_callers` | `int` | Number of parallel workers. | `4` |
| `guppy.threads_per_caller` | `int` | CPU threads allocated to each worker. | `20` |
| `guppy.total_threads` | `int` | Upper bound on total threads (used to declare resources to Snakemake). | `80` |

---

## Alignment and signal handling

| Key | Type | Description | Default |
|---|---|---|---|
| `alignment.threads` | `int` | `minimap2` threads used for both transcriptome and genome alignment. | `40` |
| `nanopolish.threads` | `int` | `nanopolish eventalign` threads. Reduce to 10–16 for very large datasets to avoid memory exhaustion. | `40` |

---

## CHEUI

| Key | Type | Description | Default |
|---|---|---|---|
| `cheui.kmer_model` | `string` | Path to the k-mer model (CSV). | `resources/CHEUI/kmer_models/model_kmer.csv` |
| `cheui.model1` | `string` | Stage-1 (pre-screen) HDF5 model file. | `resources/CHEUI/.../model1.h5` |
| `cheui.model2` | `string` | Stage-2 (refined) HDF5 model file. | `resources/CHEUI/.../model2.h5` |
| `cheui.threads` | `int` | CPU threads for inference. | `40` |
| `cheui.prob_threshold` | `float` | Minimum stage-2 probability. Default `0.999` is very strict; relax to `0.9` for more sensitivity. | `0.999` |
| `cheui.ratio_threshold` | `float` | Minimum estimated stoichiometry. | `0.1` |

---

## ELIGOS2

| Key | Type | Description | Default |
|---|---|---|---|
| `eligos2.threads` | `int` | Parallel thread count (one worker per region). | `12` |
| `eligos2.max_depth` | `int` | Per-site maximum read depth. | `2000000` |
| `eligos2.min_depth` | `int` | Per-site minimum read depth — sites below this value are skipped. | `5` |
| `eligos2.padj_threshold` | `float` | BH-adjusted p-value upper bound. | `0.0001` |
| `eligos2.oddr_threshold` | `float` | Lower bound on odds ratio (> 1.2 implies treatment is more modified than control). | `1.2` |

---

## m6Anet

| Key | Type | Description | Default |
|---|---|---|---|
| `m6anet.dataprep_threads` | `int` | Threads for the event-reading pre-processing step. | `12` |
| `m6anet.inference_threads` | `int` | Threads for the model-inference step. | `40` |
| `m6anet.readcount_max` | `int` | Upper bound on reads per transcript (memory guard). | `2000000` |
| `m6anet.prob_threshold` | `float` | Minimum probability to call a site m6A. | `0.5` |
| `m6anet.ratio_threshold` | `float` | Minimum estimated stoichiometry. | `0.1` |

---

## Nanocompore

| Key | Type | Description | Default |
|---|---|---|---|
| `nanocompore.threads` | `int` | Parallel threads for sampcomp. | `20` |
| `nanocompore.min_coverage` | `int` | Minimum reads covering a site. | `5` |
| `nanocompore.min_ref_length` | `int` | Minimum length of the reference transcript. | `10` |
| `nanocompore.pvalue_threshold` | `float` | Upper bound on the GMM p-value. | `0.05` |
| `nanocompore.lor_threshold` | `float` | Lower bound on the absolute log odds ratio (effect size). | `0.5` |

---

## DENA

| Key | Type | Description | Default |
|---|---|---|---|
| `dena.motif` | `string` | Target motif in IUPAC degenerate notation (R = A/G, H = A/C/U). | `RRACH` |
| `dena.corr_grp` | `string` | Name of the Tombo-corrected signal group. | `RawGenomeCorrected_000` |
| `dena.windows` | `string` | Upstream and downstream window sizes around the target site. | `"2 2"` |
| `dena.processes` | `int` | Number of parallel processes. | `40` |
| `dena.model` | `string` | Path to the LSTM model directory. | `resources/DENA/DENA_LSTM_Model` |
| `dena.ratio_threshold` | `float` | Minimum estimated stoichiometry. | `0.1` |
| `dena.coverage_threshold` | `int` | Minimum read coverage. | `20` |

---

## Epinano / Epinano_DiffErr

| Key | Type | Description | Default |
|---|---|---|---|
| `epinano.threads` | `int` | Parallel threads. | `12` |
| `epinano.model` | `string` | Trained SVM model file. | `resources/Epinano/models/rrach.q3.mis3.del3.linear.dump` |
| `epinano.columns` | `string` | Feature columns (1-based, comma-separated). Default `"8,13,23"` correspond to quality / mismatch / deletion columns. | `"8,13,23"` |
| `epinano.delta_threshold` | `float` | Minimum delta-error difference between treatment and control (Epinano_DiffErr only). | `0.1` |

---

## DRUMMER

| Key | Type | Description | Default |
|---|---|---|---|
| `drummer.threads` | `int` | Parallel threads. | `12` |
| `drummer.max_depth` | `int` | Upper bound on per-site read depth. | `2000000` |
| `drummer.min_depth` | `int` | Lower bound on per-site read depth. | `5` |
| `drummer.pvalue_threshold` | `float` | Upper bound on the odds-ratio-adjusted p-value. | `0.05` |
| `drummer.frac_diff_threshold` | `float` | Lower bound on treatment-vs-control proportion difference. | `0.1` |

---

## MINES

| Key | Type | Description | Default |
|---|---|---|---|
| `mines.kmer_models` | `string` | Path to the k-mer model list file (`names.txt`). | `resources/MINES/Final_Models/names.txt` |
| `mines.coverage_threshold` | `int` | Minimum read coverage — sites below this threshold are skipped. | `20` |
| `mines.ratio_threshold` | `float` | Minimum fraction of modified reads. | `0.1` |

---

## xPore

| Key | Type | Description | Default |
|---|---|---|---|
| `xpore.dataprep_threads` | `int` | Threads for the pre-processing step. | `40` |
| `xpore.readcount_max` | `int` | Upper bound on reads per transcript (memory guard). | `2000000` |
| `xpore.prob_threshold` | `float` | Minimum differential-modification probability (higher is stricter). | `0.5` |

---

## yanocomp

| Key | Type | Description | Default |
|---|---|---|---|
| `yanocomp.prep_threads` | `int` | Threads for the event-align → HDF5 pre-processing step. | `8` |
| `yanocomp.test_threads` | `int` | Threads for the GMM testing step. | `8` |
| `yanocomp.n_components` | `int` | Number of Gaussian mixture components (typically 50; larger → finer but slower). | `50` |
| `yanocomp.fdr` | `float` | Benjamini–Hochberg false-discovery-rate control. | `0.05` |
| `yanocomp.pvalue_threshold` | `float` | Additional hard p-value filter. | `0.05` |

---

## NanoSPA

| Key | Type | Description | Default |
|---|---|---|---|
| `nanospa.prob_threshold` | `float` | Minimum modification probability (output filter). | `0.5` |

---

## Utility-tool parameters

| Key | Type | Description |
|---|---|---|
| `utilities.r2d_tool` | `string` | Absolute path to the R2Dtool executable. |

---

## Reference-file paths

| Key | Type | Description |
|---|---|---|
| `reference_files.genome` | `string` | Path to the genome FASTA. |
| `reference_files.transcriptome` | `string` | Path to the transcriptome FASTA. |
| `reference_files.genes_gtf` | `string` | Path to the Ensembl-style GTF gene model. |
| `reference_files.genes_bed` | `string` | Path to the BED gene annotation (used by ELIGOS2). |

---

## Reserved parameters

The following keys are declared but **[NOT YET USED]** by the current
Snakefile. They are placeholders for planned features.

| Key | Description |
|---|---|
| `create_transcriptome` | Auto-generate transcriptome.fa from genome + GTF. |
| `convert_gtf_to_bed` | Auto-convert GTF to BED. |
| `create_gene_annotations` | Auto-build derived gene annotations. |
| `qc.min_read_length` / `min_read_quality` / `min_coverage` / `max_coverage` | Automated read-level and coverage-level QC. |
| `output.bed_format` / `liftover_to_genome` / `include_5mer_context` | Output-format control (currently hard-coded to 7-column TSV). |
| `resources.max_memory` / `max_cpus` / `gpu_available` / `gpu_devices` | Compute-resource caps (currently driven by the Snakemake `--cores` flag). |
| `logging.level` / `save_intermediates` / `cleanup_temp_files` | Logging and debugging controls. |
| `species.name` / `codon_table` / `mitochondrial_genome` / `annotation_source` | Species-aware default tuning. |
