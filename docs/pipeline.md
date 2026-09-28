# The callset pipeline (`src/harmonisation/`)

Stages that turn each tool's own output format into the harmonised callsets under
`data/callsets/`, then the metrics behind the figures. Run in order with
`bash src/harmonisation/scripts/run_all.sh`; the header of that script explains why the order is
what it is (scope pruning must happen before the audits and the evaluations, so no
table can retain a callset that is later removed).

This directory holds no figure code. Everything that computes, assembles or gates a
published figure lives in that figure's own `figures/<figure>/src/` - see
[`figure_index.md`](figure_index.md). The stages below are here because their output
is a callset or an evaluation table rather than a panel; the two working analyses
whose output a figure re-uses are listed with their tables in
[`data_availability.md`](data_availability.md).

| stage | what it does | main output |
|---|---|---|
| `00_build_registry.py` | sample/tool registry, replicate and independence assignment | `metadata/sample_registry.csv` |
| `01_extract_callsets.py` | parse every tool's native output into the common schema | `callsets/<platform>/<species>/<group>/<mod>/<tool>/<sample>.tsv` |
| `01b_liftover_missing.py`, `02b_validate_liftover.py` | coordinate lift-over for the mouse reference, and its validation | `metadata/liftover_validation*.csv` |
| `32_impute_strand.py` | infer strand where a tool reports none (Nanom6A) | `metadata/strand_imputation.csv` |
| `03_annotate_callsets.py` | reference base, 5-mer context, DRACH membership, transcript annotation (Ensembl) | callset annotation columns |
| `11_scope_split.py` | drop tool × dataset combinations that are outside the study scope | `metadata/nonm6a_scope.csv` |
| `33_center_base_filter.py --apply` | hard filter: called base must be compatible with the modification and strand | `metadata/center_base_filter.csv` |
| `05_offset_audit.py` | coordinate-offset distribution per callset versus the reference | `data/evaluation/tables/offset_*.tsv` |
| `04_build_universe.py` | candidate universe per sample (skipped by default; ~15 GB) | `harmonisation/universe/` (not deposited) |
| `06_eval_m6a_glori.py` | precision/recall/F1/MCC against GLORI per match window, with bootstrap intervals; localisation curve | `m6a_glori_confusion.tsv`, `m6a_localization_curve.tsv` |
| `07_eval_controls.py` | IVT negative-control false-positive rates, Curlcake synthetic truth, purified-site and knock-down metrics | `controls_ivt_fpr.tsv`, `curlcake_truth.tsv`, `purified_sites.tsv`, `ko_kd_metrics.tsv` |
| `08_eval_nonm6a.py` | m5C / Ψ / m1Ψ / Nm / inosine panels | `hela_nonm6a.tsv`, `nonm6a_*.tsv` |
| `09_eval_rna004.py` | RNA004 chemistry evaluations and Dorado model scans | `rna004_*.tsv` |
| `12_nonm6a_summary.py` | known-site comparison for the non-m6A tools (Figure 7) | `nonm6a_fig7_summary.tsv` |
| `13_completeness_audit.py`, `14_legacy_coverage_audit.py` | grid completeness and comparison against the earlier assembly | `metadata/completeness_audit.csv` |
| `29_anchor_audit.py`, `30_pileup_call_filter_audit.py`, `31_persite_reference_audit.py` | anchor/offset correctness, no-call leakage, per-site reference check | `anchor_audit.tsv`, `pileup_call_filter_audit.tsv`, `persite_reference_audit.tsv` |
| `10_qc_reconcile.py` | reconciliation and the QC narrative | `data/evaluation/qc_report.md` |
| `export_figure_ready.py` | per-replicate, figure-grade metrics | `figure_ready_replicates.tsv` |
| `34_export_callsets.py` | the deposited layer: deduplicate overlapping-transcript coordinates, keep the strongest call, drop bookkeeping columns | `data/callsets/` |

`20`–`27` are the replicate-aware analysis stages built on top of the callsets: they
write the region models, the Guitar BED inputs and the working figures of the
supporting analyses in `analysis/`. The numbered stages from `28` on draw published
figures and therefore live with them - `figures/<figure>/src/`, listed per figure in
[`figure_index.md`](figure_index.md).

## Shared library (`common/`)

| module | role |
|---|---|
| `src/harmonisation/common/config.py` | all roots and the locked conventions (window, coverage floor, modification vocabulary, Dorado groups) |
| `src/harmonisation/common/registry.py`, `src/harmonisation/common/rawinfo.py` | sample/tool registry and resolution of each tool's source files |
| `parsers/` | one parser per tool output format |
| `extract`/`src/harmonisation/common/annotate.py`, `src/harmonisation/common/center.py`, `src/harmonisation/common/strand.py`, `src/harmonisation/common/liftover.py` | schema harmonisation, centre-base logic, strand inference, lift-over |
| `src/harmonisation/common/match.py`, `src/harmonisation/common/evaluation.py`, `src/harmonisation/common/metrics.py`, `consensus*.py` | window matching, TP/FP/FN/TN definitions, metric computation, replicate-aware consensus |
| `src/harmonisation/common/regionmodel.py`, `src/harmonisation/common/refs.py` | transcript region models (UTR/metagene) and reference handling |
| `src/harmonisation/common/manifest.py`, `src/harmonisation/common/io_utils.py` | provenance ledger writing and table I/O |
| `src/harmonisation/common/figstyle.py`, `src/harmonisation/common/pagelayout.py`, `src/harmonisation/common/panelpage.py` | typography, page budget and collision gates, panel composition |
