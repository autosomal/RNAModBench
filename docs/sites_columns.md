# Column definitions for `data/callsets/`

Each file is one tool's retained calls on one sample, tab-separated with a header.
The layout mirrors the extraction layer: `<platform>/<species>/<dataset_group>/
<modification>/<tool>/<sample>.tsv`.

## Core interval (BED)

| column | meaning |
|---|---|
| `chrom` | native contig name of that species (`chr1`…, `chromosome` for *E. coli*, `Curlcake1…4` for the synthetic constructs) |
| `start`, `end` | 0-based, half-open; `end = start + 1`, one base per call. `start` is the modified position |
| `strand` | `+` / `-`. Nanom6A reports no strand; its own read alignment was used to infer it (99.74 % agreement with the exon annotation) |

## Score columns

`score` and `score_type` always appear together - **never interpret a score
without its `score_type`**, because the quantity differs by tool.

| tool | `score_type` | reading it |
|---|---|---|
| CHEUI_m6A, CHEUI_m5C | `Prob` | CHEUI's conversion-layer probability is ≈ 1 everywhere and carries no discrimination; the informative quantity is in `mod_ratio` (see below) |
| m6Anet | `Prob` | `probability_modified`. Its `mod_ratio` column is a *different* value, so both are kept |
| NanoSPA_m6A, NanoSPA_psU, NanoMUD_psi, NanoMUD_m1psi, NanoPsu | `Prob` | modification probability, 0.5–1 |
| Nanom6A, MINES, NanoNm, DENA (`m6a_ratio`) | modification ratio | the tool's own stoichiometry estimate |
| `Dorado_*` | `percent_modified / 100` | the `modkit` fraction |
| xPore, yanocomp | `FDR` | smaller is more significant |
| ELIGOS2_diff, ELIGOS2_solo | `adjPval` | smaller is more significant |
| DRUMMER, Nanocompore | `Pvalue` | smaller is more significant |
| EpiNano_Error | `delta_sum_err` | error-rate difference; **larger** means more modified |

| column | meaning |
|---|---|
| `mod_ratio` | modification ratio. Present only where it differs numerically from `score` (CHEUI_m6A/m5C, m6Anet); elsewhere the tool reported one value and `score` holds it |
| `frac_diff` | xPore only: differential modification rate between knock-out and wild-type |
| `mod_label` | `Dorado` only: the modification label reported by `modkit` (e.g. `m6A`) |

## Support

| column | meaning |
|---|---|
| `coverage` | number of reads supporting the call. To make coverage comparable across tools on a platform it is taken from one alignment per platform: RNA002 from the Nanom6A sorted BAM, RNA004 from the RNA004 minimap2 BAM. NanoNm keeps its own self-reported coverage (`src_coverage`). Which BAM a given sample used is recorded in `metadata/samples.csv` (`coverage_source`) |
| `src_coverage` | `Dorado` only: the effective coverage `modkit` itself reports |

## Base and motif context

| column | meaning |
|---|---|
| `ref_base` | reference base at the called position, guaranteed consistent with the modification type and strand by a hard centre-base filter (stage `33`) |
| `five_mer_raw` | the 5-mer context around the site as read from the reference |
| `drach` | m6A only: whether the site falls in a DRACH motif |

## Provenance of the call position

| column | meaning |
|---|---|
| `center_status` | `ok` for every deposited row; `no_ref_base` is possible for Nm, which has no base constraint. After the stage-33 hard filter no other value occurs |
| `offset_flag` | non-empty only in the few files whose coordinates were corrected (for example the CHEUI +4 centre fix), so a corrected call remains identifiable |

## Empty callsets

A callset that ended up with zero calls is written as a **header-only file** rather
than omitted, so the file set corresponds exactly to the tool × sample grid that was
run. `metadata/callsets_index.tsv` gives the row counts.

## What is *not* a column

Identity is carried by the path (platform, species, dataset group, modification,
tool, sample). Also omitted because they are constant within a file or guaranteed
by the pipeline: the strand and coverage *source* (see `metadata/samples.csv`),
`dist_center` / `pos_center` (zero after the hard filter), and parser bookkeeping
notes. Fingerprinted provenance for every file is in
`metadata/callsets_index.tsv` and `metadata/file_inventory.csv`.
