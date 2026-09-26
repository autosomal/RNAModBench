# `data/callsets/` — processed modification callsets

One TSV per tool × sample × modification:

```
<platform>/<species>/<dataset_group>/<modification>/<tool>/<sample>.tsv
RNA002 | RNA004
Human | Mouse | Arabidopsis | E.coli | Curlcake
m6A | m5C | Psi | m1Psi | Nm | inosine
```

406 files, 4,451,964 site rows. Read the column definitions — including what each
tool's `score` actually means — in
[`../../docs/sites_columns.md`](../../docs/sites_columns.md), and the per-file
provenance (source file, parser, fingerprint, measured row count) in
[`../../metadata/callsets_index.tsv`](../../metadata/callsets_index.tsv).

Quick look at one callset:

```bash
head -3 RNA002/Human/HeLa_WT/m6A/CHEUI_m6A/HeLa_WT1.tsv
```

Sites are BED: 0-based, half-open, one base per call, in each species' native
contig names. A tool × sample that ended with no retained calls is present as a
header-only file, so the file set covers the grid that was run without gaps.
