#!/usr/bin/env python3
"""Write Table S10 (datasets and replicate structure) from the frozen Table1.tsv."""

# --- RNAModBench path bootstrap (added when this file was deposited) ----------
import os as _rb_os, pathlib as _rb_pl


def _rb_find(start):
    for p in (start, *start.parents):
        if (p / "RNAMOD_BENCH_ROOT").exists():
            return p
    return start


_RB = _rb_pl.Path(_rb_os.environ.get("RNAMODBENCH_ROOT") or _rb_find(_rb_pl.Path(__file__).resolve().parent))
_XB = _rb_pl.Path(_rb_os.environ.get("RNAMODBENCH_LOCAL") or (_RB / "_local"))
# --------------------------------------------------------------------------- #
import csv
from pathlib import Path

SRC = Path(str(_XB / "manuscript/manuscript_rev/_alignment_20260919/Table1.tsv"))
OUT = Path(__file__).resolve().parent / "tables" / "TableS10_datasets.tsv"

rows = list(csv.DictReader(SRC.open(), delimiter="\t"))
with OUT.open("w") as fh:
    fh.write("# Table S10: Datasets used in this benchmark and their replicate structure.\n")
    fh.write("Species\tDataset\tChemistry\tStudy\tSamples\tn\tRole\n")
    for r in rows:
        fh.write("\t".join([r["species"], r["dataset"], r["chemistry"], r["source_study"],
                            r["samples"], r["n_independent_units"], r["role"]]) + "\n")
    fh.write("# note: Each row is one dataset (species x condition x chemistry). The two mouse "
             "studies are independent studies, not replicates, and are never pooled; the last "
             "column gives the number of independent sequencing units used in the site-level "
             "evaluation (biological replicates for HeLa, Arabidopsis and Curlcake; single units "
             "for mouse, E. coli IVT and the RNA004 datasets). Curlcake_IVT_rep2_partial is a "
             "depth-matched subset of the full run and is not counted as a replicate; the "
             "machine-readable sample registry deposited with the source code carries the "
             "per-sample provenance and accession of every row.\n")
print(f"TableS10: {len(rows)} rows")
