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

#: 2026-09-29 (R1-4 round): the three Human IVT runs come from a single library
#: registered under one BioSample, so the Role cell names that structure instead
#: of the bare "negative control" (the footnote carries the full declaration).
ROLE_FIX = {("Human", "HeLa, IVT", "RNA002"): "negative control (three runs of one library)"}

with OUT.open("w") as fh:
    fh.write("# Table S10: Datasets used in this benchmark and their replicate structure.\n")
    fh.write("Species\tDataset\tChemistry\tStudy\tSamples\tn\tRole\n")
    for r in rows:
        role = ROLE_FIX.get((r["species"], r["dataset"], r["chemistry"]), r["role"])
        fh.write("\t".join([r["species"], r["dataset"], r["chemistry"], r["source_study"],
                            r["samples"], r["n_independent_units"], role]) + "\n")
    fh.write("# note: Each row is one dataset (species x condition x chemistry). The two mouse "
             "studies are independent studies, not replicates, and are never pooled; the last "
             "column gives the number of independent sequencing units used in the site-level "
             "evaluation: biological replicates for Arabidopsis, the Curlcake constructs and "
             "the Human WT condition, and single units for mouse, the E. coli IVT control and "
             "the RNA004 datasets. The three Human IVT units are runs of a single library "
             "registered under one BioSample (SAMN22863220) and are therefore sequencing "
             "repeats rather than biological replicates of one another; because the Human WT "
             "and IVT libraries also come from different studies (SRP393373 and SRP428418), "
             "the Human WT-versus-IVT comparison is made across batches. Curlcake_IVT_rep2_partial is a "
             "depth-matched subset of the full run and is not counted as a replicate; the "
             "machine-readable sample registry deposited with the source code carries the "
             "per-sample provenance and accession of every row.\n")
print(f"TableS10: {len(rows)} rows")
