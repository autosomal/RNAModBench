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


#: could not tell which reading applied.  The two names are the two local splits of
#: ONE SRA run, as the sample registry says (`harmonisation/common/config.py`:
#: NON_INDEPENDENT_SAMPLE_OF["E_IVT_neg1"] = "E_IVT_neg2" with
#: NON_INDEPENDENT_REASON = "same_run_split"; SRR27228854, local halves of 882,018 and
#: 506,384 reads), so n stays 1 and the cell names the relationship.
SAMPLES_FIX = {("E.coli", "IVT", "RNA002"):
               "E_IVT_neg1, E_IVT_neg2 (two local splits of one run, SRR27228854)"}

with OUT.open("w") as fh:
    fh.write("# Table S10: Datasets used in this benchmark and their replicate structure.\n")
    fh.write("Species\tDataset\tChemistry\tStudy\tSamples\tn\tRole\n")
    for r in rows:
        role = ROLE_FIX.get((r["species"], r["dataset"], r["chemistry"]), r["role"])
        samples = SAMPLES_FIX.get((r["species"], r["dataset"], r["chemistry"]), r["samples"])
        fh.write("\t".join([r["species"], r["dataset"], r["chemistry"], r["source_study"],
                            samples, r["n_independent_units"], role]) + "\n")
    fh.write("# note: One row per dataset. The last column gives the number of independent sequencing units: "
"three biological replicates for Arabidopsis, the Curlcake constructs and the Human WT "
"condition, single units elsewhere. E_IVT_neg1 and E_IVT_neg2 are two local splits of one "
"run (SRR27228854) and count as one unit; the three Human IVT units are sequencing "
"repeats of one library (SAMN22863220), and the Human WT and IVT libraries come from "
"different studies (SRP393373, SRP428418), so that comparison is made across batches. "
"Curlcake_IVT_rep2_partial is a depth-matched subset of the full run, not a replicate.")
print(f"TableS10: {len(rows)} rows")
