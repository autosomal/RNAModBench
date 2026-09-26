#!/usr/bin/env python3
"""29 -- anchor audit report-only, filters nothing : is every callset sitting on the
right base, *on the right strand*?

Motivation (2026-09-18): the four f5c-filled Nanom6A samples shipped calls whose
anchor was 7 bp upstream of the modified A (``--legacy-binary`` converter wrote
the window anchor ``indx`` instead of the A at ``indx+7``).  Nothing in the
pipeline complained, because the wrong coordinates are still *plausible*
coordinates.  This audit makes that failure mode loud.

Strict (chain-aware) verdict -- the one that matters
----------------------------------------------------
A call is *on base* when the reference base at ``pos_raw`` is the base the
modification sits on **in transcript orientation**:

    m6A -> genome A on '+', T on '-'   m5C -> C / G
    Psi, m1Psi -> T / A                inosine -> A / T     Nm -> any

The old loose check (``ref_base in {expected, complement}``) is still reported
as ``frac_loose`` because it is the only thing that can be computed without a
strand, but it is **not** a pass criterion any more: a callset whose strand
column is wrong for half of its rows scores 1.0 on the loose check and 0.5 on
the strict one.

Verdicts
--------
``ok``                        every call on base (or Nm: no expectation)
``off_base``                  some calls sit on a base the modification cannot
                              occupy -> below ``--min-strict`` (default 0.99)
``unknown_strand``            the strand is missing, so the call cannot be
                              judged (``32_impute_strand`` should have filled it)
``suspect_axis_offset``       < 0.90 strict **and** a constant shift k explains
                              the calls (>= 0.95 at k) -> coordinate offset
``suspect_reference_mismatch``  < 0.90 and no shift explains it -> the source was
                              built against another reference frame (the
                              Arabidopsis fip37 EpiNano per-site tables were
                              built with a human FASTA; 2026-09-18)
``no_reference_base``         modification without a single-base expectation (Nm)

Diagnostics that separate the two failure modes
-----------------------------------------------
``drach_failed`` / ``drach_failed_revcomp``: the DRACH share of the failing rows
as-is and reverse-complemented.  Going from ~0.0 to the callset's normal DRACH
share proves the coordinates are fine and only the *strand* was missing.

Outputs
-------
harmonisation/evaluation/tables/anchor_audit.tsv          all (sample, tool) rows
harmonisation/evaluation/tables/nanom6a_anchor_audit.tsv  Nanom6A subset

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/29_anchor_audit.py [--fail-on-bad]

``--fail-on-bad`` exits 1 when any callset still holds an off-base or
unjudgeable row (i.e. ``33_center_base_filter`` did not do its job).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.center import (OFF_BASE, OK, UNKNOWN_STRAND, centre_base_report,  # noqa: E402
                           expected_base_pair, strict_status)
from common.config import CALLSET_ROOT, GENOMES, TABLE_DIR  # noqa: E402
from common.io_utils import read_table, write_table  # noqa: E402
from common.manifest import Inventory, log_time, setup_logger  # noqa: E402
from common.persite import Genome  # noqa: E402
from common.persite import SHIFTS as PERSITE_SHIFTS  # noqa: E402

OK_MIN = 0.95
WARN_MIN = 0.90

#: lazy per-species genome readers (only touched for sub-threshold callsets)
_GENOME_CACHE: dict[str, Genome | None] = {}


def genome_for(species: str | None) -> Genome | None:
    """Genome reader for a species (``None`` when it has none / has no .fai)."""
    if not species:
        return None
    if species not in _GENOME_CACHE:
        path = GENOMES.get(species)
        try:
            _GENOME_CACHE[species] = Genome(path) if path else None
        except (FileNotFoundError, OSError):
            _GENOME_CACHE[species] = None
    return _GENOME_CACHE[species]


def classify_low_fraction(df, mod_type: str, species: str | None,
                          fallback: str) -> tuple[str, dict, int | None, float | None]:
    """Below the expected-base bar: constant offset or wrong reference?

    Scans the callset's own positions over -2..+2 bp: a *shift* that lifts the
    compatible-base share back to >= 0.95 is a coordinate bug; when no shift
    does, the calls live in another reference frame.
    """
    want = set(expected_base_pair(mod_type))
    genome = genome_for(species)
    if genome is None or not want or "pos_raw" not in df.columns:
        return fallback, {}, None, None
    chrom = df["chrom"].astype(str).to_numpy()
    pos = pd.to_numeric(df["pos_raw"], errors="coerce").to_numpy()
    keep = ~np.isnan(pos)
    scan: dict[str, float | None] = {}
    for k in PERSITE_SHIFTS:
        hits = [genome.base(c, int(p) + k) in want
                for c, p in zip(chrom[keep], pos[keep])]
        scan[str(k)] = round(float(np.mean(hits)), 4) if hits else None
    best = max((k for k in scan if scan[k] is not None),
               key=lambda k: scan[k], default=None)
    if best is None:
        return fallback, scan, None, None
    frac_best = scan[best]
    if int(best) != 0 and frac_best is not None and frac_best >= OK_MIN:
        return "suspect_axis_offset", scan, int(best), frac_best
    return "suspect_reference_mismatch", scan, int(best), frac_best


def audit_one(path: Path, meta: dict) -> dict:
    df = read_table(path)
    rec = {
        "platform": meta["platform"], "species": meta["species"],
        "dataset_group": meta["dataset_group"], "sample": meta["sample"],
        "tool": meta["tool"], "mod_type": meta["mod_type"],
        "n_calls": len(df), "source_file": "",
        "frac_strict": np.nan, "frac_loose": np.nan,
        "n_off_base": 0, "n_unknown_strand": 0,
        "frac_drach_raw": np.nan, "frac_dist_center_0": np.nan,
        "frac_dist_drach_0": np.nan,
        "drach_failed": np.nan, "drach_failed_revcomp": np.nan,
        "best_shift": None, "frac_at_best_shift": np.nan, "shift_scan": "",
        "anchor_flag": "empty",
    }
    if df.empty:
        return rec
    if "source_file" in df.columns:
        rec["source_file"] = str(df["source_file"].iloc[0])

    rep = centre_base_report(df, meta["mod_type"])
    code = strict_status(df["strand"].astype(str).to_numpy(),
                         df["ref_base"].astype(str).to_numpy(),
                         meta["mod_type"])
    rec["frac_strict"] = rep["frac_strict"]
    rec["frac_loose"] = rep["frac_loose"]
    rec["n_off_base"] = int((code == OFF_BASE).sum())
    rec["n_unknown_strand"] = int((code == UNKNOWN_STRAND).sum())
    rec["drach_failed"] = rep["drach_failed"]
    rec["drach_failed_revcomp"] = rep["drach_failed_revcomp"]

    for col, out in (("is_drach_raw", "frac_drach_raw"),
                     ("dist_center", "frac_dist_center_0"),
                     ("dist_drach_a", "frac_dist_drach_0")):
        if col in df.columns:
            ser = pd.to_numeric(df[col], errors="coerce")
            rec[out] = (float(df[col].astype(str).eq("True").mean()) if col == "is_drach_raw"
                        else float((ser == 0).mean()))

    if (code == 3).all():                       # Nm: no single-base expectation
        rec["anchor_flag"] = "no_reference_base"
        return rec
    strict = rep["frac_strict"]
    if np.isnan(strict) or rec["n_unknown_strand"]:
        rec["anchor_flag"] = "unknown_strand" if rec["n_unknown_strand"] >= rec["n_off_base"] \
            else "off_base"
        if rec["n_unknown_strand"] and rec["n_off_base"]:
            rec["anchor_flag"] = "off_base+unknown_strand"
        return rec
    if rec["n_off_base"] == 0:
        rec["anchor_flag"] = "ok"
        return rec
    rec["anchor_flag"] = "off_base"
    if strict >= WARN_MIN:
        return rec
    flag, scan, best, frac_best = classify_low_fraction(
        df, meta["mod_type"], meta.get("species"), "suspect_axis_offset")
    rec["anchor_flag"] = flag
    rec["best_shift"] = best
    rec["frac_at_best_shift"] = np.nan if frac_best is None else frac_best
    rec["shift_scan"] = ",".join(f"{k}:{v}" for k, v in scan.items())
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", action="append", default=None)
    ap.add_argument("--fail-on-bad", action="store_true",
                    help="exit 1 when any callset still holds an off-base or "
                         "unjudgeable row (acceptance mode, run after 33)")
    args = ap.parse_args()

    logger = setup_logger("29_anchor_audit")
    inv = Inventory("29_anchor_audit")
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    with log_time(logger, "anchor audit"):
        rows = []
        for f in sorted(CALLSET_ROOT.rglob("*.tsv")):
            parts = f.relative_to(CALLSET_ROOT).parts
            if len(parts) != 6:
                logger.warning("unexpected callset path: %s", f)
                continue
            platform, species, group, mod_type, tool, fname = parts
            sample = fname[:-4]
            if args.sample and sample not in args.sample:
                continue
            try:
                rec = audit_one(f, {"platform": platform, "species": species,
                                    "dataset_group": group, "mod_type": mod_type,
                                    "tool": tool, "sample": sample})
            except Exception as exc:  # noqa: BLE001 - reported, never silent
                logger.error("audit failed: %s (%s)", f, exc)
                continue
            rows.append(rec)
            if rec["anchor_flag"] not in ("ok", "no_reference_base", "empty"):
                logger.warning("  %-22s %-16s %-8s %s (strict=%.4f off=%d unk=%d)",
                               sample, tool, rec["anchor_flag"],
                               rec["source_file"].split("/")[-1],
                               rec["frac_strict"] if not np.isnan(rec["frac_strict"]) else -1,
                               rec["n_off_base"], rec["n_unknown_strand"])

    adf = pd.DataFrame(rows)
    write_table(adf, TABLE_DIR / "anchor_audit.tsv")
    inv.record(TABLE_DIR / "anchor_audit.tsv", n_rows=len(adf))
    n6 = adf[adf["tool"] == "Nanom6A"] if len(adf) else adf
    write_table(n6, TABLE_DIR / "nanom6a_anchor_audit.tsv")
    inv.record(TABLE_DIR / "nanom6a_anchor_audit.tsv", n_rows=len(n6))
    inv.flush()

    bad = pd.DataFrame()
    if len(adf):
        logger.info("verdicts:\n%s", adf["anchor_flag"].value_counts().to_string())
        bad = adf[~adf["anchor_flag"].isin(("ok", "no_reference_base", "empty"))]
        if len(bad):
            cols = ["sample", "tool", "mod_type", "anchor_flag", "frac_strict",
                    "frac_loose", "n_off_base", "n_unknown_strand",
                    "drach_failed", "drach_failed_revcomp", "best_shift"]
            cols = [c for c in cols if c in bad.columns]
            logger.warning("%d callset(s) not fully on base:\n%s", len(bad),
                           bad[cols].to_string(index=False))
        else:
            logger.info("every callset sits on its modification's base")
    logger.info("-> %s", TABLE_DIR / "anchor_audit.tsv")
    logger.info("log: %s", logger.log_path)

    if args.fail_on_bad and len(bad):
        logger.error("--fail-on-bad: %d callset(s) still carry off-base / "
                     "unjudgeable rows", len(bad))
        sys.exit(1)


if __name__ == "__main__":
    main()
