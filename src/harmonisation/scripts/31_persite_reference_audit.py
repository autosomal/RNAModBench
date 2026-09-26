#!/usr/bin/env python3
"""31 -- per-site reference audit report-only, filters nothing : was this per-site
feature table built against the genome its BAM was aligned to?

Motivation (2026-09-18)
-----------------------
``Epinano_Variants.py`` writes ``base = refseq[refpos]`` from the FASTA handed to
``-r``.  The three Arabidopsis fip37 tables were generated with a **human** FASTA:
their ``base`` column agreed with TAIR10 at 26.5 % of positions (chance level) and
no -2..+2 bp shift explained it, while a 60 bp stretch matched GRCh38 chr2 at the
identical index.  The DiffErr chain filters on that column, so the KD callset
became "positions where the human genome carries an A" -- a plausible-looking but
meaningless callset (GLORI +/-2 enrichment 2.19 % vs 15-26 % for the other tools
on the same sample).

This script samples every ``*per.site.csv`` under ``result/EpiNano_DiffErr/`` and
reports, per table:

* agreement of ``base`` with the sample's own genome at shift 0, plus the
  -2..+2 bp shift scan (a shift = coordinate bug; no shift = wrong reference);
* the best-matching genome when the sample's own genome does not explain the
  column (diagnosis aid: it names the species the table was actually built from);
* verdict ``ok`` / ``suspect_axis_offset`` / ``suspect_reference_mismatch``.

Report-only: nothing is filtered, no callset is touched.

Outputs
-------
harmonisation/evaluation/tables/persite_reference_audit.tsv

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/harmonisation/scripts/31_persite_reference_audit.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import GENOMES, RESULT_RNA002, TABLE_DIR  # noqa: E402
from common.io_utils import write_table  # noqa: E402
from common.manifest import Inventory, log_time, setup_logger  # noqa: E402
from common.persite import Genome, audit_table, parse_persite_row, sample_rows  # noqa: E402

#: directory-name prefix -> species key of config.GENOMES.  The EpiNano result
#: tree is organised by sample, and the sample name always starts with its
#: species (Arabidopsis_*, mES*, mESCs_*, HeLa_*, E_*, Curlcake_*).
_PREFIX_SPECIES = (
    ("Arabidopsis", "Arabidopsis"),
    ("mES", "Mouse"),
    ("HeLa", "Human"),
    ("E_", "E.coli"),
    ("Curlcake", "Curlcake"),
    ("RNA0", "Curlcake"),
)
#: probe order for the "which genome was it really built from" column.
_PROBE_ORDER = ("Arabidopsis", "Mouse", "Human", "E.coli", "Curlcake")


def species_of(path: Path) -> str | None:
    name = path.parent.name
    for prefix, species in _PREFIX_SPECIES:
        if name.startswith(prefix):
            return species
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=RESULT_RNA002 / "EpiNano_DiffErr",
                    help="tree holding the per-site tables")
    ap.add_argument("--offsets", type=int, default=4,
                    help="seek positions sampled per table (default 4)")
    args = ap.parse_args()

    logger = setup_logger("31_persite_reference_audit")
    inv = Inventory("31_persite_reference_audit")
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    #: skip the salvage/work areas (``_aux/``) and the per-contig intermediates
    #: (``<bam>.tmp/<bam>.<contig>.per.site.csv``) of a run in flight -- only the
    #: finished tables that a downstream chain can consume are audited.
    tables = sorted(p for p in args.root.rglob("*per.site.csv")
                    if p.is_file()
                    and "_aux" not in p.parts
                    and not any(part.endswith(".tmp") for part in p.parts))
    logger.info("per-site tables found: %d", len(tables))
    cache: dict[str, Genome] = {}
    rows = []
    with log_time(logger, "per-site reference audit"):
        for path in tables:
            species = species_of(path)
            rec = {"sample_dir": path.parent.name, "table": path.name,
                   "species": species or "", "n_sampled": 0, "frac_shift0": None,
                   "best_shift": None, "frac_at_best_shift": None,
                   "frac_ok_strand_aware": None, "verdict": "no_rows",
                   "best_genome_guess": "", "shift_scan": ""}
            if species and species in GENOMES:
                g = cache.setdefault(species, Genome(GENOMES[species]))
                out = audit_table(path, g, offsets=args.offsets)
                rec.update({k: out.get(k) for k in
                            ("n_sampled", "frac_shift0", "best_shift",
                             "frac_at_best_shift", "frac_ok_strand_aware", "verdict")})
                rec["shift_scan"] = ",".join(f"{k}:{v}" for k, v in out["shifts"].items())
                if rec["verdict"] != "ok":
                    # which genome would explain this column? (diagnosis only)
                    parsed = [r for r in (parse_persite_row(x)
                                          for x in sample_rows(path, args.offsets)) if r]
                    best, best_frac = "", -1.0
                    for other in _PROBE_ORDER:
                        if other == species or other not in GENOMES:
                            continue
                        og = cache.setdefault(other, Genome(GENOMES[other]))
                        frac = (sum(og.base(r["ref"], r["pos"]) == r["base"] for r in parsed)
                                / len(parsed)) if parsed else 0.0
                        if frac > best_frac:
                            best, best_frac = other, frac
                    rec["best_genome_guess"] = (f"{best or 'none'} ({best_frac:.3f})"
                                                if best else "")
            else:
                logger.warning("  %s: unknown species, skipped", path)
            rows.append(rec)
            if rec["verdict"] not in ("ok",):
                logger.warning("  %-28s %-34s %-26s frac_shift0=%s best=%s",
                               rec["sample_dir"], rec["table"], rec["verdict"],
                               rec["frac_shift0"], rec["best_genome_guess"])

    df = pd.DataFrame(rows)
    out_path = TABLE_DIR / "persite_reference_audit.tsv"
    write_table(df, out_path)
    inv.record(out_path, n_rows=len(df))
    inv.flush()
    if len(df):
        logger.info("verdicts:\n%s", df["verdict"].value_counts().to_string())
        bad = df[df["verdict"].str.startswith(("suspect", "no_rows"))]
        if len(bad):
            logger.warning("%d table(s) need re-generation:\n%s", len(bad),
                           bad[["sample_dir", "table", "verdict", "frac_shift0",
                                "best_genome_guess"]].to_string(index=False))
    logger.info("-> %s", out_path)
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
