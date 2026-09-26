#!/usr/bin/env python
"""Build (and cache) the transcript region model used by the metagene analysis.

One cache file per species x class under ``sites_v2/_refs/regionmodels/``:

  Arabidopsis  TAIR10.61 (Ensembl)      mrna / ncrna
  Mouse        GRCm39.114 (Ensembl)     mrna / ncrna
  Human        GRCh38.112 .chr (Ensembl) mrna / ncrna

The annotation is taken from ``config.GTF_EXON`` (strip the trailing ``.exon``),
so the region model, the universe and the Guitar R panels all sit on the same
gene model.  Project rule 2026-09-18: Ensembl only, GENCODE is banned.

  python scripts/20_build_region_model.py [--species Human,Mouse]
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common.config import GTF_EXON                         # noqa: E402
from common.regionmodel import MODEL_DIR, parse_gtf        # noqa: E402


def gtf_for(species: str) -> Path:
    p = GTF_EXON[species]
    return Path(str(p)[:-len(".exon")] if str(p).endswith(".exon") else p)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default=",".join(GTF_EXON),
                    help="comma-separated; default all species in config")
    args = ap.parse_args()

    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    for species in args.species.split(","):
        if species not in GTF_EXON:
            print(f"skip {species}: not in config.GTF_EXON")
            continue
        gtf = gtf_for(species)
        if not gtf.exists():
            print(f"skip {species}: {gtf} not found")
            continue
        for cls in ("mrna", "ncrna"):
            out = MODEL_DIR / f"{species}.{gtf.stem}.{cls}.regionmodel.pkl"
            if out.exists() and out.stat().st_mtime >= gtf.stat().st_mtime:
                print(f"[skip] {out.name} (up to date)", flush=True)
                continue
            t0 = time.time()
            idx = parse_gtf(gtf, tx_class=cls)
            idx.save(out)
            print(f"[ok] {species} {cls}: {idx.n_tx} transcripts, "
                  f"{idx.n_segments():,} segments, "
                  f"{len(idx.tables)} chrom/strand tables, "
                  f"{time.time() - t0:.0f}s -> {out.name}", flush=True)


if __name__ == "__main__":
    main()

