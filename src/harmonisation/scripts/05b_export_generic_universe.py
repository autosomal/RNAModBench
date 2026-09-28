#!/usr/bin/env python3
"""05b -- rebuild the legacy-format all-base universe ``<sample>__universe.tsv[.gz]``.

Why this script exists
----------------------
The evaluation layer (``06_eval_m6a_glori``, ``07_eval_controls``,
``08_eval_nonm6a``, ``09_eval_rna004``, ``12_nonm6a_summary``) reads
``universe/<platform>/<species>/<sample>__universe.tsv[.gz]`` -- a
**strand-agnostic, all-base** candidate set::

    chrom  pos  base  coverage  drach  in_glori

(``drach`` = strand of the DRACH motif at that position, ``+``/``-``/``.``;
coverage floor = ``min(C_MIN_SCAN)`` = 5, the most permissive stored threshold).

``04_build_universe.py`` today writes only the *strand-aware per-mod* files
(``<sample>__<mod>.tsv``); the producer of the generic name was removed in the
2026-09-15 generation, so those readers have silently kept using files built on
2026-09-14 -- for Human that is the in-house **GENCODE** annotation, which the
the author banned project-wide on 2026-09-18.  This script rebuilds the generic file
from exactly the same primitives as ``04`` (same exon index, same ``samtools
depth`` pass, same threshold, same ``in_glori`` flag) so that both file families
stay consistent for every species.

Verification
------------
``--verify-against`` compares the freshly built frame with an existing file of
the same sample (row count, coverage range, base composition, and set equality
on ``(chrom, pos, base, coverage, drach)``); run it once on a sample whose
annotation never changed (e.g. E. coli) before trusting it on Human.  The
``drach`` derivation (genome 5-mer in both orientations) reproduces the existing
files exactly (1.000 agreement; 2026-09-19).

Usage
-----
    conda run -n benchmark-revision --no-capture-output python \
        src/harmonisation/scripts/05b_export_generic_universe.py \
        --sample HeLa_WT1 --sample HeLa_WT2 ...

    # dry validation, writes into a scratch dir instead of universe/
    ... 05b_export_generic_universe.py --sample E_ss_rd_RNA1 \
        --out-dir figures/figure6/_validate \
        --verify-against <the existing file>

Curlcake is synthetic (no BAM) and is skipped: its generic file is the construct
sequence and is unaffected by the annotation change.
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
SCRIPTS = HERE.parent
sys.path.insert(0, str(SCRIPTS.parent))

from common.annotate import _drach_mask, five_mer_matrix  # noqa: E402
from common.config import (C_MIN_SCAN, GENOMES, GLORI, SAMPLES_BY_NAME,  # noqa: E402
                           UNIVERSE_ROOT)
from common.io_utils import read_table, write_table  # noqa: E402
from common.manifest import Inventory, setup_logger  # noqa: E402
from common.match import exact_hit_mask  # noqa: E402
from common.refs import EXON_CSV  # noqa: E402


def _load_04():
    """Import the companion ``04_build_universe.py`` (module name starts with a digit)."""
    spec = importlib.util.spec_from_file_location(
        "u04_build_universe", SCRIPTS / "04_build_universe.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


U04 = _load_04()
MIN_COV = min(C_MIN_SCAN)
COLUMNS = ["chrom", "pos", "base", "coverage", "drach", "in_glori"]


# --------------------------------------------------------------------------- #
def drach_strand(arr: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """'+' / '-' / '.' per position: DRACH on the forward or reverse strand."""
    n = pos.size
    plus = _drach_mask(*five_mer_matrix(arr, pos, np.full(n, "+", dtype="<U1")))
    minus = _drach_mask(*five_mer_matrix(arr, pos, np.full(n, "-", dtype="<U1")))
    out = np.full(n, ".", dtype="<U1")
    out[minus] = "-"
    out[plus] = "+"
    return out


def build_sample(sample, samtools: str, logger) -> pd.DataFrame:
    """All-base universe frame for one sample (same primitives as 04).

    Exon-model species (Arabidopsis / Mouse / Human) keep only positions inside an
    annotated exon of known strand (``04``'s ``universe_chunk`` with the "any base"
    modification ``Nm``); species without an exon model (E. coli -- and Curlcake,
    which this script skips) use the whole sequence, exactly as the README's
 universe definition prescribes (" E. coli and Curlcake using full-length sequences ").
    """
    species = sample.species
    whole_sequence = species not in EXON_CSV
    exon_idx = {} if whole_sequence else U04.load_exon_index(species)
    if not whole_sequence and not exon_idx:
        raise RuntimeError(f"{species}: no exon index")
    glori = None
    if GLORI.get(species):
        from common.annotate import load_glori
        glori = load_glori(GLORI[species])
    bam, src = U04.coverage_bam_for(sample)
    if bam is None:
        raise RuntimeError(f"{sample.canonical}: no coverage BAM ({src})")
    logger.info("[%s] samtools depth %s%s", sample.canonical, Path(bam).name,
                "  [whole sequence]" if whole_sequence else "")
    cov = U04.covered_positions(bam, samtools, MIN_COV)

    from pyfaidx import Fasta
    fa = Fasta(str(GENOMES[species]), as_raw=True)
    fa_keys = list(fa.keys())
    targets = sorted(cov) if whole_sequence else list(exon_idx)
    frames = []
    for chrom in targets:
        if chrom not in cov:
            continue
        fa_name = U04._match_key(fa_keys, chrom)
        if fa_name is None:
            continue
        arr = U04._seq_array(str(fa[fa_name]))
        pos, depth = cov[chrom]
        if whole_sequence:
            keep = (pos >= 0) & (pos < arr.size)
            pk, dk = pos[keep], depth[keep]
            if pk.size == 0:
                continue
            base = bytes(arr[pk].tobytes()).decode("latin-1").replace("\x00", "N")
            in_gl = (exact_hit_mask(pk, glori[chrom])
                     if glori and glori.get(chrom) is not None
                     else np.zeros(pk.size, dtype=bool))
            frames.append(pd.DataFrame({
                "chrom": chrom, "pos": pk, "base": list(base), "coverage": dk,
                "drach": drach_strand(arr, pk), "in_glori": in_gl,
            }))
            continue
        # mod "Nm" == "any base": keep = strand != 0  (04: EXPECTED_BASE["Nm"] is None)
        chunk = U04.universe_chunk(arr, exon_idx[chrom], pos, depth, "Nm", chrom, glori)
        if chunk.empty:
            continue
        frames.append(pd.DataFrame({
            "chrom": chrom,
            "pos": chunk["pos"].to_numpy(),
            "base": chunk["base"].to_numpy(),
            "coverage": chunk["coverage"].to_numpy(),
            "drach": drach_strand(arr, chunk["pos"].to_numpy()),
            "in_glori": chunk["in_glori"].to_numpy(),
        }))
    fa.close()
    if not frames:
        return pd.DataFrame(columns=COLUMNS)
    out = pd.concat(frames, ignore_index=True)
    return out.sort_values(["chrom", "pos"], kind="mergesort")[COLUMNS]


# --------------------------------------------------------------------------- #
def compare_frames(new: pd.DataFrame, old: pd.DataFrame, label: str,
                   logger) -> bool:
    """Row-set comparison between a rebuilt frame and the file on disk."""
    ok = True
    if len(new) != len(old):
        logger.warning("%s: row count differs: new=%d old=%d", label, len(new), len(old))
        ok = False
    key_new = new[["chrom", "pos"]].apply(tuple, axis=1)
    key_old = old[["chrom", "pos"]].apply(tuple, axis=1)
    set_new, set_old = set(key_new), set(key_old)
    if set_new != set_old:
        logger.warning("%s: (chrom,pos) sets differ: only-new=%d only-old=%d",
                       label, len(set_new - set_old), len(set_old - set_new))
        ok = False
    common = new.merge(old, on=["chrom", "pos"], how="inner",
                       suffixes=("_new", "_old"))
    if len(common):
        for col, name in (("base", "base"), ("coverage", "coverage"),
                          ("drach", "drach")):
            nc, oc = f"{col}_new", f"{col}_old"
            if nc not in common.columns:
                continue
            if col == "base":
                same = (common[nc] == common[oc]).mean()
            else:
                same = (common[nc].to_numpy() == common[oc].to_numpy()).mean()
            if same < 1.0:
                logger.warning("%s: %s mismatch on %d/%d shared rows",
                               label, name, int((1 - same) * len(common)), len(common))
                ok = False
        logger.info("%s: %d shared rows compared (base/coverage/drach)", label, len(common))
    logger.info("%s: %s", label, "IDENTICAL" if ok else "DIFFERS")
    return ok


def _existing_generic(sample) -> Path | None:
    base = UNIVERSE_ROOT / sample.platform / sample.species
    plain = base / f"{sample.canonical}__universe.tsv"
    gz = plain.with_suffix(".tsv.gz")
    for p in (plain, gz):
        if p.exists():
            return p
    return None


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", action="append", required=True,
                    help="canonical sample name (repeatable)")
    ap.add_argument("--out-dir", default=None,
                    help="override the universe root (dry runs / validation)")
    ap.add_argument("--verify-against", default=None,
                    help="compare the rebuilt frame with this file before/after writing")
    ap.add_argument("--min-cov", type=int, default=None)
    args = ap.parse_args()

    global MIN_COV
    if args.min_cov:
        MIN_COV = args.min_cov

    logger = setup_logger("05b_export_generic_universe")
    inv = Inventory("05b_export_generic_universe")
    samtools = U04.find_samtools()
    if samtools is None:
        logger.error("samtools not found")
        sys.exit(1)

    root = Path(args.out_dir) if args.out_dir else UNIVERSE_ROOT
    verify = Path(args.verify_against) if args.verify_against else None
    rows = []
    failures = []
    for name in args.sample:
        sample = SAMPLES_BY_NAME.get(name)
        if sample is None:
            logger.error("[%s] unknown sample name -- skipped", name)
            failures.append(name)
            continue
        if sample.species == "Curlcake":
            logger.info("[%s] synthetic -- generic file not rebuilt", name)
            continue
        frame = build_sample(sample, samtools, logger)
        out = root / sample.platform / sample.species / f"{sample.canonical}__universe.tsv.gz"
        out.parent.mkdir(parents=True, exist_ok=True)
        if verify is not None:
            old = pd.read_csv(verify, sep="\t")
            if not compare_frames(frame, old, f"[{name}] vs {verify.name}", logger):
                failures.append(name)
        frame.to_csv(out, sep="\t", index=False, compression="gzip")
        rows.append({
            "sample": sample.canonical, "platform": sample.platform,
            "species": sample.species, "min_cov": MIN_COV,
            "n_rows": len(frame),
            "n_base_A": int((frame.base == "A").sum()),
            "n_base_C": int((frame.base == "C").sum()),
            "n_base_G": int((frame.base == "G").sum()),
            "n_base_T": int((frame.base == "T").sum()),
            "n_drach": int((frame.drach != ".").sum()),
            "n_in_glori": int(frame.in_glori.sum()),
            "path": str(out), "written": datetime.now().isoformat(timespec="seconds"),
        })
        logger.info("[%s] wrote %s rows=%d", name, out, len(frame))
        inv.record(out, n_rows=len(frame))

    if rows:
        summary = root / "generic_universe_summary.tsv"
        new_df = pd.DataFrame(rows)
        if summary.exists():
            prev = read_table(summary)
            prev = prev[~prev["sample"].isin(new_df["sample"])]
            new_df = pd.concat([prev, new_df], ignore_index=True)
        write_table(new_df, summary)
        logger.info("summary -> %s", summary)
    inv.flush()
    if failures:
        logger.error("failures: %s", ",".join(failures))
        sys.exit(2)


if __name__ == "__main__":
    main()
