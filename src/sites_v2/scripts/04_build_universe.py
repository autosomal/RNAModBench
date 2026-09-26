#!/usr/bin/env python3
"""04 -- build the per-sample candidate-site universe (coverage-threshold scan)
and the Curlcake synthetic-truth table.

Candidate universe U (per sample, per modification type)
-------------------------------------------------------
U = { genomic positions p :
        p lies in an exon of the correct strand (strand map from the species
        GTF)  AND  coverage(p) >= c_min  AND  ref_base(p) matches the
        modification's expected base in transcript orientation }.

Strand awareness (critical, matches the strand-aware ``GLORI_all_A`` reference):
    m6A     -> genome 'A' on '+' , genome 'T' on '-'   (a minus-strand A sits
               on a genomic T)
    m5C     -> 'C' / 'G'
    Psi/m1Psi -> 'T' / 'A'
    inosine -> 'A' / 'T'
    Nm      -> any base (still restricted to exons)
All match coordinates are 0-based single-nucleotide, identical to the callset
and GLORI spaces (see common/config.py).  Coverage comes from the sample's own
genomic BAM (``samtools depth``, every mapped read, flag-filtered) so each
sample gets its own callable denominator.

The universe is stored once at the most permissive threshold in ``C_MIN_SCAN``
(min coverage), and ``universe_summary.csv`` reports the size at every scanned
threshold by filtering on the stored ``coverage`` column.

Curlcake is synthetic: no real BAM.  Its universe is the construct sequence
itself (all A / C / T / any positions) and its m6A truth is every adenosine in
the four constructs (the library is fully in-vitro methylated).

Usage
-----
conda activate benchmark-revision
python code/sites_v2/scripts/04_build_universe.py [--sample S ...]

When to re-run
--------------
A universe is built from the coverage BAM (``samtools depth``) + the annotation
only -- it never reads a callset.  So ``11_scope_split`` deleting callsets, or
any change confined to the callsets, cannot change a universe, and the step is
commented out in ``run_all.sh`` by default.  Re-run it only when the coverage
BAMs, the annotation or the coverage threshold change (mouse/human take hours).
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.annotate import (_as_str, _drach_mask, _seq_array, five_mer_matrix)
from common.bamcov import find_samtools
from common.config import (C_MIN_SCAN, GENOMES, GLORI, GTF_EXON, MOD_REF_BASE,
                           SAMPLES, SAMPLES_BY_NAME, UNIVERSE_ROOT, in_scope,
                           resolve_coverage_bam, rna004_coverage_bam, universe_path)
from common.io_utils import read_table, write_table
from common.manifest import Inventory, log_time, setup_logger
from common.match import exact_hit_mask, fix_chromosome

# transcript-orientation expected genome base, by modification and strand
EXPECTED_BASE = {
    "m6A":    {"+": ord("A"), "-": ord("T")},
    "m5C":    {"+": ord("C"), "-": ord("G")},
    "Psi":    {"+": ord("T"), "-": ord("A")},
    "m1Psi":  {"+": ord("T"), "-": ord("A")},
    "inosine": {"+": ord("A"), "-": ord("T")},
    "Nm":     None,  # any base
}
MIN_COV = min(C_MIN_SCAN)

_DRACH = __import__("re").compile(r"[AGT][AG]AC[ACT]")
_EXON_CACHE: dict[str, dict] = {}


# --------------------------------------------------------------------------- #
# exon strand index
# --------------------------------------------------------------------------- #
def _iter_exons(gtf_path: Path):
    """Yield (contig, start_1based, end_1based, strand) for every exon feature.

    Two input shapes are supported:
      * a precomputed exon CSV (``*.gtf.exon``, header
        ``contig,source,type,start,end,...``) -- used when ``GTF_EXON`` points
        at that file;
      * a full GTF (tab-separated, 3rd column == ``exon``) -- used as a fallback
        when only the ``.gtf`` exists (e.g. Human/gencode).
    """
    p = Path(gtf_path)
    if p.suffix == ".exon" or str(p).endswith(".gtf.exon"):
        with open(p) as fh:
            rdr = csv.reader(fh)
            try:
                next(rdr)  # header
            except StopIteration:
                return
            for row in rdr:
                if len(row) < 8 or row[2] != "exon":
                    continue
                yield row[0], int(row[3]), int(row[4]), row[6]
        return
    with open(p) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] != "exon":
                continue
            yield f[0], int(f[3]), int(f[4]), f[6]


def load_exon_index(species: str) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """chrom (normalised) -> (starts 0-based, ends 0-based, strands int8)."""
    if species in _EXON_CACHE:
        return _EXON_CACHE[species]
    gtf = GTF_EXON.get(species)
    if gtf is None or not gtf.exists():
        return {}
    by_chrom: dict[str, list] = defaultdict(list)
    for contig, start, end, strand in _iter_exons(gtf):
        s = 1 if strand == "+" else (-1 if strand == "-" else 0)
        if s == 0:
            continue
        by_chrom[fix_chromosome(contig)].append((start - 1, end - 1, s))
    out = {}
    for chrom, ex in by_chrom.items():
        ex.sort()
        starts = np.array([e[0] for e in ex], dtype=np.int64)
        ends = np.array([e[1] for e in ex], dtype=np.int64)
        strands = np.array([e[2] for e in ex], dtype=np.int8)
        out[chrom] = (starts, ends, strands)
    _EXON_CACHE[species] = out
    return out


def strand_at(idx, positions: np.ndarray) -> np.ndarray:
    """int8 strand (+1/-1/0) at each 0-based position."""
    starts, ends, strands = idx
    if len(starts) == 0:
        return np.zeros(len(positions), dtype=np.int8)
    i = np.searchsorted(starts, positions, side="right") - 1
    i = np.clip(i, 0, len(starts) - 1)
    in_exon = (positions >= starts[i]) & (positions <= ends[i])
    return np.where(in_exon, strands[i], 0).astype(np.int8)


# --------------------------------------------------------------------------- #
# coverage from BAM (covered positions only)
# --------------------------------------------------------------------------- #
def covered_positions(bam: Path, samtools: str, min_cov: int = MIN_COV):
    """{chrom: (0-based positions int64, depth int64)} for covered positions."""
    out: dict[str, tuple[list, list]] = defaultdict(lambda: ([], []))
    cmd = [samtools, "depth", "-q", "0", "-Q", "0", "-d", "0", str(bam)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    for line in res.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 3:
            continue
        try:
            dep = int(f[2])
        except ValueError:
            continue
        if dep < min_cov:
            continue
        out[fix_chromosome(f[0])][0].append(int(f[1]) - 1)
        out[fix_chromosome(f[0])][1].append(dep)
    return {c: (np.array(p, dtype=np.int64), np.array(d, dtype=np.int64))
            for c, (p, d) in out.items()}


def coverage_bam_for(sample) -> tuple[Path | None, str]:
    if sample.platform == "RNA004":
        p = rna004_coverage_bam(sample)
        return (p, "rna004_minimap2_bam") if p else (None, "missing")
    return resolve_coverage_bam(sample)


# --------------------------------------------------------------------------- #
# universe chunk for one chromosome / one sample
# --------------------------------------------------------------------------- #
def universe_chunk(arr: np.ndarray, eidx, pos: np.ndarray, depth: np.ndarray,
                   mod: str, chrom: str, glori) -> pd.DataFrame:
    strand = strand_at(eidx, pos)
    base = arr[np.clip(pos, 0, arr.size - 1)]
    exp = EXPECTED_BASE.get(mod)
    if exp is None:
        keep = strand != 0
    else:
        exp_pos = np.where(strand == 1, exp["+"],
                           np.where(strand == -1, exp["-"], 0)).astype(np.uint8)
        keep = (strand != 0) & (base == exp_pos)
    if not bool(keep.any()):
        return pd.DataFrame()
    pk = pos[keep]
    dk = depth[keep]
    sk = strand[keep]
    bk = base[keep]
    strand_str = np.where(sk == 1, "+", "-")
    mat, valid = five_mer_matrix(arr, pk, strand_str)
    five_mer = _as_str(mat, valid)
    is_drach = _drach_mask(mat, valid)
    if glori and glori.get(chrom) is not None:
        in_gl = exact_hit_mask(pk, glori[chrom])
    else:
        in_gl = np.zeros(len(pk), dtype=bool)
    base_str = bytes(bk.tobytes()).decode("latin-1").replace("\x00", "N")
    return pd.DataFrame({
        "chrom": chrom, "pos": pk.astype(np.int64),
        "base": list(base_str), "five_mer": five_mer,
        "is_drach": is_drach, "coverage": dk.astype(np.int64),
        "in_glori": in_gl,
    })


def mods_of(sample) -> list[str]:
    return [m for m in MOD_REF_BASE if in_scope(sample.platform, sample.species,
                                                sample.dataset_group, m)]


# --------------------------------------------------------------------------- #
# per-species driver (reuses loaded chromosome sequence across samples)
# --------------------------------------------------------------------------- #
def process_species(species: str, samples, logger, samtools: str,
                    summaries: list, pending: list) -> None:
    if species == "Curlcake":
        for s in samples:
            build_curlcake(s, summaries, pending)
        return
    genome = GENOMES.get(species)
    if genome is None or not genome.exists():
        for s in samples:
            pending.append({"sample": s.canonical, "reason": "no genome reference"})
        return
    exon_idx = load_exon_index(species)
    if not exon_idx:
        for s in samples:
            pending.append({"sample": s.canonical, "reason": "no exon GTF"})
        return
    glori = None
    if GLORI.get(species):
        from common.annotate import load_glori
        glori = load_glori(GLORI[species])
    fa = Fasta(str(genome), as_raw=True)
    fa_keys = list(fa.keys())

    # precompute covered positions per sample
    cov_cache = {}
    for s in samples:
        bam, src = coverage_bam_for(s)
        if bam is None:
            pending.append({"sample": s.canonical, "reason": f"no coverage BAM ({src})"})
            continue
        logger.info("[%s] samtools depth %s", s.canonical, Path(bam).name)
        cov_cache[s] = (covered_positions(bam, samtools), src, str(bam))

    acc = defaultdict(list)
    for chrom, eidx in exon_idx.items():
        fa_name = _match_key(fa_keys, chrom)
        if fa_name is None:
            continue
        arr = _seq_array(str(fa[fa_name]))
        for s, (cov, src, _bam) in cov_cache.items():
            if chrom not in cov:
                continue
            pos, depth = cov[chrom]
            for mod in mods_of(s):
                df = universe_chunk(arr, eidx, pos, depth, mod, chrom, glori)
                if len(df):
                    acc[(s, mod)].append(df)
    fa.close()

    for (s, mod), dfs in acc.items():
        out = pd.concat(dfs, ignore_index=True).sort_values(["chrom", "pos"],
                                                             kind="mergesort")
        write_table(out, universe_path(s, mod))
        cov = pd.to_numeric(out["coverage"], errors="coerce")
        summaries.append({
            "sample": s.canonical, "platform": s.platform, "species": species,
            "mod_type": mod, "source": cov_cache[s][1],
            "bam": cov_cache[s][2],
            "n_cov5": int((cov >= 5).sum()), "n_cov10": int((cov >= 10).sum()),
            "n_cov20": int((cov >= 20).sum()),
            "n_in_glori_cov10": int((out["in_glori"] & (cov >= 10)).sum()),
            "n_drach_cov10": int((out["is_drach"] & (cov >= 10)).sum()),
            "n_rows": len(out),
        })
        logger.info("  %-22s %-6s universe rows=%-9d cov10=%-9d in_glori=%-7d",
                    s.canonical, mod, len(out), (cov >= 10).sum(),
                    (out["in_glori"] & (cov >= 10)).sum())


def _match_key(fa_keys, chrom):
    from common.io_utils import normalize_chrom_for_fasta
    return normalize_chrom_for_fasta(fa_keys, chrom)


# --------------------------------------------------------------------------- #
# Curlcake (synthetic) universe + truth
# --------------------------------------------------------------------------- #
def build_curlcake(sample, summaries, pending) -> None:
    from common.annotate import load_glori
    cc_fasta = GENOMES["Curlcake"]
    fa = Fasta(str(cc_fasta), as_raw=True)
    # synthetic m6A truth = every adenosine in the four constructs
    truth_rows = []
    for construct in fa.keys():
        seq = str(fa[construct])
        arr = _seq_array(seq)
        a_pos = np.flatnonzero(arr == ord("A")).astype(np.int64)
        for p in a_pos:
            mat, valid = five_mer_matrix(arr, np.array([p]), np.array(["+"]))
            fm = _as_str(mat, valid & True)[0]
            truth_rows.append({"construct": construct, "chrom": fix_chromosome(construct),
                               "pos": int(p), "ref_base": "A", "five_mer": fm,
                               "is_drach": bool(_DRACH.search(fm))})
    truth = pd.DataFrame(truth_rows)
    UNIVERSE_ROOT.mkdir(parents=True, exist_ok=True)
    (UNIVERSE_ROOT / "Curlcake").mkdir(parents=True, exist_ok=True)
    write_table(truth, UNIVERSE_ROOT / "Curlcake" / "curlcake_synthetic_truth.tsv")
    summaries.append({"sample": "Curlcake", "platform": sample.platform,
                      "species": "Curlcake", "mod_type": "m6A_truth",
                      "source": "curlcake_A", "bam": "", "n_cov5": len(truth),
                      "n_cov10": len(truth), "n_cov20": len(truth),
                      "n_in_glori_cov10": 0, "n_drach_cov10": int(truth["is_drach"].sum()),
                      "n_rows": len(truth)})

    # per-mod candidate universe from the construct sequences
    mod_bases = {"m6A": "A", "m5C": "C", "Psi": "T", "m1Psi": "T", "Nm": None}
    for mod in mods_of(sample):
        target = mod_bases.get(mod)
        rows = []
        for construct in fa.keys():
            seq = str(fa[construct])
            arr = _seq_array(seq)
            if target is None:
                pos = np.arange(arr.size, dtype=np.int64)
            else:
                pos = np.flatnonzero(arr == ord(target)).astype(np.int64)
            for p in pos:
                mat, valid = five_mer_matrix(arr, np.array([p]), np.array(["+"]))
                fm = _as_str(mat, valid & True)[0]
                rows.append({"chrom": fix_chromosome(construct), "pos": int(p),
                             "base": chr(arr[p]), "five_mer": fm,
                             "is_drach": bool(_DRACH.search(fm)),
                             "coverage": -1, "in_glori": False})
        out = pd.DataFrame(rows).sort_values(["chrom", "pos"], kind="mergesort")
        write_table(out, universe_path(sample, mod))
        summaries.append({"sample": sample.canonical, "platform": sample.platform,
                          "species": "Curlcake", "mod_type": mod,
                          "source": "curlcake_construct", "bam": "",
                          "n_cov5": len(out), "n_cov10": len(out),
                          "n_cov20": len(out), "n_in_glori_cov10": 0,
                          "n_drach_cov10": int(out["is_drach"].sum()),
                          "n_rows": len(out)})
        print(f"  Curlcake {sample.canonical} {mod}: {len(out)} candidates")
    fa.close()


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", action="append")
    args = ap.parse_args()

    logger = setup_logger("04_build_universe")
    inv = Inventory("04_build_universe")
    samtools = find_samtools()
    if samtools is None:
        logger.error("samtools not found -- cannot build coverage-based universe")
        sys.exit(1)

    if args.sample:
        samples = [SAMPLES_BY_NAME[s] for s in args.sample if s in SAMPLES_BY_NAME]
    else:
        samples = SAMPLES

    by_species: dict[str, list] = defaultdict(list)
    for s in samples:
        by_species[s.species].append(s)

    summaries, pending = [], []
    with log_time(logger, "build_universe"):
        for species, sp_samples in by_species.items():
            logger.info("=== species %s (%d samples) ===", species, len(sp_samples))
            process_species(species, sp_samples, logger, samtools, summaries, pending)

    sdf = pd.DataFrame(summaries)
    UNIVERSE_ROOT.mkdir(parents=True, exist_ok=True)
    write_table(sdf, UNIVERSE_ROOT / "universe_summary.csv")
    inv.record(UNIVERSE_ROOT / "universe_summary.csv", n_rows=len(sdf))
    if pending:
        write_table(pd.DataFrame(pending), MANIFEST_PENDING)
        inv.record(MANIFEST_PENDING, n_rows=len(pending))
    inv.flush()
    logger.info("universe summary -> %s (%d rows, %d pending)",
                UNIVERSE_ROOT / "universe_summary.csv", len(sdf), len(pending))


from common.config import MANIFEST_DIR
MANIFEST_PENDING = MANIFEST_DIR / "universe_pending.csv"


if __name__ == "__main__":
    main()
