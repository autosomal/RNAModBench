"""Reference interval helpers: exon BEDs per species (cached).

The candidate-site universe is defined inside *annotated exon regions* with a
minimum read coverage (see ``04_build_universe.py``); exons are built once from
the project's own annotation files and cached under
``harmonisation/universe/_refs/``:

* Arabidopsis -- ``Arabidopsis_thaliana.TAIR10.61.gtf.exon`` (CSV, 1-based)
* Mouse       -- ``Mus_musculus.GRCm39.114.gtf.exon`` (CSV, 1-based)
* Human       -- ``Homo_sapiens.GRCh38.112.chr.gtf.exon`` (CSV, 1-based;
*               **Ensembl only -- the author banned GENCODE on 2026-09-18**)
* E. coli     -- whole chromosome (4.6 Mb)
* Curlcake    -- whole constructs (10.1 kb)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .config import GENOMES, REF_ROOT, UNIVERSE_ROOT

REFS_DIR = UNIVERSE_ROOT / "_refs"

EXON_CSV = {
    "Arabidopsis": REF_ROOT / "arabidopsis" / "Arabidopsis_thaliana.TAIR10.61.gtf.exon",
    "Mouse": REF_ROOT / "GRCm39" / "ensembl" / "Mus_musculus.GRCm39.114.gtf.exon",
    #: Ensembl 112 (house rule 2026-09-18: GENCODE is banned project-wide)
    "Human": REF_ROOT / "GRCh38p14" / "ensembl112" / "Homo_sapiens.GRCh38.112.chr.gtf.exon",
}


def _write_bed(path: Path, chrom: str, start0: int, end0: int) -> None:
    with open(path, "a") as fh:
        fh.write(f"{chrom}\t{start0}\t{end0}\n")


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals.sort()
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def build_exon_bed(species: str, logger=None) -> Path:
    """Return the cached, merged exon BED for a species (built on first use)."""
    REFS_DIR.mkdir(parents=True, exist_ok=True)
    out = REFS_DIR / f"{species}.exon.bed"
    if out.exists():
        return out

    tmp = out.with_suffix(".bed.tmp")
    if tmp.exists():
        tmp.unlink()

    if species in EXON_CSV:
        src = EXON_CSV[species]
        per_chrom: dict[str, list[tuple[int, int]]] = {}
        with open(src, "r") as fh:
            header = fh.readline()
            for line in fh:
                f = line.rstrip("\n").split(",")
                if len(f) < 5:
                    continue
                chrom, start, end = f[0], int(f[3]), int(f[4])
                per_chrom.setdefault(chrom, []).append((start - 1, end))  # 1-based -> BED
        with open(tmp, "w") as fh:
            for chrom, ivs in per_chrom.items():
                for s, e in _merge_intervals(ivs):
                    fh.write(f"{chrom}\t{s}\t{e}\n")
        n = sum(len(_merge_intervals(v)) for v in per_chrom.values())
        if logger:
            logger.info("[%s] exon bed: %d merged intervals from %s", species, n, src.name)
    else:
        # whole-genome / whole-construct universe
        fa = Path(str(GENOMES[species]) + ".fai")
        with open(tmp, "w") as fh, open(fa) as src_fh:
            for line in src_fh:
                f = line.split("\t")
                if len(f) >= 2:
                    fh.write(f"{f[0]}\t0\t{f[1]}\n")
        if logger:
            logger.info("[%s] universe regions = whole sequences", species)

    tmp.replace(out)
    return out


def bam_style_bed(species: str, bam_reference_map: dict[str, str], logger=None) -> Path:
    """Exon BED with chromosome labels renamed to match a BAM's references.

    ``bam_reference_map`` maps normalised (``chr``-lowercase) names to the BAM's
    own reference names, e.g. ``{'chr1': '1'}`` for the HeLa BAMs.
    """
    from .match import fix_chromosome

    base = build_exon_bed(species, logger)
    out = REFS_DIR / f"{species}.{'__'.join(sorted(set(bam_reference_map.values()))[:1])}.exon.bed"
    if out.exists() and out.stat().st_mtime >= base.stat().st_mtime:
        return out
    tmp = out.with_suffix(".bed.tmp")
    skipped = 0
    with open(base) as fin, open(tmp, "w") as fout:
        for line in fin:
            f = line.rstrip("\n").split("\t")
            ref = bam_reference_map.get(fix_chromosome(f[0]))
            if ref is None:
                skipped += 1
                continue
            fout.write(f"{ref}\t{f[1]}\t{f[2]}\n")
    tmp.replace(out)
    if logger and skipped:
        logger.warning("[%s] %d exon intervals dropped (chromosome absent from BAM)",
                       species, skipped)
    return out


def reference_lengths(fasta_path: Path) -> dict[str, int]:
    fai = Path(str(fasta_path) + ".fai")
    out: dict[str, int] = {}
    with open(fai) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 2:
                out[f[0]] = int(f[1])
    return out


_DRACH_FASTA_CACHE: dict[str, dict[str, np.ndarray]] = {}


def curlcake_drach_a(fasta: Path) -> dict[str, np.ndarray]:
    """All DRACH-A positions (both strands, 0-based) of the Curlcake constructs.

    Ground truth for the fully m6A-modified constructs; also used by the offset
    audit for Curlcake samples.
    """
    key = str(fasta)
    if key in _DRACH_FASTA_CACHE:
        return _DRACH_FASTA_CACHE[key]
    seqs: dict[str, str] = {}
    name, buf = None, []
    with open(fasta) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()

    out: dict[str, np.ndarray] = {}
    for n, s in seqs.items():
        arr = np.frombuffer(s.encode("ascii"), dtype=np.uint8)
        pos = []
        for i in range(2, len(s) - 2):
            plus = (arr[i - 2] in (65, 71, 84) and arr[i - 1] in (65, 71)
                    and arr[i] == 65 and arr[i + 1] == 67 and arr[i + 2] in (65, 67, 84))
            minus = (arr[i + 2] in (65, 71, 84) and arr[i + 1] in (65, 71)
                     and arr[i] == 84 and arr[i - 1] == 67 and arr[i - 2] in (65, 67, 84))
            if plus or minus:
                pos.append(i)
        out[n] = np.array(sorted(set(pos)), dtype=np.int64)
    _DRACH_FASTA_CACHE[key] = out
    return out
