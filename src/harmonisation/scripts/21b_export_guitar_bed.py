#!/usr/bin/env python
"""Export replicate-merged BED files for the R Guitar package.

Reads ``harmonisation/callsets`` (the analysis-ready layer: BED core columns,
species-native chromosome names, one file per sample = per replicate) instead of
the wide ``callsets`` tables, because all this step needs is
chrom/start/end/strand plus the merge bookkeeping.

The published metagene was NOT replicate-aware.  Comparing
``$RNAMODBENCH_LOCAL/code/code/guitar/<Condition>_clean/<Tool>.bed`` with the per-replicate call
sets shows Arabidopsis used rep3 only, mouse used the mES_WT study only, and
HeLa used an undocumented union.  Guitar itself cannot merge replicates --
``stSampleNum`` only controls how many equidistant points are taken inside each
site interval -- so the merge happens here:

    bed/<merge>/<Condition>/<Tool>.bed            one file per tool
    bed/<merge>_pooled/<Condition>/ALL.bed        every tool of the library pooled
    bed/rep_<tag>/<Condition>/<Tool>.bed          a single replicate

``<merge>`` = union | majority | intersection over the group's *independent*
units (``common.consensus.quorum``; nested subsets dropped, one row per
sequencing unit).  Intervals are widened by 1 bp on each side, exactly as the
original ``Guitar_*.r`` scripts did.  Chromosome labels are rewritten to the
spelling of the GTF handed to Guitar.

Scope follows the manuscript: m6A on the three real transcriptomes, the
non-m6A chemistries of HeLa / Curlcake (Fig. 7D, Fig. S7), and RNA004
(Fig. 8D, Fig. S9).
"""
from __future__ import annotations

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
import argparse
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from pandas.errors import EmptyDataError

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from common.config import SITES_ROOT                       # noqa: E402
from common.consensus import quorum                        # noqa: E402
from common.match import fix_chromosome                    # noqa: E402

OUT = (_XB / "harmonisation/guitar_metagene")
CLEAN = (_RB / "data/callsets")
BED = OUT / "bed"
PAD = 1

REF_GTF = {
    "Arabidopsis": str(_XB / "reference/arabidopsis/Arabidopsis_thaliana.TAIR10.61.gtf"),
    "Mouse": str(_XB / "reference/GRCm39/ensembl/Mus_musculus.GRCm39.114.gtf"),
    "Human": str(_XB / "reference/GRCh38p14/ensembl112/Homo_sapiens.GRCh38.112.chr.gtf"),
    "Curlcake": str(_XB / "reference/curlcakes/Curlcake.gtf"),
    "E.coli": str(_XB / "reference/K_12/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.61.gtf"),
}
#: (platform, species, group, mod) blocks that feed a published GUITAR panel
TARGETS = [
    ("RNA002", "Arabidopsis", "Arabidopsis_WT", "m6A"),
    ("RNA002", "Arabidopsis", "Arabidopsis_KD", "m6A"),
    ("RNA002", "Mouse", "Mouse_WT", "m6A"),
    ("RNA002", "Mouse", "Mouse_KO", "m6A"),
    ("RNA002", "Human", "HeLa_WT", "m6A"),
    ("RNA002", "Human", "HeLa_IVT", "m6A"),
    # Fig. 7D / Fig. S7: non-m6A chemistries, HeLa WT vs IVT
    *[("RNA002", "Human", g, m) for g in ("HeLa_WT", "HeLa_IVT")
      for m in ("Psi", "m1Psi", "m5C", "Nm")],
    # Fig. 8D / Fig. S9: RNA004, HeLa WT vs IVT (+ Curlcake IVT null)
    *[("RNA004", "Human", g, m) for g in ("RNA004_HeLa_WT", "RNA004_HeLa_IVT")
      for m in ("m6A", "Psi", "m5C", "inosine")],
    ("RNA004", "Curlcake", "RNA004_Curlcake_IVT", "m6A"),
]


def gtf_chroms(gtf: str) -> dict[str, str]:
    """{normalised label: the GTF's own spelling} for one annotation."""
    out = subprocess.run(
        f"grep -v '^#' {gtf} | cut -f1 | sort -u", shell=True,
        capture_output=True, text=True, check=True).stdout.split()
    return {fix_chromosome(c).lower(): c for c in out}


def load_sites(path: Path) -> pd.DataFrame:
    """BED core of one replicate's call set."""
    try:
        d = pd.read_csv(path, sep="\t", usecols=["chrom", "start", "strand"])
    except (FileNotFoundError, EmptyDataError):
        return pd.DataFrame(columns=["chrom", "pos", "strand"])
    if d.empty:
        return pd.DataFrame(columns=["chrom", "pos", "strand"])
    return pd.DataFrame({"chrom": d.chrom.astype(str), "pos": d["start"].astype("int64"),
                         "strand": d.strand.astype(str)}).drop_duplicates()


def merge_units(frames: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    """{merge rule: sites kept by that rule, with n_units / units columns}."""
    tags = list(frames)
    stacked = pd.concat([f.assign(unit=t) for t, f in frames.items()],
                        ignore_index=True)
    g = stacked.groupby(["chrom", "pos", "strand"], sort=False)["unit"]
    counts = pd.DataFrame({"n_units": g.nunique(),
                           "units": g.apply(lambda s: "|".join(sorted(set(s))))}
                          ).reset_index()
    return {"union": counts,
            "majority": counts[counts.n_units >= quorum(len(tags))],
            "intersection": counts[counts.n_units >= len(tags)]}


def write_bed(df: pd.DataFrame, path: Path, chrom_map: dict[str, str]) -> int:
    key = df.chrom.map(fix_chromosome).str.lower()
    keep = key.isin(chrom_map)
    if not keep.any():
        return 0
    bed = pd.DataFrame({
        "chrom": key[keep].map(chrom_map).values,
        "start": (df.pos[keep] - PAD).values,
        "end": (df.pos[keep] + 1 + PAD).values,
        "name": ".", "score": 0, "strand": df.strand[keep].values})
    path.parent.mkdir(parents=True, exist_ok=True)
    bed.to_csv(path, sep="\t", header=False, index=False)
    return len(bed)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--groups", default=None, help="restrict to these dataset groups")
    ap.add_argument("--mods", default=None, help="restrict to these mod types")
    args = ap.parse_args()
    want_groups = set(args.groups.split(",")) if args.groups else None
    want_mods = set(args.mods.split(",")) if args.mods else None

    reg = pd.read_csv((_RB / "metadata/sample_registry.csv"), sep="\t",
                      dtype=str)
    reg = reg[reg.independence_class != "nested_subset"]
    reg = reg.drop_duplicates(["dataset_group", "sequencing_unit"])

    totals: dict[str, int] = defaultdict(int)
    chrom_cache: dict[str, dict[str, str]] = {}
    for platform, species, group, mod in TARGETS:
        if want_groups and group not in want_groups:
            continue
        if want_mods and mod not in want_mods:
            continue
        gdir = CLEAN / platform / species / group / mod
        if not gdir.is_dir():
            print(f"skip {group}/{mod}: no callsets directory")
            continue
        if species not in chrom_cache:
            chrom_cache[species] = gtf_chroms(REF_GTF[species])
        chrom_map = chrom_cache[species]
        units = reg[(reg.platform == platform) & (reg.dataset_group == group)]
        tools = sorted(p.name for p in gdir.iterdir() if p.is_dir())
        pooled: dict[str, list[pd.DataFrame]] = defaultdict(list)
        for tool in tools:
            frames = {}
            for r in units.to_dict("records"):
                d = load_sites(gdir / tool / f"{r['sample']}.tsv")
                if len(d):
                    frames[r["replicate_tag"]] = d
            if not frames:
                continue
            for rule, sub in merge_units(frames).items():
                totals[rule] += write_bed(sub, BED / platform / rule / group / mod
                                          / f"{tool}.bed", chrom_map)
                pooled[rule].append(sub[["chrom", "pos", "strand"]])
            for tag, d in frames.items():
                write_bed(d, BED / platform / f"rep_{tag}" / group / mod
                          / f"{tool}.bed", chrom_map)
        for rule, parts in pooled.items():
            totals[rule + "_pooled"] += write_bed(
                pd.concat(parts, ignore_index=True),
                BED / platform / f"{rule}_pooled" / group / mod / "ALL.bed",
                chrom_map)
        print(f"{platform} {group:<20} {mod:<8} tools={len(tools):>2} "
              f"units={len(units)}", flush=True)
    print("BED rows written:", dict(totals))
    print("tree:", BED)


if __name__ == "__main__":
    main()
