"""Per-site feature table helpers -- which reference does a ``base`` column follow?

Motivation (2026-09-18)
-----------------------
``Epinano_Variants.py`` builds its per-site table as ``base = refseq[refpos]``
with ``refseq`` fetched from the FASTA handed to ``-r`` (L80-84).  The DiffErr
chain then filters on that column, so a wrong ``-r`` (the three Arabidopsis
fip37 tables were built against a **human** FASTA) silently turns the callset
into "positions where the wrong genome happens to carry an A".  Nothing raised
an error -- the table simply disagreed with the genome, which is trivial to spot:

* sample a few thousand rows of the table,
* compare ``base`` with the genome base at the same 0-based position,
* if the agreement is low, scan small shifts -- a *shift* is a coordinate bug,
  no shift at all is a *wrong reference*.

These helpers implement exactly that, are read-only, and are shared by
``scripts/31_persite_reference_audit.py`` (pipeline audit) and
``code/epinano_refix/audit_per_site_reference.py`` (one-off evidence tool).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

#: shifts (bp) tried when the zero-shift agreement is poor.
SHIFTS: tuple[int, ...] = (-2, -1, 0, 1, 2)
#: agreement at the same 0-based position required to call a table healthy.
OK_MIN = 0.99
NUC = ("A", "C", "G", "T")
_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def load_fai(fasta: Path) -> dict[str, tuple[int, int, int, int]]:
    """``name -> (length, offset, linebases, linewidth)`` from an existing .fai.

    Refuses to build an index: the reference FASTAs live outside the project and
    must stay untouched.
    """
    fai = Path(str(fasta) + ".fai")
    if not fai.exists():
        raise FileNotFoundError(f"missing .fai for {fasta} (refusing to create one)")
    out: dict[str, tuple[int, int, int, int]] = {}
    with open(fai) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                out[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    return out


class Genome:
    """Byte-addressed single-base reader (read-only, no index writes)."""

    def __init__(self, fasta: Path) -> None:
        self.fasta = Path(fasta)
        self.fai = load_fai(self.fasta)
        self._fh = open(self.fasta, "rb")
        self._cache: dict[str, str] = {}
        self._alias: dict[str, str] = {}

    def seq(self, name: str) -> str:
        if name not in self._cache:
            length, offset, bpl, bpw = self.fai[name]
            self._fh.seek(offset)
            raw = self._fh.read(length + (length // bpl) + 2).replace(b"\n", b"")
            self._cache[name] = raw[:length].decode("ascii", "replace").upper()
        return self._cache[name]

    def base(self, ref: str, pos: int) -> str:
        name = self.name_for(ref)
        if name is None:
            return "N"
        seq = self.seq(name)
        return seq[pos] if 0 <= pos < len(seq) else "N"

    def name_for(self, ref: str) -> str | None:
        """Map a table/``callset`` chrom name onto a FASTA key (``1`` <-> ``chr1``)."""
        ref = str(ref)
        if ref in self.fai:
            return ref
        if ref in self._alias:
            return self._alias[ref]
        cands = [f"chr{ref}", ref[3:] if ref.lower().startswith("chr") else ref,
                 ref.upper(), ref.lower(), "Chromosome", "chromosome"]
        for c in cands:
            if c in self.fai:
                self._alias[ref] = c
                return c
        self._alias[ref] = None  # type: ignore[assignment]
        return None


def sample_rows(path: Path, offsets: int = 6, chunk: int = 1_500_000) -> list[list[str]]:
    """Row-sample a (possibly multi-GB) CSV by seeking to spread offsets + the head.

    Full scans of the multi-GB EpiNano tables are forbidden on this filesystem;
    a few MB spread over the file is plenty for a 26 %-vs-100 % decision.
    """
    size = path.stat().st_size
    rows: list[list[str]] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        fh.readline()  # header
        for i in range(offsets):
            off = int(size * i / offsets) if i else 0
            fh.seek(off)
            lines = fh.read(chunk).split("\n")
            if off:
                lines = lines[1:]  # first line is partial
            rows.extend(r.split(",") for r in lines if r.count(",") >= 9)
    return rows


def parse_persite_row(row: list[str]) -> dict | None:
    """``#Ref,pos,strand,base,cov,q_mean,q_median,q_std,mat,mis,ins,del`` -> dict."""
    if len(row) < 11:
        return None
    try:
        return {"ref": row[0], "pos": int(row[1]), "strand": row[2],
                "base": row[3].upper(), "cov": int(row[4]),
                "mat": float(row[8]), "mis": float(row[9])}
    except ValueError:
        return None


def shift_scan(rows: list[dict], genome: Genome,
               shifts: tuple[int, ...] = SHIFTS) -> dict[str, float | None]:
    """Agreement between ``base`` and ``genome[ref, pos + k]`` for every shift k."""
    out: dict[str, float | None] = {}
    for k in shifts:
        hits = [genome.base(r["ref"], r["pos"] + k) == r["base"] for r in rows]
        out[str(k)] = round(float(np.mean(hits)), 4) if hits else None
    return out


def verdict(scan: dict[str, float | None]) -> str:
    """``ok`` / ``suspect_axis_offset`` / ``suspect_reference_mismatch``."""
    if not scan or scan.get("0") is None:
        return "no_rows"
    if scan["0"] >= OK_MIN:  # type: ignore[operator]
        return "ok"
    best = max((k for k in scan if scan[k] is not None), key=lambda k: scan[k])
    if int(best) != 0 and scan[best] >= OK_MIN:  # type: ignore[operator]
        return "suspect_axis_offset"
    return "suspect_reference_mismatch"


def strand_aware_frac(rows: list[dict], genome: Genome, mod: str = "m6A") -> float | None:
    """Share of rows on the modification-compatible base (m6A: A on '+', T on '-')."""
    target = {"m6A": "A", "m5C": "C", "Psi": "T", "m1Psi": "T", "inosine": "A"}.get(mod)
    if target is None or not rows:
        return None
    ok = 0
    for r in rows:
        b = genome.base(r["ref"], r["pos"])
        ok += int((b == target) if r["strand"] == "+" else (_COMPLEMENT.get(b) == target))
    return round(ok / len(rows), 4)


def audit_table(path: Path, genome: Genome, old: Path | None = None,
                offsets: int = 6) -> dict:
    """Full record for one per-site table (see module docstring)."""
    rec: dict = {"table": str(path), "n_sampled": 0, "shifts": {},
                 "best_shift": None, "frac_at_best_shift": None,
                 "frac_shift0": None, "frac_ok_strand_aware": None,
                 "verdict": "no_rows"}
    rows = [r for r in (parse_persite_row(x) for x in sample_rows(path, offsets)) if r]
    rec["n_sampled"] = len(rows)
    if not rows:
        return rec
    scan = shift_scan(rows, genome)
    rec["shifts"] = scan
    rec["frac_shift0"] = scan.get("0")
    rec["best_shift"] = int(max((k for k in scan if scan[k] is not None),
                                key=lambda k: scan[k])) if scan.get("0") is not None else None
    rec["frac_at_best_shift"] = scan.get(str(rec["best_shift"]))
    rec["verdict"] = verdict(scan)
    rec["frac_ok_strand_aware"] = strand_aware_frac(rows, genome)
    if old is not None and Path(old).exists():
        old_rows = [r for r in (parse_persite_row(x) for x in sample_rows(old, offsets)) if r]
        n = min(len(rows), len(old_rows))
        if n:
            rec["vs_old"] = {
                "n_compared": n,
                "frac_same_pos": round(sum(rows[i]["pos"] == old_rows[i]["pos"]
                                           for i in range(n)) / n, 4),
                "frac_same_cov": round(sum(rows[i]["cov"] == old_rows[i]["cov"]
                                           for i in range(n)) / n, 4),
                "frac_same_base": round(sum(rows[i]["base"] == old_rows[i]["base"]
                                            for i in range(n)) / n, 4),
            }
    return rec
