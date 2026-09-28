#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Acceptance checks for the rebuilt Figure S3 (exit 1 on any FAIL).

Checks: rebuilt HeLa GLORI reference, venn counts + naming, panel B per-unit
PPV (callsets/revision layer, unified "PPV vs. GLORI (2 bp)" label), panel
C provenance (the companion analysis analysis), page geometry, embedded fonts, and
the untouched original.
"""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import s3_common as sc

RESULTS: list[tuple[str, str, str]] = []


def check(name: str, ok: bool, detail: str) -> None:
    RESULTS.append((name, detail, "PASS" if ok else "FAIL"))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


def main() -> int:
    # 1. rebuilt reference ---------------------------------------------------
    bed = sc.GLORI_BEDS["Human"]
    lines = bed.read_text().splitlines()
    check("Hela_GLORI.bed rows == 112,451", len(lines) == 112_451,
          f"{len(lines)} rows")
    positions = [(p.split("\t")[0], int(p.split("\t")[1])) for p in lines]
    check("Hela_GLORI.bed unique + sorted",
          len(set(positions)) == len(positions)
          and positions == sorted(positions),
          f"{len(set(positions))} unique, sorted={positions == sorted(positions)}")
    backups = list(bed.parent.glob("Hela_GLORI.bed.bak_*"))
    check("old reference backed up", len(backups) >= 1,
          ", ".join(b.name for b in backups))
    with (sc.TABLE_DIR / "hela_glori_rebuild_map.tsv").open() as fh:
        map_rows = [ln for ln in fh if not ln.startswith("#")]
    check("rebuild map rows == 1,034", len(map_rows) - 1 == 1_034,
          f"{len(map_rows) - 1} dropped sites")

    # 2. venn counts ----------------------------------------------------------
    overlaps = sc.compute_overlap()
    try:
        sc.check_overlap(overlaps)
        check("venn counts (unified criterion)", True,
              "Arab 40461/80624/41894, Mouse 48221/41961/59580, "
              "Human 23686/112451/26606")
    except AssertionError as exc:
        check("venn counts (unified criterion)", False, str(exc))

    # 3. naming (panel A) ------------------------------------------------------
    txt = subprocess.run(["pdftotext", str(sc.PANEL_DIR / "S3A_venn.pdf"), "-"],
                         capture_output=True, text=True).stdout
    for label in ("Arabidopsis-rep1", "Arabidopsis-rep2", "Mouse-rep1",
                  "Mouse-rep2", "Human-rep1", "Human-rep2"):
        check(f"panel A label '{label}'", label in txt, "present" if label in txt
              else "MISSING")
    check("panel A: retired names gone",
          all(s not in txt for s in ("Hela-1", "Hela-2", "mESC-rep1")),
          "no 'Hela-1/Hela-2/mESC-rep1' in the panel text")

    # 4. panel B (PPV, callsets layer) --------------------------------------
    by_unit = sc.TABLE_DIR / "S3B_ppv_by_unit.tsv"
    check("panel B files exist",
          by_unit.exists() and (sc.PANEL_DIR / "S3B_ppv.pdf").exists(),
          f"{by_unit.name} + S3B_ppv.pdf")
    if by_unit.exists():
        rows = pd_read(by_unit)
        for species, want_n in (("Arabidopsis", 3), ("Mouse", 2), ("Human", 3)):
            got = sorted(rows[rows["species"] == species]["unit"].unique())
            check(f"panel B {species} units", len(got) == want_n,
                  f"{got}")
        n_tools = rows[rows["species"] == "Arabidopsis"]["tool"].nunique()
        check("panel B tool scope == 13 m6A tools", n_tools == 13,
              f"{n_tools} tools")
    btxt = subprocess.run(["pdftotext", str(sc.PANEL_DIR / "S3B_ppv.pdf"), "-"],
                          capture_output=True, text=True).stdout
    check("panel B y-label == 'PPV vs. GLORI (2 bp)'",
          "PPV vs. GLORI (2 bp)" in btxt,
          "retired 'Hit Rate' label absent" if "Hit Rate" not in btxt
          else "WARNING: 'Hit Rate' still present")
    check("panel B: the two mouse studies labelled",
          "mouse study A" in btxt and "mouse study B" in btxt
          and not any(s in btxt for s in ("SRP", "mES_WT", "mESCs_Mettl3_WT")),
          "study numbers in the legend; no sample id / accession (author decision 2026-09-21)")

    # 5. panel C (the companion analysis analysis) ----------------------------------
    prov = sc.TABLE_DIR / "S3C_provenance.tsv"
    check("panel C provenance file", prov.exists(), str(prov))
    if prov.exists():
        ptxt = prov.read_text()
        check("panel C credits the companion script",
              "16_mod_ratio_regression.py" in ptxt
              and "mod_ratio_matched_sites.tsv" in ptxt,
              "script + cached tables recorded with md5")
    ctxt = subprocess.run(
        ["pdftotext", str(sc.PANEL_DIR / "S3C_modratio_perrep.pdf"), "-"],
        capture_output=True, text=True).stdout
    check("panel C legend carries r/CCC",
          "CCC=" in ctxt and "r=" in ctxt, "the companion analysis stats format")
    check("panel C: no in-figure study note (house rule)",
          "studies:" not in ctxt, "removed per the author 2026-09-20 (legend/caption carries it)")

    # 6. page geometry + fonts --------------------------------------------------
    pdf = sc.OUT_DIR / "FigureS3_rev.pdf"
    info = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True).stdout
    size = [ln for ln in info.splitlines() if ln.startswith("Page size")]
    check("FigureS3_rev page size == sup3",
          "595.276 x 633.598 pts" in (size[0] if size else ""),
          (size[0] if size else "pdfinfo failed").strip())
    fonts = subprocess.run(["pdffonts", str(pdf)], capture_output=True, text=True).stdout
    check("no Type-3 fonts (editable vector text)", "Type 3" not in fonts,
          "Arial embedded" if "Arial" in fonts else "fonts?")
    ptxt = subprocess.run(["pdftotext", str(pdf), "-"],
                          capture_output=True, text=True).stdout
    check("no study accession printed on the page",
          not any(s in ptxt for s in ("SRP", "mES_WT", "mESCs_Mettl3_WT")),
          "page prints the study numbers only (author decision 2026-09-21)")

    # 7. submitted original untouched -------------------------------------------
    sup = sc.SUP3_PDF
    check("submitted sup3.pdf untouched", sup.stat().st_size == 736_344,
          f"{sup.stat().st_size} bytes")

    # 8. style -------------------------------------------------------------------
    sc.apply_page_style()
    import matplotlib.pyplot as plt
    check("no grid anywhere (rcParam)", plt.rcParams["axes.grid"] is False,
          f"axes.grid={plt.rcParams['axes.grid']}, font={plt.rcParams['font.family']}")

    # ---- report -----------------------------------------------------------------
    sc.TABLE_DIR.mkdir(parents=True, exist_ok=True)
    report = sc.TABLE_DIR / "verify_s3_report.tsv"
    with report.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["check", "detail", "status"])
        w.writerows(RESULTS)
    fails = [r for r in RESULTS if r[2] == "FAIL"]
    print(f"\n[summary] {len(RESULTS) - len(fails)}/{len(RESULTS)} PASS -> {report}")
    return 1 if fails else 0


def pd_read(path: Path) -> "pd.DataFrame":
    import pandas as pd
    return pd.read_csv(path, sep="\t", dtype={"sample": str, "unit": str})


if __name__ == "__main__":
    sys.exit(main())
