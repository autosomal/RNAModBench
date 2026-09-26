#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Panel A of the rebuilt Figure S3: the three replicate-overlap venns.

Recomputes the GLORI replicate overlap for Arabidopsis, Mouse (mESC) and
Human (HeLa) under the project-unified criterion "modification ratio > 0.1 in
each of the two replicates" (user decision 2026-09-19), verifies the counts
against independently recomputed values, draws the panel in the submitted
sup3A style and writes a traceability table with all three historical
criteria per species.

Outputs
-------
04_revision_analysis/figS3_revision/panels/S3A_venn.pdf / .png
04_revision_analysis/figS3_revision/tables/S3A_overlap_counts.tsv
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import matplotlib.pyplot as plt

import s3_common as sc


def draw_panel_a(ax) -> None:
    """Draw the three venns into ``ax`` (top-down page coordinates)."""
    sc.label_axes_page(ax, y_top=0.0, height=185.0)
    overlaps = sc.compute_overlap()
    sc.check_overlap(overlaps)
    #: unified "species + repN" style (s3_common.REPLICATE_LABELS)
    labels = sc.REPLICATE_LABELS
    for cx, species in zip(sc.A_PANEL_CX, sc.SPECIES_ORDER):
        o = overlaps[species]
        sc.draw_venn_pair(
            ax,
            cx=cx,
            n1=o["only1"],
            n_inter=o["intersection"],
            n2=o["only2"],
            jaccard=o["jaccard"],
            label1=labels[species][0],
            label2=labels[species][1],
            title=species,
        )


def write_counts_tsv(overlaps: dict) -> Path:
    sc.TABLE_DIR.mkdir(parents=True, exist_ok=True)
    path = sc.TABLE_DIR / "S3A_overlap_counts.tsv"
    fields = [
        "species", "criterion", "rep1_sites", "rep2_sites", "only1",
        "only2", "intersection", "union", "jaccard",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for species in sc.SPECIES_ORDER:
            o = overlaps[species]
            w.writerow({
                "species": species,
                "criterion": "both replicates > 0.1 (unified 2026-09-19)",
                **{k: (f"{o[k]:.4f}" if k == "jaccard" else o[k])
                   for k in fields[2:]},
            })
        # traceability rows: the two historical criteria of the HeLa panel
        o = overlaps["Human"]
        w.writerow({
            "species": "Human",
            "criterion": "raw replicate intersection, no ratio filter (legacy figure)",
            "rep1_sites": o["rep1_sites"], "rep2_sites": o["rep2_sites"],
            "only1": o["rep1_sites"] - o["raw_intersection"],
            "only2": o["rep2_sites"] - o["raw_intersection"],
            "intersection": o["raw_intersection"],
            "union": o["rep1_sites"] + o["rep2_sites"] - o["raw_intersection"],
            "jaccard": "",
        })
        w.writerow({
            "species": "Human",
            "criterion": "mean NormeRatio > 0.1 (legacy Hela_GLORI.bed, backup "
                         "Hela_GLORI.bed.bak_20260919_063142)",
            "intersection": o["mean_gt01"],
        })
    return path


def main() -> int:
    sc.apply_page_style()
    overlaps = sc.compute_overlap()
    sc.check_overlap(overlaps)

    sc.PANEL_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(sc.PAGE_W / 72.0, 185.0 / 72.0))
    ax = fig.add_axes((0, 0, 1, 1))
    draw_panel_a(ax)
    for ext in ("pdf", "png"):
        fig.savefig(sc.PANEL_DIR / f"S3A_venn.{ext}", facecolor="white")
    plt.close(fig)
    print(f"[write] {sc.PANEL_DIR / 'S3A_venn.pdf'} (+png)")

    tsv = write_counts_tsv(overlaps)
    print(f"[write] {tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
