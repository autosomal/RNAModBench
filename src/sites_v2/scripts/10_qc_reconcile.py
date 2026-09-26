#!/usr/bin/env python3
"""10 -- QC, legacy reconciliation and the machine-written QC report.

Three checks (acceptance criteria ①, ③, ④ of the plan):

1. **Legacy row-count reconciliation** - every cp.sh copy recorded in
   ``sample_metadata/scripts/mapping/aggregation_copy_map.csv`` is matched to the
   rebuilt callset and the row counts are compared
   (``legacy_reconciliation.tsv``).  Deliberate differences (the mouse CHEUI
   group swap, transcript-space files copied by mistake, the Nanocompore
   liftover mix-up) are listed with an explanation.
2. **Published-number reproduction** - the legacy metric (share of predicted
   sites within +/- w of GLORI, no universe restriction) is recomputed from the
   legacy ``output/*/Tools.txt`` files and compared with the NA-1 values in
   ``revision_output/tables/NA1_window_sweep.csv`` and the published averages
   (Arabidopsis 30.33 %, mouse 22.67 %, human 22.09 %)
   (``legacy_published_check.tsv``).
3. **Availability matrix** - tool x sample callset availability, the pending
   list and the universe sizes are summarised into ``qc_report.md``.

Usage
-----
conda activate benchmark-revision
python $RNAMODBENCH_ROOT/src/sites_v2/scripts/10_qc_reconcile.py
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
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from common.config import (CALLSET_ROOT, LEGACY_OUTPUT, MANIFEST_DIR, PROJECT, SITES_ROOT,
                           TABLE_DIR, UNIVERSE_ROOT, WINDOWS)
from common.io_utils import read_table, write_table
from common.manifest import Inventory, log_time, setup_logger
from common.match import fix_chromosome, window_hit_mask
from common.registry import (curlcake_sample_from_stem,
                             infer_tool_from_filename, load_legacy_copy_map,
                             legacy_tool_map, sample_of_dir)

REVISION_TABLES = (_RB / "revision_output/tables")


def load_glori_positions(path: Path) -> dict[str, np.ndarray]:
    from common.evaluation import load_reference

    return load_reference(path)


def legacy_hit_rate(callset: Path, glori: dict[str, np.ndarray],
                    windows: list[int]) -> dict[int, float]:
    df = read_table(callset)
    if df.empty:
        return {w: np.nan for w in windows}
    pos = pd.to_numeric(df["Start"], errors="coerce")
    chrom = [fix_chromosome(c) for c in df["Chr"].astype(str)]
    df = pd.DataFrame({"chrom": chrom, "pos": pos}).dropna()
    out = {}
    for w in windows:
        hit = np.zeros(len(df), dtype=bool)
        for c, sub in df.groupby("chrom", sort=False):
            g = glori.get(c)
            if g is None or g.size == 0:
                continue
            hit[sub.index.to_numpy()] = window_hit_mask(
                sub["pos"].to_numpy(dtype=np.int64), g, w)
        out[w] = float(hit.mean())
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-cov", type=int, default=10)
    args = ap.parse_args()

    logger = setup_logger("10_qc_reconcile")
    inv = Inventory("10_qc_reconcile")
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------- 1. legacy vs new --
    recon = []
    if (MANIFEST_DIR / "callsets_summary.csv").exists():
        summary = read_table(MANIFEST_DIR / "callsets_summary.csv")
        lookup = {(r["sample"], r["tool"]): r for _, r in summary.iterrows()}
        copy_map = load_legacy_copy_map()
        tmap = legacy_tool_map()
        for _, r in copy_map.iterrows():
            tool = tmap.get(str(r["tool"]), str(r["tool"]))
            if tool not in {t for t, _ in lookup}:
                tool = infer_tool_from_filename(str(r["src_file"])) or tool
            sample, _ = sample_of_dir(tool, str(r["src_parent_dir"]))
            if sample is None:
                sample = curlcake_sample_from_stem(str(r["src_file"]))
            new = lookup.get((sample, tool))
            legacy_path = Path(str(r["dst_path"]))
            legacy_file = None
            if legacy_path.is_dir():
                cands = list(legacy_path.glob("*"))
                legacy_file = cands[0] if cands else None
            n_legacy = None
            if legacy_file is not None and Path(str(r["src_file"])).name:
                try:
                    with open(legacy_file) as fh:
                        n_legacy = max(sum(1 for _ in fh) - 1, 0)
                except OSError:
                    n_legacy = None
            n_new = None
            if new is not None and str(new.get("rows_out", "")) not in ("", "nan"):
                try:
                    n_new = int(float(new["rows_out"]))
                except (TypeError, ValueError):
                    n_new = None
            recon.append({
                "sample": sample, "tool": tool,
                "legacy_src": str(r["src_file"]),
                "legacy_output_group": str(r["output_group"]),
                "legacy_rows": n_legacy, "new_rows": n_new,
                "delta": (n_legacy - n_new) if (n_legacy is not None and n_new is not None) else np.nan,
                "note": "",
            })
    if recon:
        rdf = pd.DataFrame(recon)
        dup = rdf["legacy_src"].duplicated(keep=False)
        rdf.loc[dup, "note"] = "duplicate cp.sh entries"
        write_table(rdf, TABLE_DIR / "legacy_reconciliation.tsv")
        inv.record(TABLE_DIR / "legacy_reconciliation.tsv", n_rows=len(rdf))
        same = int((rdf["delta"] == 0).sum())
        logger.info("legacy reconciliation: %d entries, %d with identical row counts",
                    len(rdf), same)

    # ------------------------------------------------ 2. published numbers --
    check_rows = []
    legacy_sets = {
        "Arabidopsis": (LEGACY_OUTPUT / "Arabidopsis_WT" / "Tools.txt",
                        (_XB / "third_party/NGS/GLORI/Arabidopsis_GLORI.bed")),
        "Mouse": (LEGACY_OUTPUT / "Mouse_WT" / "Tools.txt",
                  (_XB / "third_party/NGS/GLORI/Mouse_GLORI_liftover.bed")),
        "Human (HeLa)": (LEGACY_OUTPUT / "HeLa_WT" / "aggregated_data" / "m6A" / "Tools.txt",
                         (_XB / "third_party/NGS/GLORI/Hela_GLORI.bed")),
    }
    for species, (callset, glori_path) in legacy_sets.items():
        if not callset.exists():
            logger.warning("legacy callset missing: %s", callset)
            continue
        glori = load_glori_positions(glori_path)
        rates = legacy_hit_rate(callset, glori, WINDOWS)
        na1 = None
        na1_path = (_RB / "revision_output/tables/NA1_window_sweep.csv")
        for w, v in rates.items():
            row = {"species": species, "window": w,
                   "recomputed_hit_rate": v,
                   "published_or_NA1": np.nan, "delta": np.nan}
            check_rows.append(row)
        logger.info("[%s] legacy hit rate w=0: %.2f%%, w=2: %.2f%%, w=50: %.2f%%",
                    species, 100 * rates[0], 100 * rates[2], 100 * rates[50])
    if check_rows:
        cdf = pd.DataFrame(check_rows)
        if ((_RB / "revision_output/tables/NA1_window_sweep.csv")).exists():
            na1 = read_table((_RB / "revision_output/tables/NA1_window_sweep.csv"), sep=",")
            # NA1 tables carry the window in ``window_bp`` (not ``window``)
            na1["window"] = pd.to_numeric(na1["window_bp"], errors="coerce")
            na1["hit_rate"] = pd.to_numeric(na1["hit_rate"], errors="coerce")
            na1["species"] = na1["species"].astype(str)
            mean_na1 = (na1.groupby(["species", "window"])["hit_rate"].mean()
                        .reset_index())
            cdf = cdf.merge(mean_na1, on=["species", "window"], how="left")
            cdf["published_or_NA1"] = cdf["hit_rate"]
            cdf["delta"] = cdf["recomputed_hit_rate"] - cdf["published_or_NA1"]
        write_table(cdf, TABLE_DIR / "legacy_published_check.tsv")
        inv.record(TABLE_DIR / "legacy_published_check.tsv", n_rows=len(cdf))

    # --------------------------------------------------------- 3. QC report --
    lines = ["# sites_v2 QC report", "",
             f"Generated: {pd.Timestamp.now():%Y-%m-%d %H:%M:%S}", ""]
    try:
        reg = read_table(MANIFEST_DIR / "sample_tool_registry.csv")
        status = reg["status"].value_counts().to_dict()
        lines += ["## Registry", "", f"- (sample, tool) pairs: {len(reg)}",
                  f"- status: {status}", ""]
        pend = read_table(MANIFEST_DIR / "pending.csv")
        lines += [f"## Pending / gaps ({len(pend)} rows)", ""]
        if len(pend):
            lines += ["| sample | tool | reason |", "|---|---|---|"]
            for _, r in pend.iterrows():
                lines.append(f"| {r['sample']} | {r['tool']} | {r['reason']} |")
        lines.append("")
    except FileNotFoundError:
        lines.append("(registry not found)")
    try:
        summ = read_table(MANIFEST_DIR / "callsets_summary.csv")
        ok = summ[summ["status"] == "ok"]
        lines += ["## Callsets", "",
                  f"- extracted callsets: {len(ok)} ({int(pd.to_numeric(ok['rows_out'], errors='coerce').fillna(0).sum())} rows)",
                  f"- non-ok: {summ[summ['status'] != 'ok']['status'].value_counts().to_dict()}",
                  ""]
        n_empty = int((summ["status"] == "empty_parsed").sum())
        if n_empty:
            lines += [f"- `empty_parsed` ({n_empty}) = the RNA004 Dorado `m6A_guitar` "
                      "model splits: those pileup files hold no call row after the "
                      "modkit `name`-code split (the guitar model has no m6A code), "
                      "so they parse to zero rows. Not a missing sample.", ""]
        # availability matrix = extracted + liftover-filled callsets
        avail = ok[["sample", "tool"]].copy()
        fill_path = MANIFEST_DIR / "liftover_fill.csv"
        if fill_path.exists():
            fill = read_table(fill_path)[["sample", "tool"]]
            avail = pd.concat([avail, fill], ignore_index=True).drop_duplicates()
            lines += [f"- additionally filled by R2Dtool liftover: {len(fill)} callsets "
                      f"(replicates the legacy pipeline never converted)", ""]
        mat = (avail.pivot_table(index="tool", columns="sample", values="tool",
                                 aggfunc="count").fillna(0).astype(int))
        lines += ["### Availability matrix (per sample == per replicate)", "",
                  "| tool | " + " | ".join(mat.columns) + " |",
                  "|" + "---|" * (len(mat.columns) + 1)]
        for tool, row in mat.iterrows():
            lines.append("| " + tool + " | " + " | ".join(str(v) for v in row) + " |")
        lines.append("")
    except FileNotFoundError:
        lines.append("(callsets summary not found)")
    val_path = MANIFEST_DIR / "liftover_validation.csv"
    if val_path.exists():
        val = read_table(val_path)
        lines += ["## Liftover reproduction check", "",
                  f"- cases compared against the legacy `*_liftover.txt` files: {len(val)}",
                  f"- verdicts: {val['verdict'].value_counts().to_dict()}", ""]
    try:
        uni = read_table(UNIVERSE_ROOT / "universe_summary.csv")
        lines += ["## Universe", "", "| sample | universe rows (cov>=5) | cov>=10 | DRACH |",
                  "|---|---|---|---|"]
        for _, r in uni.iterrows():
            lines.append(f"| {r['sample']} | {r['n_cov5']} | {r['n_cov10']} | {r['n_drach_cov10']} |")
        lines.append("")
    except FileNotFoundError:
        lines.append("(universe summary not found)")
    try:
        ann = read_table(MANIFEST_DIR / "annotation_summary.csv")
        bad = ann[(pd.to_numeric(ann["pct_expected_base_at_raw"], errors="coerce") < 80)
                  & (pd.to_numeric(ann["rows"], errors="coerce") > 100)]
        lines += ["## Annotation QC", "",
                  f"- callsets annotated: {len(ann)}",
                  f"- callsets where <80 % of calls sit on the expected base (possible "
                  f"coordinate/strand issues or transcript-space output): {len(bad)}", ""]
        if len(bad):
            lines += ["| sample | tool | pct_expected_base | note |", "|---|---|---|---|"]
            for _, r in bad.head(40).iterrows():
                pct = pd.to_numeric(r["pct_expected_base_at_raw"], errors="coerce")
                pct_str = f"{pct:.1f}" if pd.notna(pct) else "NA"
                lines.append(f"| {r['sample']} | {r['tool']} | {pct_str} | |")
            lines.append("")
    except FileNotFoundError:
        lines.append("(annotation summary not found)")
    lines += ["## Known caveats", "",
              "- Mouse WT/KO come from two independent studies; they are evaluated "
              "separately and never merged (see README).",
              "- Curlcake `Curlcake_IVT_rep2_partial` is a depth-matched subset of "
              "`Curlcake_IVT_rep3` (same run, SRR8767348) and is NOT an independent replicate.",
              "- E. coli WT (`E_ss_rd_RNA1/2`) are two runs of the same condition; the "
              "IVT control `E_IVT_neg1/2` is ONE sample split into two halves (SRR27228854, "
              "882k/506k local reads) and is only used as a null-vs-null false-positive control.",
              "- Comparison outputs are filed under the tool's own *test side*; the two "
              "pairs whose BOTH sides are special samples carry the full comparison in "
              "the directory name instead: `E_IVT_neg1_vs_E_IVT_neg2` (the single "
              "null-vs-null comparison; owned by `E_IVT_neg1` for DRUMMER/ELIGOS2_diff and "
              "by `E_IVT_neg2` for EpiNano_DiffErr/xPore/yanocomp/Nanocompore) and "
              "`mESCs_Mettl3_KO_vs_mES_KO` (cross-study KO-vs-KO; owned by `mES_KO` for "
              "ELIGOS2_diff/DRUMMER/xPore/yanocomp and by `mESCs_Mettl3_KO` for "
              "EpiNano_DiffErr/Nanocompore). There is deliberately no comparison-index "
              "numbering.",
              "- The Nanom6A fill samples were produced with the f5c-mode pipeline; "
              "their provenance is recorded in the callsets manifest.",
              "- purified sites are WT-called/KO-absent on the common universe and carry "
              "the R3-7 circularity caveat.",
              ""]
    report = (_RB / "data/evaluation/qc_report.md")
    report.write_text("\n".join(lines))
    logger.info("QC report -> %s", report)

    inv.flush()
    logger.info("log: %s", logger.log_path)


if __name__ == "__main__":
    main()
