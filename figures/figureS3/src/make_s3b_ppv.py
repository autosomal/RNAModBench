#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Panel B of the rebuilt Figure S3: PPV vs. GLORI (2 bp), per independent unit.

Source of record
----------------
``sites_v2/evaluation/tables/m6a_glori_confusion.tsv`` -- the revision's
per-unit confusion table (explicit universe = exonic, reference-base
compatible positions with coverage >= 10; nearest GLORI within the window;
``precision`` = TP/(TP+FP) = PPV).  This is the very table that
``export_figure_ready.py`` / Fig. 5A read, so the panel is by construction in
the same as the revised Figure 5 -- and therefore sits on the
``sites_clean`` layer, not on the retired legacy ``output/`` tree.

On top of that, ``cross_check()`` recomputes PPV at 2 bp straight from
``sites_clean/<...>/<tool>/<unit>.tsv`` + the per-unit ``__universe.tsv.gz``
+ the GLORI beds for a few unit/tool pairs and reports the delta, proving the
panel really is the sites_clean layer's numbers.

Units: Arabidopsis_WT_rep1-3, HeLa_WT1-3 and the TWO independent mouse WT mESC
samples ``mESCs_Mettl3_WT`` (SRP166020, drawn as **study A**) and ``mES_WT``
(SRP357195, **study B**) -- the two studies are drawn as separate series and are
NEVER averaged (user rule).  The in-figure legend prints only ``mouse study A``
/ ``mouse study B`` (same wording as Fig. S4/S5/S10, user 2026-09-21); the
sample ids and the accessions stay in ``tables/S3B_ppv_by_unit.tsv`` and in this
docstring.

Metric label: ``PPV vs. GLORI (2 bp)`` -- the published label "GLORI hit
rate" / "Hit Rate" is retired (fig5_revision/README.md).
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
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
#: the project's chromosome normaliser lives in the sites_v2 package
sys.path.insert(0, str(Path(str(_RB / "src/sites_v2"))))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import s3_common as sc


# --------------------------------------------------------------------------- #
# source-of-record table
# --------------------------------------------------------------------------- #
def load_confusion() -> pd.DataFrame:
    """Per-unit confusion rows at the main window for the m6A tool scope."""
    df = pd.read_csv(sc.CONFUSION_TSV, sep="\t", dtype={"sample": str})
    df = df[(df["window"] == sc.PPV_WINDOW) & df["tool"].isin(sc.M6A_TOOLS)]
    keep_units = {u for _, units in sc.UNITS_BY_SPECIES.values() for u in units}
    df = df[df["sample"].isin(keep_units)].copy()
    missing = keep_units - set(df["sample"])
    if missing:
        raise SystemExit(f"units missing from {sc.CONFUSION_TSV}: {sorted(missing)}")
    tools = sorted(df["tool"].unique())
    if tools != sorted(sc.M6A_TOOLS):
        raise SystemExit(f"tool scope mismatch: {tools}")
    df["ppv_w2"] = df["tp"] / (df["tp"] + df["fp"])
    print(f"[table] {sc.CONFUSION_TSV.name}: {len(df)} unit x tool rows, "
          f"w={sc.PPV_WINDOW}, {len(tools)} m6A tools, "
          f"{len(keep_units)} units")
    return df


# --------------------------------------------------------------------------- #
# sites_clean cross-check
# --------------------------------------------------------------------------- #
def _glori_positions(bed: Path) -> dict[str, np.ndarray]:
    """chrom (normalised) -> sorted start array of the GLORI reference."""
    df = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1],
                     names=["chrom", "start"], dtype={0: str})
    df["chrom"] = _norm(df["chrom"])
    return {c: np.sort(g["start"].to_numpy()) for c, g in df.groupby("chrom")}


def _norm(chrom) -> pd.Series:
    """Normalise chromosome labels on every side of a join.

    Uses the project's own ``common.match.fix_chromosome`` (lower-case canonical
    ``chr*``): a naive case-sensitive strip silently loses the human GLORI bed's
    lowercase ``chrx`` contig (3,119 sites -> 66 m6Anet calls in HeLa_WT1, i.e.
    0.9 % of that panel's PPV).
    """
    s = pd.Series(chrom) if not isinstance(chrom, pd.Series) else chrom
    try:
        from common.match import fix_chromosome
        return s.map(fix_chromosome)
    except Exception:  # pragma: no cover - fallback keeps the two sides aligned
        return s.astype(str).str.replace(r"^chr", "", regex=True).str.lower()


def _universe_positions(sample: str) -> dict[str, np.ndarray]:
    """pos arrays (0-based) of the per-unit candidate universe."""
    path = sc.UNIVERSE_ROOT / "RNA002"  # species folder resolved by caller
    raise NotImplementedError  # replaced below


def _universe_for(species: str, sample: str) -> dict[str, np.ndarray]:
    path = sc.UNIVERSE_ROOT / "RNA002" / species / f"{sample}__universe.tsv.gz"
    if not path.exists():
        raise FileNotFoundError(path)
    frames = []
    with gzip.open(path, "rt") as fh:
        for chunk in pd.read_csv(fh, sep="\t", chunksize=2_000_000,
                                 dtype={"chrom": str, "pos": np.int64,
                                        "coverage": np.int64}):
            frames.append(chunk[chunk["coverage"] >= 10][["chrom", "pos"]])
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if df.empty:
        return {}
    df["chrom"] = _norm(df["chrom"])
    return {c: np.sort(g["pos"].to_numpy()) for c, g in df.groupby("chrom")}


def _within(match_positions: np.ndarray, ref: np.ndarray, w: int) -> np.ndarray:
    """Boolean: each position in ``match_positions`` is within w bp of ``ref``."""
    if ref.size == 0 or match_positions.size == 0:
        return np.zeros(match_positions.size, dtype=bool)
    idx = np.searchsorted(ref, match_positions)
    left = np.clip(idx - 1, 0, ref.size - 1)
    right = np.clip(idx, 0, ref.size - 1)
    d = np.minimum(np.abs(match_positions - ref[left]),
                   np.abs(match_positions - ref[right]))
    return d <= w


def cross_check(table: pd.DataFrame, check_species: dict[str, list[str]]) -> pd.DataFrame:
    """Recompute PPV@2bp from sites_clean + universe + GLORI for sample pairs."""
    out = []
    for species, (group, units) in sc.UNITS_BY_SPECIES.items():
        picks = check_species.get(species)
        if not picks:
            continue
        refs = _glori_positions(sc.GLORI_BEDS[species])
        print(f"[cross-check] {species}: loading {len(units)} universe files ...")
        univ = {u: _universe_for(species, u) for u in units}
        for unit in units:
            for tool in picks:
                path = (sc.SITES_CLEAN / "RNA002" / species / group / "m6A"
                        / tool / f"{unit}.tsv")
                if not path.exists():
                    print(f"[cross-check] missing {path}")
                    continue
                df = pd.read_csv(path, sep="\t", dtype={"chrom": str})
                rows_total = len(df)
                # the evaluator scores SITE SETS: unique positions (chrom, start)
                df = (df.assign(chrom=_norm(df["chrom"]))
                        .drop_duplicates(subset=["chrom", "start"])
                        .reset_index(drop=True))
                n_calls_total = len(df)
                # in-universe membership per chromosome (universe already cov>=10)
                keep = np.zeros(len(df), dtype=bool)
                for chrom, sub in df.groupby(_norm(df["chrom"])):
                    u = univ[unit].get(chrom)
                    if u is None:
                        continue
                    idx = np.searchsorted(u, sub["start"].to_numpy())
                    left = np.clip(idx - 1, 0, u.size - 1)
                    right = np.clip(idx, 0, u.size - 1)
                    hit = np.minimum(np.abs(sub["start"].to_numpy() - u[left]),
                                     np.abs(sub["start"].to_numpy() - u[right])) == 0
                    keep[sub.index.to_numpy()] = hit
                calls = df[keep].reset_index(drop=True)
                tp = np.zeros(len(calls), dtype=bool)
                for chrom, sub in calls.groupby(_norm(calls["chrom"])):
                    ref = refs.get(chrom)
                    if ref is None:
                        continue
                    tp[sub.index.to_numpy()] = _within(sub["start"].to_numpy(), ref,
                                                       sc.PPV_WINDOW)
                n_in = len(calls)
                ppv = tp.sum() / n_in if n_in else float("nan")
                row = table[(table["sample"] == unit) & (table["tool"] == tool)].iloc[0]
                out.append({
                    "species": species, "unit": unit, "tool": tool,
                    "n_rows_sites_clean": rows_total,
                    "n_unique_positions_sites_clean": n_calls_total,
                    "n_in_universe_sites_clean": n_in,
                    "n_in_universe_table": int(row["n_calls_in_universe"]),
                    "ppv_sites_clean": round(float(ppv), 6),
                    "ppv_table": round(float(row["ppv_w2"]), 6),
                    "d_ppv": round(float(ppv) - float(row["ppv_w2"]), 6),
                })
                print(f"   {species}/{unit}/{tool}: sites_clean PPV={ppv:.6f} vs "
                      f"table {row['ppv_w2']:.6f} (in-universe {n_in} vs "
                      f"{int(row['n_calls_in_universe'])})")
    return pd.DataFrame(out)


# --------------------------------------------------------------------------- #
# drawing
# --------------------------------------------------------------------------- #
def _axes_rect(i: int, y_top: float, row_h: float, page: bool = False):
    x0 = sc.B_AXES_X0 + i * (sc.B_AXES_W + sc.B_AXES_GAP)
    h = sc.B_AXES_BOTTOM - sc.B_AXES_TOP
    if page:
        return (x0 / sc.PAGE_W, 1.0 - sc.B_AXES_BOTTOM / sc.PAGE_H,
                sc.B_AXES_W / sc.PAGE_W, h / sc.PAGE_H)
    y_b = sc.B_AXES_BOTTOM - y_top
    return (x0 / sc.PAGE_W, (row_h - y_b) / row_h,
            sc.B_AXES_W / sc.PAGE_W, h / row_h)


def _series(table: pd.DataFrame, species: str) -> tuple[list[str], dict]:
    """(tool order, per-unit pivots) for one species panel."""
    _, units = sc.UNITS_BY_SPECIES[species]
    sub = table[table["species"] == species]
    piv = sub.pivot_table(index="tool", columns="sample", values="ppv_w2")
    piv = piv[units]                                   # column order = units
    order = piv.mean(axis=1).sort_values(ascending=False).index.tolist()
    return order, dict(piv.loc[order])


def draw_row_b(fig: plt.Figure, table: pd.DataFrame, y_top: float, row_h: float,
               page: bool = False, titles: bool = True) -> None:
    """Three species panels of per-unit PPV vs. GLORI (2 bp)."""
    for i, species in enumerate(sc.SPECIES_ORDER):
        ax = fig.add_axes(_axes_rect(i, y_top, row_h, page=page))
        order, vals = _series(table, species)
        x = np.arange(len(order))
        color = sc.SPECIES_COLORS[species]

        if species == "Mouse":
            # two independent studies -- separate series, never averaged;
            # study A listed first so the legend reads A over B (as Fig. S4/S5)
            series = [("mESCs_Mettl3_WT", "s", "#C98A2E", 0.0),
                      ("mES_WT", "o", color, 1.0)]
        else:
            series = [(None, "o", color, 1.0)]

        if species == "Mouse":
            for unit, marker, col, face in series:
                y = vals[unit].to_numpy()
                ax.plot(x, y, color=col, marker=marker, markersize=3.4,
                        linewidth=sc.LW_MAIN, linestyle="-",
                        markerfacecolor=(col if face else "white"),
                        markeredgewidth=0.9,
                        #: legend carries the study numbers only -- sample ids and
                        #: accessions stay in the tables/provenance (user 2026-09-21)
                        label=sc.MOUSE_STUDY_LABEL[unit])
        else:
            unit_vals = [vals[u].to_numpy() for u in vals]
            ymean = np.mean(unit_vals, axis=0)
            ysd = np.std(unit_vals, axis=0, ddof=1)
            for k, unit in enumerate(vals):
                ax.plot(x, vals[unit].to_numpy(), linestyle="none", marker="o",
                        markersize=2.3, markerfacecolor="white",
                        markeredgecolor=color, markeredgewidth=0.6,
                        alpha=0.85, zorder=2)
            # mean ± SD across the independent units (user request 2026-09-20)
            ax.errorbar(x, ymean, yerr=ysd, color=color, capsize=2,
                        elinewidth=0.7, linestyle="none", alpha=0.9, zorder=2)
            ax.plot(x, ymean, color=color, marker="o", markersize=3.4,
                    linewidth=sc.LW_MAIN, markerfacecolor="white",
                    markeredgewidth=0.9, zorder=3,
                    label=f"mean ± SD (n={len(vals)} units)")

        ax.set_ylim(0, 1.02)
        ax.set_xlim(-0.6, len(order) - 0.4)
        ax.set_xticks(list(x))
        ax.set_xticklabels(order, rotation=45, ha="right", fontsize=sc.F_TICK)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        if i == 0:
            ax.set_ylabel(sc.YLABEL_PPV, fontsize=sc.F_AXIS_LABEL,
                          fontweight="bold", color=sc.INK)
            ax.set_yticklabels([f"{t:.1f}" for t in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)],
                               fontsize=sc.F_TICK)
        else:
            ax.set_yticklabels([])
        if titles:
            ax.set_title(species, fontsize=sc.F_TITLE, fontweight="bold",
                         color=sc.INK, pad=6)
        leg = ax.legend(frameon=False, loc="upper right",
                        fontsize=sc.F_LEGEND_B, handlelength=1.6,
                        borderpad=0.2, labelspacing=0.25, handletextpad=0.4)
        for txt in leg.get_texts():
            txt.set_color(sc.INK)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_color("black")
            spine.set_linewidth(sc.LW_SPINE)


# --------------------------------------------------------------------------- #
def write_tables(table: pd.DataFrame, check: pd.DataFrame) -> list[Path]:
    sc.TABLE_DIR.mkdir(parents=True, exist_ok=True)
    written = []
    by_unit = sc.TABLE_DIR / "S3B_ppv_by_unit.tsv"
    rows = table.copy()
    rows["study"] = rows["sample"].map(sc.MOUSE_STUDY).fillna("")
    cols = ["species", "dataset_group", "sample", "study", "tool", "window",
            "n_calls_total", "n_calls_in_universe", "n_calls_out_of_universe",
            "tp", "fp", "ppv_w2", "precision_ci_lo", "precision_ci_hi"]
    rows = rows[[c for c in cols if c in rows.columns]].rename(
        columns={"sample": "unit"})
    rows.to_csv(by_unit, sep="\t", index=False)
    written.append(by_unit)

    summary = []
    for species in sc.SPECIES_ORDER:
        _, units = sc.UNITS_BY_SPECIES[species]
        sub = table[table["species"] == species]
        for tool in sc.M6A_TOOLS:
            vals = sub[sub["tool"] == tool].set_index("sample")["ppv_w2"]
            vals = vals.reindex(units).dropna()
            summary.append({
                "species": species, "tool": tool, "n_units": len(vals),
                "units": ",".join(vals.index), "mean_ppv": round(vals.mean(), 6),
                "sd_ppv": round(vals.std(ddof=1), 6) if len(vals) > 1 else "",
            })
    sm = pd.DataFrame(summary).sort_values(["species", "mean_ppv"],
                                           ascending=[True, False])
    path = sc.TABLE_DIR / "S3B_ppv_mean_sd.tsv"
    sm.to_csv(path, sep="\t", index=False)
    written.append(path)

    if not check.empty:
        path = sc.TABLE_DIR / "S3B_ppv_sites_clean_crosscheck.tsv"
        check.to_csv(path, sep="\t", index=False)
        written.append(path)
    return written


def main() -> int:
    ap = argparse.ArgumentParser(description="S3 panel B: PPV vs. GLORI (2 bp)")
    ap.add_argument("--no-titles", action="store_true",
                    help="omit the species titles (house rule 2026-09-19)")
    ap.add_argument("--skip-crosscheck", action="store_true")
    args = ap.parse_args()

    sc.apply_page_style()
    table = load_confusion()

    check = pd.DataFrame()
    if not args.skip_crosscheck:
        check = cross_check(table, check_species={
            "Arabidopsis": ["m6Anet", "Nanom6A"],
            "Mouse": ["m6Anet", "Nanom6A"],
            "Human": ["m6Anet", "Nanom6A"],
        })
        if not check.empty:
            worst = check["d_ppv"].abs().max()
            print(f"[cross-check] max |delta PPV| = {worst:.2e} "
                  f"(sites_clean recomputation vs the evaluation table)")

    sc.PANEL_DIR.mkdir(parents=True, exist_ok=True)
    y_top, y_bottom = 200.0, 445.0
    row_h = y_bottom - y_top
    fig = plt.figure(figsize=(sc.PAGE_W / 72.0, row_h / 72.0))
    draw_row_b(fig, table, y_top=y_top, row_h=row_h, titles=not args.no_titles)
    for ext in ("pdf", "png"):
        fig.savefig(sc.PANEL_DIR / f"S3B_ppv.{ext}", facecolor="white")
    plt.close(fig)
    print(f"[write] {sc.PANEL_DIR / 'S3B_ppv.pdf'} (+png)")

    for p in write_tables(table, check):
        print(f"[write] {p}")

    # species-level mean across tools (context for the response letter)
    for species in sc.SPECIES_ORDER:
        sm = pd.read_csv(sc.TABLE_DIR / "S3B_ppv_mean_sd.tsv", sep="\t")
        m = sm[sm["species"] == species]["mean_ppv"].mean()
        print(f"[summary] {species}: mean PPV over {len(sc.M6A_TOOLS)} tools = {m:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
