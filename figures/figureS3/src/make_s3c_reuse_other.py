#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Panel C of the rebuilt Figure S3: per-replicate modification-ratio agreement.

The ANALYSIS behind this panel belongs to the sibling session and is NOT
re-run here:
  * ``sites_v2/scripts/15_mod_ratio_replicate_agreement.py`` -> cached tables in
    ``04_revision_analysis/mod_ratio_replicates/tables`` (read-only for us);
  * ``sites_v2/scripts/16_mod_ratio_regression_fig.py`` -> their standalone
    figure ``mod_ratio_replicates/figures/mod_ratio_regression_S3C_style.pdf``.

This module re-renders *their recipe* (per-replicate thin fit lines + bold
pooled fit, identity line, per-tool legend with r / Lin's CCC) onto this page's
panel-C box.  Rationale: their figure is laid out for a 21.6-in canvas, so
dropping it into a 2-in page panel would render its legend at ~5 pt and break
the house rule (>= 7 pt ticks/labels).  Tool colours and the tool list are
imported from their script so the visual identity is preserved; the sibling
directory is never written.

Outputs: panels/S3C_modratio_perrep.{pdf,png}, tables/S3C_provenance.tsv
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
import csv
import hashlib
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import s3_common as sc

OTHER_SCRIPT = (Path(str(_RB / "src/sites_v2/scripts"))
                / "16_mod_ratio_regression_fig.py")
OTHER_DIR = Path(str(_RB / "analysis/mod_ratio_replicates"))
MATCHED = (_RB / "analysis/mod_ratio_replicates/tables/mod_ratio_matched_sites.tsv")
SUMMARY = (_RB / "analysis/mod_ratio_replicates/tables/mod_ratio_summary_by_group.tsv")

PANELS = [("Arabidopsis", "Arabidopsis_WT"), ("Mouse", "Mouse_WT"),
          ("Human", "HeLa_WT")]
REP_ALPHA = 0.55


def load_other_recipe():
    """Import the sibling script to reuse its tool list / colours verbatim."""
    spec = importlib.util.spec_from_file_location("other_s3c", OTHER_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _axes_rect(i: int, y_top: float, row_h: float, page: bool = False):
    """Axes rect as figure fractions (``page=False``: relative to a row figure)."""
    x0 = sc.B_AXES_X0 + i * (sc.B_AXES_W + sc.B_AXES_GAP)
    h = sc.C_AXES_BOTTOM - sc.C_AXES_TOP
    if page:
        return (x0 / sc.PAGE_W, 1.0 - sc.C_AXES_BOTTOM / sc.PAGE_H,
                sc.B_AXES_W / sc.PAGE_W, h / sc.PAGE_H)
    y_b = sc.C_AXES_BOTTOM - y_top
    return (x0 / sc.PAGE_W, (row_h - y_b) / row_h,
            sc.B_AXES_W / sc.PAGE_W, h / row_h)


def _fit(xs: np.ndarray, ys: np.ndarray):
    if len(xs) < 2 or np.allclose(np.std(xs), 0):
        return None
    b, a = np.polyfit(xs, ys, 1)
    return a, b


def draw_row_c(fig: plt.Figure, sites: pd.DataFrame, summary: pd.DataFrame,
               tools: list[str], colors: dict[str, str], titles: bool = True,
               y_top: float = 420.0, row_h: float = 213.598,
               page: bool = False) -> None:
    """Three species panels at page geometry (sibling session's recipe)."""
    xx = np.linspace(0, 1, 2)
    for i, (species, group) in enumerate(PANELS):
        ax = fig.add_axes(_axes_rect(i, y_top, row_h, page=page))
        # agreement plot: equal aspect keeps the y=x identity line at 45 degrees
        ax.set_aspect("equal", adjustable="box")
        ax.plot([0, 1], [0, 1], ls="--", lw=sc.LW_AVG, color="black", zorder=1)
        handles = []
        for tool in tools:
            sub = sites[(sites.species == species) & (sites.group == group)
                        & (sites.tool == tool)]
            if sub.empty:
                continue
            col = colors[tool]
            for rep in sorted(sub["replicate_tag"].dropna().unique()):
                rs = sub[sub.replicate_tag == rep]
                f = _fit(rs["glori_ratio"].values, rs["tool_ratio"].values)
                if f is not None:
                    a, b = f
                    ax.plot(xx, np.clip(a + b * xx, 0, 1), color=col, lw=0.8,
                            alpha=REP_ALPHA, zorder=2)
            f = _fit(sub["glori_ratio"].values, sub["tool_ratio"].values)
            if f is not None:
                a, b = f
                ax.plot(xx, np.clip(a + b * xx, 0, 1), color=col, lw=sc.LW_REG,
                        zorder=3)
            row = summary[(summary.species == species) & (summary.group == group)
                          & (summary.tool == tool)]
            if len(row):
                r = row["pearson_r_mean"].values[0]
                ccc = row["ccc_mean"].values[0]
                lab = f"{tool}  r={r:.2f}, CCC={ccc:.2f}"
            else:
                lab = tool
            handles.append(plt.Line2D([], [], color=col, lw=sc.LW_REG, label=lab))
            

        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_xticklabels([f"{t:.1f}" for t in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)],
                           fontsize=sc.F_TICK)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        if i == 0:
            ax.set_yticklabels([f"{t:.1f}" for t in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)],
                               fontsize=sc.F_TICK)
            ax.set_ylabel("Tool Predicted Modification Ratio",
                          fontsize=sc.F_AXIS_LABEL, fontweight="bold", color=sc.INK)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel("GLORI Modification Ratio", fontsize=sc.F_AXIS_LABEL,
                      fontweight="bold", color=sc.INK)
        if titles:
            ax.set_title(species, fontsize=sc.F_TITLE, fontweight="bold",
                         color=sc.INK, pad=6)
        leg = ax.legend(handles=handles, loc="upper left", frameon=True,
                        fontsize=sc.F_LEGEND_C, handlelength=1.4,
                        borderpad=0.3, labelspacing=0.25)
        leg.get_frame().set_edgecolor("black")
        leg.get_frame().set_linewidth(0.6)
        for txt in leg.get_texts():
            txt.set_color(sc.INK)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_color("black")
            spine.set_linewidth(sc.LW_SPINE)


def write_provenance(summary: pd.DataFrame, tools: list[str]) -> Path:
    sc.TABLE_DIR.mkdir(parents=True, exist_ok=True)
    path = sc.TABLE_DIR / "S3C_provenance.tsv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["# panel C analysis (per-replicate fits + CCC) is the sibling "
                    "session's; this page re-renders its recipe at page geometry"])
        w.writerow(["artifact", "path", "md5", "mtime"])
        for p in (OTHER_SCRIPT, MATCHED, SUMMARY,
                  (_RB / "analysis/mod_ratio_replicates/figures/mod_ratio_regression_S3C_style.pdf")):
            if p.exists():
                md5 = hashlib.md5(p.read_bytes()).hexdigest()[:12]
                st = p.stat()
                w.writerow([p.name, str(p), md5,
                            pd.Timestamp(st.st_mtime, unit="s").isoformat()])
        w.writerow([])
        w.writerow(["species", "group", "tool", "n_replicates", "replicate_labels",
                    "pearson_r_mean", "pearson_r_sd", "ccc_mean", "ccc_sd",
                    "n_overlap_mean"])
        for sp, grp in PANELS:
            for tool in tools:
                row = summary[(summary.species == sp) & (summary.group == grp)
                              & (summary.tool == tool)]
                if row.empty:
                    continue
                r = row.iloc[0]
                w.writerow([sp, grp, tool, r.get("n_replicates", ""),
                            r.get("replicate_labels", ""),
                            f"{r['pearson_r_mean']:.4f}",
                            "" if pd.isna(r.get("pearson_r_sd")) else f"{r['pearson_r_sd']:.4f}",
                            f"{r['ccc_mean']:.4f}",
                            "" if pd.isna(r.get("ccc_sd")) else f"{r['ccc_sd']:.4f}",
                            r.get("n_overlap_mean", "")])
    return path


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description="S3 panel C (sibling session's recipe)")
    ap.add_argument("--no-titles", action="store_true")
    args = ap.parse_args()

    #: their module sets its own rcParams at import -> apply OUR page style after
    other = load_other_recipe()
    sc.apply_page_style()
    tools, colors = list(other.TOOLS), dict(other.TOOL_COLOR)
    sites = pd.read_csv(MATCHED, sep="\t")
    summary = pd.read_csv(SUMMARY, sep="\t")
    print(f"[input] {MATCHED.name}: {len(sites)} matched sites; "
          f"{SUMMARY.name}: {len(summary)} rows; tools={tools}")

    y_top, y_bottom = 420.0, 633.598
    fig = plt.figure(figsize=(sc.PAGE_W / 72.0, (y_bottom - y_top) / 72.0))
    draw_row_c(fig, sites, summary, tools, colors, titles=not args.no_titles,
               y_top=y_top, row_h=y_bottom - y_top)
    sc.PANEL_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(sc.PANEL_DIR / f"S3C_modratio_perrep.{ext}", facecolor="white")
    plt.close(fig)
    print(f"[write] {sc.PANEL_DIR / 'S3C_modratio_perrep.pdf'} (+png)")

    prov = write_provenance(summary, tools)
    print(f"[write] {prov}")
    for sp, grp in PANELS:
        sub = summary[(summary.species == sp) & (summary.group == grp)
                      & (summary.tool.isin(tools))]
        print(f"[summary] {sp}: " + ", ".join(
            f"{r.tool} r={r.pearson_r_mean:.2f}/CCC={r.ccc_mean:.2f}"
            for r in sub.itertuples()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
