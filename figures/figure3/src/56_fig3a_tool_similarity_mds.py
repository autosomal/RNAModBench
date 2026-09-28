#!/usr/bin/env python
"""Figure 3A, revision rebuild (R3-2 / E6 replicate structure, R1-6 tool scope).

The published panel A (``manuscript_figure_code/Figure3/Hierarchical_Clustering_of_Tools.ipynb``)
read one pre-aggregated ``output/<group>/Tools.txt`` per species and computed the
tool x tool Jaccard distance on **1 Mbp genomic bins** (a bin counts as shared if
*both* tools call somewhere inside it), then drew an MDS of the three species.
That panel therefore (i) carried no replicate structure and (ii) measured
regional co-occupancy rather than site-level agreement.

This script recomputes the same panel from ``harmonisation/callsets`` (the
analysis-ready layer, one file per sample = per independent unit):

* site granularity: single-nucleotide ``(chrom, start)`` sets, no binning;
* primary similarity: Jaccard index between the two tools' site sets;
* replicate structure: one Jaccard matrix per independent unit, plus the
  majority consensus over units (``common.consensus.quorum``), exactly as the
  metagene panels do;
* MDS: deterministic classical MDS (Torgerson/PCoA, eigen-decomposition of the
  double-centred distance matrix) -- no random initialisation, so the panel is
  reproducible by construction.  The legacy sklearn SMACOF solution
  (``random_state=42``) is only a rotation/reflection of the same configuration;
  the axes of an MDS are arbitrary, the configuration is what is read.
  sklearn is not installed in ``benchmark-revision``; classical MDS is the
  deterministic substitute and its Kruskal stress-1 is reported per species.
* unit solutions are overlaid on the consensus solution after a Procrustes
  similarity fit (rotation + translation + uniform scale, ``scipy.spatial``),
  which is the only meaningful way to compare two MDS configurations; the
  residual disparity is tabulated per tool and unit (the drift the reviewer's
  "consistently incorporate replicates" request asks to see).

Mouse uses its two independent studies (SRP357195 / SRP166020) as the two units
and never averages them ("majority" there = both studies agree).  Tools whose
call set is empty in any unit of a group are excluded from that group's MDS and
listed with the reason in ``fig3a_tool_inclusion.tsv`` (R1-6).

Outputs -> ``figures/figure3/{tables,figures,logs}``

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/figures/figure3/src/56_fig3a_tool_similarity_mds.py
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
import logging
import sys
from itertools import cycle
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.spatial import procrustes

HERE = _RB / "src/harmonisation"
sys.path.insert(0, str(HERE))

from common import config as C                                        # noqa: E402
from common.consensus import quorum                                   # noqa: E402
from common.figstyle import apply as apply_style                      # noqa: E402
from common.io_utils import write_table                               # noqa: E402
from common.manifest import setup_logger                              # noqa: E402

OUT = (_RB / "figures/figure3")
TAB, FIG, LOG = OUT / "tables", OUT / "figures", OUT / "logs"
CLEAN = (_RB / "data/callsets")

#: species -> (dataset group, [independent units in order], title)
GROUPS: dict[str, tuple[str, list[str], str]] = {
    "Arabidopsis": ("Arabidopsis_WT",
                    ["Arabidopsis_WT_rep1", "Arabidopsis_WT_rep2", "Arabidopsis_WT_rep3"],
                    "Arabidopsis"),
    "Mouse": ("Mouse_WT", ["mES_WT", "mESCs_Mettl3_WT"], "Mouse"),
    "Human": ("HeLa_WT", ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"], "Human"),
}
#: unit-name convention in the figure/table: replicates vs the two mouse studies
UNIT_LABEL = {
    "Arabidopsis_WT_rep1": "rep1", "Arabidopsis_WT_rep2": "rep2",
    "Arabidopsis_WT_rep3": "rep3",
    "HeLa_WT1": "rep1", "HeLa_WT2": "rep2", "HeLa_WT3": "rep3",
    "mES_WT": "studyB (SRP357195)", "mESCs_Mettl3_WT": "studyA (SRP166020)",
}

#: published panel-A identity of every tool: the legacy notebook cycled tab20
#: colours and a marker list over ``sorted(tool)``; the same rule is replayed
#: here so the redrawn panel keeps the paper's colour/marker per tool.
_MARKERS = ["o", "s", "^", "D", "v", "<", ">", "p", "*", "h"]
TOOLS_ORDER = sorted(C.ARTICLE_M6A_TOOLS)
_COLORS = plt.cm.tab20(np.linspace(0, 1, 20))
_CC = cycle(_COLORS)
_MC = cycle(_MARKERS)
TOOL_STYLE = {t: {"color": next(_CC), "marker": next(_MC)} for t in TOOLS_ORDER}

#: panel canvas (final printed size: 0.95 x \textwidth of the Wiley USG layout).
#: The row height is the one Figure 4 panel C uses for the same three-species
#: layout (1.98 in), so the MDS maps are drawn at the size the reader already
#: meets there instead of the 1.42 in that squeezed them; the extra 0.56 in
#: keeps the assembled page inside the 8.60 in budget.
PANEL_W, PANEL_H = 6.66, 1.98
FS = {"title": 9.0, "label": 8.0, "tick": 7.6, "legend": 7.2}
MDS_COMPONENTS = 2


# --------------------------------------------------------------------------- #
def sites_path(species: str, group: str, tool: str, sample: str) -> Path:
    return CLEAN / C.PLATFORM_RNA002 / species / group / "m6A" / tool / f"{sample}.tsv"


def load_sites(path: Path, logger: logging.Logger) -> pd.DataFrame:
    """Unique ``(chrom, pos)`` sites of one tool in one sample (empty if none)."""
    if not path.exists():
        logger.warning("callset missing: %s", path)
        return pd.DataFrame(columns=["chrom", "pos"])
    try:
        d = pd.read_csv(path, sep="\t", usecols=["chrom", "start"],
                        dtype={"chrom": str, "start": "int64"})
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=["chrom", "pos"])
    if d.empty:
        return pd.DataFrame(columns=["chrom", "pos"])
    return (pd.DataFrame({"chrom": d["chrom"].astype(str),
                          "pos": d["start"].astype("int64")})
            .drop_duplicates().reset_index(drop=True))


def site_keys(df: pd.DataFrame) -> np.ndarray:
    """``(chrom, pos)`` -> int64 keys, shared across the tools of one unit."""
    if df.empty:
        return np.empty(0, dtype="int64")
    cidx = pd.factorize(df["chrom"], sort=True)[0].astype("int64")
    return cidx * np.int64(10 ** 9) + df["pos"].to_numpy("int64")


def jaccard_matrix(keys: dict[str, np.ndarray], tools: list[str]) -> pd.DataFrame:
    """Pairwise Jaccard index of the tools' site-key sets."""
    sets = {t: set(keys[t].tolist()) for t in tools}
    m = np.ones((len(tools), len(tools)), dtype=float)
    for i, a in enumerate(tools):
        for j in range(i + 1, len(tools)):
            b = tools[j]
            inter = len(sets[a] & sets[b])
            union = len(sets[a]) + len(sets[b]) - inter
            jac = inter / union if union else np.nan
            m[i, j] = m[j, i] = jac
    return pd.DataFrame(m, index=tools, columns=tools)


def classical_mds(d: np.ndarray) -> np.ndarray:
    """Torgerson classical MDS (PCoA) of a symmetric distance matrix.

    Deterministic (plain eigen-decomposition), no random initialisation.
    """
    n = d.shape[0]
    a = -0.5 * d ** 2
    j = np.eye(n) - np.ones((n, n)) / n
    b = j @ a @ j
    vals, vecs = np.linalg.eigh(b)
    order = np.argsort(vals)[::-1][:MDS_COMPONENTS]
    vals, vecs = vals[order], vecs[:, order]
    vals = np.clip(vals, 0, None)
    return vecs * np.sqrt(vals)


def kruskal_stress_1(d: np.ndarray, coords: np.ndarray) -> float:
    """Kruskal stress-1: sqrt(sum((d_ij - dhat_ij)^2) / sum(d_ij^2))."""
    dd = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(-1))
    iu = np.triu_indices(d.shape[0], 1)
    obs, fit = d[iu], dd[iu]
    denom = (obs ** 2).sum()
    return float(np.sqrt(((obs - fit) ** 2).sum() / denom)) if denom else np.nan


# --------------------------------------------------------------------------- #
def build(logger: logging.Logger) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    coords_rows: list[dict] = []
    incl_rows: list[dict] = []
    sim_rows: list[dict] = []
    jac_frames: dict[str, pd.DataFrame] = {}

    for species, (group, units, title) in GROUPS.items():
        per_unit: dict[str, dict[str, np.ndarray]] = {}
        for unit in units:
            per_unit[unit] = {t: site_keys(load_sites(sites_path(species, group, t, unit),
                                                       logger))
                              for t in TOOLS_ORDER}
        sizes = pd.DataFrame({u: {t: len(v) for t, v in per_unit[u].items()}
                              for u in units})
        # a tool must be non-empty in EVERY unit of the group to enter the MDS
        usable = [t for t in TOOLS_ORDER if (sizes.loc[t] > 0).all()]
        for t in TOOLS_ORDER:
            n_empty = int((sizes.loc[t] == 0).sum())
            incl_rows.append({
                "species": species, "dataset_group": group, "tool": t,
                "n_units": len(units),
                "n_units_with_calls": int((sizes.loc[t] > 0).sum()),
                "n_sites_per_unit": "|".join(str(int(v)) for v in sizes.loc[t]),
                "included_in_mds": t in usable,
                "reason": ("all units have calls" if t in usable
                           else f"empty call set in {n_empty} of {len(units)} units"),
            })
        if len(usable) < 3:
            logger.warning("%s: only %d tools usable -- MDS skipped", species, len(usable))
            continue

        # per-unit similarity matrices
        unit_jac = {u: jaccard_matrix(per_unit[u], usable) for u in units}
        for u, jm in unit_jac.items():
            out = TAB / f"fig3a_jaccard_{species}_{UNIT_LABEL[u].split()[0]}.tsv"
            write_table(jm.rename_axis("tool").reset_index(), out)

        # majority consensus site set per tool, then its similarity matrix
        cons_keys: dict[str, np.ndarray] = {}
        for t in usable:
            stacked = np.concatenate([per_unit[u][t] for u in units])
            uniq, cnt = np.unique(stacked, return_counts=True)
            cons_keys[t] = uniq[cnt >= quorum(len(units))]
        cons_jac = jaccard_matrix(cons_keys, usable)
        write_table(cons_jac.rename_axis("tool").reset_index(),
                    TAB / f"fig3a_jaccard_{species}_consensus.tsv")

        d_cons = 1.0 - cons_jac.to_numpy(dtype=float)
        x_cons = classical_mds(d_cons)
        stress = kruskal_stress_1(d_cons, x_cons)
        logger.info("%-12s consensus MDS: %d tools, Kruskal stress-1 = %.4f",
                    species, len(usable), stress)
        # per-tool similarity summary (which tool a tool is closest to) -- the
        # quantitative basis of the "motif-trained tools cluster together" claim
        for t in usable:
            others = cons_jac.loc[t].drop(t)
            sim_rows.append({
                "species": species, "tool": t, "n_sites_consensus": len(cons_keys[t]),
                "mean_jaccard_to_others": float(others.mean()),
                "max_jaccard_tool": str(others.idxmax()),
                "max_jaccard": float(others.max()),
                "kruskal_stress_1": stress,
            })
        for i, t in enumerate(usable):
            coords_rows.append({"species": species, "tool": t, "kind": "consensus",
                                "unit": "", "x": x_cons[i, 0], "y": x_cons[i, 1],
                                "procrustes_disparity": 0.0,
                                "kruskal_stress_1": stress,
                                "n_sites": len(cons_keys[t])})

        # per-unit solutions, Procrustes-fitted to the consensus configuration
        for u in units:
            d_u = 1.0 - unit_jac[u].to_numpy(dtype=float)
            x_u = classical_mds(d_u)
            _, x_al, disparity = procrustes(x_cons, x_u)
            for i, t in enumerate(usable):
                coords_rows.append({
                    "species": species, "tool": t, "kind": "unit", "unit": UNIT_LABEL[u],
                    "x": x_al[i, 0], "y": x_al[i, 1],
                    "procrustes_disparity": float(disparity),
                    "kruskal_stress_1": kruskal_stress_1(d_u, x_u),
                    "n_sites": int(sizes.loc[t, u])})

    coords = pd.DataFrame(coords_rows)
    incl = pd.DataFrame(incl_rows)
    sim = pd.DataFrame(sim_rows)
    return coords, incl, sim


# --------------------------------------------------------------------------- #
def figure(coords: pd.DataFrame, logger: logging.Logger) -> None:
    apply_style()
    fig = plt.figure(figsize=(PANEL_W, PANEL_H))
    # panels are drawn inside a fixed area and the tool legend gets its own
    # right-hand strip in figure coordinates, so the (wide) two-column legend
    # can never spill onto the human panel or off the canvas
    gs = fig.add_gridspec(1, 3, left=0.075, right=0.715, top=0.90, bottom=0.20,
                          wspace=0.30)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]

    for k, (species, (group, units, title)) in enumerate(GROUPS.items()):
        ax = axes[k]
        sub = coords[(coords.species == species)]
        cons = sub[sub.kind == "consensus"]
        if cons.empty:
            ax.axis("off")
            continue
        xs, ys = cons["x"].to_numpy(float), cons["y"].to_numpy(float)
        pad = 0.09 * max(np.ptp(xs), np.ptp(ys))
        for r in cons.itertuples():
            st = TOOL_STYLE[r.tool]
            #: 2026-09-24: the per-unit layer (open circles + grey spokes) is
            #: gone -- the panel now matches the Figure 4C MDS style (points
            #: only, black marker edge); unit coordinates stay in the frozen
            #: table fig3a_mds_coords.tsv
            ax.plot(r.x, r.y, linestyle="none", marker=st["marker"], markersize=4.6,
                    markerfacecolor=st["color"], markeredgecolor="black",
                    markeredgewidth=0.3, zorder=3)
        ax.set_title(title, fontsize=FS["title"], fontweight="bold", pad=3)
        ax.set_xlim(xs.min() - pad, xs.max() + pad)
        ax.set_ylim(ys.min() - pad, ys.max() + pad)
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_xlabel("MDS Dimension 1", fontsize=FS["label"], labelpad=2)
        ax.tick_params(labelsize=FS["tick"], length=3, width=0.8)
        if k == 0:
            ax.set_ylabel("MDS Dimension 2", fontsize=FS["label"], labelpad=2)
        # unit spread of the drawn configuration, as a plain subtitle-free statement
        logger.info("%s: %d consensus + %d unit points", species, len(cons),
                    int((sub.kind == "unit").sum()))

    # tool legend column (layout v3: back beside the panels, two columns so the
    # 13 entries fit the 1.55 in row; the figure-wide band keeps only the
    # condition / linetype keys)
    lax = fig.add_axes([0.720, 0.0, 0.270, 1.0])
    lax.axis("off")
    handles = [Line2D([], [], linestyle="none", marker=TOOL_STYLE[t]["marker"],
                      markersize=4.2, markerfacecolor=TOOL_STYLE[t]["color"],
                      markeredgecolor=TOOL_STYLE[t]["color"], label=t)
               for t in TOOLS_ORDER]
    lax.legend(handles=handles, loc="center left", bbox_to_anchor=(0.0, 0.5),
               frameon=False,
               fontsize=FS["legend"], handletextpad=0.24, labelspacing=0.34,
               columnspacing=0.6, ncol=2, borderaxespad=0.0)
    logger.info("panel A legend: %d tools in 2 columns", len(handles))
    # bold panel letter (house rule: letters may carry the panel identity)
    fig.text(0.004, 0.965, "A", fontsize=FS["title"], fontweight="bold",
             ha="left", va="top")

    FIG.mkdir(parents=True, exist_ok=True)
    stem = FIG / "Fig3A_tool_similarity_mds"
    # exact canvas (no tight bbox): the assembly step places these panels 1:1
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=300)
    plt.close(fig)
    logger.info("wrote %s.{pdf,png} (%.2f x %.2f in)", stem.name, PANEL_W, PANEL_H)


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figures-only", action="store_true")
    args = ap.parse_args()

    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    logger = setup_logger("56_fig3a_tool_similarity_mds", log_dir=LOG)

    coords, incl, sim = build(logger)
    write_table(coords, TAB / "fig3a_mds_coords.tsv")
    write_table(incl, TAB / "fig3a_tool_inclusion.tsv")
    write_table(sim, TAB / "fig3a_similarity_summary.tsv")
    # the figure-wide legend band (61) reads the tool coding from this table so
    # marker/colour cannot drift between panel A and the legend
    from matplotlib.colors import to_hex

    style = pd.DataFrame([{"tool": t, "marker": TOOL_STYLE[t]["marker"],
                           "color": to_hex(TOOL_STYLE[t]["color"])}
                          for t in TOOLS_ORDER])
    write_table(style, TAB / "fig3a_tool_style.tsv")
    logger.info("tables: fig3a_mds_coords.tsv (%d rows), fig3a_tool_inclusion.tsv, "
                "fig3a_similarity_summary.tsv, fig3a_jaccard_*.tsv", len(coords))
    if len(sim):
        top = sim.sort_values("max_jaccard", ascending=False).head(6)
        logger.info("closest tool pairs (consensus): %s",
                    "; ".join(f"{r.species}/{r.tool}<->{r.max_jaccard_tool}"
                              f"={r.max_jaccard:.3f}" for r in top.itertuples()))
    excl = incl[~incl.included_in_mds]
    if len(excl):
        logger.warning("excluded from the MDS: %s",
                       ", ".join(f"{r.species}/{r.tool}" for r in excl.itertuples()))
    if not args.figures_only:
        figure(coords, logger)
    logger.info("output: %s", OUT)


if __name__ == "__main__":
    main()
