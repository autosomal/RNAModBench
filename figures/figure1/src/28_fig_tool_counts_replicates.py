#!/usr/bin/env python
"""Figure 1B, replicate-aware (R3-2 / E6): sites detected by each tool, per replicate.

The published panel counted the sites of every tool per dataset group from the legacy
``output/<group>/Tools.txt`` aggregates.  Those aggregates came from a chain that
collapsed the replicate number out of the directory name (``docs/pipeline.md``),
so the published dots are *one* replicate in some groups (Arabidopsis used rep3 only)
and an undocumented pile-up of replicates in others (HeLa) -- the panel carried no
biological-replicate information at all.

This script recounts the very same quantity -- the number of distinct sites a tool
reported, deduplicated on ``(chrom, pos_raw)`` exactly as the old notebook deduplicated
``(Chr, Start, End, Tools)`` -- from ``harmonisation/callsets/``, where every sample is its
own file, and plots **one marker per condition**:

  dot          mean of the group's independent sequencing units
  grey line    connects the two groups' means
  title        carries the counts (n) of both conditions

The per-replicate counts stay in the tables, not on the panel: after trying a dot cloud
and an error bar, the author asked for the plain published layout with replicate-aware
values (2026-09-17).  Mouse is the exception -- two independent studies have no
admissible mean, so both study points are drawn (never pooled).

2026-09-19 (project rule): no figure may carry small annotation text any more.  The
grey n-line that used to sit above each panel is gone; its information moved into the
panel title (``_panel_title``).  "never pooled" and "hollow marker = 0 sites" are not
drawn at all any more -- the tables and the ``replicate_display`` column carry that
wording.

``config`` semantics are honoured, not re-derived:

* ``CROSS_STUDY_GROUPS`` (Mouse WT / KO) are two independent studies: both study points
  are drawn and never averaged into one marker;
* ``NON_INDEPENDENT_SAMPLE_OF`` (``E_IVT_neg1`` is one half of ``E_IVT_neg2``) means the
  E. coli IVT control contributes exactly one point;
* a site count of 0 (a tool ran on the sample and its own threshold kept nothing) is
  drawn as a solid marker clamped to the axis floor, because a log axis cannot show zero.

Outputs
-------
figures/figure1/inputs/figures/Fig1B_tool_counts_replicates.{pdf,png}
figures/figure1/inputs/tables/per_replicate_tool_counts.tsv
figures/figure1/inputs/tables/tool_counts_group_summary.tsv
figures/figure1/inputs/tables/legacy_vs_replicates_tool_counts.tsv

Usage
-----
python $RNAMODBENCH_ROOT/figures/figure1/src/28_fig_tool_counts_replicates.py
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
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

HERE = _RB / "src/harmonisation"
sys.path.insert(0, str(HERE))

from common import config as C                                    # noqa: E402
from common.figstyle import apply as apply_style                  # noqa: E402
from common.figstyle import save                                  # noqa: E402
from common.io_utils import write_table                           # noqa: E402
from common.manifest import setup_logger                          # noqa: E402

OUT = (_RB / "figures/figure1/inputs")
TAB, FIG = OUT / "tables", OUT / "figures"

PLATFORM = "RNA002"
MOD_TYPE = "m6A"
WT_COLOR = "#3778A0"      # kept from the published panel
TREAT_COLOR = "#F5B264"   # kept from the published panel
CONNECT = "0.35"

#: publication tool list and order (same set as 24_replicate_structure.py)
TOOL_ORDER = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
              "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
              "NanoSPA_m6A", "xPore", "yanocomp"]
#: ids whose published axis label is not the bare id
DISPLAY = {"yanocomp": "Yanocomp"}
#: the legacy ``Tools.txt`` directory labels vs the canonical tool ids
LEGACY_TOOL_ALIAS = {"ELIGOS_diff": "ELIGOS2_diff", "ELIGOS_solo": "ELIGOS2_solo",
                     "Yanocomp": "yanocomp"}


@dataclass(frozen=True)
class Panel:
    species: str
    wt_group: str
    treat_group: str
    treat_label: str

    @property
    def groups(self) -> tuple[str, str]:
        return (self.wt_group, self.treat_group)


PANELS = [
    Panel("Arabidopsis", "Arabidopsis_WT", "Arabidopsis_KD", "KD"),
    Panel("Mouse", "Mouse_WT", "Mouse_KO", "KO"),
    Panel("Human", "HeLa_WT", "HeLa_IVT", "IVT"),
    Panel("E.coli", "E.coli_WT", "E.coli_IVT", "IVT"),
]

#: the four Tools.txt files the published panel actually read (now in the archive)
LEGACY_SOURCES = {
    "Arabidopsis_WT": C.LEGACY_OUTPUT / "Arabidopsis_WT" / "Tools.txt",
    "Arabidopsis_KD": C.LEGACY_OUTPUT / "Arabidopsis_KD" / "Tools.txt",
    "Mouse_WT": C.LEGACY_OUTPUT / "Mouse_WT" / "Tools.txt",
    "Mouse_KO": C.LEGACY_OUTPUT / "Mouse_KO" / "Tools.txt",
    "HeLa_WT": C.LEGACY_OUTPUT / "HeLa_WT" / "aggregated_data" / "m6A" / "Tools.txt",
    "HeLa_IVT": C.LEGACY_OUTPUT / "HeLa_IVT" / "aggregated_data" / "m6A" / "Tools.txt",
    "E.coli_WT": C.LEGACY_OUTPUT / "E.coli_WT" / "Tools.txt",
    "E.coli_IVT": C.LEGACY_OUTPUT / "E.coli_IVT" / "Tools.txt",
}


# --------------------------------------------------------------------------- #
# counting
# --------------------------------------------------------------------------- #
def count_sites(path: Path) -> int | None:
    """Distinct sites in one callset file; ``None`` when the file does not exist.

    The key ``(chrom, pos_raw)`` mirrors what the published notebook did with
    ``(Chr, Start, End, Tools)``: one row is one site, and a repeated row must not
    inflate the count.  A header-only file is a real zero, not a missing file.
    """
    if not path.exists():
        return None
    if path.stat().st_size == 0:
        return 0
    df = pd.read_csv(path, sep="\t", usecols=["chrom", "pos_raw"],
                     dtype=str, keep_default_na=False)
    return int(df.drop_duplicates().shape[0])


def count_all_callsets(logger) -> pd.DataFrame:
    """One row per (dataset_group, tool, sample) with the site count of its callset."""
    audit = pd.read_csv(C.MANIFEST_DIR / "completeness_audit.csv", sep="\t",
                        dtype=str, keep_default_na=False)
    audit = audit[(audit.platform == PLATFORM) & (audit.mod_type == MOD_TYPE)]
    akey = {(r.sample, r.tool): r for _, r in audit.iterrows()}

    rows, n_missing, n_mismatch = [], 0, 0
    root = C.CALLSET_ROOT / PLATFORM
    for species_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        for group_dir in sorted(p for p in species_dir.iterdir() if p.is_dir()):
            mod_dir = group_dir / MOD_TYPE
            if not mod_dir.is_dir():
                continue
            for tool_dir in sorted(p for p in mod_dir.iterdir() if p.is_dir()):
                tool = tool_dir.name
                for f in sorted(tool_dir.glob("*.tsv")):
                    sample = f.stem
                    s = C.SAMPLES_BY_NAME.get(sample)
                    n = count_sites(f)
                    if n is None:
                        n_missing += 1
                        continue
                    a = akey.get((sample, tool))
                    n_audit = int(a.n_sites) if a is not None and a.n_sites != "" else -1
                    if n_audit >= 0 and n_audit != n:
                        n_mismatch += 1
                        logger.warning("count mismatch %s/%s/%s: callset=%d audit=%d",
                                       group_dir.name, tool, sample, n, n_audit)
                    rows.append({
                        "platform": PLATFORM,
                        "species": species_dir.name,
                        "dataset_group": group_dir.name,
                        "condition_class": s.condition_class if s else "",
                        "tool": tool,
                        "sample": sample,
                        "replicate_tag": s.replicate_tag if s else "",
                        "study": s.study if s else "",
                        "sequencing_unit": C.sequencing_unit(s) if s else sample,
                        "independence_class": C.independence_class(s) if s else "",
                        "n_sites": n,
                        "n_sites_audit": n_audit,
                        "callset": str(f.relative_to((_RB))),
                    })
    tab = pd.DataFrame(rows).sort_values(["species", "dataset_group", "tool", "sample"])
    logger.info("callset files read: %d (missing %d, audit mismatches %d)",
                len(tab), n_missing, n_mismatch)
    return tab


def legacy_counts(logger) -> pd.DataFrame:
    """Recompute the published numbers with the published notebook's own algorithm."""
    rows = []
    for group, path in LEGACY_SOURCES.items():
        if not path.exists():
            logger.warning("legacy Tools.txt missing: %s", path)
            continue
        df = pd.read_csv(path, sep="\t", header=0,
                         usecols=["Chr", "Start", "End", "Tools"],
                         dtype={"Chr": str}, keep_default_na=False)
        df = df.drop_duplicates(subset=["Chr", "Start", "End", "Tools"])
        for tool, n in df.groupby("Tools").size().items():
            rows.append({"dataset_group": group,
                         "tool": LEGACY_TOOL_ALIAS.get(tool, tool),
                         "legacy_label": tool,
                         "legacy_count": int(n), "legacy_file": str(path)})
    tab = pd.DataFrame(rows)
    outside = sorted(set(tab.tool) - set(TOOL_ORDER))
    if outside:
        logger.info("legacy files also list non-manuscript tools: %s", ", ".join(outside))
    renamed = sorted(set(tab[tab.legacy_label != tab.tool].legacy_label))
    if renamed:
        logger.info("legacy tool labels mapped: %s", ", ".join(renamed))
    return tab


# --------------------------------------------------------------------------- #
# per-group aggregation
# --------------------------------------------------------------------------- #
def unit_members(group: str) -> dict[str, list[str]]:
    """{sequencing-unit key: members}, the representative sample first."""
    members = [s for s in C.SAMPLES
               if s.dataset_group == group and s.platform == PLATFORM]
    names = {s.canonical for s in members}
    out = {}
    for unit, mem in C.independent_units(members).items():
        rep = unit if unit in names else mem[0]
        out[unit] = [rep] + [m for m in mem if m != rep]
    return out


def units_of(group: str) -> dict[str, str]:
    """{sequencing-unit key: representative sample} for one dataset group."""
    return {u: mem[0] for u, mem in unit_members(group).items()}


def group_points(per_rep: pd.DataFrame, group: str, tool: str) -> list[tuple[str, int]]:
    """[(unit label, site count)] -- exactly one entry per independent sequencing unit.

    A unit is only ever counted once.  When the representative half has no callset
    but another member of the same unit does (E. coli IVT: DRUMMER was run on
    ``E_IVT_neg1``, the other half of ``E_IVT_neg2``'s run), the unit's value is
    taken from that member and the label records it -- dropping the point would
    hide a real zero.
    """
    sub = per_rep[(per_rep["dataset_group"] == group) & (per_rep["tool"] == tool)]
    pts = []
    for unit, members in sorted(unit_members(group).items()):
        for name in members:
            hit = sub[sub["sample"] == name]   # .sample is a method: bracket access
            if len(hit):
                label = unit if name == members[0] else f"{unit}[via {name}]"
                pts.append((label, int(hit["n_sites"].iloc[0])))
                break
    return pts


def summarise(per_rep: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for panel in PANELS:
        tools = [t for t in TOOL_ORDER
                 if not per_rep[(per_rep["dataset_group"].isin(panel.groups))
                                & (per_rep["tool"] == t)].empty]
        for group, cond in zip(panel.groups, ("WT", panel.treat_label)):
            klass = C.replicate_class(group)
            for tool in tools:
                pts = group_points(per_rep, group, tool)
                vals = [v for _, v in pts]
                if not vals:
                    continue
                mean = float(np.mean(vals))
                sd = float(np.std(vals, ddof=1)) if len(vals) > 1 else np.nan
                if klass == "cross_study":
                    disp = ("two study points, never pooled"
                            if len(vals) > 1 else "one study only for this tool")
                elif len(vals) < 2:
                    disp = "single unit, one point"
                else:
                    disp = f"mean of {len(vals)} independent units"
                rows.append({
                    "species": panel.species, "dataset_group": group, "condition_class": cond,
                    "tool": tool, "replicate_class": klass, "n_units": len(vals),
                    "units": "|".join(u for u, _ in pts),
                    "sites_each": "|".join(str(v) for _, v in pts),
                    "sites_mean": mean, "sites_sd": sd,
                    "sites_min": min(vals), "sites_max": max(vals),
                    "replicate_display": disp,
                })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# drawing
# --------------------------------------------------------------------------- #
def _order_tools(per_rep: pd.DataFrame, panel: Panel) -> list[str]:
    """Tools with a callset on *both* sides, sorted by the WT side's mean count.

    A tool that ran on one condition only cannot be drawn as a dumbbell (the
    published panel would have shown the missing side as a zero); it is left out
    of the panel and named in the log instead.
    """
    have = []
    for t in TOOL_ORDER:
        if all(group_points(per_rep, g, t) for g in panel.groups):
            have.append(t)

    def key(tool: str) -> float:
        vals = [v for _, v in group_points(per_rep, panel.wt_group, tool)]
        return float(np.mean(vals)) if vals else 0.0

    return sorted(have, key=key, reverse=True)


def _spread(pts: list[tuple[str, int]], floor: float,
            cross_study: bool) -> tuple[float, float] | None:
    """The x extent the dumbbell line has to cover for one group.

    Replicate groups are joined at their mean (the large dot); cross-study groups
    at their two study points, because there is no admissible mean there.
    """
    xs = [max(float(v), floor) for _, v in pts]
    if not xs:
        return None
    return (min(xs), max(xs)) if cross_study else (float(np.mean(xs)),) * 2


def _draw_group(ax, y: float, pts: list[tuple[str, int]], color: str,
                floor: float, cross_study: bool, ms_scale: float = 1.0) -> None:
    """One marker for the group: the mean of its independent sequencing units.

    Cross-study groups (Mouse) get their two study points instead, because there is
    no mean that may be pooled from two different studies.  Individual replicate
    values are not drawn -- they live in the two tables.

    ``ms_scale`` shrinks every marker: the printed-size page rebuild (67) draws the
    same panel in a 116 x 153 pt cell, where the 7-8.5 pt markers of this 7.4 in
    draft would touch the neighbouring row.
    """
    if not pts:
        return
    vals = [v for _, v in pts]
    if cross_study:                       # two studies: two points, no summary
        xs = [max(float(v), floor) for v in vals]
        if len(xs) > 1:
            ax.plot([min(xs), max(xs)], [y, y], color=color, lw=1.5 * ms_scale,
                    alpha=.75, zorder=2)
        for v, x in zip(vals, xs):
            ax.plot(x, y, "o", ms=7.0 * ms_scale, mfc=color, mec="white",
                    mew=.9, zorder=4)
        return

    mean = float(np.mean(vals))
    ax.plot(max(mean, floor), y, "o", ms=8.5 * ms_scale, mfc=color, mec="white",
            mew=.9, zorder=5)


def _panel_title(panel: Panel) -> str:
    """``"<species>\\nWT n = ... / <treat> n = ..."`` -- the title is the species name only.

    Since 2026-09-19 no figure of this project may carry small annotation text,
    so the grey n-line that used to sit above the axes moved into the panel
    title.  The wording of the old note is kept, minus "never pooled" (the
    legend says "two studies" already) and minus the zero-marker hint.
    """
    # species name only (revision request 2026-09-19): neither the old grey note
    # nor bracketed counts belong in a panel title; the per-condition n lives
    # in the legend ("mean of 3" / "two studies") and in the tables
    return panel.species


#: font sizes of this 7.4 x 6.74 in draft; the printed-size page rebuild (67) passes
#: its own dictionary instead, so these stay the reference rendering of this script.
FS_DEFAULT = {"tick": 10.5, "title": 14.0, "axis": 13.0, "legend": 9.5}


def draw_panel(ax, panel: Panel, per_rep: pd.DataFrame, *, fs: dict | None = None,
               legend_loc: str = "upper left", ms_scale: float = 1.0,
               title_pad: float = 11.0, order: list[str] | None = None,
               show_ylabels: bool = True) -> None:
    """Draw one species dumbbell of the 2 x 2 layout into ``ax``.

    Split out of :func:`figure` on 2026-09-21 so the printed-size page rebuild
    (``67_fig1_panels.py``) draws *this* panel rather than a copy of it; every
    default reproduces this script's own output.

    ``order`` overrides the per-species sort and ``show_ylabels`` suppresses the
    tool names: on the printed page the four sub-panels share one row order, and
    the 40 pt label column is printed once for each row of the grid (the
    layout Figure 6 uses, where the names appear in column 0 only) -- a second
    identical column would run into the neighbouring panel.
    """
    f = {**FS_DEFAULT, **(fs or {})}
    tools = list(order) if order is not None else _order_tools(per_rep, panel)
    cross = C.replicate_class(panel.wt_group) == "cross_study"
    pts_of = {t: (group_points(per_rep, panel.wt_group, t),
                  group_points(per_rep, panel.treat_group, t)) for t in tools}
    all_vals = [v for wt, tt in pts_of.values() for _, v in wt + tt]
    pos = [v for v in all_vals if v > 0]
    floor = max(0.4, min(pos) * 0.5) if pos else 0.4
    top = max(pos) * 2.2 if pos else 10

    for i, tool in enumerate(tools):
        wt, tt = pts_of[tool]
        ends = [x for pts in (wt, tt) for x in (_spread(pts, floor, cross) or ())]
        if not ends:
            # a shared row order can carry a tool this dataset never saw
            # (DENA / MINES were not run on E. coli): the row stays empty
            continue
        ax.plot([min(ends), max(ends)], [i, i], color=CONNECT,
                lw=1.1 * ms_scale, alpha=.45, zorder=1, solid_capstyle="round")
        _draw_group(ax, i, wt, WT_COLOR, floor, cross, ms_scale)
        _draw_group(ax, i, tt, TREAT_COLOR, floor, cross, ms_scale)

    ax.set_yticks(range(len(tools)))
    ax.set_yticklabels([DISPLAY.get(t, t) for t in tools] if show_ylabels else [],
                       fontsize=f["tick"])
    ax.set_ylim(len(tools) - 0.5, -0.5)
    ax.set_xscale("log")
    ax.set_xlim(floor, top)
    ax.set_title(_panel_title(panel), fontweight="bold", fontsize=f["title"],
                 pad=title_pad)
    if panel in PANELS[2:]:
        ax.set_xlabel("Counts (Log Scale)", fontsize=f["axis"])
    if cross:
        handles = [
            Line2D([], [], marker="o", ls="", ms=6.2 * ms_scale, mfc=WT_COLOR,
                   mec="white", mew=.8),
            Line2D([], [], marker="o", ls="", ms=6.2 * ms_scale, mfc=TREAT_COLOR,
                   mec="white", mew=.8),
        ]
        labels = ["WT (two studies)", f"{panel.treat_label} (two studies)"]
    else:
        handles = [
            Line2D([], [], marker="o", ls="", ms=8 * ms_scale, mfc=WT_COLOR,
                   mec="white", mew=.9),
            Line2D([], [], marker="o", ls="", ms=8 * ms_scale, mfc=TREAT_COLOR,
                   mec="white", mew=.9),
        ]
        labels = [f"WT (mean of {len(units_of(panel.wt_group))})",
                  f"{panel.treat_label} (mean of {len(units_of(panel.treat_group))})"]
    ax.legend(handles=handles, labels=labels, loc=legend_loc, frameon=False,
              fontsize=f["legend"], handletextpad=.5, labelspacing=.35,
              borderpad=.2)


def figure(per_rep: pd.DataFrame) -> None:
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 6.74), constrained_layout=True)

    for ax, panel in zip(axes.ravel(), PANELS):
        draw_panel(ax, panel, per_rep)

    FIG.mkdir(parents=True, exist_ok=True)
    stem = FIG / "Fig1B_tool_counts_replicates"
    save(fig, str(stem))
    print("wrote", stem.with_suffix(".pdf").name, "+ .png")


# --------------------------------------------------------------------------- #
def main() -> None:
    logger = setup_logger("28_fig_tool_counts_replicates")
    TAB.mkdir(parents=True, exist_ok=True)

    per_rep = count_all_callsets(logger)
    write_table(per_rep, TAB / "per_replicate_tool_counts.tsv")

    summ = summarise(per_rep)
    write_table(summ, TAB / "tool_counts_group_summary.tsv")

    for panel in PANELS:
        one_side = [t for t in TOOL_ORDER
                    if any(group_points(per_rep, g, t) for g in panel.groups)
                    and not all(group_points(per_rep, g, t) for g in panel.groups)]
        if one_side:
            logger.info("%s: one condition only, not drawn as a dumbbell: %s",
                        panel.species, ", ".join(one_side))

    legacy = legacy_counts(logger)
    if len(legacy):
        cmp = summ.merge(legacy, on=["dataset_group", "tool"], how="left")
        unmatched = cmp[cmp["legacy_count"].isna()]
        if len(unmatched):
            logger.warning("no legacy row for: %s",
                           ", ".join(f"{r.dataset_group}/{r.tool}"
                                     for _, r in unmatched.iterrows()))
        cmp["side"] = cmp["condition_class"]
        cmp["delta_mean_minus_legacy"] = cmp["sites_mean"] - cmp["legacy_count"]
        cmp["ratio_mean_over_legacy"] = cmp["sites_mean"] / cmp["legacy_count"]
        cmp["legacy_source_units"] = np.where(
            cmp["n_units"] > 1, "pooled >1 unit", "one unit")
        cols = ["species", "dataset_group", "side", "tool", "n_units", "units",
                "sites_each", "sites_mean", "sites_sd",
                "legacy_count", "delta_mean_minus_legacy", "ratio_mean_over_legacy",
                "legacy_source_units", "legacy_file"]
        cmp = cmp.sort_values(["species", "dataset_group", "tool"])
        write_table(cmp[cols], TAB / "legacy_vs_replicates_tool_counts.tsv")

    figure(per_rep)

    logger.info("--- WT side: replicate mean vs the published single-number value ---")
    for _, r in summ[summ.condition_class == "WT"].iterrows():
        logger.info("%-16s %-14s n=%d mean=%.0f sd=%s",
                    r.dataset_group, r.tool, r.n_units, r.sites_mean,
                    "-" if pd.isna(r.sites_sd) else f"{r.sites_sd:.0f}")
    logger.info("output: %s", OUT)


if __name__ == "__main__":
    main()
