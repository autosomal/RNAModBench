#!/usr/bin/env python3
"""63 -- acceptance gate for the revised Figure S5 (A-E, three columns).

The figure is *complementary by construction*: everything the rebuilt main
Figure 6 shows is out of scope here.  The gate therefore checks, in order:

1. the frozen evidence tables still carry the numbers the S5 legend quotes --
   the criterion facts (2,379 enumerated combinations per group, five selected
   optima, PPV >= p0 on every one of them, greedy path == exhaustive optima),
   the panel A anchors (intersection recall at k = 5 and its collapse to zero,
   the recall gain of the eight extra tools) and the panel B-E evidence
   (per-tool coverage / DRACH / control burden, site-set DRACH);
2. the figure rebuilds with the layout contract of ``62_figS5_figure.py`` --
   a compact **240 x 175 mm** page, five panel rows (A-E): row A in three
   species columns (Arabidopsis | Mouse with both studies in one axes | HeLa),
   rows B-E as single full-width panels (row D has no group dimension, the
   controls are not species-specific); the 13 tool names are printed **once and
   slanted 45 deg** under row D (right-aligned, no axis title) and the 5 set
   names once, horizontally, under row E -- the tool names are the *only*
   rotated text on the page, one panel letter per row, Arial only, no
   gridlines, no in-panel annotation text, no text collisions, no text off the
   canvas, legend inside the page, minimum font >= 9 pt;
3. no content of the main figure leaks back in: no dense combination cloud
   (the search-space plane lives in Fig. 6D) and no criterion wording in the
   figure itself;
4. the written PDF is exactly one A4 page with Arial embedded, the PNG is the
   same canvas at 300 dpi, and the paste-ready legend file documents panels
   A-E consistently with the byte-identical mirror in the Fig. 6 directory.

Usage: conda run -n benchmark-revision --no-capture-output python \
    src/harmonisation/scripts/63_verify_figS5.py
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
import importlib.util
import logging
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PROJECT = _RB
OUT = (_RB / "figures/figureS6")
FROZEN = (_RB / "figures/figure6")
TAB, FIG, LOG = (_RB / "figures/figure6/tables"), (_RB / "figures/figureS6/figures"), (_RB / "figures/figureS6/logs")
PDF = FIG / "FigureS5_rev.pdf"
PNG = FIG / "FigureS5_rev.png"
LEGEND_MD = FIG / "FigS5_legends.md"                       # canonical
LEGEND_MIRROR = (_RB / "figures/figure6/figures/FigS5_legends.md")    # kept identical

GROUPS = [("Arabidopsis", "Arabidopsis_WT"), ("Mouse", "studyA"),
          ("Mouse", "studyB"), ("Human", "HeLa_WT")]
P0 = {"Arabidopsis_WT": 1.50, "studyA": 0.70, "studyB": 0.70, "HeLa_WT": 0.85}
#: intersection recall at k = 5 (panel A anchor), per group
ISECT_K5 = {"Arabidopsis_WT": 0.12, "studyA": 2.24, "studyB": 1.74,
            "HeLa_WT": 1.58}
#: first k beyond 5 where the intersection is empty (isect PPV undefined)
ISECT_EMPTY_FROM = {"Arabidopsis_WT": 8, "studyA": 13, "studyB": None,
                    "HeLa_WT": 12}
#: mean DRACH fraction of the five site sets (panel E anchors)
SET_DRACH = {
    ("Arabidopsis", "Arabidopsis_WT"): (9.55, 67.25, 100.0, 100.0, 38.75),
    ("Mouse", "studyA"): (99.90, 37.23, 100.0, 100.0, 90.61),
    ("Mouse", "studyB"): (99.91, 33.35, 100.0, 100.0, 90.61),
    ("Human", "HeLa_WT"): (99.79, 14.51, 100.0, 100.0, 89.22),
}
#: configurations with no Curlcake run
NO_CURLCAKE = {"CHEUI_m6A", "DRUMMER", "ELIGOS2_diff"}

FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"[{'ok' if ok else 'FAIL'}] {name}" + (f" --- {detail}" if detail else ""),
          flush=True)
    if not ok:
        FAILURES.append(name)


def close(a: float, b: float, rtol: float = 1e-5) -> bool:
    return abs(float(a) - float(b)) <= rtol * max(abs(float(b)), 1e-12)


def rgb(colour: str) -> tuple[float, float, float]:
    """Hex colour -> the RGB tuple a scatter collection reports (4 decimals)."""
    h = colour.lstrip("#")
    return tuple(round(int(h[i:i + 2], 16) / 255.0, 4) for i in (0, 2, 4))


def printed_texts(fig) -> list[str]:
    """Every string the figure actually prints.

    ``fig.findobj(Text)`` double-counts each tick label (a Tick carries a
    ``label1`` and a ``label2`` artist that share the formatter), so the labels
    are taken from ``label1`` -- the one that is drawn -- plus the axis titles,
    the panel letters and the legend.
    """
    out = [t.get_text() for t in fig.texts]
    for ax in fig.axes:
        out += [t.label1.get_text() for t in ax.xaxis.get_major_ticks()]
        out += [t.label1.get_text() for t in ax.yaxis.get_major_ticks()]
        out += [t.get_text() for t in ax.texts]      # staggered tool names
        out += [ax.get_title(), ax.get_xlabel(), ax.get_ylabel()]
    for lg in fig.legends:
        out += [t.get_text() for t in lg.get_texts()]
    return [s for s in out if s.strip()]


def drawn_xticklabels(ax) -> list[str]:
    """The x tick labels that are printed under one panel (label1 only)."""
    return [t.label1.get_text() for t in ax.xaxis.get_major_ticks()]


def collection_colours(ax) -> set[tuple[float, float, float]]:
    """Every face and edge colour drawn inside one panel."""
    seen: set[tuple[float, float, float]] = set()
    for c in ax.collections:
        for arr in (c.get_facecolor(), c.get_edgecolor()):
            for rgba in np.atleast_2d(arr):
                seen.add(tuple(float(v) for v in np.round(rgba[:3], 4)))
    return seen


def main() -> None:
    space = pd.read_csv(TAB / "figS5_search_space.tsv", sep="\t")
    sel = pd.read_csv(TAB / "fig6_combination_selected.tsv", sep="\t")
    greedy = pd.read_csv(TAB / "figS5_greedy_1to13.tsv", sep="\t")
    toolq = pd.read_csv(TAB / "figS5_tool_quality.tsv", sep="\t")
    ctltool = pd.read_csv(TAB / "fig6_negative_control_fp_bytool.tsv", sep="\t")
    qual = pd.read_csv(TAB / "figS5_site_quality.tsv", sep="\t")

    # ---- 1a. criterion facts the legend states ------------------------------ #
    check("search space enumerates C(13,1..5)=2379 combinations per group",
          all(len(space[(space.species == sp) & (space.group == g)]) == 2379
              for sp, g in GROUPS), f"rows={len(space)}")
    check("exactly five selected optima per group",
          all(int(space[(space.species == sp) & (space.group == g)].selected.sum())
              == 5 for sp, g in GROUPS))
    check("feasible flag == (mean union PPV >= p0) on every row",
          bool((space.feasible
                == (space.union_precision_mean >= space.p0_chance_precision)).all()))
    check("p0 is a formal guardrail: every enumerated combination is feasible",
          bool(space.feasible.all()),
          "the S5 legend may only state this while the table says so")
    for sp, g in GROUPS:
        row = space[(space.species == sp) & (space.group == g)]
        check(f"chance level p0 {sp}/{g} = {P0[g]:.2f}%",
              close(100 * row.p0_chance_precision.iloc[0], P0[g], rtol=5e-3),
              f"{100 * row.p0_chance_precision.iloc[0]:.3f}%")
    sel_set = {(r.species, r.group, r.k, r.combination) for r in sel.itertuples()}
    sp_set = {(r.species, r.group, r.k, r.combination)
              for r in space[space.selected].itertuples()}
    check("selected rows identical to fig6_combination_selected.tsv",
          sel_set == sp_set, f"n={len(sel_set)}")
    check("greedy forward selection reproduces all exhaustive optima",
          bool(sel.greedy_agrees.all()))
    # the same statement, recomputed from the two frozen tables (the legend
    # makes it; the tables order the members differently, so compare the sets)
    bad: list[str] = []
    for sp, g in GROUPS:
        gp = greedy[(greedy.species == sp) & (greedy.group == g)].sort_values("k")
        check(f"greedy path has all 13 steps {sp}/{g}",
              sorted(gp.k) == list(range(1, 14)), f"n={len(gp)}")
        selg = space[(space.species == sp) & (space.group == g) & space.selected]
        for k in range(1, 6):
            want = frozenset(selg[selg.k == k].combination.iloc[0].split("+"))
            got = frozenset(gp[gp.k == k].combination.iloc[0].split("+"))
            if want != got:
                bad.append(f"{g}/k{k}")
    check("greedy k=1..5 == the enumerated optima (member sets)", not bad,
          f"{bad}")

    # ---- 1b. panel A anchors (the greedy path beyond the main figure) ------- #
    for sp, g in GROUPS:
        d = greedy[(greedy.species == sp) & (greedy.group == g)].sort_values("k")
        k5 = d[d.k == 5].iloc[0]
        check(f"intersection recall at k=5 {sp}/{g} = {ISECT_K5[g]}%",
              abs(100 * k5.isect_recall_mean - ISECT_K5[g]) <= 0.05,
              f"{100 * k5.isect_recall_mean:.2f}%")
        check(f"intersection recall reaches 0 at k=13 {sp}/{g}",
              abs(100 * d[d.k == 13].isect_recall_mean.iloc[0]) <= 1e-6,
              f"{100 * d[d.k == 13].isect_recall_mean.iloc[0]:.4f}%")
        nan_k = [int(k) for k, v in zip(d.k, d.isect_precision_mean)
                 if not np.isfinite(v)]
        exp = ISECT_EMPTY_FROM[g]
        check(f"intersection PPV undefined only where the intersection is empty "
              f"({sp}/{g})",
              (nan_k[0] if nan_k else None) == exp,
              f"first undefined k={nan_k[0] if nan_k else None} (expected {exp})")
        late = d[d.k >= 6]
        gain = 100 * (late.union_recall_mean.max() - k5.union_recall_mean)
        check(f"eight extra tools buy <= 8 pp recall {sp}/{g}",
              gain <= 8.0, f"+{gain:.1f} pp (46.7->52.4 style claim: {gain:.1f})")
        check(f"union PPV at k=13 is below k=5 {sp}/{g}",
              float(d[d.k == 13].union_precision_mean.iloc[0])
              < float(k5.union_precision_mean),
              f"{100 * float(d[d.k == 13].union_precision_mean.iloc[0]):.1f}%"
              f" < {100 * float(k5.union_precision_mean):.1f}%")

    # ---- 1c. panels B / C: per-tool site quality --------------------------- #
    for sp, g in GROUPS:
        q = toolq[(toolq.species == sp) & (toolq.group == g)]
        check(f"per-tool quality covers the 13 configurations {sp}/{g}",
              q.tool.nunique() == 13, f"n={q.tool.nunique()}")
        check(f"every configuration's calls are covered (>= 20 reads) {sp}/{g}",
              float(q.coverage_median.min()) >= 20.0,
              f"min median={q.coverage_median.min():.0f} reads")
        check(f"DRACH spread across configurations {sp}/{g}",
              float(q.frac_drach.min()) <= 0.10 and float(q.frac_drach.max()) >= 0.95,
              f"{100 * q.frac_drach.min():.1f}-{100 * q.frac_drach.max():.1f}%")

    # ---- 1d. panel D: control burden per tool ------------------------------ #
    for cname in ("Curlcake IVT", "HeLa IVT"):
        s = ctltool[ctltool.control == cname]
        check(f"{cname} per-tool rows present for all available tools",
              len(s) > 0 and set(s.tool) <= set(toolq.tool.unique())
              and s.tool.nunique() >= 10,
              f"{s.tool.nunique()} tools, samples={sorted(s['sample'].unique())}")
    check("Curlcake IVT lacks exactly the three HeLa-only configurations",
          set(toolq.tool.unique()) - set(ctltool[ctltool.control == "Curlcake IVT"]
                                         .tool.unique()) == NO_CURLCAKE,
          f"missing="
          f"{sorted(set(toolq.tool.unique()) - set(ctltool[ctltool.control == 'Curlcake IVT'].tool.unique()))}")
    med_curl = float(ctltool[ctltool.control == "Curlcake IVT"].fp_per_10kb.median())
    med_hela = float(ctltool[ctltool.control == "HeLa IVT"].fp_per_10kb.median())
    check("Curlcake IVT is the harsher control (median FP/10 kb)",
          med_curl > med_hela, f"{med_curl:.2f} vs {med_hela:.2f}")

    # ---- 1e. panel E: site-set DRACH --------------------------------------- #
    check("site-quality table carries the five plotted sets of every group",
          all({"single (k=1)", "union marginal", "intersection (k=2)",
               "intersection (k=5)", "GLORI reference"}
              <= set(qual[(qual.species == sp) & (qual.group == g)].set)
              for sp, g in GROUPS))
    order = ["single (k=1)", "union marginal", "intersection (k=2)",
             "intersection (k=5)", "GLORI reference"]
    for (sp, g), want in SET_DRACH.items():
        s = qual[(qual.species == sp) & (qual.group == g)]
        got = tuple(100 * float(s[s.set == name].frac_drach.mean())
                    for name in order)
        check(f"set DRACH means {sp}/{g} match the frozen table",
              all(abs(a - b) <= 0.05 for a, b in zip(got, want)),
              " / ".join(f"{v:.2f}" for v in got))
    inter = qual[qual.set.str.startswith("intersection")]
    check("both intersections are >= 99.9% DRACH in every unit",
          bool((inter.frac_drach.to_numpy(dtype=float) >= 0.999).all()),
          f"min={inter.frac_drach.min():.5f}")
    uni = qual[qual.set == "union (k=5)"].set_index("unit")
    sin = qual[qual.set == "single (k=1)"].set_index("unit")
    mar = qual[qual.set == "union marginal"].set_index("unit")
    check("marginal size = union - single for every unit",
          all(int(mar.loc[u].n_sites) == int(uni.loc[u].n_sites)
              - int(sin.loc[u].n_sites) for u in uni.index))
    for sp, g in GROUPS:
        f = sel[(sel.species == sp) & (sel.group == g)].set_index("k")
        q5 = qual[(qual.species == sp) & (qual.group == g)
                  & (qual.set == "union (k=5)")]
        check(f"union PPV {sp}/{g} reproduces the frozen selection table",
              close(q5.frac_glori.mean(), f.loc[5].union_precision_mean, rtol=1e-5),
              f"{q5.frac_glori.mean():.6f}")

    # ---- 2. figure geometry (rebuild in memory) ---------------------------- #
    spec = importlib.util.spec_from_file_location("figS5",
                                                  HERE.parent / "62_figS5_figure.py")
    f62 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(f62)
    tables = f62.load()
    tools = f62.tool_order(tables["members"], tables["space"])
    f62.check_tables(tables, tools, logging.getLogger("63_verify_figS5"))
    fig, meta = f62.build_figure(tables, tools)
    rep = f62.layout_report(fig)

    check("min font size >= 7.0 pt", rep["min_fontsize"] >= f62.MIN_PT,
          f"min={rep['min_fontsize']:.2f} pt")
    check("no text collisions", not rep["overlaps"],
          f"overlaps={len(rep['overlaps'])}")
    check("no text running off the canvas", not rep["clipped"],
          f"clipped={rep['clipped'][:3]}")
    legs = [lg.get_window_extent(renderer=fig.canvas.get_renderer())
            for lg in fig.legends]
    check("legend strip inside the canvas",
          bool(legs) and all(float(b.x0) >= -1.0
                             and float(b.x1) <= rep["canvas_px"][0] + 1.0
                             for b in legs), f"{len(legs)} legend(s)")
    grid_on = [ax for ax in fig.axes
               if any(t.get_visible() for t in ax.xaxis.get_gridlines())
               or any(t.get_visible() for t in ax.yaxis.get_gridlines())]
    check("no gridlines in any panel", not grid_on,
          f"axes with grid={len(grid_on)}")

    texts = [str(t.get_text()) for t in fig.findobj() if hasattr(t, "get_text")]
    letters = sorted(s for s in texts if s in list("ABCDE"))
    check("panel letters A-E exactly once each", letters == list("ABCDE"),
          f"{letters}")

    # rows B-E are single full-width panels; row A keeps three species columns
    wide = [ax for ax in fig.axes
            if close(ax.get_position().width, (f62.CANVAS_W - f62.L_LEFT
                                               - f62.L_RIGHT) / f62.CANVAS_W,
                     rtol=1e-3)]
    check("rows B-E are four full-width panels", len(wide) == 4, f"n={len(wide)}")
    check("7 axes = three column axes (row A) + four full-width rows (B-E)",
          len(fig.axes) == 7, f"axes={len(fig.axes)}")
    x0s = sorted({round(ax.get_position().x0, 6) for ax in fig.axes
                  if ax not in wide})
    check("three group columns in row A (Arabidopsis | Mouse | HeLa)",
          len(x0s) == 3, f"columns={len(x0s)}")
    row_a = sorted([ax for ax in fig.axes if ax not in wide],
                   key=lambda a: a.get_position().x0)
    check("row A carries the three species column titles",
          [a.get_title() for a in row_a] == ["Arabidopsis", "Mouse", "HeLa"],
          f"{[a.get_title() for a in row_a]}")

    # user decision 2026-09-21: exactly one kind of text is slanted -- the block
    # of 13 tool names under row D at 45 deg (parallel to one another).  Any
    # other rotated label is a regression; the k values of row A and the set
    # names of row E must stay horizontal.
    rot = rep["rotated"]
    check("the only rotated labels are the 13 tool names at 45 deg",
          len(rot) == len(tools)
          and all(abs(float(a) - 45.0) < 1e-6 and s in set(tools) for s, a in rot),
          f"{len(rot)} rotated: {sorted({a for _, a in rot})}")

    # the 13 tool names are printed once, horizontally, under row D, and the
    # five set names once, horizontally, under row E
    row_axes = sorted(wide, key=lambda a: -a.get_position().y0)   # B, C, D, E
    tool_axis, set_axis = row_axes[2], row_axes[3]
    printed = printed_texts(fig)
    tool_labels = list(tool_axis.get_xticklabels())
    check("the 13 tool names are row D's x tick labels",
          drawn_xticklabels(tool_axis) == tools, f"n={len(tool_labels)}")
    check("every tool name is slanted 45 deg and right-aligned at its tick",
          len(tool_labels) == len(tools)
          and all(abs(float(t.get_rotation()) - 45.0) < 1e-6 for t in tool_labels)
          and all(t.get_ha() == "right" for t in tool_labels),
          f"rotations={sorted({float(t.get_rotation()) for t in tool_labels})} "
          f"ha={sorted({t.get_ha() for t in tool_labels})}")
    check("row D carries no axis title (the names are the axis text)",
          not str(tool_axis.get_xlabel()).strip(), tool_axis.get_xlabel())
    check("row A's k values and row E's set names stay horizontal",
          all(abs(float(t.get_rotation())) < 1e-6
              for ax in (*row_a, set_axis) for t in ax.get_xticklabels()),
          "only the tool names may be slanted")
    check("every tool name appears exactly once in the figure",
          sorted(s for s in printed if s in set(tools)) == sorted(tools),
          f"{sum(1 for s in printed if s in set(tools))} occurrences (expect 13)")
    check("row E carries the five two-line set names, horizontally",
          drawn_xticklabels(set_axis) == f62.SET_LABELS,
          f"{drawn_xticklabels(set_axis)[:1]}")
    check("the five set names are printed once, under row E",
          sorted(s for s in printed if s in set(f62.SET_LABELS))
          == sorted(f62.SET_LABELS))
    check("rows B and C share the tool axis but print no labels",
          all(s == "" for s in drawn_xticklabels(row_axes[0]))
          and all(s == "" for s in drawn_xticklabels(row_axes[1])),
          f"{drawn_xticklabels(row_axes[0])[:2]}")

    # rows B, C and E carry the four independent groups as separate clusters
    for idx, name in ((0, "B (coverage)"), (1, "C (DRACH)"), (3, "E (site sets)")):
        colours = collection_colours(row_axes[idx])
        missing = [g for g in f62.GROUPS
                   if rgb(f62.GROUP_STYLE[g]["colour"]) not in colours]
        check(f"row {name} carries the four independent groups", not missing,
              f"missing={missing}")
        check(f"row {name} keeps mouse study B open (white fill)",
              (1.0, 1.0, 1.0) in colours, f"{len(colours)} colours drawn")

    # the mouse column of row A carries both studies: solid AND dashed data
    # series (matplotlib normalises the (0,(4,2)) pattern, hence the tuple)
    mouse_x0 = x0s[1]
    mouse_line = [ax for ax in fig.axes
                  if round(ax.get_position().x0, 6) == mouse_x0]
    dashed_ls = {"--", "-.", ":"}
    solid = dashed = 0
    for ax in mouse_line:
        for ln in ax.get_lines():
            if ln.get_color() not in (f62.UNION_C, f62.ISECT_C):
                continue
            if ln.get_linestyle() == "-":
                solid += 1
            elif isinstance(ln.get_linestyle(), tuple) or ln.get_linestyle() in dashed_ls:
                dashed += 1
    check("mouse column carries solid (study A) and dashed (study B) series",
          solid >= 4 and dashed >= 4, f"solid={solid} dashed={dashed}")

    # ---- 3. no main-figure content leaks back in --------------------------- #
    n_pts = sum(len(c.get_offsets()) for ax in fig.axes for c in ax.collections)
    check("no dense combination cloud (the search space is Fig. 6D, not S5)",
          n_pts < 1000, f"{n_pts} scatter points")
    check("the site-quality dots are actually drawn", n_pts > 200,
          f"{n_pts} scatter points")
    banned_fig = {"chance level": re.compile(r"chance level", re.I),
                  "p0": re.compile(r"p0|p_0"),
                  "search space": re.compile(r"search space", re.I),
                  "selected optimum": re.compile(r"selected optimum", re.I),
                  "all enumerated": re.compile(r"all enumerated", re.I),
                  "trajector": re.compile(r"trajector", re.I),
                  "single-tool plane": re.compile(r"single-tool plane", re.I)}
    leak = [f"{k}: {[s for s in texts if p.search(s)][:2]}"
            for k, p in banned_fig.items() if any(p.search(s) for s in texts)]
    check("no criterion / main-figure wording inside the figure", not leak,
          f"{leak}")

    # ---- 4. PDF / PNG facts ------------------------------------------------ #
    try:
        out = subprocess.run(["pdffonts", str(PDF)], capture_output=True,
                             text=True, timeout=60).stdout
        lines = [l for l in out.splitlines()[2:] if l.strip()]
        names = {l.split()[0].split("+")[-1] for l in lines}
        check("PDF fonts all Arial",
              names <= {"ArialMT", "Arial-BoldMT", "Arial-ItalicMT",
                        "Arial-BoldItalicMT"}, f"fonts={sorted(names)}")
        check("PDF fonts embedded", bool(lines)
              and all(l.split()[-5] == "yes" for l in lines))
    except FileNotFoundError:
        check("pdffonts available", False, "pdffonts not on PATH")
    try:
        out = subprocess.run(["pdfinfo", str(PDF)], capture_output=True,
                             text=True, timeout=60).stdout
        mm = re.search(r"Page size:\s+([\d.]+) x ([\d.]+) pts", out)
        pages = re.search(r"Pages:\s+(\d+)", out)
        w_mm = float(mm.group(1)) * 25.4 / 72 if mm else 0.0
        h_mm = float(mm.group(2)) * 25.4 / 72 if mm else 0.0
        check(f"PDF is one {f62.CANVAS_W * 25.4:.0f} x "
              f"{f62.CANVAS_H * 25.4:.0f} mm page, printed 1:1",
              pages and pages.group(1) == "1"
              and abs(w_mm - f62.CANVAS_W * 25.4) <= 1.0
              and abs(h_mm - f62.CANVAS_H * 25.4) <= 1.0,
              f"{w_mm:.1f} x {h_mm:.1f} mm, pages={pages.group(1) if pages else '?'}")
        check("canvas width equals the PDF page width (no rescaling)",
              abs(w_mm - f62.CANVAS_W * 25.4) <= 0.5,
              f"{w_mm:.1f} vs {f62.CANVAS_W * 25.4:.1f} mm")
    except FileNotFoundError:
        check("pdfinfo available", False, "pdfinfo not on PATH")
    if PNG.exists():
        from PIL import Image
        with Image.open(PNG) as im:
            check("300 dpi PNG on the A4 canvas",
                  abs(im.size[0] - round(f62.CANVAS_W * 300)) <= 4
                  and abs(im.size[1] - round(f62.CANVAS_H * 300)) <= 4,
                  f"{im.size}")
    else:
        check("PNG written", False, str(PNG.name))

    # ---- 5. legend file + wording residue ---------------------------------- #
    banned = {"treatment": re.compile(r"treatment", re.I),
              "T/WT": re.compile(r"(?<![A-Za-z])T/WT"),
              "hit rate": re.compile(r"hit rate", re.I),
              "GENCODE": re.compile(r"gencode", re.I),
              "ORCA (wrong reference name for Fig. 7)": re.compile(r"\bORCA\b"),
              "NanoSPA_Psu": re.compile(r"NanoSPA_Psu"),
              "NanoMUD_psi": re.compile(r"NanoMUD_psi"),
              "the retired ten-panel lettering": re.compile(r"A-J|A–J")}
    residue = [f"{k}: {[s for s in texts if p.search(s)][:2]}"
               for k, p in banned.items() if any(p.search(s) for s in texts)]
    if LEGEND_MD.exists():
        body = LEGEND_MD.read_text(encoding="utf-8")
        residue += [f"legends.md contains {k}" for k, p in banned.items()
                    if p.search(body)]
        for letter, role in (("A", "greedy"), ("B", "coverage"),
                             ("C", "DRACH"), ("D", "control"),
                             ("E", "DRACH")):
            check(f"legend documents panel {letter} ({role})",
                  re.search(rf"\*\*\(?{letter}\)?\*\*", body) is not None
                  or re.search(rf"\(?{letter}\)?[ .—-]+{role}", body, re.I)
                  is not None,
                  f"panel {letter}")
        check("legend points the criterion at the main figure and the table",
              "Fig. 6D" in body and "figS5_search_space.tsv" in body)
        check("legend states that p0 was satisfied by every combination",
              re.search(r"(every|all).{0,60}combinat", body, re.I) is not None
              and re.search(r"p0|p₀|chance", body, re.I) is not None)
        check("legend states the greedy path reproduced the optima",
              re.search(r"greedy", body, re.I) is not None
              and re.search(r"reproduc", body, re.I) is not None)
        check("legend quotes the intersection cost at k = 5",
              all(v in body for v in ("0.12", "2.24", "1.74", "1.58")),
              "0.12 / 2.24 / 1.74 / 1.58 %")
        check("legend states the two mouse studies are never averaged",
              "never averaged" in body)
        check("legend quotes the chance level p0 values",
              "1.5" in body and "0.7" in body and "0.85" in body)
        if LEGEND_MIRROR.exists():
            check("Fig. 6 directory mirror is byte-identical",
                  LEGEND_MIRROR.read_bytes() == LEGEND_MD.read_bytes(),
                  str(LEGEND_MIRROR.relative_to(PROJECT)))
        else:
            check("Fig. 6 directory mirror exists", False,
                  str(LEGEND_MIRROR.relative_to(PROJECT)))
    else:
        check("legend file written", False, str(LEGEND_MD.name))
    check("no retired naming / wording residue", not residue, f"{residue}")

    LOG.mkdir(parents=True, exist_ok=True)
    print(f"\n{'ALL CHECKS PASSED' if not FAILURES else 'FAILURES: ' + str(FAILURES)}",
          flush=True)
    (LOG / "63_verify_figS5.log").write_text(
        "failures: " + (", ".join(FAILURES) if FAILURES else "none") + "\n")
    sys.exit(1 if FAILURES else 0)


if __name__ == "__main__":
    main()
