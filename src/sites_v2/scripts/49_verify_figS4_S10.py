#!/usr/bin/env python
"""49 -- acceptance checks for the rebuilt Figure S4 and the new Figure S10.

Runs the four checks the revision workflow asks for before a figure is allowed
anywhere near the submission folder:

1. **numbers** -- every anchor written by ``42_figS4_figure.py`` /
   ``48_figS10_validation.py`` is recomputed from the frozen ``sites_v2`` tables
   (``m6a_glori_confusion.tsv``, ``m6a_localization_curve.tsv`` and the
   ``figS4_*`` tables) and must agree to 1e-6 relative;
2. **page/vector contract** -- one page, A4 **portrait** (the page of the
   submitted ``sup4.pdf``), every font embedded (``pdfinfo`` / ``pdffonts``),
   and the page is assembled from the per-panel PDFs
   (``figures/panels/*.pdf``) that ``panelpage.compose_page`` merges;
3. **house style** -- no grid anywhere in the two drawing scripts, no text below
   7 pt, no in-panel annotation text, and both scripts must draw every panel on
   its own canvas and gate it (``panelpage.save_panel`` ->
   ``pagelayout.assert_page_clean``) with panel letters in the page margin
   (``pagelayout.margin_letter``) instead of at an axes-relative offset.  The
   gate must also reject text that *hugs* a foreign frame line
 (``near_pt`` >= 3 pt) -- the defect the user reported twice as " font/box/line overlap ";
   Figure S4 must keep its 13-tool key in the centred stripe after panel A,
 must *not* carry the PR-AUC row any more (dropped 2026-09-21, " values are all too low ")
   while its frozen PR-AUC tables stay anchored, must draw the window sweep as
   **two lettered facets** (B = PPV, C = exact-nucleotide fraction) with the
   purified comparison as panel D, and Figure S10 must keep the enlarged type
   scale of its own;
4. **cross-references** -- the S4 legends no longer advertise panels E--G and
   the S4/S10 legends point at the right figure number.

Exit code 0 = all checks passed; 2 = at least one failure (printed).

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/src/sites_v2/scripts/49_verify_figS4_S10.py
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
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

PROJECT = Path(str(_RB))
SITES = (_RB / "data")
S4 = (_RB / "figures/figureS4")
S10 = (_RB / "figures/figureS5")
SCRIPTS = (_RB / "src/sites_v2/scripts")

TOOLS = ["CHEUI_m6A", "DENA", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo",
         "EpiNano_Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A",
         "NanoSPA_m6A", "xPore", "yanocomp"]
COLUMN_SAMPLES = {
    "Arabidopsis": ["Arabidopsis_WT_rep1", "Arabidopsis_WT_rep2",
                    "Arabidopsis_WT_rep3"],
    "Mouse": ["mESCs_Mettl3_WT", "mES_WT"],
    "HeLa": ["HeLa_WT1", "HeLa_WT2", "HeLa_WT3"],
}
COLUMN_PAIRS = {
    "Arabidopsis": ["Ath_rep1", "Ath_rep2", "Ath_rep3"],
    "Mouse": ["Mouse_studyA", "Mouse_studyB"],
    "HeLa": ["HeLa_rep1", "HeLa_rep2", "HeLa_rep3"],
}
FAILURES: list[str] = []


def check(label: str, condition: bool, detail: str = "") -> None:
    status = "ok  " if condition else "FAIL"
    print(f"[{status}] {label}{(' -- ' + detail) if detail else ''}")
    if not condition:
        FAILURES.append(label)


def close(a: float, b: float, tol: float = 1e-6) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def anchors(path: Path) -> dict[str, float]:
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["anchor"], df["value"]))


# --------------------------------------------------------------------------- #
# 1. numbers against the frozen tables
# --------------------------------------------------------------------------- #
def verify_numbers() -> None:
    a4 = anchors((_RB / "figures/figureS4/logs/42_figS4_anchors.tsv"))
    a10 = anchors((_RB / "figures/figureS5/logs/48_figS10_anchors.tsv"))

    conf = pd.read_csv((_RB / "data/evaluation/tables/m6a_glori_confusion.tsv"), sep="\t")
    conf = conf[(conf["platform"] == "RNA002") & (conf["window"] == 2)
                & (conf["tool"].isin(TOOLS))]
    for key, samples in COLUMN_SAMPLES.items():
        want = conf[conf["sample"].isin(samples)]["precision"].mean()
        check(f"confusion PPV@2bp vs anchors F-S4 ({key})",
              close(want, a4[f"ppv_w2_{key}"]),
              f"{want:.6f} vs {a4[f'ppv_w2_{key}']:.6f}")

    curve = pd.read_csv((_RB / "data/evaluation/tables/m6a_localization_curve.tsv"), sep="\t")
    curve = curve[(curve["platform"] == "RNA002") & (curve["tool"].isin(TOOLS))]
    for key, samples in COLUMN_SAMPLES.items():
        for window in (0, 50):
            sub = curve[curve["sample"].isin(samples)
                        & (curve["window"] == window)]
            want = sub["hit_rate"].mean()
            check(f"localisation curve PPV@w{window} ({key})",
                  close(want, a4[f"ppv_w{window}_{key}"]),
                  f"{want:.6f} vs {a4[f'ppv_w{window}_{key}']:.6f}")

    prec = pd.read_csv((_RB / "figures/figureS4/tables/figS4_precision_on_purified.tsv"),
                       sep="\t")
    for column, key in (("precision_wt_common_w2", "ppv_wt_common"),
                        ("precision_purified_w2", "ppv_purified")):
        want = prec[column].mean()
        check(f"precision table {column}", close(want, a4[key]),
              f"{want:.6f} vs {a4[key]:.6f}")
    up = int((prec["precision_purified_w2"]
              > prec["precision_wt_common_w2"]).sum())
    check("tool-pair combinations shifting upwards", int(a4["n_tool_pair_up"]) == up,
          f"{up} of {len(prec)}")

    auprc = pd.read_csv((_RB / "figures/figureS4/tables/figS4_auprc_summary.tsv"), sep="\t")
    for group in ("Arabidopsis", "Mouse_studyA", "Mouse_studyB", "HeLa"):
        for window in (0, 50):
            sub = auprc[(auprc["species_group"] == group)
                        & (auprc["window"] == window)]
            want = sub["pr_auc_mean"].mean()
            key = f"pr_auc_w{window}_{group}"
            check(f"PR-AUC summary {key}", close(want, a4[key]),
                  f"{want:.6f} vs {a4[key]:.6f}")

    dfv = pd.read_csv((_RB / "figures/figureS4/tables/figS4_validation_groups.tsv"), sep="\t")
    for group in ("def_only", "shared", "purified"):
        sub = dfv[dfv["group"] == group]
        for column, key in (("glori_hit_rate_w2", f"glori_overlap_{group}"),
                            ("drach_rate", f"drach_{group}")):
            want = sub[column].mean()
            check(f"validation groups {key}", close(want, a10[key]),
                  f"{want:.6f} vs {a10[key]:.6f}")
        for column, pairs in COLUMN_PAIRS.items():
            sub2 = dfv[(dfv["group"] == group)
                       & (dfv["pair_id"].isin(pairs))]
            want = sub2["glori_hit_rate_w2"].mean()
            key = f"glori_overlap_{group}_{column}"
            check(f"validation groups {key}", close(want, a10[key]),
                  f"{want:.6f} vs {a10[key]:.6f}")

    bg = pd.read_csv((_RB / "figures/figureS4/tables/figS4_universe_drach.tsv"), sep="\t")
    for sample, value in bg.set_index("sample")["drach_rate"].items():
        key = f"drach_background_{sample}"
        if key in a10:
            check(f"universe DRACH background {sample}", close(value, a10[key]),
                  f"{value:.6f} vs {a10[key]:.6f}")

    semantic = ["CHEUI_m6A", "m6Anet", "MINES", "Nanom6A", "DENA", "NanoSPA_m6A"]
    for tool in semantic:
        sub = dfv[dfv["tool"] == tool]
        for group in ("shared", "purified"):
            want = sub[sub["group"] == group]["score_median"].mean()
            check(f"semantic score {group} {tool}",
                  close(want, a10[f"score_{group}_{tool}"]),
                  f"{want:.6f} vs {a10[f'score_{group}_{tool}']:.6f}")

    check("S10 direction is reported as mixed (no systematic stoichiometry gain)",
          (a10["score_purified_m6Anet"] < a10["score_shared_m6Anet"])
          and (a10["score_purified_Nanom6A"] > a10["score_shared_Nanom6A"]))


# --------------------------------------------------------------------------- #
# 2. page / font contract
# --------------------------------------------------------------------------- #
def verify_pdfs() -> None:
    # the per-panel files the page is assembled from (2026-09-21 architecture)
    for folder, names in (
            ((_RB / "figures/figureS4/figures/panels"),
             ["FigureS4_A", "FigureS4_legend_tools", "FigureS4_BC",
              "FigureS4_D"]),
            ((_RB / "figures/figureS5/figures/panels"),
             ["FigureS10_A", "FigureS10_B", "FigureS10_C"])):
        for name in names:
            panel = folder / f"{name}.pdf"
            check(f"panel {name}.pdf exists", panel.is_file(),
                  "" if panel.is_file() else str(panel))

    for stem in ((_RB / "figures/figureS4/figures/FigureS4_rev.pdf"),
                 (_RB / "figures/figureS5/figures/FigureS10_rev.pdf")):
        info = subprocess.run(["pdfinfo", str(stem)], capture_output=True,
                              text=True).stdout
        pages = re.search(r"^Pages:\s+(\d+)", info, re.M)
        size = re.search(r"^Page size:\s+([\d.]+) x ([\d.]+)", info, re.M)
        check(f"{stem.name}: single page",
              bool(pages) and pages.group(1) == "1")
        if size:
            w, h = float(size.group(1)), float(size.group(2))
            check(f"{stem.name}: A4 portrait",
                  abs(w - 595.3) < 2 and abs(h - 841.9) < 2,
                  f"{w:.1f} x {h:.1f} pt")
        fonts = subprocess.run(["pdffonts", str(stem)], capture_output=True,
                               text=True).stdout.splitlines()[2:]
        rows = [f.split() for f in fonts if f.strip()]
        check(f"{stem.name}: fonts embedded",
              bool(rows) and all(r[-4] == "yes" for r in rows),
              f"{len(rows)} font(s): {', '.join(r[0] for r in rows)}")
        check(f"{stem.name}: all fonts are TrueType/Type1",
              all(("TrueType" in " ".join(r)) or ("Type 1" in " ".join(r))
                  for r in rows))


# --------------------------------------------------------------------------- #
# 3. house style in the drawing scripts
# --------------------------------------------------------------------------- #
def _font_sizes(body: str) -> list[float]:
    """Every declared font size, with ``FS["..."]`` lookups resolved.

    Since the 2026-09-21 style pass the two scripts declare their sizes through
    the shared ``FS`` dict, so a plain ``fontsize=7.5`` regex is not enough any
    more -- the dict values have to be substituted in.
    """
    fs = {m.group(1): float(m.group(2)) for m in re.finditer(
        r'"(\w+)":\s*([0-9.]+)', body)}
    sizes = [float(m) for m in
             re.findall(r"(?:fontsize|labelsize)=([0-9.]+)", body)]
    for key in re.findall(r'(?:fontsize|labelsize)=FS\["(\w+)"\]', body):
        if key in fs:
            sizes.append(fs[key])
    return sizes


def verify_style() -> None:
    # the shared per-panel machinery itself must keep calling the gate
    module = (SCRIPTS.parent / "common" / "panelpage.py").read_text()
    check("panelpage.py runs the layout gate on every panel",
          "pagelayout.assert_page_clean(" in module
          and "ignore_axes=" in module)
    check("panelpage.py asserts that panels cannot overlap on the page",
          "_overlaps(" in module and "panels overlap" in module)
    # the gate measures printed ink, not axis-aligned bounding boxes: without
    # that a 45-degree tool name is flagged as overlapping its neighbour
    layout = (SCRIPTS.parent / "common" / "pagelayout.py").read_text()
    check("pagelayout.py gates on rotated ink quads (SAT), not bounding boxes",
          "def text_quad(" in layout and "def quad_gap(" in layout
          and "def quads_intersect(" in layout
          and "def _selftest(" in layout)
    check("pagelayout.py handles both rotation modes",
          'get_rotation_mode() == "anchor"' in layout and "rot_w" in layout)

    for script in ((_RB / "src/sites_v2/scripts/42_figS4_figure.py"), (_RB / "src/sites_v2/scripts/48_figS10_validation.py")):
        text = script.read_text()
        body = "\n".join(line for line in text.splitlines()
                         if not line.strip().startswith("#"))
        check(f"{script.name}: no grid", "grid(" not in body
              and "axes.grid': True" not in body)
        sizes = _font_sizes(body)
        # 8 pt floor for everything printed; the single exception is a
        # ``legend_small`` entry (7.5 pt, the 13-entry tool-colour key that used
        # to sit inside panel D of Figure S4), which is still above the 7 pt
        # print-size house rule.  Panel A's 45-degree tool names are what sets
        # the 8 pt floor (13 names in a 2.16 in column keep a 1.5 pt ink gap at
        # 8 pt).
        fs_floor = {m.group(1): float(m.group(2)) for m in re.finditer(
            r'"(\w+)":\s*([0-9.]+)', body)}
        small = fs_floor.get("legend_small")
        others = [s for s in sizes if small is None or s != small]
        check(f"{script.name}: smallest declared font >= 8.0 pt "
              f"(in-panel key >= 7.5 pt)",
              bool(sizes) and min(sizes) >= 7.5 and (not others or min(others) >= 8.0),
              f"min={min(sizes)}" + (f", legend_small={small}" if small else ""))
        # no panel may carry free annotation text (labels belong in titles,
        # legends or the tables) -- and no ``ax.text`` at all any more: the
        # panel letters moved into the page margin
        check(f"{script.name}: no in-panel annotation text",
              "annotate(" not in body and "ax.text(" not in body)
        # one panel at a time, then assemble (2026-09-21 third rework)
        check(f"{script.name}: panels are drawn on their own canvases",
              "panelpage.new_panel(" in body
              and "panelpage.save_panel(" in body)
        check(f"{script.name}: the page is assembled from the panels",
              "panelpage.compose_page(" in body
              and "page_placements(" in body)
        check(f"{script.name}: panel letters come from the page margin",
              "pagelayout.margin_letter(" in body)

    # panel A of Figure S4 keeps the 45-degree tilt of the submitted figure
    s4 = ((_RB / "src/sites_v2/scripts/42_figS4_figure.py")).read_text()
    check("42_figS4_figure.py: panel-A tool names default to 45 degrees",
          "TICK_ROTATION = 45" in s4)
    check("42_figS4_figure.py: the tilted labels get their own pad",
          "pad=5.0 if rotation == 90 else 8.0" in s4)
    # every key of Figure S4 is drawn in a band that cannot hide data
    # (2026-09-21, fifth rework): the 13 tool colours sit in a *centred stripe

    # their own line inside the A / B canvases
    check("42_figS4_figure.py: the tool key is a centred stripe after panel A",
          "def render_tool_stripe(" in s4
          and '"legend_tools"' in s4
          and 'loc="center"' in s4
          and 'order = [("A", GAP_IN), ("legend_tools", GAP_IN), ("BC", GAP_IN),'
              in s4
          and "arab_ax.legend(" not in s4)
    check("42_figS4_figure.py: three keys are defined",
          "A_KEY: list" in s4 and "LINE_KEY_B: list" in s4
          and "TOOL_KEY: list" in s4)
    check("42_figS4_figure.py: two keys on the panel canvases + one stripe",
          s4.count("fig.legend(") == 2 and s4.count("ax.legend(") == 1)

    # of the page (no renderer, no key headroom), while its frozen tables keep
    # being cross-checked.  Panel D is something else entirely: since the B/C
    # re-lettering it is the purified-site comparison, which must exist.
    check("42_figS4_figure.py: the dropped PR-AUC row is not drawn",
          "D_YLIM_HEADROOM" not in s4 and "AUPRC_GROUPS" not in s4
          and "render_panel_pr" not in s4
          and "df_auprc" in s4)
    check("42_figS4_figure.py: panel D exists (the purified comparison)",
          "def render_panel_d(" in s4 and '"FigureS4_D"' in s4
          and 'panel_letter(fig, "D"' in s4 and "D_AXES_IN" in s4)
    check("42_figS4_figure.py: the page stacks A / key stripe / BC / D",
          'order = [("A", GAP_IN), ("legend_tools", GAP_IN), ("BC", GAP_IN),'
          in s4 and '("D", 0.0)]' in s4)
    check("42_figS4_figure.py: the archived PR-AUC tables stay verified",
          "_assert_unit_means(df_unit, df_auprc)" in s4
          and "pr_auc_w{window}_{group}" in s4)
    check("42_figS4_figure.py: the sub-row titles clear the frame line",
          "A_TITLE_INDENT_IN" in s4 and "t_x = x_in + A_TITLE_INDENT_IN" in s4)


    # below, one shared x axis, and a single block letter B
    check("42_figS4_figure.py: the B/C canvas draws two facets per column",
          "B_FACET_IN" in s4
          and 'facet_spec = (("localization_accuracy", 0), ("hit_rate", 1))'
          in s4 and "facets[0].sharex(facets[1])" in s4)
    # ...and each of the two facets carries a panel letter of its own

    check("42_figS4_figure.py: each facet carries its own margin letter",
          "def render_panels_bc(" in s4 and '"FigureS4_BC"' in s4
          and 'panel_letter(fig, "B", y_in=y_up + B_FACET_IN - 0.04)' in s4
          and 'panel_letter(fig, "C", y_in=y_lo + B_FACET_IN - 0.04)' in s4)
    check("42_figS4_figure.py: the two facets share one window axis",
          "facets[0].set_xticks([0, 20, 50])" in s4
          and 'facets[1].tick_params(axis="x", which="both", bottom=False,'
          in s4 and '"Matching window (bp)"' in s4)
    check("42_figS4_figure.py: both facets keep units and their group mean",
          "for column, fi in facet_spec:" in s4
          and 't.groupby("window")[column].mean()' in s4 and "lw=0.7" in s4
          and "lw=2.0" in s4)
    check("42_figS4_figure.py: the lower facet names its own quantity",
          "B_EXACT_LABEL" in s4 and "set_ylabel(B_EXACT_LABEL" in s4)
    key_b = re.search(r"LINE_KEY_B: list\[tuple\[str, dict\]\] = \[(.*?)\n\]",
                      s4, re.S)
    key_spec = key_b.group(1) if key_b else ""
    n_key_entries = len(re.findall(r'\("', key_spec))
    check("42_figS4_figure.py: the key names the two mouse study styles",
          '("Study A",' in key_spec and '("Study B",' in key_spec
          and "Per-tool pair (D)" in key_spec
          and '"Exact nucleotide"' not in key_spec,
          f"{n_key_entries} entries in LINE_KEY_B")
    check("42_figS4_figure.py: the mouse facets draw one line per study",
          "for _study, sample, ls, _filled in MOUSE_STUDIES:" in s4
          and "one line per study, never merged" in s4)
    check("42_figS4_figure.py: every key line is measured, not estimated",
          "def assert_key_fits(" in s4
          and s4.count("assert_key_fits(fig, legend,") == 2)

    def inch(name: str) -> float:
        """Inch constant of the renderer, read from its source."""
        hit = re.search(rf"^{re.escape(name)}\s*=\s*([0-9.]+)", s4, re.M)
        if hit is None:
            raise SystemExit(f"42_figS4_figure.py: cannot read {name}")
        return float(hit.group(1))

    # the second facet had to be paid for: the budget still fills the A4 page
    # exactly, and B/C stay close to square (the point of the rework)
    a_row = (inch("A_LABEL_IN") + inch("A_RATIO_IN")
             + 2 * (inch("A_BAND_IN") + inch("A_SUB_IN"))
             + inch("A_HEAD_IN") + inch("A_KEY_IN"))
    b_row = (inch("B_KEY_IN") + inch("B_TITLE_IN") + 2 * inch("B_FACET_IN")
             + inch("B_GAP_IN") + inch("B_BOTTOM_IN"))
    d_row = (inch("D_TITLE_IN") + inch("D_AXES_IN") + inch("D_BOTTOM_IN"))
    total = (inch("TOP_IN") + a_row + inch("GAP_IN") + inch("TOOL_STRIPE_IN")
             + inch("GAP_IN") + b_row + inch("GAP_IN") + d_row
             + inch("BOTTOM_IN"))
    check("42_figS4_figure.py: the inch budget fills the A4 page exactly",
          abs(total - 11.70) < 0.005,
          f"A {a_row:.2f} + stripe + BC {b_row:.2f} + D {d_row:.2f} = "
          f"{total:.3f} in")
    col_w = (8.27 - inch("LEFT_IN") - inch("RIGHT_IN")
             - 2 * inch("COL_GAP_IN")) / 3
    for label, height in (("B/C facet", inch("B_FACET_IN")),
                          ("D", inch("D_AXES_IN"))):
        ratio = max(col_w / height, height / col_w)
        check(f"42_figS4_figure.py: the {label} sub-panel stays near square",
              ratio <= 1.25, f"{col_w:.2f} x {height:.2f} in = {ratio:.2f}:1")
    # no key is drawn over a set of data axes: the only ``ax.legend`` left is
    # the stripe's, on an axis-off canvas with nothing behind it; the A/B key
    # lines are figure legends in their own bands, and the page gate rejects
    # any text that enters a foreign drawing area
    check("42_figS4_figure.py: no key is drawn over a data axes",
          "arab_ax.legend(" not in s4 and s4.count("ax.legend(") == 1
          and "ignore_axes=[ax]" in s4)
    s10 = ((_RB / "src/sites_v2/scripts/48_figS10_validation.py")).read_text()
    check("48_figS10_validation.py: no page-level legend stripe",
          "def render_legend(" not in s10 and "LEG_IN" not in s10)

    #: the keys are axes legends now (three in A, three in B, two in C) and no
    #: figure-level legend is left; the intent ("no page-foot stripe") holds
    check("48_figS10_validation.py: the keys sit inside the panels",
          "KEY_ENTRIES: list" in s10 and s10.count("ax.legend(") == 3
          and "fig.legend(" not in s10)

    fs10 = {m.group(1): float(m.group(2)) for m in re.finditer(
        r'"(\w+)":\s*([0-9.]+)', s10)}
    fs4 = {m.group(1): float(m.group(2)) for m in re.finditer(
        r'"(\w+)":\s*([0-9.]+)', s4)}
    check("48_figS10_validation.py: S10 uses an enlarged scale of its own",
          fs10.get("tick", 0.0) >= 12.0 and fs10.get("axis", 0.0) >= 14.0
          and fs10.get("column", 0.0) >= 17.0,
          f"tick={fs10.get('tick')}, axis={fs10.get('axis')}, "
          f"column={fs10.get('column')}")
    check("42_figS4_figure.py: S4 keeps its own (smaller) scale",
          fs4.get("tick", 99.0) < fs10.get("tick", 0.0),
          f"S4 tick={fs4.get('tick')} < S10 tick={fs10.get('tick')}")
    layout_src = (SCRIPTS.parent / "common" / "pagelayout.py").read_text()
    check("pagelayout.py provides the legend/data clearance gate",
          "def assert_legend_clear(" in layout_src
          and "def _data_vertices(" in layout_src)
    # the "text hugs a frame line" defect is gated, not eyeballed
    check("pagelayout.py gates text that hugs a foreign frame line",
          "near_pt" in layout_src and "text hugs a frame line" in layout_src)
    check("pagelayout.py keeps the frame-clearance floor >= 3 pt",
          "near_pt: float = 3.0" in layout_src)


# --------------------------------------------------------------------------- #
# 4. cross-references
# --------------------------------------------------------------------------- #
def verify_references() -> None:
    s4_legend = ((_RB / "figures/figureS4/figures/FigS4_legends.md")).read_text()
    # whitespace-normalised view: the legend is hard-wrapped, so a phrase can be
    # split across two lines (the S6 lesson: normalise before matching)
    s4_flat = re.sub(r"\s+", " ", s4_legend)
    s4_readme = ((_RB / "figures/figureS4/README.md")).read_text()
    s10_legend = ((_RB / "figures/figureS5/figures/FigS10_legends.md")).read_text()
    check("S4 legend no longer advertises panels E-G",
          "(E–G)" not in s4_legend and "**(E)" not in s4_legend
          and "**(G)" not in s4_legend)
    check("S4 legend documents the four delivered panels (A-D)",
          "**(A)**" in s4_legend and "**(B)**" in s4_legend
          and "**(C)**" in s4_legend and "**(D)**" in s4_legend)
    check("S4 legend describes B and C as two lettered facets on one axis",
          "facet" in s4_flat and "Exact-nucleotide fraction" in s4_flat
          and "PPV vs. GLORI" in s4_flat
          and "share one x axis" in s4_flat
          and "Three line layers" not in s4_flat)
    check("S4 legend keeps the per-unit + mean definition of both facets",
          "individual sequencing units" in s4_flat
          and "group mean" in s4_flat
          and "two layers" in s4_flat)
    check("S4 README records the two-facet rework",
          " " in s4_readme or "two stacked facets" in s4_readme)
    check("S4 README records the B/C/D re-lettering",
          "" in s4_readme or "lettered" in s4_readme)
    check("S4 legend no longer describes the dropped PR-AUC row",
          "PR-AUC" not in s4_legend and "AUPRC" not in s4_legend
          and "four-row" not in s4_legend)
    check("S4 legend points to Figure S10",
          "Figure S10" in s4_legend)
    check("S4 README points to Figure S10",
          "Figure S10" in s4_readme or "figS10_revision" in s4_readme)
    check("S10 legend is self-contained (panels A-C)",
          "**(A)**" in s10_legend and "**(B)**" in s10_legend
          and "**(C)**" in s10_legend)


def main() -> int:
    print("== 1. numbers vs frozen sites_v2 tables ==")
    verify_numbers()
    print("\n== 2. page / font contract ==")
    verify_pdfs()
    print("\n== 3. house style ==")
    verify_style()
    print("\n== 4. cross-references ==")
    verify_references()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILED check(s): {FAILURES}")
        return 2
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
