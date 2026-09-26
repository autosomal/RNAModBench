#!/usr/bin/env python
"""S8 data layer -- frozen tables shaped into the panels of the rebuilt S8.

Nothing is recomputed here except the bootstrap interval of panel E and the
least-squares fits of panel F, which run on frozen matched-pair tables.  Every
other number is read from ``tables/source/`` (see ``SOURCE.md`` there for
provenance).

The page is **RNA004 only** (user decision, 2026-09-21): the RNA002 Curlcake
runs of the original benchmark are not drawn anywhere on this figure.  The
unmodified Curlcake control therefore carries the RNA004 entries only -- the
eight Dorado m6A models at the 50 % operating point plus the two m6A tools that
were run on that library (m6Anet, NanoSPA_m6A).

Panels
------
A  HeLa RNA004 detection counts as a dot matrix (rows = models, columns =
   WT / IVT) in three blocks: the eight Dorado m6A models, the five other m6A
   tools and the ten other-modification models.
B  the eight ORCA channels on the same HeLa RNA004 library, in their own panel
   (ORCA calls several modification types from one run).
C  the unmodified RNA004 Curlcake control: the calls per entry (dot column, the
   same area law as A and B) and the false-positive rate per 10^6 candidate
   sites on a log axis with the matching specificity on the top axis; both
   blocks share one row list, so the labels are printed once.
D  matching window (0-50 bp) against PPV and against the exact-nucleotide
   fraction, all 13 tools.
E  Pearson r against GLORI (Fisher-z 95 % CI), one lollipop per ratio tool.
F  the effect size: the OLS slope of the tool ratio on the GLORI ratio with its
   bootstrap 95 % CI (reference 1.0 = proportional).
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

import numpy as np
import pandas as pd

from s8_style import PROJECT, SRC_TABLES, TABLES

#: frozen window-sweep table of the sites_v2 evaluation (read-only)
LOCALIZATION = ((_RB / "data/evaluation/tables/m6a_localization_curve.tsv"))

#: frozen confusion table of the sites_v2 evaluation (read-only)
CONFUSION = ((_RB / "data/evaluation/tables/m6a_glori_confusion.tsv"))

#: the only chemistry on this page (RNA004); sample labels of the frozen tables
HELA_UNITS = ["HeLa_RNA004_WT", "HeLa_RNA004_IVT"]
CURLAKE_RNA004_UNIT = "Curlcake_RNA004_IVT"

#: Curlcake Dorado call sets are evaluated at this operating point
DORADO_CURLAKE_PCT = 50

#: row blocks of panel A, in print order
BLOCKS_A = ["Dorado m6A models", "other m6A tools", "other modification models"]

#: the single block of panel B and the blocks of panel C, in print order
BLOCK_ORCA = "ORCA channels"
BLOCKS_C = ["m6A DRACH", "m6A (non-DRACH)", "inosine + m6A", "other m6A tools"]

#: Dorado family of the Curlcake run -> (printed block, colour family).  The run
#: names the two v5.1.0 families differently: "all" is that caller's DRACH
#: equivalent (the same key the old HeLa link map used).
CURLAKE_FAMILY = {
    "m6A_DRACH": ("m6A DRACH", "m6A_DRACH"),
    "all": ("m6A DRACH", "m6A_DRACH"),
    "pseU_m6A": ("m6A (non-DRACH)", "m6A"),
    "inosine_m6A": ("inosine + m6A", "inosine_m6A"),
}

#: tools whose score is a p-value, not a modification ratio -- excluded from the
#: correlation and effect-size panels (their matched sets are tiny)
NON_RATIO_TOOLS = {"ELIGOS2_solo", "NanoSPA_m6A"}


#: the v5.1.0 run that loads the whole model set at once is named ``..._all`` in
#: the call tables, and each suffixed channel (``_all_Psi``, ``_all_m5C``) is that
#: run's output for one modification.  The code name never reaches the page: the
#: m6A channel is printed as the DRACH-equivalent model of that release (the same
#: reading as ``CURLAKE_FAMILY`` above), the other channels under their own
#: modification.  This is the single label authority for every S8 panel.
ALL_RUN = {"all": "m6A DRACH", "all_Psi": "\u03a8", "all_m5C": "m5C"}


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(SRC_TABLES / name, sep="\t")


def pretty_tool(tool: str) -> str:
    """Printed label of a tool/model entry (caller version kept for Dorado)."""
    if tool.startswith("ORCA_"):
        return tool.replace("ORCA_", "ORCA ")
    if not tool.startswith("Dorado_"):
        return tool
    s = tool[len("Dorado_"):].replace("_otherMod", "")
    m = re.match(r"(?P<mode>hac|sup)@(?P<ver>v[\d.]+)_?(?P<model>.*)", s)
    if not m:
        return tool.replace("Dorado_", "")
    raw = re.sub(r"@v\d+$", "", m.group("model"))
    if raw in ALL_RUN:
        model = ALL_RUN[raw]
    else:
        model = raw.replace("_", " ").strip()
        model = model.replace("inosine m6A", "inosine+m6A") or "m6A"
    return re.sub(r"\s+", " ", f"{m.group('mode')} {m.group('ver')} {model}")


def family_of(tool: str, mod_type: str) -> str:
    """Colour family key of an entry."""
    if "_DRACH" in tool:
        return "m6A_DRACH"
    if mod_type == "m6A":
        return "inosine_m6A" if "inosine" in tool else "m6A"
    if mod_type == "Psi":
        return "pseU"
    if mod_type == "m5C":
        return "m5C"
    return "tool"


def _block_of_a(tool: str, mod_type: str) -> str:
    """Row block of panel A."""
    if mod_type == "m6A":
        return "Dorado m6A models" if tool.startswith("Dorado_") else "other m6A tools"
    return "other modification models"


# --------------------------------------------------------------------------- #
# A -- HeLa RNA004 detection counts (dot matrix, WT | IVT)
# --------------------------------------------------------------------------- #
def panel_a() -> pd.DataFrame:
    """Rows = HeLa RNA004 entries, columns = WT / IVT, value = detected sites.

    The ORCA channels are **not** rows of this panel (they have their own panel
    B); the remaining entries are grouped into the three ``BLOCKS_A`` and sorted
    inside each block by WT count, so the label column is printed once and the
    blocks are separated by a blank line.
    """
    counts = _read("s8in_counts.tsv")
    he = counts[counts["sample"].isin(HELA_UNITS)]
    rows: list[dict] = []
    for tool, g in he.groupby("tool", sort=False):
        tool = str(tool)
        if tool.startswith("ORCA_"):
            continue
        mt = str(g["mod_type"].iloc[0])
        w = g[g["sample"] == "HeLa_RNA004_WT"]
        v = g[g["sample"] == "HeLa_RNA004_IVT"]
        rows.append({
            "tool": tool, "label": pretty_tool(tool),
            "family": family_of(tool, mt), "mod_type": mt,
            "block": _block_of_a(tool, mt),
            
            #: differential tools have no callset on the unmodified control, and
            #: the counts table has no row for them), so it stays NaN and the
            #: bar panel draws nothing -- it used to be forced to 0.0, which the
            #: panel then printed as its "measured zero" open circle.
            "WT": float(w["n_sites"].iloc[0]) if len(w) else float("nan"),
            "IVT": float(v["n_sites"].iloc[0]) if len(v) else float("nan"),
        })
    df = pd.DataFrame(rows)
    # two entries can share a printed label (the inosine+m6A model is run once on
    # the m6A channel and once on the inosine channel): disambiguate the copies
    dup = df["label"].duplicated(keep=False)
    df.loc[dup, "label"] = [f"{l} [{m}]" for l, m in
                            zip(df.loc[dup, "label"], df.loc[dup, "mod_type"])]
    df["block"] = pd.Categorical(df["block"], BLOCKS_A, ordered=True)
    df = df.sort_values(["block", "WT"], ascending=[True, False], kind="stable")
    df["log10WT"] = np.log10(df["WT"] + 1.0)
    df["log10IVT"] = np.log10(df["IVT"] + 1.0)
    return df.reset_index(drop=True)


# --------------------------------------------------------------------------- #
# B -- the ORCA channels on the same HeLa RNA004 library
# --------------------------------------------------------------------------- #
def panel_orca() -> pd.DataFrame:
    """One row per ORCA channel: the calls on HeLa RNA004 WT and IVT.

    ORCA reports several modification types from one run, so its eight channels
    are a facet of their own instead of rows of the m6A count matrix.
    """
    df = _read("figS8_orca_counts.tsv").copy()
    df["tool"] = df["tool"].astype(str)
    df["label"] = [pretty_tool(t) for t in df["tool"]]
    df["family"] = "tool"
    df["block"] = BLOCK_ORCA
    df["WT"] = df["WT"].astype(float)
    df["IVT"] = df["IVT"].astype(float)
    df["log10WT"] = np.log10(df["WT"] + 1.0)
    df["log10IVT"] = np.log10(df["IVT"] + 1.0)
    return df.sort_values("WT", ascending=False,
                          kind="stable").reset_index(drop=True)


# --------------------------------------------------------------------------- #
# C -- the unmodified RNA004 Curlcake control
# --------------------------------------------------------------------------- #
def panel_curlcake() -> pd.DataFrame:
    """One row per entry measured on the unmodified RNA004 Curlcake control.

    ``calls`` is the number of sites the entry reports on the synthetic control
    (a false positive by construction), ``per1e6`` the same count per 10^6
    candidate adenosines and ``specificity`` = 1 - per1e6 / 10^6.  The Dorado
    models are read at ``DORADO_CURLAKE_PCT``; the two m6A tools were run once
    on this library and carry no threshold.
    """
    scan = _read("s8in_fpr_curlcake_scan.tsv")
    r4 = scan[(scan["sample"] == CURLAKE_RNA004_UNIT) & (scan["mod_type"] == "m6A")]
    r4 = r4[(r4["threshold_pct"].isna())
            | (r4["threshold_pct"] == DORADO_CURLAKE_PCT)]
    rows: list[dict] = []
    for _, r in r4.iterrows():
        fam = str(r["family"])
        block, colour = CURLAKE_FAMILY.get(fam, ("other m6A tools", "tool"))
        per1e6 = float(r["fp_per_1e6_candidates"])
        rows.append({
            "tool": str(r["tool"]), "label": pretty_tool(str(r["tool"])),
            "block": block, "family": colour, "run_family": fam,
            "calls": float(r["n_fp"]), "per_10kb": float(r["fp_per_10kb"]),
            "per1e6": per1e6, "specificity": 1.0 - per1e6 / 1e6,
        })
    df = pd.DataFrame(rows)
    df["block"] = pd.Categorical(df["block"], BLOCKS_C, ordered=True)
    df["log10calls"] = np.log10(df["calls"] + 1.0)
    df = df.sort_values(["block", "calls"], ascending=[True, False], kind="stable")
    return df.reset_index(drop=True)


# --------------------------------------------------------------------------- #
# D -- matching window against detection quality
# --------------------------------------------------------------------------- #
#: the two tools the legend quotes by name when reading the window sweep
HIGHLIGHT = {"m6Anet": "m6Anet", "ELIGOS2_solo": "ELIGOS2_solo"}

#: the best DRACH-specific Dorado model by PPV (numbers quoted in the legend)
BEST_DRACH = "Dorado_sup@v5.0.0_m6A_DRACH@v1"


def panel_window() -> pd.DataFrame:
    """Window (0-50 bp) against overlap hit rate and single-nucleotide accuracy.

    ``hit_rate`` is the share of calls that overlap a GLORI site within the
    window (it can only grow with the window), ``localization_accuracy`` is the
    share of those overlaps that are exact single-nucleotide matches -- the
    resolution/sensitivity trade-off the manuscript reports.
    """
    d = pd.read_csv(LOCALIZATION, sep="\t")
    d = d[(d["species"] == "Human") & (d["sample"] == "HeLa_RNA004_WT")]
    d = d[["tool", "window", "hit_rate", "localization_accuracy",
           "n_calls_in_universe"]].copy()
    d["label"] = [pretty_tool(t) for t in d["tool"]]
    return d.sort_values(["tool", "window"]).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# E -- Pearson r and effect size per ratio tool
# --------------------------------------------------------------------------- #
def _fisher_ci(r: float, n: int) -> tuple[float, float]:
    """95 % CI of a Pearson r via the Fisher z transform."""
    if not np.isfinite(r) or n < 4:
        return (float("nan"), float("nan"))
    z = np.arctanh(min(max(r, -0.999999), 0.999999))
    se = 1.0 / np.sqrt(n - 3.0)
    return (float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se)))


#: bootstrap settings of the panel E effect size (fixed seed, documented)
BOOT = 1000
BOOT_SEED = 20260921


def _slope_ci(pairs: pd.DataFrame, tool: str, rng: np.random.Generator,
              ) -> tuple[float, float, float]:
    """OLS slope of the tool ratio on the GLORI ratio, with a bootstrap CI.

    A slope of 1.0 means the tool tracks GLORI proportionally; > 1 means it
    systematically overestimates the modification ratio.
    """
    g = pairs[pairs["tool"] == tool][["ratio", "glori_ratio"]]
    g = g.apply(pd.to_numeric, errors="coerce").dropna()
    if len(g) < 30:
        return (np.nan, np.nan, np.nan)
    x = g["glori_ratio"].to_numpy(float)
    y = g["ratio"].to_numpy(float)
    slope = float(np.polyfit(x, y, 1)[0])
    draws = np.empty(BOOT)
    n = len(x)
    for i in range(BOOT):
        idx = rng.integers(0, n, n)
        draws[i] = np.polyfit(x[idx], y[idx], 1)[0]
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return (slope, float(lo), float(hi))


def panel_effect() -> pd.DataFrame:
    """One row per ratio tool: Pearson r and the effect size, both with 95 % CI.

    r and the matched-site count come from the frozen agreement table (Fisher-z
    CI computed here), the slope from the frozen matched-pair table.  Tools
    whose score is not a modification ratio (ELIGOS2, NanoSPA) have no row.
    """
    agree = _read("figS8_modratio_agreement.tsv").set_index("tool")
    
    #: panel G carries Lin's concordance correlation coefficient next to the
    #: calibration slope.  CCC and its interval come from the frozen delivered
    #: table `tables/S8F_glori_agreement.tsv` (the same numbers the reply letter
    #: and the manuscript quote), never recomputed here.
    ccc = pd.read_csv(TABLES / "S8F_glori_agreement.tsv", sep="\t").set_index("tool")
    pairs = pd.read_csv(SRC_TABLES / "figS8F_pairs.tsv", sep="\t")
    rng = np.random.default_rng(BOOT_SEED)
    rows: list[dict] = []
    for tool, a in agree.iterrows():
        tool = str(tool)
        if tool in NON_RATIO_TOOLS:
            continue
        n_m = int(a["n_matched"])
        rr = float(a["pearson_r"])
        lo, hi = _fisher_ci(rr, n_m)
        sl, sl_lo, sl_hi = _slope_ci(pairs, tool, rng)
        if not np.isfinite(sl):
            continue
        if tool not in ccc.index:
            raise SystemExit(f"no Lin's CCC for ratio tool {tool} in S8F_glori_agreement.tsv")
        c = ccc.loc[tool]
        rows.append({
            "tool": tool, "label": pretty_tool(tool),
            "family": family_of(tool, "m6A"),
            "r": rr, "r_lo": lo, "r_hi": hi,
            "slope": sl, "slope_lo": sl_lo, "slope_hi": sl_hi,
            "ccc": float(c["ccc"]), "ccc_lo": float(c["ccc_lo95"]),
            "ccc_hi": float(c["ccc_hi95"]),
            "n_matched": n_m,
        })
    df = pd.DataFrame(rows)
    return df.sort_values("r", ascending=False,
                          na_position="last").reset_index(drop=True)


# --------------------------------------------------------------------------- #
# anchors written to tables/S8_anchors.tsv (acceptance + reply-letter numbers)
# --------------------------------------------------------------------------- #
def _span(loc: pd.DataFrame, tool: str, col: str) -> str:
    """``first -> last`` value of one curve across the window sweep."""
    g = loc[loc["tool"] == tool].sort_values("window")
    return f"{g[col].iloc[0]:.4f} -> {g[col].iloc[-1]:.4f}"


def anchors() -> list[tuple[str, str, str]]:
    """(panel, quantity, value) rows for ``tables/S8_anchors.tsv``."""
    a = panel_a()
    a_idx = a.set_index("label")
    b = panel_orca()
    c = panel_curlcake()
    e = panel_effect()
    e_idx = e.set_index("label")
    loc = panel_window()

    dr = a[a["label"].str.contains("DRACH")]
    cur_dr = c[c["block"] == "m6A DRACH"]
    cur_nd = c[c["block"] == "m6A (non-DRACH)"]
    c_idx = c.set_index("label")
    othermod = int((a["block"] == "other modification models").sum())
    return [
        ("A", "rows on the HeLa count matrix (no ORCA)",
         f"{len(a)} (Dorado {int((a['block'] == 'Dorado m6A models').sum())}, "
         f"other m6A tools {int((a['block'] == 'other m6A tools').sum())}, "
         f"other modification models {othermod})"),
        ("A", "m6Anet detected sites, HeLa RNA004 WT",
         f"{int(a_idx.loc['m6Anet', 'WT'])}"),
        ("A", "DRACH models detected sites, HeLa RNA004 IVT",
         f"{int(dr['IVT'].min())}-{int(dr['IVT'].max())}"),
        ("B", "ORCA channels drawn (WT calls)",
         f"{len(b)} channels, {int(b['WT'].min())}-{int(b['WT'].max())} calls"),
        ("C", "Curcake false-positive calls, DRACH models (50 %)",
         f"{int(cur_dr['calls'].min())}-{int(cur_dr['calls'].max())}"),
        ("C", "Curcake false-positive calls, non-DRACH m6A models (50 %)",
         f"{int(cur_nd['calls'].min())}-{int(cur_nd['calls'].max())}"),
        ("C", "Curcake false positives per 10 kb, m6Anet",
         f"{float(c_idx.loc['m6Anet', 'per_10kb']):.1f}"),
        ("C", "Curcake false positives per 10 kb, NanoSPA_m6A",
         f"{float(c_idx.loc['NanoSPA_m6A', 'per_10kb']):.1f}"),
        ("C", "entries on the Curlcake control", f"{len(c)}"),
        ("C", "specificity range on the Curlcake control",
         f"{c['specificity'].min() * 100:.2f}-"
         f"{c['specificity'].max() * 100:.2f} %"),
        ("D", "ELIGOS2_solo hit rate 0 -> 50 bp",
         _span(loc, "ELIGOS2_solo", "hit_rate")),
        ("D", "ELIGOS2_solo localization accuracy 0 -> 50 bp",
         _span(loc, "ELIGOS2_solo", "localization_accuracy")),
        ("D", "m6Anet hit rate / accuracy 0 -> 50 bp",
         f"{_span(loc, 'm6Anet', 'hit_rate')} / "
         f"{_span(loc, 'm6Anet', 'localization_accuracy')}"),
        ("D", "best DRACH hit rate / accuracy 0 -> 50 bp",
         f"{_span(loc, BEST_DRACH, 'hit_rate')} / "
         f"{_span(loc, BEST_DRACH, 'localization_accuracy')}"),
        ("E", "m6Anet r vs GLORI (95 % CI)",
         f"{e_idx.loc['m6Anet', 'r']:.3f} ({e_idx.loc['m6Anet', 'r_lo']:.3f}-"
         f"{e_idx.loc['m6Anet', 'r_hi']:.3f}), n = "
         f"{int(e_idx.loc['m6Anet', 'n_matched'])}"),
        ("E", "m6Anet effect size (OLS slope of tool on GLORI ratio)",
         f"{e_idx.loc['m6Anet', 'slope']:.3f} "
         f"({e_idx.loc['m6Anet', 'slope_lo']:.3f}-"
         f"{e_idx.loc['m6Anet', 'slope_hi']:.3f})"),
        #: 2026-09-26: panel G now carries the calibration slope and Lin's CCC
        ("G", "m6Anet Lin's CCC (95 % CI)",
         f"{e_idx.loc['m6Anet', 'ccc']:.3f} ({e_idx.loc['m6Anet', 'ccc_lo']:.3f}-"
         f"{e_idx.loc['m6Anet', 'ccc_hi']:.3f})"),
        ("G", "ratio tools, Lin's CCC range",
         f"{len(e)} ({e['ccc'].min():.3f}-{e['ccc'].max():.3f})"),
        ("E", "ratio tools on the r / slope blocks",
         f"{int(e['r'].notna().sum())}"),
        ("F", "ratio tools with an effect size (slope, bootstrap CI)",
         f"{len(e)} (slopes {e['slope'].min():.3f}-{e['slope'].max():.3f})"),
    ]
