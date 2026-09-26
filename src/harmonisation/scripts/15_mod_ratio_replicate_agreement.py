#!/usr/bin/env python
"""
Revision analysis NA-6 (redo): modification-ratio agreement with GLORI,
accounting for BIOLOGICAL REPLICATES.

Motivation
----------
The original manuscript reported a single Pearson r between each tool's
predicted modification ratio and the GLORI measured ratio (Fig. S3C), computed
on collapsed / replicate-3-only call sets, and interpreted it as "quantitative
accuracy".  Reviewers flagged two issues:
  * R3-2 / E6 : biological-replicate structure was not incorporated.
  * R3-5 / R1-9 / E7 : Pearson r measures ASSOCIATION, not AGREEMENT.

This script re-derives the analysis from the reconstructed per-replicate
call sets in harmonisation/callsets, and for every (species x group x tool x
replicate) reports association (Pearson r, Spearman rho, both with bootstrap
95% CI) AND agreement (Lin's CCC with CI, Bland-Altman bias + limits of
agreement, MAE, RMSE).  It also reproduces the OLD "pooled" (replicate-collapsed)
single number so the two views can be contrasted directly.

Outputs -> analysis/mod_ratio_replicates/{tables,figures}
"""

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
import os
import glob
import json
import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
BENCH = str(_RB)
SITES = f"{_RB}/data"
# Read from the curated `callsets` layer (ref_base hard-filtered to be
# compatible with mod-type + strand; rebuild already baked in) rather than the
# raw `callsets`.  Layout is identical: <species>/<group>/<mod>/<tool>/<sample>.tsv
CALLSETS = f"{SITES}/callsets/RNA002"
REGISTRY = f"{SITES}/manifest/sample_registry.csv"
OUT = f"{_XB}/mod_ratio_replicates"
TABDIR = f"{OUT}/tables"
FIGDIR = f"{OUT}/figures"
os.makedirs(TABDIR, exist_ok=True)
os.makedirs(FIGDIR, exist_ok=True)

GLORI = {
    "Arabidopsis": f"{_XB}/third_party/NGS/GLORI/Arabidopsis_GLORI.bed",
    "Mouse":       f"{_XB}/third_party/NGS/GLORI/Mouse_GLORI_liftover.bed",
    "Human":       f"{_XB}/third_party/NGS/GLORI/Hela_GLORI.bed",
}
# m6A tools that report a genuine per-site modification ratio (0-1), comparable
# to the GLORI ratio.  Restricted to the FOUR tools used in the original
# manuscript mod-ratio analysis (Total_new/mod_ratio_correlation_change.ipynb):
# m6Anet, MINES, DENA, Nanom6A.  Others report probability / p-value / FDR /
# error, not a stoichiometry, so they are excluded from the ratio analysis.
TOOLS = ["m6Anet", "MINES", "DENA", "Nanom6A"]
TOOL_LABEL = {
    "m6Anet": "m6Anet", "MINES": "MINES", "Nanom6A": "Nanom6A", "DENA": "DENA",
}
MIN_OVERLAP = 30          # minimum GLORI-matched sites to report a metric
BOOT = 1000               # bootstrap replicates
RNG = np.random.default_rng(20260917)

# Figure style: match the manuscript's house style — all-Arial, large fonts,
# a closed thick black box (all 4 spines, linewidth 2), and NO in-plot gridlines.
rcParams.update({
    "font.family": "Arial",
    "font.size": 16,
    "axes.titlesize": 20,
    "axes.labelsize": 18,
    "axes.linewidth": 2.0,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 13,
    "axes.grid": False,
    "figure.dpi": 120,
    "savefig.dpi": 300,
})
SP_COLORS = {"Arabidopsis": "#1b7837", "Mouse": "#7570b8", "Human": "#d95f02"}
REP_COLORS = ["#2166ac", "#b2182b", "#4d9221"]   # rep1 / rep2 / rep3
# global replicate-tag -> colour so the same replicate is the same colour in
# every panel (single-replicate tools must not silently inherit "rep1" colour).
TAG_COLORS = {
    "rep1": "#2166ac", "rep2": "#b2182b", "rep3": "#4d9221",
    "studyA": "#7570b8", "studyB": "#e6a817",
}

# --------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------
def lin_ccc(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    vx, vy = x.var(ddof=0), y.var(ddof=0)
    cova = np.mean((x - x.mean()) * (y - y.mean()))
    denom = vx + vy + (x.mean() - y.mean()) ** 2
    return (2 * cova / denom) if denom > 0 else np.nan


def boot_ci(fn, x, y, B=BOOT):
    n = len(x)
    if n < 5:
        return (np.nan, np.nan)
    idx = np.arange(n)
    vals = np.empty(B)
    for b in range(B):
        s = RNG.choice(idx, size=n, replace=True)
        try:
            vals[b] = fn(x[s], y[s])
        except Exception:
            vals[b] = np.nan
    vals = vals[~np.isnan(vals)]
    if len(vals) < 10:
        return (np.nan, np.nan)
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def eval_pair(glori, tool):
    """Return metric dict comparing one replicate's tool vs GLORI ratios."""
    glor = np.asarray(glori, float)
    tol = np.asarray(tool, float)
    n = len(tol)
    if n < MIN_OVERLAP:
        return None
    r, p_r = stats.pearsonr(glor, tol)
    rho, p_rho = stats.spearmanr(glor, tol)
    ccc = lin_ccc(glor, tol)
    diff = tol - glor
    bias = diff.mean()
    sd = diff.std(ddof=1)
    mae = np.abs(diff).mean()
    rmse = np.sqrt((diff ** 2).mean())
    r_lo, r_hi = boot_ci(lambda a, b: stats.pearsonr(a, b)[0], glor, tol)
    c_lo, c_hi = boot_ci(lin_ccc, glor, tol)
    return dict(
        n_overlap=n, pearson_r=r, pearson_p=p_r, pearson_lo=r_lo, pearson_hi=r_hi,
        spearman_rho=rho, spearman_p=p_rho, ccc=ccc, ccc_lo=c_lo, ccc_hi=c_hi,
        ba_bias=bias, ba_sd=sd, ba_loa_low=bias - 1.96 * sd, ba_loa_hi=bias + 1.96 * sd,
        mae=mae, rmse=rmse,
        mean_glori=glor.mean(), mean_tool=tol.mean(),
    )


# --------------------------------------------------------------------------
# Data loading
# --------------------------------------------------------------------------
def load_glori(path):
    d = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 3],
                    names=["chrom", "pos_raw", "glori_ratio"])
    d["chrom"] = "chr" + d["chrom"].astype(str).str.replace("^chr", "", regex=True)
    d["pos_raw"] = pd.to_numeric(d["pos_raw"], errors="coerce")
    d = d.dropna(subset=["pos_raw"])
    d["pos_raw"] = d["pos_raw"].astype(int)
    d = d.groupby(["chrom", "pos_raw"], as_index=False)["glori_ratio"].mean()
    return d


def extract_ratio(d):
    """Return a per-site modification ratio (0-1) Series aligned to d.index.

    The reconstructed call sets store the ratio in different columns per tool:
      * m6Anet / Nanom6A  -> dedicated `mod_ratio` column (score is a probability)
      * DENA / MINES      -> `score` column, score_type in {m6a_ratio, mod_ratio}
                             (the `mod_ratio` column is empty in the current build)
    Rule: prefer a populated `mod_ratio`; otherwise fall back to `score` when its
    score_type denotes a ratio.  Values are clipped to [0, 1].
    """
    mr = (pd.to_numeric(d["mod_ratio"], errors="coerce")
          if "mod_ratio" in d.columns else pd.Series(np.nan, index=d.index))
    if mr.notna().sum() >= MIN_OVERLAP:
        ratio = mr
    else:
        st = " ".join(map(str, d["score_type"].dropna().unique())).lower()
        if "ratio" in st:
            ratio = pd.to_numeric(d["score"], errors="coerce")
        else:
            ratio = mr  # all-NaN -> caller filters to empty
    return ratio


def load_registry():
    reg = pd.read_csv(REGISTRY, sep="\t")
    m = {}
    for _, row in reg.iterrows():
        m[row["sample"]] = dict(
            replicate_tag=row.get("replicate_tag", ""),
            sequencing_unit=row.get("sequencing_unit", row["sample"]),
            independence_class=row.get("independence_class", "independent"),
            dataset_group=row.get("dataset_group", ""),
        )
    return m


def iter_group_tool_files():
    """Yield (species, group, tool, replicate_sample, tsv_path)."""
    for species in ["Arabidopsis", "Mouse", "Human"]:
        for gdir in sorted(glob.glob(f"{CALLSETS}/{species}/*/m6A")):
            # gdir = .../<species>/<group>/m6A  -> group is the parent of "m6A"
            group = os.path.basename(os.path.dirname(gdir))
            for tool in TOOLS:
                for f in sorted(glob.glob(f"{gdir}/{tool}/*.tsv")):
                    sample = os.path.basename(f)[:-4]
                    yield species, group, tool, sample, f


# --------------------------------------------------------------------------
# Build tidy per-site overlap table + per-replicate metrics
# --------------------------------------------------------------------------
def main():
    gloris = {sp: load_glori(p) for sp, p in GLORI.items()}
    reg = load_registry()

    site_rows = []      # tidy per matched site (for scatter / BA / bins)
    rep_rows = []       # per-replicate metric summaries
    skipped = []

    for species, group, tool, sample, f in iter_group_tool_files():
        try:
            d = pd.read_csv(f, sep="\t", low_memory=False)
        except pd.errors.EmptyDataError:
            continue
        # callsets uses BED `start` (0-based) as the position column
        if "pos_raw" not in d.columns and "start" in d.columns:
            d["pos_raw"] = pd.to_numeric(d["start"], errors="coerce")
        if "mod_ratio" not in d.columns and "score" not in d.columns:
            continue
        d["tool_ratio"] = extract_ratio(d)
        d = d.dropna(subset=["tool_ratio", "chrom", "pos_raw"])
        d = d[(d["tool_ratio"] >= 0) & (d["tool_ratio"] <= 1)]
        if len(d) == 0:
            skipped.append((species, group, tool, sample, "no_mod_ratio"))
            continue
        # reference-consistency gate: m6A sites must sit on an A/T reference base.
        # A callset whose positions are off the current reference (wrong genome
        # build / bad run) drops to ~random A/T.  Flag & exclude these so the
        # agreement numbers are never computed on mis-mapped coordinates.
        if "ref_base" in d.columns:
            rb = d["ref_base"].astype(str).str.upper()
            known = rb.isin(list("ACGT"))
            if known.sum() >= 50:
                at = float((rb[known] == "A").sum() + (rb[known] == "T").sum()) / int(known.sum())
                if at < 0.90:
                    skipped.append((species, group, tool, sample,
                                    f"off_reference_AT={100*at:.0f}%"))
                    continue
        d["pos_raw"] = d["pos_raw"].astype(int)
        g = gloris[species]
        # collapse duplicate (chrom,pos) within a replicate (mean tool ratio)
        d = d.groupby(["chrom", "pos_raw"], as_index=False).agg(
            tool_ratio=("tool_ratio", "mean"),
        )
        # exact single-nucleotide genomic match (chrom + pos), w=0
        d["_key"] = d["chrom"] + ":" + d["pos_raw"].astype(str)
        gg = g.copy()
        gg["_key"] = gg["chrom"] + ":" + gg["pos_raw"].astype(str)
        gm = dict(zip(gg["_key"], gg["glori_ratio"]))
        d["glori_ratio"] = d["_key"].map(gm)
        m = d.dropna(subset=["glori_ratio"])

        rmeta = reg.get(sample, {})
        n_rep = rmeta.get("independence_class", "independent")
        for _, rr in m.iterrows():
            site_rows.append(dict(
                species=species, group=group, tool=tool, sample=sample,
                replicate_tag=rmeta.get("replicate_tag", ""),
                independence_class=n_rep,
                sequencing_unit=rmeta.get("sequencing_unit", sample),
                chrom=rr["chrom"], pos_raw=rr["pos_raw"],
                glori_ratio=rr["glori_ratio"], tool_ratio=rr["tool_ratio"],
            ))

        met = eval_pair(m["glori_ratio"].values, m["tool_ratio"].values)
        if met is None:
            skipped.append((species, group, tool, sample, f"overlap<{MIN_OVERLAP}({len(m)})"))
            continue
        rep_rows.append(dict(
            species=species, group=group, tool=tool, sample=sample,
            replicate_tag=rmeta.get("replicate_tag", ""),
            independence_class=n_rep,
            sequencing_unit=rmeta.get("sequencing_unit", sample),
            n_calls=len(d), **met,
        ))

    sites = pd.DataFrame(site_rows)
    per_rep = pd.DataFrame(rep_rows)
    sites.to_csv(f"{TABDIR}/mod_ratio_matched_sites.tsv", sep="\t", index=False)
    per_rep.to_csv(f"{TABDIR}/mod_ratio_per_replicate.tsv", sep="\t", index=False)
    with open(f"{TABDIR}/skipped_pairs.json", "w") as fh:
        json.dump([dict(species=a, group=b, tool=c, sample=d, reason=e)
                   for a, b, c, d, e in skipped], fh, indent=2)

    print(f"[sites] {len(sites)} matched sites; {len(per_rep)} per-replicate rows; "
          f"{len(skipped)} skipped (sample,tool) pairs")
    return sites, per_rep


# --------------------------------------------------------------------------
# Aggregation: per-replicate mean +/- SD vs pooled (legacy) contrast
# --------------------------------------------------------------------------
def aggregate(sites, per_rep):
    # independence: how many independent sequencing units per group
    units = (per_rep.groupby(["species", "group"])["sequencing_unit"]
             .nunique().rename("n_sequencing_units"))
    # how many are genuinely INDEPENDENT biological/sequencing replicates
    # (cross_study / nested_subset / same_run_split do NOT support inference)
    indep = (per_rep[per_rep.independence_class == "independent"]
             .groupby(["species", "group"])["sequencing_unit"].nunique()
             .rename("n_independent_units"))

    # per-replicate summary (mean +/- SD of each metric across replicates)
    summ = []
    for (sp, grp, tool), sub in per_rep.groupby(["species", "group", "tool"]):
        row = dict(species=sp, group=grp, tool=tool,
                   n_replicates=len(sub),
                   replicate_labels=",".join(sorted(sub["replicate_tag"].astype(str))))
        for metric in ["pearson_r", "ccc", "ba_bias", "mae", "rmse", "spearman_rho"]:
            vals = sub[metric].dropna().values
            row[f"{metric}_mean"] = vals.mean() if len(vals) else np.nan
            row[f"{metric}_sd"] = vals.std(ddof=1) if len(vals) > 1 else np.nan
        row["n_overlap_mean"] = sub["n_overlap"].mean()
        row["n_overlap_min"] = sub["n_overlap"].min()
        row["n_overlap_max"] = sub["n_overlap"].max()
        # POOLED / LEGACY: concatenate all replicate sites, ignore replicate id
        pool = sites[(sites.species == sp) & (sites.group == grp) & (sites.tool == tool)]
        pm = eval_pair(pool["glori_ratio"].values, pool["tool_ratio"].values)
        if pm:
            row["pooled_pearson_r"] = pm["pearson_r"]
            row["pooled_ccc"] = pm["ccc"]
            row["pooled_bias"] = pm["ba_bias"]
            row["pooled_mae"] = pm["mae"]
            row["pooled_n"] = pm["n_overlap"]
        summ.append(row)
    summary = pd.DataFrame(summ)
    if len(units):
        summary = summary.merge(units.reset_index(), on=["species", "group"], how="left")
    summary["n_independent_units"] = (
        summary.merge(indep.reset_index(), on=["species", "group"], how="left")
        ["n_independent_units"].fillna(0).astype(int).values)
    # is inference allowed (>=2 independent biological/sequencing replicates)?
    # cross-study / nested / split units count as descriptive only.
    summary["replicate_inference"] = np.where(
        summary["n_independent_units"] >= 2, "yes", "descriptive_only")
    summary.to_csv(f"{TABDIR}/mod_ratio_summary_by_group.tsv", sep="\t", index=False)

    # headline table: pooled vs per-replicate r / CCC  (legacy-vs-new contrast)
    head = summary[[
        "species", "group", "tool", "n_replicates", "n_sequencing_units",
        "n_independent_units", "replicate_inference", "replicate_labels",
        "pooled_pearson_r", "pearson_r_mean", "pearson_r_sd",
        "pooled_ccc", "ccc_mean", "ccc_sd",
        "pooled_bias", "ba_bias_mean", "ba_bias_sd", "mae_mean",
        "n_overlap_mean",
    ]].copy()
    head.to_csv(f"{TABDIR}/pooled_vs_per_replicate.tsv", sep="\t", index=False)
    return summary


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------
def strip_spines(ax):
    # Manuscript house style: keep the full closed box, bold black borders, no grid.
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(2.0)
        s.set_color("black")
    ax.grid(False)
    ax.set_axisbelow(True)


def fig_scatter_per_replicate(sites, per_rep):
    """One row per species; columns = tools with data. Colour = replicate.
    Shows replicate-to-replicate spread of the ratio relationship vs GLORI."""
    # pick the 3 species with most coverage
    species = ["Arabidopsis", "Human", "Mouse"]
    tools_by_sp = {}
    for sp in species:
        ts = sorted(per_rep[per_rep.species == sp]["tool"].unique())
        tools_by_sp[sp] = ts
    ncols = max(len(t) for t in tools_by_sp.values())
    fig, axes = plt.subplots(len(species), ncols,
                             figsize=(4.3 * ncols, 4.5 * len(species)))
    axes = np.atleast_2d(axes)
    present_tags = []
    for ri, sp in enumerate(species):
        # one representative group per species with replicates (WT preferred)
        grp = "Arabidopsis_WT" if sp == "Arabidopsis" else (
            "HeLa_WT" if sp == "Human" else "Mouse_WT")
        for ci in range(ncols):
            ax = axes[ri, ci]
            strip_spines(ax)
            ax.set_xlim(0, 1); ax.set_ylim(0, 1)
            ax.set_box_aspect(1)
            ax.plot([0, 1], [0, 1], ls="--", lw=1, color="grey")
            if ci >= len(tools_by_sp[sp]):
                ax.axis("off"); continue
            tool = tools_by_sp[sp][ci]
            sub = sites[(sites.species == sp) & (sites.group == grp) & (sites.tool == tool)]
            for rep in sorted(sub["replicate_tag"].unique()):
                rs = sub[sub.replicate_tag == rep]
                if len(rs) == 0:
                    continue
                col = TAG_COLORS.get(str(rep), "#666666")
                if str(rep) not in present_tags:
                    present_tags.append(str(rep))
                ax.scatter(rs["glori_ratio"], rs["tool_ratio"], s=8, alpha=0.35,
                           linewidths=0, color=col, label=str(rep))
            pr = per_rep[(per_rep.species == sp) & (per_rep.group == grp) & (per_rep.tool == tool)]
            if len(pr):
                ttl = f"{TOOL_LABEL[tool]}  (r={pr['pearson_r'].mean():.2f}, CCC={pr['ccc'].mean():.2f})"
            else:
                ttl = TOOL_LABEL[tool]
            ax.set_title(f"{grp}\n{ttl}", fontsize=11)
            if ci == 0:
                ax.set_ylabel("Tool predicted ratio")
            if ri == len(species) - 1:
                ax.set_xlabel("GLORI measured ratio")
    # combined legend over all replicate tags actually present
    handles = [plt.Line2D([], [], marker="o", ls="", color=TAG_COLORS[t],
                          markersize=8, label=t) for t in present_tags]
    fig.legend(handles=handles, loc="lower center", ncol=len(present_tags),
               frameon=False, bbox_to_anchor=(0.5, -0.012))
    fig.suptitle("Modification ratio vs GLORI — per biological replicate", y=1.0)
    fig.tight_layout(rect=[0, 0.025, 1, 0.99])
    fig.savefig(f"{FIGDIR}/mod_ratio_scatter_per_replicate.png", bbox_inches="tight")
    fig.savefig(f"{FIGDIR}/mod_ratio_scatter_per_replicate.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_forest(summary):
    """Forest / dot plot of Pearson r and Lin's CCC per (group,tool),
    showing per-replicate points + mean +/- SD. Highlights that inference is
    only valid where >=2 independent replicates exist."""
    df = summary.copy()
    # group label: "Species / WT / Tool"  (strip the leading "<Species>_" from group)
    df["lab"] = (df["species"] + " / "
                 + df["group"].str.replace(r"^[A-Za-z.]+_", "", regex=True)
                 + " / " + df["tool"].map(TOOL_LABEL))
    # order: species then group then tool
    df = df.sort_values(["species", "group", "tool"]).reset_index(drop=True)

    per_rep = pd.read_csv(f"{TABDIR}/mod_ratio_per_replicate.tsv", sep="\t")
    per_rep["lab"] = (per_rep["species"] + " / "
                      + per_rep["group"].str.replace(r"^[A-Za-z.]+_", "", regex=True)
                      + " / " + per_rep["tool"].map(TOOL_LABEL))

    n = len(df)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 0.34 * n + 2.6))
    for ax, metric, pooled_col, lo_col, hi_col, title, xlabel in [
        (axes[0], "pearson_r", "pooled_pearson_r", "pearson_lo", "pearson_hi",
         "Association  (Pearson r)", "Pearson r"),
        (axes[1], "ccc", "pooled_ccc", None, None,
         "Agreement  (Lin's CCC)", "Lin's concordance correlation"),
    ]:
        strip_spines(ax)
        ymap = {lab: i for i, lab in enumerate(df["lab"])}
        # per-replicate raw points
        for lab, i in ymap.items():
            pr = per_rep[per_rep.lab == lab]
            ax.scatter(pr[metric].dropna().values,
                       [i] * len(pr[metric].dropna()), s=26, color="#999999",
                       zorder=2, linewidths=0)
        # mean +/- SD
        for _, row in df.iterrows():
            i = ymap[row["lab"]]
            mean = row[f"{metric}_mean"]; sd = row[f"{metric}_sd"]
            if np.isnan(mean):
                continue
            infer = row["replicate_inference"] == "yes"
            col = SP_COLORS[row["species"]]
            mk = "o" if infer else "s"
            if not np.isnan(sd):
                ax.plot([mean - sd, mean + sd], [i, i], color=col, lw=2.2,
                        alpha=0.75, zorder=3, solid_capstyle="round")
            ax.scatter([mean], [i], marker=mk, s=70, color=col, zorder=4,
                       edgecolor="black", linewidth=0.5)
            p = row[pooled_col]
            if not np.isnan(p):
                ax.scatter([p], [i], marker="x", s=60, color="#d95f02", zorder=5,
                           linewidths=1.6)
        ax.set_yticks(range(n)); ax.set_yticklabels(df["lab"], fontsize=10)
        ax.set_xlabel(xlabel)
        ax.set_title(title)
        ax.axvline(0 if metric == "ccc" else 0, color="#cccccc", lw=0.8)
        ax.set_ylim(-1, n)
        ax.set_xlim(-0.6, 1.05)
    # shared legend
    from matplotlib.lines import Line2D
    leg = [
        Line2D([], [], marker="o", ls="", color="#333333", markersize=8,
               label="Mean, replicate-inference group"),
        Line2D([], [], marker="s", ls="", color="#333333", markersize=8,
               label="Mean, single/cross-study (descriptive)"),
        Line2D([], [], marker="x", ls="", color="#d95f02", markersize=9,
               label="Pooled (original, replicate-collapsed)"),
        Line2D([], [], ls="-", color="#333333", lw=2.2, label="+/- SD across replicates"),
        Line2D([], [], marker="o", ls="", color="#999999", markersize=6,
               label="Individual replicates"),
    ]
    axes[1].legend(handles=leg, loc="lower right", frameon=False, fontsize=10)
    fig.suptitle("GLORI-vs-tool modification-ratio: association and agreement, "
                 "per biological replicate", y=1.0)
    fig.tight_layout(rect=[0, 0, 0.86, 0.97])
    fig.savefig(f"{FIGDIR}/mod_ratio_forest_association_agreement.png", bbox_inches="tight")
    fig.savefig(f"{FIGDIR}/mod_ratio_forest_association_agreement.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_bland_altman(sites, summary):
    """Bland-Altman (tool - GLORI vs mean) for the top tool per species WT group,
    one panel per replicate, with bias + limits of agreement."""
    cases = [("Arabidopsis", "Arabidopsis_WT"), ("Human", "HeLa_WT"), ("Mouse", "Mouse_WT")]
    # best tool by mean CCC per species
    best = {}
    for sp, grp in cases:
        sub = summary[(summary.species == sp) & (summary.group == grp)]
        if len(sub) == 0:
            continue
        best[(sp, grp)] = sub.sort_values("ccc_mean", ascending=False)["tool"].iloc[0]

    # grid: rows = cases, cols = replicates
    ncols = 3
    fig, axes = plt.subplots(len(cases), ncols,
                             figsize=(5.0 * ncols, 4.7 * len(cases)))
    axes = np.atleast_2d(axes)
    for ri, (sp, grp) in enumerate(cases):
        tool = best[(sp, grp)]
        sub = sites[(sites.species == sp) & (sites.group == grp) & (sites.tool == tool)]
        reps = sorted(sub["replicate_tag"].unique())[:ncols]
        for ci in range(ncols):
            ax = axes[ri, ci]
            strip_spines(ax)
            ax.set_xlim(0, 1)
            ax.set_box_aspect(1)
            if ci >= len(reps):
                ax.axis("off"); continue
            rep = reps[ci]
            rs = sub[sub.replicate_tag == rep]
            mean = (rs["glori_ratio"] + rs["tool_ratio"]) / 2
            diff = rs["tool_ratio"] - rs["glori_ratio"]
            if len(rs) < MIN_OVERLAP:
                ax.set_title(f"{sp} {TOOL_LABEL[tool]} {rep}\n(insufficient sites)", fontsize=11)
                ax.set_xlabel("Mean of two methods"); ax.set_ylabel("Tool − GLORI")
                continue
            bias = diff.mean(); sd = diff.std(ddof=1)
            lo, hi = bias - 1.96 * sd, bias + 1.96 * sd
            ax.scatter(mean, diff, s=9, alpha=0.35, color=SP_COLORS[sp], linewidths=0)
            ax.axhline(bias, color="#333333", lw=1.6)
            ax.axhline(lo, color="#b2182b", ls="--", lw=1.3)
            ax.axhline(hi, color="#b2182b", ls="--", lw=1.3)
            # bias / LoA values live in the title: no small annotation text in
            # the panel (project rule 2026-09-19)
            ax.set_title(f"{grp}: {TOOL_LABEL[tool]}  {rep}  (n={len(rs)})\n"
                         f"bias = {bias:+.2f}  ·  95% LoA = {lo:+.2f} … {hi:+.2f}",
                         fontsize=11)
            ax.set_xlabel("Mean of two methods"); ax.set_ylabel("Tool − GLORI ratio")
    fig.suptitle("Bland–Altman agreement, per biological replicate "
                 "(best tool per species by Lin's CCC)", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(f"{FIGDIR}/bland_altman_per_replicate.png", bbox_inches="tight")
    fig.savefig(f"{FIGDIR}/bland_altman_per_replicate.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_error_bins(sites, summary):
    """MAE vs GLORI ratio bins — shows tools work better at high modification
    ratio; per-replicate curves so replicate structure is visible."""
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0))
    cases = [("Arabidopsis", "Arabidopsis_WT"), ("Human", "HeLa_WT"), ("Mouse", "Mouse_WT")]
    bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    blab = ["0.0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"]
    for ax, (sp, grp) in zip(axes, cases):
        strip_spines(ax)
        ax.set_box_aspect(1)
        sub = sites[(sites.species == sp) & (sites.group == grp)]
        tools = sorted(sub["tool"].unique())
        sub = sub.copy()
        sub["bin"] = pd.cut(sub["glori_ratio"], bins=bins, labels=blab, include_lowest=True)
        for tool in tools:
            t = sub[sub.tool == tool]
            means, xs = [], []
            for j, bl in enumerate(blab):
                bt = t[t.bin == bl]
                if len(bt) >= MIN_OVERLAP:
                    xs.append(j); means.append((bt["tool_ratio"] - bt["glori_ratio"]).abs().mean())
            if xs:
                ax.plot(xs, means, marker="o", ms=6, lw=1.8, label=TOOL_LABEL[tool])
        ax.set_xticks(range(len(blab))); ax.set_xticklabels(blab, rotation=30, fontsize=10)
        ax.set_xlabel("GLORI measured ratio (bin)")
        ax.set_ylabel("Mean absolute error")
        ax.set_title(f"{sp} ({grp})")
    axes[0].legend(frameon=False, title="Tool", loc="upper right")
    axes[0].set_ylim(0, 1); axes[1].set_ylim(0, 1); axes[2].set_ylim(0, 1)
    fig.suptitle("Absolute ratio error vs modification level, per tool", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{FIGDIR}/mod_ratio_error_by_level.png", bbox_inches="tight")
    fig.savefig(f"{FIGDIR}/mod_ratio_error_by_level.pdf", bbox_inches="tight")
    plt.close(fig)


def write_readme(summary, per_rep, sites):
    """Concise markdown summary of headline numbers for the response letter."""
    md = ["# Modification-ratio agreement with GLORI — re-analysis with biological replicates\n",
          "Addresses R3-2/E6 (replicate structure), R3-5/R1-9/E7 (association vs agreement).\n",
          "## What changed vs the original Fig. S3C",
          "- Metrics now computed **per biological replicate** (each tool call set kept separate), "
          "not on a replicate-collapsed / rep3-only set.",
          "- Added **agreement** metrics (Lin's CCC, Bland–Altman bias + limits of agreement, MAE) "
          "alongside association (Pearson r, Spearman rho).",
          "- Pearson r and CCC carry **site-level bootstrap 95% CI**.",
          "- Inference restricted to groups with ≥ 2 independent sequencing units; "
          "mouse / E. coli / Curlcake are single-replicate or cross-study (descriptive only).\n",
          "## Headline agreement (WT groups, mean +/- SD across replicates)\n",
          "| Species | Group | Tool | n reps | Pearson r (pooled) | Pearson r (per-rep mean±SD) | CCC (pooled) | CCC (per-rep mean±SD) | Bias | MAE |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    wt = summary[summary.group.str.contains("_WT")]
    for _, row in wt.sort_values(["species", "tool"]).iterrows():
        md.append("| {} | {} | {} | {} | {:.3f} | {:.3f} ± {:.3f} | {:.3f} | {:.3f} ± {:.3f} | {:.3f} | {:.3f} |".format(
            row["species"], row["group"], TOOL_LABEL[row["tool"]],
            int(row["n_replicates"]),
            row.get("pooled_pearson_r", np.nan), row["pearson_r_mean"], row["pearson_r_sd"],
            row.get("pooled_ccc", np.nan), row["ccc_mean"], row["ccc_sd"],
            row["ba_bias_mean"], row["mae_mean"]))
    md.append("\n## Figures (in `figures/`)")
    md.append("- `mod_ratio_scatter_per_replicate` — tool-vs-GLORI ratio scatter, coloured by replicate.")
    md.append("- `mod_ratio_forest_association_agreement` — Pearson r and CCC per group/tool with replicate spread vs the pooled (original) value.")
    md.append("- `bland_altman_per_replicate` — agreement limits per replicate for the best tool per species.")
    md.append("- `mod_ratio_error_by_level` — absolute error across GLORI ratio bins.\n")
    md.append("## Tables (in `tables/`)")
    md.append("- `mod_ratio_matched_sites.tsv` — every GLORI-matched site (tidy, with replicate metadata).")
    md.append("- `mod_ratio_per_replicate.tsv` — per-replicate metrics.")
    md.append("- `mod_ratio_summary_by_group.tsv` — mean±SD across replicates + pooled contrast + inference flag.")
    md.append("- `pooled_vs_per_replicate.tsv` — headline pooled-vs-per-replicate table.")
    md.append("- `skipped_pairs.json` — (sample,tool) pairs without a usable ratio (e.g. MINES/DENA rep1-2 lack a `mod_ratio` column).")
    with open(f"{OUT}/README.md", "w") as fh:
        fh.write("\n".join(md) + "\n")


if __name__ == "__main__":
    import sys
    cache = "--figures-only" in sys.argv
    if cache:
        sites = pd.read_csv(f"{TABDIR}/mod_ratio_matched_sites.tsv", sep="\t")
        per_rep = pd.read_csv(f"{TABDIR}/mod_ratio_per_replicate.tsv", sep="\t")
        summary = pd.read_csv(f"{TABDIR}/mod_ratio_summary_by_group.tsv", sep="\t")
        print("[figures-only] loaded cached tables:", len(sites), "sites,",
              len(per_rep), "replicate rows,", len(summary), "summary rows")
    else:
        sites, per_rep = main()
        summary = aggregate(sites, per_rep)
    print("[summary] rows:", len(summary))
    fig_scatter_per_replicate(sites, per_rep)
    fig_forest(summary)
    fig_bland_altman(sites, summary)
    fig_error_bins(sites, summary)
    write_readme(summary, per_rep, sites)
    print("[done] outputs in", OUT)
