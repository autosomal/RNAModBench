#!/usr/bin/env python3
"""73 -- the published NanoNm call set beside our own re-analysis.

Why this figure exists
----------------------
Every non-m6A panel of this study characterises the tools on **unmodified**
controls only, because no m6A-centred reference can validate a tool that
predicts another modification (R2-2).  For NanoNm a published positive call set
does exist: the authors of the model deposited their own per-site Nm
stoichiometry for three cell lines (HEK293T, HeLa, HepG2) in the supplement of
their Molecular Cell paper.  This stage turns that table into the evidence the
figure needs, and measures how far our own re-analysis of the same model lands
from it.

What it writes
--------------
``sources/mmc6_fingerprint.tsv``      the source file's bytes/md5, its sheets, the
                                     parser that was used, and the md5 of every
                                     input the tables are derived from
``sources/ensembl112_gene_spans.tsv`` cached gene_name -> GRCh38 span, built from
                                     the Ensembl 112 GTF (gene lines only)
``tables/s11_published_sites.tsv``    one row per published site, with a
                                     ``parse_ok`` flag rather than silent drops
``tables/s11_coordinate_check.tsv``   per cell line: how many parsed sites fall
                                     inside the gene the sheet itself names
``tables/s11_overlap.tsv``            our HeLa WT / IVT calls against the published
                                     HeLa set, at the analysis ratio threshold
``tables/s11_ratio_agreement.tsv``    per-site agreement of the two mod_ratio
                                     values at the shared sites
``tables/s11_ratio_distribution.tsv`` the distributions panel C draws
``logs/73_tables_anchors.tsv``        every number the caption quotes, so the
                                     verifier can recompute rather than trust

Two limits are structural, not fixable by a better script, and the caption has to
carry both: the published table has **no unmodified control** (no IVT, no
Curlcake), so it cannot support any claim about false positives; and it was
produced from wild-type cells yet overlaps our *unmodified IVT* calls slightly
more than our wild-type calls, so it is a superset of both and cannot be read as
evidence that the shared sites are modified.

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    analysis/figS11_revision/scripts/73_figS11_tables.py
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
import hashlib
import re
import sys
import time
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PROJECT = _RB
OUT = (_RB / "analysis/figS11_revision")
TAB, SRC, LOG = (_RB / "analysis/figS11_revision/tables"), (_RB / "analysis/figS11_revision/sources"), (_RB / "analysis/figS11_revision/logs")

XLSX = (_XB / "third_party/NanoNm/mmc6.xlsx")
GTF = PROJECT.parent / "reference" / "GRCh38p14" / "ensembl112" / "Homo_sapiens.GRCh38.112.chr.gtf"
FAI = PROJECT.parent / "reference" / "GRCh38p13" / "GRCh38.primary_assembly.genome.fa.fai"
CALLSETS = (_RB / "data/callsets/RNA002/Human")

#: sheet name -> label used on the page and in the caption
CELL_LINES = {"Table.HEK293.ratio": "HEK293T", "HeLa": "HeLa", "HepG2": "HepG2"}
#: our own re-analysis: (folder, label, units)
CONDITIONS = (
    ("HeLa_WT", "HeLa WT", ("HeLa_WT1", "HeLa_WT2", "HeLa_WT3")),
    ("HeLa_IVT", "HeLa IVT", ("HeLa_IVT_rep1", "HeLa_IVT_rep2", "HeLa_IVT_rep3")),
)
#: the study's own analysis threshold for these call sets (10 % modified reads);
#: the delivered calls carry a lower bound already (support >= 20 plus the model
#: probability), this is the ratio cut the results section quotes
RATIO_MIN = 0.1
#: a site must land inside the gene its own sheet names, otherwise the published
#: coordinates are not in the assembly we analyse
MIN_IN_NAMED_GENE = 0.90
MAX_UNPARSED = 0.05

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
POS_RE = re.compile(r"^(chr[0-9XYM]+)_(\d+)$")


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


# --------------------------------------------------------------------------- #
# xlsx reading: openpyxl when the environment has it, stdlib otherwise
# --------------------------------------------------------------------------- #
def read_xlsx(path: Path) -> tuple[dict[str, pd.DataFrame], str]:
    """Return {sheet: DataFrame} plus the parser that was used.

    The repository's analysis environment has no openpyxl, and reading one
    supplementary table is not a reason to change an environment, so the stdlib
    path is the default: an xlsx is a zip of XML, and the two sheets patterns
    here (inline/shared strings plus numeric cells) are parsed directly.
    """
    try:
        import openpyxl                                        # noqa: F401
    except ImportError:
        pass
    else:
        book = pd.read_excel(path, sheet_name=None, engine="openpyxl")
        return {name: df for name, df in book.items()}, "openpyxl"
    return _read_xlsx_stdlib(path), "stdlib(zipfile+ElementTree)"


def _read_xlsx_stdlib(path: Path) -> dict[str, pd.DataFrame]:
    ns = NS
    out: dict[str, pd.DataFrame] = {}
    with zipfile.ZipFile(path) as z:
        strings: list[str] = []
        if "xl/sharedStrings.xml" in z.namelist():
            with z.open("xl/sharedStrings.xml") as fh:
                for _ev, el in ET.iterparse(fh, events=("end",)):
                    if el.tag == ns + "si":
                        strings.append("".join(t.text or "" for t in el.iter(ns + "t")))
                        el.clear()
        # sheet order follows the workbook relationships, not the file names
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        target = {r.get("Id"): r.get("Target") for r in rels}
        sheets = []
        for sh in wb.iter(ns + "sheet"):
            tgt = target[sh.get("{http://schemas.openxmlformats.org/officeDocument/"
                                "2006/relationships}id")]
            sheets.append((sh.get("name"), "xl/" + tgt.lstrip("/")))
        for name, member in sheets:
            rows: list[list[str | None]] = []
            with z.open(member) as fh:
                for _ev, el in ET.iterparse(fh, events=("end",)):
                    if el.tag != ns + "row":
                        continue
                    cells: list[str | None] = []
                    for c in el.findall(ns + "c"):
                        v = c.find(ns + "v")
                        txt = v.text if v is not None else None
                        if txt is not None and c.get("t") == "s":
                            txt = strings[int(txt)]
                        cells.append(txt)
                    rows.append(cells)
                    el.clear()
            width = max(len(r) for r in rows)
            out[name] = pd.DataFrame([r + [None] * (width - len(r)) for r in rows[1:]],
                                     columns=[str(i) for i in range(width)])
    return out


# --------------------------------------------------------------------------- #
# gene spans (cached; the GTF is 1.5 GB so the cache is keyed on its mtime)
# --------------------------------------------------------------------------- #
def gene_spans() -> pd.DataFrame:
    cache = SRC / "ensembl112_gene_spans.tsv"
    stamp = SRC / "ensembl112_gene_spans.mtime"
    if cache.exists() and stamp.exists() and stamp.read_text().strip() == str(GTF.stat().st_mtime_ns):
        return pd.read_csv(cache, sep="\t")
    rows: list[tuple[str, int, int, str]] = []
    with GTF.open() as fh:
        for line in fh:
            if "\tgene\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            m = re.search(r'gene_name "([^"]+)"', f[8])
            if m:
                rows.append((f[0], int(f[3]), int(f[4]), m.group(1)))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "gene"])
    df.sort_values(["gene", "chrom", "start"]).to_csv(cache, sep="\t", index=False)
    stamp.write_text(str(GTF.stat().st_mtime_ns))
    return df


# --------------------------------------------------------------------------- #
def main() -> None:
    t0 = time.time()
    for d in (TAB, SRC, LOG):
        d.mkdir(parents=True, exist_ok=True)
    if not XLSX.exists():
        raise SystemExit(f"source table missing: {XLSX}")

    # ---- our own re-analysis, one row per site (max ratio over the units) ----
    ours: dict[str, pd.DataFrame] = {}
    prov: list[tuple[str, str, str, str]] = []
    for folder, label, units in CONDITIONS:
        parts = []
        for u in units:
            path = CALLSETS / folder / "Nm" / "NanoNm" / f"{u}.tsv"
            d = pd.read_csv(path, sep="\t", usecols=["chrom", "start", "score"])
            d["key"] = d.chrom + ":" + d.start.astype("int64").astype(str)
            d = d[d.score >= RATIO_MIN]
            parts.append(d)
            prov.append((f"our callset {label}", str(path), md5(path),
                         f"{len(d)} rows at ratio >= {RATIO_MIN}"))
        merged = pd.concat(parts).groupby("key", as_index=False).score.max()
        ours[label] = merged.rename(columns={"score": "our_ratio"})
        print(f"[{label}] unit-level rows {sum(len(p) for p in parts)} -> "
              f"{len(merged)} unique sites at ratio >= {RATIO_MIN}")

    # ---- the published call set --------------------------------------------- #
    book, parser = read_xlsx(XLSX)
    missing = [s for s in CELL_LINES if s not in book]
    if missing:
        raise SystemExit(f"sheets missing from {XLSX.name}: {missing}")
    #: the .fai already carries the ``chr`` prefix on this assembly; index it both
    #: ways so the lookup cannot depend on which name a source happens to use
    chrom_len: dict[str, int] = {}
    for line in FAI.read_text().splitlines():
        f = line.split("\t")
        chrom_len[f[0]] = int(f[1])
        chrom_len.setdefault("chr" + f[0], int(f[1]))
    spans = gene_spans()
    span_ix = spans.groupby("gene").agg(gchrom=("chrom", lambda s: s.mode().iat[0]),
                                        gstart=("start", "min"), gend=("end", "max"))

    site_rows, coord_rows = [], []
    for sheet, cell in CELL_LINES.items():
        raw = book[sheet]
        raw.columns = [str(c) for c in raw.columns]
        rec = pd.DataFrame({
            "pos_id": raw["0"].astype(str),
            "gene": raw["1"].astype(str),
            "ratio": pd.to_numeric(raw["2"], errors="coerce"),
            "mod_count": pd.to_numeric(raw["3"], errors="coerce"),
            "total_count": pd.to_numeric(raw["4"], errors="coerce"),
        })
        rec.insert(0, "cell_line", cell)
        ex = rec.pos_id.str.extract(POS_RE)
        rec["chrom"] = ex[0]
        rec["pos"] = pd.to_numeric(ex[1], errors="coerce")
        rec["parse_ok"] = rec.chrom.notna() & rec.pos.notna()
        rec["pos0"] = (rec.pos - 1).round()
        rec["in_chr_length"] = [
            bool(p) and c in chrom_len and p <= chrom_len[c]
            for c, p in zip(rec.chrom, rec.pos)
        ]
        j = rec.join(span_ix, on="gene")
        rec["gene_named"] = j.gstart.notna().to_numpy()
        #: the GTF writes bare chromosome numbers and the deposited table writes
        #: ``chr1``; the check is on the position inside the span of the gene the
        #: sheet itself names, on the chromosome it names
        rec["in_named_gene"] = [
            bool(n) and bool(ok) and c is not None and gc is not None
            and str(c).removeprefix("chr") == str(gc).removeprefix("chr")
            and p is not None and (p - 1) >= gs and (p - 1) <= ge
            for c, p, ok, n, gc, gs, ge in zip(rec.chrom, rec.pos, rec.parse_ok,
                                               rec.gene_named, j.gchrom, j.gstart, j.gend)
        ]
        site_rows.append(rec)
        parsed = rec[rec.parse_ok]
        coord_rows.append({
            "cell_line": cell, "sheet": sheet,
            "n_rows": len(rec), "n_parsed": int(rec.parse_ok.sum()),
            "n_unparsed": int((~rec.parse_ok).sum()),
            "frac_unparsed": round(float((~rec.parse_ok).mean()), 4),
            "n_in_chromosome_length": int(rec.in_chr_length.sum()),
            "n_gene_named": int(rec.gene_named.sum()),
            "n_in_named_gene": int(rec.in_named_gene.sum()),
            "frac_in_named_gene": round(
                float(rec.in_named_gene.sum() / max(int(rec.parse_ok.sum()), 1)), 4),
            "n_genes": int(parsed.gene.nunique()),
        })
        print(f"[{cell}] rows {len(rec)} parsed {int(rec.parse_ok.sum())} "
              f"({rec.in_named_gene.mean():.1%} inside the named gene)")

    sites = pd.concat(site_rows, ignore_index=True)
    coord = pd.DataFrame(coord_rows)
    ok = coord[coord.cell_line == "HeLa"].iloc[0]
    if ok.frac_in_named_gene < MIN_IN_NAMED_GENE:
        raise SystemExit(f"published coordinates do not land in the named genes "
                         f"({ok.frac_in_named_gene:.1%} < {MIN_IN_NAMED_GENE:.1%}); "
                         f"the assembly assumption is wrong -- stop")
    if coord.frac_unparsed.max() > MAX_UNPARSED:
        raise SystemExit(f"too many unparseable rows "
                         f"({coord.frac_unparsed.max():.1%} > {MAX_UNPARSED:.1%})")

    sites = sites[["cell_line", "chrom", "pos", "pos0", "gene", "ratio", "mod_count",
                   "total_count", "parse_ok", "in_chr_length", "in_named_gene"]]
    sites.to_csv(TAB / "s11_published_sites.tsv", sep="\t", index=False)
    coord.to_csv(TAB / "s11_coordinate_check.tsv", sep="\t", index=False)

    # ---- overlap and agreement, site level ---------------------------------- #
    hela = sites[(sites.cell_line == "HeLa") & sites.parse_ok].copy()
    hela["key"] = hela.chrom + ":" + hela.pos0.astype("int64").astype(str)
    hela = hela.drop_duplicates("key").rename(columns={
        "ratio": "published_ratio", "total_count": "published_total_reads",
        "gene": "published_gene"})
    published = set(hela.key)
    ov_rows, ag_rows = [], []
    for _folder, label, _units in CONDITIONS:
        mine = ours[label]
        shared = mine.merge(hela[["key", "published_ratio", "published_total_reads",
                                  "published_gene"]], on="key")
        inter = len(shared)
        n_ours = len(mine)
        ov_rows.append({
            "condition": label, "ratio_min": RATIO_MIN,
            "n_our_sites": n_ours, "n_shared": inter,
            "frac_ours_shared": round(inter / n_ours, 4),
            "n_published_unique": len(published),
            "frac_published_shared": round(inter / len(published), 4),
            "n_published_not_ours": len(published) - inter,
            "jaccard": round(inter / (n_ours + len(published) - inter), 4),
        })
        ag_rows.append({
            "condition": label, "n_shared": inter,
            "pearson_log10_ratio": round(float(np.corrcoef(
                np.log10(shared.our_ratio),
                np.log10(shared.published_ratio.clip(lower=1e-3)))[0, 1]), 4),
            "spearman_ratio": round(float(shared.our_ratio.corr(
                shared.published_ratio, method="spearman")), 4),
            "median_our_ratio": round(float(shared.our_ratio.median()), 4),
            "median_published_ratio": round(float(shared.published_ratio.median()), 4),
            "q1_our_ratio": round(float(shared.our_ratio.quantile(.25)), 4),
            "q3_our_ratio": round(float(shared.our_ratio.quantile(.75)), 4),
            "q1_published_ratio": round(float(shared.published_ratio.quantile(.25)), 4),
            "q3_published_ratio": round(float(shared.published_ratio.quantile(.75)), 4),
            "median_published_total_reads": int(shared.published_total_reads.median()),
        })
        print(f"[{label}] ours {n_ours}, shared {inter} "
              f"({inter / n_ours:.1%} of ours, {inter / len(published):.1%} of published), "
              f"r={ag_rows[-1]['pearson_log10_ratio']}")
    ov = pd.DataFrame(ov_rows)
    ag = pd.DataFrame(ag_rows)
    ov.to_csv(TAB / "s11_overlap.tsv", sep="\t", index=False)
    ag.to_csv(TAB / "s11_ratio_agreement.tsv", sep="\t", index=False)
    #: the site-level pairs panel B draws, kept as a table so the figure and the
    #: verifier read the same rows rather than each rebuilding the merge
    shared_all = pd.concat([shared.assign(condition=label)
                            for _f, label, _u in CONDITIONS
                            for shared in [ours[label].merge(
                                hela[["key", "published_ratio", "published_total_reads",
                                      "published_gene"]], on="key")]])
    shared_all[["condition", "key", "published_gene", "our_ratio", "published_ratio",
                "published_total_reads"]].to_csv(
        TAB / "s11_shared_sites.tsv", sep="\t", index=False)

    # ---- the distributions panel C of the standalone figure drew, and the
    
    # distribution moved into Figure S8D as NanoNm's wild-type side, so the 5th and
    # 95th percentiles are needed there (they are the band's whisker ends) and the
    # published HeLa row is written out on its own for that panel to read.
    dist_rows = []
    for _sheet, cell in CELL_LINES.items():
        d = sites[(sites.cell_line == cell) & sites.parse_ok].ratio.dropna()
        dist_rows.append({
            "source": "published", "label": cell, "n": len(d),
            "min": round(float(d.min()), 4),
            "q05": round(float(d.quantile(.05)), 4),
            "q1": round(float(d.quantile(.25)), 4),
            "median": round(float(d.median()), 4),
            "q3": round(float(d.quantile(.75)), 4),
            "q95": round(float(d.quantile(.95)), 4),
            "max": round(float(d.max()), 4),
            "median_total_reads": int(sites[(sites.cell_line == cell)
                                            & sites.parse_ok].total_count.median()),
        })
    for _folder, label, _units in CONDITIONS:
        d = ours[label].our_ratio
        dist_rows.append({
            "source": "this study", "label": label, "n": len(d),
            "min": round(float(d.min()), 4),
            "q05": round(float(d.quantile(.05)), 4),
            "q1": round(float(d.quantile(.25)), 4),
            "median": round(float(d.median()), 4),
            "q3": round(float(d.quantile(.75)), 4),
            "q95": round(float(d.quantile(.95)), 4),
            "max": round(float(d.max()), 4), "median_total_reads": pd.NA,
        })
    dist = pd.DataFrame(dist_rows)
    dist.to_csv(TAB / "s11_ratio_distribution.tsv", sep="\t", index=False)
    hela_pub = dist[(dist.source == "published") & (dist.label == "HeLa")]
    if len(hela_pub) != 1:
        raise SystemExit("the published HeLa distribution row is missing or duplicated")
    hela_pub.to_csv(TAB / "s11_published_hela_score.tsv", sep="\t", index=False)

    # ---- the anchors the caption quotes ------------------------------------- #
    anchors = [
        ("published_hek293t_sites", int(coord.loc[coord.cell_line == "HEK293T",
                                                 "n_parsed"].iloc[0])),
        ("published_hela_sites", int(ok.n_parsed)),
        ("published_hepG2_sites", int(coord.loc[coord.cell_line == "HepG2",
                                                "n_parsed"].iloc[0])),
        ("published_hela_unique_sites", len(published)),
        ("published_hela_genes", int(ok.n_genes)),
        ("unparsed_rows_total", int(coord.n_unparsed.sum())),
        ("frac_in_named_gene_min", float(coord.frac_in_named_gene.min())),
        ("our_hela_wt_sites", int(ov.loc[ov.condition == "HeLa WT", "n_our_sites"].iloc[0])),
        ("our_hela_ivt_sites", int(ov.loc[ov.condition == "HeLa IVT", "n_our_sites"].iloc[0])),
        ("shared_wt", int(ov.loc[ov.condition == "HeLa WT", "n_shared"].iloc[0])),
        ("shared_ivt", int(ov.loc[ov.condition == "HeLa IVT", "n_shared"].iloc[0])),
        ("frac_ours_shared_wt", float(ov.loc[ov.condition == "HeLa WT",
                                             "frac_ours_shared"].iloc[0])),
        ("frac_ours_shared_ivt", float(ov.loc[ov.condition == "HeLa IVT",
                                              "frac_ours_shared"].iloc[0])),
        ("frac_published_shared_wt", float(ov.loc[ov.condition == "HeLa WT",
                                                  "frac_published_shared"].iloc[0])),
        ("frac_published_shared_ivt", float(ov.loc[ov.condition == "HeLa IVT",
                                                   "frac_published_shared"].iloc[0])),
        ("pearson_log10_wt", float(ag.loc[ag.condition == "HeLa WT",
                                          "pearson_log10_ratio"].iloc[0])),
        ("pearson_log10_ivt", float(ag.loc[ag.condition == "HeLa IVT",
                                           "pearson_log10_ratio"].iloc[0])),
        ("spearman_wt", float(ag.loc[ag.condition == "HeLa WT", "spearman_ratio"].iloc[0])),
        ("spearman_ivt", float(ag.loc[ag.condition == "HeLa IVT", "spearman_ratio"].iloc[0])),
        ("published_hela_median_ratio", float(dist.loc[(dist.source == "published")
                                                       & (dist.label == "HeLa"),
                                                       "median"].iloc[0])),
        ("published_hek293t_median_ratio", float(dist.loc[(dist.source == "published")
                                                         & (dist.label == "HEK293T"),
                                                         "median"].iloc[0])),
        ("published_hepG2_median_ratio", float(dist.loc[(dist.source == "published")
                                                        & (dist.label == "HepG2"),
                                                        "median"].iloc[0])),
        ("published_hek293t_median_reads", int(dist.loc[(dist.source == "published")
                                                        & (dist.label == "HEK293T"),
                                                        "median_total_reads"].iloc[0])),
        ("published_hela_median_reads", int(dist.loc[(dist.source == "published")
                                                     & (dist.label == "HeLa"),
                                                     "median_total_reads"].iloc[0])),
        ("published_hepG2_median_reads", int(dist.loc[(dist.source == "published")
                                                      & (dist.label == "HepG2"),
                                                      "median_total_reads"].iloc[0])),
        ("shared_published_median_wt", float(ag.loc[ag.condition == "HeLa WT",
                                                    "median_published_ratio"].iloc[0])),
        ("shared_published_median_ivt", float(ag.loc[ag.condition == "HeLa IVT",
                                                     "median_published_ratio"].iloc[0])),
        ("frac_unparsed_hela", float(coord.loc[coord.cell_line == "HeLa",
                                                "frac_unparsed"].iloc[0])),
        ("frac_in_named_gene_hela", float(coord.loc[coord.cell_line == "HeLa",
                                                    "frac_in_named_gene"].iloc[0])),
        ("our_wt_median_ratio", float(ag.loc[ag.condition == "HeLa WT",
                                             "median_our_ratio"].iloc[0])),
        ("our_ivt_median_ratio", float(ag.loc[ag.condition == "HeLa IVT",
                                              "median_our_ratio"].iloc[0])),
    ]
    pd.DataFrame(anchors, columns=["key", "value"]).to_csv(
        LOG / "73_tables_anchors.tsv", sep="\t", index=False)

    # ---- provenance ---------------------------------------------------------- #
    with (SRC / "mmc6_fingerprint.tsv").open("w", encoding="utf-8") as fh:
        fh.write("role\tpath\tbytes\tmd5\tdetail\n")
        fh.write(f"source\tthe deposited table\t{XLSX.stat().st_size}\t{md5(XLSX)}\t"
                 f"sheets {len(CELL_LINES)}: {', '.join(CELL_LINES.values())}; "
                 f"parser {parser}\n")
        for role, path, digest, detail in prov:
            fh.write(f"{role}\t{path}\t{Path(path).stat().st_size}\t{digest}\t{detail}\n")
        fh.write(f"annotation\t{GTF}\t{GTF.stat().st_size}\t-\tgene_name spans "
                 f"({len(spans)} named genes)\n")
    print(f"[done] tables + anchors in {time.time() - t0:.1f} s "
          f"(parser: {parser})")
    print(f"[log] {LOG / '73_tables_anchors.tsv'}")


if __name__ == "__main__":
    sys.exit(main())
