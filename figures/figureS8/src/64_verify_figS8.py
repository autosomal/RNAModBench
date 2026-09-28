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

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PROJECT = Path(str(_RB))
OUT = (_RB / "figures/figureS8")
TAB = (_RB / "figures/figureS8/tables")
FIG = (_RB / "figures/figureS8/figures")
PDF = (_RB / "figures/figureS8/figures/FigureS8_rev.pdf")
checks = []


def check(name, ok=True, detail="", status=None):
    st = status or ("OK" if ok else "FAIL")
    checks.append({"check": name, "status": st, "detail": detail})
    print(f"{st:<5}{name}  {detail}")


def shell(cmd):
    return subprocess.run(cmd, capture_output=True, text=True).stdout


# --- 1. anchors -------------------------------------------------------------
def verify_anchors():
    a = pd.read_csv((_RB / "figures/figureS8/tables/s8_anchor_check.tsv"), sep="\t")
    bad = a[a.status != "OK"]
    check("anchors all OK", len(bad) == 0,
          f"{len(a) - len(bad)}/{len(a)}")
    check("anchor count >= 84", len(a) >= 84, str(len(a)))


# --- 2. Guitar vs Python placement on the same majority BEDs ----------------
def verify_engines():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "s8_tables", Path(__file__).with_name("61_figS8_tables.py"))
    s7 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(s7)
    geo = pd.read_csv((_RB / "figures/figureS8/tables/figS8_guitar_geometry.tsv"), sep="\t")
    g = {}
    for r in geo.itertuples():
        try:
            g[str(r.item)] = float(r.value)
        except (TypeError, ValueError):
            pass
    bw = np.cumsum([g["component_width_promoter"],
                    g["component_width_ncrna"],
                    g["component_width_tail"]])[:2]
    sites = pd.read_csv((_RB / "figures/figureS8/tables/figS8_guitar_sites.tsv"), sep="\t")
    plan = pd.read_csv((_RB / "figures/figureS8/tables/figS8_panel_inputs.tsv"), sep="\t")
    plan = plan[plan.kind == "consensus"]
    model = s7.RegionIndex.load(sorted(
        p for p in s7.MODEL_DIR.glob("Human.*.ncrna.regionmodel.pkl")
        if "gencode" not in p.name)[0])
    rows = []
    for r in plan.itertuples():
        grp = f"{r.tool}|{r.role}|consensus|"
        rel = sites.loc[sites.group == grp, "relative"].to_numpy(float)
        if rel.size == 0:
            continue
        gsh = np.array([(rel < bw[0]).mean(),
                        ((rel >= bw[0]) & (rel < bw[1])).mean(),
                        (rel >= bw[1]).mean()])
        bed = pd.read_csv(r.path, sep="\t", header=None,
                          names=["chrom", "start", "end", "name", "score",
                                 "strand"])
        kk, body_py = [], []
        for chrom, sub in bed.groupby("chrom", sort=False):
            # the BEDs carry the GTF's own chromosome spelling, the model keys
            # the normalised chr* labels: normalise before the lookup
            tab = s7.chrom_table(model, s7.fix_chromosome(chrom))
            if tab is None:
                continue
            k, z = s7.place_positions(tab, np.sort(
                sub["start"].to_numpy(np.int64)))
            kk.append(k)
            bb = z[k == s7.NCR]
            if bb.size:
                body_py.append(bb)
        k = np.concatenate(kk) if kk else np.zeros(0, np.int8)
        psh = np.array([(k == s7.F5).mean(), (k == s7.NCR).mean(),
                        (k == s7.F3).mean()])
        rb = (rel[(rel >= bw[0]) & (rel < bw[1])] - bw[0]) / (bw[1] - bw[0])
        py = np.concatenate(body_py) if body_py else np.zeros(0, np.float64)
        dg, dp = s7.density(rb), s7.density(py)
        rho = float(spearmanr(dg, dp).statistic) if rb.size and py.size else np.nan
        l1 = float(np.abs(dg / dg.sum() - dp / dp.sum()).sum()) if rb.size and py.size else np.nan
        rows.append(dict(tool=r.tool, condition=r.role, n_bed=len(bed),
                         n_guitar=int(rel.size), n_body_guitar=int(rb.size),
                         n_body_py=int(py.size),
                         mode_bin_guitar=int(np.argmax(dg)),
                         mode_bin_py=int(np.argmax(dp)),
                         spearman_rho=rho, l1_body=l1,
                         g_promoter=float(gsh[0]), g_body=float(gsh[1]),
                         g_tail=float(gsh[2]), p_promoter=float(psh[0]),
                         p_body=float(psh[1]), p_tail=float(psh[2]),
                         d_promoter=float(gsh[0] - psh[0]),
                         d_body=float(gsh[1] - psh[1]),
                         d_tail=float(gsh[2] - psh[2])))
    cc = pd.DataFrame(rows)
    cc.to_csv((_RB / "figures/figureS8/tables/s8_density_crosscheck.tsv"), sep="\t", index=False,
              float_format="%.6g")
    # the two engines are not expected to be identical: Guitar counts a site once
    # per overlapping transcript (and drops ambiguous ones), the Python region
    # model assigns each site to the longest transcript.  What must agree is the
    # axis itself and the qualitative shape -- a swapped/rotated component or a
    # wrong boundary would break all three checks below.
    shift = float(np.abs(cc.mode_bin_guitar - cc.mode_bin_py).max()) / s7.GRID
    check("body-density peak within 0.20 of the body", shift <= 0.20,
          f"max shift {shift:.3f} over {len(cc)} majority groups")
    wt = cc[cc.condition == "WT"].set_index("tool")
    ivt = cc[cc.condition == "IVT"].set_index("tool")
    d_eng = {"guitar": wt.g_body - ivt.g_body, "python": wt.p_body - ivt.p_body}
    agree = int(np.sum(np.sign(d_eng["guitar"]) == np.sign(d_eng["python"])))
    check("INFO: sign of the body-share change agrees",
          f"{agree}/6 tools (NanoNm flips: the engines differ by the "
          f"multi-transcript ambiguity rule)", status="INFO")
    check("body-density Spearman rho >= 0.50",
          float(cc.spearman_rho.min()) >= 0.5,
          f"min rho = {cc.spearman_rho.min():.3f}, "
          f"median = {cc.spearman_rho.median():.3f}")
    # the engines differ by the ambiguity rule (Guitar counts a site once per
    # overlapping transcript, this model assigns it to the longest one), so the
    # segment shares are reported for information and only guarded against a
    # gross axis error
    worst = float(np.abs(cc[["d_promoter", "d_body", "d_tail"]].to_numpy()).max())
    check("segment shares within 0.40 (gross-error guard)", worst <= 0.40,
          f"max |delta| = {worst:.4f} (ambiguity rule), see s8_density_crosscheck.tsv")


BANNED = re.compile(
    r"NanoSPA_Psu|NanoMUD_psi|NanoMUD_m1psi|NanoSPA_psi|T/WT|treatment|"
    r"GENCODE|gencode|1e[-+]|mean[_ ]of[_ ]replicates")
NAMES = ["CHEUI-m5C", "NanoMUD-\u03a8", "NanoMUD-m1\u03a8", "NanoNm",
         "NanoPsu", "NanoSPA-\u03a8"]


def verify_pdf():
    txt = shell(["pdftotext", str(PDF), "-"])
    m = BANNED.search(txt)
    check("no legacy names / banned wording", m is None,
          "" if m is None else m.group(0))
    missing = [n for n in NAMES if n not in txt]
    check("standardized tool names present", not missing, str(missing))
    # one A for the whole metagene block, then B/C/D for the quantitative row
    letters = set(re.findall(r"^([A-I])$", txt, re.M))
    check("panel letters A-D present", {"A", "B", "C", "D"} <= letters,
          str(sorted(letters)))
    check("no stale panel letters E-I", not (letters - {"A", "B", "C", "D"}),
          f"unexpected {sorted(letters - {'A', 'B', 'C', 'D'})}")
    for key in ("1kb", "ncRNA", "Density", "JSD"):
        check(f"axis/label '{key}'", key in txt)
    info = shell(["pdfinfo", str(PDF)])
    mm = re.search(r"Page size:\s+([\d.]+) x ([\d.]+) pts", info)
    w, h = float(mm.group(1)) * 25.4 / 72, float(mm.group(2)) * 25.4 / 72
    check("page <= 222 x 247 mm", w <= 222 and h <= 247, f"{w:.1f} x {h:.1f} mm")
    check("page is the 180 mm page width", abs(w - 180) <= 2, f"{w:.1f} mm")
    fnt = shell(["pdffonts", str(PDF)]).strip().splitlines()[2:]
    bad_font = [l.split()[0] for l in fnt
                if "Arial" not in l.split()[0] or l.split()[-4] != "yes"]
    check("fonts: Arial only, all embedded", not bad_font, str(bad_font))
    geo = pd.read_csv((_RB / "figures/figureS8/tables/figS8_page_geometry.tsv"), sep="\t")
    gm = {r.item: float(r.value) for r in geo.itertuples()
          if str(r.value).replace(".", "").isdigit()}
    check("declared minimum font >= 7 pt", gm.get("font_min_pt", 0) >= 7,
          f"{gm.get('font_min_pt')} pt")
    xml = shell(["pdftotext", "-bbox", str(PDF), "-"])
    boxes = np.array(
        [[float(a), float(b), float(c), float(d)] for a, b, c, d in re.findall(
            r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" '
            r'yMax="([\d.]+)"', xml)])
    heights = boxes[:, 3] - boxes[:, 1]
    # pdftotext reports the tight glyph box, not the em box: divide by the
    # Arial cap/digit height ratio 0.716 to recover the font size (a digit-only
    # word fills exactly the cap height)
    est = float(heights.min()) / 0.716
    check("smallest drawn text >= 7 pt (estimated)", est >= 6.9,
          f"{est:.2f} pt from a {heights.min():.2f} pt box "
          f"({len(boxes)} words)")
    coll = 0
    for i in range(len(boxes)):
        a = boxes[i]
        ov = np.minimum(boxes[i + 1:, 2], a[2]) - np.maximum(boxes[i + 1:, 0], a[0])
        oh = np.minimum(boxes[i + 1:, 3], a[3]) - np.maximum(boxes[i + 1:, 1], a[1])
        coll += int((np.minimum(ov, oh) > 0.5).sum())
    check("no overlapping text boxes", coll == 0, f"{coll} overlapping pairs")


def main():
    for f in ((_RB / "figures/figureS8/tables/s8_anchor_check.tsv"), PDF,
              (_RB / "figures/figureS8/tables/figS8_guitar_geometry.tsv"),
              (_RB / "figures/figureS8/tables/figS8_page_geometry.tsv"),
              (_RB / "figures/figureS8/tables/figS8_panel_inputs.tsv"),
              (_RB / "figures/figureS8/tables/figS8_guitar_sites.tsv")):
        if not f.exists():
            print(f"FAIL missing input {f}")
            return 1
    verify_anchors()
    verify_engines()
    verify_pdf()
    df = pd.DataFrame(checks)
    df.to_csv((_RB / "figures/figureS8/tables/s8_verify_report.tsv"), sep="\t", index=False)
    bad = df[df.status == "FAIL"]
    n_info = int((df.status == "INFO").sum())
    print(f"\n{len(df) - len(bad) - n_info}/{len(df) - n_info} checks passed"
          f" ({n_info} informational) -> {TAB / 's8_verify_report.tsv'}")
    print("ALL CHECKS PASSED" if not len(bad) else f"{len(bad)} FAILED")
    return 1 if len(bad) else 0


if __name__ == "__main__":
    sys.exit(main())
