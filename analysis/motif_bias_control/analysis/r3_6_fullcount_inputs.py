#!/usr/bin/env python
"""Derive figure-input tables from the full-count callsets analysis
(r3_6_callsets_kl.py outputs), in the same schemas as the legacy
Top-5-based inputs, so r3_6_figures.py can render either dataset.

Outputs: fullcount_kl_pooled.tsv, fullcount_top5.tsv, fullcount_ggac.tsv
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
import glob
import os
from collections import defaultdict

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CLEAN = str(_RB / "data/callsets/RNA002")
SP = {"Arabidopsis": "Arabidopsis", "Mouse": "Mouse", "Human (HeLa)": "Human"}
COMP = str.maketrans("ACUG", "UGAC")


def norm5(s):
    if s[2] == "U":
        return s.translate(COMP)[::-1]
    return s


def main():
    kl = pd.read_csv(os.path.join(HERE, "fullcount_kl.tsv"), sep="\t")
    pooled = kl.groupby(["Tool", "Species_Name"], as_index=False)["KL_Divergence"].mean()
    pooled.to_csv(os.path.join(HERE, "fullcount_kl_pooled.tsv"), sep="\t", index=False)

    counts = defaultdict(lambda: defaultdict(int))
    for spname, spdir in SP.items():
        for path in sorted(glob.glob(f"{CLEAN}/{spdir}/*_WT/m6A/*/*.tsv")):
            tool = path.split("/m6A/")[1].split("/")[0]
            df = pd.read_csv(path, sep="\t", usecols=["five_mer_raw"], dtype=str)
            for s in df["five_mer_raw"].str.replace("T", "U", regex=False):
                if len(s) == 5:
                    k = norm5(s)
                    if k[2] == "A":
                        counts[(spname, tool)][k] += 1
    rows = []
    for (sp, tool), c in counts.items():
        tot = sum(c.values())
        for rank, (km, n) in enumerate(sorted(c.items(), key=lambda x: -x[1])[:5], 1):
            rows.append({"Species": sp, "Tool": tool, "Rank": rank, "5mer": km,
                         "Frequency": n, "Relative_Frequency": n / tot})
    top5 = pd.DataFrame(rows)
    top5.to_csv(os.path.join(HERE, "fullcount_top5.tsv"), sep="\t", index=False)

    # AGAC/GGAC fractions over ALL normalized 5-mers (no truncation)
    full = []
    for (sp, tool), c in counts.items():
        tot = sum(c.values())
        ag = sum(n for k, n in c.items() if k[:4] in ("AGAC",)) / tot
        gg = sum(n for k, n in c.items() if k[:4] == "GGAC") / tot
        full.append({"Tool": tool, "Species": sp, "AGAC": ag, "GGAC": gg})
    ff = pd.DataFrame(full)
    out = ff.pivot_table(index="Tool", columns=["Species"], values=["AGAC", "GGAC"])
    out.columns = [f"{c2}_{c1}" for c1, c2 in out.columns]
    out = out.fillna(0.0)
    for sp in SP:
        out[f"{sp}_AGAC_minus_GGAC"] = out[f"{sp}_AGAC"] - out[f"{sp}_GGAC"]
    out.to_csv(os.path.join(HERE, "fullcount_ggac.tsv"), sep="\t")
    print("wrote fullcount_kl_pooled.tsv, fullcount_top5.tsv, fullcount_ggac.tsv")


if __name__ == "__main__":
    main()
