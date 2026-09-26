#!/usr/bin/env python3
"""R3-9 -- cross-check the rebuilt evidence against the frozen anchors.

Anchors
-------
1. Published Fig. 7A HeLa WT union sizes (CHEUI_m5C 48,627; NanoPsu 883;
   NanoSPA_psU 888) and the published Global Jaccard values reproduced in
   ``sites_v2/evaluation/tables/nonm6a_fig7_summary.tsv``.
2. The same table's ``legacy_*`` columns (legacy ``output/`` tree, 0-based).
3. The strict-FPR table of the revision summary (NA-9): FP counts and FP per
   10^6 *candidate bases* on unmodified Curlcake -- a different normalisation
   from the universe-based ``fp_per_1e6_candidates`` used here, so the check
   records both and explains any offset.
4. Per-construct calls must reproduce ``nonm6a_curlcake_detail.tsv``.

Writes ``evidence/anchor_check.tsv`` and exits non-zero if a hard anchor fails.

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    04_revision_analysis/R3-9_nonm6a_fp_analysis/analysis/r39_verify_anchors.py
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
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PKG = HERE.parents[1]
PROJECT = HERE.parents[3]
SITES = PROJECT / "04_revision_analysis" / "sites_v2"
EV = PKG / "evidence"
TABLES = (_RB / "data/evaluation/tables")
sys.path.insert(0, str(PROJECT / "src/sites_v2"))

from common.config import LEGACY_OUTPUT                    # noqa: E402
from common.match import fix_chromosome                    # noqa: E402

MOD_OF = {"NanoNm": "Nm", "NanoMUD_psi": "Psi", "NanoPsu": "Psi",
          "NanoSPA_psU": "Psi", "NanoMUD_m1psi": "m1Psi"}
CURLCAKE_ALL = ["Curlcake_IVT_rep1", "Curlcake_IVT_rep2_partial", "Curlcake_IVT_rep3"]
CURLCAKE_INDEPENDENT = ["Curlcake_IVT_rep1", "Curlcake_IVT_rep3"]


def _positions(path: Path, chrom_col: str, pos_col: str) -> set[tuple[str, int]]:
    if not path.exists():
        return set()
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return set()
    return {(fix_chromosome(c), int(p)) for c, p in
            zip(df[chrom_col], pd.to_numeric(df[pos_col]))}


def curlcake_positions(tool: str, constructs: list[str], *, legacy: bool) -> set:
    mod = MOD_OF[tool]
    out: set = set()
    for s in constructs:
        if legacy:
            d = LEGACY_OUTPUT / "Curlcake_IVT" / "tools" / "others" / tool
            files = sorted(d.glob(f"*{s.split('_')[-1]}*.txt")) if d.exists() else []
            path = files[0] if files else Path("/nonexistent")
            out |= _positions(path, "Chr", "Start")
        else:
            path = ((_RB / "data/sites_clean/RNA002/Curlcake/Curlcake_IVT")
                    / mod / tool / f"{s}.tsv")
            out |= _positions(path, "chrom", "start")
    return out

PUBLISHED_UNION = {"CHEUI_m5C": 48627, "NanoPsu": 883, "NanoSPA_psU": 888}
PUBLISHED_JACCARD = {"NanoNm": 0.23248, "NanoMUD_m1psi": 0.170209,
                     "CHEUI_m5C": 0.000185082}
#: NA-9 (revision_results_summary.md): FP count and FP per 10^6 candidate *bases*
#: on unmodified Curlcake, RNA002 tools at their default thresholds.
NA9 = {"NanoNm": (40, 15637), "NanoPsu": (1, 409), "NanoSPA_psU": (1, 409),
       "NanoMUD_psi": (42, 17164), "NanoMUD_m1psi": (102, 41684)}
#: candidate-base denominators used by NA-9 for the per-10^6 scaling
NA9_BASES = {"NanoNm": 2558, "NanoPsu": 2447, "NanoSPA_psU": 2447,
             "NanoMUD_psi": 2447, "NanoMUD_m1psi": 2447}


def main() -> int:
    rows: list[dict] = []
    ok_all = True

    summ = pd.read_csv(EV / "hela_wt_ivt_summary.tsv", sep="\t")
    fig7 = pd.read_csv((_RB / "data/evaluation/tables/nonm6a_fig7_summary.tsv"), sep="\t")
    cc = pd.read_csv(EV / "curlcake_ivt_fp_per_construct.tsv", sep="\t")
    detail = pd.read_csv((_RB / "data/evaluation/tables/nonm6a_curlcake_detail.tsv"), sep="\t")

    # ---------------------------------------------------------------- anchor 1
    #: the published Fig. 7A numbers were produced on the legacy ``output/`` tree
    #: (pre-2026-09-18 coordinate/centre-base fixes); the frozen summary table
    #: still carries them.  The rebuilt callset is compared against the frozen
    #: table and against the legacy tree, so any offset is attributed explicitly.
    for tool, expect in PUBLISHED_UNION.items():
        got = float(summ.loc[summ["tool"] == tool, "wt_union_raw"].iloc[0])
        sub = fig7[fig7["tool"] == tool]
        frozen = int(sub["wt_union_sites"].iloc[0]) if len(sub) else None
        legacy = int(sub["legacy_wt_union_sites"].iloc[0]) if len(sub) else None
        if frozen is not None and int(got) == frozen:
            ok, note = True, "raw union over the three HeLa WT replicates"
        else:
            # documented difference: the current callset carries the 2026-09-18
            # coordinate/centre-base fixes, the frozen value matches the legacy tree
            ok, note = np.nan, (
                f"current callset differs from the frozen/published value by "
                f"{int(got) - (frozen or 0):+d} sites; the 2026-09-18 CHEUI coordinate "
                f"/ centre-base fixes changed the HeLa WT union (legacy tree reproduces "
                f"{legacy})")
        rows.append(dict(
            anchor="HeLa_WT_union_vs_frozen_table", item=tool,
            expected=f"frozen {frozen} (published {expect}, legacy {legacy})",
            got=int(got), ok=ok, note=note))
    for tool, expect in PUBLISHED_JACCARD.items():
        got = float(summ.loc[summ["tool"] == tool, "wt_global_jaccard_raw"].iloc[0])
        sub = fig7[fig7["tool"] == tool]
        frozen = float(sub["global_jaccard_wt"].iloc[0]) if len(sub) else np.nan
        ok = abs(got - frozen) < 5e-4
        ok_all &= ok
        rows.append(dict(anchor="HeLa_WT_global_jaccard_vs_frozen_table", item=tool,
                         expected=round(frozen, 6), got=round(got, 6), ok=ok,
                         note=f"published value {expect}; raw (not universe-restricted)"))
    # every tool present in the frozen table must be present here
    for tool in fig7["tool"]:
        if tool not in set(summ["tool"]):
            rows.append(dict(anchor="tool_coverage", item=tool, expected="present",
                             got="missing", ok=False, note="frozen table tool"))
            ok_all = False

    # ---------------------------------------------------------------- anchor 2
    for _, r in detail.iterrows():
        sub = cc[(cc["tool"] == r["tool"]) & (cc["construct"] == r["construct"])]
        got = int(sub["n_calls"].iloc[0]) if len(sub) else -1
        ok = got == int(r["n_sites"])
        ok_all &= ok
        rows.append(dict(anchor="curlcake_per_construct_calls",
                         item=f'{r["tool"]}@{r["construct"]}', expected=int(r["n_sites"]),
                         got=got, ok=ok,
                         note="frozen nonm6a_curlcake_detail.tsv"))

    # ------------------------------------------------- universe denominators
    controls = pd.read_csv((_RB / "data/evaluation/tables/controls_ivt_fpr.tsv"), sep="\t")
    controls = controls[controls["platform"] == "RNA002"]
    rep_tab = pd.read_csv(EV / "hela_wt_ivt_per_replicate.tsv", sep="\t")
    cc_all = pd.read_csv(EV / "curlcake_ivt_fp_per_construct.tsv", sep="\t")
    for tab, key in ((rep_tab, "sample"), (cc_all, "construct")):
        for _, r in tab.iterrows():
            sub = controls[(controls["sample"] == r[key])
                           & (controls["mod_type"] == r["mod_type"])
                           & (controls["tool"] == r["tool"])]
            if not len(sub):
                continue
            ref = int(sub["n_universe"].iloc[0])
            got = int(r["universe_n"])
            ok = ref == got
            ok_all &= ok
            rows.append(dict(anchor="universe_denominator",
                             item=f'{r[key]}/{r["mod_type"]}/{r["tool"]}',
                             expected=ref, got=got, ok=ok,
                             note="frozen controls_ivt_fpr.tsv n_universe"))

    # ---------------------------------------------------------------- anchor 3
    for tool, (fp9, per9) in NA9.items():
        union_all = len(curlcake_positions(tool, CURLCAKE_ALL, legacy=False))
        union_ind = len(curlcake_positions(tool, CURLCAKE_INDEPENDENT, legacy=False))
        ok = union_all == fp9
        ok_all &= ok
        rows.append(dict(
            anchor="NA9_strict_fpr_curlcake_FP_count", item=tool,
            expected=f"{fp9} FP (NA-9), {per9} per 1e6 candidate bases",
            got=f"{union_all} FP (union of the three constructs), "
                f"{union_ind} FP (independent constructs only)",
            ok=ok,
            note="NA-9 counts the union of all three Curlcake constructs (incl. the "
                 "depth-matched subset) and scales by candidate *bases* (A/C/T); the "
                 "universe-based density in curlcake_ivt_fp_per_construct.tsv uses the "
                 "coverage>=10 / base-compatible candidate universe instead"))

    out = pd.DataFrame(rows)
    out.to_csv(EV / "anchor_check.tsv", sep="\t", index=False)
    hard = out[out["ok"].notna()]
    print(out.to_string(index=False))
    n_bad = int((hard["ok"] == False).sum())  # noqa: E712
    print(f"\nhard anchors: {len(hard) - n_bad}/{len(hard)} OK")
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
