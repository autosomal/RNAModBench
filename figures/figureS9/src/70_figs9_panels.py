#!/usr/bin/env python
"""70 -- render the nine Figure S9 panels (A-I), one file at a time, at print size.

Each panel is drawn on its own canvas of exactly the size it occupies on the
assembled page and has to pass the layout gate (``pagelayout.assert_page_clean``:
no text collisions, no foreign-panel intrusion, no text on a spine, no grid, no
font below 7 pt) before its PDF is written.

Usage
-----
conda run -n benchmark-revision --no-capture-output python \
    $RNAMODBENCH_ROOT/analysis/figS9_rebuild/scripts/70_figs9_panels.py [A B ...]
"""
from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import s9_data as D           # noqa: E402
import s9_panels as P         # noqa: E402
import s9_style as S          # noqa: E402
from common import pagelayout  # noqa: E402
from common.panelpage import new_panel, save_panel  # noqa: E402

_lines: list[str] = []


def log(msg: str) -> None:
    line = f"[{datetime.now():%H:%M:%S}] {msg}"
    print(line, flush=True)
    _lines.append(line)


def assert_on_canvas(fig, key: str) -> None:
    """Abort when any drawn ink sticks out of the panel canvas.

    A title placed above an axis that already reaches the top of the canvas is
    silently clipped at export -- the defect class the Figure 8 revision hit with
    a lost x title.  ``get_tightbbox`` reports the real ink box in inches, so
    the check is arithmetic instead of a look at the preview.
    """
    fig.canvas.draw()
    box = fig.get_tightbbox(fig.canvas.get_renderer())
    over = []
    if box.x0 < -0.004:
        over.append(f"left {-box.x0:.3f} in")
    if box.y0 < -0.004:
        over.append(f"bottom {-box.y0:.3f} in")
    if box.x1 > fig.get_figwidth() + 0.004:
        over.append(f"right {box.x1 - fig.get_figwidth():.3f} in")
    if box.y1 > fig.get_figheight() + 0.004:
        over.append(f"top {box.y1 - fig.get_figheight():.3f} in")
    if over:
        raise SystemExit(f"panel {key}: ink leaves the canvas "
                         f"({', '.join(over)}) -- the print would clip it")


def column_left() -> dict[str, float]:
    """Label-column width shared by all nine panels of the 3 x 3 grid.

    Every cell is 4.027 x 3.480/3.550 in (2026-09-30) and the nine boxes begin on
    one vertical line, so one column serves them all.  Two needs set its width:
    A's 23 model names at ``FS["row_a"]`` (locked at 8 pt by their row pitch) and
    F's longest tool name at ``FS["row"]`` plus the clearance the panel letter
    needs.  Measured on 2026-09-30: A needs 1.856 in, F needs 1.861 in at 9.5 pt
    -- F bound the column that day by 0.005 in, so the nine boxes lost 0.005 in
    of width (accepted); the assertion below turns a runaway into a hard failure
    instead of a silent squeeze.
    """
    
    #: column, so every box starts and ends on the same two vertical lines and the
    #: page reads as one grid.  A's longest model name at FS["row_a"] (+0.09 in
    #: pad) vs F's longest tool name at FS["row"] (+0.32 in for the panel letter,
    #: 0.30 in of which is the 28 pt glyph itself).
    a_need = P.label_width_in(P.bar_labels(D.panel_a()), S.FS["row_a"]) + 0.09
    f_need = P.label_width_in(P.effect_labels(D.panel_effect()),
                              S.FS["row"]) + 0.32
    common = max(a_need, f_need)
    
    #: type change pushes it past this ceiling, the boxes would silently shrink;
    #: fail loudly and rebalance the scale instead
    assert common <= 1.90, (f"label column grew to {common:.4f} in "
                            f"(a_need {a_need:.4f}, f_need {f_need:.4f})")
    return {key: common for key in "ABCDEFGHI"}


def render(key: str, left_in: float) -> None:
    w, h = S.PIECE[key]
    S.apply()
    fig = new_panel(w, h)
    ignore: list = []

    if key == "A":
        P.draw_a(fig, D.panel_a(), left_in)
    elif key == "B":
        P.draw_orca(fig, D.panel_orca(), left_in)
    elif key == "C":
        P.draw_curlcake(fig, D.panel_curlcake(), left_in)
    elif key == "D":
        P.draw_window_ppv(fig, D.panel_window(), left_in)
        
        #: D draws the first half, E the second
        P.draw_sweep_key(fig, D.panel_window(), half="left", left_in=left_in)
    elif key == "E":
        P.draw_window_exact(fig, D.panel_window(), left_in)
        P.draw_sweep_key(fig, D.panel_window(), half="right", left_in=left_in)
    elif key == "F":
        P.draw_effect_r(fig, D.panel_effect(), left_in)
    elif key == "G":
        P.draw_effect_slope(fig, D.panel_effect(), left_in)
    elif key == "H":
        P.draw_stability_coverage(fig, D.panel_stability("H"), left_in)
    elif key == "I":
        P.draw_stability_ratio(fig, D.panel_stability("I"), left_in)
    else:
        raise SystemExit(f"unknown panel {key!r}")

    for ax in fig.axes:
        if ax.get_legend() is not None:
            try:
                pagelayout.assert_legend_clear(ax, verbose=False)
            except SystemExit as exc:
                log(f"panel {key}: {exc}")
    assert_on_canvas(fig, key)
    out = S.PANELS / f"figS9{key}.pdf"
    save_panel(fig, out, ignore_axes=list(ignore))
    log(f"panel {key}: {w:.3f} x {h:.3f} in -> {out.relative_to(S.PROJECT)}")


def main() -> None:
    keys = [a.upper() for a in sys.argv[1:]] or list("ABCDEFGHI")
    S.PANELS.mkdir(parents=True, exist_ok=True)
    S.LOGS.mkdir(parents=True, exist_ok=True)
    S.apply()
    cols = column_left()
    log("label columns: " + ", ".join(f"{k} {v:.3f} in" for k, v in cols.items()))
    for key in keys:
        render(key, cols[key])
    (S.TABLES / "S9_anchors.tsv").write_text(
        "panel\tquantity\tvalue\n"
        + "\n".join("\t".join(r) for r in D.anchors()) + "\n")
    log(f"dots drawn in total: {P.DOTS_DRAWN}")
    (S.LOGS / "70_panels.log").write_text("\n".join(_lines) + "\n")


if __name__ == "__main__":
    main()
