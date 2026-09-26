import numpy as np
from matplotlib.collections import PathCollection
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox


def qa_tick_overlaps(fig):
    """Pairwise overlap between the x tick labels of every axis.

    A 45 deg slant only stays readable while each label's horizontal footprint
    stays inside its slot, so the count of colliding pairs (and the worst
    overlap area) is the acceptance number for the rotation.
    """
    fig.canvas.draw()
    ren = fig.canvas.get_renderer()
    rows = []
    for i, ax in enumerate(fig.axes):
        labels = [t for t in ax.get_xticklabels()
                  if t.get_visible() and t.get_text()]
        boxes = [t.get_window_extent(renderer=ren) for t in labels]
        hits = []
        for a in range(len(boxes)):
            for b in range(a + 1, len(boxes)):
                ov = Bbox.intersection(boxes[a], boxes[b])
                if ov is not None:
                    hits.append((labels[a].get_text(), labels[b].get_text(),
                                 round(ov.width * ov.height, 1)))
        if hits:
            rows.append((i, len(hits), max(h[2] for h in hits), hits[:3]))
    print("[qa] x tick labels: axes, colliding pairs, worst overlap area (pt^2)")
    if not rows:
        print("  none")
    for i, n, worst, sample in rows:
        print(f"  ax{i:<2} pairs={n:<4} worst={worst:<7} eg {sample}")
    return rows


def qa_overlaps(fig, path=None):
    fig.canvas.draw()
    ren = fig.canvas.get_renderer()
    rows = []
    for i, ax in enumerate(fig.axes):
        leg = ax.get_legend()
        if leg is None:
            continue
        lb = leg.get_window_extent(renderer=ren)
        pts = []
        for c in ax.collections:
            if isinstance(c, PathCollection) and c.get_offsets().size:
                pts.append(c.get_offsets())
        for ln in ax.lines:
            #: 2026-09-24: coerce to arrays -- some backends hand back plain
            #: lists here, and .size then raises (fig2 G panel, hlines era)
            xd = np.asarray(ln.get_xdata(), dtype=float)
            yd = np.asarray(ln.get_ydata(), dtype=float)
            if xd.size > 1:
                pts.append(np.column_stack([xd, yd]))
        for pt in ax.patches:
            if isinstance(pt, Rectangle):
                pts.append(np.array([[pt.get_x(), pt.get_y()],
                                     [pt.get_x() + pt.get_width(),
                                      pt.get_y() + pt.get_height()]]))
        n = 0
        for p in pts:
            xy = ax.transData.transform(np.asarray(p, dtype=float))
            inside = ((xy[:, 0] > lb.x0) & (xy[:, 0] < lb.x1) &
                      (xy[:, 1] > lb.y0) & (xy[:, 1] < lb.y1))
            n += int(inside.sum())
        rows.append((i, n, round(100 * lb.width * lb.height /
                                 (ax.get_window_extent(renderer=ren).width *
                                  ax.get_window_extent(renderer=ren).height), 1)))
    print("[qa] legend: axes, data points covered, legend area % of panel")
    for i, n, frac in rows:
        flag = "  <-- OVERLAP" if n else ""
        print(f"  ax{i:<2} points_covered={n:<4} area={frac}%{flag}")
    return rows
