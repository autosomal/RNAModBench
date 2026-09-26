"""Confusion-matrix metrics and bootstrap confidence intervals.

Definitions (shared by every evaluation script; written into the README):

For a candidate universe ``U`` (see ``04_build_universe.py``), a reference site
set ``G`` (GLORI for m6A; Curlcake ground truth for the synthetic controls) and
a tool call set ``C``, all restricted to ``U``:

* ``TP(w) = |{u in U : u in C and dist(u, G) <= w}|``
* ``FP(w) = |{u in U : u in C and dist(u, G) >  w}|``
* ``FN(w) = |{u in U : u not in C and dist(u, G) <= w}|``  (reference sites that
  were testable but not called)
* ``TN(w) = |{u in U : u not in C and dist(u, G) >  w}|``

``precision = TP/(TP+FP)`` (this is the manuscript's "GLORI hit rate"),
``recall = TP/(TP+FN)``, ``F1``, ``MCC`` and ``specificity = TN/(TN+FP)`` follow
the usual formulas.  Calls outside ``U`` are never TP/FP; they are reported
separately as ``n_calls_out_of_universe``.

Bootstrap CIs resample *sites* (B = 1000, fixed seed) and are reported for
precision and recall.
"""

from __future__ import annotations

import math

import numpy as np


def window_hit_mask(query: np.ndarray, reference_sorted: np.ndarray, w: int) -> np.ndarray:
    """Does each query position have a reference position within +/- w?"""
    if reference_sorted.size == 0 or query.size == 0:
        return np.zeros(query.shape, dtype=bool)
    left = np.searchsorted(reference_sorted, query - w, side="left")
    right = np.searchsorted(reference_sorted, query + w, side="right")
    return right > left


def confusion_from_universe(call_pos: np.ndarray, glori_pos: np.ndarray,
                            universe_pos: np.ndarray, w: int) -> dict:
    """Full TP/FP/FN/TN against an explicit universe position array."""
    tp, fp, fn, _ = _counts(call_pos, glori_pos, universe_pos, w)
    tn = int(universe_pos.size) - tp - fp - fn
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "n_calls": int(call_pos.size), "n_glori": int(glori_pos.size),
            "universe": int(universe_pos.size)}


def _counts(call_pos: np.ndarray, glori_pos: np.ndarray,
            universe_pos: np.ndarray, w: int) -> tuple[int, int, int, int]:
    """TP/FP/FN/(unused) on sorted position arrays; all inputs must be subsets of U."""
    n_calls = call_pos.size
    hit_call = window_hit_mask(call_pos, glori_pos, w)
    tp = int(hit_call.sum())
    fp = int(n_calls - tp)
    glori_hit = window_hit_mask(glori_pos, call_pos, w)
    fn = int(glori_pos.size - int(glori_hit.sum()))
    return tp, fp, fn, 0


def precision(tp: int, fp: int) -> float:
    return tp / (tp + fp) if (tp + fp) > 0 else float("nan")


def recall(tp: int, fn: int) -> float:
    return tp / (tp + fn) if (tp + fn) > 0 else float("nan")


def f1(p: float, r: float) -> float:
    if not np.isfinite(p) or not np.isfinite(r) or (p + r) == 0:
        return float("nan")
    return 2.0 * p * r / (p + r)


def specificity(tn: int, fp: int) -> float:
    return tn / (tn + fp) if (tn + fp) > 0 else float("nan")


def mcc(tp: int, fp: int, fn: int, tn: int) -> float:
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if den == 0:
        return float("nan")
    return (tp * tn - fp * fn) / den


def bootstrap_ci(values: np.ndarray, n_bootstrap: int = 1000, confidence: float = 0.95,
                 seed: int = 20260914) -> tuple[float, float]:
    """Percentile bootstrap CI of the mean for a 0/1 vector (site-level)."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    n = values.size
    means = np.empty(n_bootstrap)
    for i in range(n_bootstrap):
        means[i] = values[rng.integers(0, n, n)].mean()
    alpha = 1.0 - confidence
    return float(np.percentile(means, 100 * alpha / 2)), float(
        np.percentile(means, 100 * (1 - alpha / 2)))


def precision_recall_ci(tp: int, fp: int, fn: int, n_bootstrap: int = 1000,
                        seed: int = 20260914) -> dict:
    """Site-level bootstrap CIs for precision and recall (0/1 vectors)."""
    prec_values = np.concatenate([np.ones(tp), np.zeros(fp)]) if (tp + fp) else np.array([])
    rec_values = np.concatenate([np.ones(tp), np.zeros(fn)]) if (tp + fn) else np.array([])
    p_lo, p_hi = bootstrap_ci(prec_values, n_bootstrap, seed=seed)
    r_lo, r_hi = bootstrap_ci(rec_values, n_bootstrap, seed=seed + 1)
    return {"precision_ci_lo": p_lo, "precision_ci_hi": p_hi,
            "recall_ci_lo": r_lo, "recall_ci_hi": r_hi}


def summarize(tp: int, fp: int, fn: int, tn: int) -> dict:
    p = precision(tp, fp)
    r = recall(tp, fn)
    return {"precision": p, "recall": r, "f1": f1(p, r),
            "specificity": specificity(tn, fp), "mcc": mcc(tp, fp, fn, tn)}
