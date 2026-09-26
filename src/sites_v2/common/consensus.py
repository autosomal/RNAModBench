"""Consensus rule shared by every replicate-aware revision analysis.

``quorum(n) = n // 2 + 1`` -- a site must be called by *strictly more than half*
of the group's independent units.  ``ceil(n/2)`` would make the "majority" of a
two-unit group identical to its union, which is exactly the collapse the
revision is supposed to avoid, so the strict form is used everywhere instead.
"""

from __future__ import annotations


def quorum(n_units: int) -> int:
    return n_units // 2 + 1
