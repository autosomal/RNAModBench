"""Logging + file inventory for the harmonisation rebuild.

Every script opens a :class:`RunLog` (stdout + file) and records produced files
in the shared ``harmonisation/manifest/file_inventory.csv``.  WARN/ERROR lines are
never silent: the log path is printed at the end of every run.
"""

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .config import LOG_DIR, MANIFEST_DIR
from .io_utils import file_fingerprint, now_stamp, rel

INVENTORY = MANIFEST_DIR / "file_inventory.csv"
_INVENTORY_COLUMNS = ["file", "rows", "bytes", "fingerprint", "source", "generated_by",
                      "generated_at"]


def setup_logger(name: str, *, log_dir: Path | None = None,
                 stamp: str | None = None) -> logging.Logger:
    """Logger writing to stdout and ``logs/<name>_<stamp>.log``."""
    log_dir = log_dir or LOG_DIR
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = stamp or now_stamp().replace("-", "").replace(":", "").replace(" ", "_")
    path = log_dir / f"{name}_{stamp}.log"

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%H:%M:%S")

    fh = logging.FileHandler(path, encoding="utf-8")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    logger.log_path = str(path)  # type: ignore[attr-defined]
    return logger


@dataclass
class Inventory:
    """Accumulates produced files and appends them to file_inventory.csv."""

    script: str
    rows: list[dict] = field(default_factory=list)

    def record(self, path: Path, *, n_rows: int | None = None,
               source: str = "") -> None:
        try:
            n_bytes = Path(path).stat().st_size
        except FileNotFoundError:
            return
        self.rows.append({
            "file": rel(Path(path)),
            "rows": "" if n_rows is None else int(n_rows),
            "bytes": int(n_bytes),
            "fingerprint": file_fingerprint(Path(path)),
            "source": source,
            "generated_by": self.script,
            "generated_at": now_stamp(),
        })

    def flush(self) -> None:
        if not self.rows:
            return
        df = pd.DataFrame(self.rows, columns=_INVENTORY_COLUMNS)
        MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
        if INVENTORY.exists():
            old = pd.read_csv(INVENTORY, sep="\t", dtype=str, keep_default_na=False)
            # replace previous records of the same script + file
            key = old["generated_by"].astype(str) + "|" + old["file"].astype(str)
            new_key = df["generated_by"].astype(str) + "|" + df["file"].astype(str)
            old = old[~key.isin(set(new_key))]
            df = pd.concat([old, df], ignore_index=True)
        df.to_csv(INVENTORY, sep="\t", index=False)


def log_time(logger: logging.Logger, label: str):
    """Context manager logging elapsed seconds (tiny helper, no deps)."""
    import contextlib
    import time

    @contextlib.contextmanager
    def _cm():
        t0 = time.time()
        yield
        logger.info("%s done in %.1fs", label, time.time() - t0)

    return _cm()
