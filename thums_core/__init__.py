"""THUMS — thermal storage in repurposed wells with phase-change materials."""
from __future__ import annotations
import datetime as _dt
from . import errors, kpi, doe, results          # noqa: F401

__version__ = "0.2.0"

def _stamp():
    c = results.git_commit()
    return f"THUMS v{__version__} · commit {c}"

VERSION_STAMP = _stamp()
