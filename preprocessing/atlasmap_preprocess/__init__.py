"""AtlasMap Preprocessing Pipeline.

Convert H5AD single-cell data to multi-resolution Zarr bins for high-performance visualization.
"""

from __future__ import annotations

import os
from pathlib import Path

__version__ = "0.2.0"

# Ensure Numba has a writable cache directory. Some environments (e.g. read-only
# site-packages) can cause Scanpy's `@njit(cache=True)` imports to fail unless
# NUMBA_CACHE_DIR is set.
if "NUMBA_CACHE_DIR" not in os.environ:
    candidates: list[Path] = []
    xdg_cache = os.environ.get("XDG_CACHE_HOME")
    if xdg_cache:
        candidates.append(Path(xdg_cache) / "numba")
    candidates.append(Path.home() / ".cache" / "numba")
    candidates.append(Path.cwd() / ".numba_cache")

    for candidate in candidates:
        try:
            candidate.mkdir(parents=True, exist_ok=True)
        except Exception:
            continue
        os.environ["NUMBA_CACHE_DIR"] = str(candidate)
        break

from .pipeline import PreprocessingPipeline
from .config import PreprocessConfig

__all__ = ["PreprocessingPipeline", "PreprocessConfig", "__version__"]
