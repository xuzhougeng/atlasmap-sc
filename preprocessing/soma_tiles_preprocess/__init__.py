"""SOMA-Tiles Preprocessing Pipeline.

Convert H5AD single-cell data to multi-resolution Zarr bins for high-performance visualization.
"""

__version__ = "0.1.0"

from .pipeline import PreprocessingPipeline
from .config import PreprocessConfig

__all__ = ["PreprocessingPipeline", "PreprocessConfig", "__version__"]
