"""Zarr writer for multi-resolution bin data."""

import logging
from pathlib import Path
from typing import Any

import numpy as np
import zarr

logger = logging.getLogger(__name__)


class ZarrBinWriter:
    """Writer for multi-resolution bin data in Zarr format."""

    def __init__(
        self,
        path: Path,
        zoom_levels: int,
        tile_size: int,
        n_genes: int,
        n_categories: dict[str, int],
        numeric_columns: list[str] | None = None,
        compressor: str = "zstd",
        compression_level: int = 3,
    ):
        """Initialize the Zarr writer.

        Args:
            path: Output path for Zarr store
            zoom_levels: Number of zoom levels
            tile_size: Tile size in pixels
            n_genes: Number of pre-aggregated genes
            n_categories: Dict mapping category column names to number of categories
            numeric_columns: List of numeric column names
            compressor: Compression algorithm ("zstd" or "blosc")
            compression_level: Compression level
        """
        self.path = Path(path)
        self.zoom_levels = zoom_levels
        self.tile_size = tile_size
        self.n_genes = n_genes
        self.n_categories = n_categories
        self.numeric_columns = numeric_columns or []
        self.compression_level = compression_level

        # Create root directory
        self.path.mkdir(parents=True, exist_ok=True)

        # Create root group using zarr v3 API
        self.root = zarr.open_group(str(self.path), mode="w")

        # Add metadata
        self.root.attrs["zoom_levels"] = zoom_levels
        self.root.attrs["tile_size"] = tile_size
        self.root.attrs["n_genes"] = n_genes
        self.root.attrs["n_categories"] = dict(n_categories)
        self.root.attrs["numeric_columns"] = self.numeric_columns
        self.root.attrs["format_version"] = "1.0"

        # Create groups for each zoom level
        for z in range(zoom_levels):
            self.root.create_group(f"zoom_{z}")

    def write_zoom_level(self, zoom: int, bins_data: list[dict[str, Any]]) -> None:
        """Write bin data for a single zoom level.

        Args:
            zoom: Zoom level
            bins_data: List of bin data dictionaries with keys:
                - bin_x: int
                - bin_y: int
                - cell_count: int
                - cell_ids: np.ndarray of uint32
                - expression_mean: np.ndarray of float32
                - expression_max: np.ndarray of float32
                - category_counts: dict[str, np.ndarray]
                - numeric_medians: dict[str, float] (optional)
        """
        if not bins_data:
            logger.warning(f"No bins for zoom level {zoom}")
            return

        n_bins = len(bins_data)
        group = self.root[f"zoom_{zoom}"]

        # Prepare arrays
        bin_coords = np.array(
            [(b["bin_x"], b["bin_y"]) for b in bins_data],
            dtype=np.int32,
        )
        cell_counts = np.array(
            [b["cell_count"] for b in bins_data],
            dtype=np.uint32,
        )
        expression_mean = np.array(
            [b["expression_mean"] for b in bins_data],
            dtype=np.float32,
        )
        expression_max = np.array(
            [b["expression_max"] for b in bins_data],
            dtype=np.float32,
        )

        # Write coordinate and count arrays using zarr v3 API
        group.create_array(
            "bin_coords",
            data=bin_coords,
            chunks=(min(1000, n_bins), 2),
        )
        group.create_array(
            "cell_count",
            data=cell_counts,
            chunks=(min(1000, n_bins),),
        )

        # Write expression arrays
        group.create_array(
            "expression_mean",
            data=expression_mean,
            chunks=(min(1000, n_bins), min(100, self.n_genes)),
        )
        group.create_array(
            "expression_max",
            data=expression_max,
            chunks=(min(1000, n_bins), min(100, self.n_genes)),
        )

        # Write category counts
        cat_group = group.create_group("category_counts")
        for col_name, n_cats in self.n_categories.items():
            cat_data = np.array(
                [b["category_counts"].get(col_name, np.zeros(n_cats, dtype=np.uint32))
                 for b in bins_data],
                dtype=np.uint32,
            )
            cat_group.create_array(
                col_name,
                data=cat_data,
                chunks=(min(1000, n_bins), n_cats),
            )

        # Write numeric medians
        if self.numeric_columns:
            num_group = group.create_group("numeric_medians")
            for col_name in self.numeric_columns:
                num_data = np.array(
                    [b.get("numeric_medians", {}).get(col_name, np.nan)
                     for b in bins_data],
                    dtype=np.float32,
                )
                num_group.create_array(
                    col_name,
                    data=num_data,
                    chunks=(min(1000, n_bins),),
                )

        # Write cell IDs as variable-length data
        # Store as object array with ragged arrays
        cell_ids_group = group.create_group("cell_ids")

        # For memory efficiency, store cell IDs in chunks
        # Each chunk contains cell IDs for a range of bins
        chunk_size = 100
        for chunk_start in range(0, n_bins, chunk_size):
            chunk_end = min(chunk_start + chunk_size, n_bins)
            chunk_cell_ids = [bins_data[i]["cell_ids"] for i in range(chunk_start, chunk_end)]

            # Flatten and store with offsets
            offsets = np.array([0] + [len(ids) for ids in chunk_cell_ids], dtype=np.uint64)
            offsets = np.cumsum(offsets)
            flat_ids = np.concatenate(chunk_cell_ids) if chunk_cell_ids else np.array([], dtype=np.uint32)

            chunk_group = cell_ids_group.create_group(f"chunk_{chunk_start}")
            chunk_group.create_array("offsets", data=offsets)
            chunk_group.create_array("cell_ids", data=flat_ids)

        # Store zoom-level metadata
        n_bins_per_axis = 2 ** zoom
        group.attrs["n_bins"] = n_bins
        group.attrs["n_bins_per_axis"] = n_bins_per_axis
        group.attrs["bin_size"] = 256.0 / n_bins_per_axis

        logger.info(f"Wrote {n_bins:,} bins for zoom level {zoom}")

    def finalize(self) -> None:
        """Finalize the Zarr store."""
        logger.info(f"Finalized Zarr store at: {self.path}")

    def create_spatial_index(self, zoom: int) -> None:
        """Create a spatial index for fast tile-based queries.

        Args:
            zoom: Zoom level to index
        """
        group = self.root[f"zoom_{zoom}"]
        bin_coords = group["bin_coords"][:]

        # Create a mapping from tile coordinates to bin indices
        # This enables fast lookup of bins within a tile
        n_bins_per_axis = 2 ** zoom
        bins_per_tile = max(1, n_bins_per_axis // (256 // self.tile_size))

        tile_to_bins: dict[tuple[int, int], list[int]] = {}
        for bin_idx, (bin_x, bin_y) in enumerate(bin_coords):
            tile_x = bin_x // bins_per_tile
            tile_y = bin_y // bins_per_tile
            key = (int(tile_x), int(tile_y))
            if key not in tile_to_bins:
                tile_to_bins[key] = []
            tile_to_bins[key].append(bin_idx)

        # Store as JSON-compatible structure in attrs
        group.attrs["tile_index"] = {
            f"{tx},{ty}": indices for (tx, ty), indices in tile_to_bins.items()
        }
