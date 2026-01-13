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
        write_cell_ids: bool = False,
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
            write_cell_ids: Whether to write per-bin cell id lists (off by default)
            compressor: Compression algorithm ("zstd" or "blosc")
            compression_level: Compression level
        """
        self.path = Path(path)
        self.zoom_levels = zoom_levels
        self.tile_size = tile_size
        self.n_genes = n_genes
        self.n_categories = n_categories
        self.numeric_columns = numeric_columns or []
        self.write_cell_ids_enabled = write_cell_ids
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

    def write_zoom_level_base(
        self,
        zoom: int,
        *,
        bin_coords: np.ndarray,
        cell_count: np.ndarray,
        expression_mean: np.ndarray,
        expression_max: np.ndarray,
    ) -> None:
        """Write non-category arrays for a zoom level."""
        if bin_coords.size == 0:
            logger.warning(f"No bins for zoom level {zoom}")
            return

        n_bins = bin_coords.shape[0]
        group = self.root[f"zoom_{zoom}"]

        group.create_array(
            "bin_coords",
            data=bin_coords.astype(np.int32, copy=False),
            chunks=(min(1000, n_bins), 2),
        )
        group.create_array(
            "cell_count",
            data=cell_count.astype(np.uint32, copy=False),
            chunks=(min(1000, n_bins),),
        )
        group.create_array(
            "expression_mean",
            data=expression_mean.astype(np.float32, copy=False),
            chunks=(min(1000, n_bins), min(100, self.n_genes)),
        )
        group.create_array(
            "expression_max",
            data=expression_max.astype(np.float32, copy=False),
            chunks=(min(1000, n_bins), min(100, self.n_genes)),
        )

        # Ensure groups exist for optional products
        group.create_group("category_counts")
        if self.numeric_columns:
            group.create_group("numeric_medians")

        # Store zoom-level metadata
        n_bins_per_axis = 2 ** zoom
        group.attrs["n_bins"] = n_bins
        group.attrs["n_bins_per_axis"] = n_bins_per_axis
        group.attrs["bin_size"] = 256.0 / n_bins_per_axis

        logger.info(f"Wrote base arrays for {n_bins:,} bins (zoom {zoom})")

    def write_category_counts_column(
        self,
        zoom: int,
        *,
        column: str,
        data: np.ndarray,
    ) -> None:
        """Write category counts for a single column at a zoom level."""
        group = self.root[f"zoom_{zoom}"]
        cat_group = group["category_counts"]
        cat_group.create_array(
            column,
            data=data.astype(np.uint32, copy=False),
            chunks=(min(1000, data.shape[0]), data.shape[1]),
        )

    def write_numeric_medians_column(
        self,
        zoom: int,
        *,
        column: str,
        data: np.ndarray,
    ) -> None:
        """Write numeric aggregate for a single column at a zoom level."""
        group = self.root[f"zoom_{zoom}"]
        num_group = group["numeric_medians"]
        num_group.create_array(
            column,
            data=data.astype(np.float32, copy=False),
            chunks=(min(1000, data.shape[0]),),
        )

    def write_cell_ids(
        self,
        zoom: int,
        *,
        cell_ids: list[np.ndarray],
    ) -> None:
        """Write per-bin cell ids (ragged) for a zoom level."""
        if not self.write_cell_ids_enabled:
            return
        if not cell_ids:
            logger.warning(f"No cell_ids for zoom level {zoom}")
            return

        group = self.root[f"zoom_{zoom}"]
        cell_ids_group = group.create_group("cell_ids")
        n_bins = len(cell_ids)

        chunk_size = 100
        for chunk_start in range(0, n_bins, chunk_size):
            chunk_end = min(chunk_start + chunk_size, n_bins)
            chunk_cell_ids = cell_ids[chunk_start:chunk_end]

            offsets = np.array([0] + [len(ids) for ids in chunk_cell_ids], dtype=np.uint64)
            offsets = np.cumsum(offsets)
            flat_ids = (
                np.concatenate(chunk_cell_ids) if chunk_cell_ids else np.array([], dtype=np.uint32)
            )

            chunk_group = cell_ids_group.create_group(f"chunk_{chunk_start}")
            chunk_group.create_array("offsets", data=offsets)
            chunk_group.create_array("cell_ids", data=flat_ids)

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

        self.write_zoom_level_base(
            zoom,
            bin_coords=bin_coords,
            cell_count=cell_counts,
            expression_mean=expression_mean,
            expression_max=expression_max,
        )

        # Write category counts
        for col_name, n_cats in self.n_categories.items():
            cat_data = np.array(
                [
                    b["category_counts"].get(col_name, np.zeros(n_cats, dtype=np.uint32))
                    for b in bins_data
                ],
                dtype=np.uint32,
            )
            self.write_category_counts_column(zoom, column=col_name, data=cat_data)

        # Write numeric medians
        if self.numeric_columns:
            for col_name in self.numeric_columns:
                num_data = np.array(
                    [b.get("numeric_medians", {}).get(col_name, np.nan) for b in bins_data],
                    dtype=np.float32,
                )
                self.write_numeric_medians_column(zoom, column=col_name, data=num_data)

        # Optional per-bin cell id lists
        if self.write_cell_ids_enabled:
            self.write_cell_ids(
                zoom,
                cell_ids=[b["cell_ids"] for b in bins_data],
            )

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
