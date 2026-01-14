"""Quadtree-based spatial binning for multi-resolution visualization."""

import numpy as np
from numba import jit, prange


class QuadtreeBinner:
    """Quadtree-based spatial binner for multi-resolution cell binning.

    Uses power-of-2 bin sizes for efficient quadtree subdivision:
    - Zoom 0: 256 units per bin (1 bin covers entire space)
    - Zoom 1: 128 units per bin (4 bins)
    - Zoom 2: 64 units per bin (16 bins)
    - ...
    - Zoom 7: 2 units per bin (16384 bins per axis)
    """

    def __init__(self, coordinate_range: float = 256.0, zoom_levels: int = 8):
        """Initialize the quadtree binner.

        Args:
            coordinate_range: Maximum coordinate value (coordinates in [0, coordinate_range))
            zoom_levels: Number of zoom levels to support (0 to zoom_levels-1)
        """
        self.coordinate_range = coordinate_range
        self.zoom_levels = zoom_levels

        # Precompute bin sizes for each zoom level
        # At zoom 0, bin_size = coordinate_range (1 bin)
        # At zoom z, bin_size = coordinate_range / 2^z
        self.bin_sizes = np.array([
            coordinate_range / (2 ** z) for z in range(zoom_levels)
        ])

    def get_bin_size(self, zoom: int) -> float:
        """Get bin size for a given zoom level."""
        if zoom < 0 or zoom >= self.zoom_levels:
            raise ValueError(f"Zoom {zoom} out of range [0, {self.zoom_levels})")
        return self.bin_sizes[zoom]

    def get_n_bins_per_axis(self, zoom: int) -> int:
        """Get number of bins per axis at a given zoom level."""
        return 2 ** zoom

    def assign_bins(self, coords: np.ndarray, zoom: int) -> np.ndarray:
        """Assign cells to bins at a given zoom level.

        Args:
            coords: Cell coordinates, shape (n_cells, 2)
            zoom: Zoom level

        Returns:
            Bin assignments, shape (n_cells, 2) with (bin_x, bin_y) per cell
        """
        bin_size = self.get_bin_size(zoom)
        n_bins = self.get_n_bins_per_axis(zoom)

        # Compute bin indices
        bin_assignments = np.floor(coords / bin_size).astype(np.int32)

        # Clamp to valid range [0, n_bins - 1]
        bin_assignments = np.clip(bin_assignments, 0, n_bins - 1)

        return bin_assignments

    def get_bin_bounds(self, bin_x: int, bin_y: int, zoom: int) -> tuple[float, float, float, float]:
        """Get coordinate bounds for a bin.

        Args:
            bin_x: Bin X index
            bin_y: Bin Y index
            zoom: Zoom level

        Returns:
            Tuple of (min_x, min_y, max_x, max_y)
        """
        bin_size = self.get_bin_size(zoom)
        min_x = bin_x * bin_size
        min_y = bin_y * bin_size
        max_x = min_x + bin_size
        max_y = min_y + bin_size
        return (min_x, min_y, max_x, max_y)

    def get_bin_center(self, bin_x: int, bin_y: int, zoom: int) -> tuple[float, float]:
        """Get center coordinates of a bin.

        Args:
            bin_x: Bin X index
            bin_y: Bin Y index
            zoom: Zoom level

        Returns:
            Tuple of (center_x, center_y)
        """
        min_x, min_y, max_x, max_y = self.get_bin_bounds(bin_x, bin_y, zoom)
        return ((min_x + max_x) / 2, (min_y + max_y) / 2)

    def get_tile_bins(self, tile_x: int, tile_y: int, zoom: int, tile_size: int = 256) -> list[tuple[int, int]]:
        """Get all bin indices that fall within a tile.

        Args:
            tile_x: Tile X index
            tile_y: Tile Y index
            zoom: Zoom level
            tile_size: Tile size in pixels

        Returns:
            List of (bin_x, bin_y) tuples for bins in this tile
        """
        # At each zoom level, there's one bin per pixel-equivalent in the tile
        # For zoom z, each tile covers tile_size bins (since we align bins to tiles)
        bins_per_tile = self.get_n_bins_per_axis(zoom) // max(1, self.get_n_bins_per_axis(zoom) // tile_size)

        # But for simplicity, we just return all bins whose center falls in this tile
        bins = []
        n_bins = self.get_n_bins_per_axis(zoom)
        bin_size = self.get_bin_size(zoom)

        # Tile bounds in coordinate space
        tile_coord_size = self.coordinate_range / (2 ** zoom) * tile_size / 256
        tile_min_x = tile_x * tile_coord_size
        tile_max_x = tile_min_x + tile_coord_size
        tile_min_y = tile_y * tile_coord_size
        tile_max_y = tile_min_y + tile_coord_size

        # Find bins that overlap with tile
        min_bin_x = int(tile_min_x / bin_size)
        max_bin_x = min(n_bins - 1, int(tile_max_x / bin_size))
        min_bin_y = int(tile_min_y / bin_size)
        max_bin_y = min(n_bins - 1, int(tile_max_y / bin_size))

        for bx in range(min_bin_x, max_bin_x + 1):
            for by in range(min_bin_y, max_bin_y + 1):
                bins.append((bx, by))

        return bins


@jit(nopython=True, parallel=True, cache=True)
def _fast_bin_assignment(coords: np.ndarray, bin_size: float, n_bins: int) -> np.ndarray:
    """Numba-accelerated bin assignment.

    Args:
        coords: Cell coordinates, shape (n_cells, 2)
        bin_size: Size of each bin
        n_bins: Number of bins per axis

    Returns:
        Bin assignments, shape (n_cells, 2)
    """
    n_cells = coords.shape[0]
    result = np.empty((n_cells, 2), dtype=np.int32)

    for i in prange(n_cells):
        bx = int(coords[i, 0] / bin_size)
        by = int(coords[i, 1] / bin_size)

        # Clamp to valid range
        result[i, 0] = max(0, min(n_bins - 1, bx))
        result[i, 1] = max(0, min(n_bins - 1, by))

    return result
