"""Coordinate normalization utilities."""

import numpy as np


class CoordinateNormalizer:
    """Normalizes coordinates to a fixed range for consistent binning."""

    def __init__(self, target_range: float = 256.0, padding: float = 0.02):
        """Initialize the normalizer.

        Args:
            target_range: Target coordinate range [0, target_range)
            padding: Padding fraction to add to bounds (to avoid edge cases)
        """
        self.target_range = target_range
        self.padding = padding

        # Will be set during normalization
        self.original_min: np.ndarray | None = None
        self.original_max: np.ndarray | None = None
        self.scale: float | None = None

    def normalize(self, coords: np.ndarray) -> np.ndarray:
        """Normalize coordinates to [0, target_range).

        Args:
            coords: Input coordinates, shape (n_cells, 2)

        Returns:
            Normalized coordinates, shape (n_cells, 2)
        """
        # Compute bounds
        self.original_min = coords.min(axis=0)
        self.original_max = coords.max(axis=0)

        # Add padding
        range_vals = self.original_max - self.original_min
        pad = range_vals * self.padding
        self.original_min = self.original_min - pad
        self.original_max = self.original_max + pad

        # Use the larger dimension to maintain aspect ratio
        range_vals = self.original_max - self.original_min
        max_range = max(range_vals[0], range_vals[1])

        # Scale to [0, target_range)
        self.scale = (self.target_range - 0.001) / max_range  # Slight offset to stay in range

        # Center the smaller dimension
        centered_min = self.original_min.copy()
        if range_vals[0] < range_vals[1]:
            offset = (max_range - range_vals[0]) / 2
            centered_min[0] -= offset
        elif range_vals[1] < range_vals[0]:
            offset = (max_range - range_vals[1]) / 2
            centered_min[1] -= offset

        # Normalize
        normalized = (coords - centered_min) * self.scale

        # Clip to valid range
        normalized = np.clip(normalized, 0, self.target_range - 0.001)

        return normalized.astype(np.float64)

    def denormalize(self, coords: np.ndarray) -> np.ndarray:
        """Convert normalized coordinates back to original space.

        Args:
            coords: Normalized coordinates, shape (n_cells, 2)

        Returns:
            Original coordinates, shape (n_cells, 2)
        """
        if self.scale is None:
            raise ValueError("Must call normalize() before denormalize()")

        # Reverse the normalization
        range_vals = self.original_max - self.original_min
        max_range = max(range_vals[0], range_vals[1])

        centered_min = self.original_min.copy()
        if range_vals[0] < range_vals[1]:
            offset = (max_range - range_vals[0]) / 2
            centered_min[0] -= offset
        elif range_vals[1] < range_vals[0]:
            offset = (max_range - range_vals[1]) / 2
            centered_min[1] -= offset

        return (coords / self.scale) + centered_min

    def get_transform_params(self) -> dict:
        """Get transformation parameters for metadata storage.

        Returns:
            Dictionary with transformation parameters
        """
        if self.scale is None:
            raise ValueError("Must call normalize() before get_transform_params()")

        return {
            "original_min": self.original_min.tolist(),
            "original_max": self.original_max.tolist(),
            "scale": float(self.scale),
            "target_range": float(self.target_range),
            "padding": float(self.padding),
        }
