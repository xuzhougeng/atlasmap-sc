"""Configuration schema for preprocessing pipeline."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import yaml


@dataclass
class PreprocessConfig:
    """Configuration for the preprocessing pipeline."""

    # Input/Output
    input_path: Path
    output_dir: Path

    # Coordinate settings
    # NOTE: `umap_key` is kept for backward compatibility; it now represents the
    # default coordinate system key (not necessarily UMAP).
    umap_key: str = "X_umap"
    # Optional list of coordinate keys to preprocess (supports multiple).
    # If empty, falls back to [`umap_key`].
    coordinate_keys: list[str] = field(default_factory=list)
    coordinate_range: float = 256.0  # Normalize to [0, coordinate_range)

    # Binning settings
    zoom_levels: int = 12  # Recommend 12 for clear tiles at low zoom levels
    tile_size: int = 256  # Pixels per tile

    # Gene selection
    n_genes: int = 500  # Number of genes to pre-aggregate
    hvg_n_top: int = 300  # Top highly variable genes
    marker_genes: list[str] = field(default_factory=list)
    use_all_expressed: bool = False  # If True, use all genes with expression > 0
    min_cells_expressed: int = 3  # Minimum cells where gene must be expressed

    # Category settings
    category_columns: list[str] = field(default_factory=list)  # obs columns to include
    default_category: str = "cell_type"
    exclude_category_columns: list[str] = field(default_factory=list)  # 显式排除的类别列
    max_category_cardinality: int = 2000  # 超过此唯一值数的列自动跳过（避免内存爆炸）
    max_category_matrix_elems: int = 50_000_000  # n_bins * n_cats 的上限

    # Numeric column settings
    numeric_columns: list[str] = field(default_factory=list)  # 显式指定的数值列
    exclude_numeric_columns: list[str] = field(default_factory=list)  # 排除自动检测的列

    # Zarr settings
    zarr_compressor: str = "zstd"
    zarr_compression_level: int = 3
    chunk_size: int = 256
    write_cell_ids: bool = False  # Optional: store per-bin cell id lists (expensive)

    # Memory settings
    batch_size: int = 100000  # Cells to process at once

    # SOMA storage settings
    enable_soma: bool = True  # Enable TileDBSOMA storage for all genes
    soma_compression: str = "zstd"
    soma_compression_level: int = 9

    # Metadata
    dataset_name: Optional[str] = None
    dataset_description: Optional[str] = None

    @classmethod
    def from_yaml(cls, path: Path) -> "PreprocessConfig":
        """Load configuration from a YAML file."""
        with open(path) as f:
            data = yaml.safe_load(f)

        # Convert path strings to Path objects
        if "input_path" in data:
            data["input_path"] = Path(data["input_path"])
        if "output_dir" in data:
            data["output_dir"] = Path(data["output_dir"])

        return cls(**data)

    def to_yaml(self, path: Path) -> None:
        """Save configuration to a YAML file."""
        data = {
            "input_path": str(self.input_path),
            "output_dir": str(self.output_dir),
            "umap_key": self.umap_key,
            "coordinate_keys": self.coordinate_keys or [self.umap_key],
            "coordinate_range": self.coordinate_range,
            "zoom_levels": self.zoom_levels,
            "tile_size": self.tile_size,
            "n_genes": self.n_genes,
            "hvg_n_top": self.hvg_n_top,
            "marker_genes": self.marker_genes,
            "use_all_expressed": self.use_all_expressed,
            "min_cells_expressed": self.min_cells_expressed,
            "category_columns": self.category_columns,
            "default_category": self.default_category,
            "exclude_category_columns": self.exclude_category_columns,
            "max_category_cardinality": self.max_category_cardinality,
            "max_category_matrix_elems": self.max_category_matrix_elems,
            "numeric_columns": self.numeric_columns,
            "exclude_numeric_columns": self.exclude_numeric_columns,
            "zarr_compressor": self.zarr_compressor,
            "zarr_compression_level": self.zarr_compression_level,
            "chunk_size": self.chunk_size,
            "write_cell_ids": self.write_cell_ids,
            "batch_size": self.batch_size,
            "enable_soma": self.enable_soma,
            "soma_compression": self.soma_compression,
            "soma_compression_level": self.soma_compression_level,
            "dataset_name": self.dataset_name,
            "dataset_description": self.dataset_description,
        }
        with open(path, "w") as f:
            yaml.dump(data, f, default_flow_style=False)

    def validate(self) -> None:
        """Validate configuration parameters."""
        if not self.input_path.exists():
            raise FileNotFoundError(f"Input file not found: {self.input_path}")

        # Coordinate keys
        keys = self.coordinate_keys if self.coordinate_keys else [self.umap_key]
        keys = [str(k).strip() for k in keys if str(k).strip()]
        if not keys:
            raise ValueError("No coordinate keys provided (set coordinate_keys or umap_key)")
        # De-duplicate while preserving order
        seen: set[str] = set()
        uniq: list[str] = []
        for k in keys:
            if k in seen:
                continue
            seen.add(k)
            uniq.append(k)
        self.coordinate_keys = uniq
        # Ensure default key exists in the list
        if self.umap_key not in self.coordinate_keys:
            self.umap_key = self.coordinate_keys[0]

        if self.zoom_levels < 1 or self.zoom_levels > 12:
            raise ValueError(f"zoom_levels must be between 1 and 12, got {self.zoom_levels}")

        if not self.use_all_expressed and self.n_genes < 1:
            raise ValueError(f"n_genes must be positive, got {self.n_genes}")

        if self.min_cells_expressed < 1:
            raise ValueError(f"min_cells_expressed must be >= 1, got {self.min_cells_expressed}")

        if self.tile_size not in [128, 256, 512]:
            raise ValueError(f"tile_size must be 128, 256, or 512, got {self.tile_size}")
