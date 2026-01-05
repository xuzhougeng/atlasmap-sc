"""Main preprocessing pipeline for SOMA-Tiles."""

import logging
from pathlib import Path
from typing import Optional
import json

import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm

from .config import PreprocessConfig
from .binning.quadtree import QuadtreeBinner
from .binning.aggregator import ExpressionAggregator
from .binning.normalizer import CoordinateNormalizer
from .io.zarr_writer import ZarrBinWriter
from .io.soma_writer import SomaWriter
from .genes.selector import GeneSelector

logger = logging.getLogger(__name__)


class PreprocessingPipeline:
    """Main preprocessing pipeline for converting H5AD to multi-resolution Zarr bins."""

    def __init__(self, config: PreprocessConfig):
        """Initialize the preprocessing pipeline.

        Args:
            config: Preprocessing configuration
        """
        self.config = config
        self.adata: Optional[sc.AnnData] = None
        self.normalized_coords: Optional[np.ndarray] = None
        self.selected_genes: Optional[list[str]] = None
        self.category_mapping: dict[str, dict[str, int]] = {}
        self.numeric_columns: list[str] = []  # 存储检测到的数值列
        self.column_name_mapping: dict[str, str] = {}  # 原始列名 -> 清理后列名

    def run(self) -> None:
        """Execute the full preprocessing pipeline."""
        logger.info("Starting preprocessing pipeline")

        # Step 1: Load H5AD
        self._load_data()

        # Step 2: Extract and normalize coordinates
        self._process_coordinates()

        # Step 3: Select genes for pre-aggregation
        self._select_genes()

        # Step 4: Build category mappings
        self._build_category_mappings()

        # Step 5: Build multi-resolution bins and write to Zarr
        self._build_and_write_bins()

        # Step 6: Write SOMA store (if enabled)
        if self.config.enable_soma:
            self._write_soma_store()

        # Step 7: Write metadata
        self._write_metadata()

        logger.info("Pipeline completed successfully")

    def _load_data(self) -> None:
        """Load H5AD file."""
        logger.info(f"Loading H5AD from: {self.config.input_path}")
        self.adata = sc.read_h5ad(self.config.input_path)
        logger.info(f"Loaded {self.adata.n_obs:,} cells, {self.adata.n_vars:,} genes")

    def _process_coordinates(self) -> None:
        """Extract and normalize UMAP coordinates."""
        logger.info(f"Extracting coordinates from obsm['{self.config.umap_key}']")

        if self.config.umap_key not in self.adata.obsm:
            available_keys = list(self.adata.obsm.keys())
            raise KeyError(
                f"UMAP key '{self.config.umap_key}' not found. "
                f"Available keys: {available_keys}"
            )

        coords = self.adata.obsm[self.config.umap_key][:, :2].astype(np.float64)

        # Normalize to [0, coordinate_range)
        normalizer = CoordinateNormalizer(self.config.coordinate_range)
        self.normalized_coords = normalizer.normalize(coords)

        logger.info(
            f"Normalized coordinates to [{0}, {self.config.coordinate_range}), "
            f"shape: {self.normalized_coords.shape}"
        )

    def _select_genes(self) -> None:
        """Select genes for pre-aggregation."""
        logger.info("Selecting genes for pre-aggregation")

        selector = GeneSelector(self.adata)
        self.selected_genes = selector.select(
            n_genes=self.config.n_genes,
            hvg_n_top=self.config.hvg_n_top,
            marker_genes=self.config.marker_genes,
            use_all_expressed=self.config.use_all_expressed,
            min_cells_expressed=self.config.min_cells_expressed,
        )

        logger.info(f"Selected {len(self.selected_genes)} genes for pre-aggregation")

    def _sanitize_column_name(self, col: str) -> str:
        """清理列名，替换不适合 URL 的字符。

        Args:
            col: 原始列名

        Returns:
            清理后的列名（. 替换为 _）
        """
        return col.replace(".", "_")

    def _is_numeric_column(self, col: str) -> bool:
        """判断列是否为数值型。

        Args:
            col: obs 列名

        Returns:
            True 如果是数值型列
        """
        # 显式排除
        if col in self.config.exclude_numeric_columns:
            return False
        # 显式指定
        if col in self.config.numeric_columns:
            return True
        # 自动检测
        return pd.api.types.is_numeric_dtype(self.adata.obs[col])

    def _build_category_mappings(self) -> None:
        """Build mappings from category names to integer indices."""
        logger.info("Building category mappings")

        # Determine which columns to include
        columns = self.config.category_columns.copy()
        if not columns and self.config.default_category in self.adata.obs.columns:
            columns = [self.config.default_category]

        for col in columns:
            if col not in self.adata.obs.columns:
                logger.warning(f"Category column '{col}' not found in obs, skipping")
                continue

            # 清理列名（替换不适合 URL 的字符）
            sanitized_col = self._sanitize_column_name(col)
            self.column_name_mapping[sanitized_col] = col  # 清理后 -> 原始

            if self._is_numeric_column(col):
                # 数值列：跳过分类映射
                self.numeric_columns.append(sanitized_col)
                if sanitized_col != col:
                    logger.info(f"  {col} -> {sanitized_col}: numeric (will aggregate with median)")
                else:
                    logger.info(f"  {col}: numeric (will aggregate with median)")
            else:
                # 分类列：使用清理后的列名作为 key
                categories = self.adata.obs[col].astype(str).unique()
                self.category_mapping[sanitized_col] = {cat: idx for idx, cat in enumerate(sorted(categories))}
                if sanitized_col != col:
                    logger.info(f"  {col} -> {sanitized_col}: {len(categories)} categories")
                else:
                    logger.info(f"  {col}: {len(categories)} categories")

    def _build_and_write_bins(self) -> None:
        """Build multi-resolution bins and write to Zarr."""
        logger.info(f"Building {self.config.zoom_levels} zoom levels")

        # Create output directory
        zarr_path = self.config.output_dir / "zarr" / "bins.zarr"
        zarr_path.parent.mkdir(parents=True, exist_ok=True)

        # Initialize Zarr writer
        writer = ZarrBinWriter(
            zarr_path,
            zoom_levels=self.config.zoom_levels,
            tile_size=self.config.tile_size,
            n_genes=len(self.selected_genes),
            n_categories={col: len(cats) for col, cats in self.category_mapping.items()},
            numeric_columns=self.numeric_columns,
            compressor=self.config.zarr_compressor,
            compression_level=self.config.zarr_compression_level,
        )

        # Initialize binner and aggregator
        binner = QuadtreeBinner(
            coordinate_range=self.config.coordinate_range,
            zoom_levels=self.config.zoom_levels,
        )

        # Get gene indices for expression extraction
        gene_indices = [
            self.adata.var_names.get_loc(gene)
            for gene in self.selected_genes
            if gene in self.adata.var_names
        ]

        aggregator = ExpressionAggregator(
            adata=self.adata,
            gene_indices=gene_indices,
            batch_size=self.config.batch_size,
        )

        # Process each zoom level
        for zoom in tqdm(range(self.config.zoom_levels), desc="Zoom levels"):
            logger.info(f"Processing zoom level {zoom}")

            # Compute bin assignments for all cells
            bin_assignments = binner.assign_bins(self.normalized_coords, zoom)

            # Get unique bins
            unique_bins = np.unique(bin_assignments, axis=0)
            logger.info(f"  Zoom {zoom}: {len(unique_bins):,} bins")

            # Aggregate data per bin
            bins_data = []
            for bin_x, bin_y in tqdm(unique_bins, desc=f"  Aggregating bins", leave=False):
                mask = (bin_assignments[:, 0] == bin_x) & (bin_assignments[:, 1] == bin_y)
                cell_indices = np.where(mask)[0]

                if len(cell_indices) == 0:
                    continue

                # Aggregate expression
                expr_mean, expr_max = aggregator.aggregate(cell_indices)

                # Aggregate categories
                category_counts = {}
                for sanitized_col, mapping in self.category_mapping.items():
                    orig_col = self.column_name_mapping[sanitized_col]
                    cats = self.adata.obs[orig_col].iloc[cell_indices].astype(str)
                    counts = np.zeros(len(mapping), dtype=np.uint32)
                    for cat, idx in mapping.items():
                        counts[idx] = (cats == cat).sum()
                    category_counts[sanitized_col] = counts

                # Aggregate numeric columns (median)
                numeric_medians = {}
                for sanitized_col in self.numeric_columns:
                    orig_col = self.column_name_mapping[sanitized_col]
                    values = self.adata.obs[orig_col].iloc[cell_indices].values
                    numeric_medians[sanitized_col] = float(np.median(values))

                bins_data.append({
                    "bin_x": int(bin_x),
                    "bin_y": int(bin_y),
                    "cell_count": len(cell_indices),
                    "cell_ids": cell_indices.astype(np.uint32),
                    "expression_mean": expr_mean,
                    "expression_max": expr_max,
                    "category_counts": category_counts,
                    "numeric_medians": numeric_medians,
                })

            # Write bins for this zoom level
            writer.write_zoom_level(zoom, bins_data)

        writer.finalize()
        logger.info(f"Zarr data written to: {zarr_path}")

    def _write_soma_store(self) -> None:
        """Write complete expression data to SOMA store."""
        logger.info("Writing TileDBSOMA store")

        soma_path = self.config.output_dir / "soma" / "experiment.soma"

        writer = SomaWriter(
            path=soma_path,
            zoom_levels=self.config.zoom_levels,
            coordinate_range=self.config.coordinate_range,
            compression=self.config.soma_compression,
            compression_level=self.config.soma_compression_level,
        )

        writer.write_from_adata(
            adata=self.adata,
            normalized_coords=self.normalized_coords,
            preaggregated_genes=self.selected_genes,
            category_columns=list(self.category_mapping.keys()),
            numeric_columns=self.numeric_columns,
            column_name_mapping=self.column_name_mapping,
        )

        writer.finalize()
        logger.info(f"SOMA store written to: {soma_path}")

    def _write_metadata(self) -> None:
        """Write metadata files."""
        logger.info("Writing metadata")

        metadata_dir = self.config.output_dir / "zarr"

        # Gene index mapping
        gene_index = {gene: idx for idx, gene in enumerate(self.selected_genes)}
        with open(metadata_dir / "gene_index.json", "w") as f:
            json.dump(gene_index, f, indent=2)

        # Main metadata
        metadata = {
            "dataset_name": self.config.dataset_name or self.config.input_path.stem,
            "dataset_description": self.config.dataset_description,
            "n_cells": int(self.adata.n_obs),
            "n_genes_total": int(self.adata.n_vars),
            "n_genes_preaggregated": len(self.selected_genes),
            "zoom_levels": self.config.zoom_levels,
            "tile_size": self.config.tile_size,
            "coordinate_range": self.config.coordinate_range,
            "preaggregated_genes": self.selected_genes,
            "categories": {
                col: {
                    "values": list(mapping.keys()),
                    "mapping": mapping,
                }
                for col, mapping in self.category_mapping.items()
            },
            "numeric_columns": self.numeric_columns,
            "umap_key": self.config.umap_key,
            "bounds": {
                "min_x": 0.0,
                "max_x": float(self.config.coordinate_range),
                "min_y": 0.0,
                "max_y": float(self.config.coordinate_range),
            },
            "soma_enabled": self.config.enable_soma,
            "soma_path": "soma/experiment.soma" if self.config.enable_soma else None,
            "all_genes_queryable": self.config.enable_soma,
        }

        with open(metadata_dir / "metadata.json", "w") as f:
            json.dump(metadata, f, indent=2)

        logger.info(f"Metadata written to: {metadata_dir}")
