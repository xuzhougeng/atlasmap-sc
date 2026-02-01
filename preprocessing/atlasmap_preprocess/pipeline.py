"""Main preprocessing pipeline for AtlasMap."""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import json
import re

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


@dataclass(frozen=True)
class _CoarsenPlan:
    """Coarsen mapping from zoom z -> z-1."""

    order: np.ndarray  # indices to sort child bins by parent id
    starts: np.ndarray  # start indices (in sorted order) for each parent bin


@dataclass(frozen=True)
class _BinPyramid:
    """Bin pyramid derived from max zoom assignments."""

    max_zoom: int
    bin_index_per_cell_max: np.ndarray  # int32[n_cells] -> [0..n_bins_max)
    bin_coords_by_zoom: list[np.ndarray]  # per zoom: int32[n_bins,2]
    cell_count_by_zoom: list[np.ndarray]  # per zoom: uint32[n_bins]
    coarsen_plans: list[Optional[_CoarsenPlan]]  # index z: plan from z -> z-1 (z>0)


class PreprocessingPipeline:
    """Main preprocessing pipeline for converting H5AD to multi-resolution Zarr bins."""

    def __init__(self, config: PreprocessConfig):
        """Initialize the preprocessing pipeline.

        Args:
            config: Preprocessing configuration
        """
        self.config = config
        self.adata: Optional[sc.AnnData] = None
        self.normalized_coords: dict[str, np.ndarray] = {}
        self.coordinate_store_names: dict[str, str] = {}
        self.selected_genes: Optional[list[str]] = None
        self.category_mapping: dict[str, dict[str, int]] = {}
        self.numeric_columns: list[str] = []  # 存储检测到的数值列
        self.column_name_mapping: dict[str, str] = {}  # 原始列名 -> 清理后列名
        self.reusing_zarr: bool = False  # 标记是否在复用 zarr 输出

    def run(self) -> None:
        """Execute the full preprocessing pipeline."""
        logger.info("Starting preprocessing pipeline")

        # Step 1: Load H5AD
        self._load_data()

        # Check if Zarr output already exists
        if self._zarr_exists():
            logger.info("=" * 60)
            logger.info(f"Detected existing Zarr metadata at: {self._zarr_metadata_path()}")
            logger.info("Skipping Zarr generation and reusing existing metadata.")
            logger.info("To regenerate Zarr, delete the output_dir/zarr/ directory.")
            logger.info("=" * 60)
            
            self.reusing_zarr = True
            
            # Load existing metadata and align config
            self._load_existing_zarr_metadata()
            
            # Still need to compute normalized coordinates for SOMA
            self._process_coordinates()
            
            # Skip gene selection, category mapping, and bin building
            # (metadata provides gene list and category info)
        else:
            # Full pipeline: generate everything
            self.reusing_zarr = False
            
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

        # Step 7: Write/update metadata
        if self.reusing_zarr:
            self._update_existing_metadata_soma_flags()
        else:
            self._write_metadata()

        logger.info("Pipeline completed successfully")

    def _zarr_metadata_path(self) -> Path:
        """Return the path to zarr/metadata.json."""
        return self.config.output_dir / "zarr" / "metadata.json"

    def _zarr_exists(self) -> bool:
        """Check if Zarr output already exists."""
        return self._zarr_metadata_path().exists()

    def _load_existing_zarr_metadata(self) -> None:
        """Load existing Zarr metadata and configure pipeline state for reuse.
        
        This method:
        - Reads metadata.json
        - Sets selected_genes, category_mapping, numeric_columns from metadata
        - Aligns config settings (zoom_levels, coordinate_range, coord keys) with metadata
        - Builds column_name_mapping for consistency
        """
        metadata_path = self._zarr_metadata_path()
        logger.info(f"Loading existing Zarr metadata from: {metadata_path}")
        
        with open(metadata_path) as f:
            metadata = json.load(f)
        
        # Extract gene list
        self.selected_genes = metadata.get("preaggregated_genes", [])
        logger.info(f"Reusing {len(self.selected_genes)} pre-aggregated genes")
        
        # Extract numeric columns
        self.numeric_columns = metadata.get("numeric_columns", [])
        if self.numeric_columns:
            logger.info(f"Reusing {len(self.numeric_columns)} numeric columns")
        
        # Extract category columns (sanitized names)
        categories_meta = metadata.get("categories", {})
        category_columns = list(categories_meta.keys())
        
        # Rebuild category_mapping from metadata
        for col, cat_info in categories_meta.items():
            self.category_mapping[col] = cat_info.get("mapping", {})
        
        if category_columns:
            logger.info(f"Reusing {len(category_columns)} category columns: {category_columns}")
        
        # Override config settings to match existing metadata
        self.config.zoom_levels = metadata.get("zoom_levels", self.config.zoom_levels)
        self.config.coordinate_range = metadata.get("coordinate_range", self.config.coordinate_range)
        
        # Use the default coordinate system from metadata
        default_coord = metadata.get("default_coordinate_system") or metadata.get("umap_key")
        if default_coord:
            self.config.umap_key = default_coord
        
        # Only process the default coordinate system (needed for SOMA)
        self.config.coordinate_keys = [self.config.umap_key]
        
        logger.info(f"Config aligned: zoom_levels={self.config.zoom_levels}, "
                   f"coordinate_range={self.config.coordinate_range}, "
                   f"default_coord={self.config.umap_key}")
        
        # Build column_name_mapping using the same sanitize rule
        # Map sanitized -> original column name
        all_columns = list(category_columns) + self.numeric_columns
        for sanitized_col in all_columns:
            # Try to find the original column in adata.obs
            # Reverse the sanitize rule: sanitized has _ where original had .
            # Simple approach: iterate obs columns and check if sanitized matches
            for orig_col in self.adata.obs.columns:
                if self._sanitize_column_name(orig_col) == sanitized_col:
                    self.column_name_mapping[sanitized_col] = orig_col
                    break
            else:
                # Fallback: assume sanitized == original
                self.column_name_mapping[sanitized_col] = sanitized_col
        
        logger.info(f"Column name mapping built with {len(self.column_name_mapping)} entries")

    def _load_data(self) -> None:
        """Load H5AD file."""
        logger.info(f"Loading H5AD from: {self.config.input_path}")
        self.adata = sc.read_h5ad(self.config.input_path)
        logger.info(f"Loaded {self.adata.n_obs:,} cells, {self.adata.n_vars:,} genes")

    def _process_coordinates(self) -> None:
        """Extract and normalize coordinates (supports multiple coordinate systems)."""
        self.coordinate_store_names = self._build_coordinate_store_names()

        for key in self.config.coordinate_keys:
            logger.info(f"Extracting coordinates from obsm['{key}']")

            if key not in self.adata.obsm:
                available_keys = list(self.adata.obsm.keys())
                raise KeyError(
                    f"Coordinate key '{key}' not found. "
                    f"Available keys: {available_keys}"
                )

            coords = self.adata.obsm[key][:, :2].astype(np.float64)

            # Normalize to [0, coordinate_range)
            normalizer = CoordinateNormalizer(self.config.coordinate_range)
            self.normalized_coords[key] = normalizer.normalize(coords)

            logger.info(
                f"Normalized '{key}' to [{0}, {self.config.coordinate_range}), "
                f"shape: {self.normalized_coords[key].shape}"
            )

    def _sanitize_coordinate_id(self, key: str) -> str:
        """Convert an obsm key into a filesystem-safe identifier."""
        safe = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_")
        return safe or "coord"

    def _build_coordinate_store_names(self) -> dict[str, str]:
        """Build mapping: coordinate key -> Zarr store directory name."""
        out: dict[str, str] = {}
        used: set[str] = set()

        default_key = self.config.umap_key
        for key in self.config.coordinate_keys:
            if key == default_key:
                name = "bins.zarr"
            else:
                suffix = self._sanitize_coordinate_id(key)
                name = f"bins.{suffix}.zarr"
                i = 2
                while name in used or name == "bins.zarr":
                    name = f"bins.{suffix}_{i}.zarr"
                    i += 1

            out[key] = name
            used.add(name)

        return out

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
        # 默认不做自动推断：仅当用户显式指定时才作为数值列处理
        return col in self.config.numeric_columns

    def _build_category_mappings(self) -> None:
        """Build mappings from category names to integer indices."""
        logger.info("Building category mappings")

        # Build set of excluded columns (from config)
        exclude_set = set(self.config.exclude_category_columns)
        # Also sanitize exclusion names for matching
        exclude_sanitized = {self._sanitize_column_name(c) for c in exclude_set}

        # Determine which columns to include
        user_specified = bool(self.config.category_columns)
        columns = self.config.category_columns.copy()
        # If user didn't specify any category columns, default to categorical-only columns.
        # Numeric columns are ignored unless explicitly listed in config.numeric_columns.
        if not columns:
            columns = []
            for col in self.adata.obs.columns:
                series = self.adata.obs[col]
                if pd.api.types.is_bool_dtype(series) or not pd.api.types.is_numeric_dtype(series):
                    columns.append(col)

        # Opt-in numeric columns (median aggregation)
        for col in self.config.numeric_columns:
            if col in self.adata.obs.columns and col not in columns:
                columns.append(col)

        for col in columns:
            if col not in self.adata.obs.columns:
                logger.warning(f"Category column '{col}' not found in obs, skipping")
                continue

            # 清理列名（替换不适合 URL 的字符）
            sanitized_col = self._sanitize_column_name(col)

            # Check exclusion list (both original and sanitized names)
            if col in exclude_set or sanitized_col in exclude_sanitized:
                logger.info(f"  {col}: excluded by --exclude-category")
                continue

            self.column_name_mapping[sanitized_col] = col  # 清理后 -> 原始

            if self._is_numeric_column(col):
                # 数值列：跳过分类映射
                self.numeric_columns.append(sanitized_col)
                if sanitized_col != col:
                    logger.info(f"  {col} -> {sanitized_col}: numeric (will aggregate with mean)")
                else:
                    logger.info(f"  {col}: numeric (will aggregate with mean)")
            else:
                # 分类列：先检查基数（唯一值数量）
                n_unique = self.adata.obs[col].nunique()

                if n_unique > self.config.max_category_cardinality:
                    if user_specified:
                        # 用户显式指定了这个列，报错以提醒
                        raise ValueError(
                            f"Category column '{col}' has {n_unique:,} unique values, "
                            f"exceeding max_category_cardinality={self.config.max_category_cardinality:,}. "
                            f"High-cardinality columns (like cell IDs) are not suitable for per-bin category statistics. "
                            f"Either use --exclude-category {col} or increase --max-category-cardinality."
                        )
                    else:
                        # 自动选择的列，警告并跳过
                        logger.warning(
                            f"  {col}: skipped (cardinality {n_unique:,} > max {self.config.max_category_cardinality:,})"
                        )
                        continue

                # 分类列：使用清理后的列名作为 key
                categories = self.adata.obs[col].astype(str).unique()
                self.category_mapping[sanitized_col] = {cat: idx for idx, cat in enumerate(sorted(categories))}
                if sanitized_col != col:
                    logger.info(f"  {col} -> {sanitized_col}: {len(categories)} categories")
                else:
                    logger.info(f"  {col}: {len(categories)} categories")

    def _build_bin_pyramid(self, coords: np.ndarray, binner: QuadtreeBinner) -> _BinPyramid:
        """Compute bin assignments once (at max zoom) and coarsen to all zooms."""
        zoom_levels = self.config.zoom_levels
        max_zoom = zoom_levels - 1

        bin_assignments = binner.assign_bins(coords, max_zoom)
        n_bins_per_axis = binner.get_n_bins_per_axis(max_zoom)
        total_bins = n_bins_per_axis * n_bins_per_axis

        bin_id = (
            bin_assignments[:, 0].astype(np.int64) * n_bins_per_axis
            + bin_assignments[:, 1].astype(np.int64)
        )

        counts_full = np.bincount(bin_id, minlength=total_bins)
        nonzero_bin_ids = np.flatnonzero(counts_full)

        cell_count_max = counts_full[nonzero_bin_ids].astype(np.uint32, copy=False)
        bin_coords_max = np.column_stack(
            (nonzero_bin_ids // n_bins_per_axis, nonzero_bin_ids % n_bins_per_axis),
        ).astype(np.int32, copy=False)

        # Map each cell to its non-empty bin index at max zoom.
        bin_id_to_idx = np.full(total_bins, -1, dtype=np.int32)
        bin_id_to_idx[nonzero_bin_ids] = np.arange(nonzero_bin_ids.size, dtype=np.int32)
        bin_index_per_cell_max = bin_id_to_idx[bin_id].astype(np.int32, copy=False)

        bin_coords_by_zoom: list[np.ndarray] = [np.empty((0, 2), dtype=np.int32) for _ in range(zoom_levels)]
        cell_count_by_zoom: list[np.ndarray] = [np.empty((0,), dtype=np.uint32) for _ in range(zoom_levels)]
        coarsen_plans: list[Optional[_CoarsenPlan]] = [None for _ in range(zoom_levels)]

        bin_coords_by_zoom[max_zoom] = bin_coords_max
        cell_count_by_zoom[max_zoom] = cell_count_max

        cur_coords = bin_coords_max
        cur_counts = cell_count_max

        for zoom in range(max_zoom, 0, -1):
            parent_zoom = zoom - 1
            parent_bins_per_axis = binner.get_n_bins_per_axis(parent_zoom)

            parent_coords = (cur_coords // 2).astype(np.int32, copy=False)
            parent_bin_id = (
                parent_coords[:, 0].astype(np.int64) * parent_bins_per_axis
                + parent_coords[:, 1].astype(np.int64)
            )

            order = np.argsort(parent_bin_id, kind="stable")
            parent_bin_id_sorted = parent_bin_id[order]
            boundaries = np.flatnonzero(np.diff(parent_bin_id_sorted)) + 1
            starts = np.concatenate(([0], boundaries)).astype(np.int64, copy=False)
            unique_parent_ids = parent_bin_id_sorted[starts]

            parent_coords_unique = np.column_stack(
                (unique_parent_ids // parent_bins_per_axis, unique_parent_ids % parent_bins_per_axis),
            ).astype(np.int32, copy=False)

            counts_sorted = cur_counts[order].astype(np.uint64, copy=False)
            parent_counts = np.add.reduceat(counts_sorted, starts).astype(np.uint32, copy=False)

            coarsen_plans[zoom] = _CoarsenPlan(order=order, starts=starts)
            bin_coords_by_zoom[parent_zoom] = parent_coords_unique
            cell_count_by_zoom[parent_zoom] = parent_counts

            cur_coords = parent_coords_unique
            cur_counts = parent_counts

        return _BinPyramid(
            max_zoom=max_zoom,
            bin_index_per_cell_max=bin_index_per_cell_max,
            bin_coords_by_zoom=bin_coords_by_zoom,
            cell_count_by_zoom=cell_count_by_zoom,
            coarsen_plans=coarsen_plans,
        )

    def _build_and_write_bins(self) -> None:
        """Build multi-resolution bins and write to Zarr (per coordinate system)."""
        logger.info(f"Building {self.config.zoom_levels} zoom levels")

        zarr_dir = self.config.output_dir / "zarr"
        zarr_dir.mkdir(parents=True, exist_ok=True)

        # Initialize binner and aggregator (shared across coordinate systems)
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

        for coord_key in self.config.coordinate_keys:
            store_name = self.coordinate_store_names.get(coord_key, "bins.zarr")
            zarr_path = zarr_dir / store_name
            logger.info(f"Writing Zarr bins for '{coord_key}' -> {zarr_path}")

            # Initialize Zarr writer
            writer = ZarrBinWriter(
                zarr_path,
                zoom_levels=self.config.zoom_levels,
                tile_size=self.config.tile_size,
                n_genes=len(self.selected_genes),
                n_categories={col: len(cats) for col, cats in self.category_mapping.items()},
                numeric_columns=self.numeric_columns,
                write_cell_ids=getattr(self.config, "write_cell_ids", False),
                compressor=self.config.zarr_compressor,
                compression_level=self.config.zarr_compression_level,
            )

            coords = self.normalized_coords[coord_key]

            pyramid = self._build_bin_pyramid(coords, binner)
            max_zoom = pyramid.max_zoom
            n_bins_max = pyramid.bin_coords_by_zoom[max_zoom].shape[0]

            logger.info(
                f"Computed bin pyramid for '{coord_key}': "
                f"zoom {max_zoom} has {n_bins_max:,} non-empty bins"
            )

            # Expression aggregation: compute once at max zoom, then coarsen.
            expr_sum, expr_max = aggregator.aggregate_sum_max_by_bin(
                pyramid.bin_index_per_cell_max,
                n_bins=n_bins_max,
            )

            cur_sum = expr_sum
            cur_max = expr_max
            for zoom in tqdm(range(max_zoom, -1, -1), desc=f"Writing zoom levels ({coord_key})"):
                bin_coords_zoom = pyramid.bin_coords_by_zoom[zoom]
                cell_count_zoom = pyramid.cell_count_by_zoom[zoom]
                if bin_coords_zoom.size == 0:
                    continue

                expr_mean = (cur_sum / cell_count_zoom[:, None]).astype(np.float32, copy=False)
                writer.write_zoom_level_base(
                    zoom,
                    bin_coords=bin_coords_zoom,
                    cell_count=cell_count_zoom,
                    expression_mean=expr_mean,
                    expression_max=cur_max,
                )

                if zoom > 0:
                    plan = pyramid.coarsen_plans[zoom]
                    assert plan is not None
                    cur_sum = np.add.reduceat(cur_sum[plan.order], plan.starts, axis=0)
                    cur_max = np.maximum.reduceat(cur_max[plan.order], plan.starts, axis=0)

            # Category counts: encode once per column, then coarsen using the same plans.
            # Track columns that are skipped due to matrix size limits (for metadata consistency)
            skipped_categories: set[str] = set()

            for sanitized_col, mapping in tqdm(
                list(self.category_mapping.items()),
                desc=f"Writing categories ({coord_key})",
                leave=False,
            ):
                orig_col = self.column_name_mapping[sanitized_col]
                categories = list(mapping.keys())
                n_cats = len(categories)
                if n_cats == 0:
                    continue

                # Second-layer check: n_bins_max * n_cats must not exceed max_category_matrix_elems
                matrix_elems = n_bins_max * n_cats
                if matrix_elems > self.config.max_category_matrix_elems:
                    logger.warning(
                        f"Skipping category '{sanitized_col}': matrix size {n_bins_max:,} bins x {n_cats:,} cats = "
                        f"{matrix_elems:,} elements exceeds max_category_matrix_elems={self.config.max_category_matrix_elems:,}"
                    )
                    skipped_categories.add(sanitized_col)
                    continue

                codes = pd.Categorical(
                    self.adata.obs[orig_col].astype(str),
                    categories=categories,
                ).codes.astype(np.int32, copy=False)

                valid = codes >= 0
                if not np.any(valid):
                    continue

                bin_idx = pyramid.bin_index_per_cell_max[valid].astype(np.int64, copy=False)
                idx = bin_idx * n_cats + codes[valid].astype(np.int64, copy=False)
                counts_flat = np.bincount(idx, minlength=n_bins_max * n_cats)
                cur_counts = counts_flat.reshape((n_bins_max, n_cats)).astype(np.uint32, copy=False)

                for zoom in range(max_zoom, -1, -1):
                    writer.write_category_counts_column(zoom, column=sanitized_col, data=cur_counts)
                    if zoom > 0:
                        plan = pyramid.coarsen_plans[zoom]
                        assert plan is not None
                        cur_counts = np.add.reduceat(
                            cur_counts[plan.order].astype(np.uint64, copy=False),
                            plan.starts,
                            axis=0,
                        ).astype(np.uint32, copy=False)

            # Remove skipped categories from category_mapping to keep metadata consistent
            for col in skipped_categories:
                if col in self.category_mapping:
                    del self.category_mapping[col]

            # Numeric columns: aggregate by mean (faster than median) and coarsen.
            for sanitized_col in tqdm(
                list(self.numeric_columns),
                desc=f"Writing numeric ({coord_key})",
                leave=False,
            ):
                orig_col = self.column_name_mapping[sanitized_col]
                values = np.asarray(self.adata.obs[orig_col].values, dtype=np.float64)
                finite = np.isfinite(values)

                if not np.any(finite):
                    for zoom in range(max_zoom, -1, -1):
                        writer.write_numeric_medians_column(
                            zoom,
                            column=sanitized_col,
                            data=np.full(
                                (pyramid.bin_coords_by_zoom[zoom].shape[0],),
                                np.nan,
                                dtype=np.float32,
                            ),
                        )
                    continue

                bin_idx_finite = pyramid.bin_index_per_cell_max[finite]
                sum_per_bin = np.bincount(
                    bin_idx_finite,
                    weights=values[finite],
                    minlength=n_bins_max,
                ).astype(np.float64, copy=False)
                n_per_bin = np.bincount(bin_idx_finite, minlength=n_bins_max).astype(np.uint32, copy=False)

                cur_sum = sum_per_bin
                cur_n = n_per_bin
                for zoom in range(max_zoom, -1, -1):
                    with np.errstate(invalid="ignore", divide="ignore"):
                        mean = cur_sum / cur_n
                    mean = np.where(cur_n == 0, np.nan, mean).astype(np.float32, copy=False)
                    writer.write_numeric_medians_column(zoom, column=sanitized_col, data=mean)

                    if zoom > 0:
                        plan = pyramid.coarsen_plans[zoom]
                        assert plan is not None
                        cur_sum = np.add.reduceat(cur_sum[plan.order], plan.starts)
                        cur_n = np.add.reduceat(cur_n[plan.order].astype(np.uint64), plan.starts).astype(
                            np.uint32,
                            copy=False,
                        )

            # Optional per-bin cell ids (expensive; off by default).
            if getattr(self.config, "write_cell_ids", False):
                for zoom in tqdm(
                    range(max_zoom, -1, -1),
                    desc=f"Writing cell_ids ({coord_key})",
                    leave=False,
                ):
                    bin_assignments = binner.assign_bins(coords, zoom)
                    n_bins_axis = binner.get_n_bins_per_axis(zoom)
                    bin_id = (
                        bin_assignments[:, 0].astype(np.int64) * n_bins_axis
                        + bin_assignments[:, 1].astype(np.int64)
                    )
                    order = np.argsort(bin_id, kind="stable")
                    bin_id_sorted = bin_id[order]
                    boundaries = np.flatnonzero(np.diff(bin_id_sorted)) + 1
                    starts = np.concatenate(([0], boundaries))
                    ends = np.concatenate((boundaries, [order.size]))
                    cell_ids = [order[s:e].astype(np.uint32, copy=False) for s, e in zip(starts, ends)]
                    writer.write_cell_ids(zoom, cell_ids=cell_ids)

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
            normalized_coords=self.normalized_coords[self.config.umap_key],
            preaggregated_genes=self.selected_genes,
            category_columns=list(self.category_mapping.keys()),
            numeric_columns=self.numeric_columns,
            column_name_mapping=self.column_name_mapping,
        )

        writer.finalize()
        logger.info(f"SOMA store written to: {soma_path}")

    def _update_existing_metadata_soma_flags(self) -> None:
        """Update only SOMA-related fields in existing metadata.json.
        
        This is called when reusing Zarr output to update only:
        - soma_enabled
        - soma_path
        - all_genes_queryable
        """
        metadata_path = self._zarr_metadata_path()
        logger.info("Updating SOMA flags in existing metadata")
        
        # Load existing metadata
        with open(metadata_path) as f:
            metadata = json.load(f)
        
        # Update only SOMA-related fields
        metadata["soma_enabled"] = self.config.enable_soma
        metadata["soma_path"] = "soma/experiment.soma" if self.config.enable_soma else None
        metadata["all_genes_queryable"] = self.config.enable_soma
        
        # Write back
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=2)
        
        logger.info(f"Updated SOMA flags in: {metadata_path}")

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
            "coordinate_systems": [
                {"key": key, "zarr_path": self.coordinate_store_names.get(key, "bins.zarr")}
                for key in self.config.coordinate_keys
            ],
            "default_coordinate_system": self.config.umap_key,
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
