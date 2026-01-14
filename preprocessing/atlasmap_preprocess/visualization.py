"""Visualization module for validating preprocessed Zarr output."""

import json
import logging
import random
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import zarr

logger = logging.getLogger(__name__)


class ZarrVisualizer:
    """Visualize preprocessed Zarr bin data."""

    def __init__(self, zarr_dir: Path, coord_key: str | None = None):
        """Initialize visualizer with Zarr output directory.

        Args:
            zarr_dir: Path to the Zarr output directory (contains bins.zarr, metadata.json)
            coord_key: Optional coordinate key to visualize (e.g., X_umap, X_tsne)
        """
        self.zarr_dir = Path(zarr_dir)
        self.metadata_path = self.zarr_dir / "metadata.json"
        self.gene_index_path = self.zarr_dir / "gene_index.json"

        # Load metadata
        with open(self.metadata_path) as f:
            self.metadata = json.load(f)

        # Resolve which Zarr store to open
        bins_rel = "bins.zarr"
        if coord_key:
            found = None
            for item in (self.metadata.get("coordinate_systems") or []):
                if item.get("key") == coord_key:
                    found = item.get("zarr_path")
                    break
            if found:
                bins_rel = str(found)
            else:
                default_key = self.metadata.get("default_coordinate_system") or self.metadata.get("umap_key")
                if coord_key != default_key:
                    available = [i.get("key") for i in (self.metadata.get("coordinate_systems") or [])]
                    raise ValueError(f"Coordinate '{coord_key}' not found. Available: {available}")

        self.bins_path = self.zarr_dir / bins_rel

        # Load gene index
        with open(self.gene_index_path) as f:
            self.gene_index = json.load(f)

        # Inverse gene index (index -> gene name)
        self.index_to_gene = {v: k for k, v in self.gene_index.items()}

        # Open Zarr store
        self.store = zarr.open(str(self.bins_path), mode="r")

        logger.info(f"Loaded Zarr data from {self.zarr_dir}")
        logger.info(f"  Store: {self.bins_path.name}")
        logger.info(f"  Cells: {self.metadata['n_cells']:,}")
        logger.info(f"  Genes: {self.metadata['n_genes_preaggregated']}")
        logger.info(f"  Zoom levels: {self.metadata['zoom_levels']}")
        logger.info(f"  Categories: {list(self.metadata['categories'].keys())}")

    @property
    def available_zoom_levels(self) -> list[int]:
        """Return list of available zoom levels."""
        return list(range(self.metadata["zoom_levels"]))

    @property
    def available_categories(self) -> list[str]:
        """Return list of available category columns."""
        return list(self.metadata["categories"].keys())

    @property
    def available_genes(self) -> list[str]:
        """Return list of pre-aggregated genes."""
        return self.metadata["preaggregated_genes"]

    def get_random_genes(self, n: int = 3) -> list[str]:
        """Get n random genes from pre-aggregated genes."""
        genes = self.available_genes
        return random.sample(genes, min(n, len(genes)))

    def _load_zoom_data(self, zoom: int) -> dict:
        """Load data for a specific zoom level."""
        zoom_group = self.store[f"zoom_{zoom}"]
        return {
            "bin_coords": zoom_group["bin_coords"][:],
            "cell_count": zoom_group["cell_count"][:],
            "expression_mean": zoom_group["expression_mean"][:],
            "expression_max": zoom_group["expression_max"][:],
        }

    def _load_category_counts(self, zoom: int, category: str) -> np.ndarray:
        """Load category counts for a specific zoom level and category."""
        zoom_group = self.store[f"zoom_{zoom}"]
        return zoom_group[f"category_counts/{category}"][:]

    def _coords_to_center(self, bin_coords: np.ndarray, zoom: int) -> tuple[np.ndarray, np.ndarray]:
        """Convert bin coordinates to center positions for plotting."""
        n_bins_per_axis = 2**zoom
        bin_size = self.metadata["coordinate_range"] / n_bins_per_axis
        x = bin_coords[:, 0] * bin_size + bin_size / 2
        y = bin_coords[:, 1] * bin_size + bin_size / 2
        return x, y

    def plot_category(
        self,
        zoom: int,
        category: str,
        output_path: Optional[Path] = None,
        figsize: tuple[float, float] = (10, 10),
        dpi: int = 150,
        point_size: Optional[float] = None,
    ) -> plt.Figure:
        """Plot bins colored by dominant category.

        Args:
            zoom: Zoom level to visualize
            category: Category column name
            output_path: Optional path to save figure
            figsize: Figure size in inches
            dpi: Figure resolution
            point_size: Point size (auto-calculated if None)

        Returns:
            Matplotlib figure
        """
        if category not in self.available_categories:
            raise ValueError(f"Category '{category}' not found. Available: {self.available_categories}")

        # Load data
        data = self._load_zoom_data(zoom)
        category_counts = self._load_category_counts(zoom, category)

        # Get dominant category for each bin
        dominant_idx = np.argmax(category_counts, axis=1)
        category_values = self.metadata["categories"][category]["values"]

        # Convert to center coordinates
        x, y = self._coords_to_center(data["bin_coords"], zoom)

        # Auto-calculate point size based on zoom level
        if point_size is None:
            point_size = max(1, 500 / (2**zoom))

        # Create figure
        fig, ax = plt.subplots(figsize=figsize)

        # Create color map
        n_categories = len(category_values)
        cmap = plt.cm.get_cmap("tab20" if n_categories <= 20 else "nipy_spectral")
        colors = [cmap(i / n_categories) for i in dominant_idx]

        # Plot
        scatter = ax.scatter(x, y, c=dominant_idx, cmap=cmap, s=point_size, alpha=0.8)

        # Add legend (limit to top 10 categories by count)
        category_totals = category_counts.sum(axis=0)
        top_indices = np.argsort(category_totals)[::-1][:10]
        handles = [
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=cmap(i / n_categories), markersize=8)
            for i in top_indices
        ]
        labels = [category_values[i] for i in top_indices]
        ax.legend(handles, labels, loc="upper right", fontsize=8)

        # Style
        ax.set_xlim(0, self.metadata["coordinate_range"])
        ax.set_ylim(0, self.metadata["coordinate_range"])
        ax.set_aspect("equal")
        ax.set_title(f"{category} - Zoom {zoom} ({len(x):,} bins)", fontsize=12)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        plt.tight_layout()

        if output_path:
            fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
            logger.info(f"Saved: {output_path}")

        return fig

    def plot_expression(
        self,
        zoom: int,
        gene: str,
        output_path: Optional[Path] = None,
        use_max: bool = False,
        figsize: tuple[float, float] = (10, 10),
        dpi: int = 150,
        point_size: Optional[float] = None,
        cmap: str = "viridis",
    ) -> plt.Figure:
        """Plot bins colored by gene expression.

        Args:
            zoom: Zoom level to visualize
            gene: Gene name
            output_path: Optional path to save figure
            use_max: If True, use max expression; otherwise use mean
            figsize: Figure size in inches
            dpi: Figure resolution
            point_size: Point size (auto-calculated if None)
            cmap: Colormap name

        Returns:
            Matplotlib figure
        """
        if gene not in self.gene_index:
            raise ValueError(f"Gene '{gene}' not found in pre-aggregated genes")

        gene_idx = self.gene_index[gene]

        # Load data
        data = self._load_zoom_data(zoom)
        expression = data["expression_max" if use_max else "expression_mean"][:, gene_idx]

        # Convert to center coordinates
        x, y = self._coords_to_center(data["bin_coords"], zoom)

        # Auto-calculate point size based on zoom level
        if point_size is None:
            point_size = max(1, 500 / (2**zoom))

        # Create figure
        fig, ax = plt.subplots(figsize=figsize)

        # Plot
        scatter = ax.scatter(x, y, c=expression, cmap=cmap, s=point_size, alpha=0.8)
        plt.colorbar(scatter, ax=ax, label="Expression" + (" (max)" if use_max else " (mean)"))

        # Style
        ax.set_xlim(0, self.metadata["coordinate_range"])
        ax.set_ylim(0, self.metadata["coordinate_range"])
        ax.set_aspect("equal")
        ax.set_title(f"{gene} - Zoom {zoom} ({len(x):,} bins)", fontsize=12)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        plt.tight_layout()

        if output_path:
            fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
            logger.info(f"Saved: {output_path}")

        return fig

    def plot_multi_zoom_category(
        self,
        category: str,
        zoom_levels: list[int],
        output_path: Optional[Path] = None,
        figsize_per_plot: tuple[float, float] = (5, 5),
        dpi: int = 150,
    ) -> plt.Figure:
        """Plot category across multiple zoom levels.

        Args:
            category: Category column name
            zoom_levels: List of zoom levels to show
            output_path: Optional path to save figure
            figsize_per_plot: Size per subplot
            dpi: Figure resolution

        Returns:
            Matplotlib figure
        """
        n_zooms = len(zoom_levels)
        fig, axes = plt.subplots(1, n_zooms, figsize=(figsize_per_plot[0] * n_zooms, figsize_per_plot[1]))

        if n_zooms == 1:
            axes = [axes]

        category_values = self.metadata["categories"][category]["values"]
        n_categories = len(category_values)
        cmap = plt.cm.get_cmap("tab20" if n_categories <= 20 else "nipy_spectral")

        for ax, zoom in zip(axes, zoom_levels):
            data = self._load_zoom_data(zoom)
            category_counts = self._load_category_counts(zoom, category)
            dominant_idx = np.argmax(category_counts, axis=1)

            x, y = self._coords_to_center(data["bin_coords"], zoom)
            point_size = max(1, 300 / (2**zoom))

            ax.scatter(x, y, c=dominant_idx, cmap=cmap, s=point_size, alpha=0.8)
            ax.set_xlim(0, self.metadata["coordinate_range"])
            ax.set_ylim(0, self.metadata["coordinate_range"])
            ax.set_aspect("equal")
            ax.set_title(f"Zoom {zoom}\n({len(x):,} bins)", fontsize=10)

        fig.suptitle(f"Category: {category}", fontsize=14, y=1.02)
        plt.tight_layout()

        if output_path:
            fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
            logger.info(f"Saved: {output_path}")

        return fig

    def plot_multi_zoom_expression(
        self,
        gene: str,
        zoom_levels: list[int],
        output_path: Optional[Path] = None,
        use_max: bool = False,
        figsize_per_plot: tuple[float, float] = (5, 5),
        dpi: int = 150,
        cmap: str = "viridis",
    ) -> plt.Figure:
        """Plot gene expression across multiple zoom levels.

        Args:
            gene: Gene name
            zoom_levels: List of zoom levels to show
            output_path: Optional path to save figure
            use_max: If True, use max expression; otherwise use mean
            figsize_per_plot: Size per subplot
            dpi: Figure resolution
            cmap: Colormap name

        Returns:
            Matplotlib figure
        """
        gene_idx = self.gene_index[gene]
        n_zooms = len(zoom_levels)
        fig, axes = plt.subplots(1, n_zooms, figsize=(figsize_per_plot[0] * n_zooms, figsize_per_plot[1]))

        if n_zooms == 1:
            axes = [axes]

        # Find global expression range for consistent colorbar
        all_expr = []
        for zoom in zoom_levels:
            data = self._load_zoom_data(zoom)
            expr = data["expression_max" if use_max else "expression_mean"][:, gene_idx]
            all_expr.extend(expr)
        vmin, vmax = np.percentile(all_expr, [0, 99])

        for ax, zoom in zip(axes, zoom_levels):
            data = self._load_zoom_data(zoom)
            expression = data["expression_max" if use_max else "expression_mean"][:, gene_idx]

            x, y = self._coords_to_center(data["bin_coords"], zoom)
            point_size = max(1, 300 / (2**zoom))

            scatter = ax.scatter(x, y, c=expression, cmap=cmap, s=point_size, alpha=0.8, vmin=vmin, vmax=vmax)
            ax.set_xlim(0, self.metadata["coordinate_range"])
            ax.set_ylim(0, self.metadata["coordinate_range"])
            ax.set_aspect("equal")
            ax.set_title(f"Zoom {zoom}\n({len(x):,} bins)", fontsize=10)

        # Add shared colorbar
        fig.colorbar(scatter, ax=axes, label="Expression" + (" (max)" if use_max else " (mean)"), shrink=0.8)

        fig.suptitle(f"Gene: {gene}", fontsize=14, y=1.02)
        plt.tight_layout()

        if output_path:
            fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
            logger.info(f"Saved: {output_path}")

        return fig

    def generate_report(
        self,
        output_dir: Path,
        zoom_levels: list[int] = [3, 5, 7],
        categories: Optional[list[str]] = None,
        genes: Optional[list[str]] = None,
        n_random_genes: int = 3,
        fmt: str = "png",
        dpi: int = 150,
    ) -> None:
        """Generate a complete visualization report.

        Args:
            output_dir: Directory to save figures
            zoom_levels: Zoom levels to visualize
            categories: Categories to visualize (None for all)
            genes: Genes to visualize (None for random selection)
            n_random_genes: Number of random genes if genes is None
            fmt: Output format (png, svg, pdf)
            dpi: Figure resolution
        """
        output_dir = Path(output_dir)

        # Filter zoom levels to available ones
        zoom_levels = [z for z in zoom_levels if z < self.metadata["zoom_levels"]]
        if not zoom_levels:
            zoom_levels = [self.metadata["zoom_levels"] - 1]

        # Categories
        categories = categories or self.available_categories
        category_dir = output_dir / "category"
        category_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Generating category visualizations for: {categories}")
        for cat in categories:
            # Individual zoom levels
            for zoom in zoom_levels:
                self.plot_category(
                    zoom=zoom,
                    category=cat,
                    output_path=category_dir / f"{cat}_zoom_{zoom}.{fmt}",
                    dpi=dpi,
                )
            # Multi-zoom comparison
            self.plot_multi_zoom_category(
                category=cat,
                zoom_levels=zoom_levels,
                output_path=category_dir / f"{cat}_multi_zoom.{fmt}",
                dpi=dpi,
            )

        # Genes
        genes = genes or self.get_random_genes(n_random_genes)
        expression_dir = output_dir / "expression"
        expression_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Generating expression visualizations for: {genes}")
        for gene in genes:
            # Individual zoom levels
            for zoom in zoom_levels:
                self.plot_expression(
                    zoom=zoom,
                    gene=gene,
                    output_path=expression_dir / f"{gene}_zoom_{zoom}.{fmt}",
                    dpi=dpi,
                )
            # Multi-zoom comparison
            self.plot_multi_zoom_expression(
                gene=gene,
                zoom_levels=zoom_levels,
                output_path=expression_dir / f"{gene}_multi_zoom.{fmt}",
                dpi=dpi,
            )

        logger.info(f"Report generated in: {output_dir}")
