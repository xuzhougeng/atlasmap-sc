"""Command-line interface for SOMA-Tiles preprocessing."""

import click
from pathlib import Path
import logging
import sys

from .config import PreprocessConfig
from .pipeline import PreprocessingPipeline


def setup_logging(verbose: bool) -> None:
    """Configure logging based on verbosity level."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


@click.group()
@click.version_option(version="0.1.0")
def main():
    """SOMA-Tiles Preprocessing CLI.

    Convert H5AD single-cell data to multi-resolution Zarr bins
    for high-performance visualization.
    """
    pass


@main.command()
@click.option(
    "--input", "-i",
    "input_path",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to input H5AD file",
)
@click.option(
    "--output", "-o",
    "output_dir",
    required=True,
    type=click.Path(path_type=Path),
    help="Output directory for preprocessed data",
)
@click.option(
    "--umap-key",
    default="X_umap",
    help="Default coordinate key in adata.obsm (legacy name; default: X_umap)",
)
@click.option(
    "--coord-key",
    "coordinate_keys",
    multiple=True,
    help="Coordinate key(s) in adata.obsm (can specify multiple); overrides --umap-key if provided",
)
@click.option(
    "--zoom-levels", "-z",
    default=8,
    type=int,
    help="Number of zoom levels (1-12)",
)
@click.option(
    "--n-genes", "-g",
    default=500,
    type=int,
    help="Number of genes to pre-aggregate",
)
@click.option(
    "--all-expressed", "-a",
    "use_all_expressed",
    is_flag=True,
    default=False,
    help="Use all genes with non-zero expression instead of top N genes",
)
@click.option(
    "--min-cells",
    default=3,
    type=int,
    help="Minimum cells where gene must be expressed (only with --all-expressed)",
)
@click.option(
    "--category", "-c",
    "categories",
    multiple=True,
    help="Obs columns to include for aggregation (can specify multiple). If omitted, defaults to non-numeric obs columns.",
)
@click.option(
    "--numeric",
    "numeric_columns",
    multiple=True,
    help="Obs columns to treat as numeric (per-bin mean aggregation). Opt-in; ignored by default.",
)
@click.option(
    "--exclude-numeric",
    "exclude_numeric_columns",
    multiple=True,
    help="Obs columns to exclude from numeric aggregation (overrides --numeric).",
)
@click.option(
    "--chunk-size",
    default=256,
    type=int,
    help="Zarr chunk size",
)
@click.option(
    "--write-cell-ids",
    is_flag=True,
    default=False,
    help="Write per-bin cell id lists to Zarr (large; usually not needed).",
)
@click.option(
    "--batch-size",
    default=100000,
    type=int,
    help="Batch size for processing cells",
)
@click.option(
    "--no-soma",
    is_flag=True,
    default=False,
    help="Disable TileDBSOMA storage (skip full expression matrix)",
)
@click.option(
    "--name",
    default=None,
    help="Dataset name for metadata",
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    help="Enable verbose logging",
)
def run(
    input_path: Path,
    output_dir: Path,
    umap_key: str,
    coordinate_keys: tuple[str, ...],
    zoom_levels: int,
    n_genes: int,
    use_all_expressed: bool,
    min_cells: int,
    categories: tuple[str, ...],
    numeric_columns: tuple[str, ...],
    exclude_numeric_columns: tuple[str, ...],
    chunk_size: int,
    write_cell_ids: bool,
    batch_size: int,
    no_soma: bool,
    name: str | None,
    verbose: bool,
):
    """Run the preprocessing pipeline.

    Example:
        soma-preprocess run -i data.h5ad -o ./preprocessed -g 500 -z 8
    """
    setup_logging(verbose)
    logger = logging.getLogger(__name__)

    logger.info(f"Starting preprocessing pipeline")
    logger.info(f"Input: {input_path}")
    logger.info(f"Output: {output_dir}")

    # Build configuration
    config = PreprocessConfig(
        input_path=input_path,
        output_dir=output_dir,
        umap_key=umap_key,
        coordinate_keys=list(coordinate_keys),
        zoom_levels=zoom_levels,
        n_genes=n_genes,
        use_all_expressed=use_all_expressed,
        min_cells_expressed=min_cells,
        category_columns=list(categories),
        numeric_columns=list(numeric_columns),
        exclude_numeric_columns=list(exclude_numeric_columns),
        chunk_size=chunk_size,
        write_cell_ids=write_cell_ids,
        batch_size=batch_size,
        enable_soma=not no_soma,
        dataset_name=name,
    )

    # Validate configuration
    try:
        config.validate()
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Configuration error: {e}")
        raise click.ClickException(str(e))

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run pipeline
    pipeline = PreprocessingPipeline(config)
    try:
        pipeline.run()
        logger.info("Preprocessing completed successfully!")
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        raise click.ClickException(str(e))


@main.command()
@click.option(
    "--config", "-c",
    "config_path",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to YAML configuration file",
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    help="Enable verbose logging",
)
def from_config(config_path: Path, verbose: bool):
    """Run preprocessing from a YAML configuration file.

    Example:
        soma-preprocess from-config -c config.yaml
    """
    setup_logging(verbose)
    logger = logging.getLogger(__name__)

    logger.info(f"Loading configuration from: {config_path}")

    config = PreprocessConfig.from_yaml(config_path)
    config.validate()

    # Create output directory
    config.output_dir.mkdir(parents=True, exist_ok=True)

    # Run pipeline
    pipeline = PreprocessingPipeline(config)
    try:
        pipeline.run()
        logger.info("Preprocessing completed successfully!")
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        raise click.ClickException(str(e))


@main.command()
@click.option(
    "--output", "-o",
    "output_path",
    required=True,
    type=click.Path(path_type=Path),
    help="Output path for example configuration",
)
def init_config(output_path: Path):
    """Generate an example configuration file.

    Example:
        soma-preprocess init-config -o config.yaml
    """
    config = PreprocessConfig(
        input_path=Path("data.h5ad"),
        output_dir=Path("./preprocessed"),
        category_columns=["cell_type", "leiden"],
        marker_genes=["CD3D", "CD8A", "MS4A1", "CD14"],
        dataset_name="Example Dataset",
        dataset_description="Example single-cell dataset",
    )
    config.to_yaml(output_path)
    click.echo(f"Configuration saved to: {output_path}")


@main.command()
@click.option(
    "--input", "-i",
    "input_dir",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to Zarr output directory (contains bins.zarr, metadata.json)",
)
@click.option(
    "--coord",
    "coord_key",
    default=None,
    help="Coordinate key to visualize (e.g., X_umap, X_tsne); default uses bins.zarr",
)
@click.option(
    "--output", "-o",
    "output_dir",
    required=True,
    type=click.Path(path_type=Path),
    help="Output directory for figures",
)
@click.option(
    "--zoom-levels", "-z",
    "zoom_levels_str",
    default="3,5,7",
    help="Comma-separated zoom levels to visualize (default: 3,5,7)",
)
@click.option(
    "--category", "-c",
    "categories",
    multiple=True,
    help="Category columns to visualize (default: all)",
)
@click.option(
    "--gene", "-g",
    "genes",
    multiple=True,
    help="Genes to visualize (default: random selection)",
)
@click.option(
    "--n-genes",
    default=3,
    type=int,
    help="Number of random genes to visualize if --gene not specified",
)
@click.option(
    "--format", "-f",
    "fmt",
    default="png",
    type=click.Choice(["png", "svg", "pdf"]),
    help="Output format (default: png)",
)
@click.option(
    "--dpi",
    default=150,
    type=int,
    help="Figure resolution (default: 150)",
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    help="Enable verbose logging",
)
def visualize(
    input_dir: Path,
    coord_key: str | None,
    output_dir: Path,
    zoom_levels_str: str,
    categories: tuple[str, ...],
    genes: tuple[str, ...],
    n_genes: int,
    fmt: str,
    dpi: int,
    verbose: bool,
):
    """Visualize preprocessed Zarr output for validation.

    Generate static images showing category distributions and gene expression
    across multiple zoom levels.

    Example:
        soma-preprocess visualize -i ./output/zarr -o ./figures
        soma-preprocess visualize -i ./output/zarr -o ./figures -z 3,5,7 -g CD3D -g CD8A
    """
    setup_logging(verbose)
    logger = logging.getLogger(__name__)

    from .visualization import ZarrVisualizer

    # Parse zoom levels
    try:
        zoom_levels = [int(z.strip()) for z in zoom_levels_str.split(",")]
    except ValueError:
        raise click.ClickException(f"Invalid zoom levels: {zoom_levels_str}")

    logger.info(f"Visualizing Zarr data from: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Zoom levels: {zoom_levels}")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize visualizer
    try:
        visualizer = ZarrVisualizer(input_dir, coord_key=coord_key)
    except Exception as e:
        logger.error(f"Failed to load Zarr data: {e}")
        raise click.ClickException(str(e))

    # Generate report
    try:
        visualizer.generate_report(
            output_dir=output_dir,
            zoom_levels=zoom_levels,
            categories=list(categories) if categories else None,
            genes=list(genes) if genes else None,
            n_random_genes=n_genes,
            fmt=fmt,
            dpi=dpi,
        )
        logger.info("Visualization completed successfully!")
    except Exception as e:
        logger.error(f"Visualization failed: {e}")
        raise click.ClickException(str(e))


if __name__ == "__main__":
    main()
