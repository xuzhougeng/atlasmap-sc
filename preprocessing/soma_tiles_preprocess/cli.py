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
    help="Key for UMAP coordinates in adata.obsm",
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
    default=1,
    type=int,
    help="Minimum cells where gene must be expressed (only with --all-expressed)",
)
@click.option(
    "--category", "-c",
    "categories",
    multiple=True,
    help="Category columns to include (can specify multiple)",
)
@click.option(
    "--chunk-size",
    default=256,
    type=int,
    help="Zarr chunk size",
)
@click.option(
    "--batch-size",
    default=100000,
    type=int,
    help="Batch size for processing cells",
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
    zoom_levels: int,
    n_genes: int,
    use_all_expressed: bool,
    min_cells: int,
    categories: tuple[str, ...],
    chunk_size: int,
    batch_size: int,
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
        zoom_levels=zoom_levels,
        n_genes=n_genes,
        use_all_expressed=use_all_expressed,
        min_cells_expressed=min_cells,
        category_columns=list(categories),
        chunk_size=chunk_size,
        batch_size=batch_size,
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


if __name__ == "__main__":
    main()
