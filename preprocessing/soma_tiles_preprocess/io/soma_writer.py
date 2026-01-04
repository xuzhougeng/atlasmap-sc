"""TileDBSOMA writer for complete single-cell data."""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp
import tiledbsoma as soma
from anndata import AnnData

logger = logging.getLogger(__name__)


class SomaWriter:
    """Writer for complete single-cell data in TileDBSOMA format."""

    def __init__(
        self,
        path: Path,
        zoom_levels: int,
        coordinate_range: float,
        compression: str = "zstd",
        compression_level: int = 9,
    ):
        """Initialize the SOMA writer.

        Args:
            path: Output path for SOMA Experiment
            zoom_levels: Number of zoom levels (for bin coordinate columns)
            coordinate_range: Coordinate normalization range
            compression: Compression algorithm
            compression_level: Compression level
        """
        self.path = Path(path)
        self.zoom_levels = zoom_levels
        self.coordinate_range = coordinate_range
        self.compression = compression
        self.compression_level = compression_level

        # Ensure parent directory exists
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def write_from_adata(
        self,
        adata: AnnData,
        normalized_coords: np.ndarray,
        preaggregated_genes: list[str],
        category_columns: list[str],
    ) -> None:
        """Write complete AnnData to SOMA Experiment.

        Args:
            adata: AnnData object with expression data
            normalized_coords: Normalized 2D coordinates [n_cells, 2]
            preaggregated_genes: List of genes in the Zarr pre-aggregation
            category_columns: Category columns to include in obs
        """
        logger.info(f"Writing SOMA Experiment to: {self.path}")

        # Compute bin coordinates for all zoom levels
        bin_coords = self._compute_bin_coordinates(normalized_coords)

        # Build DataFrames
        obs_df = self._build_obs_dataframe(adata, normalized_coords, bin_coords, category_columns)
        var_df = self._build_var_dataframe(adata, preaggregated_genes)

        # Get expression matrix
        X = adata.X
        if not sp.issparse(X):
            X = sp.csr_matrix(X)
        elif not sp.isspmatrix_csr(X):
            X = X.tocsr()

        n_cells = adata.n_obs
        n_genes = adata.n_vars

        # Create SOMA Experiment using context manager
        with soma.Experiment.create(str(self.path)) as exp:
            # Write obs DataFrame with proper domain
            logger.info(f"Writing obs DataFrame: {n_cells} cells")
            obs_schema = pa.Schema.from_pandas(obs_df)
            obs_arrow = pa.Table.from_pandas(obs_df, preserve_index=False)

            # Create obs with domain specification
            obs_uri = str(self.path / "obs")
            with soma.DataFrame.create(
                obs_uri,
                schema=obs_schema,
                index_column_names=["soma_joinid"],
                domain=[[0, n_cells - 1]],
            ) as obs_tbl:
                obs_tbl.write(obs_arrow)

            # Add obs to experiment
            exp.set("obs", soma.DataFrame.open(obs_uri))

            # Create measurement collection
            ms = exp.add_new_collection("ms")
            rna = ms.add_new_collection("RNA")

            # Write var DataFrame with proper domain
            logger.info(f"Writing var DataFrame: {n_genes} genes")
            var_schema = pa.Schema.from_pandas(var_df)
            var_arrow = pa.Table.from_pandas(var_df, preserve_index=False)

            var_uri = str(self.path / "ms" / "RNA" / "var")
            with soma.DataFrame.create(
                var_uri,
                schema=var_schema,
                index_column_names=["soma_joinid"],
                domain=[[0, n_genes - 1]],
            ) as var_tbl:
                var_tbl.write(var_arrow)

            # Add var to RNA measurement
            rna.set("var", soma.DataFrame.open(var_uri))

            # Create X collection and write expression matrix
            X_collection = rna.add_new_collection("X")

            logger.info(f"Writing expression matrix: {X.shape}")
            X_uri = str(self.path / "ms" / "RNA" / "X" / "data")
            self._write_sparse_matrix(X, X_uri)

            # Add X/data to collection
            X_collection.set("data", soma.SparseNDArray.open(X_uri))

        logger.info("SOMA Experiment written successfully")

    def _compute_bin_coordinates(
        self,
        normalized_coords: np.ndarray,
    ) -> dict[int, np.ndarray]:
        """Compute bin coordinates for all zoom levels.

        Args:
            normalized_coords: Normalized coordinates [n_cells, 2]

        Returns:
            Dict mapping zoom level to bin coordinates [n_cells, 2]
        """
        bin_coords = {}

        for zoom in range(self.zoom_levels):
            n_bins_per_axis = 2 ** zoom
            bin_size = self.coordinate_range / n_bins_per_axis

            # Compute bin indices
            bin_x = np.floor(normalized_coords[:, 0] / bin_size).astype(np.int32)
            bin_y = np.floor(normalized_coords[:, 1] / bin_size).astype(np.int32)

            # Clamp to valid range
            bin_x = np.clip(bin_x, 0, n_bins_per_axis - 1)
            bin_y = np.clip(bin_y, 0, n_bins_per_axis - 1)

            bin_coords[zoom] = np.column_stack([bin_x, bin_y])

        return bin_coords

    def _build_obs_dataframe(
        self,
        adata: AnnData,
        normalized_coords: np.ndarray,
        bin_coords: dict[int, np.ndarray],
        category_columns: list[str],
    ) -> pd.DataFrame:
        """Build obs DataFrame with spatial indices.

        Args:
            adata: AnnData object
            normalized_coords: Normalized coordinates
            bin_coords: Bin coordinates per zoom level
            category_columns: Category columns to include

        Returns:
            DataFrame with cell metadata and spatial indices
        """
        n_cells = adata.n_obs

        # Start with soma_joinid
        obs_data = {
            "soma_joinid": np.arange(n_cells, dtype=np.int64),
            "cell_id": adata.obs_names.astype(str).values,
            "x_normalized": normalized_coords[:, 0].astype(np.float32),
            "y_normalized": normalized_coords[:, 1].astype(np.float32),
        }

        # Add bin coordinates for each zoom level
        for zoom in range(self.zoom_levels):
            coords = bin_coords[zoom]
            obs_data[f"bin_z{zoom}_x"] = coords[:, 0]
            obs_data[f"bin_z{zoom}_y"] = coords[:, 1]

        # Add category columns from original obs
        for col in category_columns:
            if col in adata.obs.columns:
                obs_data[col] = adata.obs[col].astype(str).values

        return pd.DataFrame(obs_data)

    def _build_var_dataframe(
        self,
        adata: AnnData,
        preaggregated_genes: list[str],
    ) -> pd.DataFrame:
        """Build var DataFrame with gene metadata.

        Args:
            adata: AnnData object
            preaggregated_genes: List of pre-aggregated genes

        Returns:
            DataFrame with gene metadata
        """
        n_genes = adata.n_vars
        preagg_set = set(preaggregated_genes)

        # Build preagg_index: -1 if not in list, else index in list
        preagg_indices = []
        for gene in adata.var_names:
            if gene in preagg_set:
                preagg_indices.append(preaggregated_genes.index(gene))
            else:
                preagg_indices.append(-1)

        var_data = {
            "soma_joinid": np.arange(n_genes, dtype=np.int64),
            "gene_id": adata.var_names.astype(str).values,
            "gene_name": adata.var_names.astype(str).values,
            "is_preaggregated": np.array([g in preagg_set for g in adata.var_names]),
            "preagg_index": np.array(preagg_indices, dtype=np.int32),
        }

        return pd.DataFrame(var_data)

    def _write_sparse_matrix(
        self,
        X: sp.csr_matrix,
        path: str,
    ) -> None:
        """Write sparse matrix to SOMA SparseNDArray.

        Args:
            X: CSR sparse matrix
            path: Output path
        """
        # Convert to COO format for SOMA
        X_coo = X.tocoo()

        # Create SparseNDArray
        arr = soma.SparseNDArray.create(
            path,
            type=pa.float32(),
            shape=X.shape,
        )

        # Write data in chunks to manage memory
        chunk_size = 10_000_000  # 10M non-zero values per chunk
        n_nonzero = X_coo.nnz

        logger.info(f"Writing {n_nonzero:,} non-zero values in chunks")

        for start in range(0, n_nonzero, chunk_size):
            end = min(start + chunk_size, n_nonzero)

            # Create PyArrow table for this chunk
            chunk_table = pa.table({
                "soma_dim_0": pa.array(X_coo.row[start:end].astype(np.int64)),
                "soma_dim_1": pa.array(X_coo.col[start:end].astype(np.int64)),
                "soma_data": pa.array(X_coo.data[start:end].astype(np.float32)),
            })

            arr.write(chunk_table)

            if (end - start) == chunk_size:
                logger.info(f"  Written {end:,}/{n_nonzero:,} values")

    def finalize(self) -> None:
        """Finalize the SOMA store."""
        logger.info(f"Finalized SOMA store at: {self.path}")
