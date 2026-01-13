"""Expression aggregation for binned cells."""

import numpy as np
from numba import njit
from scipy import sparse
import scanpy as sc


@njit(cache=True)
def _aggregate_csr_by_bin(
    indptr: np.ndarray,
    indices: np.ndarray,
    data: np.ndarray,
    bin_index_per_cell: np.ndarray,
    n_bins: int,
    n_genes: int,
) -> tuple[np.ndarray, np.ndarray]:
    sums = np.zeros((n_bins, n_genes), dtype=np.float32)
    maxs = np.zeros((n_bins, n_genes), dtype=np.float32)

    n_cells = bin_index_per_cell.shape[0]
    for cell in range(n_cells):
        b = bin_index_per_cell[cell]
        start = indptr[cell]
        end = indptr[cell + 1]
        for k in range(start, end):
            g = indices[k]
            v = data[k]
            sums[b, g] += v
            if v > maxs[b, g]:
                maxs[b, g] = v

    return sums, maxs


@njit(cache=True)
def _aggregate_dense_by_bin(
    X: np.ndarray,
    bin_index_per_cell: np.ndarray,
    n_bins: int,
    n_genes: int,
) -> tuple[np.ndarray, np.ndarray]:
    sums = np.zeros((n_bins, n_genes), dtype=np.float32)
    maxs = np.zeros((n_bins, n_genes), dtype=np.float32)

    n_cells = bin_index_per_cell.shape[0]
    for cell in range(n_cells):
        b = bin_index_per_cell[cell]
        for g in range(n_genes):
            v = X[cell, g]
            sums[b, g] += v
            if v > maxs[b, g]:
                maxs[b, g] = v

    return sums, maxs


class ExpressionAggregator:
    """Aggregates gene expression values for groups of cells."""

    def __init__(
        self,
        adata: sc.AnnData,
        gene_indices: list[int],
        batch_size: int = 100000,
    ):
        """Initialize the expression aggregator.

        Args:
            adata: AnnData object with expression data
            gene_indices: Indices of genes to aggregate
            batch_size: Batch size for processing (for memory efficiency)
        """
        self.adata = adata
        self.gene_indices = np.array(gene_indices, dtype=np.int32)
        self.batch_size = batch_size
        self.n_genes = len(gene_indices)

        # Check if X is sparse
        self.is_sparse = sparse.issparse(adata.X)

    def aggregate(
        self,
        cell_indices: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Aggregate expression for a group of cells.

        Args:
            cell_indices: Indices of cells to aggregate

        Returns:
            Tuple of (mean_expression, max_expression), each shape (n_genes,)
        """
        if len(cell_indices) == 0:
            return (
                np.zeros(self.n_genes, dtype=np.float32),
                np.zeros(self.n_genes, dtype=np.float32),
            )

        # Extract expression matrix for selected cells and genes
        if self.is_sparse:
            expr = self.adata.X[cell_indices][:, self.gene_indices]
            if sparse.issparse(expr):
                expr = expr.toarray()
        else:
            expr = self.adata.X[np.ix_(cell_indices, self.gene_indices)]

        expr = np.asarray(expr, dtype=np.float32)

        # Compute aggregations
        mean_expr = np.mean(expr, axis=0).astype(np.float32)
        max_expr = np.max(expr, axis=0).astype(np.float32)

        return mean_expr, max_expr

    def aggregate_batched(
        self,
        cell_indices: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Aggregate expression using batched processing for memory efficiency.

        Args:
            cell_indices: Indices of cells to aggregate

        Returns:
            Tuple of (mean_expression, max_expression), each shape (n_genes,)
        """
        if len(cell_indices) == 0:
            return (
                np.zeros(self.n_genes, dtype=np.float32),
                np.zeros(self.n_genes, dtype=np.float32),
            )

        # For small groups, use direct aggregation
        if len(cell_indices) <= self.batch_size:
            return self.aggregate(cell_indices)

        # Batched aggregation
        sum_expr = np.zeros(self.n_genes, dtype=np.float64)
        max_expr = np.full(self.n_genes, -np.inf, dtype=np.float64)
        n_cells = len(cell_indices)

        for start in range(0, n_cells, self.batch_size):
            end = min(start + self.batch_size, n_cells)
            batch_indices = cell_indices[start:end]

            if self.is_sparse:
                batch_expr = self.adata.X[batch_indices][:, self.gene_indices]
                if sparse.issparse(batch_expr):
                    batch_expr = batch_expr.toarray()
            else:
                batch_expr = self.adata.X[np.ix_(batch_indices, self.gene_indices)]

            batch_expr = np.asarray(batch_expr, dtype=np.float64)

            sum_expr += np.sum(batch_expr, axis=0)
            max_expr = np.maximum(max_expr, np.max(batch_expr, axis=0))

        mean_expr = (sum_expr / n_cells).astype(np.float32)
        max_expr = max_expr.astype(np.float32)

        return mean_expr, max_expr

    def aggregate_sum_max_by_bin(
        self,
        bin_index_per_cell: np.ndarray,
        n_bins: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Aggregate expression for all bins in one pass.

        This is the fast path used by the preprocessing pipeline. It avoids per-bin
        slicing and sparse densification by aggregating over all cells at once.

        Args:
            bin_index_per_cell: Per-cell bin index (0..n_bins-1), shape (n_cells,)
            n_bins: Number of non-empty bins at the zoom level

        Returns:
            Tuple of (sum_expression, max_expression), each shape (n_bins, n_genes)
        """
        if n_bins <= 0 or self.n_genes == 0:
            return (
                np.zeros((max(n_bins, 0), self.n_genes), dtype=np.float32),
                np.zeros((max(n_bins, 0), self.n_genes), dtype=np.float32),
            )

        if bin_index_per_cell.dtype != np.int32:
            bin_index_per_cell = bin_index_per_cell.astype(np.int32, copy=False)

        if self.is_sparse:
            X = self.adata.X[:, self.gene_indices]
            if not sparse.isspmatrix_csr(X):
                X = X.tocsr()
            indptr = X.indptr
            indices = X.indices.astype(np.int32, copy=False)
            data = X.data.astype(np.float32, copy=False)
            sums, maxs = _aggregate_csr_by_bin(
                indptr=indptr,
                indices=indices,
                data=data,
                bin_index_per_cell=bin_index_per_cell,
                n_bins=n_bins,
                n_genes=self.n_genes,
            )
            return sums, maxs

        X = np.asarray(self.adata.X[:, self.gene_indices], dtype=np.float32)
        sums, maxs = _aggregate_dense_by_bin(
            X=X,
            bin_index_per_cell=bin_index_per_cell,
            n_bins=n_bins,
            n_genes=self.n_genes,
        )
        return sums, maxs
