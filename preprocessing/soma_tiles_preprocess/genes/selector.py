"""Gene selection for pre-aggregation."""

import logging
from typing import Optional

import numpy as np
import scanpy as sc

logger = logging.getLogger(__name__)


class GeneSelector:
    """Selects genes for pre-aggregation based on various criteria."""

    def __init__(self, adata: sc.AnnData):
        """Initialize the gene selector.

        Args:
            adata: AnnData object with expression data
        """
        self.adata = adata
        self.gene_names = list(adata.var_names)

    def select(
        self,
        n_genes: int = 500,
        hvg_n_top: int = 300,
        marker_genes: Optional[list[str]] = None,
    ) -> list[str]:
        """Select genes for pre-aggregation.

        Selection strategy:
        1. Highly variable genes (top hvg_n_top)
        2. User-specified marker genes
        3. Fill remaining with most expressed genes

        Args:
            n_genes: Total number of genes to select
            hvg_n_top: Number of highly variable genes to include
            marker_genes: User-specified marker genes to include

        Returns:
            List of selected gene names
        """
        selected: set[str] = set()

        # 1. Add highly variable genes
        hvg = self._get_highly_variable_genes(hvg_n_top)
        selected.update(hvg)
        logger.info(f"Added {len(hvg)} highly variable genes")

        # 2. Add user-specified marker genes
        if marker_genes:
            valid_markers = [g for g in marker_genes if g in self.gene_names]
            selected.update(valid_markers)
            logger.info(f"Added {len(valid_markers)} marker genes (of {len(marker_genes)} specified)")

        # 3. Fill remaining with most expressed genes
        if len(selected) < n_genes:
            n_remaining = n_genes - len(selected)
            most_expressed = self._get_most_expressed_genes(n_remaining * 2)  # Get extra for filtering
            for gene in most_expressed:
                if len(selected) >= n_genes:
                    break
                if gene not in selected:
                    selected.add(gene)
            logger.info(f"Added {len(selected) - len(hvg) - len(valid_markers if marker_genes else [])} additional genes")

        # Sort and limit
        final_genes = sorted(list(selected)[:n_genes])
        logger.info(f"Final selection: {len(final_genes)} genes")

        return final_genes

    def _get_highly_variable_genes(self, n_top: int) -> list[str]:
        """Get highly variable genes.

        Args:
            n_top: Number of top HVGs to return

        Returns:
            List of highly variable gene names
        """
        # Check if HVG already computed
        if "highly_variable" in self.adata.var.columns:
            hvg_mask = self.adata.var["highly_variable"]
            hvg_genes = self.adata.var_names[hvg_mask].tolist()

            # If there's a ranking, use it
            if "highly_variable_rank" in self.adata.var.columns:
                hvg_df = self.adata.var[hvg_mask].copy()
                hvg_df = hvg_df.sort_values("highly_variable_rank")
                return hvg_df.index.tolist()[:n_top]

            return hvg_genes[:n_top]

        # Compute HVG if not present
        logger.info("Computing highly variable genes...")
        adata_copy = self.adata.copy()

        try:
            # Normalize if needed
            if adata_copy.X.max() > 100:  # Likely raw counts
                sc.pp.normalize_total(adata_copy, target_sum=1e4)
                sc.pp.log1p(adata_copy)

            sc.pp.highly_variable_genes(adata_copy, n_top_genes=n_top, flavor="seurat_v3")
            hvg_mask = adata_copy.var["highly_variable"]
            return adata_copy.var_names[hvg_mask].tolist()
        except Exception as e:
            logger.warning(f"Failed to compute HVG: {e}")
            # Fallback to most variable by simple variance
            return self._get_most_variable_genes(n_top)

    def _get_most_variable_genes(self, n_top: int) -> list[str]:
        """Get genes with highest variance (simple fallback).

        Args:
            n_top: Number of genes to return

        Returns:
            List of gene names with highest variance
        """
        from scipy import sparse

        X = self.adata.X
        if sparse.issparse(X):
            # Compute variance for sparse matrix
            mean = np.array(X.mean(axis=0)).flatten()
            sq_mean = np.array(X.power(2).mean(axis=0)).flatten()
            var = sq_mean - mean ** 2
        else:
            var = np.var(X, axis=0)

        top_indices = np.argsort(var)[::-1][:n_top]
        return [self.gene_names[i] for i in top_indices]

    def _get_most_expressed_genes(self, n_top: int) -> list[str]:
        """Get genes with highest total expression.

        Args:
            n_top: Number of genes to return

        Returns:
            List of gene names with highest total expression
        """
        from scipy import sparse

        X = self.adata.X
        if sparse.issparse(X):
            total_expr = np.array(X.sum(axis=0)).flatten()
        else:
            total_expr = np.sum(X, axis=0)

        top_indices = np.argsort(total_expr)[::-1][:n_top]
        return [self.gene_names[i] for i in top_indices]
