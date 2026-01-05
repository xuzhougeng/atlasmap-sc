#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Seurat RDS to h5ad converter.

This script embeds the R extraction code and handles the full pipeline:
1. Write embedded R script to temp file
2. Run R script to export MTX/rows/cols/metadata from Seurat object
3. Read intermediate files and build AnnData
4. Write h5ad output
"""

import os
import sys
import argparse
import gzip
import tempfile
import subprocess
import shutil
import glob
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io as spio
from scipy import sparse

# Embedded R script for extracting data from Seurat RDS
R_SCRIPT = r'''
suppressPackageStartupMessages({
  library(Seurat)   # loads SeuratObject too
  library(Matrix)   # for writeMM/coercion
})

args <- commandArgs(trailingOnly = TRUE)
rds_file <- args[1]
prefix   <- args[2]

seu <- readRDS(rds_file)

# Always get the RNA assay
rna <- GetAssay(seu, "RNA")

# Fetch counts robustly across Seurat v4/v5:
# - v4 Assay: use slot = "counts"
# - v5 layered AssayX: use layer = "counts"
# - v5 with multiple layers: join layers first
counts <- tryCatch(
  GetAssayData(rna, slot = "counts"),
  error = function(e) {
    # Try with layer = "counts" first
    tryCatch(
      GetAssayData(rna, layer = "counts"),
      error = function(e2) {
        # If multiple layers exist, join them first
        message("Multiple layers detected, joining layers...")
        seu <- JoinLayers(seu)
        rna <- GetAssay(seu, "RNA")
        GetAssayData(rna, layer = "counts")
      }
    )
  }
)

# Ensure sparse dgCMatrix for writeMM
if (!inherits(counts, "dgCMatrix")) {
  counts <- Matrix(counts, sparse = TRUE)
}

mtx_file      <- paste(prefix, "mtx", sep = ".")
row_file      <- paste(prefix, "rows", sep = ".")
col_file      <- paste(prefix, "cols", sep = ".")
metadata_file <- paste(prefix, "metadata.csv", sep = ".")

# Write Matrix Market and dimnames
Matrix::writeMM(counts, mtx_file)
write.table(rownames(counts), row_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(colnames(counts), col_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
write.csv(seu[[]], metadata_file)

# Export dimensional reductions (PCA, UMAP, tSNE, etc.)
reductions <- Reductions(seu)
if (length(reductions) > 0) {
  message("[R] Found reductions: ", paste(reductions, collapse = ", "))
  for (red_name in reductions) {
    emb <- Embeddings(seu, reduction = red_name)
    out_file <- paste(prefix, paste0("reduction_", red_name, ".csv"), sep = ".")
    write.csv(emb, out_file, row.names = TRUE)
    message("[R] Exported reduction: ", red_name, " -> ", out_file)
  }
}

# Export variable features if available
var_features <- tryCatch(
  VariableFeatures(seu),
  error = function(e) NULL
)
if (!is.null(var_features) && length(var_features) > 0) {
  var_file <- paste(prefix, "var_features.txt", sep = ".")
  write.table(var_features, var_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
  message("[R] Exported ", length(var_features), " variable features")
}

message("[R] Export complete: ", prefix)
'''


def _open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def _smart_path(base, exts):
    """返回第一个存在的路径，如 base+'.mtx' 或 base+'.mtx.gz'。"""
    for ext in exts:
        p = f"{base}.{ext}"
        if os.path.exists(p):
            return p
    return None

def read_names(path):
    with _open_text(path) as fh:
        return [line.rstrip("\n\r") for line in fh]

def attach_metadata(adata, meta_path, specified_index=None):
    """Attach metadata to AnnData object.

    Returns:
        tuple: (adata, list of merged column names)
    """
    if not meta_path or not os.path.exists(meta_path):
        sys.stderr.write(f"[info] metadata 文件未提供或不存在：{meta_path}\n")
        return adata, []

    df = pd.read_csv(meta_path)
    # 优先使用用户指定索引列
    if specified_index and specified_index in df.columns:
        df = df.set_index(specified_index)
    else:
        # 常见列名自动识别
        candidates = ["barcode", "cell", "cell_id", "Cell", "cells", "obs", "obs_id"]
        ok = False
        for c in candidates:
            if c in df.columns:
                df = df.set_index(c)
                ok = True
                break
        if not ok:
            # 自动寻找与 obs_names 高重合度的列
            for c in df.columns:
                overlap = df[c].isin(adata.obs_names).mean()
                if overlap > 0.9:  # 90% 以上匹配就认定
                    df = df.set_index(c)
                    ok = True
                    break
        if not ok:
            sys.stderr.write("[warn] 未找到能与细胞条形码匹配的列，跳过合并元数据。\n")
            return adata, []

    # 对齐并合并
    df = df.loc[df.index.intersection(adata.obs_names)].copy()
    df = df.reindex(adata.obs_names)  # 与 adata.obs 顺序一致
    merged_columns = list(df.columns)
    adata.obs = pd.concat([adata.obs, df], axis=1)
    sys.stderr.write(f"[info] 已合并元数据列：{merged_columns}\n")
    return adata, merged_columns


def attach_reductions(adata, prefix):
    """Load dimensional reductions (PCA, UMAP, tSNE, etc.) from exported CSV files."""
    # Standard mapping: Seurat name -> AnnData obsm key
    name_mapping = {
        'pca': 'X_pca',
        'umap': 'X_umap',
        'tsne': 'X_tsne',
    }

    # Find all reduction files
    pattern = f"{prefix}.reduction_*.csv"
    reduction_files = glob.glob(pattern)

    for red_file in reduction_files:
        # Extract reduction name from filename: prefix.reduction_pca.csv -> pca
        basename = os.path.basename(red_file)
        # Handle prefix with dots
        red_name = basename.split('.reduction_')[-1].replace('.csv', '')

        df = pd.read_csv(red_file, index_col=0)

        # Align with adata.obs_names
        if not df.index.equals(adata.obs_names):
            # Try to reindex
            common = df.index.intersection(adata.obs_names)
            if len(common) < len(adata.obs_names) * 0.9:
                sys.stderr.write(f"[warn] reduction {red_name}: 仅 {len(common)}/{len(adata.obs_names)} 细胞匹配，跳过\n")
                continue
            df = df.reindex(adata.obs_names)

        # Convert to numpy array
        emb = df.values.astype(np.float32)

        # Determine obsm key
        obsm_key = name_mapping.get(red_name.lower(), f'X_{red_name}')
        adata.obsm[obsm_key] = emb
        sys.stderr.write(f"[info] 已添加降维: {obsm_key} (shape: {emb.shape})\n")

    return adata


def attach_var_features(adata, prefix):
    """Load variable features and mark them in adata.var."""
    var_file = f"{prefix}.var_features.txt"
    if not os.path.exists(var_file):
        return adata

    with open(var_file, 'r') as f:
        var_features = [line.strip() for line in f if line.strip()]

    if len(var_features) > 0:
        # Mark highly variable genes
        adata.var['highly_variable'] = adata.var_names.isin(var_features)
        n_hvg = adata.var['highly_variable'].sum()
        sys.stderr.write(f"[info] 已标记 {n_hvg} 个高变基因 (highly_variable)\n")

    return adata


def generate_soma_preprocess_hint(h5ad_path, columns):
    """Generate a suggested soma-preprocess command."""
    if not columns:
        return None

    # Build the command
    base_name = os.path.splitext(os.path.basename(h5ad_path))[0]
    output_dir = f"data/{base_name}"

    c_flags = ' '.join(f"-c '{col}'" for col in columns)
    cmd = f"soma-preprocess run --input {h5ad_path} --output {output_dir} --zoom-levels 11 -a --name \"{base_name}\" {c_flags}"

    return cmd

def run_r_extraction(rds_path, prefix, rscript_cmd="Rscript"):
    """Write embedded R script and execute it to extract data from RDS."""
    # Write R script to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
        f.write(R_SCRIPT)
        r_script_path = f.name

    try:
        sys.stderr.write(f"[info] Running R extraction: {rds_path} -> {prefix}.*\n")
        result = subprocess.run(
            [rscript_cmd, r_script_path, rds_path, prefix],
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            sys.stderr.write(f"[error] R script failed:\n{result.stderr}\n")
            sys.exit(1)
        if result.stderr:
            sys.stderr.write(result.stderr)
        if result.stdout:
            sys.stderr.write(result.stdout)
    finally:
        os.unlink(r_script_path)


def main():
    ap = argparse.ArgumentParser(
        description="Convert Seurat RDS to h5ad. Can also read pre-exported MTX files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Direct RDS conversion (recommended)
  %(prog)s input.rds -o output.h5ad

  # From pre-exported MTX files
  %(prog)s --prefix /path/to/prefix -o output.h5ad

  # Keep intermediate files
  %(prog)s input.rds -o output.h5ad --keep-intermediate
"""
    )
    ap.add_argument("input", nargs="?", help="Seurat RDS 文件路径")
    ap.add_argument("-o", "--output", help="输出 h5ad 文件名")
    ap.add_argument("--prefix", help="已导出的 MTX 文件前缀（跳过 R 提取步骤）")
    ap.add_argument("--metadata", help="可选：细胞元数据 CSV 路径", default=None)
    ap.add_argument("--metadata-index", help="可选：元数据中用于匹配细胞的列名（如 'barcode'）", default=None)
    ap.add_argument("--make-genes-unique", action="store_true", help="将重复基因名去重")
    ap.add_argument("--keep-intermediate", action="store_true", help="保留中间文件（mtx/rows/cols）")
    ap.add_argument("--rscript", default="Rscript", help="Rscript 可执行文件路径（默认：Rscript）")
    args = ap.parse_args()

    # Determine mode: RDS input or pre-exported prefix
    if args.prefix:
        # Use pre-exported files
        prefix = args.prefix
        out_h5ad = args.output or f"{prefix}.h5ad"
        cleanup_intermediate = False
    elif args.input:
        # Convert from RDS
        if not os.path.exists(args.input):
            sys.stderr.write(f"[error] RDS 文件不存在: {args.input}\n")
            sys.exit(1)

        # Determine output path and prefix
        base = os.path.splitext(args.input)[0]
        out_h5ad = args.output or f"{base}.h5ad"

        # Use temp directory for intermediate files
        if args.keep_intermediate:
            prefix = base
            cleanup_intermediate = False
        else:
            tmpdir = tempfile.mkdtemp(prefix="seurat2h5ad_")
            prefix = os.path.join(tmpdir, "data")
            cleanup_intermediate = True

        # Run R extraction
        run_r_extraction(args.input, prefix, args.rscript)
    else:
        ap.print_help()
        sys.stderr.write("\n[error] 请提供 RDS 文件路径或 --prefix 参数\n")
        sys.exit(1)

    mtx_path  = _smart_path(prefix, ["mtx", "mtx.gz"])
    rows_path = _smart_path(prefix, ["rows", "rows.gz"])
    cols_path = _smart_path(prefix, ["cols", "cols.gz"])

    if not (mtx_path and rows_path and cols_path):
        sys.stderr.write("找不到必要文件，请确认以下三者存在（允许 .gz）：\n"
                         f"  {prefix}.mtx   / {prefix}.mtx.gz\n"
                         f"  {prefix}.rows  / {prefix}.rows.gz\n"
                         f"  {prefix}.cols  / {prefix}.cols.gz\n")
        sys.exit(1)

    # 读取行列名
    genes = read_names(rows_path)
    cells = read_names(cols_path)

    # 读 mtx（Seurat 导出为 genes × cells），转置为 cells × genes
    # 这里用 scipy.io.mmread 更快更稳；之后仍然构建 AnnData（供 scanpy 使用）
    X = spio.mmread(mtx_path)  # COO
    if not sparse.issparse(X):
        X = sparse.coo_matrix(X)
    # 转置并转 CSR（scanpy/anndata 常用）
    X = X.transpose().tocsr()

    # 形状检查
    if X.shape[0] != len(cells) or X.shape[1] != len(genes):
        sys.stderr.write(f"[error] 维度不匹配：X(cells×genes)={X.shape}, cells={len(cells)}, genes={len(genes)}\n")
        sys.exit(2)

    # 构建 AnnData
    adata = sc.AnnData(X=X)
    adata.obs_names = pd.Index(cells, name="barcode")
    adata.var_names = pd.Index(genes, name="gene")

    if args.make_genes_unique:
        adata.var_names_make_unique()

    # 把原始计数也放到 layers['counts']
    adata.layers["counts"] = adata.X.copy()

    # 处理元数据
    meta_path = args.metadata if args.metadata else f"{prefix}.metadata.csv"
    merged_columns = []
    if os.path.exists(meta_path):
        adata, merged_columns = attach_metadata(adata, meta_path, args.metadata_index)

    # 加载降维数据 (PCA, UMAP, tSNE, etc.)
    adata = attach_reductions(adata, prefix)

    # 加载高变基因标记
    adata = attach_var_features(adata, prefix)

    # 写出 h5ad（lzf 压缩更快；也可用 'gzip' 更小）
    adata.write(out_h5ad, compression="lzf")
    print(f"[ok] 写出：{out_h5ad}")
    print(f"[ok] 形状：cells × genes = {adata.n_obs} × {adata.n_vars}")
    if len(adata.obsm) > 0:
        print(f"[ok] 降维：{list(adata.obsm.keys())}")

    # Cleanup intermediate files if needed
    if cleanup_intermediate:
        tmpdir = os.path.dirname(prefix)
        shutil.rmtree(tmpdir, ignore_errors=True)
        sys.stderr.write(f"[info] 已清理临时目录: {tmpdir}\n")

    # Print soma-preprocess hint
    if merged_columns:
        hint_cmd = generate_soma_preprocess_hint(out_h5ad, merged_columns)
        if hint_cmd:
            print(f"\n[hint] 下一步可运行 soma-preprocess：")
            print(f"  {hint_cmd}")


if __name__ == "__main__":
    main()
