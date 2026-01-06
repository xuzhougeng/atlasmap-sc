#!/usr/bin/env python3
"""Print summary information about an H5AD file."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

try:
    import anndata
except ImportError:
    sys.exit(
        "Missing dependency: anndata.\n"
        "Install with: cd preprocessing && pip install -e .\n"
        "Or: make install-python"
    )


def truncate_list(items: list[str], max_show: int = 10) -> list[str]:
    """Return list with truncation notice if too long."""
    if len(items) <= max_show:
        return items
    return items[:max_show] + [f"... and {len(items) - max_show} more"]


def get_info(path: Path) -> dict[str, Any]:
    """Extract summary info from h5ad file."""
    adata = anndata.read_h5ad(path, backed="r")

    info: dict[str, Any] = {
        "file": str(path),
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
    }

    # X matrix info
    if adata.X is not None:
        x_type = type(adata.X).__name__
        info["X"] = {"shape": list(adata.shape), "type": x_type}
    else:
        info["X"] = None

    # layers
    info["layers"] = list(adata.layers.keys()) if adata.layers else []

    # obsm (embeddings/coordinates)
    obsm_info = {}
    for k in adata.obsm.keys():
        try:
            obsm_info[k] = list(adata.obsm[k].shape)
        except Exception:
            obsm_info[k] = "unknown"
    info["obsm"] = obsm_info
    info["has_X_umap"] = "X_umap" in obsm_info

    # obs columns
    obs_cols = list(adata.obs.columns)
    info["obs"] = {"n_columns": len(obs_cols), "columns": obs_cols}

    # var columns
    var_cols = list(adata.var.columns)
    info["var"] = {"n_columns": len(var_cols), "columns": var_cols}

    # raw
    info["has_raw"] = adata.raw is not None

    # uns keys
    info["uns_keys"] = list(adata.uns.keys()) if adata.uns else []

    # obsp / varp / varm
    info["obsp_keys"] = list(adata.obsp.keys()) if adata.obsp else []
    info["varp_keys"] = list(adata.varp.keys()) if adata.varp else []
    info["varm_keys"] = list(adata.varm.keys()) if adata.varm else []

    adata.file.close()
    return info


def print_human(info: dict[str, Any]) -> None:
    """Print info in human-readable format."""
    print(f"File: {info['file']}")
    print(f"Cells (n_obs): {info['n_obs']:,}")
    print(f"Genes (n_vars): {info['n_vars']:,}")
    print()

    # X
    if info["X"]:
        print(f"X matrix: {info['X']['shape']} ({info['X']['type']})")
    else:
        print("X matrix: None")

    # layers
    layers = info["layers"]
    if layers:
        print(f"Layers ({len(layers)}): {', '.join(truncate_list(layers))}")
    else:
        print("Layers: None")

    # obsm
    print()
    obsm = info["obsm"]
    if obsm:
        print(f"Embeddings/Coordinates (obsm) [{len(obsm)}]:")
        for k, shape in obsm.items():
            marker = " <-- X_umap" if k == "X_umap" else ""
            print(f"  {k}: {shape}{marker}")
    else:
        print("Embeddings/Coordinates (obsm): None")
    if not info["has_X_umap"]:
        print("  WARNING: X_umap not found (default for preprocessing; use --coord-key/--umap-key to choose another)")

    # obs
    print()
    obs = info["obs"]
    cols_display = truncate_list(obs["columns"], 15)
    print(f"obs columns ({obs['n_columns']}): {', '.join(cols_display)}")

    # var
    var = info["var"]
    cols_display = truncate_list(var["columns"], 15)
    print(f"var columns ({var['n_columns']}): {', '.join(cols_display)}")

    # raw
    print()
    print(f"Has raw: {info['has_raw']}")

    # uns
    uns = info["uns_keys"]
    if uns:
        print(f"uns keys ({len(uns)}): {', '.join(truncate_list(uns))}")

    # obsp/varp/varm
    for name in ("obsp_keys", "varp_keys", "varm_keys"):
        keys = info[name]
        if keys:
            print(f"{name.replace('_keys', '')} ({len(keys)}): {', '.join(truncate_list(keys))}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Print summary information about an H5AD file.")
    parser.add_argument("-i", "--input", required=True, type=Path, help="Path to H5AD file")
    parser.add_argument("--json", action="store_true", help="Output as JSON")
    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: file not found: {args.input}", file=sys.stderr)
        return 1

    info = get_info(args.input)

    if args.json:
        print(json.dumps(info, indent=2))
    else:
        print_human(info)

    return 0


if __name__ == "__main__":
    sys.exit(main())
