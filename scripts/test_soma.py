#!/usr/bin/env python3
"""
Test script: validate a TileDB-SOMA Experiment produced by soma-tiles preprocessing.

What it does:
  - Opens an Experiment at --soma (experiment.soma directory)
  - Reads a small slice of obs/var (schema + a few rows)
  - Resolves a gene (specified or random)
  - Reads a sparse subset of X for (cells subset) x (that gene) and prints basic stats

Notes:
  - This script requires the preprocessing Python deps (tiledbsoma, pyarrow, numpy, pandas).
    Recommended install:
      cd preprocessing
      uv venv && source .venv/bin/activate
      uv pip install -e .
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
from typing import Any, NoReturn, Sequence, Tuple


def _die(msg: str, code: int = 2) -> NoReturn:
    print(msg, file=sys.stderr)
    raise SystemExit(code)


def _import_deps() -> Tuple[Any, Any, Any, Any]:
    try:
        import numpy as np  # type: ignore
        import pandas as pd  # type: ignore
        import pyarrow as pa  # type: ignore
        import tiledbsoma as soma  # type: ignore
    except Exception as e:
        _die(
            "Missing dependencies for SOMA test.\n"
            "Install/activate preprocessing env first, e.g.:\n"
            "  cd preprocessing\n"
            "  uv venv && source .venv/bin/activate\n"
            "  uv pip install -e .\n"
            f"Import error: {e}"
        )
    return np, pd, pa, soma


def _resolve_experiment_uri(path: str) -> str:
    p = os.path.abspath(os.path.expanduser(path))
    if os.path.isdir(p) and p.endswith(".soma"):
        return p
    # allow passing parent soma dir: .../soma -> .../soma/experiment.soma
    if os.path.isdir(p) and os.path.isdir(os.path.join(p, "experiment.soma")):
        return os.path.join(p, "experiment.soma")
    _die(f"--soma must point to an Experiment directory (*.soma). Got: {p}")
    raise AssertionError("unreachable")


def _read_df_head(df: Any, limit: int) -> "pa.Table":  # type: ignore[name-defined]
    """
    Read up to `limit` rows from a SOMA DataFrame without accidentally materializing the full table.

    Strategy:
      - try value_filter (fast)
      - otherwise iterate read() batches until we have enough rows
    """
    import pyarrow as pa  # type: ignore

    limit = int(limit)
    if limit <= 0:
        return pa.table({})

    # Prefer server-side filter if supported.
    try:
        reader = df.read(value_filter=f"soma_joinid < {limit}")
        if hasattr(reader, "concat"):
            return reader.concat()
    except Exception:
        pass

    # Batch-limited fallback.
    try:
        reader = df.read()
    except Exception as e:
        _die(f"Failed to read DataFrame: {e}")

    tables = []
    n = 0
    for chunk in reader:
        tbl = chunk if isinstance(chunk, pa.Table) else pa.Table.from_batches([chunk])
        if tbl.num_rows <= 0:
            continue
        tables.append(tbl)
        n += int(tbl.num_rows)
        if n >= limit:
            break
    if not tables:
        return pa.table({})
    out = pa.concat_tables(tables, promote=True)
    return out.slice(0, limit)


def _get_ms_rna(exp: Any) -> Any:
    # Common access patterns for tiledbsoma Experiment.
    # exp.ms is a collection with measurement(s)
    ms = getattr(exp, "ms", None)
    if ms is None:
        _die("Experiment has no .ms; unexpected SOMA layout.")
    try:
        return ms["RNA"]
    except Exception:
        # Try first measurement if name differs
        try:
            keys = list(ms.keys())  # type: ignore[attr-defined]
            if keys:
                return ms[keys[0]]
        except Exception:
            pass
    _die("Failed to access measurement: exp.ms['RNA'] (and no fallback found).")
    raise AssertionError("unreachable")


def _get_var_df(rna: Any) -> Any:
    for attr in ("var",):
        if hasattr(rna, attr):
            return getattr(rna, attr)
    # some versions expose mapping access
    try:
        return rna["var"]
    except Exception:
        pass
    try:
        return rna.get("var")  # type: ignore[attr-defined]
    except Exception:
        pass
    _die("Failed to access var DataFrame under measurement (expected ms/RNA/var).")
    raise AssertionError("unreachable")


def _get_obs_df(exp: Any) -> Any:
    if hasattr(exp, "obs"):
        return exp.obs
    try:
        return exp["obs"]
    except Exception:
        pass
    try:
        return exp.get("obs")  # type: ignore[attr-defined]
    except Exception:
        pass
    _die("Failed to access obs DataFrame (expected /obs).")
    raise AssertionError("unreachable")


def _get_x_sparse(rna: Any) -> Any:
    # Expect ms/RNA/X/data
    X = getattr(rna, "X", None)
    if X is not None:
        try:
            return X["data"]
        except Exception:
            pass
    # Fallback mapping access
    try:
        X = rna["X"]
        return X["data"]
    except Exception:
        pass
    _die("Failed to access sparse matrix (expected ms/RNA/X/data).")
    raise AssertionError("unreachable")


def _resolve_gene(var_tbl: "pa.Table", gene: str) -> Tuple[int, str]:  # type: ignore[name-defined]
    # Returns: (gene_joinid, gene_id)
    import pyarrow.compute as pc  # type: ignore

    if gene and gene != "random":
        mask = pc.equal(var_tbl["gene_id"], gene)
        idx = pc.indices_nonzero(mask)
        if len(idx) > 0:
            i = int(idx[0].as_py())
            return int(var_tbl["soma_joinid"][i].as_py()), str(var_tbl["gene_id"][i].as_py())

        # try gene_name
        if "gene_name" in var_tbl.column_names:
            mask = pc.equal(var_tbl["gene_name"], gene)
            idx = pc.indices_nonzero(mask)
            if len(idx) > 0:
                i = int(idx[0].as_py())
                return int(var_tbl["soma_joinid"][i].as_py()), str(var_tbl["gene_id"][i].as_py())

        raise ValueError(f"Gene not found in var head: {gene}")

    # random: prefer preaggregated genes if available
    if "is_preaggregated" in var_tbl.column_names:
        pre = pc.equal(var_tbl["is_preaggregated"], True)  # noqa: E712
        idx = pc.indices_nonzero(pre)
        if len(idx) > 0:
            i = int(idx[random.randrange(len(idx))].as_py())
            return int(var_tbl["soma_joinid"][i].as_py()), str(var_tbl["gene_id"][i].as_py())

    i = random.randrange(var_tbl.num_rows)
    return int(var_tbl["soma_joinid"][i].as_py()), str(var_tbl["gene_id"][i].as_py())


def _first_table(reader: Any) -> "pa.Table | None":  # type: ignore[name-defined]
    """Extract first table from a SOMA reader result."""
    import pyarrow as pa  # type: ignore

    if reader is None:
        return None
    # If it's already a table, return it
    if isinstance(reader, pa.Table):
        return reader
    # SparseNDArrayRead has .tables() iterator
    if hasattr(reader, "tables"):
        for tbl in reader.tables():
            if tbl is not None and tbl.num_rows > 0:
                return tbl
        return None
    # Try concat (some readers support this)
    if hasattr(reader, "concat"):
        try:
            return reader.concat()
        except Exception:
            pass
    # Try iterating directly
    try:
        for chunk in reader:
            if isinstance(chunk, pa.Table):
                return chunk
            if hasattr(chunk, "to_table"):
                return chunk.to_table()
    except TypeError:
        pass
    return None


def _pick_cells(obs_limit_tbl: "pa.Table", n_cells: int) -> Sequence[int]:  # type: ignore[name-defined]
    # Prefer using soma_joinid values from the table slice.
    joinids = [int(x.as_py()) for x in obs_limit_tbl["soma_joinid"]]
    if not joinids:
        _die("No cells in obs slice; cannot pick cells.")
    if n_cells >= len(joinids):
        return joinids
    return random.sample(joinids, k=int(n_cells))


def _read_sparse_gene_subset(
    np: Any,
    arr: Any,
    cell_ids: Sequence[int],
    gene_joinid: int,
    max_rows: int,
) -> "pa.Table":  # type: ignore[name-defined]
    # Best-effort: try coordinate-restricted reads first (fast), fallback to scanning.
    coords = (np.asarray(cell_ids, dtype=np.int64), np.asarray([gene_joinid], dtype=np.int64))

    # Attempt common signatures
    for call in (
        lambda: arr.read(coords=coords),
        lambda: arr.read(coords=list(coords)),
        lambda: arr.read(coords={"soma_dim_0": coords[0], "soma_dim_1": coords[1]}),
        lambda: arr.read(coords),
    ):
        try:
            t = _first_table(call())
            if t is None:
                continue
            return t.slice(0, int(max_rows))
        except Exception:
            continue

    # Fallback: scan chunks until enough hits (may be slow on large datasets).
    import pyarrow as pa  # type: ignore
    import pyarrow.compute as pc  # type: ignore

    want_cells = set(int(x) for x in cell_ids)
    batches = []
    seen = 0
    try:
        reader = arr.read()
    except Exception as e:
        _die(f"Failed to read sparse matrix: {e}")

    # SparseNDArrayRead returns tables via .tables() iterator
    table_iter = reader.tables() if hasattr(reader, "tables") else reader
    for chunk in table_iter:
        tbl = chunk if isinstance(chunk, pa.Table) else pa.Table.from_batches([chunk])
        if not tbl.num_rows:
            continue
        # Filter by gene + cell subset
        m1 = pc.equal(tbl["soma_dim_1"], gene_joinid)
        tbl = tbl.filter(m1)
        if not tbl.num_rows:
            continue
        m0 = pc.is_in(tbl["soma_dim_0"], value_set=pa.array(sorted(want_cells), type=pa.int64()))
        tbl = tbl.filter(m0)
        if not tbl.num_rows:
            continue
        batches.append(tbl)
        seen += tbl.num_rows
        if seen >= int(max_rows):
            break
    if not batches:
        return pa.table({"soma_dim_0": [], "soma_dim_1": [], "soma_data": []})
    out = pa.concat_tables(batches, promote=True)
    return out.slice(0, int(max_rows))


def truncate_list(items: list[str], max_show: int = 10) -> list[str]:
    """Return list with truncation notice if too long."""
    if len(items) <= max_show:
        return items
    return items[:max_show] + [f"... and {len(items) - max_show} more"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--soma",
        required=True,
        help="Path to TileDB-SOMA Experiment directory (e.g. /.../soma/experiment.soma or /.../soma).",
    )
    p.add_argument("--gene", default="random", help="Gene ID/name or 'random'")
    p.add_argument("--n-cells", type=int, default=256, help="Number of cells to sample (from obs head)")
    p.add_argument("--head", type=int, default=5000, help="Read this many obs/var rows (max) for sampling")
    p.add_argument("--max-nnz", type=int, default=200000, help="Max nnz rows to print/analyze from X subset")
    p.add_argument("--seed", type=int, default=0, help="Random seed (0 means deterministic default)")
    p.add_argument("--json", action="store_true", help="Output as JSON")
    return p.parse_args()


def _get_cell_type_stats(obs_tbl: "pa.Table") -> dict[str, Any] | None:  # type: ignore[name-defined]
    """Get cell_type statistics from obs table. Returns None if no cell_type column found."""
    import pyarrow as pa  # type: ignore
    import pyarrow.compute as pc  # type: ignore

    col_name = None
    if "cell_type" in obs_tbl.column_names:
        col_name = "cell_type"
    elif "celltype" in obs_tbl.column_names:
        col_name = "celltype"

    if col_name is None:
        return None

    col = obs_tbl[col_name]
    # Count occurrences of each value using value_counts
    try:
        value_counts = pc.value_counts(col)
        counts = {}
        for i in range(value_counts.num_rows):
            val = value_counts["values"][i]
            count = value_counts["counts"][i]
            if val.is_valid:
                counts[str(val.as_py())] = int(count.as_py())
        return {"column": col_name, "counts": counts}
    except Exception:
        # Fallback: manual counting
        unique_values = pc.unique(col)
        counts = {}
        for val in unique_values:
            if val.is_valid:
                val_str = val.as_py()
                mask = pc.equal(col, val)
                count = int(pc.sum(pc.cast(mask, pa.int64())).as_py())
                counts[str(val_str)] = count
        return {"column": col_name, "counts": counts}


def get_info(
    exp: Any,
    obs_df: Any,
    var_df: Any,
    X: Any,
    obs_tbl: "pa.Table",  # type: ignore[name-defined]
    var_tbl: "pa.Table",  # type: ignore[name-defined]
    gene_joinid: int,
    gene_id: str,
    cell_ids: Sequence[int],
    x_tbl: "pa.Table | None",  # type: ignore[name-defined]
    uri: str,
    args: argparse.Namespace,
) -> dict[str, Any]:
    """Collect all information into a dictionary."""
    import pyarrow as pa  # type: ignore

    info: dict[str, Any] = {
        "experiment": uri,
        "obs": {
            "n_rows": int(obs_tbl.num_rows),
            "n_columns": len(obs_tbl.column_names),
            "columns": list(obs_tbl.column_names),
        },
        "var": {
            "n_rows": int(var_tbl.num_rows),
            "n_columns": len(var_tbl.column_names),
            "columns": list(var_tbl.column_names),
        },
        "gene": {
            "gene_id": gene_id,
            "soma_joinid": gene_joinid,
            "requested": args.gene,
        },
        "cells": {
            "sampled": len(cell_ids),
            "from_obs_head": int(obs_tbl.num_rows),
            "sample_size": int(args.n_cells),
        },
    }

    # Cell type statistics
    cell_type_stats = _get_cell_type_stats(obs_tbl)
    if cell_type_stats:
        info["cell_type"] = cell_type_stats

    if x_tbl is not None and x_tbl.num_rows > 0:
        x = pa.compute.cast(x_tbl["soma_data"], pa.float64())
        mean_v = float(pa.compute.mean(x).as_py())
        max_v = float(pa.compute.max(x).as_py())
        min_v = float(pa.compute.min(x).as_py())
        nnz = int(x_tbl.num_rows)

        show = min(10, nnz)
        dim0 = x_tbl["soma_dim_0"].to_pylist()[:show]
        data = x_tbl["soma_data"].to_pylist()[:show]
        entries = [{"cell_joinid": int(dim0[i]), "value": float(data[i])} for i in range(show)]

        info["X_subset"] = {
            "nnz": nnz,
            "min": min_v,
            "mean": mean_v,
            "max": max_v,
            "first_entries": entries,
        }
    else:
        info["X_subset"] = {"nnz": 0, "warning": "Gene may be unexpressed in sampled cells"}

    return info


def print_human(info: dict[str, Any]) -> None:
    """Print info in human-readable format."""
    print(f"Experiment: {info['experiment']}")
    print()

    obs = info["obs"]
    print(f"obs: rows≈{obs['n_rows']:,} cols={obs['n_columns']}")
    cols_display = truncate_list(obs["columns"], 15)
    print(f"  columns: {', '.join(cols_display)}")
    print()

    var = info["var"]
    print(f"var: rows≈{var['n_rows']:,} cols={var['n_columns']}")
    cols_display = truncate_list(var["columns"], 15)
    print(f"  columns: {', '.join(cols_display)}")
    print()

    gene = info["gene"]
    print(f"Gene: {gene['gene_id']} (soma_joinid={gene['soma_joinid']}, requested={gene['requested']})")
    print()

    cells = info["cells"]
    print(f"Cells: sampled={cells['sampled']} (from obs head={cells['from_obs_head']}, sample_size={cells['sample_size']})")
    print()

    # Cell type statistics
    if "cell_type" in info:
        ct_stats = info["cell_type"]
        print(f"Cell type statistics ({ct_stats['column']}):")
        sorted_counts = sorted(ct_stats["counts"].items(), key=lambda x: x[1], reverse=True)
        for ct_name, count in sorted_counts:
            print(f"  {ct_name}: {count:,}")
        print()

    x_subset = info["X_subset"]
    if x_subset["nnz"] == 0:
        print(f"X subset: nnz=0")
        if "warning" in x_subset:
            print(f"  WARNING: {x_subset['warning']}")
    else:
        print(f"X subset: nnz={x_subset['nnz']:,} min={x_subset['min']:.6g} mean={x_subset['mean']:.6g} max={x_subset['max']:.6g}")
        print("  First entries: (cell_joinid, value)")
        for entry in x_subset["first_entries"]:
            print(f"    {entry['cell_joinid']}\t{entry['value']}")


def main() -> int:
    args = parse_args()
    random.seed(int(args.seed))

    np, pd, pa, soma = _import_deps()

    uri = _resolve_experiment_uri(args.soma)

    try:
        exp = soma.Experiment.open(uri)
    except Exception as e:
        _die(f"Failed to open experiment: {e}")

    with exp:
        obs_df = _get_obs_df(exp)
        rna = _get_ms_rna(exp)
        var_df = _get_var_df(rna)
        X = _get_x_sparse(rna)

        obs_tbl = _read_df_head(obs_df, limit=int(args.head))
        var_tbl = _read_df_head(var_df, limit=int(args.head))

        try:
            gene_joinid, gene_id = _resolve_gene(var_tbl, args.gene)
        except ValueError:
            if args.gene and args.gene != "random":
                # For explicit gene, fall back to reading full var (usually small, e.g. ~20k genes).
                try:
                    reader = var_df.read()
                    var_full = reader.concat() if hasattr(reader, "concat") else _read_df_head(var_df, limit=10**9)
                except Exception as e:
                    _die(f"Failed to read full var for gene lookup: {e}")
                gene_joinid, gene_id = _resolve_gene(var_full, args.gene)
            else:
                raise

        cell_ids = _pick_cells(obs_tbl, n_cells=int(args.n_cells))

        x_tbl = _read_sparse_gene_subset(
            np=np,
            arr=X,
            cell_ids=cell_ids,
            gene_joinid=gene_joinid,
            max_rows=int(args.max_nnz),
        )

        info = get_info(
            exp=exp,
            obs_df=obs_df,
            var_df=var_df,
            X=X,
            obs_tbl=obs_tbl,
            var_tbl=var_tbl,
            gene_joinid=gene_joinid,
            gene_id=gene_id,
            cell_ids=cell_ids,
            x_tbl=x_tbl,
            uri=uri,
            args=args,
        )

        if args.json:
            print(json.dumps(info, indent=2))
        else:
            print_human(info)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


