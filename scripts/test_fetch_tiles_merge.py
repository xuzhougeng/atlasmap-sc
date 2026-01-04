#!/usr/bin/env python3
"""
Test script: fetch server tiles at a given zoom and stitch into one image.

Supports:
  - category coloring (e.g. cell_type): /d/{dataset}/tiles/{z}/{x}/{y}/category/{column}.png
  - expression coloring (random gene): /d/{dataset}/tiles/{z}/{x}/{y}/expression/{gene}.png?colormap=viridis

Notes:
  - For large zoom/bounds, the stitched image can be huge. This script defaults to downsampling
    so the output's max dimension is <= --max-dim (default: 8192).
  - Uses /api/datasets to auto-detect the default dataset if --dataset is not specified.
"""

from __future__ import annotations

import argparse
import io
import json
import math
import random
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple


try:
    from PIL import Image
except Exception as e:  # pragma: no cover
    raise SystemExit(
        "Missing dependency: Pillow. Install with: pip install pillow\n"
        f"Import error: {e}"
    )


@dataclass(frozen=True)
class Meta:
    tile_size: int
    zoom_levels: int
    bounds_max: float


def http_get_bytes(url: str, timeout_s: float) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "soma-tiles-test/1.0"})
    with urllib.request.urlopen(req, timeout=timeout_s) as resp:
        return resp.read()


def http_get_json(url: str, timeout_s: float) -> Any:
    data = http_get_bytes(url, timeout_s=timeout_s)
    return json.loads(data.decode("utf-8"))


def fetch_default_dataset(base_url: str, timeout_s: float) -> str:
    """Fetch the default dataset ID from /api/datasets."""
    resp = http_get_json(f"{base_url}/api/datasets", timeout_s=timeout_s)
    default = resp.get("default")
    if not default:
        datasets = resp.get("datasets") or []
        if datasets:
            default = datasets[0].get("id") if isinstance(datasets[0], dict) else datasets[0]
    if not default:
        raise RuntimeError("No datasets available from /api/datasets")
    return default


def fetch_metadata(base_url: str, dataset: str, timeout_s: float) -> Meta:
    md = http_get_json(f"{base_url}/d/{dataset}/api/metadata", timeout_s=timeout_s)
    tile_size = int(md.get("tile_size", 256))
    zoom_levels = int(md["zoom_levels"])
    bounds = md.get("bounds") or {}
    bounds_max = float(bounds.get("max_x", 256))
    return Meta(tile_size=tile_size, zoom_levels=zoom_levels, bounds_max=bounds_max)


def pick_gene(base_url: str, dataset: str, timeout_s: float, gene: Optional[str]) -> str:
    if gene and gene != "random":
        return gene
    genes_resp = http_get_json(f"{base_url}/d/{dataset}/api/genes", timeout_s=timeout_s)
    genes = genes_resp.get("genes") or []
    if not genes:
        raise RuntimeError(f"No genes returned by /d/{dataset}/api/genes; cannot pick random gene.")
    return random.choice(list(genes))


def tiles_per_axis(bounds_max: float, tile_size: int, zoom: int) -> int:
    # CRS.Simple in Leaflet: pixel coords scale by 2^zoom; tile index is floor(pixel / tile_size).
    # If coordinates are in [0, bounds_max], then pixel span is bounds_max * 2^zoom.
    span_px = bounds_max * (2**zoom)
    return max(1, int(math.ceil(span_px / tile_size)))


def build_tile_url(
    base_url: str,
    dataset: str,
    zoom: int,
    x: int,
    y: int,
    mode: str,
    category: str,
    gene: str,
    colormap: str,
) -> str:
    if mode == "category":
        return f"{base_url}/d/{dataset}/tiles/{zoom}/{x}/{y}/category/{urllib.parse.quote(category)}.png"
    if mode == "expression":
        qp = urllib.parse.urlencode({"colormap": colormap})
        return f"{base_url}/d/{dataset}/tiles/{zoom}/{x}/{y}/expression/{urllib.parse.quote(gene)}.png?{qp}"
    raise ValueError(f"Unknown mode: {mode}")


def decode_png(png_bytes: bytes) -> Image.Image:
    img = Image.open(io.BytesIO(png_bytes))
    img.load()
    return img.convert("RGBA")


def is_fully_transparent(img: Image.Image) -> bool:
    if img.mode != "RGBA":
        img = img.convert("RGBA")
    alpha = img.getchannel("A")
    return alpha.getextrema() == (0, 0)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--server", default="http://localhost:8080", help="Base URL of Go server")
    p.add_argument("--dataset", default="", help="Dataset ID (auto-detect default if empty)")
    p.add_argument("--zoom", type=int, default=7, help="Zoom level to fetch (e.g. 7)")
    p.add_argument(
        "--mode",
        choices=["category", "expression"],
        default="category",
        help="Coloring mode",
    )
    p.add_argument("--category", default="cell_type", help="Category column name (for --mode category)")
    p.add_argument("--gene", default="random", help="Gene name or 'random' (for --mode expression)")
    p.add_argument("--colormap", default="viridis", help="Colormap (for --mode expression)")
    p.add_argument("--timeout", type=float, default=20.0, help="HTTP timeout seconds")
    p.add_argument("--threads", type=int, default=16, help="Parallel download threads")
    p.add_argument(
        "--max-dim",
        type=int,
        default=8192,
        help="Downsample so stitched image max(width,height) <= this. Set 0 to disable.",
    )
    p.add_argument(
        "--skip-empty",
        action="store_true",
        help="Skip fully-transparent tiles (keeps the output black/transparent there).",
    )
    p.add_argument("--x0", type=int, default=0, help="Start tile x (inclusive)")
    p.add_argument("--x1", type=int, default=0, help="End tile x (inclusive); -1 means full range")
    p.add_argument("--y0", type=int, default=0, help="Start tile y (inclusive)")
    p.add_argument("--y1", type=int, default=0, help="End tile y (inclusive); -1 means full range")
    p.add_argument(
        "--full-grid",
        action="store_true",
        help="Fetch and stitch the full tile grid for the zoom (sets x1=y1=-1).",
    )
    p.add_argument(
        "--out",
        default="stitched.png",
        help="Output file path (PNG). If relative, relative to repo root.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    base_url = args.server.rstrip("/")

    # Resolve dataset ID
    if args.dataset:
        dataset = args.dataset
    else:
        dataset = fetch_default_dataset(base_url, timeout_s=args.timeout)
    print(f"[info] dataset={dataset}")

    meta = fetch_metadata(base_url, dataset, timeout_s=args.timeout)
    if args.zoom < 0 or args.zoom >= meta.zoom_levels:
        raise SystemExit(f"--zoom out of range: {args.zoom} (zoom_levels={meta.zoom_levels})")

    chosen_gene = ""
    if args.mode == "expression":
        chosen_gene = pick_gene(base_url, dataset, timeout_s=args.timeout, gene=args.gene)
        print(f"[info] mode=expression gene={chosen_gene} colormap={args.colormap}")
    else:
        print(f"[info] mode=category column={args.category}")

    n = tiles_per_axis(meta.bounds_max, meta.tile_size, args.zoom)
    if args.full_grid:
        args.x1 = -1
        args.y1 = -1
    x0 = max(0, args.x0)
    y0 = max(0, args.y0)
    x1 = n - 1 if args.x1 < 0 else min(args.x1, n - 1)
    y1 = n - 1 if args.y1 < 0 else min(args.y1, n - 1)

    tiles_x = x1 - x0 + 1
    tiles_y = y1 - y0 + 1
    if tiles_x <= 0 or tiles_y <= 0:
        raise SystemExit("Empty tile range; check --x0/--x1/--y0/--y1")

    full_w = tiles_x * meta.tile_size
    full_h = tiles_y * meta.tile_size

    scale = 1.0
    if args.max_dim and args.max_dim > 0:
        scale = min(1.0, float(args.max_dim) / float(max(full_w, full_h)))
    out_w = max(1, int(round(full_w * scale)))
    out_h = max(1, int(round(full_h * scale)))
    tile_out = max(1, int(round(meta.tile_size * scale)))

    print(
        "[info] metadata:",
        f"tile_size={meta.tile_size} zoom_levels={meta.zoom_levels} bounds_max={meta.bounds_max}",
    )
    print(
        "[info] grid:",
        f"zoom={args.zoom} tiles_per_axis≈{n} range_x=[{x0},{x1}] range_y=[{y0},{y1}]",
    )
    print(f"[info] output: {args.out} size={out_w}x{out_h} scale={scale:.4f} tile_out={tile_out}")

    canvas = Image.new("RGBA", (out_w, out_h), (0, 0, 0, 0))

    def fetch_one(x: int, y: int) -> Tuple[int, int, bytes]:
        url = build_tile_url(
            base_url=base_url,
            dataset=dataset,
            zoom=args.zoom,
            x=x,
            y=y,
            mode=args.mode,
            category=args.category,
            gene=chosen_gene,
            colormap=args.colormap,
        )
        try:
            data = http_get_bytes(url, timeout_s=args.timeout)
            return x, y, data
        except urllib.error.HTTPError as e:
            # For debugging: keep a marker tile if server ever returns non-200.
            raise RuntimeError(f"HTTP {e.code} for {url}: {e.read().decode('utf-8', 'ignore')}") from e
        except Exception as e:
            raise RuntimeError(f"Failed to fetch {url}: {e}") from e

    total = tiles_x * tiles_y
    done = 0
    t0 = time.time()

    with ThreadPoolExecutor(max_workers=max(1, int(args.threads))) as ex:
        futures = []
        for yy in range(y0, y1 + 1):
            for xx in range(x0, x1 + 1):
                futures.append(ex.submit(fetch_one, xx, yy))

        for fut in as_completed(futures):
            x, y, png = fut.result()
            img = decode_png(png)
            if args.skip_empty and is_fully_transparent(img):
                pass
            else:
                if scale != 1.0:
                    img = img.resize((tile_out, tile_out), resample=Image.Resampling.BILINEAR)
                ox = (x - x0) * tile_out
                oy = (y - y0) * tile_out
                canvas.paste(img, (ox, oy), img)

            done += 1
            if done % 200 == 0 or done == total:
                dt = time.time() - t0
                rate = done / dt if dt > 0 else 0.0
                print(f"[progress] {done}/{total} tiles ({rate:.1f} tiles/s)")

    canvas.save(args.out)
    print("[done] wrote", args.out)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\n[abort] interrupted", file=sys.stderr)
        raise SystemExit(130)


