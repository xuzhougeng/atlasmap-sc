# **AtlasMap: enabling low-cost, map-style exploration of million-cell single-cell atlases**

High-performance single-cell visualization system for 10M+ cells.

![](doc/images/Architecture.png)

Please refer to the following guides:

- [Quick Start](./doc/quick-start.md)
- [Detailed Tutorial](./doc/how-to-use.md)
- [Deployment](./doc/deploy.md)

You can view the AtlasMap browser in action at <https://explorer.plantcellatlas.com/>.

If you use **AtlasMap** in your work, please cite:

> **AtlasMap: enabling low-cost, map-style exploration of million-cell single-cell atlases**  
> Zhou-Geng Xu, Hao-Chen Xue, Guihui Qin, Jia-Wei Wang  
> doi: [10.64898/2026.01.14.699595](https://doi.org/10.64898/2026.01.14.699595)


## Project Structure

```
atlasmap-sc/
├── preprocessing/          # Python preprocessing pipeline
│   ├── atlasmap_preprocess/
│   │   ├── pipeline.py     # Main pipeline
│   │   ├── binning/        # Quadtree binning
│   │   └── io/             # Zarr & SOMA I/O
│   └── pyproject.toml
│
├── server/                 # Go backend
│   ├── cmd/server/         # Entry point
│   └── internal/
│       ├── api/            # REST API
│       ├── render/         # Tile rendering
│       └── data/zarr/      # Zarr reader
│
├── frontend/               # Vanilla JS/TS frontend
│   ├── src/
│   │   ├── map/            # Leaflet integration
│   │   ├── components/     # UI components
│   │   └── api/            # API client
│   └── index.html
│
├── config/                 # Configuration files
├── data/                   # Data directory (gitignored)
│   └── preprocessed/
│       ├── zarr/           # Bin aggregations
│       └── soma/           # Full expression (TileDBSOMA)
└── docker-compose.yml
```

## Architecture

AtlasMap is a two-stage pipeline: **offline preprocessing** turns an `.h5ad` into a tile-friendly,
multi-zoom bin store, and an **online Go server** renders Leaflet-compatible PNG tiles on demand.
A static TypeScript frontend consumes the JSON + tile endpoints.

```text
┌──────────────────────────────────────────────┐
│               Offline (Python)               │
│  preprocessing/ (atlasmap-preprocess CLI)    │
│  - normalize UMAP coords to [0,256)          │
│  - select preaggregated genes                │
│  - quadtree binning for zoom=0..Z-1          │
│  - per-bin aggregation (sparse non-empty)    │
└───────────────────────────┬──────────────────┘
                            │ writes
                            ▼
┌────────────────────────────────────────────────────────────────────────┐
│ Data artifacts (per dataset)                                            │
│                                                                        │
│  zarr/metadata.json + zarr/gene_index.json                              │
│  zarr/bins.zarr  (Zarr v3, zstd chunks)                                 │
│    zoom_{z}/bin_coords        int32  [N,2]   (bin_x,bin_y)              │
│    zoom_{z}/cell_count        uint32 [N]                                  │
│    zoom_{z}/expression_mean   float32 [N,G] (preaggregated genes)       │
│    zoom_{z}/expression_max    float32 [N,G]                              │
│    zoom_{z}/category_counts/* uint32 [N,C] (per column)                 │
│    zoom_{z}/cell_ids/*        ragged cell indices (optional; off by default) │
│                                                                        │
│  soma/experiment.soma (optional; TileDB-SOMA full matrix + obs/var)     │
└───────────────────────────┬────────────────────────────────────────────┘
                            │ read
                            ▼
┌────────────────────────────────────────────────────────────────────────┐
│ Online (Go)                                                           │
│  server/ (chi router)                                                 │
│  - multi-dataset registry (/api/datasets, /d/{dataset}/...)            │
│  - Zarr reader (loads sparse bins + per-gene/category vectors)         │
│  - tile service: render from highest zoom and slice per {z,x,y}        │
│  - PNG renderer (256x256; base/category/expression)                    │
│  - in-memory tile cache (BigCache; dataset-aware keys)                 │
│  - (optional, -tags soma) SOMA endpoints + DE job runner (SQLite)      │
└───────────────────────────┬────────────────────────────────────────────┘
                            │ HTTP (PNG tiles + JSON)
                            ▼
┌────────────────────────────────────────────────────────────────────────┐
│ Web Frontend (TS + Vite + Leaflet)                                     │
│  frontend/                                                            │
│  - dataset selector + metadata bootstrap                               │
│  - Leaflet tile layers: base/category/expression                        │
│  - gene/category controls + filters + (optional) DE UI                 │
└────────────────────────────────────────────────────────────────────────┘
```

## License

MIT
