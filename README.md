# SOMA-Tiles

High-performance single-cell visualization system for 10M+ cells.

Detailed usage/deployment guide: `doc/how-to-use.md`.

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    Offline Preprocessing                        │
│  (Run on machine with >16GB RAM)                                │
│                                                                 │
│  H5AD ──► Multi-resolution UMAP Bins ──► TileDBSOMA Store      │
│           (Zarr, 8 zoom levels)          (Full expression)      │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Online Service (<2GB RAM)                    │
│                                                                 │
│  Go Server ──► Zarr Column Slices ──► PNG Tile Rendering       │
│            ──► Bin Aggregation      ──► Legend JSON             │
│  (SOMA for arbitrary gene queries)                              │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Web Frontend                                 │
│                                                                 │
│  Leaflet Tiles + Annotation/Gene Controls + Filter UI          │
└─────────────────────────────────────────────────────────────────┘
```

## Quick Start

### Prerequisites

- Go 1.22+
- Node.js 20+
- Python 3.9+
- Docker (optional)

### Installation

```bash
# Clone the repository
git clone https://github.com/your-org/soma-tiles.git
cd soma-tiles

# Install all dependencies
make install
```

### Preprocessing

```bash
# Run preprocessing on your H5AD file
make preprocess INPUT=path/to/data.h5ad

# Or use the CLI directly
cd preprocessing
pip install -e .
soma-preprocess run -i data.h5ad -o ../data/preprocessed -g 500 -z 8
```

### Development

```bash
# Start development servers (Go backend + Vite frontend)
make dev

# Or start separately:
# Terminal 1: Go server
cd server && go run ./cmd/server -config ../config/server.yaml

# Terminal 2: Frontend
cd frontend && npm run dev
```

### Docker

```bash
# Build and start all services
docker-compose up -d

# Run preprocessing in Docker
docker-compose run --rm preprocess run \
    --input /data/raw/input.h5ad \
    --output /data/preprocessed
```

## Project Structure

```
soma-tiles/
├── preprocessing/          # Python preprocessing pipeline
│   ├── soma_tiles_preprocess/
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

## API Endpoints

### Tiles
- `GET /tiles/{z}/{x}/{y}.png` - Base tile
- `GET /tiles/{z}/{x}/{y}/expression/{gene}.png` - Expression colored tile

### Metadata
- `GET /api/metadata` - Dataset metadata
- `GET /api/genes` - List of pre-aggregated genes
- `GET /api/categories` - Category information

## Configuration

See `config/server.yaml` for server configuration options.

### Single Dataset (Legacy Format)

```yaml
server:
  port: 8080
  cors_origins:
    - "http://localhost:3000"

data:
  zarr_path: "/data/preprocessed/zarr/bins.zarr"
  soma_path: "/data/preprocessed/soma"

cache:
  tile_size_mb: 512
  tile_ttl_minutes: 10

render:
  tile_size: 256
  default_colormap: viridis
```

### Multiple Datasets

You can configure multiple datasets. The first dataset becomes the default:

```yaml
server:
  port: 8080
  cors_origins:
    - "http://localhost:3000"

data:
  pbmc:
    zarr_path: "/data/pbmc/zarr/bins.zarr"
    soma_path: "/data/pbmc/soma"
  liver:
    zarr_path: "/data/liver/zarr/bins.zarr"
    soma_path: "/data/liver/soma"

cache:
  tile_size_mb: 512
  tile_ttl_minutes: 10

render:
  tile_size: 256
  default_colormap: viridis
```

When multiple datasets are configured:
- The frontend displays a dataset selector dropdown
- Each dataset is accessible via `/d/{dataset}/api/...` and `/d/{dataset}/tiles/...`
- Legacy routes (`/api/...`, `/tiles/...`) point to the default dataset

## License

MIT
