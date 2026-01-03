# SOMA-Tiles

High-performance single-cell visualization system for 10M+ cells.

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    Offline Preprocessing                        │
│  (Run on machine with >16GB RAM)                                │
│                                                                 │
│  H5AD ──► TileDB-SOMA ──► Multi-resolution UMAP Bins           │
│           (X, obs, var)   (Zarr, 8 zoom levels)                │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Online Service (<2GB RAM)                    │
│                                                                 │
│  Go Server ──► SOMA Column Slices ──► PNG Tile Rendering       │
│            ──► Bin Aggregation       ──► Legend JSON            │
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
cd server && go run ./cmd/server

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
│   │   └── io/             # Zarr I/O
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

See `config/server.yaml` for server configuration options:

```yaml
server:
  port: 8080
  cors_origins:
    - "http://localhost:3000"

data:
  zarr_path: "/data/preprocessed/zarr/bins.zarr"

cache:
  tile_size_mb: 512
  tile_ttl_minutes: 10

render:
  tile_size: 256
  default_colormap: viridis
```

## License

MIT
