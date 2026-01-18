# Quick Start

Two common ways to run this project:

- **Without Docker**: Recommended for local development and debugging
- **With Docker**: Recommended for deployment/demo — consistent environment, minimal dependencies

For more detailed usage, see [how-to-use.md](./how-to-use.md).  
For project deployment and server launch, please refer to [deploy.md](./deploy.md).

## Without Docker (Local Development)

### 0) Prerequisites

- Python >= 3.9 (recommend using `uv` for management)
- Node.js 20+
- Go (optional; `make` will attempt to bootstrap if missing)

### 1) Install Dependencies

```bash
make install
```

### 2) Preprocess Data

```bash
mkdir -p data/raw data/dataset_a
wget -O data/raw/input.h5ad https://datasets.cellxgene.cziscience.com/12fbed11-bd9d-43c2-97c1-6a79a2bb1be2.h5ad
make preprocess INPUT=data/raw/input.h5ad OUTPUT=data/dataset_a
```

### 3) Start Development Mode

```bash
make dev SERVER_CONFIG=../config/server.yaml
```

Defaults:

- Frontend: `http://localhost:3000`
- Backend: `http://localhost:8080`


## With Docker (Recommended)

### 0) Prerequisites

- Docker + `docker compose`

### 1) Clone the Repository

```bash
git clone https://github.com/xuzhougeng/atlasmap-sc.git
cd atlasmap-sc
```

### 2) Prepare Data Directory (Host Machine)

`docker-compose.yml` mounts the host's `./data` to `/data` inside containers.

Recommended multi-dataset layout (host paths):

```text
data/
  raw/
    input.h5ad
  dataset_a/
    zarr/metadata.json
    zarr/gene_index.json
    zarr/bins.zarr/...
    soma/... (optional)
```

### 3) Run Preprocessing (Choose One)

#### Option A: Run Preprocessing on Host (Recommended, Faster Iteration)

```bash
mkdir -p data/raw data/dataset_a
wget -O data/raw/input.h5ad https://datasets.cellxgene.cziscience.com/12fbed11-bd9d-43c2-97c1-6a79a2bb1be2.h5ad

# Create venv with uv and install preprocessing dependencies
make preprocess-venv

# Run preprocessing, output to data/dataset_a
make preprocess INPUT=data/raw/input.h5ad OUTPUT=data/dataset_a
```

#### Option B: Run Preprocessing in Container (Most Consistent Environment)

```bash
mkdir -p data/raw data/dataset_a
# Place input.h5ad in data/raw/input.h5ad

docker compose --profile tools run --rm preprocess run \
  -i /data/raw/input.h5ad \
  -o /data/dataset_a \
  -z 8 -g 500
```

### 4) Configure Backend Data Paths (Container Paths)

Since containers see `/data/...`, `config/server.yaml` must use container-internal paths, e.g.:

```yaml
data:
  dataset_a:
    zarr_path: "/data/dataset_a/zarr/bins.zarr"
    soma_path: "/data/dataset_a/soma"  # optional
```

### 5) Start Services

```bash
docker compose up -d --build
```

Access:

- Frontend: `http://localhost:3000`
- Backend health check: `http://localhost:8080/health`

View logs:

```bash
docker compose logs -f server
docker compose logs -f frontend
```

Stop:

```bash
docker compose down
```