# AtlasMap User Guide (how-to-use)

AtlasMap's typical workflow consists of two steps:

1. **Offline Preprocessing (Python)**: Convert `.h5ad` to multi-resolution `Zarr` (`bins.zarr`) + metadata (`metadata.json`/`gene_index.json`).
2. **Online Service (Go + Frontend)**: Go server reads `Zarr` on-demand and renders PNG tiles; frontend displays interactive UMAP via `/api` and `/tiles`.

---

## Directory and Data Conventions

The `data/` directory in the repository is ignored by `.gitignore` by default. We recommend organizing locally as follows:

```text
data/
  raw/                      # Raw input (optional)
    input.h5ad
  preprocessed/             # Preprocessing output
    zarr/
      bins.zarr/            # Zarr store (note: it's a directory)
      bins.X_tsne.zarr/     # (Optional) Zarr store for other coordinate systems (when using --coord-key multiple times)
      metadata.json
      gene_index.json
    soma/                   # TileDBSOMA storage (optional, enabled by default)
      experiment.soma/      # Complete single-cell data
```

**Important**: The backend reads the `bins.zarr` directory specified by `data.zarr_path`, and looks for `metadata.json` and `gene_index.json` in its **parent directory** (i.e., `.../zarr/metadata.json`).

---

## Method A: Local Run (Development/Debugging)

### 0) Dependencies

- Go 1.22+
- Node.js 20+
- Python 3.9+

Optional: Use `make install` to install all dependencies at once (installs Python/Go/frontend dependencies separately).

### 0.5) View H5AD Information (Optional)

Before running preprocessing, you can view basic information about the `.h5ad` (cells/genes count, layers, coordinate keys, etc.):

```bash
python scripts/h5ad_info.py -i data.h5ad
```

Example output:

```
File: data.h5ad
Cells (n_obs): 10,000
Genes (n_vars): 2,000

X matrix: [10000, 2000] (...)
Layers: None

Embeddings/Coordinates (obsm) [2]:
  X_pca: [10000, 50]
  X_umap: [10000, 2] <-- X_umap

obs columns (5): cell_type, batch, ...
var columns (2): gene_name, highly_variable
...
```

Add `--json` for JSON format output.

### 1) Preprocess Data (Python)

Simplest (using Makefile):

```bash
make preprocess INPUT=path/to/input.h5ad
```

Or use CLI directly (more flexible):

```bash
cd preprocessing
pip install -e .
atlasmap-preprocess run -i path/to/input.h5ad -o ../data/preprocessed -z 10 -g 500
```

Common parameters (consistent with code, see `preprocessing/atlasmap_preprocess/cli.py`):

- `--coord-key`: Coordinate key in `.h5ad` `adata.obsm[...]` (can be specified multiple times to generate multiple coordinate system tiles)
- `--umap-key`: Default coordinate key (legacy parameter name; used to specify default/main coordinate system, defaults to `X_umap`)
- `--category/-c`: Category fields to write (from `adata.obs`, can be specified multiple times)
- `--zoom-levels/-z`: Number of zoom levels (default 8, starting from 0 up to the specified value)
- `--n-genes/-g`: Number of pre-aggregated genes (default 500)
- `--all-expressed/-a`: Use all expressed genes instead of top N genes
- `--min-cells`: Minimum cells a gene must be expressed in (default 3, only used with `--all-expressed`)
- `--no-soma`: Disable TileDBSOMA storage (enabled by default, stores complete expression matrix for arbitrary gene queries)

If you need to specify marker genes/more advanced parameters: generate config file first then run:

```bash
cd preprocessing
atlasmap-preprocess init-config -o config.yaml
# Edit config.yaml (e.g., marker_genes / category_columns etc.)
atlasmap-preprocess from-config -c config.yaml
```

### 1.5) Validate Preprocessing Results (Optional)

Use the `visualize` command to generate static images and validate preprocessing output:

> Tip: Here `zoom` represents "binning resolution" (higher = finer). With default `--zoom-levels 8`:
> - `zoom=0`: 1×1 bin, usually shows only 1 point/block (very coarse)
> - `zoom=7`: 128×128 bins, shows complete UMAP structure (recommended as default check/display level)

```bash
# Basic usage (defaults to drawing zoom 3,5,7 for comparison; randomly selects 3 genes)
atlasmap-preprocess visualize -i ../data/preprocessed/zarr -o ../data/figures

# Recommended: only view highest resolution (closest to global effect when opening frontend)
atlasmap-preprocess visualize -i ../data/preprocessed/zarr -o ../data/figures -z 7

# Specify zoom levels and genes
atlasmap-preprocess visualize -i ../data/preprocessed/zarr -o ../data/figures -z 3,5,7 -g CD3D -g CD8A
```

Output structure:
```
figures/
├── category/          # Colored by category
│   ├── cell_type_zoom_7.png
│   └── cell_type_multi_zoom.png
└── expression/        # Colored by gene expression
    ├── GENE1_zoom_7.png
    └── GENE1_multi_zoom.png
```

### 1.6) Build BLASTP Database (Optional, for Multi-species Search)

If you need to use BLASTP for cross-dataset homologous gene search, you need to build a protein sequence database for each dataset.

#### Prerequisites

- Install NCBI BLAST+ toolkit (includes `makeblastdb` and `blastp` commands)

```bash
# Ubuntu/Debian
sudo apt-get install ncbi-blast+

# macOS
brew install blast

# Or use conda
conda install -c bioconda blast
```

#### Prepare Protein Sequence File

Protein sequence file must be in FASTA format:

```fasta
>Afi_001234
MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPT
LVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDT
>Afi_005678
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNP
KVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHH
>Afi_009999
MALTVRIQAACLLLLLLAAALCPTHGAAAGGSTKEVVEEAENGRDAPANGNAENEENRQ
```

**Important Notes**:
- Sequence ID **must** match the `gene_id` in the dataset
- Sequence ID (part after `>`) will be used as the `sseqid` field in BLASTP results
- This ID must be found in the dataset's `gene_index.json`, otherwise the frontend cannot navigate to expression view
- Recommend using the gene ID from the dataset as sequence ID, not gene name or other identifiers

#### Build Database

```bash
# Basic usage
makeblastdb -in proteins.fasta -dbtype prot -out /path/to/blast/db/proteins

# Example: build database for human dataset
mkdir -p data/preprocessed/blast
makeblastdb -in data/raw/human_proteins.fasta \
            -dbtype prot \
            -out data/preprocessed/blast/human_proteins
```

This generates the following files:
```
data/preprocessed/blast/
├── human_proteins.phr    # Header info
├── human_proteins.pin    # Index
├── human_proteins.psq    # Sequence data
└── human_proteins.pal    # (Optional) Alias
```

#### Specify Database in Configuration File

Edit `config/server.yaml`, add `blastp_path` for each dataset (**note: path prefix without extension**):

```yaml
data:
  blood:
    zarr_path: "/data/preprocessed/zarr/bins.zarr"
    soma_path: "/data/preprocessed/soma"
    blastp_path: "/data/preprocessed/blast/blood_proteins"  # Without .phr/.pin/.psq
  
  liver:
    zarr_path: "/data/liver/zarr/bins.zarr"
    soma_path: "/data/liver/soma"
    blastp_path: "/data/preprocessed/blast/liver_proteins"
```

#### Multi-dataset Shared Database

If multiple datasets' genes come from the same species or gene set, they can share the same BLASTP database:

```yaml
data:
  tissue_a:
    zarr_path: "/data/tissue_a/zarr/bins.zarr"
    blastp_path: "/data/shared/blast/species_proteins"
  
  tissue_b:
    zarr_path: "/data/tissue_b/zarr/bins.zarr"
    blastp_path: "/data/shared/blast/species_proteins"  # Shared database
```

The backend automatically deduplicates execution (runs `blastp` only once), but results are expanded to all datasets using that database.

#### Verify Database

```bash
# Check if database was built successfully
blastdbcmd -db /data/preprocessed/blast/blood_proteins -info

# Test query (using example sequence)
echo ">test
MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGK" | \
  blastp -db /data/preprocessed/blast/blood_proteins \
         -outfmt "6 sseqid pident length evalue bitscore" \
         -evalue 1e-5 \
         -max_target_seqs 5
```

After configuration, a search button will appear in the frontend toolbar. Click to enter the BLASTP search page.

### 2) Start Backend (Go)

The backend Go module is in `server/`. We recommend starting from the `server/` directory with explicit config file:

```bash
cd server
go run ./cmd/server -config ../config/server.yaml
```

If you need to enable **TileDB-SOMA (arbitrary gene × cell expression queries)**:

Step 1: Install TileDB C library (`libtiledb` / headers)

```bash
conda install -c conda-forge -y tiledb
```

> Version compatibility: TileDB-Go and TileDB core need to match. This repository's `server/go.mod` uses `github.com/TileDB-Inc/TileDB-Go v0.38.0`, corresponding to TileDB `2.29.x` (e.g., `2.29.1`).

Step 2: Build/run with build tag (requires CGO):

```bash
cd server

export CGO_ENABLED=1
export CGO_CFLAGS="-I$CONDA_PREFIX/include"
export CGO_LDFLAGS="-L$CONDA_PREFIX/lib -ltiledb -Wl,-rpath,$CONDA_PREFIX/lib"

go run -tags soma ./cmd/server -config ../config/server.yaml
```

If you see many `deprecated` compiler warnings: they don't affect operation. To silence them:

```bash
export CGO_CFLAGS="-I$CONDA_PREFIX/include -Wno-deprecated-declarations"
```

Key configuration items in `config/server.yaml`:

- `server.port`: Default `8080`
- `data.zarr_path`: Points to `.../zarr/bins.zarr` (recommend using relative or absolute paths locally)

**Multi-dataset Configuration**: You can define multiple datasets under `data`, the first one is the default:

```yaml
data:
  pbmc:
    zarr_path: "/data/pbmc/zarr/bins.zarr"
    soma_path: "/data/pbmc/soma"
  liver:
    zarr_path: "/data/liver/zarr/bins.zarr"
    soma_path: "/data/liver/soma"
```

The frontend will display a dataset selection dropdown, and API can access specific datasets via `/d/{dataset}/api/...`.

If preprocessing generated multiple coordinate systems via `--coord-key`, the frontend will also display a `Coord` dropdown; you can also append `?coord=X_tsne` to any tiles/API request to switch coordinate systems.

Self-check after startup (note: first get dataset ID via `/api/datasets`):

```bash
curl http://localhost:8080/health
curl http://localhost:8080/api/datasets
curl http://localhost:8080/d/{dataset}/api/metadata  # Replace {dataset} with actual ID
```

### 3) Start Frontend (Vite)

```bash
cd frontend
npm ci
npm run dev
```

Default opens: `http://localhost:3000`  
During local development, Vite proxies `/api` and `/d` to `http://localhost:8080` (see `frontend/vite.config.ts`).

---

## Method B: Docker Compose (Deployment/Demo)

### 1) Prepare Configuration and Data Paths

`docker-compose.yml` mounts host's `./data` to container's `/data`, so container should use:

- `data.zarr_path: "/data/preprocessed/zarr/bins.zarr"`

If you're using `config/server.yaml` directly, ensure `data.zarr_path` **is valid for the container path**; otherwise the backend container will fail to start because it can't find the Zarr data.

### 2) Build and Start

```bash
docker compose up -d --build
```

Access:

- Frontend: `http://localhost:3000`
- Backend health check: `http://localhost:8080/health`

### 3) (Optional) Run Preprocessing in Docker

The `preprocess` service in `docker-compose.yml` is under the `tools` profile by default, needs explicit enabling:

```bash
mkdir -p data/raw data/preprocessed
# Put input.h5ad in data/raw/input.h5ad

docker compose --profile tools run --rm preprocess run \
  -i /data/raw/input.h5ad \
  -o /data/preprocessed \
  -z 8 -g 500
```

After preprocessing completes, ensure `data/preprocessed/zarr/bins.zarr` and sibling `metadata.json`/`gene_index.json` are generated, then start/restart the backend service.

---

## Troubleshooting

### Backend startup error: `Failed to initialize Zarr reader`

Usually the `data.zarr_path` in `config/server.yaml` is incorrect or preprocessing artifacts are incomplete:

- `data.zarr_path` must point to `.../zarr/bins.zarr` (directory)
- `.../zarr/metadata.json` and `.../zarr/gene_index.json` must exist

### Preprocessing error: `Coordinate key 'X_umap' not found`

The input `.h5ad` does not contain the corresponding `adata.obsm[...]` coordinate key. Use `--coord-key` (or legacy parameter `--umap-key`) to specify the correct key name.

### Gene expression tile is always blank/gene not found

Expression tiles can only render **pre-aggregated genes**; first check available genes via `GET /d/{dataset}/api/genes`, or increase `--n-genes` during preprocessing, or use `from-config` to add marker genes if needed.
