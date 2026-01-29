## API Reference (Backend)

Code location: `server/internal/api/routes.go`

General notes:
- All API and Tiles endpoints require the `/d/{dataset}` prefix, e.g., `/d/pbmc/api/metadata`, `/d/pbmc/tiles/0/0/0.png`
- `{dataset}` is the dataset ID, which can be obtained via `GET /api/datasets`

### Health

- `GET /health`: Health check endpoint

### Datasets (Global)

- `GET /api/datasets`: List available datasets and default dataset (this endpoint does not require `/d/{dataset}` prefix)
  - Response includes `has_blastp` field indicating whether BLASTP database is configured

### Gene Lookup (Global)

- `GET /api/gene_lookup?gene_id={gene_id}`: Find datasets containing the specified gene (for cross-dataset navigation)
  - Response: `{ "gene_id": "...", "datasets": ["dataset1", "dataset2", ...] }`
  - When the URL only provides a `gene` parameter, the frontend automatically calls this endpoint to locate the dataset

### API (JSON)

- `GET /d/{dataset}/api/metadata`: Dataset metadata (required for frontend initialization)
- `GET /d/{dataset}/api/stats`: Statistics (e.g., n_cells, n_genes, zoom_levels, dataset_name)
- `GET /d/{dataset}/api/genes`: Pre-aggregated gene list
- `GET /d/{dataset}/api/genes/{gene}`: Gene information (e.g., index/preaggregated)
- `GET /d/{dataset}/api/genes/{gene}/stats?zoom={z}`: Gene statistics (`zoom` is optional, defaults to 0)
- `GET /d/{dataset}/api/genes/{gene}/bins?threshold=&offset=&limit=`: Query bins expressing the gene (limit defaults to 100, max 1000)
- `GET /d/{dataset}/api/genes/{gene}/category/{column}/means`: Mean expression of the gene by category (for category comparison)
- `GET /d/{dataset}/api/categories`: Available category columns and values
- `GET /d/{dataset}/api/categories/{column}/colors`: Category color mapping
- `GET /d/{dataset}/api/categories/{column}/centroids`: Category centroids (weighted by `cell_count` `(x,y)`)
- `GET /d/{dataset}/api/categories/{column}/legend`: Category legend

### SOMA (Experiment-level Expression Query, TileDB-SOMA)

> Note: This endpoint requires the backend to be built with `-tags soma`, and the runtime environment must have the TileDB C library installed (`libtiledb` + headers); otherwise, it returns `501` (not implemented).

- `GET /d/{dataset}/api/soma/expression?gene={gene}&cells={c1,c2,...}&mode={sparse|dense}`
  - `gene`: Gene name (`ms/RNA/var.gene_id`; note: not `gene_name`)
  - `cells`: Cell `soma_joinid` (integers, comma-separated; not `cell_id` strings)
  - `mode`:
    - `sparse` (default): Returns non-zero items `items=[{cell_joinid,value}, ...]`
    - `dense`: Returns values aligned with input `cells` `values=[...]` (missing values are 0)
  - Limit: Maximum 10000 `cells` (exceeding returns 400)

Example:

```bash
curl -sS 'http://localhost:8080/d/retina/api/soma/expression?gene=ENSMICG00000038116&cells=95,118,143&mode=dense'
```

Status codes (current implementation):

- 404: Dataset does not have `soma_path` configured or experiment does not exist
- 501: Server not built with `-tags soma`
- 400: Parameter error / gene does not exist / TileDB query failed (error body is plain text)

### SOMA Differential Expression Analysis (DE Job, TileDB-SOMA)

> Also requires backend built with `-tags soma`.

Differential expression analysis is an asynchronous task: after submission, a `job_id` is returned; the frontend polls for status and retrieves results when complete.

#### Submit Job

- `POST /d/{dataset}/api/soma/de/jobs`
  - Body (JSON):
    - `groupby`: string (obs column name, e.g., `cell_type`)
    - `group1`: string[] (group 1 categories, at least one)
    - `group2`: string[] (group 2 categories; empty means one-vs-rest)
    - `tests`: string[] (`["ttest","ranksum"]`; defaults to both)
    - `max_cells_per_group`: int (max cells sampled per group; default 2000, max 20000)
    - `seed`: int (sampling random seed; default 0)
    - `limit`: int (return top N genes; default 50, max 500)
  - Response: `{ "job_id": "...", "status": "queued" }`

Example:

```bash
curl -X POST 'http://localhost:8080/d/retina/api/soma/de/jobs' \
  -H 'Content-Type: application/json' \
  -d '{"groupby":"cell_type","group1":["retinal rod cell"],"group2":[]}'
```

#### Query Status

- `GET /d/{dataset}/api/soma/de/jobs/{job_id}`
  - Response: `{ job_id, status, created_at, started_at, finished_at, progress:{phase,done,total}, error }`
  - `status`: `queued|running|completed|failed|cancelled`

Example:

```bash
curl -sS 'http://localhost:8080/d/retina/api/soma/de/jobs/8d72e092e14fc785' | jq
```

#### Retrieve Results

- `GET /d/{dataset}/api/soma/de/jobs/{job_id}/result?offset=&limit=&order_by=`
  - Only available when `status=completed`
  - Query parameters:
    - `offset`: Starting position (default 0)
    - `limit`: Number of results (default 50, max 500)
    - `order_by`: Sort order, options: `fdr_ranksum` (default), `fdr_ttest`, `p_ranksum`, `p_ttest`, `abs_log2fc`
  - Response:
    ```json
    {
      "params": {...},
      "n1": 500, "n2": 500,
      "total": 20000,
      "offset": 0, "limit": 50,
      "order_by": "fdr_ranksum",
      "items": [
        {"gene":"...", "gene_joinid":123, "mean1":0.5, "mean2":0.1, "pct1":0.8, "pct2":0.2,
         "log2fc":2.3, "p_ttest":0.001, "fdr_ttest":0.01, "p_ranksum":0.002, "fdr_ranksum":0.02},
        ...
      ]
    }
    ```

Example:

```bash
curl 'http://localhost:8080/d/retina/api/soma/de/jobs/8d72e092e14fc785/result | jq
```

#### Cancel Job

- `DELETE /d/{dataset}/api/soma/de/jobs/{job_id}`

#### obs Metadata Query (for Frontend Dropdowns)

- `GET /d/{dataset}/api/soma/obs/columns`: Returns available obs column names
- `GET /d/{dataset}/api/soma/obs/{column}/values`: Returns unique values for the column

#### Statistical Methods

- **t-test**: Welch t-test (does not assume equal variance), p-value calculated using Student-t distribution CDF
- **ranksum**: Mann-Whitney U test (Wilcoxon rank-sum), with tie correction, normal approximation p-value
- **FDR**: Benjamini-Hochberg multiple testing correction

#### Persistence and Limits

- **SQLite Persistence**: Job status and results are stored in SQLite database (`de.sqlite_path`), historical results can be queried after server restart
- **Restart Recovery**: On server restart, running jobs are marked as `failed`; queued jobs are re-queued
- **TTL Cleanup**: Completed jobs are retained for **7 days** by default (`de.retention_days`), then automatically cleaned up
- Max 20000 cells per group (`max_cells_per_group`)
- Default max 1 concurrent DE job (configurable via `de.max_concurrent`), others are queued

### BLASTP Protein Sequence Search (Global)

> Used for cross-dataset homologous gene search. BLASTP searches the protein databases configured for each dataset and returns a hit list.

#### Submit BLASTP Job

- `POST /api/blastp/jobs`
  - Body (JSON):
    - `sequence`: string (protein sequence, FASTA format or plain sequence)
    - `max_hits`: int (max hits per database, default 10, max 100)
    - `evalue`: float (E-value threshold, default 1e-5)
    - `datasets`: string[] (optional, limit search to specific datasets; empty searches all datasets with blastp_path configured)
    - `num_threads`: int (threads per blastp process, default 1, max 4)
  - Response: `{ "job_id": "...", "status": "queued" }`

Example:

```bash
curl -X POST 'http://localhost:8080/api/blastp/jobs' \
  -H 'Content-Type: application/json' \
  -d '{"sequence":">query\nMSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGK","max_hits":10}'
```

#### Query Status

- `GET /api/blastp/jobs/{job_id}`
  - Response: `{ job_id, status, created_at, started_at, finished_at, progress:{phase,done,total}, error }`
  - `status`: `queued|running|completed|failed|cancelled`

#### Retrieve Results

- `GET /api/blastp/jobs/{job_id}/result?offset=&limit=&order_by=`
  - Only available when `status=completed`
  - Query parameters:
    - `offset`: Starting position (default 0)
    - `limit`: Number of results (default 100, max 500)
    - `order_by`: Sort order, options: `bitscore` (default), `evalue`, `pident`, `length`
  - Response:
    ```json
    {
      "params": {...},
      "total": 25,
      "offset": 0, "limit": 100,
      "order_by": "bitscore",
      "items": [
        {"dataset_id":"blood", "gene_id":"Afi_001234", "pident":95.5, "length":250, "evalue":1e-50, "bitscore":450},
        ...
      ]
    }
    ```

#### Cancel Job

- `DELETE /api/blastp/jobs/{job_id}`

#### Shared Database Handling

When multiple datasets are configured with the same `blastp_path`:
- Backend runs only one `blastp` command per unique `blastp_path` (deduplicated execution)
- Results are expanded to include a row for each dataset using that database (same gene appears multiple times, once per dataset)

### Tiles (PNG)

- `GET /d/{dataset}/tiles/{z}/{x}/{y}.png`: Base tile
- `GET /d/{dataset}/tiles/{z}/{x}/{y}/expression/{gene}.png?colormap=viridis`: Gene expression coloring (`colormap` is optional, defaults to viridis)
- `GET /d/{dataset}/tiles/{z}/{x}/{y}/category/{column}.png?categories=...`: Category coloring (optional filtering)
  - Filter (GET): `categories` can be repeated (`?categories=T&categories=B`), or JSON array (`?categories=["T","B"]`), or comma-separated (`?categories=T,B`)
  - Filter (POST): `POST /d/{dataset}/tiles/{z}/{x}/{y}/category/{column}.png`, body supports `["T","B"]` or `{"categories":["T","B"]}` (also supports form-encoded)
