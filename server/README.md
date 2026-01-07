# SOMA-Tiles Server（Go）API 文档

本目录是 SOMA-Tiles 的 Go 后端：读取预处理产物（Zarr v3 store）并按需渲染 PNG tiles，同时提供数据集元信息、基因与分类相关的 JSON API。

## 快速启动

在仓库根目录执行（推荐）：

```bash
cd server
go run ./cmd/server -config ../config/server.yaml
```

如果你需要启用 **TileDB-SOMA（任意基因 × 任意细胞表达查询）**（需要 CGO + TileDB C 库）：

```bash
cd server

# 安装 TileDB core（conda-forge）
conda install -c conda-forge -y tiledb

export CGO_ENABLED=1
export CGO_CFLAGS="-I$CONDA_PREFIX/include"
export CGO_LDFLAGS="-L$CONDA_PREFIX/lib -ltiledb -Wl,-rpath,$CONDA_PREFIX/lib"

go run -tags soma ./cmd/server -config ../config/server.yaml
```

默认监听端口来自配置文件 `server.port`（默认 `8080`）。健康检查：

```bash
curl -sS http://localhost:8080/health
```

## 配置文件（`config/server.yaml`）

后端启动时会读取 YAML 配置（见 `../config/server.yaml`）。

### 单数据集配置（旧格式，仍然支持）

- `data.zarr_path`: 指向预处理输出的 `bins.zarr` **目录**
- 注意：服务会在 `bins.zarr` 的**上一级目录**读取 `metadata.json` 与 `gene_index.json`

```yaml
server:
  port: 8080
  cors_origins:
    - "http://localhost:3000"

data:
  zarr_path: "/abs/path/to/data/preprocessed/zarr/bins.zarr"
  soma_path: "/abs/path/to/data/preprocessed/soma"

cache:
  tile_size_mb: 512
  tile_ttl_minutes: 10

render:
  tile_size: 256
  default_colormap: viridis

de:
  max_concurrent: 1                      # 同时运行的最大 DE 任务数（默认 1）
  sqlite_path: "./data/de_jobs.sqlite"   # SQLite 数据库路径（持久化任务状态与结果）
  retention_days: 7                      # 已完成任务保留天数（默认 7 天）
```

### 多数据集配置（新格式）

可以在 `data` 下定义多个数据集，每个数据集有自己的 `zarr_path` 和 `soma_path`。**第一个数据集会作为默认数据集**。

```yaml
server:
  port: 8080
  cors_origins:
    - "http://localhost:3000"

data:
  pbmc:
    zarr_path: "/data/pbmc/zarr/bins.zarr"
    soma_path: "/data/pbmc/soma"
    blastp_path: "/data/pbmc/blast/proteins"  # 可选：BLASTP 数据库路径前缀
  liver:
    zarr_path: "/data/liver/zarr/bins.zarr"
    soma_path: "/data/liver/soma"
    blastp_path: "/data/liver/blast/proteins"

cache:
  tile_size_mb: 512
  tile_ttl_minutes: 10

render:
  tile_size: 256
  default_colormap: viridis

de:
  max_concurrent: 1
  sqlite_path: "./data/de_jobs.sqlite"
  retention_days: 7
```

多数据集时，前端会显示数据集选择下拉框，用户可以切换不同的数据集。

### BLASTP 数据库配置

`blastp_path` 指定 BLASTP 数据库的路径前缀（不含 `.phr`、`.pin`、`.psq` 等扩展名）。该路径直接用于 `blastp -db` 参数。

**要求**：
- 系统 PATH 中需有 `blastp` 可执行文件（NCBI BLAST+）
- 数据库需使用 `makeblastdb` 预先构建
- BLAST 输出的 `sseqid` 字段应与数据集中的 `gene_id` 一致

**示例：构建 BLASTP 数据库**

```bash
# 假设 proteins.fasta 中序列 ID 格式为 >Afi_001234
makeblastdb -in proteins.fasta -dbtype prot -out /data/pbmc/blast/proteins
```

**多数据集共享同一数据库**

如果多个数据集配置了相同的 `blastp_path`，后端会：
1. 对该数据库只执行一次 `blastp` 命令（去重执行，节省资源）
2. 将结果展开到所有关联的数据集（结果表中会出现多行，每个数据集一行）

## 多坐标系（coord，可选）

如果预处理时通过 `--coord-key` 生成了多套坐标系（`metadata.json` 里包含 `coordinate_systems`），后端会在启动时自动加载同一 `zarr/` 目录下的多套 `bins.*.zarr`。

所有数据集 scoped 的 Tiles/API 接口都支持可选的查询参数：

- `coord={key}`：坐标系 key（通常对应 `.h5ad` 的 `adata.obsm[key]`，例如 `X_umap`、`X_tsne`）

示例：

```bash
curl "http://localhost:8080/d/pbmc/api/metadata?coord=X_tsne"
curl "http://localhost:8080/d/pbmc/tiles/2/1/1/category/cell_type.png?coord=X_tsne" -o tile.png
```

不带 `coord` 时会使用默认坐标系（`default_coordinate_system`）。

## 通用约定

- Base URL：`http://localhost:8080`
- JSON 接口：`Content-Type: application/json`
- Tile 接口：`Content-Type: image/png`，并带有 `Cache-Control: public, max-age=3600`
- 错误响应（JSON 接口）：多数情况下返回纯文本错误体（`http.Error`），并使用：
  - 400：参数非法
  - 404：资源不存在（如 gene/category 不存在）
  - 500：内部错误
- Tile 接口容错：如果渲染/读取数据失败，通常会返回**空白（透明）PNG tile**，HTTP 状态仍为 `200`（用于前端地图不因局部错误中断）。
- 示例说明：本文档的 JSON 示例会用小规模数据展示字段结构；实际生产数据里数组可能很大。

## 接口总览

### Health

- `GET /health`

### 数据集列表（全局）

- `GET /api/datasets` — 返回可用数据集列表及默认数据集

### Tiles（PNG）

所有 Tiles 接口均需在路径前加 `/d/{dataset}`：
- `GET /d/{dataset}/tiles/{z}/{x}/{y}.png`
- `GET /d/{dataset}/tiles/{z}/{x}/{y}/category/{column}.png`
- `GET /d/{dataset}/tiles/{z}/{x}/{y}/expression/{gene}.png?colormap={name}`

### API（JSON）

所有 API 接口均需在路径前加 `/d/{dataset}`：
- `GET /d/{dataset}/api/metadata`
- `GET /d/{dataset}/api/stats`
- `GET /d/{dataset}/api/genes`
- `GET /d/{dataset}/api/genes/{gene}`
- `GET /d/{dataset}/api/genes/{gene}/stats`
- `GET /d/{dataset}/api/genes/{gene}/bins?threshold=&offset=&limit=`
- `GET /d/{dataset}/api/soma/expression?gene={gene}&cells={c1,c2,...}&mode={sparse|dense}`（需要 `-tags soma`）
- `GET /d/{dataset}/api/categories`
- `GET /d/{dataset}/api/categories/{column}/colors`
- `GET /d/{dataset}/api/categories/{column}/legend`

## Health

### `GET /health`

用于探活。

```bash
curl -sS http://localhost:8080/health
```

响应（纯文本）：

```text
OK
```

## Tiles（PNG）

> `{z}/{x}/{y}` 使用 Leaflet 标准 tile 坐标：`z` 为缩放级别，`x/y` 为 tile 索引。
> 所有 Tiles 接口均需在路径前加 `/d/{dataset}`。

### `GET /d/{dataset}/tiles/{z}/{x}/{y}.png`

基础 tile（当前按 bin 的 cell_count 强度渲染）。

```bash
curl -sS -o tile.png "http://localhost:8080/d/default/tiles/0/0/0.png"
file tile.png
```

### `GET /d/{dataset}/tiles/{z}/{x}/{y}/category/{column}.png`

按分类字段着色（例如 `cell_type`）。`column` 来自 `/d/{dataset}/api/categories` 的 key。

```bash
curl -sS -o cell_type.png "http://localhost:8080/d/default/tiles/0/0/0/category/cell_type.png"
```

### `GET /d/{dataset}/tiles/{z}/{x}/{y}/expression/{gene}.png?colormap={name}`

按基因表达着色。

- `gene`：基因名（通常需要在预聚合基因列表内，见 `/d/{dataset}/api/genes`）
- `colormap`：可选，默认 `viridis`。支持：`viridis`、`plasma`、`inferno`、`magma`

```bash
curl -sS -o expr.png "http://localhost:8080/d/default/tiles/0/0/0/expression/CD3D.png?colormap=viridis"
```

## 元数据与统计（JSON）

> 所有 API 接口均需在路径前加 `/d/{dataset}`。

### `GET /d/{dataset}/api/metadata`

返回 Zarr 元数据（由 `metadata.json` 加载，字段名以 Go struct 的 json tag 为准）。

```bash
curl -sS "http://localhost:8080/d/default/api/metadata"
```

代表性响应（示例）：

```json
{
  "zoom_levels": 8,
  "tile_size": 256,
  "n_genes_preaggregated": 3,
  "n_categories": {
    "cell_type": 3
  },
  "format_version": "1",
  "dataset_name": "example_dataset",
  "n_cells": 100000,
  "preaggregated_genes": ["CD3D", "CD8A", "MS4A1"],
  "gene_index": {
    "CD3D": 0,
    "CD8A": 1,
    "MS4A1": 2
  },
  "categories": {
    "cell_type": {
      "values": ["T", "B", "Myeloid"],
      "mapping": {
        "T": 0,
        "B": 1,
        "Myeloid": 2
      }
    }
  },
  "bounds": {
    "min_x": 0,
    "max_x": 256,
    "min_y": 0,
    "max_y": 256
  }
}
```

### `GET /d/{dataset}/api/stats`

返回轻量统计信息：

```bash
curl -sS "http://localhost:8080/d/default/api/stats"
```

字段说明：

- `n_genes`：**预聚合基因数**（与 `/d/{dataset}/api/metadata` 的 `n_genes_preaggregated` 一致）

代表性响应（示例）：

```json
{
  "n_cells": 100000,
  "n_genes": 3,
  "zoom_levels": 8,
  "dataset_name": "example_dataset"
}
```

## 基因相关（JSON）

### `GET /d/{dataset}/api/genes`

返回预聚合基因列表（用于表达渲染/查询）。

```bash
curl -sS "http://localhost:8080/d/default/api/genes"
```

代表性响应（示例）：

```json
{
  "genes": ["CD3D", "CD8A", "MS4A1"],
  "total": 3,
  "preaggregated_count": 3
}
```

### `GET /d/{dataset}/api/genes/{gene}`

查询单个基因的索引信息（仅对预聚合基因生效；不存在则 404）。

```bash
curl -sS "http://localhost:8080/d/default/api/genes/CD3D"
```

代表性响应（示例）：

```json
{
  "name": "CD3D",
  "index": 0,
  "preaggregated": true
}
```

### `GET /d/{dataset}/api/genes/{gene}/stats`

基因表达统计（使用 zoom=0 的表达数组计算）。

```bash
curl -sS "http://localhost:8080/d/default/api/genes/CD3D/stats"
```

代表性响应（示例）：

```json
{
  "gene": "CD3D",
  "index": 0,
  "expressing_bins": 1234,
  "total_bins": 65536,
  "total_cells": 100000,
  "mean_expression": 0.42,
  "max_expression": 8.9
}
```

### `GET /d/{dataset}/api/genes/{gene}/bins?threshold=&offset=&limit=`

查询表达量大于等于阈值的 bins（使用 zoom=0 的表达数组，返回分页结果）。

Query 参数：

- `threshold`：浮点数，默认 `0`
- `offset`：整数，默认 `0`
- `limit`：整数，默认 `100`，最大 `1000`

状态码说明：

- `200`：成功
- `500`：后端查询失败（**当前实现里 gene 不存在也会走 500**，错误文本类似 `gene not found: XXX`）

```bash
curl -sS "http://localhost:8080/d/default/api/genes/CD3D/bins?threshold=1.0&offset=0&limit=2"
```

代表性响应（示例）：

```json
{
  "bins": [
    {
      "bin_index": 10,
      "bin_x": 120,
      "bin_y": 80,
      "cell_count": 15,
      "expression": 1.23,
      "categories": {
        "cell_type": "T"
      }
    },
    {
      "bin_index": 42,
      "bin_x": 121,
      "bin_y": 81,
      "cell_count": 9,
      "expression": 1.05,
      "categories": {
        "cell_type": "T"
      }
    }
  ],
  "total_count": 987,
  "offset": 0,
  "limit": 2,
  "gene": "CD3D",
  "threshold": 1.0
}
```

## 分类相关（JSON）

### `GET /d/{dataset}/api/categories`

返回所有可用分类字段及其取值（来自 metadata）。

```bash
curl -sS "http://localhost:8080/d/default/api/categories"
```

代表性响应（示例）：

```json
{
  "cell_type": {
    "values": ["T", "B", "Myeloid"],
    "mapping": {
      "T": 0,
      "B": 1,
      "Myeloid": 2
    }
  }
}
```

### `GET /d/{dataset}/api/categories/{column}/colors`

返回该分类字段的颜色映射（value -> hex color）。

```bash
curl -sS "http://localhost:8080/d/default/api/categories/cell_type/colors"
```

代表性响应（示例）：

```json
{
  "T": "#1f77b4",
  "B": "#ff7f0e",
  "Myeloid": "#2ca02c"
}
```

### `GET /d/{dataset}/api/categories/{column}/legend`

返回 legend 数组（用于前端渲染图例）。

```bash
curl -sS "http://localhost:8080/d/default/api/categories/cell_type/legend"
```

代表性响应（示例）：

```json
[
  { "value": "T", "color": "#1f77b4", "index": 0 },
  { "value": "B", "color": "#ff7f0e", "index": 1 },
  { "value": "Myeloid", "color": "#2ca02c", "index": 2 }
]
```

## SOMA 表达查询（TileDB-SOMA，实验级）

该接口用于从 `soma/experiment.soma` 的完整表达矩阵里查询 **任意基因 × 任意细胞** 的表达值。

前提：

- 配置里该数据集设置了 `soma_path`（指向 `.../soma` 目录；内部会寻找 `experiment.soma`）
- 后端以 `-tags soma` 构建/运行，并且系统能找到 TileDB core（`libtiledb` + headers）
- 版本兼容：本仓库 `server/go.mod` 固定 `TileDB-Go v0.38.0`，对应 TileDB core `2.29.x`

### `GET /d/{dataset}/api/soma/expression?gene=&cells=&mode=`

- `gene`：基因名（`ms/RNA/var.gene_id`）
- `cells`：细胞 `soma_joinid`（整数，逗号分隔；不是 `cell_id`）
- `mode`：
  - `dense`：返回与输入 `cells` 对齐的 `values`（缺失视为 0）
  - `sparse`：仅返回非零项 `items=[{cell_joinid,value}, ...]`

示例（retina）：

```bash
curl -sS 'http://localhost:8080/d/retina/api/soma/expression?gene=ENSMICG00000038116&cells=95,118,143,206,276,367,373,385,386,497&mode=dense'
```

返回（示例）：

```json
{
  "cells": [95,118,143,206,276,367,373,385,386,497],
  "gene": "ENSMICG00000038116",
  "gene_joinid": 2315,
  "mode": "dense",
  "nnz": 10,
  "values": [2.0556595,0.8455493,2.0095484,1.3307723,2.1727157,2.3340511,2.3182132,2.1789465,2.3348932,2.1009095]
}
```

状态码说明（当前实现）：

- 404：该 dataset 未配置 `soma_path` 或 experiment 不存在
- 501：服务未以 `-tags soma` 构建
- 400：参数错误 / gene 不存在 / TileDB 查询失败（错误体为纯文本）

> Docker 注意：当前 `server/Dockerfile` 使用 `CGO_ENABLED=0` 的静态构建，默认不支持 `-tags soma`。如需在容器内启用 SOMA，需要自定义镜像并安装 TileDB core。

## SOMA 差异表达分析（DE Job）

该接口用于在线计算两组细胞之间的差异表达基因，支持 **全基因** t-test 和秩和检验。

前提与 SOMA 表达查询相同：需要 `-tags soma` 构建并安装 TileDB core。

### 工作流程

1. **提交任务**：`POST /d/{dataset}/api/soma/de/jobs` 返回 `job_id`
2. **轮询状态**：`GET /d/{dataset}/api/soma/de/jobs/{job_id}` 直到 `status=completed`
3. **拉取结果**：`GET /d/{dataset}/api/soma/de/jobs/{job_id}/result?offset=&limit=`

### `POST /d/{dataset}/api/soma/de/jobs`

提交差异表达分析任务。

请求体（JSON）：

```json
{
  "groupby": "cell_type",
  "group1": ["T"],
  "group2": ["B", "Myeloid"],
  "tests": ["ttest", "ranksum"],
  "max_cells_per_group": 2000,
  "seed": 0,
  "limit": 50
}
```

参数说明：

- `groupby`：obs 列名（用于分组，如 `cell_type`）
- `group1`：组1 类别（至少一个）
- `group2`：组2 类别（可选；空表示 one-vs-rest，即 group1 vs 其他所有细胞）
- `tests`：要运行的检验（默认 `["ttest", "ranksum"]`）
- `max_cells_per_group`：每组最多采样细胞数（默认 2000，上限 20000）
- `seed`：采样随机种子（默认 0，保证可复现）
- `limit`：返回 top N 基因（默认 50，上限 500）

响应：

```json
{
  "job_id": "a1b2c3d4e5f6g7h8",
  "status": "queued"
}
```

### `GET /d/{dataset}/api/soma/de/jobs/{job_id}`

查询任务状态。

响应：

```json
{
  "job_id": "a1b2c3d4e5f6g7h8",
  "status": "running",
  "created_at": "2025-01-05T10:00:00Z",
  "started_at": "2025-01-05T10:00:01Z",
  "finished_at": null,
  "progress": {
    "phase": "scanning_group1",
    "done": 5000,
    "total": 20000
  },
  "error": ""
}
```

状态值：`queued` | `running` | `completed` | `failed` | `cancelled`

### `GET /d/{dataset}/api/soma/de/jobs/{job_id}/result`

拉取完成任务的结果（分页，支持排序）。

Query 参数：

- `offset`：起始位置（默认 0）
- `limit`：返回条数（默认 50，上限 500）
- `order_by`：排序方式，可选值：
  - `fdr_ranksum`（默认）：按秩和检验 FDR 升序，然后 |log2fc| 降序
  - `fdr_ttest`：按 t-test FDR 升序
  - `p_ranksum`：按秩和检验 p-value 升序
  - `p_ttest`：按 t-test p-value 升序
  - `abs_log2fc`：按 |log2fc| 降序

响应：

```json
{
  "params": { ... },
  "n1": 500,
  "n2": 1500,
  "total": 20000,
  "offset": 0,
  "limit": 50,
  "order_by": "fdr_ranksum",
  "items": [
    {
      "gene": "CD3D",
      "gene_joinid": 123,
      "mean1": 2.5,
      "mean2": 0.3,
      "pct1": 0.85,
      "pct2": 0.12,
      "log2fc": 3.06,
      "p_ttest": 1.2e-10,
      "fdr_ttest": 2.4e-6,
      "p_ranksum": 5.6e-12,
      "fdr_ranksum": 1.1e-7
    }
  ]
}
```

### `DELETE /d/{dataset}/api/soma/de/jobs/{job_id}`

取消任务（best-effort）。

### obs 元数据查询

用于前端获取可用的分组列和取值：

- `GET /d/{dataset}/api/soma/obs/columns`：返回 obs 列名列表
- `GET /d/{dataset}/api/soma/obs/{column}/values`：返回该列的唯一值

### 统计方法

- **Welch t-test**：不假设等方差的双样本 t 检验，p-value 用 Student-t 分布 CDF
- **Mann-Whitney U**：秩和检验（Wilcoxon rank-sum），带 tie 校正，使用正态近似计算 p-value
- **BH-FDR**：Benjamini-Hochberg 多重检验校正

### 持久化与恢复

- **SQLite 持久化**：任务状态与全基因结果存储在 SQLite 数据库（`de.sqlite_path`），不占用内存，服务重启后可继续查询历史结果
- **重启恢复**：
  - 正在运行的任务会被标记为 `failed`（error: "server restarted"）
  - 排队中的任务会自动重新入队
- **TTL 清理**：完成的任务默认保留 **7 天**（`de.retention_days`），每小时检查并自动清理过期任务

### 限制与注意事项

- **采样限制**：为控制计算成本，每组最多采样 20000 细胞（超出会随机采样）
- **并发限制**：默认同时最多运行 1 个 DE 任务（可通过 `de.max_concurrent` 配置），其余排队等待
- **全基因计算**：会对所有基因计算统计量，大数据集可能需要数分钟
