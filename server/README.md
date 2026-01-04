# SOMA-Tiles Server（Go）API 文档

本目录是 SOMA-Tiles 的 Go 后端：读取预处理产物（Zarr v3 store）并按需渲染 PNG tiles，同时提供数据集元信息、基因与分类相关的 JSON API。

## 快速启动

在仓库根目录执行（推荐）：

```bash
cd server
go run ./cmd/server -config ../config/server.yaml
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

多数据集时，前端会显示数据集选择下拉框，用户可以切换不同的数据集。

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

默认数据集（旧路由，向后兼容）：
- `GET /tiles/{z}/{x}/{y}.png`
- `GET /tiles/{z}/{x}/{y}/category/{column}.png`
- `GET /tiles/{z}/{x}/{y}/expression/{gene}.png?colormap={name}`

按数据集访问（多数据集路由）：
- `GET /d/{dataset}/tiles/{z}/{x}/{y}.png`
- `GET /d/{dataset}/tiles/{z}/{x}/{y}/category/{column}.png`
- `GET /d/{dataset}/tiles/{z}/{x}/{y}/expression/{gene}.png?colormap={name}`

### API（JSON）

默认数据集（旧路由，向后兼容）：
- `GET /api/metadata`
- `GET /api/stats`
- `GET /api/genes`
- `GET /api/genes/{gene}`
- `GET /api/genes/{gene}/stats`
- `GET /api/genes/{gene}/bins?threshold=&offset=&limit=`
- `GET /api/categories`
- `GET /api/categories/{column}/colors`
- `GET /api/categories/{column}/legend`

按数据集访问（多数据集路由）：
- `GET /d/{dataset}/api/metadata`
- `GET /d/{dataset}/api/stats`
- `GET /d/{dataset}/api/genes`
- 等等...

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

### `GET /tiles/{z}/{x}/{y}.png`

基础 tile（当前按 bin 的 cell_count 强度渲染）。

```bash
curl -sS -o tile.png "http://localhost:8080/tiles/0/0/0.png"
file tile.png
```

### `GET /tiles/{z}/{x}/{y}/category/{column}.png`

按分类字段着色（例如 `cell_type`）。`column` 来自 `/api/categories` 的 key。

```bash
curl -sS -o cell_type.png "http://localhost:8080/tiles/0/0/0/category/cell_type.png"
```

### `GET /tiles/{z}/{x}/{y}/expression/{gene}.png?colormap={name}`

按基因表达着色。

- `gene`：基因名（通常需要在预聚合基因列表内，见 `/api/genes`）
- `colormap`：可选，默认 `viridis`。支持：`viridis`、`plasma`、`inferno`、`magma`

```bash
curl -sS -o expr.png "http://localhost:8080/tiles/0/0/0/expression/CD3D.png?colormap=viridis"
```

## 元数据与统计（JSON）

### `GET /api/metadata`

返回 Zarr 元数据（由 `metadata.json` 加载，字段名以 Go struct 的 json tag 为准）。

```bash
curl -sS "http://localhost:8080/api/metadata"
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

### `GET /api/stats`

返回轻量统计信息：

```bash
curl -sS "http://localhost:8080/api/stats"
```

字段说明：

- `n_genes`：**预聚合基因数**（与 `/api/metadata` 的 `n_genes_preaggregated` 一致）

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

### `GET /api/genes`

返回预聚合基因列表（用于表达渲染/查询）。

```bash
curl -sS "http://localhost:8080/api/genes"
```

代表性响应（示例）：

```json
{
  "genes": ["CD3D", "CD8A", "MS4A1"],
  "total": 3,
  "preaggregated_count": 3
}
```

### `GET /api/genes/{gene}`

查询单个基因的索引信息（仅对预聚合基因生效；不存在则 404）。

```bash
curl -sS "http://localhost:8080/api/genes/CD3D"
```

代表性响应（示例）：

```json
{
  "name": "CD3D",
  "index": 0,
  "preaggregated": true
}
```

### `GET /api/genes/{gene}/stats`

基因表达统计（使用 zoom=0 的表达数组计算）。

```bash
curl -sS "http://localhost:8080/api/genes/CD3D/stats"
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

### `GET /api/genes/{gene}/bins?threshold=&offset=&limit=`

查询表达量大于等于阈值的 bins（使用 zoom=0 的表达数组，返回分页结果）。

Query 参数：

- `threshold`：浮点数，默认 `0`
- `offset`：整数，默认 `0`
- `limit`：整数，默认 `100`，最大 `1000`

状态码说明：

- `200`：成功
- `500`：后端查询失败（**当前实现里 gene 不存在也会走 500**，错误文本类似 `gene not found: XXX`）

```bash
curl -sS "http://localhost:8080/api/genes/CD3D/bins?threshold=1.0&offset=0&limit=2"
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

### `GET /api/categories`

返回所有可用分类字段及其取值（来自 metadata）。

```bash
curl -sS "http://localhost:8080/api/categories"
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

### `GET /api/categories/{column}/colors`

返回该分类字段的颜色映射（value -> hex color）。

```bash
curl -sS "http://localhost:8080/api/categories/cell_type/colors"
```

代表性响应（示例）：

```json
{
  "T": "#1f77b4",
  "B": "#ff7f0e",
  "Myeloid": "#2ca02c"
}
```

### `GET /api/categories/{column}/legend`

返回 legend 数组（用于前端渲染图例）。

```bash
curl -sS "http://localhost:8080/api/categories/cell_type/legend"
```

代表性响应（示例）：

```json
[
  { "value": "T", "color": "#1f77b4", "index": 0 },
  { "value": "B", "color": "#ff7f0e", "index": 1 },
  { "value": "Myeloid", "color": "#2ca02c", "index": 2 }
]
```

