## 常用接口（后端）

代码位置：`server/internal/api/routes.go`

通用说明：
- 所有 API 与 Tiles 接口均需在路径前加 `/d/{dataset}`，例如 `/d/pbmc/api/metadata`、`/d/pbmc/tiles/0/0/0.png`
- `{dataset}` 为数据集 ID，可通过 `GET /api/datasets` 获取可用列表及默认值

### Health

- `GET /health`：健康检查

### 数据集（全局）

- `GET /api/datasets`：列出可用数据集与默认数据集（该接口不加 `/d/{dataset}`）

### API（JSON）

- `GET /d/{dataset}/api/metadata`：数据集元数据（前端启动依赖）
- `GET /d/{dataset}/api/stats`：统计信息（如 n_cells、n_genes、zoom_levels、dataset_name）
- `GET /d/{dataset}/api/genes`：预聚合基因列表
- `GET /d/{dataset}/api/genes/{gene}`：基因信息（如 index/preaggregated）
- `GET /d/{dataset}/api/genes/{gene}/stats?zoom={z}`：基因统计（`zoom` 可选，默认 0）
- `GET /d/{dataset}/api/genes/{gene}/bins?threshold=&offset=&limit=`：查询表达该基因的 bins（limit 默认 100，最大 1000）
- `GET /d/{dataset}/api/genes/{gene}/category/{column}/means`：该基因在指定分类上的均值（用于类别对比）
- `GET /d/{dataset}/api/categories`：可用分类字段与取值
- `GET /d/{dataset}/api/categories/{column}/colors`：分类颜色映射
- `GET /d/{dataset}/api/categories/{column}/legend`：分类图例

### Tiles（PNG）

- `GET /d/{dataset}/tiles/{z}/{x}/{y}.png`：基础 tile
- `GET /d/{dataset}/tiles/{z}/{x}/{y}/expression/{gene}.png?colormap=viridis`：基因表达着色（`colormap` 可选，默认 viridis）
- `GET /d/{dataset}/tiles/{z}/{x}/{y}/category/{column}.png?categories=...`：分类着色（可选过滤）
  - 过滤（GET）：`categories` 可重复（`?categories=T&categories=B`），或 JSON 数组（`?categories=["T","B"]`），或逗号分隔（`?categories=T,B`）
  - 过滤（POST）：`POST /d/{dataset}/tiles/{z}/{x}/{y}/category/{column}.png`，body 支持 `["T","B"]` 或 `{"categories":["T","B"]}`（也兼容 form-encoded）
