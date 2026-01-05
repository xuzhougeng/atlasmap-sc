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

### SOMA（实验级表达查询，TileDB-SOMA）

> 注意：该接口需要后端以 `-tags soma` 构建，并且运行环境安装 TileDB C 库（`libtiledb` + headers）；否则会返回 `501`（not implemented）。

- `GET /d/{dataset}/api/soma/expression?gene={gene}&cells={c1,c2,...}&mode={sparse|dense}`
  - `gene`：基因名（`ms/RNA/var.gene_id`；注意不是 `gene_name`）
  - `cells`：细胞 `soma_joinid`（整数，逗号分隔；不是 `cell_id` 字符串）
  - `mode`：
    - `sparse`（默认）：返回非零项 `items=[{cell_joinid,value}, ...]`
    - `dense`：返回与输入 `cells` 对齐的 `values=[...]`（缺失视为 0）
  - 限制：`cells` 最多 10000 个（超出会 400）

示例：

```bash
curl -sS 'http://localhost:8080/d/retina/api/soma/expression?gene=ENSMICG00000038116&cells=95,118,143&mode=dense'
```

状态码说明（当前实现）：

- 404：该 dataset 未配置 `soma_path` 或 experiment 不存在
- 501：服务未以 `-tags soma` 构建
- 400：参数错误 / gene 不存在 / TileDB 查询失败（错误体为纯文本）

### SOMA 差异表达分析（DE Job，TileDB-SOMA）

> 同样需要后端以 `-tags soma` 构建。

差异表达分析为异步任务：提交后返回 `job_id`，前端轮询状态，完成后拉取结果。

#### 提交任务

- `POST /d/{dataset}/api/soma/de/jobs`
  - Body（JSON）：
    - `groupby`: string（obs 列名，如 `cell_type`）
    - `group1`: string[]（组1 类别，至少一个）
    - `group2`: string[]（组2 类别；空表示 one-vs-rest）
    - `tests`: string[]（`["ttest","ranksum"]`；默认两者都算）
    - `max_cells_per_group`: int（每组最多采样细胞数；默认 2000，上限 20000）
    - `seed`: int（采样随机种子；默认 0）
    - `limit`: int（返回 top N 基因；默认 50，上限 500）
  - Response: `{ "job_id": "...", "status": "queued" }`

示例：

```bash
curl -X POST 'http://localhost:8080/d/retina/api/soma/de/jobs' \
  -H 'Content-Type: application/json' \
  -d '{"groupby":"cell_type","group1":["T"],"group2":["B"]}'
```

#### 查询状态

- `GET /d/{dataset}/api/soma/de/jobs/{job_id}`
  - Response: `{ job_id, status, created_at, started_at, finished_at, progress:{phase,done,total}, error }`
  - `status`: `queued|running|completed|failed|cancelled`

#### 拉取结果

- `GET /d/{dataset}/api/soma/de/jobs/{job_id}/result?offset=&limit=&order_by=`
  - 仅当 `status=completed` 时可用
  - Query 参数：
    - `offset`：起始位置（默认 0）
    - `limit`：返回条数（默认 50，上限 500）
    - `order_by`：排序方式，可选值：`fdr_ranksum`（默认）、`fdr_ttest`、`p_ranksum`、`p_ttest`、`abs_log2fc`
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

#### 取消任务

- `DELETE /d/{dataset}/api/soma/de/jobs/{job_id}`

#### obs 元数据查询（用于前端下拉）

- `GET /d/{dataset}/api/soma/obs/columns`：返回 obs 可用列名
- `GET /d/{dataset}/api/soma/obs/{column}/values`：返回该列的唯一值

#### 统计方法说明

- **t-test**：Welch t-test（不假设等方差），p-value 用 Student-t 分布 CDF
- **ranksum**：Mann-Whitney U 检验（Wilcoxon rank-sum），带 tie 校正，正态近似 p-value
- **FDR**：Benjamini-Hochberg 多重校正

#### 持久化与限制

- **SQLite 持久化**：任务状态与结果存储在 SQLite 数据库（`de.sqlite_path`），服务重启后可继续查询历史结果
- **重启恢复**：服务重启时，正在运行的任务会被标记为 `failed`；排队中的任务会重新入队
- **TTL 清理**：完成的任务默认保留 **7 天**（`de.retention_days`），到期自动清理
- 每组最多采样 20000 细胞（`max_cells_per_group`）
- 默认同时运行最多 1 个 DE 任务（可通过配置 `de.max_concurrent` 调整），其余排队

### Tiles（PNG）

- `GET /d/{dataset}/tiles/{z}/{x}/{y}.png`：基础 tile
- `GET /d/{dataset}/tiles/{z}/{x}/{y}/expression/{gene}.png?colormap=viridis`：基因表达着色（`colormap` 可选，默认 viridis）
- `GET /d/{dataset}/tiles/{z}/{x}/{y}/category/{column}.png?categories=...`：分类着色（可选过滤）
  - 过滤（GET）：`categories` 可重复（`?categories=T&categories=B`），或 JSON 数组（`?categories=["T","B"]`），或逗号分隔（`?categories=T,B`）
  - 过滤（POST）：`POST /d/{dataset}/tiles/{z}/{x}/{y}/category/{column}.png`，body 支持 `["T","B"]` 或 `{"categories":["T","B"]}`（也兼容 form-encoded）
