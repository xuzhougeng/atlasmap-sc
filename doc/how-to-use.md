# SOMA-Tiles 使用指南（how-to-use）

SOMA-Tiles 的典型工作流分两步：

1. **离线预处理（Python）**：把 `.h5ad` 转成多分辨率的 `Zarr`（`bins.zarr`）+ 元数据（`metadata.json`/`gene_index.json`）。
2. **在线服务（Go + 前端）**：Go 服务器按需读取 `Zarr` 并渲染 PNG tiles；前端通过 `/api` 和 `/tiles` 展示交互式 UMAP。

---

## 目录与数据约定

仓库里的 `data/` 默认被 `.gitignore` 忽略，建议本地按以下结构组织：

```text
data/
  raw/                      # 原始输入（可选）
    input.h5ad
  preprocessed/             # 预处理输出
    zarr/
      bins.zarr/            # Zarr store（注意：是目录）
      metadata.json
      gene_index.json
```

**重要**：后端读取 `data.zarr_path` 指向的 `bins.zarr` 目录，并且会在其**上一级目录**查找 `metadata.json` 与 `gene_index.json`（即 `.../zarr/metadata.json`）。

---

## 方式 A：本地运行（开发/调试）

### 0) 依赖

- Go 1.22+
- Node.js 20+
- Python 3.9+

可选：用 `make install` 一次性安装依赖（会分别安装 Python/Go/前端依赖）。

### 1) 预处理数据（Python）

最简单（使用 Makefile）：

```bash
make preprocess INPUT=path/to/input.h5ad
```

或直接使用 CLI（更灵活）：

```bash
cd preprocessing
pip install -e .
soma-preprocess run -i path/to/input.h5ad -o ../data/preprocessed -z 8 -g 500
```

常用参数（与代码一致，见 `preprocessing/soma_tiles_preprocess/cli.py`）：

- `--umap-key`：`.h5ad` 里 `adata.obsm[...]` 的坐标键（默认 `X_umap`）
- `--category/-c`：要写入的分类字段（来自 `adata.obs`，可多次指定）
- `--zoom-levels/-z`：缩放层级数（默认 8，范围 1–12）
- `--n-genes/-g`：预聚合基因数（默认 500）
- `--all-expressed/-a`：使用所有表达的基因而非 top N 基因
- `--min-cells`：基因至少在多少个细胞中表达（默认 3，仅与 `--all-expressed` 一起使用）

如果需要指定 marker genes/更多高级参数：先生成配置文件再运行：

```bash
cd preprocessing
soma-preprocess init-config -o config.yaml
# 编辑 config.yaml（如 marker_genes / category_columns 等）
soma-preprocess from-config -c config.yaml
```

### 1.5) 验证预处理结果（可选）

使用 `visualize` 命令生成静态图片，验证预处理输出是否正确：

```bash
# 基本用法（随机选择 3 个基因）
soma-preprocess visualize -i ../data/preprocessed/zarr -o ../data/figures

# 指定 zoom 级别和基因
soma-preprocess visualize -i ../data/preprocessed/zarr -o ../data/figures -z 3,5,7 -g CD3D -g CD8A
```

输出结构：
```
figures/
├── category/          # 按类别着色
│   ├── cell_type_zoom_3.png
│   └── cell_type_multi_zoom.png
└── expression/        # 按基因表达着色
    ├── GENE1_zoom_3.png
    └── GENE1_multi_zoom.png
```

### 2) 启动后端（Go）

后端 Go 模块位于 `server/`，推荐在 `server/` 目录启动，并显式指定配置文件：

```bash
cd server
go run ./cmd/server -config ../config/server.yaml
```

关键配置项在 `config/server.yaml`：

- `server.port`：默认 `8080`
- `data.zarr_path`：指向 `.../zarr/bins.zarr`（本地建议用相对路径或绝对路径）

启动后自检：

```bash
curl http://localhost:8080/health
curl http://localhost:8080/api/metadata
```

### 3) 启动前端（Vite）

```bash
cd frontend
npm ci
npm run dev
```

默认打开：`http://localhost:3000`  
本地开发时，Vite 会把 `/api` 与 `/tiles` 代理到 `http://localhost:8080`（见 `frontend/vite.config.ts`）。

---

## 方式 B：Docker Compose（部署/演示）

### 1) 准备配置与数据路径

`docker-compose.yml` 会把宿主机 `./data` 挂载到容器内的 `/data`，因此容器内应使用：

- `data.zarr_path: "/data/preprocessed/zarr/bins.zarr"`

如果你直接使用 `config/server.yaml`，请确认其中的 `data.zarr_path` **对容器内路径有效**；否则后端容器会因找不到 Zarr 数据而启动失败。

### 2) 构建并启动

```bash
docker compose up -d --build
```

访问：

- 前端：`http://localhost:3000`
- 后端健康检查：`http://localhost:8080/health`

### 3) （可选）在 Docker 里跑预处理

`docker-compose.yml` 里 `preprocess` 服务默认在 `tools` profile 下，需要显式启用：

```bash
mkdir -p data/raw data/preprocessed
# 把 input.h5ad 放到 data/raw/input.h5ad

docker compose --profile tools run --rm preprocess run \
  -i /data/raw/input.h5ad \
  -o /data/preprocessed \
  -z 8 -g 500
```

预处理完成后，确保生成 `data/preprocessed/zarr/bins.zarr` 与同级的 `metadata.json`/`gene_index.json`，再启动/重启后端服务即可。

---

## 常用接口（后端）

代码位置：`server/internal/api/routes.go`

- `GET /health`：健康检查
- `GET /api/metadata`：数据集元数据（前端启动依赖）
- `GET /api/genes`：预聚合基因列表
- `GET /api/categories`：可用分类字段与取值
- `GET /tiles/{z}/{x}/{y}.png`：基础 tile
- `GET /tiles/{z}/{x}/{y}/expression/{gene}.png?colormap=viridis`：基因表达着色
- `GET /tiles/{z}/{x}/{y}/category/{column}.png`：分类着色

---

## 常见问题（Troubleshooting）

### 后端启动报错：`Failed to initialize Zarr reader`

通常是 `config/server.yaml` 里的 `data.zarr_path` 不对或预处理产物不完整：

- `data.zarr_path` 必须指向 `.../zarr/bins.zarr`（目录）
- `.../zarr/metadata.json` 与 `.../zarr/gene_index.json` 必须存在

### 预处理报错：`UMAP key 'X_umap' not found`

输入 `.h5ad` 不包含 `adata.obsm['X_umap']`，请用 `--umap-key` 指定正确键名。

### 基因表达 tile 一直空白/找不到 gene

表达 tile 只能对**预聚合基因**渲染；请先从 `GET /api/genes` 查看可用基因，或在预处理时提高 `--n-genes`，必要时使用 `from-config` 加入 marker genes。
