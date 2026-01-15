# AtlasMap 使用指南（how-to-use）

AtlasMap 的典型工作流分两步：

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
      bins.X_tsne.zarr/     # （可选）其它坐标系的 Zarr store（当使用 --coord-key 多次指定时）
      metadata.json
      gene_index.json
    soma/                   # TileDBSOMA 存储（可选，默认启用）
      experiment.soma/      # 完整单细胞数据
```

**重要**：后端读取 `data.zarr_path` 指向的 `bins.zarr` 目录，并且会在其**上一级目录**查找 `metadata.json` 与 `gene_index.json`（即 `.../zarr/metadata.json`）。

---

## 方式 A：本地运行（开发/调试）

### 0) 依赖

- Go 1.22+
- Node.js 20+
- Python 3.9+

可选：用 `make install` 一次性安装依赖（会分别安装 Python/Go/前端依赖）。

### 0.5) 查看 H5AD 信息（可选）

在运行预处理前，可以先查看 `.h5ad` 的基本信息（细胞/基因数、layers、坐标键等）：

```bash
python scripts/h5ad_info.py -i data.h5ad
```

输出示例：

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

如需 JSON 格式输出，可加 `--json`。

### 1) 预处理数据（Python）

最简单（使用 Makefile）：

```bash
make preprocess INPUT=path/to/input.h5ad
```

或直接使用 CLI（更灵活）：

```bash
cd preprocessing
pip install -e .
atlasmap-preprocess run -i path/to/input.h5ad -o ../data/preprocessed -z 10 -g 500
```

常用参数（与代码一致，见 `preprocessing/atlasmap_preprocess/cli.py`）：

- `--coord-key`：`.h5ad` 里 `adata.obsm[...]` 的坐标键（可多次指定，生成多套坐标系 tiles）
- `--umap-key`：默认坐标键（旧参数名；用于指定默认/主坐标系，默认 `X_umap`）
- `--category/-c`：要写入的分类字段（来自 `adata.obs`，可多次指定）
- `--zoom-levels/-z`：缩放层级数（默认 8，从0开始, 一直到设置的参数为止)
- `--n-genes/-g`：预聚合基因数（默认 500）
- `--all-expressed/-a`：使用所有表达的基因而非 top N 基因
- `--min-cells`：基因至少在多少个细胞中表达（默认 3，仅与 `--all-expressed` 一起使用）
- `--no-soma`：禁用 TileDBSOMA 存储（默认启用，存储完整表达矩阵供任意基因查询）

如果需要指定 marker genes/更多高级参数：先生成配置文件再运行：

```bash
cd preprocessing
atlasmap-preprocess init-config -o config.yaml
# 编辑 config.yaml（如 marker_genes / category_columns 等）
atlasmap-preprocess from-config -c config.yaml
```

### 1.5) 验证预处理结果（可选）

使用 `visualize` 命令生成静态图片，验证预处理输出是否正确：

> 提示：这里的 `zoom` 表示"分箱分辨率"（越大越细）。以默认 `--zoom-levels 8` 为例：
> - `zoom=0`：1×1 个 bin，画出来通常只有 1 个点/块（非常粗）
> - `zoom=7`：128×128 个 bin，能看到完整的 UMAP 结构（建议作为默认检查/展示层）

```bash
# 基本用法（默认会画 3,5,7 三个 zoom，便于对比；随机选择 3 个基因）
atlasmap-preprocess visualize -i ../data/preprocessed/zarr -o ../data/figures

# 推荐：只看最高分辨率（最接近前端打开时的全局效果）
atlasmap-preprocess visualize -i ../data/preprocessed/zarr -o ../data/figures -z 7

# 指定 zoom 级别和基因
atlasmap-preprocess visualize -i ../data/preprocessed/zarr -o ../data/figures -z 3,5,7 -g CD3D -g CD8A
```

输出结构：
```
figures/
├── category/          # 按类别着色
│   ├── cell_type_zoom_7.png
│   └── cell_type_multi_zoom.png
└── expression/        # 按基因表达着色
    ├── GENE1_zoom_7.png
    └── GENE1_multi_zoom.png
```

### 1.6) 构建 BLASTP 数据库（可选，用于多物种检索）

如果需要使用 BLASTP 进行跨数据集的同源基因检索，需要为每个数据集构建蛋白序列数据库。

#### 前置要求

- 安装 NCBI BLAST+ 工具包（包含 `makeblastdb` 和 `blastp` 命令）

```bash
# Ubuntu/Debian
sudo apt-get install ncbi-blast+

# macOS
brew install blast

# 或使用 conda
conda install -c bioconda blast
```

#### 准备蛋白序列文件

蛋白序列文件需要是 FASTA 格式，

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

**重要提示**：
- 序列 ID **必须**与数据集中的 `gene_id` 一致
- 序列 ID（`>` 后面的部分）会作为 BLASTP 结果的 `sseqid` 字段
- 这个 ID 必须能在数据集的 `gene_index.json` 中找到，否则前端无法跳转到表达视图
- 建议使用数据集中的基因 ID 作为序列 ID，而非基因名称或其他标识符

#### 构建数据库

```bash
# 基本用法
makeblastdb -in proteins.fasta -dbtype prot -out /path/to/blast/db/proteins

# 示例：为 human 数据集构建数据库
mkdir -p data/preprocessed/blast
makeblastdb -in data/raw/human_proteins.fasta \
            -dbtype prot \
            -out data/preprocessed/blast/human_proteins
```

这会生成以下文件：
```
data/preprocessed/blast/
├── human_proteins.phr    # 头信息
├── human_proteins.pin    # 索引
├── human_proteins.psq    # 序列数据
└── human_proteins.pal    # (可选) 别名
```

#### 在配置文件中指定数据库

编辑 `config/server.yaml`，为每个数据集添加 `blastp_path`（**注意：路径前缀不含扩展名**）：

```yaml
data:
  blood:
    zarr_path: "/data/preprocessed/zarr/bins.zarr"
    soma_path: "/data/preprocessed/soma"
    blastp_path: "/data/preprocessed/blast/blood_proteins"  # 不含 .phr/.pin/.psq
  
  liver:
    zarr_path: "/data/liver/zarr/bins.zarr"
    soma_path: "/data/liver/soma"
    blastp_path: "/data/preprocessed/blast/liver_proteins"
```

#### 多数据集共享数据库

如果多个数据集的基因来自同一物种或同一基因集合，可以共享同一个 BLASTP 数据库：

```yaml
data:
  tissue_a:
    zarr_path: "/data/tissue_a/zarr/bins.zarr"
    blastp_path: "/data/shared/blast/species_proteins"
  
  tissue_b:
    zarr_path: "/data/tissue_b/zarr/bins.zarr"
    blastp_path: "/data/shared/blast/species_proteins"  # 共享同一数据库
```

后端会自动去重执行（只运行一次 `blastp`），但结果会展开到所有使用该数据库的数据集。

#### 验证数据库

```bash
# 检查数据库是否构建成功
blastdbcmd -db /data/preprocessed/blast/blood_proteins -info

# 测试查询（使用示例序列）
echo ">test
MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGK" | \
  blastp -db /data/preprocessed/blast/blood_proteins \
         -outfmt "6 sseqid pident length evalue bitscore" \
         -evalue 1e-5 \
         -max_target_seqs 5
```

配置完成后，前端工具栏会出现 🔍 按钮，点击即可进入 BLASTP 搜索页面。

### 2) 启动后端（Go）

后端 Go 模块位于 `server/`，推荐在 `server/` 目录启动，并显式指定配置文件：

```bash
cd server
go run ./cmd/server -config ../config/server.yaml
```

如果你需要启用 **TileDB-SOMA（任意基因 × 任意细胞表达查询）**：

第一步, 安装 TileDB C 库（`libtiledb` / headers）

```bash
conda install -c conda-forge -y tiledb
```

> 版本兼容：TileDB-Go 与 TileDB core 需要匹配。本仓库 `server/go.mod` 固定使用 `github.com/TileDB-Inc/TileDB-Go v0.38.0`，对应 TileDB `2.29.x`（例如 `2.29.1`）。

第二步, 使用 build tag 构建/运行（需要 CGO）：

```bash
cd server

export CGO_ENABLED=1
export CGO_CFLAGS="-I$CONDA_PREFIX/include"
export CGO_LDFLAGS="-L$CONDA_PREFIX/lib -ltiledb -Wl,-rpath,$CONDA_PREFIX/lib"

go run -tags soma ./cmd/server -config ../config/server.yaml
```

如果你看到一堆 `deprecated` 的编译 warning：不影响运行。想静默可加：

```bash
export CGO_CFLAGS="-I$CONDA_PREFIX/include -Wno-deprecated-declarations"
```

关键配置项在 `config/server.yaml`：

- `server.port`：默认 `8080`
- `data.zarr_path`：指向 `.../zarr/bins.zarr`（本地建议用相对路径或绝对路径）

**多数据集配置**：可以在 `data` 下定义多个数据集，第一个为默认数据集：

```yaml
data:
  pbmc:
    zarr_path: "/data/pbmc/zarr/bins.zarr"
    soma_path: "/data/pbmc/soma"
  liver:
    zarr_path: "/data/liver/zarr/bins.zarr"
    soma_path: "/data/liver/soma"
```

前端会显示数据集选择下拉框，API 可通过 `/d/{dataset}/api/...` 访问特定数据集。

如果预处理时通过 `--coord-key` 生成了多套坐标系，前端会额外显示 `Coord` 下拉框；也可以在任意 tiles/API 请求后追加 `?coord=X_tsne` 来切换坐标系。

启动后自检（注意：需先通过 `/api/datasets` 获取数据集 ID）：

```bash
curl http://localhost:8080/health
curl http://localhost:8080/api/datasets
curl http://localhost:8080/d/{dataset}/api/metadata  # 将 {dataset} 替换为实际 ID
```

### 3) 启动前端（Vite）

```bash
cd frontend
npm ci
npm run dev
```

默认打开：`http://localhost:3000`  
本地开发时，Vite 会把 `/api` 与 `/d` 代理到 `http://localhost:8080`（见 `frontend/vite.config.ts`）。

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

## 常见问题（Troubleshooting）

### 后端启动报错：`Failed to initialize Zarr reader`

通常是 `config/server.yaml` 里的 `data.zarr_path` 不对或预处理产物不完整：

- `data.zarr_path` 必须指向 `.../zarr/bins.zarr`（目录）
- `.../zarr/metadata.json` 与 `.../zarr/gene_index.json` 必须存在

### 预处理报错：`Coordinate key 'X_umap' not found`

输入 `.h5ad` 不包含对应的 `adata.obsm[...]` 坐标键，请用 `--coord-key`（或旧参数 `--umap-key`）指定正确键名。

### 基因表达 tile 一直空白/找不到 gene

表达 tile 只能对**预聚合基因**渲染；请先从 `GET /d/{dataset}/api/genes` 查看可用基因，或在预处理时提高 `--n-genes`，必要时使用 `from-config` 加入 marker genes。
