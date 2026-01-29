# 快速开始

运行本项目有两种常见方式：

- **不使用 Docker**：推荐用于本地开发和调试
- **使用 Docker**：推荐用于部署/演示——环境一致，依赖最少

更详细的使用说明，请参阅 [how-to-use_zh.md](./how-to-use_zh.md)。  
项目部署和服务器启动，请参阅 [deploy_zh.md](./deploy_zh.md)。

## 不使用 Docker（本地开发）

### 0) 前置依赖

- Python >= 3.9（推荐使用 `uv` 管理）
- Node.js 20+
- Go（可选；`make` 会在缺失时尝试自动安装）

### 1) 安装依赖

```bash
make install
```

### 2) 预处理数据

```bash
mkdir -p data/raw data/dataset_a
wget -O data/raw/input.h5ad https://datasets.cellxgene.cziscience.com/12fbed11-bd9d-43c2-97c1-6a79a2bb1be2.h5ad
make preprocess INPUT=data/raw/input.h5ad OUTPUT=data/dataset_a
```

### 3) 启动开发模式

```bash
make dev SERVER_CONFIG=../config/server.yaml
```

默认地址：

- 前端：`http://localhost:3000`
- 后端：`http://localhost:8080`


## 使用 Docker（推荐）

### 0) 前置依赖

- Docker + `docker compose`

### 1) 克隆仓库

```bash
git clone https://github.com/xuzhougeng/atlasmap-sc.git
cd atlasmap-sc
```

### 2) 准备数据目录（宿主机）

`docker-compose.yml` 会将宿主机的 `./data` 挂载到容器内的 `/data`。

推荐的多数据集目录结构（宿主机路径）：

```text
data/
  raw/
    input.h5ad
  dataset_a/
    zarr/metadata.json
    zarr/gene_index.json
    zarr/bins.zarr/...
    soma/...（可选）
```

### 3) 运行预处理（二选一）

#### 方式 A：在宿主机上运行预处理（推荐，迭代更快）

```bash
mkdir -p data/raw data/dataset_a
wget -O data/raw/input.h5ad https://datasets.cellxgene.cziscience.com/12fbed11-bd9d-43c2-97c1-6a79a2bb1be2.h5ad

# 使用 uv 创建虚拟环境并安装预处理依赖
make preprocess-venv

# 运行预处理，输出到 data/dataset_a
make preprocess INPUT=data/raw/input.h5ad OUTPUT=data/dataset_a
```

#### 方式 B：在容器中运行预处理（环境最一致）

```bash
mkdir -p data/raw data/dataset_a
# 将 input.h5ad 放到 data/raw/input.h5ad

docker compose --profile tools run --rm preprocess run \
  -i /data/raw/input.h5ad \
  -o /data/dataset_a \
  -z 8 -g 500
```

### 4) 配置后端数据路径（容器路径）

由于容器内看到的是 `/data/...`，`config/server.yaml` 必须使用容器内部路径，例如：

```yaml
data:
  dataset_a:
    zarr_path: "/data/dataset_a/zarr/bins.zarr"
    soma_path: "/data/dataset_a/soma"  # 可选
```

### 5) 启动服务

```bash
docker compose up -d --build
```

访问地址：

- 前端：`http://localhost:3000`
- 后端健康检查：`http://localhost:8080/health`

查看日志：

```bash
docker compose logs -f server
docker compose logs -f frontend
```

停止服务：

```bash
docker compose down
```
