# 容器部署（Docker Compose）

本项目提供 `docker-compose.yml`，包含：

- `server`：Go 后端（`8080`）
- `frontend`：Nginx 静态站点 + 反向代理到后端（宿主机 `3000` -> 容器 `80`）
- `preprocess`：可选的预处理容器（默认在 `tools` profile 下，手动运行）

### 0) 前置条件

- 已安装 Docker 与 Docker Compose（`docker compose`）
- 已准备好预处理产物（`bins.zarr` + `metadata.json` + `gene_index.json`），或计划在容器内运行预处理

### 1) 准备数据目录

`docker-compose.yml` 会把宿主机的 `./data` 挂载到容器内的 `/data`：

- 后端：`./data:/data:ro`
- 预处理：`./data:/data`（可写）

推荐多数据集结构：

```text
data/
  dataset_a/
    zarr/metadata.json
    zarr/gene_index.json
    zarr/bins.zarr/...
    soma/...(可选)
  dataset_b/
    ...
```

### 2) 配置后端路径（容器内路径）

`server` 容器通过 `./config:/app/config:ro` 挂载配置文件，因此需要保证 `config/server.yaml` 里使用的是**容器内路径**（`/data/...`），例如：

```yaml
data:
  dataset_a:
    zarr_path: "/data/dataset_a/zarr/bins.zarr"
    soma_path: "/data/dataset_a/soma"
  dataset_b:
    zarr_path: "/data/dataset_b/zarr/bins.zarr"
    soma_path: "/data/dataset_b/soma"
```

### 3) 启动服务

在仓库根目录执行：

```bash
docker compose up -d --build
```

访问：

- 前端：`http://localhost:3000`
- 后端健康检查：`http://localhost:8080/health`

说明：生产环境通常只需要对外暴露 `frontend`（`3000` 或 `80/443`）；`frontend` 会把 `/api` 与 `/d/...` 反向代理到后端。

查看日志：

```bash
docker compose logs -f server
docker compose logs -f frontend
```

停止并清理：

```bash
docker compose down
```

### 4) （可选）在容器内运行预处理

`preprocess` 服务在 `tools` profile 下，示例：

```bash
mkdir -p data/raw data/dataset_a
# 把 input.h5ad 放到 data/raw/input.h5ad

docker compose --profile tools run --rm preprocess run \
  -i /data/raw/input.h5ad \
  -o /data/dataset_a \
  -z 8 -g 500
```

预处理完成后，更新 `config/server.yaml` 指向对应 `/data/...` 路径，并重启后端：

```bash
docker compose restart server
```

### 5) 对外提供服务（可选）

- 修改 `docker-compose.yml` 的端口映射，将 `frontend` 暴露到 `80/443`，或在外部反向代理（Nginx/Traefik）前置。
- 若浏览器直接访问后端（不经 `frontend` 代理），需要在 `config/server.yaml` 的 `server.cors_origins` 中加入实际前端域名来源。
