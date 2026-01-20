# 前后端分离部署

鉴于高性能的云服务器价格高昂，我们推荐通过前后端分离的方式进行项目的部署。

部署的时候需要注意如下两个问题:

1. 当页面通过 HTTPS 加载时（https://），前端 JavaScript 不能再去请求不安全的 HTTP 资源（http://），包括 API 接口，浏览器会直接拦截请求，表现为 fetch/XHR 报错。
1. 跨域访问仍可能受到 CORS、证书有效性、Cookie Secure 属性等因素影响

## NIGNX 托管静态页面，转发API到HTTP后台 (推荐)

前端服务器，通过docker进行部署。

```bash
docker run -d --name hoptoper/atlasmap-frontend:latest \
  -p 3000:80 \
  -e BACKEND_HOST=<YOUR_BACKEND_SERVER_NAME> \
  -e BACKEND_PORT=<YOUR_PORT> \
  --restart unless-stopped \
  hoptoper/atlasmap-frontend:latest
```

后端处理好数据，配置server.yaml，然后启动服务

```bash
server/bin/server -config "data/server.yaml"
```

设置后台服务器nginx转发

```yaml
server {
    listen <YOUR_PORT>;
    server_name <YOUR_BACKEND_SERVER_NAME>;

    # 可选：如果 tile/结果较大
    client_max_body_size 50m;

    # /api (no trailing slash) -> /api/
    location = /api { return 301 /api/; }

    location /api/ {
        proxy_pass http://127.0.0.1:8080/api/;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 300s;
    }

    # /d (no trailing slash) -> /d/
    location = /d { return 301 /d/; }

    location /d/ {
        proxy_pass http://127.0.0.1:8080/d/;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 300s;
    }

    # 可选：健康检查直通
    location /health {
        proxy_pass http://127.0.0.1:8080/health;
    }
}
```

## Cloudflare Pages + Workers

在frontend生成用于发布的静态文件

```bash
UPSTREAM_ORIGIN=http://your.backend.server npm run build:cloudflare
```

此时会在frontend目录下生成两个子目录:

- dist: 静态资源站点
- functions : 基于 worker 实现反向代理

使用 wrangler 一键发布

```bash
# 1) 登录
npx wrangler login

# 2) 创建 Pages 项目
npx wrangler pages project create <PROJECT_NAME>

# 3) 部署：在项目根目录执行（根目录要有 functions/）
npx wrangler pages deploy dist --project-name <PROJECT_NAME>
```

注: 后台服务器 nginx 转发设置保持不变，同 [NIGNX 托管静态页面，转发API到HTTP后台 (推荐)](#nignx-托管静态页面转发api到http后台推荐).

可选: 自定义域名

## Cloudflare Pages + Cloudflare Tunnel

目标：**无需公网服务器/公网 IP**（后端在内网/NAT/校园网），实现**前后端分离部署**：

- 前端：Cloudflare Pages 托管静态站点
- 后端：任意一台可出网的机器上运行 `server/bin/server`（只需本机可监听端口）
- 连接：用 Cloudflare Tunnel 把后端暴露到一个 Cloudflare 子域名；前端通过 Pages Functions 同源代理 `/api/*`、`/d/*` 到该子域名（浏览器侧无 CORS 问题）

注: 后端主机必须能连通到Cloudflare tunnel服务, 参考[Connectivity pre-checks](https://developers.cloudflare.com/cloudflare-one/networks/connectors/cloudflare-tunnel/troubleshoot-tunnels/connectivity-prechecks/), 即，如果内网服务器被屏蔽这类服务，或者无法访问外网，无法使用改方案。

### 后端：启动服务（内网机器即可）

后端数据准备与 `server.yaml` 配置同其它方式一致，然后在后端机器上启动：

```bash
server/bin/server -config "data/server.yaml"
```

默认后端监听 `:8080`（以实际配置为准）。**不需要**给这台机器配置公网 IP，也不需要额外的 nginx 反代（Tunnel 会负责入口）。

### 后端：创建并运行 Cloudflare Tunnel

从Github上下载cloudflared

```bash
wget https://github.com/cloudflare/cloudflared/releases/download/2026.1.1/cloudflared-linux-amd64
mv cloudflared-linux-amd64 cloudflared
chmod +x cloudflared
```

参考[Create a locally-managed tunnel](https://developers.cloudflare.com/cloudflare-one/networks/connectors/cloudflare-tunnel/do-more-with-tunnels/local-management/create-local-tunnel/)

安装完成后，登录并创建 Tunnel, 此时会弹出一个网页，要求你选择一个 Cloudflare 托管域名下的子域名作为上游入口

```bash
cloudflared tunnel login
```

创建一个隧道，名字自取，这里我写的是atlasmap

```bash
cloudflared tunnel create atlasmap
```

后续步骤中的`<USER>`, `TUNNEL_ID` 和  `HOST_DOMAIN` 是需要根据实际结果进行修改。

创建 `cloudflared` 配置文件（示例路径 `~/.cloudflared/config.yml`）：

```yaml
tunnel: atlasmap
credentials-file: /home/<USER>/.cloudflared/<TUNNEL_ID>.json

ingress:
  - hostname: <HOST_DOMAIN>
    service: http://127.0.0.1:8080
  - service: http_status:404
```

建立DNS记录（非常重要，也可以在Cloudflare dashboard中DNS里面手动创建）

```bash
cloudflared tunnel route dns atlasmap <HOST_DOMAIN>
```

启动 Tunnel：

```bash
cloudflared tunnel run atlasmap
```

日志信息只要不出现ERR都算运行正常。

可选：将 Tunnel 配置为 systemd 服务，保证重启后自启动（按 cloudflared 官方文档即可）。

### 前端：构建 Cloudflare Pages 产物（包含 Functions 代理）

本项目已提供 Pages Functions 的同源代理实现：部署到 Pages 后，`/api/*` 与 `/d/*` 会被转发到 `UPSTREAM_ORIGIN` 指定的上游。

在 `frontend/` 目录构建（把上游指向刚才的 Tunnel 子域名）：

```bash
cd frontend
UPSTREAM_ORIGIN=<HOST_DOMAIN> npm run build:cloudflare
```

构建完成后会生成：

- `frontend/dist/`：静态站点
- `frontend/functions/`：Pages Functions（同源反代 `/api/*`、`/d/*`）

### 部署到 Cloudflare Pages

用 wrangler 一键发布

```bash
cd frontend
npx wrangler login
npx wrangler pages project create <PROJECT_NAME>
npx wrangler pages deploy dist --project-name <PROJECT_NAME>
```

### 验证与排错

- 访问 Pages 域名（例如 `https://<PROJECT>.pages.dev`）
- 打开浏览器 Network：确认请求打到同源的 `/api/...`、`/d/...`（由 Pages Functions 转发到 `https://api.example.com`）
- 如果接口返回 502/超时：
  - 确认后端机器上 `server` 正在运行，且本机可访问 `http://127.0.0.1:8080/health`（或你配置的健康检查）
  - 确认 `cloudflared tunnel run ...` 正在运行，且 `api.example.com` 已在 Cloudflare DNS 中生效
  - 确认 `UPSTREAM_ORIGIN` 使用 `https://api.example.com`（Tunnel 默认提供 HTTPS 入口）