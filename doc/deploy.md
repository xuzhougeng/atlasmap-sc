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

```
server/bin/server -config "data/server.yaml"
```

设置后台服务器nginx转发

```
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
- functionns : 基于worker实现方向代理

使用cloudflare wrangler进行发布

```bash
# 1) 登录
npx wrangler login

# 2) （可选）创建 Pages 项目
npx wrangler pages project create <PROJECT_NAME>

# 3) 部署：在项目根目录执行（根目录要有 functions/）
npx wrangler pages deploy dist --project-name <PROJECT_NAME>
```

注: 后台服务器 nginx 转发设置保持不变，同 [NIGNX 托管静态页面，转发API到HTTP后台 (推荐)](#nignx-托管静态页面转发api到http后台推荐).

可选: 自定义域名