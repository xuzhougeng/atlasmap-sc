# Frontend-Backend Separation Deployment

Given the high cost of high-performance cloud servers, we recommend deploying the project with frontend-backend separation.

When deploying, pay attention to the following two issues:

1. When pages are loaded via HTTPS (https://), frontend JavaScript cannot request insecure HTTP resources (http://), including API endpoints. The browser will block requests directly, manifesting as fetch/XHR errors.
2. Cross-origin access may still be affected by CORS, certificate validity, Cookie Secure attributes, and other factors

## NGINX Hosting Static Pages, Forwarding API to HTTP Backend (Recommended)

Deploy the frontend server using Docker.

```bash
docker run -d --name hoptoper/atlasmap-frontend:latest \
  -p 3000:80 \
  -e BACKEND_HOST=<YOUR_BACKEND_SERVER_NAME> \
  -e BACKEND_PORT=<YOUR_PORT> \
  --restart unless-stopped \
  hoptoper/atlasmap-frontend:latest
```

Prepare the backend data, configure server.yaml, then start the service:

```bash
server/bin/server -config "data/server.yaml"
```

Configure nginx forwarding on the backend server:

```yaml
server {
    listen <YOUR_PORT>;
    server_name <YOUR_BACKEND_SERVER_NAME>;

    # Optional: if tiles/results are large
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

    # Optional: health check passthrough
    location /health {
        proxy_pass http://127.0.0.1:8080/health;
    }
}
```

## Cloudflare Pages + Workers

Generate static files for deployment in frontend:

```bash
UPSTREAM_ORIGIN=http://your.backend.server npm run build:cloudflare
```

This will generate two subdirectories in the frontend directory:

- dist: Static resource site
- functions: Reverse proxy implemented with workers

Deploy with wrangler:

```bash
# 1) Login
npx wrangler login

# 2) Create Pages project
npx wrangler pages project create <PROJECT_NAME>

# 3) Deploy: execute in project root (must have functions/)
npx wrangler pages deploy dist --project-name <PROJECT_NAME>
```

Note: Backend server nginx forwarding settings remain unchanged, same as [NGINX Hosting Static Pages, Forwarding API to HTTP Backend (Recommended)](#nginx-hosting-static-pages-forwarding-api-to-http-backend-recommended).

Optional: Custom domain

## Cloudflare Pages + Cloudflare Tunnel

Goal: **No public server/public IP required** (backend on intranet/NAT/campus network), achieve **frontend-backend separation deployment**:

- Frontend: Cloudflare Pages hosts static site
- Backend: Run `server/bin/server` on any machine with outbound internet access (only needs to listen on local port)
- Connection: Use Cloudflare Tunnel to expose the backend to a Cloudflare subdomain; frontend uses Pages Functions for same-origin proxy of `/api/*`, `/d/*` to that subdomain (no CORS issues on browser side)

Note: The backend host must be able to connect to Cloudflare tunnel service. See [Connectivity pre-checks](https://developers.cloudflare.com/cloudflare-one/networks/connectors/cloudflare-tunnel/troubleshoot-tunnels/connectivity-prechecks/). If the intranet server is blocked from these services or cannot access the internet, this solution cannot be used.

### Backend: Start Service (Intranet Machine is Fine)

Backend data preparation and `server.yaml` configuration are the same as other methods, then start on the backend machine:

```bash
server/bin/server -config "data/server.yaml"
```

Default backend listens on `:8080` (or as configured). **No need** to configure a public IP for this machine, and no need for additional nginx reverse proxy (Tunnel handles the entry point).

### Backend: Create and Run Cloudflare Tunnel

Download cloudflared from GitHub:

```bash
wget https://github.com/cloudflare/cloudflared/releases/download/2026.1.1/cloudflared-linux-amd64
mv cloudflared-linux-amd64 cloudflared
chmod +x cloudflared
```

See [Create a locally-managed tunnel](https://developers.cloudflare.com/cloudflare-one/networks/connectors/cloudflare-tunnel/do-more-with-tunnels/local-management/create-local-tunnel/)

After installation, login and create a Tunnel. A webpage will pop up asking you to select a subdomain under a Cloudflare-hosted domain as the upstream entry:

```bash
cloudflared tunnel login
```

Create a tunnel (name it as you like, here I'm using atlasmap):

```bash
cloudflared tunnel create atlasmap
```

In the following steps, `<USER>`, `TUNNEL_ID`, and `HOST_DOMAIN` need to be replaced with actual values.

Create `cloudflared` configuration file (example path `~/.cloudflared/config.yml`):

```yaml
tunnel: atlasmap
credentials-file: /home/<USER>/.cloudflared/<TUNNEL_ID>.json

ingress:
  - hostname: <HOST_DOMAIN>
    service: http://127.0.0.1:8080
  - service: http_status:404
```

Create DNS record (very important, can also be created manually in Cloudflare dashboard DNS):

```bash
cloudflared tunnel route dns atlasmap <HOST_DOMAIN>
```

Start the Tunnel:

```bash
cloudflared tunnel run atlasmap
```

As long as no ERR appears in the logs, it's running normally.

Optional: Configure Tunnel as a systemd service to ensure auto-start after reboot (follow cloudflared official documentation).

### Frontend: Build Cloudflare Pages Artifacts (Including Functions Proxy)

This project provides a Pages Functions same-origin proxy implementation: after deploying to Pages, `/api/*` and `/d/*` will be forwarded to the upstream specified by `UPSTREAM_ORIGIN`.

Build in the `frontend/` directory (point upstream to the Tunnel subdomain just created):

```bash
cd frontend
UPSTREAM_ORIGIN=<HOST_DOMAIN> npm run build:cloudflare
```

After building:

- `frontend/dist/`: Static site
- `frontend/functions/`: Pages Functions (same-origin reverse proxy for `/api/*`, `/d/*`)

### Deploy to Cloudflare Pages

Deploy with wrangler:

```bash
cd frontend
npx wrangler login
npx wrangler pages project create <PROJECT_NAME>
npx wrangler pages deploy dist --project-name <PROJECT_NAME>
```

### Verification and Troubleshooting

- Visit Pages domain (e.g., `https://<PROJECT>.pages.dev`)
- Open browser Network: confirm requests hit same-origin `/api/...`, `/d/...` (forwarded by Pages Functions to `https://api.example.com`)
- If API returns 502/timeout:
  - Confirm `server` is running on the backend machine, and can access `http://127.0.0.1:8080/health` locally (or your configured health check)
  - Confirm `cloudflared tunnel run ...` is running, and `api.example.com` is active in Cloudflare DNS
  - Confirm `UPSTREAM_ORIGIN` uses `https://api.example.com` (Tunnel provides HTTPS entry by default)
