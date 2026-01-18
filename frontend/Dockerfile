# AtlasMap Frontend Dockerfile
# Multi-stage build

# Build stage
FROM node:20-alpine AS builder

WORKDIR /app

# Copy package files
COPY package.json package-lock.json* ./

# Install dependencies
RUN npm ci

# Copy source
COPY . .

# Build
RUN npm run build

# Production stage
FROM nginx:alpine

# Copy built assets
COPY --from=builder /app/dist /usr/share/nginx/html

# Default backend target (can override at `docker run -e ...`)
ENV BACKEND_HOST=server
ENV BACKEND_PORT=8080

# Copy nginx config template (rendered by nginx docker entrypoint)
COPY nginx.conf.template /etc/nginx/templates/default.conf.template

EXPOSE 80

CMD ["nginx", "-g", "daemon off;"]
