

容器构建

```bash
docker build -t atlasmap-frontend:latest ./frontend
docker tag atlasmap-frontend:latest hoptoper/atlasmap-frontend:latest

docker login
docker push hoptoper/atlasmap-frontend:latest
```

启动

```bash
docker run -d --name atlasmap-frontend \
  -p 3000:80 \
  -e BACKEND_HOST=xxx.xxx.xxx \
  -e BACKEND_PORT=80 \
  --restart unless-stopped \
  atlasmap-frontend:latest
```