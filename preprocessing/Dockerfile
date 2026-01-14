# AtlasMap Preprocessing Dockerfile
FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy package files
COPY pyproject.toml ./
COPY atlasmap_preprocess/ ./atlasmap_preprocess/

# Install Python package
RUN pip install --no-cache-dir -e .

# Set entrypoint
ENTRYPOINT ["atlasmap-preprocess"]
CMD ["--help"]
