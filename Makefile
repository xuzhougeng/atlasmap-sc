.PHONY: all build run dev dev-soma test clean preprocess

# Variables
GO_SERVER_DIR := server
FRONTEND_DIR := frontend
PREPROCESS_DIR := preprocessing
DATA_DIR := data
SERVER_CONFIG ?= ../config/server.yaml

# Default target
all: build

# Build all components
build: build-server build-frontend

# Build Go server
build-server:
	cd $(GO_SERVER_DIR) && go build -o bin/server ./cmd/server

# Build frontend
build-frontend:
	cd $(FRONTEND_DIR) && npm ci && npm run build

# Run development servers
dev:
	@echo "Starting development servers..."
	@$(MAKE) -j2 dev-server dev-frontend

dev-server:
	cd $(GO_SERVER_DIR) && go run ./cmd/server -config "$(SERVER_CONFIG)"

dev-soma:
	@echo "Starting SOMA-enabled development servers..."
	@$(MAKE) -j2 dev-soma-server dev-frontend

dev-soma-server:
	cd $(GO_SERVER_DIR) && \
		CGO_ENABLED=1 \
		CGO_CFLAGS="-I$$CONDA_PREFIX/include" \
		CGO_LDFLAGS="-L$$CONDA_PREFIX/lib -ltiledb -Wl,-rpath,$$CONDA_PREFIX/lib" \
		go run -tags soma ./cmd/server -config "$(SERVER_CONFIG)"

dev-frontend:
	cd $(FRONTEND_DIR) && npm run dev

# Run tests
test: test-server test-frontend

test-server:
	cd $(GO_SERVER_DIR) && go test ./...

test-frontend:
	cd $(FRONTEND_DIR) && npm test

# Run preprocessing
preprocess:
	@echo "Running preprocessing pipeline..."
	@if [ -z "$(INPUT)" ]; then \
		echo "Usage: make preprocess INPUT=path/to/data.h5ad"; \
		exit 1; \
	fi
	cd $(PREPROCESS_DIR) && python -m atlasmap_preprocess.cli run \
		--input $(INPUT) \
		--output ../$(DATA_DIR)/preprocessed \
		--zoom-levels 8 \
		--n-genes 500

# Install Python dependencies
install-python:
	cd $(PREPROCESS_DIR) && pip install -e ".[dev]"

# Install Go dependencies
install-go:
	cd $(GO_SERVER_DIR) && go mod download

# Install frontend dependencies
install-frontend:
	cd $(FRONTEND_DIR) && npm ci

# Install all dependencies
install: install-python install-go install-frontend

# Docker commands
docker-build:
	docker-compose build

docker-up:
	docker-compose up -d

docker-down:
	docker-compose down

docker-logs:
	docker-compose logs -f

docker-preprocess:
	docker-compose run --rm preprocess run \
		--input /data/raw/input.h5ad \
		--output /data/preprocessed

# Clean build artifacts
clean:
	rm -rf $(GO_SERVER_DIR)/bin
	rm -rf $(FRONTEND_DIR)/dist
	rm -rf $(FRONTEND_DIR)/node_modules
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type d -name "*.egg-info" -exec rm -rf {} +

# Format code
fmt:
	cd $(GO_SERVER_DIR) && go fmt ./...
	cd $(PREPROCESS_DIR) && black .
	cd $(FRONTEND_DIR) && npm run lint -- --fix

# Help
help:
	@echo "AtlasMap Makefile"
	@echo ""
	@echo "Usage:"
	@echo "  make build          - Build all components"
	@echo "  make dev [SERVER_CONFIG=../config/server.yaml]      - Start development servers"
	@echo "  make dev-soma [SERVER_CONFIG=../config/server.yaml] - Start SOMA-enabled dev server (requires conda TileDB)"
	@echo "  make test           - Run all tests"
	@echo "  make install        - Install all dependencies"
	@echo "  make preprocess INPUT=file.h5ad - Run preprocessing"
	@echo "  make docker-build   - Build Docker images"
	@echo "  make docker-up      - Start Docker containers"
	@echo "  make clean          - Clean build artifacts"
	@echo "  make help           - Show this help"
