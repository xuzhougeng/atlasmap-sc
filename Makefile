.PHONY: all build run dev dev-soma test clean preprocess preprocess-venv ensure-go

# Variables
GO_SERVER_DIR := server
FRONTEND_DIR := frontend
PREPROCESS_DIR := preprocessing
DATA_DIR := data
SERVER_CONFIG ?= ../config/server.yaml

# Python (preprocessing) environment (uv)
UV ?= uv
PREPROCESS_VENV_DIR := $(PREPROCESS_DIR)/.venv
PREPROCESS_PYTHON ?= 3.11
UV_LINK_MODE ?= copy

# Preprocess input/output paths (relative to repo root, unless absolute path)
# When running from preprocessing/, we need to prepend ../ for relative paths
OUTPUT ?= $(DATA_DIR)/preprocessed
PREPROCESS_INPUT := $(if $(filter /%,$(INPUT)),$(INPUT),../$(INPUT))
PREPROCESS_OUTPUT := $(if $(filter /%,$(OUTPUT)),$(OUTPUT),../$(OUTPUT))

# Go bootstrap (only used when system Go is missing)
# NOTE: server/go.mod specifies:
#   go 1.24.0
#   toolchain go1.24.11
GO_TOOLCHAIN_VERSION ?= 1.24.11
GO_TOOLS_DIR ?= .tools
GO_BOOTSTRAP_DIR := $(abspath $(GO_TOOLS_DIR)/go$(GO_TOOLCHAIN_VERSION))
GO_ENV := PATH="$(GO_BOOTSTRAP_DIR)/bin:$$PATH"

# Default target
all: build

# Build all components
build: build-server build-frontend

# Build Go server
build-server:
	cd $(GO_SERVER_DIR) && $(GO_ENV) go build -o bin/server ./cmd/server

# Build frontend
build-frontend:
	cd $(FRONTEND_DIR) && npm ci && npm run build

# Run development servers
dev:
	@echo "Starting development servers..."
	@$(MAKE) -j2 dev-server dev-frontend

dev-server:
	cd $(GO_SERVER_DIR) && $(GO_ENV) go run ./cmd/server -config "$(SERVER_CONFIG)"

dev-soma:
	@echo "Starting SOMA-enabled development servers..."
	@$(MAKE) -j2 dev-soma-server dev-frontend

dev-soma-server:
	cd $(GO_SERVER_DIR) && \
		CGO_ENABLED=1 \
		CGO_CFLAGS="-I$$CONDA_PREFIX/include" \
		CGO_LDFLAGS="-L$$CONDA_PREFIX/lib -ltiledb -Wl,-rpath,$$CONDA_PREFIX/lib" \
		$(GO_ENV) go run -tags soma ./cmd/server -config "$(SERVER_CONFIG)"

dev-frontend:
	cd $(FRONTEND_DIR) && npm run dev

# Run tests
test: test-server test-frontend

test-server:
	cd $(GO_SERVER_DIR) && $(GO_ENV) go test ./...

test-frontend:
	cd $(FRONTEND_DIR) && npm test

# Run preprocessing
preprocess: preprocess-venv
	@echo "Running preprocessing pipeline..."
	@if [ -z "$(INPUT)" ]; then \
		echo "Usage: make preprocess INPUT=path/to/data.h5ad [OUTPUT=data/dataset_x]"; \
		exit 1; \
	fi
	cd $(PREPROCESS_DIR) && UV_LINK_MODE=$(UV_LINK_MODE) $(UV) run -m atlasmap_preprocess.cli run \
		--input $(PREPROCESS_INPUT) \
		--output $(PREPROCESS_OUTPUT) \
		--zoom-levels 8 \
		--n-genes 500

# Create preprocessing venv (Python>=3.9 required by preprocessing/pyproject.toml)
preprocess-venv:
	@command -v "$(UV)" >/dev/null 2>&1 || { \
		echo "uv not found. Install it first: https://docs.astral.sh/uv/"; \
		exit 1; \
	}
	@cd $(PREPROCESS_DIR) && \
		[ -d ".venv" ] || ( \
			echo "Creating preprocessing venv with Python $(PREPROCESS_PYTHON)..." && \
			$(UV) python install $(PREPROCESS_PYTHON) && \
			$(UV) venv --python $(PREPROCESS_PYTHON) \
		)
	@cd $(PREPROCESS_DIR) && \
		echo "Installing preprocessing deps into .venv..." && \
		UV_LINK_MODE=$(UV_LINK_MODE) $(UV) pip install -e ".[dev]" --link-mode=$(UV_LINK_MODE)

# Install Python dependencies
install-python:
	@$(MAKE) preprocess-venv

# Install Go dependencies
ensure-go:
	@set -e; \
	if command -v go >/dev/null 2>&1; then \
		echo "Found Go: $$(go version)"; \
	else \
		echo "Go not found. Bootstrapping Go $(GO_TOOLCHAIN_VERSION) into $(GO_BOOTSTRAP_DIR)"; \
		OS=$$(uname -s | tr '[:upper:]' '[:lower:]'); \
		ARCH=$$(uname -m); \
		case "$$ARCH" in \
			x86_64|amd64) ARCH=amd64 ;; \
			aarch64|arm64) ARCH=arm64 ;; \
			*) echo "Unsupported architecture: $$ARCH"; exit 1 ;; \
		esac; \
		URL="https://go.dev/dl/go$(GO_TOOLCHAIN_VERSION).$$OS-$$ARCH.tar.gz"; \
		TMPDIR=$$(mktemp -d); \
		echo "Downloading $$URL"; \
		if command -v curl >/dev/null 2>&1; then \
			curl -fsSL "$$URL" -o "$$TMPDIR/go.tgz"; \
		elif command -v wget >/dev/null 2>&1; then \
			wget -qO "$$TMPDIR/go.tgz" "$$URL"; \
		else \
			echo "Neither curl nor wget is available. Please install one of them to bootstrap Go."; \
			exit 1; \
		fi; \
		tar -C "$$TMPDIR" -xzf "$$TMPDIR/go.tgz"; \
		rm -rf "$(GO_BOOTSTRAP_DIR)"; \
		mkdir -p "$(GO_BOOTSTRAP_DIR)"; \
		mv "$$TMPDIR/go"/* "$(GO_BOOTSTRAP_DIR)/"; \
		rm -rf "$$TMPDIR"; \
		echo "Bootstrapped Go: $(GO_BOOTSTRAP_DIR)/bin/go"; \
		"$(GO_BOOTSTRAP_DIR)/bin/go" version; \
	fi

install-go: ensure-go
	cd $(GO_SERVER_DIR) && $(GO_ENV) go mod download

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
	cd $(GO_SERVER_DIR) && $(GO_ENV) go fmt ./...
	cd $(PREPROCESS_DIR) && UV_LINK_MODE=$(UV_LINK_MODE) $(UV) run black .
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
