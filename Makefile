.PHONY: all build run dev dev-soma test clean preprocess preprocess-venv ensure-go \
	check-node check-npm prepare-go prepare-python prepare-frontend prepare

# Variables
GO_SERVER_DIR := server
FRONTEND_DIR := frontend
PREPROCESS_DIR := preprocessing
DATA_DIR := data
SERVER_CONFIG ?= ../config/server.yaml
NODE ?= node
NPM ?= npm
TRASH_DIR ?= .Trash

# Python (preprocessing) environment (uv)
UV_BIN ?= $(HOME)/.local/bin/uv
UV ?= uv
UV_INSTALL_URL ?= https://astral.sh/uv/install.sh
PREPROCESS_VENV_DIR := $(PREPROCESS_DIR)/.venv
PREPROCESS_PYTHON ?= 3.11
UV_LINK_MODE ?= copy

# Preprocessing defaults (can be overridden, e.g. `make preprocess PREPROCESS_NO_SOMA=0 ZOOM_LEVEL=8`)
# - PREPROCESS_NO_SOMA=1: skip writing TileDBSOMA store (faster, smaller output)
PREPROCESS_NO_SOMA ?= 1
# - ZOOM_LEVEL: required. Number of zoom levels to generate (passed to preprocessing as --zoom-levels).
ZOOM_LEVEL ?=
PREPROCESS_N_GENES ?= 500
# - PREPROCESS_ALL_EXPRESSED=1: use all expressed genes (skip HVG/top-N)
# - PREPROCESS_MIN_CELLS=1: include genes expressed in >= N cells (only with all-expressed)
PREPROCESS_ALL_EXPRESSED ?= 1
PREPROCESS_MIN_CELLS ?= 1
# - PREPROCESS_EXCLUDE_CATEGORIES: space-separated list of category columns to exclude
#   Example: PREPROCESS_EXCLUDE_CATEGORIES="observation_joinid barcode"
PREPROCESS_EXCLUDE_CATEGORIES ?=
# - PREPROCESS_MAX_CATEGORY_CARDINALITY: max unique values for category columns (default 2000)
PREPROCESS_MAX_CATEGORY_CARDINALITY ?= 2000

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

# Node bootstrap (used when system Node is missing or unsupported)
# Frontend requires Node ^18.0.0 or >=20.0.0 (Vite 5)
NODE_TOOLCHAIN_VERSION ?= 20.11.1
NODE_BOOTSTRAP_DIR := $(abspath $(GO_TOOLS_DIR)/node-v$(NODE_TOOLCHAIN_VERSION))
NODE_ENV := PATH="$(NODE_BOOTSTRAP_DIR)/bin:$$PATH"

# Default target
all: build

# Build all components
build: build-server build-frontend

# Build Go server
build-server: ensure-go
	cd $(GO_SERVER_DIR) && $(GO_ENV) go build -o bin/server ./cmd/server

# Build frontend
build-frontend: prepare-frontend
	cd $(FRONTEND_DIR) && $(NODE_ENV) npm run build

# Run development servers
dev:
	@echo "Starting development servers..."
	@$(MAKE) -j2 dev-server dev-frontend

dev-server: ensure-go
	cd $(GO_SERVER_DIR) && $(GO_ENV) go run ./cmd/server -config "$(SERVER_CONFIG)"

dev-soma:
	@echo "Starting SOMA-enabled development servers..."
	@$(MAKE) -j2 dev-soma-server dev-frontend

dev-soma-server: ensure-go
	cd $(GO_SERVER_DIR) && \
		CGO_ENABLED=1 \
		CGO_CFLAGS="-I$$CONDA_PREFIX/include" \
		CGO_LDFLAGS="-L$$CONDA_PREFIX/lib -ltiledb -Wl,-rpath,$$CONDA_PREFIX/lib" \
		$(GO_ENV) go run -tags soma ./cmd/server -config "$(SERVER_CONFIG)"

dev-frontend: prepare-frontend
	cd $(FRONTEND_DIR) && $(NODE_ENV) npm run dev

# Run tests
test: test-server test-frontend

test-server: ensure-go
	cd $(GO_SERVER_DIR) && $(GO_ENV) go test ./...

test-frontend: prepare-frontend
	cd $(FRONTEND_DIR) && $(NODE_ENV) npm test

# Run preprocessing
preprocess: preprocess-venv
	@echo "Running preprocessing pipeline..."
	@if [ -z "$(INPUT)" ]; then \
		echo "Usage: make preprocess INPUT=path/to/data.h5ad OUTPUT=data/dataset_x ZOOM_LEVEL=8 [PREPROCESS_NO_SOMA=0]"; \
		exit 1; \
	fi
	@if [ -z "$(ZOOM_LEVEL)" ]; then \
		echo "Missing required var: ZOOM_LEVEL (example: ZOOM_LEVEL=8)"; \
		exit 1; \
	fi
	@set -e; \
	ZOOM="$(ZOOM_LEVEL)"; \
	NO_SOMA_FLAG=""; \
	case "$(PREPROCESS_NO_SOMA)" in \
		1|true|TRUE|yes|YES) NO_SOMA_FLAG="--no-soma" ;; \
	esac; \
	ALL_EXPRESSED_FLAG=""; \
	MIN_CELLS_ARGS=""; \
	case "$(PREPROCESS_ALL_EXPRESSED)" in \
		1|true|TRUE|yes|YES) ALL_EXPRESSED_FLAG="--all-expressed"; MIN_CELLS_ARGS="--min-cells $(PREPROCESS_MIN_CELLS)" ;; \
	esac; \
	EXCLUDE_CAT_ARGS=""; \
	for col in $(PREPROCESS_EXCLUDE_CATEGORIES); do \
		EXCLUDE_CAT_ARGS="$$EXCLUDE_CAT_ARGS --exclude-category $$col"; \
	done; \
	MAX_CAT_CARD_ARG="--max-category-cardinality $(PREPROCESS_MAX_CATEGORY_CARDINALITY)"; \
	echo "Preprocess: zoom_levels=$$ZOOM, all_expressed=$(PREPROCESS_ALL_EXPRESSED), min_cells=$(PREPROCESS_MIN_CELLS), n_genes=$(PREPROCESS_N_GENES), no_soma=$(PREPROCESS_NO_SOMA), max_cat_card=$(PREPROCESS_MAX_CATEGORY_CARDINALITY)"; \
	cd $(PREPROCESS_DIR) && PATH="$$HOME/.local/bin:$$PATH" UV_LINK_MODE=$(UV_LINK_MODE) $(UV) run -m atlasmap_preprocess.cli run \
		--input $(PREPROCESS_INPUT) \
		--output $(PREPROCESS_OUTPUT) \
		--zoom-levels "$$ZOOM" \
		--n-genes "$(PREPROCESS_N_GENES)" \
		$$ALL_EXPRESSED_FLAG $$MIN_CELLS_ARGS \
		$$EXCLUDE_CAT_ARGS $$MAX_CAT_CARD_ARG \
		$$NO_SOMA_FLAG

# Create preprocessing venv (Python>=3.9 required by preprocessing/pyproject.toml)
ensure-uv:
	@set -e; \
	if command -v "$(UV)" >/dev/null 2>&1; then \
		echo "Found uv: $$("$(UV)" --version)"; \
	elif [ -x "$(UV_BIN)" ]; then \
		echo "Found uv: $$("$(UV_BIN)" --version)"; \
	else \
		echo "uv not found. Installing via Astral installer..."; \
		if command -v curl >/dev/null 2>&1; then \
			curl -LsSf "$(UV_INSTALL_URL)" | sh; \
		elif command -v wget >/dev/null 2>&1; then \
			wget -qO- "$(UV_INSTALL_URL)" | sh; \
		else \
			echo "Neither curl nor wget is available to install uv."; \
			exit 1; \
		fi; \
		if [ -x "$(UV_BIN)" ]; then \
			echo "Installed uv: $$("$(UV_BIN)" --version)"; \
		else \
			echo "uv install finished, but $(UV_BIN) not found."; \
			echo "Try: export PATH=\"$$HOME/.local/bin:$$PATH\""; \
			exit 1; \
		fi; \
	fi

preprocess-venv: ensure-uv
	@cd $(PREPROCESS_DIR) && \
		[ -d ".venv" ] || ( \
			echo "Creating preprocessing venv with Python $(PREPROCESS_PYTHON)..." && \
			PATH="$$HOME/.local/bin:$$PATH" $(UV) python install $(PREPROCESS_PYTHON) && \
			PATH="$$HOME/.local/bin:$$PATH" $(UV) venv --python $(PREPROCESS_PYTHON) \
		)
	@cd $(PREPROCESS_DIR) && \
		echo "Installing preprocessing deps into .venv..." && \
		PATH="$$HOME/.local/bin:$$PATH" UV_LINK_MODE=$(UV_LINK_MODE) $(UV) pip install -e ".[dev]" --link-mode=$(UV_LINK_MODE)

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
		if [ -e "$(GO_BOOTSTRAP_DIR)" ]; then \
			mkdir -p "$(TRASH_DIR)"; \
			TS=$$(date +%Y%m%d-%H%M%S); \
			mv "$(GO_BOOTSTRAP_DIR)" "$(TRASH_DIR)/go$(GO_TOOLCHAIN_VERSION)-$$TS"; \
		fi; \
		mkdir -p "$(GO_BOOTSTRAP_DIR)"; \
		mv "$$TMPDIR/go"/* "$(GO_BOOTSTRAP_DIR)/"; \
		rm -rf "$$TMPDIR"; \
		echo "Bootstrapped Go: $(GO_BOOTSTRAP_DIR)/bin/go"; \
		"$(GO_BOOTSTRAP_DIR)/bin/go" version; \
	fi

install-go: ensure-go
	cd $(GO_SERVER_DIR) && $(GO_ENV) go mod download

# Frontend runtime checks
check-node:
	@set -e; \
	NEED_BOOTSTRAP=1; \
	if command -v "$(NODE)" >/dev/null 2>&1; then \
		NODEV=$$($(NODE) -v | sed 's/^v//'); \
		MAJOR=$${NODEV%%.*}; \
		if [ "$$MAJOR" -eq 18 ] || [ "$$MAJOR" -ge 20 ]; then \
			echo "Found Node: v$$NODEV"; \
			NEED_BOOTSTRAP=0; \
		else \
			echo "Unsupported Node version: v$$NODEV"; \
			echo "Frontend requires Node ^18.0.0 or >=20.0.0 (Vite 5)."; \
			echo "Will bootstrap Node $(NODE_TOOLCHAIN_VERSION) into $(NODE_BOOTSTRAP_DIR)"; \
		fi; \
	else \
		echo "node not found. Will bootstrap Node $(NODE_TOOLCHAIN_VERSION) into $(NODE_BOOTSTRAP_DIR)"; \
	fi; \
	if [ "$$NEED_BOOTSTRAP" -eq 1 ]; then \
		OS=$$(uname -s | tr '[:upper:]' '[:lower:]'); \
		ARCH=$$(uname -m); \
		case "$$OS" in \
			linux|darwin) ;; \
			*) echo "Unsupported OS for Node bootstrap: $$OS"; exit 1 ;; \
		esac; \
		case "$$ARCH" in \
			x86_64|amd64) ARCH=x64 ;; \
			aarch64|arm64) ARCH=arm64 ;; \
			*) echo "Unsupported architecture for Node bootstrap: $$ARCH"; exit 1 ;; \
		esac; \
		URL="https://nodejs.org/dist/v$(NODE_TOOLCHAIN_VERSION)/node-v$(NODE_TOOLCHAIN_VERSION)-$$OS-$$ARCH.tar.xz"; \
		TMPDIR=$$(mktemp -d); \
		echo "Downloading $$URL"; \
		if command -v curl >/dev/null 2>&1; then \
			curl -fsSL "$$URL" -o "$$TMPDIR/node.tar.xz"; \
		elif command -v wget >/dev/null 2>&1; then \
			wget -qO "$$TMPDIR/node.tar.xz" "$$URL"; \
		else \
			echo "Neither curl nor wget is available. Please install one of them to bootstrap Node."; \
			exit 1; \
		fi; \
		tar -C "$$TMPDIR" -xJf "$$TMPDIR/node.tar.xz"; \
		if [ -e "$(NODE_BOOTSTRAP_DIR)" ]; then \
			mkdir -p "$(TRASH_DIR)"; \
			TS=$$(date +%Y%m%d-%H%M%S); \
			mv "$(NODE_BOOTSTRAP_DIR)" "$(TRASH_DIR)/node-v$(NODE_TOOLCHAIN_VERSION)-$$TS"; \
		fi; \
		mkdir -p "$(NODE_BOOTSTRAP_DIR)"; \
		mv "$$TMPDIR/node-v$(NODE_TOOLCHAIN_VERSION)-$$OS-$$ARCH"/* "$(NODE_BOOTSTRAP_DIR)/"; \
		rm -rf "$$TMPDIR"; \
		echo "Bootstrapped Node: $(NODE_BOOTSTRAP_DIR)/bin/node"; \
		"$(NODE_BOOTSTRAP_DIR)/bin/node" -v; \
		"$(NODE_BOOTSTRAP_DIR)/bin/npm" -v; \
	fi

check-npm:
	@echo "Found npm: $$($(NODE_ENV) npm -v)"

# Install frontend dependencies
install-frontend: check-node check-npm
	cd $(FRONTEND_DIR) && $(NODE_ENV) npm ci

# Prepare targets (recommended entrypoint for fresh clones / CI)
prepare-go: install-go

prepare-python: preprocess-venv

prepare-frontend: install-frontend

prepare: prepare-go prepare-python prepare-frontend
	@echo "Prepared: Go toolchain + server deps, Python preprocessing venv, frontend deps."

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
fmt: ensure-go
	cd $(GO_SERVER_DIR) && $(GO_ENV) go fmt ./...
	cd $(PREPROCESS_DIR) && UV_LINK_MODE=$(UV_LINK_MODE) $(UV) run black .
	cd $(FRONTEND_DIR) && $(NODE_ENV) npm run lint -- --fix

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
	@echo "  make prepare        - Prepare Go/Python/Node dependencies (recommended)"
	@echo "  make preprocess INPUT=file.h5ad - Run preprocessing"
	@echo "  make docker-build   - Build Docker images"
	@echo "  make docker-up      - Start Docker containers"
	@echo "  make clean          - Clean build artifacts"
	@echo "  make help           - Show this help"
