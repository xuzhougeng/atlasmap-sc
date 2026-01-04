package config

import (
	"os"
	"path/filepath"
	"testing"
)

func TestLoad_LegacyFormat(t *testing.T) {
	content := `
server:
  port: 9000
data:
  zarr_path: "/data/legacy/bins.zarr"
  soma_path: "/data/legacy/soma"
cache:
  tile_size_mb: 256
`
	cfg := loadFromString(t, content)

	if cfg.Server.Port != 9000 {
		t.Errorf("expected port 9000, got %d", cfg.Server.Port)
	}
	if cfg.Data.DefaultDataset != "default" {
		t.Errorf("expected default dataset 'default', got %q", cfg.Data.DefaultDataset)
	}
	ds, ok := cfg.Data.Datasets["default"]
	if !ok {
		t.Fatal("expected 'default' dataset")
	}
	if ds.ZarrPath != "/data/legacy/bins.zarr" {
		t.Errorf("unexpected zarr_path: %s", ds.ZarrPath)
	}
	if ds.SomaPath != "/data/legacy/soma" {
		t.Errorf("unexpected soma_path: %s", ds.SomaPath)
	}
}

func TestLoad_MultiDatasetFormat(t *testing.T) {
	content := `
server:
  port: 8080
data:
  pbmc:
    zarr_path: "/data/pbmc/bins.zarr"
    soma_path: "/data/pbmc/soma"
  liver:
    zarr_path: "/data/liver/bins.zarr"
    soma_path: "/data/liver/soma"
`
	cfg := loadFromString(t, content)

	if len(cfg.Data.Datasets) != 2 {
		t.Fatalf("expected 2 datasets, got %d", len(cfg.Data.Datasets))
	}

	// First dataset in YAML order should be default
	if cfg.Data.DefaultDataset != "pbmc" {
		t.Errorf("expected default dataset 'pbmc', got %q", cfg.Data.DefaultDataset)
	}

	pbmc, ok := cfg.Data.Datasets["pbmc"]
	if !ok {
		t.Fatal("expected 'pbmc' dataset")
	}
	if pbmc.ZarrPath != "/data/pbmc/bins.zarr" {
		t.Errorf("unexpected pbmc zarr_path: %s", pbmc.ZarrPath)
	}

	liver, ok := cfg.Data.Datasets["liver"]
	if !ok {
		t.Fatal("expected 'liver' dataset")
	}
	if liver.ZarrPath != "/data/liver/bins.zarr" {
		t.Errorf("unexpected liver zarr_path: %s", liver.ZarrPath)
	}

	// Check order preserved
	ids := cfg.Data.DatasetIDs()
	if len(ids) != 2 || ids[0] != "pbmc" || ids[1] != "liver" {
		t.Errorf("unexpected dataset order: %v", ids)
	}
}

func TestLoad_DefaultsApplied(t *testing.T) {
	content := `
server:
  port: 0
data:
  test:
    zarr_path: "/test/bins.zarr"
`
	cfg := loadFromString(t, content)

	if cfg.Server.Port != 8080 {
		t.Errorf("expected default port 8080, got %d", cfg.Server.Port)
	}
	if cfg.Cache.TileSizeMB != 512 {
		t.Errorf("expected default cache size 512, got %d", cfg.Cache.TileSizeMB)
	}
	if cfg.Render.TileSize != 256 {
		t.Errorf("expected default tile size 256, got %d", cfg.Render.TileSize)
	}
}

func TestLoad_NoDataSection(t *testing.T) {
	content := `
server:
  port: 8080
`
	cfg := loadFromString(t, content)

	if cfg.Data.DefaultDataset != "default" {
		t.Errorf("expected default dataset, got %q", cfg.Data.DefaultDataset)
	}
	if len(cfg.Data.Datasets) != 1 {
		t.Errorf("expected 1 default dataset, got %d", len(cfg.Data.Datasets))
	}
}

func loadFromString(t *testing.T, content string) *Config {
	t.Helper()

	dir := t.TempDir()
	path := filepath.Join(dir, "config.yaml")
	if err := os.WriteFile(path, []byte(content), 0644); err != nil {
		t.Fatalf("failed to write temp config: %v", err)
	}

	cfg, err := Load(path)
	if err != nil {
		t.Fatalf("failed to load config: %v", err)
	}
	return cfg
}

