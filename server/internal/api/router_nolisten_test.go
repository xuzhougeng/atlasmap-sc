package api

import (
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/soma-tiles/server/internal/cache"
	"github.com/soma-tiles/server/internal/data/zarr"
	"github.com/soma-tiles/server/internal/render"
	"github.com/soma-tiles/server/internal/service"
)

func localTestZarrPath(t *testing.T) string {
	t.Helper()

	wd, err := os.Getwd()
	if err != nil {
		t.Fatalf("failed to get working dir: %v", err)
	}

	// Package directory: server/internal/api
	// Repo root: ../../../
	p := filepath.Clean(filepath.Join(wd, "../../../data/preprocessed/zarr/bins.zarr"))
	if _, err := os.Stat(p); os.IsNotExist(err) {
		t.Skipf("test data not found at %s, skipping", p)
	}
	return p
}

func TestGeneStatsEndpoint_NoListen(t *testing.T) {
	zarrReader, err := zarr.NewReader(localTestZarrPath(t))
	if err != nil {
		t.Fatalf("Failed to initialize Zarr reader: %v", err)
	}
	defer zarrReader.Close()

	cacheManager, err := cache.NewManager(cache.Config{
		TileCacheSizeMB: 16,
		TileTTL:         1 * time.Minute,
		QueryCacheSize:  10,
	})
	if err != nil {
		t.Fatalf("Failed to initialize cache: %v", err)
	}
	defer cacheManager.Close()

	tileRenderer := render.NewTileRenderer(render.Config{
		TileSize:        256,
		DefaultColormap: "viridis",
	})

	tileService := service.NewTileService(service.TileServiceConfig{
		DatasetID:  "default",
		ZarrReader: zarrReader,
		Cache:      cacheManager,
		Renderer:   tileRenderer,
	})

	// Create registry with single dataset
	registry := NewDatasetRegistry("default", []string{"default"}, "")
	registry.Register("default", tileService)

	router := NewRouter(RouterConfig{
		Registry:    registry,
		CORSOrigins: []string{"http://localhost:3000"},
	})

	metadata := zarrReader.Metadata()
	if len(metadata.Genes) == 0 {
		t.Skip("No genes available in test data")
	}
	gene := metadata.Genes[0]
	if _, ok := metadata.GeneIndex["ENSG00000167286"]; ok {
		gene = "ENSG00000167286"
	}

	req := httptest.NewRequest(http.MethodGet, "/api/genes/"+gene+"/stats", nil)
	rec := httptest.NewRecorder()
	router.ServeHTTP(rec, req)

	if rec.Code != http.StatusOK {
		t.Fatalf("expected %d, got %d: %s", http.StatusOK, rec.Code, rec.Body.String())
	}

	var payload map[string]any
	if err := json.Unmarshal(rec.Body.Bytes(), &payload); err != nil {
		t.Fatalf("failed to decode JSON: %v", err)
	}
	if got, _ := payload["gene"].(string); got != gene {
		t.Fatalf("unexpected gene: got %q want %q", got, gene)
	}
}
