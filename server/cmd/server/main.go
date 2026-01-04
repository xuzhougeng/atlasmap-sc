// Package main is the entry point for the SOMA-Tiles server.
package main

import (
	"context"
	"flag"
	"fmt"
	"log"
	"net/http"
	"os"
	"os/signal"
	"syscall"
	"time"

	"github.com/soma-tiles/server/internal/api"
	"github.com/soma-tiles/server/internal/cache"
	"github.com/soma-tiles/server/internal/config"
	"github.com/soma-tiles/server/internal/data/zarr"
	"github.com/soma-tiles/server/internal/render"
	"github.com/soma-tiles/server/internal/service"
)

func main() {
	// Parse command line flags
	configPath := flag.String("config", "config/server.yaml", "Path to configuration file")
	flag.Parse()

	// Load configuration
	cfg, err := config.Load(*configPath)
	if err != nil {
		log.Fatalf("Failed to load configuration: %v", err)
	}

	log.Printf("Starting SOMA-Tiles server on port %d", cfg.Server.Port)

	// Initialize components
	ctx := context.Background()

	// Initialize cache manager (shared across all datasets)
	cacheManager, err := cache.NewManager(cache.Config{
		TileCacheSizeMB: cfg.Cache.TileSizeMB,
		TileTTL:         time.Duration(cfg.Cache.TileTTLMinutes) * time.Minute,
		QueryCacheSize:  1000,
	})
	if err != nil {
		log.Fatalf("Failed to initialize cache: %v", err)
	}

	// Initialize tile renderer (shared across all datasets)
	tileRenderer := render.NewTileRenderer(render.Config{
		TileSize:        cfg.Render.TileSize,
		DefaultColormap: cfg.Render.DefaultColormap,
	})

	// Initialize dataset registry
	datasetIDs := cfg.Data.DatasetIDs()
	registry := api.NewDatasetRegistry(cfg.Data.DefaultDataset, datasetIDs, cfg.Server.Title)

	log.Printf("Initializing %d dataset(s), default: %s", len(datasetIDs), cfg.Data.DefaultDataset)

	// Initialize each dataset
	for _, datasetID := range datasetIDs {
		ds := cfg.Data.Datasets[datasetID]

		zarrReader, err := zarr.NewReader(ds.ZarrPath)
		if err != nil {
			log.Fatalf("Failed to initialize Zarr reader for dataset %q: %v", datasetID, err)
		}

		log.Printf("  [%s] Loaded from: %s", datasetID, ds.ZarrPath)
		log.Printf("    Zoom levels: %d, Genes: %d", zarrReader.Metadata().ZoomLevels, zarrReader.Metadata().NGenes)

		tileService := service.NewTileService(service.TileServiceConfig{
			DatasetID:  datasetID,
			ZarrReader: zarrReader,
			Cache:      cacheManager,
			Renderer:   tileRenderer,
		})

		registry.Register(datasetID, tileService)
	}

	// Set up HTTP router
	router := api.NewRouter(api.RouterConfig{
		Registry:    registry,
		CORSOrigins: cfg.Server.CORSOrigins,
	})

	// Create HTTP server
	server := &http.Server{
		Addr:         fmt.Sprintf(":%d", cfg.Server.Port),
		Handler:      router,
		ReadTimeout:  30 * time.Second,
		WriteTimeout: 60 * time.Second,
		IdleTimeout:  120 * time.Second,
	}

	// Start server in goroutine
	go func() {
		log.Printf("Server listening on http://localhost:%d", cfg.Server.Port)
		if err := server.ListenAndServe(); err != nil && err != http.ErrServerClosed {
			log.Fatalf("Server failed: %v", err)
		}
	}()

	// Wait for interrupt signal
	quit := make(chan os.Signal, 1)
	signal.Notify(quit, syscall.SIGINT, syscall.SIGTERM)
	<-quit

	log.Println("Shutting down server...")

	// Graceful shutdown with timeout
	shutdownCtx, cancel := context.WithTimeout(ctx, 30*time.Second)
	defer cancel()

	if err := server.Shutdown(shutdownCtx); err != nil {
		log.Printf("Server forced to shutdown: %v", err)
	}

	log.Println("Server stopped")
}
