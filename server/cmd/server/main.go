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
	"path/filepath"
	"syscall"
	"time"

	"github.com/soma-tiles/server/internal/api"
	"github.com/soma-tiles/server/internal/cache"
	"github.com/soma-tiles/server/internal/config"
	"github.com/soma-tiles/server/internal/data/soma"
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

		var somaReader *soma.Reader
		if ds.SomaPath != "" {
			r, err := soma.NewReader(ds.SomaPath)
			if err != nil {
				log.Printf("  [%s] SOMA not initialized: %v", datasetID, err)
			} else {
				somaReader = r
				log.Printf("  [%s] SOMA experiment: %s (supported=%v)", datasetID, somaReader.ExperimentURI(), somaReader.Supported())
			}
		}

		// Initialize tile services for each coordinate system (if present).
		md := zarrReader.Metadata()
		coordSystems := md.CoordinateSystems
		defaultCoord := md.DefaultCoordinateSystem
		if defaultCoord == "" {
			defaultCoord = md.UMAPKey
		}

		metadataDir := filepath.Dir(ds.ZarrPath)
		dsZarrPath := filepath.Clean(ds.ZarrPath)

		// Backward-compat: if metadata doesn't list coordinate systems, register a single one.
		if len(coordSystems) == 0 {
			key := defaultCoord
			if key == "" {
				key = "X_umap"
			}
			rel, err := filepath.Rel(metadataDir, dsZarrPath)
			if err != nil || rel == "." || rel == "" || rel == ".." || rel[:1] == string(filepath.Separator) {
				rel = filepath.Base(dsZarrPath)
			}
			coordSystems = []zarr.CoordinateSystem{{Key: key, ZarrPath: rel}}
			defaultCoord = key
		}
		if defaultCoord == "" && len(coordSystems) > 0 {
			defaultCoord = coordSystems[0].Key
		}

		for _, cs := range coordSystems {
			if cs.Key == "" {
				log.Fatalf("  [%s] Invalid coordinate system entry with empty key", datasetID)
			}
			if cs.ZarrPath == "" {
				log.Fatalf("  [%s] Invalid coordinate system entry %q with empty zarr_path", datasetID, cs.Key)
			}

			storePath := cs.ZarrPath
			if !filepath.IsAbs(storePath) {
				storePath = filepath.Join(metadataDir, storePath)
			}
			storePath = filepath.Clean(storePath)

			reader := zarrReader
			if storePath != dsZarrPath {
				r, err := zarr.NewReader(storePath)
				if err != nil {
					log.Fatalf("Failed to initialize Zarr reader for dataset %q coord %q: %v", datasetID, cs.Key, err)
				}
				reader = r
				log.Printf("  [%s] Coord %q loaded from: %s", datasetID, cs.Key, storePath)
			} else {
				log.Printf("  [%s] Coord %q uses primary Zarr store", datasetID, cs.Key)
			}

			tileService := service.NewTileService(service.TileServiceConfig{
				DatasetID:  datasetID + "|" + cs.Key, // isolate cache namespace across coords
				CoordKey:   cs.Key,
				ZarrReader: reader,
				SomaReader: somaReader,
				Cache:      cacheManager,
				Renderer:   tileRenderer,
			})

			registry.RegisterCoord(datasetID, cs.Key, tileService)
		}
		registry.SetDefaultCoord(datasetID, defaultCoord)

		if ds.H5ADPath != "" {
			registry.SetH5ADPath(datasetID, ds.H5ADPath)
		}
	}

	// Initialize job manager for DE jobs (SQLite persistence)
	jobManager, err := api.NewJobManager(api.JobManagerConfig{
		MaxConcurrent: cfg.DE.MaxConcurrent,
		SQLitePath:    cfg.DE.SQLitePath,
		RetentionDays: cfg.DE.RetentionDays,
		CleanupPeriod: 1 * time.Hour,
	})
	if err != nil {
		log.Fatalf("Failed to initialize job manager: %v", err)
	}
	log.Printf("DE job manager: max_concurrent=%d, retention_days=%d, sqlite=%s",
		cfg.DE.MaxConcurrent, cfg.DE.RetentionDays, cfg.DE.SQLitePath)

	// Wire up DE service as job executor
	deService := service.NewDEService(registry)
	jobManager.Executor = deService.ExecuteDEJob

	jobManager.Start()
	defer jobManager.Stop()

	// Set up HTTP router
	router := api.NewRouter(api.RouterConfig{
		Registry:    registry,
		CORSOrigins: cfg.Server.CORSOrigins,
		JobManager:  jobManager,
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
