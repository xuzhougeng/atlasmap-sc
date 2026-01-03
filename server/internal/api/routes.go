// Package api provides HTTP handlers for the SOMA-Tiles server.
package api

import (
	"encoding/json"
	"net/http"
	"strconv"

	"github.com/go-chi/chi/v5"
	"github.com/go-chi/chi/v5/middleware"
	"github.com/go-chi/cors"
	"github.com/soma-tiles/server/internal/data/zarr"
	"github.com/soma-tiles/server/internal/service"
)

// RouterConfig contains router configuration.
type RouterConfig struct {
	TileService *service.TileService
	ZarrReader  *zarr.Reader
	CORSOrigins []string
}

// NewRouter creates a new HTTP router.
func NewRouter(cfg RouterConfig) *chi.Mux {
	r := chi.NewRouter()

	// Middleware
	r.Use(middleware.Logger)
	r.Use(middleware.Recoverer)
	r.Use(middleware.Compress(5))

	// CORS
	r.Use(cors.Handler(cors.Options{
		AllowedOrigins:   cfg.CORSOrigins,
		AllowedMethods:   []string{"GET", "POST", "PUT", "DELETE", "OPTIONS"},
		AllowedHeaders:   []string{"Accept", "Authorization", "Content-Type"},
		ExposedHeaders:   []string{"Link"},
		AllowCredentials: true,
		MaxAge:           300,
	}))

	// Health check
	r.Get("/health", func(w http.ResponseWriter, r *http.Request) {
		w.WriteHeader(http.StatusOK)
		w.Write([]byte("OK"))
	})

	// Tile endpoints
	r.Get("/tiles/{z}/{x}/{y}.png", tileHandler(cfg.TileService))
	r.Get("/tiles/{z}/{x}/{y}/expression/{gene}.png", expressionTileHandler(cfg.TileService))
	r.Get("/tiles/{z}/{x}/{y}/category/{column}.png", categoryTileHandler(cfg.TileService))

	// API endpoints
	r.Route("/api", func(r chi.Router) {
		r.Get("/metadata", metadataHandler(cfg.TileService))
		r.Get("/genes", genesHandler(cfg.TileService))
		r.Get("/genes/{gene}", geneInfoHandler(cfg.TileService))
		r.Get("/genes/{gene}/bins", geneBinsHandler(cfg.TileService))
		r.Get("/genes/{gene}/stats", geneStatsHandler(cfg.TileService))
		r.Get("/categories", categoriesHandler(cfg.TileService))
		r.Get("/categories/{column}/colors", categoryColorsHandler(cfg.TileService))
		r.Get("/categories/{column}/legend", categoryLegendHandler(cfg.TileService))
		r.Get("/stats", statsHandler(cfg.TileService))
	})

	return r
}

func tileHandler(svc *service.TileService) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		z, err := strconv.Atoi(chi.URLParam(r, "z"))
		if err != nil {
			http.Error(w, "invalid z", http.StatusBadRequest)
			return
		}
		x, err := strconv.Atoi(chi.URLParam(r, "x"))
		if err != nil {
			http.Error(w, "invalid x", http.StatusBadRequest)
			return
		}
		y, err := strconv.Atoi(chi.URLParam(r, "y"))
		if err != nil {
			http.Error(w, "invalid y", http.StatusBadRequest)
			return
		}

		data, err := svc.GetTile(z, x, y)
		if err != nil {
			// Return empty tile on error
			data, _ = svc.GetEmptyTile()
		}

		w.Header().Set("Content-Type", "image/png")
		w.Header().Set("Cache-Control", "public, max-age=3600")
		w.Write(data)
	}
}

func expressionTileHandler(svc *service.TileService) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		z, _ := strconv.Atoi(chi.URLParam(r, "z"))
		x, _ := strconv.Atoi(chi.URLParam(r, "x"))
		y, _ := strconv.Atoi(chi.URLParam(r, "y"))
		gene := chi.URLParam(r, "gene")
		colormap := r.URL.Query().Get("colormap")
		if colormap == "" {
			colormap = "viridis"
		}

		data, err := svc.GetExpressionTile(z, x, y, gene, colormap)
		if err != nil {
			data, _ = svc.GetEmptyTile()
		}

		w.Header().Set("Content-Type", "image/png")
		w.Header().Set("Cache-Control", "public, max-age=3600")
		w.Write(data)
	}
}

func categoryTileHandler(svc *service.TileService) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		z, _ := strconv.Atoi(chi.URLParam(r, "z"))
		x, _ := strconv.Atoi(chi.URLParam(r, "x"))
		y, _ := strconv.Atoi(chi.URLParam(r, "y"))
		column := chi.URLParam(r, "column")

		data, err := svc.GetCategoryTile(z, x, y, column)
		if err != nil {
			data, _ = svc.GetEmptyTile()
		}

		w.Header().Set("Content-Type", "image/png")
		w.Header().Set("Cache-Control", "public, max-age=3600")
		w.Write(data)
	}
}

func metadataHandler(svc *service.TileService) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		metadata := svc.Metadata()
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(metadata)
	}
}

func genesHandler(svc *service.TileService) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		metadata := svc.Metadata()
		response := map[string]interface{}{
			"genes":              metadata.Genes,
			"total":              len(metadata.Genes),
			"preaggregated_count": metadata.NGenes,
		}
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(response)
	}
}

func geneInfoHandler(svc *service.TileService) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		gene := chi.URLParam(r, "gene")
		metadata := svc.Metadata()

		idx, ok := metadata.GeneIndex[gene]
		if !ok {
			http.Error(w, "gene not found", http.StatusNotFound)
			return
		}

		response := map[string]interface{}{
			"name":          gene,
			"index":         idx,
			"preaggregated": true,
		}
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(response)
	}
}

func geneBinsHandler(svc *service.TileService) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		gene := chi.URLParam(r, "gene")

		// Parse query parameters
		threshold, _ := strconv.ParseFloat(r.URL.Query().Get("threshold"), 32)
		offset, _ := strconv.Atoi(r.URL.Query().Get("offset"))
		limit, _ := strconv.Atoi(r.URL.Query().Get("limit"))
		if limit <= 0 || limit > 1000 {
			limit = 100
		}

		result, err := svc.QueryBinsExpressingGene(gene, float32(threshold), offset, limit)
		if err != nil {
			http.Error(w, err.Error(), http.StatusInternalServerError)
			return
		}

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(result)
	}
}

func geneStatsHandler(svc *service.TileService) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		gene := chi.URLParam(r, "gene")

		stats, err := svc.GetGeneStats(gene)
		if err != nil {
			http.Error(w, err.Error(), http.StatusNotFound)
			return
		}

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(stats)
	}
}

func categoriesHandler(svc *service.TileService) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		metadata := svc.Metadata()
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(metadata.Categories)
	}
}

func categoryColorsHandler(svc *service.TileService) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		column := chi.URLParam(r, "column")

		colors, err := svc.GetCategoryColors(column)
		if err != nil {
			http.Error(w, err.Error(), http.StatusNotFound)
			return
		}

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(colors)
	}
}

func categoryLegendHandler(svc *service.TileService) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		column := chi.URLParam(r, "column")

		legend, err := svc.GetCategoryLegend(column)
		if err != nil {
			http.Error(w, err.Error(), http.StatusNotFound)
			return
		}

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(legend)
	}
}

func statsHandler(svc *service.TileService) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		metadata := svc.Metadata()
		response := map[string]interface{}{
			"n_cells":       metadata.NCells,
			"n_genes":       metadata.NGenes,
			"zoom_levels":   metadata.ZoomLevels,
			"dataset_name":  metadata.DatasetName,
		}
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(response)
	}
}
