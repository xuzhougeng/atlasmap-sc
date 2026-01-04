// Package api provides HTTP handlers for the SOMA-Tiles server.
package api

import (
	"bytes"
	"context"
	"encoding/json"
	"errors"
	"io"
	"net/http"
	"net/url"
	"strconv"
	"strings"

	"github.com/go-chi/chi/v5"
	"github.com/go-chi/chi/v5/middleware"
	"github.com/go-chi/cors"
	"github.com/soma-tiles/server/internal/service"
)

// RouterConfig contains router configuration.
type RouterConfig struct {
	Registry    *DatasetRegistry
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

	// Global datasets endpoint (not dataset-scoped)
	r.Get("/api/datasets", datasetsHandler(cfg.Registry))

	// Default dataset routes (legacy compatibility)
	defaultSvc := cfg.Registry.Default()
	if defaultSvc != nil {
		mountDatasetRoutes(r, defaultSvc)
	}

	// Dataset-scoped routes: /d/{dataset}/...
	r.Route("/d/{dataset}", func(r chi.Router) {
		r.Use(datasetMiddleware(cfg.Registry))

		// Tile endpoints
		r.Get("/tiles/{z}/{x}/{y}.png", datasetTileHandler)
		r.Get("/tiles/{z}/{x}/{y}/expression/{gene}.png", datasetExpressionTileHandler)
		r.Get("/tiles/{z}/{x}/{y}/category/{column}.png", datasetCategoryTileHandler)
		r.Post("/tiles/{z}/{x}/{y}/category/{column}.png", datasetCategoryTileHandler)

		// API endpoints
		r.Route("/api", func(r chi.Router) {
			r.Get("/metadata", datasetMetadataHandler)
			r.Get("/genes", datasetGenesHandler)
			r.Get("/genes/{gene}", datasetGeneInfoHandler)
			r.Get("/genes/{gene}/bins", datasetGeneBinsHandler)
			r.Get("/genes/{gene}/stats", datasetGeneStatsHandler)
			r.Get("/genes/{gene}/category/{column}/means", datasetGeneCategoryMeansHandler)
			r.Get("/categories", datasetCategoriesHandler)
			r.Get("/categories/{column}/colors", datasetCategoryColorsHandler)
			r.Get("/categories/{column}/legend", datasetCategoryLegendHandler)
			r.Get("/stats", datasetStatsHandler)
		})
	})

	return r
}

// mountDatasetRoutes mounts routes for a specific tile service (used for default dataset).
func mountDatasetRoutes(r chi.Router, svc *service.TileService) {
	// Tile endpoints
	r.Get("/tiles/{z}/{x}/{y}.png", tileHandler(svc))
	r.Get("/tiles/{z}/{x}/{y}/expression/{gene}.png", expressionTileHandler(svc))
	r.Get("/tiles/{z}/{x}/{y}/category/{column}.png", categoryTileHandler(svc))
	r.Post("/tiles/{z}/{x}/{y}/category/{column}.png", categoryTileHandler(svc))

	// API endpoints
	r.Route("/api", func(r chi.Router) {
		r.Get("/metadata", metadataHandler(svc))
		r.Get("/genes", genesHandler(svc))
		r.Get("/genes/{gene}", geneInfoHandler(svc))
		r.Get("/genes/{gene}/bins", geneBinsHandler(svc))
		r.Get("/genes/{gene}/stats", geneStatsHandler(svc))
		r.Get("/genes/{gene}/category/{column}/means", geneCategoryMeansHandler(svc))
		r.Get("/categories", categoriesHandler(svc))
		r.Get("/categories/{column}/colors", categoryColorsHandler(svc))
		r.Get("/categories/{column}/legend", categoryLegendHandler(svc))
		r.Get("/stats", statsHandler(svc))
	})
}

// Context key for dataset service
type ctxKey string

const datasetServiceKey ctxKey = "datasetService"

// datasetMiddleware resolves the dataset from URL and injects the tile service into context.
func datasetMiddleware(registry *DatasetRegistry) func(http.Handler) http.Handler {
	return func(next http.Handler) http.Handler {
		return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
			datasetID := chi.URLParam(r, "dataset")
			svc := registry.Get(datasetID)
			if svc == nil {
				http.Error(w, "dataset not found: "+datasetID, http.StatusNotFound)
				return
			}
			ctx := context.WithValue(r.Context(), datasetServiceKey, svc)
			next.ServeHTTP(w, r.WithContext(ctx))
		})
	}
}

func getDatasetService(r *http.Request) *service.TileService {
	if svc, ok := r.Context().Value(datasetServiceKey).(*service.TileService); ok {
		return svc
	}
	return nil
}

// datasetsHandler returns the list of available datasets.
func datasetsHandler(registry *DatasetRegistry) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		response := map[string]interface{}{
			"default":  registry.DefaultDatasetID(),
			"datasets": registry.Datasets(),
		}
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(response)
	}
}

// Dataset-scoped handlers (get service from context)
func datasetTileHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not found", http.StatusInternalServerError)
		return
	}
	tileHandler(svc)(w, r)
}

func datasetExpressionTileHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not found", http.StatusInternalServerError)
		return
	}
	expressionTileHandler(svc)(w, r)
}

func datasetCategoryTileHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not found", http.StatusInternalServerError)
		return
	}
	categoryTileHandler(svc)(w, r)
}

func datasetMetadataHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not found", http.StatusInternalServerError)
		return
	}
	metadataHandler(svc)(w, r)
}

func datasetGenesHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not found", http.StatusInternalServerError)
		return
	}
	genesHandler(svc)(w, r)
}

func datasetGeneInfoHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not found", http.StatusInternalServerError)
		return
	}
	geneInfoHandler(svc)(w, r)
}

func datasetGeneBinsHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not found", http.StatusInternalServerError)
		return
	}
	geneBinsHandler(svc)(w, r)
}

func datasetGeneStatsHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not found", http.StatusInternalServerError)
		return
	}
	geneStatsHandler(svc)(w, r)
}

func datasetGeneCategoryMeansHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not found", http.StatusInternalServerError)
		return
	}
	geneCategoryMeansHandler(svc)(w, r)
}

func datasetCategoriesHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not found", http.StatusInternalServerError)
		return
	}
	categoriesHandler(svc)(w, r)
}

func datasetCategoryColorsHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not found", http.StatusInternalServerError)
		return
	}
	categoryColorsHandler(svc)(w, r)
}

func datasetCategoryLegendHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not found", http.StatusInternalServerError)
		return
	}
	categoryLegendHandler(svc)(w, r)
}

func datasetStatsHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not found", http.StatusInternalServerError)
		return
	}
	statsHandler(svc)(w, r)
}

// Original handlers (take service as parameter)
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

		// Parse optional category filter
		var (
			categoryFilter []string
			hasFilter      bool
			err            error
		)
		if r.Method == http.MethodPost {
			categoryFilter, hasFilter, err = parseCategoryFilterBody(r)
			if err != nil {
				http.Error(w, err.Error(), http.StatusBadRequest)
				return
			}
		} else {
			categoryFilter, hasFilter = parseCategoryFilter(r.URL.Query())
		}
		if !hasFilter {
			categoryFilter = nil
		}

		data, err := svc.GetCategoryTile(z, x, y, column, categoryFilter)
		if err != nil {
			data, _ = svc.GetEmptyTile()
		}

		w.Header().Set("Content-Type", "image/png")
		w.Header().Set("Cache-Control", "public, max-age=3600")
		w.Write(data)
	}
}

func parseCategoryFilter(query url.Values) ([]string, bool) {
	rawValues, present := query["categories"]
	if !present {
		return nil, false
	}

	// Support repeated query parameters:
	//   ?categories=T&categories=B
	if len(rawValues) > 1 {
		out := make([]string, 0, len(rawValues))
		for _, v := range rawValues {
			v = strings.TrimSpace(v)
			if v != "" {
				out = append(out, v)
			}
		}
		return out, true
	}

	raw := strings.TrimSpace(rawValues[0])
	if raw == "" {
		// Explicit "filter to none".
		return make([]string, 0), true
	}

	// Preferred format (frontend): JSON array, e.g. ["T","B"] (allows commas in values).
	if strings.HasPrefix(raw, "[") {
		var categories []string
		if err := json.Unmarshal([]byte(raw), &categories); err == nil {
			if categories == nil {
				return make([]string, 0), true
			}
			return categories, true
		}
		// Fall through to comma-separated parsing for tolerance.
	}

	// Legacy format: comma-separated list, e.g. T,B
	parts := strings.Split(raw, ",")
	out := make([]string, 0, len(parts))
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p != "" {
			out = append(out, p)
		}
	}
	return out, true
}

const maxCategoryFilterBodyBytes = 10 << 20 // 10 MiB

func parseCategoryFilterBody(r *http.Request) ([]string, bool, error) {
	if r.Body == nil {
		return nil, false, nil
	}

	body, err := io.ReadAll(io.LimitReader(r.Body, maxCategoryFilterBodyBytes+1))
	if err != nil {
		return nil, false, err
	}
	if len(body) > maxCategoryFilterBodyBytes {
		return nil, false, errors.New("category filter body too large")
	}

	raw := bytes.TrimSpace(body)
	if len(raw) == 0 {
		return nil, false, nil
	}

	// Preferred POST payload: a JSON array, e.g. ["T","B"].
	// Also supports an object payload: {"categories":[...]}.
	if len(raw) > 0 && raw[0] == '{' {
		var payload map[string]json.RawMessage
		if err := json.Unmarshal(raw, &payload); err == nil {
			rawCategories, ok := payload["categories"]
			if !ok {
				return nil, false, nil
			}

			rawCategories = bytes.TrimSpace(rawCategories)
			if len(rawCategories) == 0 || bytes.Equal(rawCategories, []byte("null")) {
				return nil, false, nil
			}

			var categories []string
			if err := json.Unmarshal(rawCategories, &categories); err == nil {
				if categories == nil {
					return make([]string, 0), true, nil
				}
				return categories, true, nil
			}

			var categoriesString string
			if err := json.Unmarshal(rawCategories, &categoriesString); err == nil {
				filter, ok := parseCategoryFilter(url.Values{"categories": {categoriesString}})
				return filter, ok, nil
			}
		}
		// Fall through to tolerant parsing of the raw body.
	}

	// Support form-encoded bodies:
	//   categories=T&categories=B
	//   categories=["T","B"]
	if bytes.Contains(raw, []byte("=")) && raw[0] != '[' {
		if q, err := url.ParseQuery(string(raw)); err == nil {
			filter, ok := parseCategoryFilter(q)
			return filter, ok, nil
		}
	}

	filter, ok := parseCategoryFilter(url.Values{"categories": {string(raw)}})
	return filter, ok, nil
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
			"genes":               metadata.Genes,
			"total":               len(metadata.Genes),
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

		// Parse zoom parameter (default to 0 for backward compatibility)
		zoom := 0
		if zoomStr := r.URL.Query().Get("zoom"); zoomStr != "" {
			z, err := strconv.Atoi(zoomStr)
			if err != nil {
				http.Error(w, "invalid zoom parameter", http.StatusBadRequest)
				return
			}
			zoom = z
		}

		stats, err := svc.GetGeneStats(gene, zoom)
		if err != nil {
			status := http.StatusNotFound
			if strings.Contains(err.Error(), "invalid zoom") {
				status = http.StatusBadRequest
			}
			http.Error(w, err.Error(), status)
			return
		}

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(stats)
	}
}

func geneCategoryMeansHandler(svc *service.TileService) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		gene := chi.URLParam(r, "gene")
		column := chi.URLParam(r, "column")

		items, err := svc.GetGeneCategoryMeans(gene, column)
		if err != nil {
			status := http.StatusInternalServerError
			if strings.Contains(err.Error(), "not found") {
				status = http.StatusNotFound
			}
			http.Error(w, err.Error(), status)
			return
		}

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(map[string]interface{}{
			"gene":   gene,
			"column": column,
			"items":  items,
		})
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
			"n_cells":      metadata.NCells,
			"n_genes":      metadata.NGenes,
			"zoom_levels":  metadata.ZoomLevels,
			"dataset_name": metadata.DatasetName,
		}
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(response)
	}
}
