// Package api provides HTTP handlers for the SOMA-Tiles server.
package api

import (
	"bytes"
	"context"
	"encoding/json"
	"errors"
	"io"
	"math"
	"mime"
	"net/http"
	"net/url"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"github.com/go-chi/chi/v5"
	"github.com/go-chi/chi/v5/middleware"
	"github.com/go-chi/cors"
	"github.com/soma-tiles/server/internal/blaststore"
	"github.com/soma-tiles/server/internal/data/soma"
	"github.com/soma-tiles/server/internal/destore"
	"github.com/soma-tiles/server/internal/service"
)

// RouterConfig contains router configuration.
type RouterConfig struct {
	Registry        *DatasetRegistry
	CORSOrigins     []string
	JobManager      *JobManager
	BlastJobManager *BlastJobManager
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

	// Global gene_lookup endpoint (resolves gene_id -> matching datasets)
	r.Get("/api/gene_lookup", geneLookupHandler(cfg.Registry))

	// Global BLASTP job endpoints (not dataset-scoped)
	r.Route("/api/blastp/jobs", func(r chi.Router) {
		r.Post("/", blastJobSubmitHandler(cfg.BlastJobManager))
		r.Get("/{job_id}", blastJobStatusHandler(cfg.BlastJobManager))
		r.Get("/{job_id}/result", blastJobResultHandler(cfg.BlastJobManager))
		r.Delete("/{job_id}", blastJobCancelHandler(cfg.BlastJobManager))
	})

	// Dataset-scoped routes: /d/{dataset}/...
	r.Route("/d/{dataset}", func(r chi.Router) {
		r.Use(datasetMiddleware(cfg.Registry))

		// Tile endpoints
		r.Get("/tiles/{z}/{x}/{y}.png", datasetTileHandler)
		r.Get("/tiles/{z}/{x}/{y}/expression/{gene}.png", datasetExpressionTileHandler)
		// NOTE: chi treats '.' as a param delimiter when the route pattern is `{gene}.png`,
		// which breaks genes/columns containing '.' (e.g. "Azfi-s0217.g058558").
		// Add a fallback route that captures the full segment (including ".png") and strip
		// the extension in the handler.
		r.Get("/tiles/{z}/{x}/{y}/expression/{gene}", datasetExpressionTileHandler)
		r.Get("/tiles/{z}/{x}/{y}/category/{column}.png", datasetCategoryTileHandler)
		r.Post("/tiles/{z}/{x}/{y}/category/{column}.png", datasetCategoryTileHandler)
		r.Get("/tiles/{z}/{x}/{y}/category/{column}", datasetCategoryTileHandler)
		r.Post("/tiles/{z}/{x}/{y}/category/{column}", datasetCategoryTileHandler)

		// API endpoints
		r.Route("/api", func(r chi.Router) {
			r.Get("/metadata", datasetMetadataHandler)
			r.Get("/h5ad", datasetH5ADDownloadHandler(cfg.Registry))
			r.Get("/genes", datasetGenesHandler)
			r.Get("/genes/{gene}", datasetGeneInfoHandler)
			r.Get("/genes/{gene}/bins", datasetGeneBinsHandler)
			r.Get("/genes/{gene}/stats", datasetGeneStatsHandler)
			r.Get("/genes/{gene}/category/{column}/means", datasetGeneCategoryMeansHandler)
			r.Get("/categories", datasetCategoriesHandler)
			r.Get("/categories/{column}/colors", datasetCategoryColorsHandler)
			r.Get("/categories/{column}/centroids", datasetCategoryCentroidsHandler)
			r.Get("/categories/{column}/legend", datasetCategoryLegendHandler)
			r.Get("/stats", datasetStatsHandler)
			r.Get("/soma/expression", datasetSomaExpressionHandler)

			// SOMA DE job endpoints
			r.Route("/soma/de/jobs", func(r chi.Router) {
				r.Post("/", deJobSubmitHandler(cfg.JobManager))
				r.Get("/{job_id}", deJobStatusHandler(cfg.JobManager))
				r.Get("/{job_id}/result", deJobResultHandler(cfg.JobManager))
				r.Delete("/{job_id}", deJobCancelHandler(cfg.JobManager))
			})

			// SOMA obs metadata endpoints (for groupby column discovery)
			r.Get("/soma/obs/columns", datasetSomaObsColumnsHandler)
			r.Get("/soma/obs/{column}/values", datasetSomaObsValuesHandler)
		})
	})

	return r
}

// Context key for dataset service
type ctxKey string

const datasetServiceKey ctxKey = "datasetService"

// datasetMiddleware resolves the dataset from URL and injects the tile service into context.
func datasetMiddleware(registry *DatasetRegistry) func(http.Handler) http.Handler {
	return func(next http.Handler) http.Handler {
		return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
			datasetID := chi.URLParam(r, "dataset")
			coordKey := strings.TrimSpace(r.URL.Query().Get("coord"))
			svc := registry.GetCoord(datasetID, coordKey)
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
			"title":    registry.Title(),
		}
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(response)
	}
}

// geneLookupHandler resolves a gene_id to the list of datasets containing it.
func geneLookupHandler(registry *DatasetRegistry) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		geneID := strings.TrimSpace(r.URL.Query().Get("gene_id"))
		if geneID == "" {
			http.Error(w, "missing required query param: gene_id", http.StatusBadRequest)
			return
		}

		var matchingDatasets []string
		for _, dsID := range registry.DatasetIDs() {
			svc := registry.Get(dsID)
			if svc == nil {
				continue
			}
			md := svc.Metadata()
			if md.GeneIndex == nil {
				continue
			}
			if _, ok := md.GeneIndex[geneID]; ok {
				matchingDatasets = append(matchingDatasets, dsID)
			}
		}

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(map[string]interface{}{
			"gene_id":  geneID,
			"datasets": matchingDatasets,
		})
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

func datasetH5ADDownloadHandler(registry *DatasetRegistry) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		datasetID := chi.URLParam(r, "dataset")
		path := registry.H5ADPath(datasetID)
		if path == "" {
			http.NotFound(w, r)
			return
		}

		info, err := os.Stat(path)
		if err != nil || info.IsDir() {
			http.NotFound(w, r)
			return
		}

		filename := filepath.Base(path)
		disposition := mime.FormatMediaType("attachment", map[string]string{"filename": filename})
		if disposition != "" {
			w.Header().Set("Content-Disposition", disposition)
		} else {
			w.Header().Set("Content-Disposition", "attachment")
		}
		w.Header().Set("Content-Type", "application/octet-stream")

		http.ServeFile(w, r, path)
	}
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

func datasetCategoryCentroidsHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not found", http.StatusInternalServerError)
		return
	}
	categoryCentroidsHandler(svc)(w, r)
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
		gene = strings.TrimSuffix(gene, ".png")
		colormap := r.URL.Query().Get("colormap")
		if colormap == "" {
			colormap = "viridis"
		}
		pointSize := parsePointSize(r.URL.Query())

		var (
			minPtr *float32
			maxPtr *float32
		)
		if minStr := strings.TrimSpace(r.URL.Query().Get("min")); minStr != "" {
			if minF64, err := strconv.ParseFloat(minStr, 32); err == nil {
				if !math.IsNaN(minF64) && !math.IsInf(minF64, 0) {
					minF32 := float32(minF64)
					minPtr = &minF32
				}
			}
		}
		if maxStr := strings.TrimSpace(r.URL.Query().Get("max")); maxStr != "" {
			if maxF64, err := strconv.ParseFloat(maxStr, 32); err == nil {
				if !math.IsNaN(maxF64) && !math.IsInf(maxF64, 0) {
					maxF32 := float32(maxF64)
					maxPtr = &maxF32
				}
			}
		}

		data, err := svc.GetExpressionTile(z, x, y, gene, colormap, minPtr, maxPtr, pointSize)
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
		column = strings.TrimSuffix(column, ".png")
		pointSize := parsePointSize(r.URL.Query())

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

		data, err := svc.GetCategoryTile(z, x, y, column, categoryFilter, pointSize)
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

func parsePointSize(query url.Values) float64 {
	const defaultPointSize = 1.0
	raw := strings.TrimSpace(query.Get("point_size"))
	if raw == "" {
		return defaultPointSize
	}
	v, err := strconv.ParseFloat(raw, 64)
	if err != nil || math.IsNaN(v) || math.IsInf(v, 0) {
		return defaultPointSize
	}
	// Clamp to a sane range. 1.0 means "fill the full bin".
	if v < 0.1 {
		v = 0.1
	}
	if v > 5.0 {
		v = 5.0
	}
	// Quantize for stable caching.
	return math.Round(v*1000) / 1000
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

func categoryCentroidsHandler(svc *service.TileService) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		column := chi.URLParam(r, "column")

		centroids, err := svc.GetCategoryCentroids(column)
		if err != nil {
			status := http.StatusInternalServerError
			if strings.Contains(err.Error(), "not found") {
				status = http.StatusNotFound
			}
			http.Error(w, err.Error(), status)
			return
		}

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(centroids)
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

type somaExprItem struct {
	CellJoinID int64   `json:"cell_joinid"`
	Value      float32 `json:"value"`
}

func datasetSomaExpressionHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not available", http.StatusInternalServerError)
		return
	}
	sr := svc.Soma()
	if sr == nil {
		http.Error(w, "soma not configured for this dataset", http.StatusNotFound)
		return
	}
	if !sr.Supported() {
		http.Error(w, soma.ErrUnsupported.Error(), http.StatusNotImplemented)
		return
	}

	gene := strings.TrimSpace(r.URL.Query().Get("gene"))
	if gene == "" {
		http.Error(w, "missing required query param: gene", http.StatusBadRequest)
		return
	}
	cellsParam := strings.TrimSpace(r.URL.Query().Get("cells"))
	if cellsParam == "" {
		http.Error(w, "missing required query param: cells (comma-separated cell soma_joinid list)", http.StatusBadRequest)
		return
	}
	cells, err := parseCSVInt64(cellsParam, 10000)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}
	mode := strings.TrimSpace(r.URL.Query().Get("mode"))
	if mode == "" {
		mode = "sparse"
	}
	if mode != "sparse" && mode != "dense" {
		http.Error(w, "invalid mode (expected sparse or dense)", http.StatusBadRequest)
		return
	}

	geneJoinID, vals, err := sr.ExpressionByCellJoinID(gene, cells)
	if err != nil {
		// Most common errors: gene not found, TileDB not available, dataset missing.
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}

	w.Header().Set("Content-Type", "application/json")

	if mode == "dense" {
		out := make([]float32, len(cells))
		nnz := 0
		for i, cid := range cells {
			if v, ok := vals[cid]; ok {
				out[i] = v
				nnz++
			}
		}
		json.NewEncoder(w).Encode(map[string]interface{}{
			"gene":        gene,
			"gene_joinid": geneJoinID,
			"mode":        "dense",
			"cells":       cells,
			"values":      out,
			"nnz":         nnz,
		})
		return
	}

	// sparse: return only non-zeros
	items := make([]somaExprItem, 0, len(vals))
	for cid, v := range vals {
		items = append(items, somaExprItem{CellJoinID: cid, Value: v})
	}
	sort.Slice(items, func(i, j int) bool { return items[i].CellJoinID < items[j].CellJoinID })
	json.NewEncoder(w).Encode(map[string]interface{}{
		"gene":        gene,
		"gene_joinid": geneJoinID,
		"mode":        "sparse",
		"items":       items,
		"nnz":         len(items),
	})
}

func parseCSVInt64(s string, maxItems int) ([]int64, error) {
	parts := strings.Split(s, ",")
	out := make([]int64, 0, len(parts))
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p == "" {
			continue
		}
		if maxItems > 0 && len(out) >= maxItems {
			return nil, errors.New("too many cells (increase limit or batch requests)")
		}
		v, err := strconv.ParseInt(p, 10, 64)
		if err != nil {
			return nil, errors.New("invalid cell joinid: " + p)
		}
		out = append(out, v)
	}
	if len(out) == 0 {
		return nil, errors.New("empty cells list")
	}
	return out, nil
}

// DE Job handlers

type deJobSubmitRequest struct {
	Groupby          string   `json:"groupby"`
	Group1           []string `json:"group1"`
	Group2           []string `json:"group2"`
	Tests            []string `json:"tests"`
	MaxCellsPerGroup int      `json:"max_cells_per_group"`
	Seed             int      `json:"seed"`
	Limit            int      `json:"limit"`
}

func deJobSubmitHandler(jm *JobManager) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		if jm == nil {
			http.Error(w, "job manager not configured", http.StatusNotImplemented)
			return
		}

		svc := getDatasetService(r)
		if svc == nil {
			http.Error(w, "dataset service not available", http.StatusInternalServerError)
			return
		}
		sr := svc.Soma()
		if sr == nil {
			http.Error(w, "soma not configured for this dataset", http.StatusNotFound)
			return
		}
		if !sr.Supported() {
			http.Error(w, soma.ErrUnsupported.Error(), http.StatusNotImplemented)
			return
		}

		var req deJobSubmitRequest
		if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
			http.Error(w, "invalid request body: "+err.Error(), http.StatusBadRequest)
			return
		}

		// Validate required fields
		if req.Groupby == "" {
			http.Error(w, "groupby is required", http.StatusBadRequest)
			return
		}
		if len(req.Group1) == 0 {
			http.Error(w, "group1 is required (at least one category value)", http.StatusBadRequest)
			return
		}

		// Apply defaults
		if len(req.Tests) == 0 {
			req.Tests = []string{"ttest", "ranksum"}
		}
		if req.MaxCellsPerGroup <= 0 {
			req.MaxCellsPerGroup = 2000
		}
		if req.MaxCellsPerGroup > 20000 {
			req.MaxCellsPerGroup = 20000
		}
		if req.Limit <= 0 {
			req.Limit = 50
		}
		if req.Limit > 500 {
			req.Limit = 500
		}

		datasetID := chi.URLParam(r, "dataset")
		params := destore.DEJobParams{
			DatasetID:        datasetID,
			Groupby:          req.Groupby,
			Group1:           req.Group1,
			Group2:           req.Group2,
			Tests:            req.Tests,
			MaxCellsPerGroup: req.MaxCellsPerGroup,
			Seed:             req.Seed,
			Limit:            req.Limit,
		}

		job, err := jm.Submit(params)
		if err != nil {
			http.Error(w, "failed to submit job: "+err.Error(), http.StatusInternalServerError)
			return
		}

		w.Header().Set("Content-Type", "application/json")
		w.WriteHeader(http.StatusAccepted)
		json.NewEncoder(w).Encode(map[string]interface{}{
			"job_id": job.ID,
			"status": job.Status,
		})
	}
}

func deJobStatusHandler(jm *JobManager) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		if jm == nil {
			http.Error(w, "job manager not configured", http.StatusNotImplemented)
			return
		}

		jobID := chi.URLParam(r, "job_id")
		job := jm.Get(jobID)
		if job == nil {
			http.Error(w, "job not found", http.StatusNotFound)
			return
		}

		// Check dataset matches
		datasetID := chi.URLParam(r, "dataset")
		if job.Params.DatasetID != datasetID {
			http.Error(w, "job not found", http.StatusNotFound)
			return
		}

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(map[string]interface{}{
			"job_id":      job.ID,
			"status":      job.Status,
			"created_at":  job.CreatedAt,
			"started_at":  job.StartedAt,
			"finished_at": job.FinishedAt,
			"progress":    job.Progress,
			"n1":          job.N1,
			"n2":          job.N2,
			"error":       job.Error,
		})
	}
}

func deJobResultHandler(jm *JobManager) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		if jm == nil {
			http.Error(w, "job manager not configured", http.StatusNotImplemented)
			return
		}

		jobID := chi.URLParam(r, "job_id")
		job := jm.Get(jobID)
		if job == nil {
			http.Error(w, "job not found", http.StatusNotFound)
			return
		}

		datasetID := chi.URLParam(r, "dataset")
		if job.Params.DatasetID != datasetID {
			http.Error(w, "job not found", http.StatusNotFound)
			return
		}

		if job.Status != destore.JobStatusCompleted {
			http.Error(w, "job not completed (status: "+string(job.Status)+")", http.StatusBadRequest)
			return
		}

		// Parse pagination and order params
		offset := 0
		limit := job.Params.Limit
		if limit <= 0 {
			limit = 50
		}
		orderBy := r.URL.Query().Get("order_by")
		if orderBy == "" {
			orderBy = "fdr_ranksum"
		}
		if offsetStr := r.URL.Query().Get("offset"); offsetStr != "" {
			if v, err := strconv.Atoi(offsetStr); err == nil && v >= 0 {
				offset = v
			}
		}
		if limitStr := r.URL.Query().Get("limit"); limitStr != "" {
			if v, err := strconv.Atoi(limitStr); err == nil && v > 0 {
				limit = v
				if limit > 500 {
					limit = 500
				}
			}
		}

		// Query results from SQLite
		items, total, err := jm.Store().QueryResults(jobID, orderBy, offset, limit)
		if err != nil {
			http.Error(w, "failed to query results: "+err.Error(), http.StatusInternalServerError)
			return
		}

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(map[string]interface{}{
			"params":   job.Params,
			"n1":       job.N1,
			"n2":       job.N2,
			"total":    total,
			"offset":   offset,
			"limit":    limit,
			"order_by": orderBy,
			"items":    items,
		})
	}
}

func deJobCancelHandler(jm *JobManager) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		if jm == nil {
			http.Error(w, "job manager not configured", http.StatusNotImplemented)
			return
		}

		jobID := chi.URLParam(r, "job_id")
		job := jm.Get(jobID)
		if job == nil {
			http.Error(w, "job not found", http.StatusNotFound)
			return
		}

		datasetID := chi.URLParam(r, "dataset")
		if job.Params.DatasetID != datasetID {
			http.Error(w, "job not found", http.StatusNotFound)
			return
		}

		jm.Cancel(jobID)

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(map[string]interface{}{
			"job_id":    jobID,
			"cancelled": true,
		})
	}
}

// SOMA obs metadata handlers

func datasetSomaObsColumnsHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not available", http.StatusInternalServerError)
		return
	}
	sr := svc.Soma()
	if sr == nil {
		http.Error(w, "soma not configured for this dataset", http.StatusNotFound)
		return
	}
	if !sr.Supported() {
		http.Error(w, soma.ErrUnsupported.Error(), http.StatusNotImplemented)
		return
	}

	columns, err := sr.ObsColumns()
	if err != nil {
		http.Error(w, "failed to get obs columns: "+err.Error(), http.StatusInternalServerError)
		return
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"columns": columns,
	})
}

func datasetSomaObsValuesHandler(w http.ResponseWriter, r *http.Request) {
	svc := getDatasetService(r)
	if svc == nil {
		http.Error(w, "dataset service not available", http.StatusInternalServerError)
		return
	}
	sr := svc.Soma()
	if sr == nil {
		http.Error(w, "soma not configured for this dataset", http.StatusNotFound)
		return
	}
	if !sr.Supported() {
		http.Error(w, soma.ErrUnsupported.Error(), http.StatusNotImplemented)
		return
	}

	column := chi.URLParam(r, "column")
	values, err := sr.ObsColumnValues(column)
	if err != nil {
		http.Error(w, "failed to get column values: "+err.Error(), http.StatusBadRequest)
		return
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"column": column,
		"values": values,
	})
}

// BLAST Job handlers

type blastJobSubmitRequest struct {
	Sequence   string   `json:"sequence"`
	MaxHits    int      `json:"max_hits"`
	Evalue     float64  `json:"evalue"`
	Datasets   []string `json:"datasets"`
	NumThreads int      `json:"num_threads"`
}

func blastJobSubmitHandler(jm *BlastJobManager) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		if jm == nil {
			http.Error(w, "blast job manager not configured", http.StatusNotImplemented)
			return
		}

		var req blastJobSubmitRequest
		if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
			http.Error(w, "invalid request body: "+err.Error(), http.StatusBadRequest)
			return
		}

		// Validate required fields
		if strings.TrimSpace(req.Sequence) == "" {
			http.Error(w, "sequence is required", http.StatusBadRequest)
			return
		}

		// Apply defaults
		if req.MaxHits <= 0 {
			req.MaxHits = 10
		}
		if req.MaxHits > 100 {
			req.MaxHits = 100
		}
		if req.Evalue <= 0 {
			req.Evalue = 1e-5
		}
		if req.NumThreads <= 0 {
			req.NumThreads = 1
		}
		if req.NumThreads > 4 {
			req.NumThreads = 4
		}

		params := blaststore.BlastJobParams{
			Sequence:   req.Sequence,
			MaxHits:    req.MaxHits,
			Evalue:     req.Evalue,
			Datasets:   req.Datasets,
			NumThreads: req.NumThreads,
		}

		job, err := jm.Submit(params)
		if err != nil {
			http.Error(w, "failed to submit job: "+err.Error(), http.StatusInternalServerError)
			return
		}

		w.Header().Set("Content-Type", "application/json")
		w.WriteHeader(http.StatusAccepted)
		json.NewEncoder(w).Encode(map[string]interface{}{
			"job_id": job.ID,
			"status": job.Status,
		})
	}
}

func blastJobStatusHandler(jm *BlastJobManager) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		if jm == nil {
			http.Error(w, "blast job manager not configured", http.StatusNotImplemented)
			return
		}

		jobID := chi.URLParam(r, "job_id")
		job := jm.Get(jobID)
		if job == nil {
			http.Error(w, "job not found", http.StatusNotFound)
			return
		}

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(map[string]interface{}{
			"job_id":      job.ID,
			"status":      job.Status,
			"created_at":  job.CreatedAt,
			"started_at":  job.StartedAt,
			"finished_at": job.FinishedAt,
			"progress":    job.Progress,
			"error":       job.Error,
		})
	}
}

func blastJobResultHandler(jm *BlastJobManager) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		if jm == nil {
			http.Error(w, "blast job manager not configured", http.StatusNotImplemented)
			return
		}

		jobID := chi.URLParam(r, "job_id")
		job := jm.Get(jobID)
		if job == nil {
			http.Error(w, "job not found", http.StatusNotFound)
			return
		}

		if job.Status != blaststore.JobStatusCompleted {
			http.Error(w, "job not completed (status: "+string(job.Status)+")", http.StatusBadRequest)
			return
		}

		// Parse pagination and order params
		offset := 0
		limit := 100
		orderBy := r.URL.Query().Get("order_by")
		if orderBy == "" {
			orderBy = "bitscore"
		}
		if offsetStr := r.URL.Query().Get("offset"); offsetStr != "" {
			if v, err := strconv.Atoi(offsetStr); err == nil && v >= 0 {
				offset = v
			}
		}
		if limitStr := r.URL.Query().Get("limit"); limitStr != "" {
			if v, err := strconv.Atoi(limitStr); err == nil && v > 0 {
				limit = v
				if limit > 500 {
					limit = 500
				}
			}
		}

		// Query results from SQLite
		items, total, err := jm.Store().QueryResults(jobID, orderBy, offset, limit)
		if err != nil {
			http.Error(w, "failed to query results: "+err.Error(), http.StatusInternalServerError)
			return
		}

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(map[string]interface{}{
			"params":   job.Params,
			"total":    total,
			"offset":   offset,
			"limit":    limit,
			"order_by": orderBy,
			"items":    items,
		})
	}
}

func blastJobCancelHandler(jm *BlastJobManager) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		if jm == nil {
			http.Error(w, "blast job manager not configured", http.StatusNotImplemented)
			return
		}

		jobID := chi.URLParam(r, "job_id")
		job := jm.Get(jobID)
		if job == nil {
			http.Error(w, "job not found", http.StatusNotFound)
			return
		}

		jm.Cancel(jobID)

		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(map[string]interface{}{
			"job_id":    jobID,
			"cancelled": true,
		})
	}
}
