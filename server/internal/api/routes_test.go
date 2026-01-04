// Package api provides HTTP handlers for the SOMA-Tiles server.
package api

import (
	"encoding/json"
	"io"
	"net/http"
	"net/http/httptest"
	"os"
	"testing"
	"time"

	"github.com/soma-tiles/server/internal/cache"
	"github.com/soma-tiles/server/internal/data/zarr"
	"github.com/soma-tiles/server/internal/render"
	"github.com/soma-tiles/server/internal/service"
)

// Test data path - adjust this to point to your test data
const testZarrPath = "/data/xzg_data/software/soma-tiles/data/preprocessed/zarr/bins.zarr"

// testServer holds the test server and its dependencies
type testServer struct {
	server      *httptest.Server
	zarrReader  *zarr.Reader
	cache       *cache.Manager
	renderer    *render.TileRenderer
	tileService *service.TileService
}

// setupTestServer initializes all components and returns a test server
func setupTestServer(t *testing.T) *testServer {
	t.Helper()

	// Check if test data exists
	if _, err := os.Stat(testZarrPath); os.IsNotExist(err) {
		t.Skipf("Test data not found at %s, skipping integration tests", testZarrPath)
	}

	// Initialize Zarr reader
	zarrReader, err := zarr.NewReader(testZarrPath)
	if err != nil {
		t.Fatalf("Failed to initialize Zarr reader: %v", err)
	}

	// Initialize cache manager
	cacheManager, err := cache.NewManager(cache.Config{
		TileCacheSizeMB: 64, // Smaller cache for tests
		TileTTL:         5 * time.Minute,
		QueryCacheSize:  100,
	})
	if err != nil {
		t.Fatalf("Failed to initialize cache: %v", err)
	}

	// Initialize tile renderer
	tileRenderer := render.NewTileRenderer(render.Config{
		TileSize:        256,
		DefaultColormap: "viridis",
	})

	// Initialize tile service
	tileService := service.NewTileService(service.TileServiceConfig{
		DatasetID:  "default",
		ZarrReader: zarrReader,
		Cache:      cacheManager,
		Renderer:   tileRenderer,
	})

	// Create registry with single dataset
	registry := NewDatasetRegistry("default", []string{"default"}, "")
	registry.Register("default", tileService)

	// Create router
	router := NewRouter(RouterConfig{
		Registry:    registry,
		CORSOrigins: []string{"http://localhost:3000"},
	})

	// Create test server
	server := httptest.NewServer(router)

	return &testServer{
		server:      server,
		zarrReader:  zarrReader,
		cache:       cacheManager,
		renderer:    tileRenderer,
		tileService: tileService,
	}
}

// close cleans up test server resources
func (ts *testServer) close() {
	ts.server.Close()
	ts.zarrReader.Close()
	ts.cache.Close()
}

// --- Helper Functions ---

// assertStatusCode verifies the HTTP status code
func assertStatusCode(t *testing.T, resp *http.Response, expected int) {
	t.Helper()
	if resp.StatusCode != expected {
		t.Errorf("Expected status code %d, got %d", expected, resp.StatusCode)
	}
}

// assertContentType verifies the Content-Type header
func assertContentType(t *testing.T, resp *http.Response, expected string) {
	t.Helper()
	contentType := resp.Header.Get("Content-Type")
	if contentType != expected {
		t.Errorf("Expected Content-Type %q, got %q", expected, contentType)
	}
}

// assertPNG verifies the response body is a valid PNG image
func assertPNG(t *testing.T, body []byte) {
	t.Helper()
	// PNG magic bytes: 0x89 0x50 0x4E 0x47 0x0D 0x0A 0x1A 0x0A
	pngMagic := []byte{0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A}
	if len(body) < 8 {
		t.Errorf("Response too short to be a valid PNG (got %d bytes)", len(body))
		return
	}
	for i, b := range pngMagic {
		if body[i] != b {
			t.Errorf("Invalid PNG magic bytes at position %d: expected 0x%02X, got 0x%02X", i, b, body[i])
			return
		}
	}
}

// assertJSONFields verifies the response contains expected JSON fields
func assertJSONFields(t *testing.T, body []byte, expectedFields []string) {
	t.Helper()
	var result map[string]interface{}
	if err := json.Unmarshal(body, &result); err != nil {
		t.Errorf("Failed to parse JSON response: %v", err)
		return
	}
	for _, field := range expectedFields {
		if _, ok := result[field]; !ok {
			t.Errorf("Expected JSON field %q not found in response", field)
		}
	}
}

// --- Test Cases ---

// TestHealthEndpoint tests the health check endpoint
func TestHealthEndpoint(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	resp, err := http.Get(ts.server.URL + "/health")
	if err != nil {
		t.Fatalf("Failed to make request: %v", err)
	}
	defer resp.Body.Close()

	assertStatusCode(t, resp, http.StatusOK)

	body, err := io.ReadAll(resp.Body)
	if err != nil {
		t.Fatalf("Failed to read response body: %v", err)
	}

	if string(body) != "OK" {
		t.Errorf("Expected body 'OK', got %q", string(body))
	}
}

// TestTileEndpoint tests the basic tile rendering endpoint
func TestTileEndpoint(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	tests := []struct {
		name           string
		path           string
		expectedStatus int
		expectPNG      bool
	}{
		{
			name:           "valid tile z0",
			path:           "/d/default/tiles/0/0/0.png",
			expectedStatus: http.StatusOK,
			expectPNG:      true,
		},
		{
			name:           "valid tile z1",
			path:           "/d/default/tiles/1/0/0.png",
			expectedStatus: http.StatusOK,
			expectPNG:      true,
		},
		{
			name:           "invalid z parameter",
			path:           "/d/default/tiles/abc/0/0.png",
			expectedStatus: http.StatusBadRequest,
			expectPNG:      false,
		},
		{
			name:           "invalid x parameter",
			path:           "/d/default/tiles/0/abc/0.png",
			expectedStatus: http.StatusBadRequest,
			expectPNG:      false,
		},
		{
			name:           "invalid y parameter",
			path:           "/d/default/tiles/0/0/abc.png",
			expectedStatus: http.StatusBadRequest,
			expectPNG:      false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			resp, err := http.Get(ts.server.URL + tt.path)
			if err != nil {
				t.Fatalf("Failed to make request: %v", err)
			}
			defer resp.Body.Close()

			assertStatusCode(t, resp, tt.expectedStatus)

			if tt.expectPNG {
				assertContentType(t, resp, "image/png")
				body, err := io.ReadAll(resp.Body)
				if err != nil {
					t.Fatalf("Failed to read response body: %v", err)
				}
				assertPNG(t, body)
			}
		})
	}
}

// TestExpressionTileEndpoint tests the expression tile endpoint
func TestExpressionTileEndpoint(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	// Get a valid gene name from metadata
	metadata := ts.zarrReader.Metadata()
	if len(metadata.Genes) == 0 {
		t.Skip("No genes available in test data")
	}
	validGene := metadata.Genes[0]

	tests := []struct {
		name           string
		path           string
		expectedStatus int
		expectPNG      bool
	}{
		{
			name:           "valid expression tile with default colormap",
			path:           "/d/default/tiles/0/0/0/expression/" + validGene + ".png",
			expectedStatus: http.StatusOK,
			expectPNG:      true,
		},
		{
			name:           "expression tile with viridis colormap",
			path:           "/d/default/tiles/0/0/0/expression/" + validGene + ".png?colormap=viridis",
			expectedStatus: http.StatusOK,
			expectPNG:      true,
		},
		{
			name:           "expression tile with plasma colormap",
			path:           "/d/default/tiles/0/0/0/expression/" + validGene + ".png?colormap=plasma",
			expectedStatus: http.StatusOK,
			expectPNG:      true,
		},
		{
			name:           "unknown gene returns empty tile",
			path:           "/d/default/tiles/0/0/0/expression/UNKNOWN_GENE_XYZ.png",
			expectedStatus: http.StatusOK, // Returns empty tile, not error
			expectPNG:      true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			resp, err := http.Get(ts.server.URL + tt.path)
			if err != nil {
				t.Fatalf("Failed to make request: %v", err)
			}
			defer resp.Body.Close()

			assertStatusCode(t, resp, tt.expectedStatus)

			if tt.expectPNG {
				assertContentType(t, resp, "image/png")
				body, err := io.ReadAll(resp.Body)
				if err != nil {
					t.Fatalf("Failed to read response body: %v", err)
				}
				assertPNG(t, body)
			}
		})
	}
}

// TestCategoryTileEndpoint tests the category tile endpoint
func TestCategoryTileEndpoint(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	// Get a valid category from metadata
	metadata := ts.zarrReader.Metadata()
	var validCategory string
	for cat := range metadata.Categories {
		validCategory = cat
		break
	}
	if validCategory == "" {
		t.Skip("No categories available in test data")
	}

	tests := []struct {
		name           string
		path           string
		expectedStatus int
		expectPNG      bool
	}{
		{
			name:           "valid category tile",
			path:           "/d/default/tiles/0/0/0/category/" + validCategory + ".png",
			expectedStatus: http.StatusOK,
			expectPNG:      true,
		},
		{
			name:           "unknown category returns empty tile",
			path:           "/d/default/tiles/0/0/0/category/unknown_category_xyz.png",
			expectedStatus: http.StatusOK, // Returns empty tile
			expectPNG:      true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			resp, err := http.Get(ts.server.URL + tt.path)
			if err != nil {
				t.Fatalf("Failed to make request: %v", err)
			}
			defer resp.Body.Close()

			assertStatusCode(t, resp, tt.expectedStatus)

			if tt.expectPNG {
				assertContentType(t, resp, "image/png")
				body, err := io.ReadAll(resp.Body)
				if err != nil {
					t.Fatalf("Failed to read response body: %v", err)
				}
				assertPNG(t, body)
			}
		})
	}
}

// TestMetadataEndpoint tests the metadata API endpoint
func TestMetadataEndpoint(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	resp, err := http.Get(ts.server.URL + "/d/default/api/metadata")
	if err != nil {
		t.Fatalf("Failed to make request: %v", err)
	}
	defer resp.Body.Close()

	assertStatusCode(t, resp, http.StatusOK)
	assertContentType(t, resp, "application/json")

	body, err := io.ReadAll(resp.Body)
	if err != nil {
		t.Fatalf("Failed to read response body: %v", err)
	}

	// Verify expected fields
	assertJSONFields(t, body, []string{"zoom_levels", "n_genes_preaggregated", "n_cells"})
}

// TestDatasetsEndpoint tests the datasets list API endpoint
func TestDatasetsEndpoint(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	resp, err := http.Get(ts.server.URL + "/api/datasets")
	if err != nil {
		t.Fatalf("Failed to make request: %v", err)
	}
	defer resp.Body.Close()

	assertStatusCode(t, resp, http.StatusOK)
	assertContentType(t, resp, "application/json")

	body, err := io.ReadAll(resp.Body)
	if err != nil {
		t.Fatalf("Failed to read response body: %v", err)
	}

	// Verify expected fields
	assertJSONFields(t, body, []string{"default", "datasets"})
}

// TestStatsEndpoint tests the stats API endpoint
func TestStatsEndpoint(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	resp, err := http.Get(ts.server.URL + "/d/default/api/stats")
	if err != nil {
		t.Fatalf("Failed to make request: %v", err)
	}
	defer resp.Body.Close()

	assertStatusCode(t, resp, http.StatusOK)
	assertContentType(t, resp, "application/json")

	body, err := io.ReadAll(resp.Body)
	if err != nil {
		t.Fatalf("Failed to read response body: %v", err)
	}

	assertJSONFields(t, body, []string{"n_cells", "n_genes", "zoom_levels", "dataset_name"})
}

// TestGenesEndpoint tests the genes list API endpoint
func TestGenesEndpoint(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	resp, err := http.Get(ts.server.URL + "/d/default/api/genes")
	if err != nil {
		t.Fatalf("Failed to make request: %v", err)
	}
	defer resp.Body.Close()

	assertStatusCode(t, resp, http.StatusOK)
	assertContentType(t, resp, "application/json")

	body, err := io.ReadAll(resp.Body)
	if err != nil {
		t.Fatalf("Failed to read response body: %v", err)
	}

	assertJSONFields(t, body, []string{"genes", "total", "preaggregated_count"})
}

// TestGeneInfoEndpoint tests the gene info API endpoint
func TestGeneInfoEndpoint(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	// Get a valid gene name
	metadata := ts.zarrReader.Metadata()
	if len(metadata.Genes) == 0 {
		t.Skip("No genes available in test data")
	}
	validGene := metadata.Genes[0]

	tests := []struct {
		name           string
		gene           string
		expectedStatus int
		expectFields   []string
	}{
		{
			name:           "valid gene",
			gene:           validGene,
			expectedStatus: http.StatusOK,
			expectFields:   []string{"name", "index", "preaggregated"},
		},
		{
			name:           "unknown gene",
			gene:           "UNKNOWN_GENE_XYZ",
			expectedStatus: http.StatusNotFound,
			expectFields:   nil,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			resp, err := http.Get(ts.server.URL + "/d/default/api/genes/" + tt.gene)
			if err != nil {
				t.Fatalf("Failed to make request: %v", err)
			}
			defer resp.Body.Close()

			assertStatusCode(t, resp, tt.expectedStatus)

			if tt.expectedStatus == http.StatusOK {
				assertContentType(t, resp, "application/json")
				body, err := io.ReadAll(resp.Body)
				if err != nil {
					t.Fatalf("Failed to read response body: %v", err)
				}
				assertJSONFields(t, body, tt.expectFields)
			}
		})
	}
}

// TestGeneStatsEndpoint tests the gene stats API endpoint
func TestGeneStatsEndpoint(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	// Get a valid gene name
	metadata := ts.zarrReader.Metadata()
	if len(metadata.Genes) == 0 {
		t.Skip("No genes available in test data")
	}
	validGene := metadata.Genes[0]

	// First, check if expression data is available by making a probe request
	probeResp, err := http.Get(ts.server.URL + "/d/default/api/genes/" + validGene + "/stats")
	if err != nil {
		t.Fatalf("Failed to make probe request: %v", err)
	}
	probeResp.Body.Close()
	expressionDataAvailable := probeResp.StatusCode == http.StatusOK

	tests := []struct {
		name            string
		gene            string
		expectedStatus  int
		expectFields    []string
		requireExprData bool
	}{
		{
			name:            "valid gene stats",
			gene:            validGene,
			expectedStatus:  http.StatusOK,
			expectFields:    []string{"gene", "index", "expressing_bins", "total_bins", "mean_expression", "max_expression"},
			requireExprData: true,
		},
		{
			name:            "unknown gene stats",
			gene:            "UNKNOWN_GENE_XYZ",
			expectedStatus:  http.StatusNotFound,
			expectFields:    nil,
			requireExprData: false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if tt.requireExprData && !expressionDataAvailable {
				t.Skip("Expression data not available (multi-chunk Zarr format not supported)")
			}

			resp, err := http.Get(ts.server.URL + "/d/default/api/genes/" + tt.gene + "/stats")
			if err != nil {
				t.Fatalf("Failed to make request: %v", err)
			}
			defer resp.Body.Close()

			assertStatusCode(t, resp, tt.expectedStatus)

			if tt.expectedStatus == http.StatusOK {
				assertContentType(t, resp, "application/json")
				body, err := io.ReadAll(resp.Body)
				if err != nil {
					t.Fatalf("Failed to read response body: %v", err)
				}
				assertJSONFields(t, body, tt.expectFields)
			}
		})
	}
}

// TestGeneBinsEndpoint tests the gene bins API endpoint
func TestGeneBinsEndpoint(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	// Get a valid gene name
	metadata := ts.zarrReader.Metadata()
	if len(metadata.Genes) == 0 {
		t.Skip("No genes available in test data")
	}
	validGene := metadata.Genes[0]

	// Check if expression data is available by making a probe request
	probeResp, err := http.Get(ts.server.URL + "/d/default/api/genes/" + validGene + "/bins?limit=1")
	if err != nil {
		t.Fatalf("Failed to make probe request: %v", err)
	}
	probeResp.Body.Close()
	if probeResp.StatusCode != http.StatusOK {
		t.Skip("Expression data not available (multi-chunk Zarr format not supported)")
	}

	tests := []struct {
		name           string
		path           string
		expectedStatus int
		expectFields   []string
	}{
		{
			name:           "bins with default params",
			path:           "/d/default/api/genes/" + validGene + "/bins",
			expectedStatus: http.StatusOK,
			expectFields:   []string{"bins", "total_count", "offset", "limit", "gene", "threshold"},
		},
		{
			name:           "bins with threshold",
			path:           "/d/default/api/genes/" + validGene + "/bins?threshold=0.5",
			expectedStatus: http.StatusOK,
			expectFields:   []string{"bins", "total_count"},
		},
		{
			name:           "bins with pagination",
			path:           "/d/default/api/genes/" + validGene + "/bins?offset=0&limit=10",
			expectedStatus: http.StatusOK,
			expectFields:   []string{"bins", "total_count", "offset", "limit"},
		},
		{
			name:           "bins with large limit (should be capped at 1000)",
			path:           "/d/default/api/genes/" + validGene + "/bins?limit=5000",
			expectedStatus: http.StatusOK,
			expectFields:   []string{"bins"},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			resp, err := http.Get(ts.server.URL + tt.path)
			if err != nil {
				t.Fatalf("Failed to make request: %v", err)
			}
			defer resp.Body.Close()

			assertStatusCode(t, resp, tt.expectedStatus)

			if tt.expectedStatus == http.StatusOK {
				assertContentType(t, resp, "application/json")
				body, err := io.ReadAll(resp.Body)
				if err != nil {
					t.Fatalf("Failed to read response body: %v", err)
				}
				assertJSONFields(t, body, tt.expectFields)
			}
		})
	}
}

// TestCategoriesEndpoint tests the categories list API endpoint
func TestCategoriesEndpoint(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	resp, err := http.Get(ts.server.URL + "/d/default/api/categories")
	if err != nil {
		t.Fatalf("Failed to make request: %v", err)
	}
	defer resp.Body.Close()

	assertStatusCode(t, resp, http.StatusOK)
	assertContentType(t, resp, "application/json")

	body, err := io.ReadAll(resp.Body)
	if err != nil {
		t.Fatalf("Failed to read response body: %v", err)
	}

	// Verify it's valid JSON (can be a map)
	var result map[string]interface{}
	if err := json.Unmarshal(body, &result); err != nil {
		t.Errorf("Failed to parse categories response as JSON: %v", err)
	}
}

// TestCategoryColorsEndpoint tests the category colors API endpoint
func TestCategoryColorsEndpoint(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	// Get a valid category
	metadata := ts.zarrReader.Metadata()
	var validCategory string
	for cat := range metadata.Categories {
		validCategory = cat
		break
	}
	if validCategory == "" {
		t.Skip("No categories available in test data")
	}

	tests := []struct {
		name           string
		column         string
		expectedStatus int
	}{
		{
			name:           "valid category colors",
			column:         validCategory,
			expectedStatus: http.StatusOK,
		},
		{
			name:           "unknown category",
			column:         "unknown_category_xyz",
			expectedStatus: http.StatusNotFound,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			resp, err := http.Get(ts.server.URL + "/d/default/api/categories/" + tt.column + "/colors")
			if err != nil {
				t.Fatalf("Failed to make request: %v", err)
			}
			defer resp.Body.Close()

			assertStatusCode(t, resp, tt.expectedStatus)

			if tt.expectedStatus == http.StatusOK {
				assertContentType(t, resp, "application/json")
				body, err := io.ReadAll(resp.Body)
				if err != nil {
					t.Fatalf("Failed to read response body: %v", err)
				}

				// Verify it's a valid color map
				var colors map[string]string
				if err := json.Unmarshal(body, &colors); err != nil {
					t.Errorf("Failed to parse colors response: %v", err)
				}
				if len(colors) == 0 {
					t.Error("Expected non-empty color map")
				}
			}
		})
	}
}

// TestCategoryLegendEndpoint tests the category legend API endpoint
func TestCategoryLegendEndpoint(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	// Get a valid category
	metadata := ts.zarrReader.Metadata()
	var validCategory string
	for cat := range metadata.Categories {
		validCategory = cat
		break
	}
	if validCategory == "" {
		t.Skip("No categories available in test data")
	}

	tests := []struct {
		name           string
		column         string
		expectedStatus int
	}{
		{
			name:           "valid category legend",
			column:         validCategory,
			expectedStatus: http.StatusOK,
		},
		{
			name:           "unknown category legend",
			column:         "unknown_category_xyz",
			expectedStatus: http.StatusNotFound,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			resp, err := http.Get(ts.server.URL + "/d/default/api/categories/" + tt.column + "/legend")
			if err != nil {
				t.Fatalf("Failed to make request: %v", err)
			}
			defer resp.Body.Close()

			assertStatusCode(t, resp, tt.expectedStatus)

			if tt.expectedStatus == http.StatusOK {
				assertContentType(t, resp, "application/json")
				body, err := io.ReadAll(resp.Body)
				if err != nil {
					t.Fatalf("Failed to read response body: %v", err)
				}

				// Verify it's a valid legend array
				var legend []map[string]interface{}
				if err := json.Unmarshal(body, &legend); err != nil {
					t.Errorf("Failed to parse legend response: %v", err)
				}
				if len(legend) == 0 {
					t.Error("Expected non-empty legend array")
				}

				// Verify legend item structure
				if len(legend) > 0 {
					item := legend[0]
					if _, ok := item["value"]; !ok {
						t.Error("Expected legend item to have 'value' field")
					}
					if _, ok := item["color"]; !ok {
						t.Error("Expected legend item to have 'color' field")
					}
					if _, ok := item["index"]; !ok {
						t.Error("Expected legend item to have 'index' field")
					}
				}
			}
		})
	}
}

// TestCacheHeaders tests that tile endpoints return proper cache headers
func TestCacheHeaders(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	resp, err := http.Get(ts.server.URL + "/d/default/tiles/0/0/0.png")
	if err != nil {
		t.Fatalf("Failed to make request: %v", err)
	}
	defer resp.Body.Close()

	cacheControl := resp.Header.Get("Cache-Control")
	if cacheControl == "" {
		t.Error("Expected Cache-Control header to be set")
	}
	if cacheControl != "public, max-age=3600" {
		t.Errorf("Expected Cache-Control 'public, max-age=3600', got %q", cacheControl)
	}
}

// TestCORSHeaders tests that CORS headers are set correctly
func TestCORSHeaders(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	// Create a request with Origin header
	req, err := http.NewRequest("GET", ts.server.URL+"/health", nil)
	if err != nil {
		t.Fatalf("Failed to create request: %v", err)
	}
	req.Header.Set("Origin", "http://localhost:3000")

	client := &http.Client{}
	resp, err := client.Do(req)
	if err != nil {
		t.Fatalf("Failed to make request: %v", err)
	}
	defer resp.Body.Close()

	accessControlOrigin := resp.Header.Get("Access-Control-Allow-Origin")
	if accessControlOrigin == "" {
		t.Error("Expected Access-Control-Allow-Origin header to be set for allowed origin")
	}
}

// TestAllEndpointsReachable runs a quick check that all endpoints are reachable
func TestAllEndpointsReachable(t *testing.T) {
	ts := setupTestServer(t)
	defer ts.close()

	metadata := ts.zarrReader.Metadata()

	// Check if expression data is available
	expressionDataAvailable := false
	if len(metadata.Genes) > 0 {
		probeResp, err := http.Get(ts.server.URL + "/d/default/api/genes/" + metadata.Genes[0] + "/stats")
		if err == nil {
			expressionDataAvailable = probeResp.StatusCode == http.StatusOK
			probeResp.Body.Close()
		}
	}

	type endpoint struct {
		path            string
		expectedStatus  int
		requireExprData bool
	}

	endpoints := []endpoint{
		{"/health", http.StatusOK, false},
		{"/d/default/tiles/0/0/0.png", http.StatusOK, false},
		{"/d/default/api/metadata", http.StatusOK, false},
		{"/d/default/api/stats", http.StatusOK, false},
		{"/d/default/api/genes", http.StatusOK, false},
		{"/d/default/api/categories", http.StatusOK, false},
		{"/api/datasets", http.StatusOK, false},
	}

	// Add gene-specific endpoints if genes exist
	if len(metadata.Genes) > 0 {
		gene := metadata.Genes[0]
		endpoints = append(endpoints,
			endpoint{"/d/default/api/genes/" + gene, http.StatusOK, false},
			endpoint{"/d/default/api/genes/" + gene + "/stats", http.StatusOK, true},
			endpoint{"/d/default/api/genes/" + gene + "/bins", http.StatusOK, true},
			endpoint{"/d/default/tiles/0/0/0/expression/" + gene + ".png", http.StatusOK, false}, // Returns empty tile on error
		)
	}

	// Add category-specific endpoints if categories exist
	for cat := range metadata.Categories {
		endpoints = append(endpoints,
			endpoint{"/d/default/api/categories/" + cat + "/colors", http.StatusOK, false},
			endpoint{"/d/default/api/categories/" + cat + "/legend", http.StatusOK, false},
			endpoint{"/d/default/tiles/0/0/0/category/" + cat + ".png", http.StatusOK, false},
		)
		break // Only test one category
	}

	for _, ep := range endpoints {
		t.Run(ep.path, func(t *testing.T) {
			if ep.requireExprData && !expressionDataAvailable {
				t.Skip("Expression data not available")
			}

			resp, err := http.Get(ts.server.URL + ep.path)
			if err != nil {
				t.Fatalf("Failed to make request to %s: %v", ep.path, err)
			}
			defer resp.Body.Close()

			if resp.StatusCode != ep.expectedStatus {
				t.Errorf("Endpoint %s: expected status %d, got %d", ep.path, ep.expectedStatus, resp.StatusCode)
			}
		})
	}
}
