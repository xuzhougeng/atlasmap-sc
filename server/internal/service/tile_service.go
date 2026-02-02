// Package service provides business logic for the tile server.
package service

import (
	"fmt"
	"sort"
	"sync"

	"github.com/atlasmap-sc/server/internal/cache"
	"github.com/atlasmap-sc/server/internal/data/soma"
	"github.com/atlasmap-sc/server/internal/data/zarr"
	"github.com/atlasmap-sc/server/internal/render"
)

// TileServiceConfig contains tile service configuration.
type TileServiceConfig struct {
	DatasetID  string
	CoordKey   string
	ZarrReader *zarr.Reader
	SomaReader *soma.Reader
	Cache      *cache.Manager
	Renderer   *render.TileRenderer
}

// TileService handles tile rendering and serving.
type TileService struct {
	datasetID string
	coordKey  string
	zarr      *zarr.Reader
	soma      *soma.Reader
	cache     *cache.Manager
	renderer  *render.TileRenderer

	renderZoom int
	tileSize   int // tile size in pixels (typically 256)

	// Per-zoom bins cache
	binsZoomMu    sync.Mutex
	binsZoomCache map[int][]zarr.Bin
	binsZoomErr   map[int]error

	// Legacy single-zoom bins cache (for backward compatibility with functions that use renderZoom)
	binsOnce sync.Once
	bins     []zarr.Bin
	binsErr  error

	// Per-zoom expression cache: key = "zoom:gene"
	exprMu    sync.Mutex
	exprCache map[string]expressionCacheEntry

	// Per-zoom category cache: key = "zoom:column"
	catMu    sync.Mutex
	catCache map[string][]int

	geneCategoryMeansMu    sync.Mutex
	geneCategoryMeansCache map[string][]GeneCategoryMeanItem

	categoryCentroidsMu    sync.Mutex
	categoryCentroidsCache map[string][]CategoryCentroidItem

	// Spatial index for cell-level queries (lazy loaded)
	spatialIndexOnce sync.Once
	spatialIndex     *soma.SpatialIndex
	spatialIndexErr  error
}

// NewTileService creates a new tile service.
func NewTileService(cfg TileServiceConfig) *TileService {
	renderZoom := 0
	tileSize := 256 // default tile size
	if md := cfg.ZarrReader.Metadata(); md != nil {
		if md.ZoomLevels > 0 {
			renderZoom = md.ZoomLevels - 1
		}
		if md.TileSize > 0 {
			tileSize = md.TileSize
		}
	}

	datasetID := cfg.DatasetID
	if datasetID == "" {
		datasetID = "default"
	}

	return &TileService{
		datasetID:              datasetID,
		coordKey:               cfg.CoordKey,
		zarr:                   cfg.ZarrReader,
		soma:                   cfg.SomaReader,
		cache:                  cfg.Cache,
		renderer:               cfg.Renderer,
		renderZoom:             renderZoom,
		tileSize:               tileSize,
		binsZoomCache:          make(map[int][]zarr.Bin),
		binsZoomErr:            make(map[int]error),
		exprCache:              make(map[string]expressionCacheEntry),
		catCache:               make(map[string][]int),
		geneCategoryMeansCache: make(map[string][]GeneCategoryMeanItem),
		categoryCentroidsCache: make(map[string][]CategoryCentroidItem),
	}
}

type expressionCacheEntry struct {
	values  []float32
	min     float32
	max     float32
	autoMax float32
}

type GeneCategoryMeanItem struct {
	Value          string  `json:"value"`
	Color          string  `json:"color"`
	Index          int     `json:"index"`
	BinCount       int     `json:"bin_count"`
	MeanExpression float32 `json:"mean_expression"`
}

func (s *TileService) loadBins() error {
	s.binsOnce.Do(func() {
		s.bins, s.binsErr = s.zarr.GetBins(s.renderZoom, 0, 0)
	})
	return s.binsErr
}

// loadBinsForZoom lazily loads bins for a specific zoom level.
func (s *TileService) loadBinsForZoom(zoom int) ([]zarr.Bin, error) {
	s.binsZoomMu.Lock()
	defer s.binsZoomMu.Unlock()

	// Check if already cached
	if bins, ok := s.binsZoomCache[zoom]; ok {
		if err, ok := s.binsZoomErr[zoom]; ok && err != nil {
			return nil, err
		}
		return bins, nil
	}

	// Load from Zarr
	bins, err := s.zarr.GetBins(zoom, 0, 0)
	s.binsZoomCache[zoom] = bins
	s.binsZoomErr[zoom] = err
	return bins, err
}

// sourceZoomForTile computes the appropriate source zoom level for rendering a tile.
// This ensures binSize >= 1px by limiting sourceZoom so that binsPerTileAxis <= tileSize.
// Formula: sourceZoom = min(renderZoom, mapZoom + tilePixelsLog2)
// where tilePixelsLog2 = log2(tileSize).
func (s *TileService) sourceZoomForTile(mapZoom int) int {
	// Compute log2 of tile size (for 256, this is 8)
	tilePixelsLog2 := 0
	ts := s.tileSize
	for ts > 1 {
		ts >>= 1
		tilePixelsLog2++
	}

	// sourceZoom = min(renderZoom, mapZoom + tilePixelsLog2)
	sourceZoom := mapZoom + tilePixelsLog2
	if sourceZoom > s.renderZoom {
		sourceZoom = s.renderZoom
	}
	return sourceZoom
}

// binsForTile returns bins for a tile, along with their indices, binsPerTileAxis, flippedTileY, sourceZoom, and any error.
// sourceZoom is the zoom level used to load the bins (may be < renderZoom to ensure binSize >= 1px).
func (s *TileService) binsForTile(mapZoom, tileX, tileY int) ([]zarr.Bin, []int, int, int, int, error) {
	if mapZoom < 0 || mapZoom > s.renderZoom {
		return nil, nil, 0, 0, 0, fmt.Errorf("invalid zoom level: %d", mapZoom)
	}
	if tileX < 0 || tileY < 0 {
		return nil, nil, 0, 0, 0, fmt.Errorf("invalid tile coordinates: %d/%d", tileX, tileY)
	}

	tilesPerAxis := 1 << mapZoom
	if tileX >= tilesPerAxis || tileY >= tilesPerAxis {
		return nil, nil, 0, 0, 0, fmt.Errorf("tile out of range: %d/%d (tiles_per_axis=%d)", tileX, tileY, tilesPerAxis)
	}

	// Compute sourceZoom: use a zoom level that ensures binSize >= 1px
	sourceZoom := s.sourceZoomForTile(mapZoom)

	// Load bins for the source zoom level
	bins, err := s.loadBinsForZoom(sourceZoom)
	if err != nil {
		return nil, nil, 0, 0, 0, err
	}

	// binsPerTileAxis at the source zoom level relative to the requested mapZoom
	binsPerTileAxis := 1 << (sourceZoom - mapZoom)
	startX := int32(tileX * binsPerTileAxis)
	endX := int32((tileX + 1) * binsPerTileAxis)
	// Flip Y so that tileY=0 corresponds to the top (max-Y) of the dataset,
	// matching the conventional Cartesian orientation for UMAP plots.
	flippedTileY := (tilesPerAxis - 1) - tileY
	startY := int32(flippedTileY * binsPerTileAxis)
	endY := int32((flippedTileY + 1) * binsPerTileAxis)

	outBins := make([]zarr.Bin, 0, 64)
	outIdx := make([]int, 0, 64)
	for i, bin := range bins {
		if bin.X < startX || bin.X >= endX {
			continue
		}
		if bin.Y < startY || bin.Y >= endY {
			continue
		}
		outBins = append(outBins, bin)
		outIdx = append(outIdx, i)
	}

	return outBins, outIdx, binsPerTileAxis, flippedTileY, sourceZoom, nil
}

// expressionForGene loads expression for a gene at renderZoom (backward compatible).
func (s *TileService) expressionForGene(gene string) (expressionCacheEntry, error) {
	return s.expressionForGeneAtZoom(gene, s.renderZoom)
}

// expressionForGeneAtZoom loads expression for a gene at a specific zoom level.
func (s *TileService) expressionForGeneAtZoom(gene string, zoom int) (expressionCacheEntry, error) {
	s.exprMu.Lock()
	defer s.exprMu.Unlock()

	// Cache key includes zoom level
	cacheKey := fmt.Sprintf("%d:%s", zoom, gene)
	if cached, ok := s.exprCache[cacheKey]; ok {
		return cached, nil
	}

	geneIdx, ok := s.zarr.Metadata().GeneIndex[gene]
	if !ok {
		return expressionCacheEntry{}, fmt.Errorf("gene not found: %s", gene)
	}

	values, err := s.zarr.GetExpression(zoom, geneIdx, "mean")
	if err != nil {
		return expressionCacheEntry{}, err
	}

	if len(values) == 0 {
		cached := expressionCacheEntry{values: values, min: 0, max: 0, autoMax: 0}
		s.exprCache[cacheKey] = cached
		return cached, nil
	}

	minV := values[0]
	maxV := values[0]
	expressing := 0
	if values[0] > 0 {
		expressing = 1
	}
	for _, v := range values[1:] {
		if v < minV {
			minV = v
		}
		if v > maxV {
			maxV = v
		}
		if v > 0 {
			expressing++
		}
	}

	// Auto max is a robust upper bound to avoid outliers compressing the colormap.
	// Use the 80th percentile of non-zero values (nearest-rank).
	autoMax := float32(0)
	if expressing > 0 {
		expressingValues := make([]float32, 0, expressing)
		for _, v := range values {
			if v > 0 {
				expressingValues = append(expressingValues, v)
			}
		}
		sort.Slice(expressingValues, func(i, j int) bool { return expressingValues[i] < expressingValues[j] })
		n := len(expressingValues)
		if n > 0 {
			// idx = ceil(0.80*n) - 1, computed with integers.
			idx := (80*n+99)/100 - 1
			if idx < 0 {
				idx = 0
			} else if idx >= n {
				idx = n - 1
			}
			autoMax = expressingValues[idx]
		}
	}
	if autoMax <= 0 {
		autoMax = maxV
	}

	cached := expressionCacheEntry{values: values, min: minV, max: maxV, autoMax: autoMax}
	s.exprCache[cacheKey] = cached
	return cached, nil
}

// categoryForColumn loads category data at renderZoom (backward compatible).
func (s *TileService) categoryForColumn(column string) ([]int, error) {
	return s.categoryForColumnAtZoom(column, s.renderZoom)
}

// categoryForColumnAtZoom loads category data at a specific zoom level.
func (s *TileService) categoryForColumnAtZoom(column string, zoom int) ([]int, error) {
	s.catMu.Lock()
	defer s.catMu.Unlock()

	// Cache key includes zoom level
	cacheKey := fmt.Sprintf("%d:%s", zoom, column)
	if cached, ok := s.catCache[cacheKey]; ok {
		return cached, nil
	}

	values, err := s.zarr.GetCategoryDataForBins(zoom, column)
	if err != nil {
		return nil, err
	}
	s.catCache[cacheKey] = values
	return values, nil
}

// GetTile returns a rendered tile PNG.
func (s *TileService) GetTile(z, x, y int) ([]byte, error) {
	// Check cache (prefix with dataset ID)
	cacheKey := s.datasetID + ":" + cache.TileKey(z, x, y, nil)
	if data, ok := s.cache.GetTile(cacheKey); ok {
		return data, nil
	}

	// Load bins for this tile (sourceZoom ensures binSize >= 1px)
	bins, _, binsPerTileAxis, renderTileY, _, err := s.binsForTile(z, x, y)
	if err != nil {
		return nil, fmt.Errorf("failed to load bins: %w", err)
	}

	// Render tile
	data, err := s.renderer.RenderTile(bins, binsPerTileAxis, x, renderTileY)
	if err != nil {
		return nil, fmt.Errorf("failed to render tile: %w", err)
	}

	// Cache result
	s.cache.SetTile(cacheKey, data)

	return data, nil
}

// GetExpressionTile returns a tile colored by gene expression.
func (s *TileService) GetExpressionTile(
	z, x, y int,
	gene string,
	colormap string,
	exprMin *float32,
	exprMax *float32,
	pointSize float64,
) ([]byte, error) {
	// Check cache (prefix with dataset ID)
	var (
		cacheMin    *float32
		cacheMax    *float32
		cacheMinVal float32
		cacheMaxVal float32
	)
	if exprMin != nil {
		cacheMinVal = *exprMin
		cacheMin = &cacheMinVal
	}
	if exprMax != nil {
		cacheMaxVal = *exprMax
		cacheMax = &cacheMaxVal
	}
	if cacheMin != nil && cacheMax != nil && cacheMaxVal < cacheMinVal {
		cacheMinVal, cacheMaxVal = cacheMaxVal, cacheMinVal
	}
	cacheKey := s.datasetID + ":" + cache.ExpressionTileKey(z, x, y, gene, colormap, cacheMin, cacheMax, pointSize)
	if data, ok := s.cache.GetTile(cacheKey); ok {
		return data, nil
	}

	// Load bins (sourceZoom ensures binSize >= 1px)
	bins, indices, binsPerTileAxis, renderTileY, sourceZoom, err := s.binsForTile(z, x, y)
	if err != nil {
		return nil, fmt.Errorf("failed to load bins: %w", err)
	}

	// Load expression at the same zoom level as bins
	exprAll, err := s.expressionForGeneAtZoom(gene, sourceZoom)
	if err != nil {
		return nil, fmt.Errorf("failed to load expression: %w", err)
	}

	exprTile := make([]float32, len(indices))
	for i, idx := range indices {
		if idx >= 0 && idx < len(exprAll.values) {
			exprTile[i] = exprAll.values[idx]
		}
	}

	minV := float32(0)
	maxV := exprAll.autoMax
	if maxV <= 0 {
		maxV = exprAll.max
	}
	if exprMin != nil {
		minV = *exprMin
	}
	if exprMax != nil {
		maxV = *exprMax
	}
	if maxV < minV {
		minV, maxV = maxV, minV
	}

	// Render tile
	data, err := s.renderer.RenderExpressionTile(
		bins,
		exprTile,
		minV,
		maxV,
		binsPerTileAxis,
		x, renderTileY,
		colormap,
		pointSize,
	)
	if err != nil {
		return nil, fmt.Errorf("failed to render tile: %w", err)
	}

	// Cache result
	s.cache.SetTile(cacheKey, data)

	return data, nil
}

// GetEmptyTile returns an empty tile.
func (s *TileService) GetEmptyTile() ([]byte, error) {
	return s.renderer.CreateEmptyTile()
}

// Metadata returns dataset metadata.
func (s *TileService) Metadata() *zarr.ZarrMetadata {
	md := s.zarr.Metadata()
	if md == nil {
		return nil
	}
	out := *md
	if s.coordKey != "" {
		out.CoordinateSystem = s.coordKey
	}
	return &out
}

func (s *TileService) Soma() *soma.Reader {
	return s.soma
}

// GetCategoryTile returns a tile colored by category.
func (s *TileService) GetCategoryTile(z, x, y int, column string, categoryFilter []string, pointSize float64) ([]byte, error) {
	// Check cache (prefix with dataset ID)
	cacheKey := s.datasetID + ":" + cache.CategoryTileKey(z, x, y, column, categoryFilter, pointSize)
	if data, ok := s.cache.GetTile(cacheKey); ok {
		return data, nil
	}

	// Load bins (sourceZoom ensures binSize >= 1px)
	bins, indices, binsPerTileAxis, renderTileY, sourceZoom, err := s.binsForTile(z, x, y)
	if err != nil {
		return nil, fmt.Errorf("failed to load bins: %w", err)
	}

	// Load category data at the same zoom level as bins
	categoryAll, err := s.categoryForColumnAtZoom(column, sourceZoom)
	if err != nil {
		return nil, fmt.Errorf("failed to load category data: %w", err)
	}

	// Build allowed category indices if filter is provided
	var allowedIndices map[int]bool
	if categoryFilter != nil {
		allowedIndices = make(map[int]bool)
		catInfo, ok := s.zarr.Metadata().Categories[column]
		if ok {
			for _, filterValue := range categoryFilter {
				if idx, exists := catInfo.Mapping[filterValue]; exists {
					allowedIndices[idx] = true
				}
			}
		}
	}

	categoryTile := make([]int, len(indices))
	for i, idx := range indices {
		if idx >= 0 && idx < len(categoryAll) {
			catValue := categoryAll[idx]
			// Apply filter: set to -1 if not in allowed list
			if allowedIndices != nil && !allowedIndices[catValue] {
				categoryTile[i] = -1
			} else {
				categoryTile[i] = catValue
			}
		}
	}

	// Render tile
	data, err := s.renderer.RenderCategoryTile(bins, categoryTile, binsPerTileAxis, x, renderTileY, pointSize)
	if err != nil {
		return nil, fmt.Errorf("failed to render tile: %w", err)
	}

	// Cache result
	s.cache.SetTile(cacheKey, data)

	return data, nil
}

// CategoryLegendItem represents a legend item for a category.
type CategoryLegendItem struct {
	Value     string `json:"value"`
	Color     string `json:"color"`
	Index     int    `json:"index"`
	CellCount int    `json:"cell_count"`
}

type CategoryCentroidItem struct {
	Value     string   `json:"value"`
	Color     string   `json:"color"`
	Index     int      `json:"index"`
	BinCount  int      `json:"bin_count"`
	CellCount int      `json:"cell_count"`
	X         *float64 `json:"x"`
	Y         *float64 `json:"y"`
}

// GetCategoryLegend returns legend data for a category column.
func (s *TileService) GetCategoryLegend(column string) ([]CategoryLegendItem, error) {
	catInfo, ok := s.zarr.Metadata().Categories[column]
	if !ok {
		return nil, fmt.Errorf("category not found: %s", column)
	}

	colors, err := s.zarr.GetCategoryColors(column)
	if err != nil {
		return nil, err
	}

	// Load bins and category data to calculate cell counts
	if err := s.loadBins(); err != nil {
		return nil, err
	}

	catAll, err := s.categoryForColumn(column)
	if err != nil {
		return nil, err
	}

	// Calculate cell counts per category
	nCats := len(catInfo.Values)
	cellCounts := make([]int, nCats)
	n := len(s.bins)
	if len(catAll) < n {
		n = len(catAll)
	}
	for i := 0; i < n; i++ {
		catIdx := catAll[i]
		if catIdx >= 0 && catIdx < nCats {
			cellCounts[catIdx] += int(s.bins[i].CellCount)
		}
	}

	legend := make([]CategoryLegendItem, len(catInfo.Values))
	for i, value := range catInfo.Values {
		legend[i] = CategoryLegendItem{
			Value:     value,
			Color:     colors[value],
			Index:     i,
			CellCount: cellCounts[i],
		}
	}

	return legend, nil
}

func (s *TileService) GetCategoryCentroids(column string) ([]CategoryCentroidItem, error) {
	s.categoryCentroidsMu.Lock()
	if cached, ok := s.categoryCentroidsCache[column]; ok {
		s.categoryCentroidsMu.Unlock()
		return cached, nil
	}
	s.categoryCentroidsMu.Unlock()

	md := s.zarr.Metadata()
	if md == nil {
		return nil, fmt.Errorf("metadata not available")
	}

	catInfo, ok := md.Categories[column]
	if !ok {
		return nil, fmt.Errorf("category not found: %s", column)
	}

	colors, err := s.GetCategoryColors(column)
	if err != nil {
		return nil, err
	}

	if err := s.loadBins(); err != nil {
		return nil, err
	}
	catAll, err := s.categoryForColumn(column)
	if err != nil {
		return nil, err
	}

	nCats := len(catInfo.Values)
	sumX := make([]float64, nCats)
	sumY := make([]float64, nCats)
	sumW := make([]float64, nCats)
	binCounts := make([]int, nCats)
	cellCounts := make([]int, nCats)

	binsPerAxis := 1 << s.renderZoom
	binSizeX := (md.Bounds.MaxX - md.Bounds.MinX) / float64(binsPerAxis)
	binSizeY := (md.Bounds.MaxY - md.Bounds.MinY) / float64(binsPerAxis)

	n := len(s.bins)
	if len(catAll) < n {
		n = len(catAll)
	}
	for i := 0; i < n; i++ {
		catIdx := catAll[i]
		if catIdx < 0 || catIdx >= nCats {
			continue
		}

		w := float64(s.bins[i].CellCount)
		if w <= 0 {
			continue
		}

		binCounts[catIdx]++
		cellCounts[catIdx] += int(s.bins[i].CellCount)

		x := md.Bounds.MinX + (float64(s.bins[i].X)+0.5)*binSizeX
		y := md.Bounds.MinY + (float64(s.bins[i].Y)+0.5)*binSizeY

		sumX[catIdx] += x * w
		sumY[catIdx] += y * w
		sumW[catIdx] += w
	}

	out := make([]CategoryCentroidItem, nCats)
	for idx, value := range catInfo.Values {
		var cx *float64
		var cy *float64
		if sumW[idx] > 0 {
			x := sumX[idx] / sumW[idx]
			y := sumY[idx] / sumW[idx]
			cx = &x
			cy = &y
		}
		out[idx] = CategoryCentroidItem{
			Value:     value,
			Color:     colors[value],
			Index:     idx,
			BinCount:  binCounts[idx],
			CellCount: cellCounts[idx],
			X:         cx,
			Y:         cy,
		}
	}

	s.categoryCentroidsMu.Lock()
	s.categoryCentroidsCache[column] = out
	s.categoryCentroidsMu.Unlock()

	return out, nil
}

// GetCategoryColors returns the color mapping for a category column.
func (s *TileService) GetCategoryColors(column string) (map[string]string, error) {
	return s.zarr.GetCategoryColors(column)
}

// GetGeneCategoryMeans returns mean gene expression per category value.
// Expression values are per-bin (from preprocessing) at the server's render zoom level.
func (s *TileService) GetGeneCategoryMeans(gene, column string) ([]GeneCategoryMeanItem, error) {
	cacheKey := gene + "\x00" + column

	s.geneCategoryMeansMu.Lock()
	if cached, ok := s.geneCategoryMeansCache[cacheKey]; ok {
		s.geneCategoryMeansMu.Unlock()
		return cached, nil
	}
	s.geneCategoryMeansMu.Unlock()

	md := s.zarr.Metadata()
	if md == nil {
		return nil, fmt.Errorf("metadata not available")
	}

	catInfo, ok := md.Categories[column]
	if !ok {
		return nil, fmt.Errorf("category not found: %s", column)
	}

	exprAll, err := s.expressionForGene(gene)
	if err != nil {
		return nil, err
	}

	catAll, err := s.categoryForColumn(column)
	if err != nil {
		return nil, err
	}

	nCats := len(catInfo.Values)
	sums := make([]float64, nCats)
	counts := make([]int, nCats)

	n := len(exprAll.values)
	if len(catAll) < n {
		n = len(catAll)
	}

	for i := 0; i < n; i++ {
		catIdx := catAll[i]
		if catIdx < 0 || catIdx >= nCats {
			continue
		}
		sums[catIdx] += float64(exprAll.values[i])
		counts[catIdx]++
	}

	colors, err := s.GetCategoryColors(column)
	if err != nil {
		return nil, err
	}

	out := make([]GeneCategoryMeanItem, nCats)
	for idx, value := range catInfo.Values {
		mean := float32(0)
		if counts[idx] > 0 {
			mean = float32(sums[idx] / float64(counts[idx]))
		}
		out[idx] = GeneCategoryMeanItem{
			Value:          value,
			Color:          colors[value],
			Index:          idx,
			BinCount:       counts[idx],
			MeanExpression: mean,
		}
	}

	s.geneCategoryMeansMu.Lock()
	s.geneCategoryMeansCache[cacheKey] = out
	s.geneCategoryMeansMu.Unlock()

	return out, nil
}

// QueryBinsExpressingGene returns bins where gene expression is above threshold.
func (s *TileService) QueryBinsExpressingGene(
	gene string,
	threshold float32,
	offset, limit int,
) (*zarr.BinQueryResult, error) {
	return s.zarr.GetBinsExpressingGene(gene, threshold, offset, limit)
}

// GetGeneStats returns statistics for a gene at a specific zoom level.
func (s *TileService) GetGeneStats(gene string, zoom int) (*zarr.GeneStats, error) {
	return s.zarr.GetGeneStats(gene, zoom)
}

// CellInfo represents a single cell with its position and optional attributes.
type CellInfo struct {
	JoinID     int64   `json:"joinid"`
	X          float32 `json:"x"`
	Y          float32 `json:"y"`
	Expression float32 `json:"expression,omitempty"`
	Category   string  `json:"category,omitempty"`
}

// CellQueryResult is the response for GetCellsInBounds.
type CellQueryResult struct {
	Cells      []CellInfo `json:"cells"`
	TotalCount int        `json:"total_count"`
	Truncated  bool       `json:"truncated"`
}

// loadSpatialIndex lazily loads the spatial index.
func (s *TileService) loadSpatialIndex() (*soma.SpatialIndex, error) {
	s.spatialIndexOnce.Do(func() {
		if s.soma == nil || !s.soma.Supported() {
			s.spatialIndexErr = fmt.Errorf("soma not configured or not supported")
			return
		}

		md := s.zarr.Metadata()
		if md == nil {
			s.spatialIndexErr = fmt.Errorf("metadata not available")
			return
		}

		coordRange := md.Bounds.MaxX - md.Bounds.MinX
		if coordRange <= 0 {
			coordRange = 256.0
		}

		s.spatialIndex, s.spatialIndexErr = s.soma.LoadSpatialIndex(s.renderZoom, coordRange)
	})
	return s.spatialIndex, s.spatialIndexErr
}

// GetCellsInBounds returns cells within the given bounding box.
// Options:
//   - gene: if set, include expression values for this gene
//   - category: if set, include category values for this column
//   - categoryFilter: if set, only return cells with category values in this list
//   - limit: maximum number of cells to return (default 5000)
func (s *TileService) GetCellsInBounds(
	minX, minY, maxX, maxY float64,
	gene string,
	category string,
	categoryFilter []string,
	limit int,
) (*CellQueryResult, error) {
	if limit <= 0 {
		limit = 5000
	}
	if limit > 50000 {
		limit = 50000
	}

	// Match tile filtering semantics:
	// - nil  => no filter (show all)
	// - []   => filter-to-none (show none)
	if categoryFilter != nil && len(categoryFilter) == 0 {
		return &CellQueryResult{
			Cells:      []CellInfo{},
			TotalCount: 0,
			Truncated:  false,
		}, nil
	}

	idx, err := s.loadSpatialIndex()
	if err != nil {
		return nil, err
	}

	// Query cells in bbox (with limit + 1 to detect truncation)
	cellJoinIDs := idx.QueryCellsInBounds(minX, minY, maxX, maxY, limit+1)
	truncated := len(cellJoinIDs) > limit
	if truncated {
		cellJoinIDs = cellJoinIDs[:limit]
	}

	if len(cellJoinIDs) == 0 {
		return &CellQueryResult{
			Cells:      []CellInfo{},
			TotalCount: 0,
			Truncated:  false,
		}, nil
	}

	// Get coordinates
	xMap, yMap := idx.GetCellCoordinates(cellJoinIDs)

	// Optionally load expression
	var exprMap map[int64]float32
	if gene != "" && s.soma != nil && s.soma.Supported() {
		_, exprMap, err = s.soma.ExpressionByCellJoinID(gene, cellJoinIDs)
		if err != nil {
			// Gene not found or other error - continue without expression
			exprMap = nil
		}
	}

	// Optionally load category
	var catMap map[int64]string
	if category != "" && s.soma != nil && s.soma.Supported() {
		catMap, err = s.soma.ReadObsCategoryColumn(category, cellJoinIDs)
		if err != nil {
			// Column not found or other error - continue without category
			catMap = nil
		}
	}

	// Build category filter set
	var catFilterSet map[string]struct{}
	if len(categoryFilter) > 0 && catMap != nil {
		catFilterSet = make(map[string]struct{}, len(categoryFilter))
		for _, v := range categoryFilter {
			catFilterSet[v] = struct{}{}
		}
	}

	// Build result
	cells := make([]CellInfo, 0, len(cellJoinIDs))
	for _, cid := range cellJoinIDs {
		x, okX := xMap[cid]
		y, okY := yMap[cid]
		if !okX || !okY {
			continue
		}

		// Apply category filter
		catVal := ""
		if catMap != nil {
			catVal = catMap[cid]
			if catFilterSet != nil {
				if _, ok := catFilterSet[catVal]; !ok {
					continue
				}
			}
		}

		cell := CellInfo{
			JoinID: cid,
			X:      x,
			Y:      y,
		}
		if exprMap != nil {
			cell.Expression = exprMap[cid] // defaults to 0 if not in map
		}
		if catMap != nil {
			cell.Category = catVal
		}

		cells = append(cells, cell)
	}

	return &CellQueryResult{
		Cells:      cells,
		TotalCount: len(cells),
		Truncated:  truncated,
	}, nil
}
