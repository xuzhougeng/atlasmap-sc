// Package service provides business logic for the tile server.
package service

import (
	"fmt"
	"sync"

	"github.com/soma-tiles/server/internal/cache"
	"github.com/soma-tiles/server/internal/data/soma"
	"github.com/soma-tiles/server/internal/data/zarr"
	"github.com/soma-tiles/server/internal/render"
)

// TileServiceConfig contains tile service configuration.
type TileServiceConfig struct {
	DatasetID  string
	ZarrReader *zarr.Reader
	SomaReader *soma.Reader
	Cache      *cache.Manager
	Renderer   *render.TileRenderer
}

// TileService handles tile rendering and serving.
type TileService struct {
	datasetID string
	zarr      *zarr.Reader
	soma      *soma.Reader
	cache     *cache.Manager
	renderer  *render.TileRenderer

	renderZoom int

	binsOnce sync.Once
	bins     []zarr.Bin
	binsErr  error

	exprMu    sync.Mutex
	exprCache map[string]expressionCacheEntry

	catMu    sync.Mutex
	catCache map[string][]int

	geneCategoryMeansMu    sync.Mutex
	geneCategoryMeansCache map[string][]GeneCategoryMeanItem

	categoryCentroidsMu    sync.Mutex
	categoryCentroidsCache map[string][]CategoryCentroidItem
}

// NewTileService creates a new tile service.
func NewTileService(cfg TileServiceConfig) *TileService {
	renderZoom := 0
	if md := cfg.ZarrReader.Metadata(); md != nil && md.ZoomLevels > 0 {
		renderZoom = md.ZoomLevels - 1
	}

	datasetID := cfg.DatasetID
	if datasetID == "" {
		datasetID = "default"
	}

	return &TileService{
		datasetID:              datasetID,
		zarr:                   cfg.ZarrReader,
		soma:                   cfg.SomaReader,
		cache:                  cfg.Cache,
		renderer:               cfg.Renderer,
		renderZoom:             renderZoom,
		exprCache:              make(map[string]expressionCacheEntry),
		catCache:               make(map[string][]int),
		geneCategoryMeansCache: make(map[string][]GeneCategoryMeanItem),
		categoryCentroidsCache: make(map[string][]CategoryCentroidItem),
	}
}

type expressionCacheEntry struct {
	values []float32
	min    float32
	max    float32
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

func (s *TileService) binsForTile(mapZoom, tileX, tileY int) ([]zarr.Bin, []int, int, int, error) {
	if mapZoom < 0 || mapZoom > s.renderZoom {
		return nil, nil, 0, 0, fmt.Errorf("invalid zoom level: %d", mapZoom)
	}
	if tileX < 0 || tileY < 0 {
		return nil, nil, 0, 0, fmt.Errorf("invalid tile coordinates: %d/%d", tileX, tileY)
	}

	if err := s.loadBins(); err != nil {
		return nil, nil, 0, 0, err
	}

	tilesPerAxis := 1 << mapZoom
	if tileX >= tilesPerAxis || tileY >= tilesPerAxis {
		return nil, nil, 0, 0, fmt.Errorf("tile out of range: %d/%d (tiles_per_axis=%d)", tileX, tileY, tilesPerAxis)
	}

	binsPerTileAxis := 1 << (s.renderZoom - mapZoom)
	startX := int32(tileX * binsPerTileAxis)
	endX := int32((tileX + 1) * binsPerTileAxis)
	// Flip Y so that tileY=0 corresponds to the top (max-Y) of the dataset,
	// matching the conventional Cartesian orientation for UMAP plots.
	flippedTileY := (tilesPerAxis - 1) - tileY
	startY := int32(flippedTileY * binsPerTileAxis)
	endY := int32((flippedTileY + 1) * binsPerTileAxis)

	outBins := make([]zarr.Bin, 0, 64)
	outIdx := make([]int, 0, 64)
	for i, bin := range s.bins {
		if bin.X < startX || bin.X >= endX {
			continue
		}
		if bin.Y < startY || bin.Y >= endY {
			continue
		}
		outBins = append(outBins, bin)
		outIdx = append(outIdx, i)
	}

	return outBins, outIdx, binsPerTileAxis, flippedTileY, nil
}

func (s *TileService) expressionForGene(gene string) (expressionCacheEntry, error) {
	s.exprMu.Lock()
	defer s.exprMu.Unlock()

	if cached, ok := s.exprCache[gene]; ok {
		return cached, nil
	}

	geneIdx, ok := s.zarr.Metadata().GeneIndex[gene]
	if !ok {
		return expressionCacheEntry{}, fmt.Errorf("gene not found: %s", gene)
	}

	values, err := s.zarr.GetExpression(s.renderZoom, geneIdx, "mean")
	if err != nil {
		return expressionCacheEntry{}, err
	}

	if len(values) == 0 {
		cached := expressionCacheEntry{values: values, min: 0, max: 0}
		s.exprCache[gene] = cached
		return cached, nil
	}

	minV := values[0]
	maxV := values[0]
	for _, v := range values[1:] {
		if v < minV {
			minV = v
		}
		if v > maxV {
			maxV = v
		}
	}

	cached := expressionCacheEntry{values: values, min: minV, max: maxV}
	s.exprCache[gene] = cached
	return cached, nil
}

func (s *TileService) categoryForColumn(column string) ([]int, error) {
	s.catMu.Lock()
	defer s.catMu.Unlock()

	if cached, ok := s.catCache[column]; ok {
		return cached, nil
	}

	values, err := s.zarr.GetCategoryDataForBins(s.renderZoom, column)
	if err != nil {
		return nil, err
	}
	s.catCache[column] = values
	return values, nil
}

// GetTile returns a rendered tile PNG.
func (s *TileService) GetTile(z, x, y int) ([]byte, error) {
	// Check cache (prefix with dataset ID)
	cacheKey := s.datasetID + ":" + cache.TileKey(z, x, y, nil)
	if data, ok := s.cache.GetTile(cacheKey); ok {
		return data, nil
	}

	// Load bins for this tile (rendered from highest-available zoom level)
	bins, _, binsPerTileAxis, renderTileY, err := s.binsForTile(z, x, y)
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
	cacheKey := s.datasetID + ":" + cache.ExpressionTileKey(z, x, y, gene, colormap, cacheMin, cacheMax)
	if data, ok := s.cache.GetTile(cacheKey); ok {
		return data, nil
	}

	// Load bins (rendered from highest-available zoom level)
	bins, indices, binsPerTileAxis, renderTileY, err := s.binsForTile(z, x, y)
	if err != nil {
		return nil, fmt.Errorf("failed to load bins: %w", err)
	}

	exprAll, err := s.expressionForGene(gene)
	if err != nil {
		return nil, fmt.Errorf("failed to load expression: %w", err)
	}

	exprTile := make([]float32, len(indices))
	for i, idx := range indices {
		if idx >= 0 && idx < len(exprAll.values) {
			exprTile[i] = exprAll.values[idx]
		}
	}

	minV := exprAll.min
	maxV := exprAll.max
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
	return s.zarr.Metadata()
}

func (s *TileService) Soma() *soma.Reader {
	return s.soma
}

// GetCategoryTile returns a tile colored by category.
func (s *TileService) GetCategoryTile(z, x, y int, column string, categoryFilter []string) ([]byte, error) {
	// Check cache (prefix with dataset ID)
	cacheKey := s.datasetID + ":" + cache.CategoryTileKey(z, x, y, column, categoryFilter)
	if data, ok := s.cache.GetTile(cacheKey); ok {
		return data, nil
	}

	// Load bins (rendered from highest-available zoom level)
	bins, indices, binsPerTileAxis, renderTileY, err := s.binsForTile(z, x, y)
	if err != nil {
		return nil, fmt.Errorf("failed to load bins: %w", err)
	}

	// Load category data
	categoryAll, err := s.categoryForColumn(column)
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
	data, err := s.renderer.RenderCategoryTile(bins, categoryTile, binsPerTileAxis, x, renderTileY)
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
