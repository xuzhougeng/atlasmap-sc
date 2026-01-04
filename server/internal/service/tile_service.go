// Package service provides business logic for the tile server.
package service

import (
	"fmt"
	"sync"

	"github.com/soma-tiles/server/internal/cache"
	"github.com/soma-tiles/server/internal/data/zarr"
	"github.com/soma-tiles/server/internal/render"
)

// TileServiceConfig contains tile service configuration.
type TileServiceConfig struct {
	ZarrReader *zarr.Reader
	Cache      *cache.Manager
	Renderer   *render.TileRenderer
}

// TileService handles tile rendering and serving.
type TileService struct {
	zarr     *zarr.Reader
	cache    *cache.Manager
	renderer *render.TileRenderer

	renderZoom int

	binsOnce sync.Once
	bins     []zarr.Bin
	binsErr  error

	exprMu    sync.Mutex
	exprCache map[string]expressionCacheEntry

	catMu    sync.Mutex
	catCache map[string][]int
}

// NewTileService creates a new tile service.
func NewTileService(cfg TileServiceConfig) *TileService {
	renderZoom := 0
	if md := cfg.ZarrReader.Metadata(); md != nil && md.ZoomLevels > 0 {
		renderZoom = md.ZoomLevels - 1
	}

	return &TileService{
		zarr:       cfg.ZarrReader,
		cache:      cfg.Cache,
		renderer:   cfg.Renderer,
		renderZoom: renderZoom,
		exprCache:  make(map[string]expressionCacheEntry),
		catCache:   make(map[string][]int),
	}
}

type expressionCacheEntry struct {
	values []float32
	min    float32
	max    float32
}

func (s *TileService) loadBins() error {
	s.binsOnce.Do(func() {
		s.bins, s.binsErr = s.zarr.GetBins(s.renderZoom, 0, 0)
	})
	return s.binsErr
}

func (s *TileService) binsForTile(mapZoom, tileX, tileY int) ([]zarr.Bin, []int, int, error) {
	if mapZoom < 0 || mapZoom > s.renderZoom {
		return nil, nil, 0, fmt.Errorf("invalid zoom level: %d", mapZoom)
	}
	if tileX < 0 || tileY < 0 {
		return nil, nil, 0, fmt.Errorf("invalid tile coordinates: %d/%d", tileX, tileY)
	}

	if err := s.loadBins(); err != nil {
		return nil, nil, 0, err
	}

	tilesPerAxis := 1 << mapZoom
	if tileX >= tilesPerAxis || tileY >= tilesPerAxis {
		return nil, nil, 0, fmt.Errorf("tile out of range: %d/%d (tiles_per_axis=%d)", tileX, tileY, tilesPerAxis)
	}

	binsPerTileAxis := 1 << (s.renderZoom - mapZoom)
	startX := int32(tileX * binsPerTileAxis)
	endX := int32((tileX + 1) * binsPerTileAxis)
	startY := int32(tileY * binsPerTileAxis)
	endY := int32((tileY + 1) * binsPerTileAxis)

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

	return outBins, outIdx, binsPerTileAxis, nil
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
	// Check cache
	cacheKey := cache.TileKey(z, x, y, nil)
	if data, ok := s.cache.GetTile(cacheKey); ok {
		return data, nil
	}

	// Load bins for this tile (rendered from highest-available zoom level)
	bins, _, binsPerTileAxis, err := s.binsForTile(z, x, y)
	if err != nil {
		return nil, fmt.Errorf("failed to load bins: %w", err)
	}

	// Render tile
	data, err := s.renderer.RenderTile(bins, binsPerTileAxis, x, y)
	if err != nil {
		return nil, fmt.Errorf("failed to render tile: %w", err)
	}

	// Cache result
	s.cache.SetTile(cacheKey, data)

	return data, nil
}

// GetExpressionTile returns a tile colored by gene expression.
func (s *TileService) GetExpressionTile(z, x, y int, gene string, colormap string) ([]byte, error) {
	// Check cache
	cacheKey := cache.ExpressionTileKey(z, x, y, gene, colormap)
	if data, ok := s.cache.GetTile(cacheKey); ok {
		return data, nil
	}

	// Load bins (rendered from highest-available zoom level)
	bins, indices, binsPerTileAxis, err := s.binsForTile(z, x, y)
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

	// Render tile
	data, err := s.renderer.RenderExpressionTile(
		bins,
		exprTile,
		exprAll.min,
		exprAll.max,
		binsPerTileAxis,
		x, y,
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

// GetCategoryTile returns a tile colored by category.
func (s *TileService) GetCategoryTile(z, x, y int, column string) ([]byte, error) {
	// Check cache
	cacheKey := cache.CategoryTileKey(z, x, y, column)
	if data, ok := s.cache.GetTile(cacheKey); ok {
		return data, nil
	}

	// Load bins (rendered from highest-available zoom level)
	bins, indices, binsPerTileAxis, err := s.binsForTile(z, x, y)
	if err != nil {
		return nil, fmt.Errorf("failed to load bins: %w", err)
	}

	// Load category data
	categoryAll, err := s.categoryForColumn(column)
	if err != nil {
		return nil, fmt.Errorf("failed to load category data: %w", err)
	}

	categoryTile := make([]int, len(indices))
	for i, idx := range indices {
		if idx >= 0 && idx < len(categoryAll) {
			categoryTile[i] = categoryAll[idx]
		}
	}

	// Render tile
	data, err := s.renderer.RenderCategoryTile(bins, categoryTile, binsPerTileAxis, x, y)
	if err != nil {
		return nil, fmt.Errorf("failed to render tile: %w", err)
	}

	// Cache result
	s.cache.SetTile(cacheKey, data)

	return data, nil
}

// CategoryLegendItem represents a legend item for a category.
type CategoryLegendItem struct {
	Value string `json:"value"`
	Color string `json:"color"`
	Index int    `json:"index"`
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

	legend := make([]CategoryLegendItem, len(catInfo.Values))
	for i, value := range catInfo.Values {
		legend[i] = CategoryLegendItem{
			Value: value,
			Color: colors[value],
			Index: i,
		}
	}

	return legend, nil
}

// GetCategoryColors returns the color mapping for a category column.
func (s *TileService) GetCategoryColors(column string) (map[string]string, error) {
	return s.zarr.GetCategoryColors(column)
}

// QueryBinsExpressingGene returns bins where gene expression is above threshold.
func (s *TileService) QueryBinsExpressingGene(
	gene string,
	threshold float32,
	offset, limit int,
) (*zarr.BinQueryResult, error) {
	return s.zarr.GetBinsExpressingGene(gene, threshold, offset, limit)
}

// GetGeneStats returns statistics for a gene.
func (s *TileService) GetGeneStats(gene string) (*zarr.GeneStats, error) {
	return s.zarr.GetGeneStats(gene)
}
