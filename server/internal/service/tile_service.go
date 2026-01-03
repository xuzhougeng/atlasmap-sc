// Package service provides business logic for the tile server.
package service

import (
	"fmt"

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
}

// NewTileService creates a new tile service.
func NewTileService(cfg TileServiceConfig) *TileService {
	return &TileService{
		zarr:     cfg.ZarrReader,
		cache:    cfg.Cache,
		renderer: cfg.Renderer,
	}
}

// GetTile returns a rendered tile PNG.
func (s *TileService) GetTile(z, x, y int) ([]byte, error) {
	// Check cache
	cacheKey := cache.TileKey(z, x, y, nil)
	if data, ok := s.cache.GetTile(cacheKey); ok {
		return data, nil
	}

	// Load bins for this tile
	bins, err := s.zarr.GetBins(z, x, y)
	if err != nil {
		return nil, fmt.Errorf("failed to load bins: %w", err)
	}

	// Render tile
	data, err := s.renderer.RenderTile(bins, z, x, y)
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

	// Get gene index
	geneIdx, ok := s.zarr.Metadata().GeneIndex[gene]
	if !ok {
		return nil, fmt.Errorf("gene not found: %s", gene)
	}

	// Load bins
	bins, err := s.zarr.GetBins(z, x, y)
	if err != nil {
		return nil, fmt.Errorf("failed to load bins: %w", err)
	}

	// Load expression
	expression, err := s.zarr.GetExpression(z, geneIdx, "mean")
	if err != nil {
		return nil, fmt.Errorf("failed to load expression: %w", err)
	}

	// Render tile
	data, err := s.renderer.RenderExpressionTile(bins, expression, z, x, y, colormap)
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

	// Load bins for this tile
	bins, err := s.zarr.GetBins(z, x, y)
	if err != nil {
		return nil, fmt.Errorf("failed to load bins: %w", err)
	}

	// Load category data
	categoryIdx, err := s.zarr.GetCategoryDataForBins(z, column)
	if err != nil {
		return nil, fmt.Errorf("failed to load category data: %w", err)
	}

	// Render tile
	data, err := s.renderer.RenderCategoryTile(bins, categoryIdx, z, x, y)
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
