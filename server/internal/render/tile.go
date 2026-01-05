// Package render provides tile rendering using fogleman/gg.
package render

import (
	"bytes"
	"image"
	"image/color"
	"image/png"
	"sync"

	"github.com/fogleman/gg"
	"github.com/soma-tiles/server/internal/data/zarr"
	"github.com/soma-tiles/server/pkg/colormap"
)

// Config contains renderer configuration.
type Config struct {
	TileSize        int
	DefaultColormap string
}

// TileRenderer renders tiles from bin data.
type TileRenderer struct {
	config      Config
	contextPool sync.Pool
	bufferPool  sync.Pool
	colormaps   map[string]colormap.Colormap
}

// NewTileRenderer creates a new tile renderer.
func NewTileRenderer(cfg Config) *TileRenderer {
	r := &TileRenderer{
		config: cfg,
		contextPool: sync.Pool{
			New: func() interface{} {
				return gg.NewContext(cfg.TileSize, cfg.TileSize)
			},
		},
		bufferPool: sync.Pool{
			New: func() interface{} {
				return bytes.NewBuffer(make([]byte, 0, 32*1024))
			},
		},
		colormaps: make(map[string]colormap.Colormap),
	}

	// Initialize colormaps
	r.colormaps["viridis"] = colormap.Viridis
	r.colormaps["plasma"] = colormap.Plasma
	r.colormaps["inferno"] = colormap.Inferno
	r.colormaps["magma"] = colormap.Magma
	r.colormaps["seurat"] = colormap.Seurat
	r.colormaps["categorical"] = colormap.Categorical

	return r
}

// RenderTile renders a tile from bins using category coloring.
func (r *TileRenderer) RenderTile(bins []zarr.Bin, binsPerTileAxis int, tileX, tileY int) ([]byte, error) {
	// Get context from pool
	dc := r.contextPool.Get().(*gg.Context)
	defer r.contextPool.Put(dc)

	// Clear canvas with white background
	dc.SetColor(color.White)
	dc.Clear()

	if len(bins) == 0 {
		return r.encodeContext(dc)
	}

	// Calculate rendering parameters
	tileSize := float64(r.config.TileSize)
	if binsPerTileAxis <= 0 {
		binsPerTileAxis = 1
	}
	nBinsPerTile := float64(binsPerTileAxis)
	binSize := tileSize / nBinsPerTile

	// Render each bin
	for _, bin := range bins {
		// Calculate pixel position within tile
		// Bin coords are global, need to convert to tile-local
		localX := float64(bin.X) - float64(tileX)*nBinsPerTile
		localY := float64(bin.Y) - float64(tileY)*nBinsPerTile

		px := localX * binSize
		py := localY * binSize

		// Skip if outside tile
		if px < 0 || px >= tileSize || py < 0 || py >= tileSize {
			continue
		}

		// Color based on cell count (for now)
		intensity := float64(bin.CellCount) / 1000.0 // Normalize
		if intensity > 1 {
			intensity = 1
		}

		c := r.colormaps["viridis"].At(intensity)
		dc.SetColor(c)

		// Fill the entire bin area (tile-like rendering).
		dc.DrawRectangle(px, py, binSize, binSize)
		dc.Fill()
	}

	return r.encodeContext(dc)
}

// RenderExpressionTile renders a tile colored by gene expression.
func (r *TileRenderer) RenderExpressionTile(
	bins []zarr.Bin,
	expression []float32,
	exprMin float32,
	exprMax float32,
	binsPerTileAxis int,
	tileX, tileY int,
	colormapName string,
) ([]byte, error) {
	// Get context from pool
	dc := r.contextPool.Get().(*gg.Context)
	defer r.contextPool.Put(dc)

	// Clear canvas
	dc.SetColor(color.White)
	dc.Clear()

	if len(bins) == 0 || len(expression) == 0 {
		return r.encodeContext(dc)
	}

	// Get colormap
	cmap, ok := r.colormaps[colormapName]
	if !ok {
		cmap = r.colormaps[r.config.DefaultColormap]
	}

	exprRange := exprMax - exprMin
	if exprRange == 0 {
		exprRange = 1
	}

	// Calculate rendering parameters
	tileSize := float64(r.config.TileSize)
	if binsPerTileAxis <= 0 {
		binsPerTileAxis = 1
	}
	nBinsPerTile := float64(binsPerTileAxis)
	binSize := tileSize / nBinsPerTile

	// Render each bin
	for i, bin := range bins {
		if i >= len(expression) {
			break
		}

		// Calculate pixel position
		localX := float64(bin.X) - float64(tileX)*nBinsPerTile
		localY := float64(bin.Y) - float64(tileY)*nBinsPerTile

		px := localX * binSize
		py := localY * binSize

		if px < 0 || px >= tileSize || py < 0 || py >= tileSize {
			continue
		}

		// Normalize expression value
		normalized := float64(expression[i]-exprMin) / float64(exprRange)

		c := cmap.At(normalized)
		dc.SetColor(c)

		// Fill the entire bin area (tile-like rendering).
		dc.DrawRectangle(px, py, binSize, binSize)
		dc.Fill()
	}

	return r.encodeContext(dc)
}

// RenderCategoryTile renders a tile colored by category.
func (r *TileRenderer) RenderCategoryTile(
	bins []zarr.Bin,
	categoryIdx []int,
	binsPerTileAxis int,
	tileX, tileY int,
) ([]byte, error) {
	dc := r.contextPool.Get().(*gg.Context)
	defer r.contextPool.Put(dc)

	dc.SetColor(color.White)
	dc.Clear()

	if len(bins) == 0 {
		return r.encodeContext(dc)
	}

	cmap := r.colormaps["categorical"]

	tileSize := float64(r.config.TileSize)
	if binsPerTileAxis <= 0 {
		binsPerTileAxis = 1
	}
	nBinsPerTile := float64(binsPerTileAxis)
	binSize := tileSize / nBinsPerTile

	for i, bin := range bins {
		localX := float64(bin.X) - float64(tileX)*nBinsPerTile
		localY := float64(bin.Y) - float64(tileY)*nBinsPerTile

		px := localX * binSize
		py := localY * binSize

		if px < 0 || px >= tileSize || py < 0 || py >= tileSize {
			continue
		}

		catIdx := 0
		if i < len(categoryIdx) {
			catIdx = categoryIdx[i]
		}

		// Skip rendering if category is filtered out (-1)
		if catIdx < 0 {
			continue
		}

		c := cmap.AtIndex(catIdx)
		dc.SetColor(c)

		// Fill the entire bin area (tile-like rendering).
		dc.DrawRectangle(px, py, binSize, binSize)
		dc.Fill()
	}

	return r.encodeContext(dc)
}

func (r *TileRenderer) encodeContext(dc *gg.Context) ([]byte, error) {
	buf := r.bufferPool.Get().(*bytes.Buffer)
	defer func() {
		buf.Reset()
		r.bufferPool.Put(buf)
	}()

	// Use fast PNG encoder
	encoder := png.Encoder{CompressionLevel: png.BestSpeed}
	if err := encoder.Encode(buf, dc.Image()); err != nil {
		return nil, err
	}

	// Copy buffer contents (buffer will be reused)
	result := make([]byte, buf.Len())
	copy(result, buf.Bytes())
	return result, nil
}

// CreateEmptyTile creates an empty transparent tile.
func (r *TileRenderer) CreateEmptyTile() ([]byte, error) {
	img := image.NewRGBA(image.Rect(0, 0, r.config.TileSize, r.config.TileSize))
	// Fill with transparent white
	for i := 0; i < len(img.Pix); i += 4 {
		img.Pix[i] = 255   // R
		img.Pix[i+1] = 255 // G
		img.Pix[i+2] = 255 // B
		img.Pix[i+3] = 0   // A (transparent)
	}

	buf := bytes.NewBuffer(nil)
	if err := png.Encode(buf, img); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}
