// Package zarr provides a reader for Zarr v3 stores.
package zarr

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"sync"
	"unsafe"

	"github.com/klauspost/compress/zstd"
)

// Reader provides access to Zarr bin data.
type Reader struct {
	basePath string
	metadata *ZarrMetadata
	mu       sync.RWMutex
	decoder  *zstd.Decoder

	// Cached data per zoom level
	zoomData map[int]*ZoomLevelData
}

// ZarrMetadata contains metadata about the Zarr store.
type ZarrMetadata struct {
	ZoomLevels    int                     `json:"zoom_levels"`
	TileSize      int                     `json:"tile_size"`
	NGenes        int                     `json:"n_genes_preaggregated"`
	NCategories   map[string]int          `json:"n_categories"`
	FormatVersion string                  `json:"format_version"`
	DatasetName   string                  `json:"dataset_name"`
	NCells        int                     `json:"n_cells"`
	Genes         []string                `json:"preaggregated_genes"`
	GeneIndex     map[string]int          `json:"gene_index"`
	Categories    map[string]CategoryInfo `json:"categories"`
	Bounds        Bounds                  `json:"bounds"`
	// Optional metadata for multiple coordinate systems (preprocessed as multiple Zarr stores).
	CoordinateSystems       []CoordinateSystem `json:"coordinate_systems,omitempty"`
	DefaultCoordinateSystem string             `json:"default_coordinate_system,omitempty"`
	// Backward-compatible field written by older preprocessors.
	UMAPKey string `json:"umap_key,omitempty"`
	// Filled by the server when serving a specific coordinate system.
	CoordinateSystem string `json:"coordinate_system,omitempty"`
}

// CoordinateSystem describes one coordinate system's Zarr store location.
type CoordinateSystem struct {
	Key      string `json:"key"`
	ZarrPath string `json:"zarr_path"`
}

// CategoryInfo contains information about a category column.
type CategoryInfo struct {
	Values  []string       `json:"values"`
	Mapping map[string]int `json:"mapping"`
}

// Bounds represents coordinate bounds.
type Bounds struct {
	MinX float64 `json:"min_x"`
	MaxX float64 `json:"max_x"`
	MinY float64 `json:"min_y"`
	MaxY float64 `json:"max_y"`
}

// ZoomLevelData contains cached data for a zoom level.
type ZoomLevelData struct {
	NBins        int
	NBinsPerAxis int
	BinSize      float64
	BinCoords    [][2]int32
	CellCounts   []uint32
	// Expression data is loaded on demand
}

// Bin represents a single spatial bin.
type Bin struct {
	X              int32
	Y              int32
	CellCount      uint32
	ExpressionMean []float32
	ExpressionMax  []float32
	CategoryCounts map[string][]uint32
}

// ZarrV3ArrayMeta represents Zarr v3 array metadata (zarr.json).
type ZarrV3ArrayMeta struct {
	Shape     []int  `json:"shape"`
	DataType  string `json:"data_type"`
	ChunkGrid struct {
		Name          string `json:"name"`
		Configuration struct {
			ChunkShape []int `json:"chunk_shape"`
		} `json:"configuration"`
	} `json:"chunk_grid"`
	ChunkKeyEncoding struct {
		Name          string `json:"name"`
		Configuration struct {
			Separator string `json:"separator"`
		} `json:"configuration"`
	} `json:"chunk_key_encoding"`
	FillValue interface{} `json:"fill_value"`
	Codecs    []struct {
		Name          string                 `json:"name"`
		Configuration map[string]interface{} `json:"configuration"`
	} `json:"codecs"`
	ZarrFormat int    `json:"zarr_format"`
	NodeType   string `json:"node_type"`
}

// NewReader creates a new Zarr reader.
func NewReader(basePath string) (*Reader, error) {
	decoder, err := zstd.NewReader(nil)
	if err != nil {
		return nil, fmt.Errorf("failed to create zstd decoder: %w", err)
	}

	r := &Reader{
		basePath: basePath,
		zoomData: make(map[int]*ZoomLevelData),
		decoder:  decoder,
	}

	// Load metadata
	if err := r.loadMetadata(); err != nil {
		return nil, fmt.Errorf("failed to load metadata: %w", err)
	}

	// Preload zoom level metadata
	for z := 0; z < r.metadata.ZoomLevels; z++ {
		if err := r.loadZoomLevel(z); err != nil {
			return nil, fmt.Errorf("failed to load zoom level %d: %w", z, err)
		}
	}

	return r, nil
}

// Metadata returns the Zarr metadata.
func (r *Reader) Metadata() *ZarrMetadata {
	return r.metadata
}

func (r *Reader) loadMetadata() error {
	// Load main metadata from parent directory
	metadataPath := filepath.Join(r.basePath, "..", "metadata.json")
	data, err := os.ReadFile(metadataPath)
	if err != nil {
		return fmt.Errorf("failed to read metadata.json: %w", err)
	}

	var metadata ZarrMetadata
	if err := json.Unmarshal(data, &metadata); err != nil {
		return fmt.Errorf("failed to parse metadata.json: %w", err)
	}

	// Build gene index from gene list if not present
	if metadata.GeneIndex == nil {
		metadata.GeneIndex = make(map[string]int)
		for i, gene := range metadata.Genes {
			metadata.GeneIndex[gene] = i
		}
	}

	// Load gene index file if exists
	geneIndexPath := filepath.Join(r.basePath, "..", "gene_index.json")
	if data, err := os.ReadFile(geneIndexPath); err == nil {
		if err := json.Unmarshal(data, &metadata.GeneIndex); err != nil {
			return fmt.Errorf("failed to parse gene_index.json: %w", err)
		}
	}

	r.metadata = &metadata
	return nil
}

func (r *Reader) loadZoomLevel(zoom int) error {
	zoomPath := filepath.Join(r.basePath, fmt.Sprintf("zoom_%d", zoom))

	// Check if directory exists
	if _, err := os.Stat(zoomPath); os.IsNotExist(err) {
		return fmt.Errorf("zoom level %d not found", zoom)
	}

	nBinsAxis := 1 << zoom
	zd := &ZoomLevelData{
		NBinsPerAxis: nBinsAxis,
		BinSize:      256.0 / float64(nBinsAxis),
	}

	// Try to load cell_count to get bin count
	cellCountPath := filepath.Join(zoomPath, "cell_count")
	if meta, err := r.loadArrayMeta(cellCountPath); err == nil {
		if len(meta.Shape) > 0 {
			zd.NBins = meta.Shape[0]
		}
	}

	r.zoomData[zoom] = zd
	return nil
}

// loadArrayMeta loads Zarr v3 array metadata.
func (r *Reader) loadArrayMeta(arrayPath string) (*ZarrV3ArrayMeta, error) {
	metaPath := filepath.Join(arrayPath, "zarr.json")
	data, err := os.ReadFile(metaPath)
	if err != nil {
		return nil, err
	}

	var meta ZarrV3ArrayMeta
	if err := json.Unmarshal(data, &meta); err != nil {
		return nil, err
	}

	return &meta, nil
}

// readChunk reads and decompresses a chunk from Zarr v3 format.
func (r *Reader) readChunk(arrayPath string, chunkKey string) ([]byte, error) {
	// Zarr v3 stores chunks in c/ directory
	chunkPath := filepath.Join(arrayPath, "c", chunkKey)

	compressedData, err := os.ReadFile(chunkPath)
	if err != nil {
		return nil, err
	}

	// Decompress with zstd
	decompressed, err := r.decoder.DecodeAll(compressedData, nil)
	if err != nil {
		return nil, fmt.Errorf("zstd decompress failed: %w", err)
	}

	return decompressed, nil
}

func (r *Reader) encodeChunkKey(meta *ZarrV3ArrayMeta, chunkIndices []int) string {
	sep := meta.ChunkKeyEncoding.Configuration.Separator
	if sep == "" {
		sep = "/"
	}
	parts := make([]string, len(chunkIndices))
	for i, idx := range chunkIndices {
		parts[i] = strconv.Itoa(idx)
	}
	return strings.Join(parts, sep)
}

func (r *Reader) chunkShapeAt(meta *ZarrV3ArrayMeta, chunkIndices []int) ([]int, error) {
	if len(meta.Shape) == 0 || len(meta.ChunkGrid.Configuration.ChunkShape) == 0 {
		return nil, fmt.Errorf("invalid zarr metadata: missing shape/chunk_shape")
	}
	if len(meta.Shape) != len(meta.ChunkGrid.Configuration.ChunkShape) {
		return nil, fmt.Errorf("invalid zarr metadata: shape dims (%d) != chunk dims (%d)", len(meta.Shape), len(meta.ChunkGrid.Configuration.ChunkShape))
	}
	if len(chunkIndices) != len(meta.Shape) {
		return nil, fmt.Errorf("invalid chunk indices: got %d dims, expected %d", len(chunkIndices), len(meta.Shape))
	}

	actual := make([]int, len(meta.Shape))
	for d := range meta.Shape {
		chunkLen := meta.ChunkGrid.Configuration.ChunkShape[d]
		if chunkLen <= 0 {
			return nil, fmt.Errorf("invalid chunk shape at dim %d: %d", d, chunkLen)
		}
		start := chunkIndices[d] * chunkLen
		if start < 0 || start >= meta.Shape[d] {
			return nil, fmt.Errorf("chunk index out of range at dim %d: start=%d shape=%d", d, start, meta.Shape[d])
		}
		remaining := meta.Shape[d] - start
		if remaining < chunkLen {
			chunkLen = remaining
		}
		actual[d] = chunkLen
	}

	return actual, nil
}

func zarrDTypeSize(dataType string) (int, error) {
	switch dataType {
	case "float32", "int32", "uint32":
		return 4, nil
	case "uint64":
		return 8, nil
	default:
		return 0, fmt.Errorf("unsupported zarr data_type: %s", dataType)
	}
}

func zarrFillValueBytes(meta *ZarrV3ArrayMeta) ([]byte, error) {
	size, err := zarrDTypeSize(meta.DataType)
	if err != nil {
		return nil, err
	}

	// Default fill to 0 if unspecified.
	fill := meta.FillValue
	if fill == nil {
		return make([]byte, size), nil
	}

	switch meta.DataType {
	case "float32":
		var v float32
		switch t := fill.(type) {
		case float64:
			v = float32(t)
		case float32:
			v = t
		case int:
			v = float32(t)
		default:
			return nil, fmt.Errorf("unsupported fill_value type for float32: %T", fill)
		}
		bits := mathFloat32bits(v)
		return []byte{byte(bits), byte(bits >> 8), byte(bits >> 16), byte(bits >> 24)}, nil
	case "int32":
		var v int32
		switch t := fill.(type) {
		case float64:
			v = int32(t)
		case int:
			v = int32(t)
		case int32:
			v = t
		default:
			return nil, fmt.Errorf("unsupported fill_value type for int32: %T", fill)
		}
		u := uint32(v)
		return []byte{byte(u), byte(u >> 8), byte(u >> 16), byte(u >> 24)}, nil
	case "uint32":
		var v uint32
		switch t := fill.(type) {
		case float64:
			v = uint32(t)
		case int:
			v = uint32(t)
		case uint32:
			v = t
		default:
			return nil, fmt.Errorf("unsupported fill_value type for uint32: %T", fill)
		}
		return []byte{byte(v), byte(v >> 8), byte(v >> 16), byte(v >> 24)}, nil
	case "uint64":
		var v uint64
		switch t := fill.(type) {
		case float64:
			v = uint64(t)
		case int:
			v = uint64(t)
		case uint64:
			v = t
		default:
			return nil, fmt.Errorf("unsupported fill_value type for uint64: %T", fill)
		}
		return []byte{
			byte(v), byte(v >> 8), byte(v >> 16), byte(v >> 24),
			byte(v >> 32), byte(v >> 40), byte(v >> 48), byte(v >> 56),
		}, nil
	default:
		return nil, fmt.Errorf("unsupported zarr data_type: %s", meta.DataType)
	}
}

func repeatFillBytes(fill []byte, n int) []byte {
	if n <= 0 {
		return nil
	}
	if len(fill) == 0 {
		return make([]byte, n)
	}
	// Fast path: fill is all zeros; make() already zero-initializes.
	allZero := true
	for _, b := range fill {
		if b != 0 {
			allZero = false
			break
		}
	}
	if allZero {
		return make([]byte, len(fill)*n)
	}

	out := make([]byte, len(fill)*n)
	for i := 0; i < n; i++ {
		copy(out[i*len(fill):(i+1)*len(fill)], fill)
	}
	return out
}

func product(ints []int) int {
	p := 1
	for _, v := range ints {
		p *= v
	}
	return p
}

func (r *Reader) readChunkAt(arrayPath string, meta *ZarrV3ArrayMeta, chunkIndices []int) ([]byte, error) {
	key := r.encodeChunkKey(meta, chunkIndices)
	data, err := r.readChunk(arrayPath, key)
	if err == nil {
		return data, nil
	}

	// Backwards-compatible fallback: some writers may drop trailing singleton chunk dims
	// (e.g. store [N,2] chunks as c/<rowChunk> instead of c/<rowChunk>/0).
	var altErr error
	if len(chunkIndices) > 1 {
		trailingAllZero := true
		for _, v := range chunkIndices[1:] {
			if v != 0 {
				trailingAllZero = false
				break
			}
		}
		if trailingAllZero {
			altKey := strconv.Itoa(chunkIndices[0])
			altData, altReadErr := r.readChunk(arrayPath, altKey)
			if altReadErr == nil {
				return altData, nil
			}
			altErr = altReadErr
		}
	}

	// If the chunk is not present on disk, it represents an all-fill-value chunk.
	if os.IsNotExist(err) && (altErr == nil || os.IsNotExist(altErr)) {
		shape, shapeErr := r.chunkShapeAt(meta, chunkIndices)
		if shapeErr != nil {
			return nil, shapeErr
		}
		elementCount := product(shape)
		fillBytes, fillErr := zarrFillValueBytes(meta)
		if fillErr != nil {
			return nil, fillErr
		}
		return repeatFillBytes(fillBytes, elementCount), nil
	}

	return nil, err
}

// GetBins returns all bins for a tile at the specified zoom level.
func (r *Reader) GetBins(zoom, tileX, tileY int) ([]Bin, error) {
	r.mu.RLock()
	defer r.mu.RUnlock()

	if zoom < 0 || zoom >= r.metadata.ZoomLevels {
		return nil, fmt.Errorf("invalid zoom level: %d", zoom)
	}

	zoomPath := filepath.Join(r.basePath, fmt.Sprintf("zoom_%d", zoom))

	// Load bin coordinates (shape: [N, 2], type: int32)
	binCoordsPath := filepath.Join(zoomPath, "bin_coords")
	binCoordsMeta, err := r.loadArrayMeta(binCoordsPath)
	if err != nil {
		return nil, fmt.Errorf("failed to load bin_coords metadata: %w", err)
	}

	// Load cell counts (shape: [N], type: uint32)
	cellCountPath := filepath.Join(zoomPath, "cell_count")
	cellCountMeta, err := r.loadArrayMeta(cellCountPath)
	if err != nil {
		return nil, fmt.Errorf("failed to load cell_count metadata: %w", err)
	}

	// Parse cell counts across all chunks.
	if len(cellCountMeta.Shape) != 1 {
		return nil, fmt.Errorf("unexpected cell_count shape: %v", cellCountMeta.Shape)
	}
	if len(cellCountMeta.ChunkGrid.Configuration.ChunkShape) != 1 {
		return nil, fmt.Errorf("unexpected cell_count chunk shape: %v", cellCountMeta.ChunkGrid.Configuration.ChunkShape)
	}

	nBins := cellCountMeta.Shape[0]
	cellCounts := make([]uint32, nBins)
	cellChunk := cellCountMeta.ChunkGrid.Configuration.ChunkShape[0]
	nCellChunks := ceilDiv(nBins, cellChunk)

	for chunk := 0; chunk < nCellChunks; chunk++ {
		chunkStart := chunk * cellChunk
		chunkLen := min(cellChunk, nBins-chunkStart)

		chunkData, err := r.readChunkAt(cellCountPath, cellCountMeta, []int{chunk})
		if err != nil {
			return nil, fmt.Errorf("failed to load cell_count chunk %d: %w", chunk, err)
		}
		if len(chunkData) < chunkLen*4 {
			return nil, fmt.Errorf("cell_count chunk %d too short: got %d bytes, expected %d", chunk, len(chunkData), chunkLen*4)
		}

		for i := 0; i < chunkLen; i++ {
			off := i * 4
			cellCounts[chunkStart+i] = uint32(chunkData[off]) |
				uint32(chunkData[off+1])<<8 |
				uint32(chunkData[off+2])<<16 |
				uint32(chunkData[off+3])<<24
		}
	}

	// Parse bin coordinates across all chunks.
	if len(binCoordsMeta.Shape) != 2 || binCoordsMeta.Shape[0] != nBins || binCoordsMeta.Shape[1] != 2 {
		return nil, fmt.Errorf("unexpected bin_coords shape: %v (expected [%d,2])", binCoordsMeta.Shape, nBins)
	}
	if len(binCoordsMeta.ChunkGrid.Configuration.ChunkShape) != 2 {
		return nil, fmt.Errorf("unexpected bin_coords chunk shape: %v", binCoordsMeta.ChunkGrid.Configuration.ChunkShape)
	}

	binCoords := make([][2]int32, nBins)
	coordChunkRows := binCoordsMeta.ChunkGrid.Configuration.ChunkShape[0]
	coordChunkCols := binCoordsMeta.ChunkGrid.Configuration.ChunkShape[1]
	nCoordRowChunks := ceilDiv(nBins, coordChunkRows)
	nCoordColChunks := ceilDiv(2, coordChunkCols)

	for rowChunk := 0; rowChunk < nCoordRowChunks; rowChunk++ {
		rowStart := rowChunk * coordChunkRows
		rowLen := min(coordChunkRows, nBins-rowStart)

		for colChunk := 0; colChunk < nCoordColChunks; colChunk++ {
			colStart := colChunk * coordChunkCols
			colLen := min(coordChunkCols, 2-colStart)

			chunkData, err := r.readChunkAt(binCoordsPath, binCoordsMeta, []int{rowChunk, colChunk})
			if err != nil {
				return nil, fmt.Errorf("failed to load bin_coords chunk %d/%d: %w", rowChunk, colChunk, err)
			}
			if len(chunkData) < rowLen*colLen*4 {
				return nil, fmt.Errorf("bin_coords chunk %d/%d too short: got %d bytes, expected %d", rowChunk, colChunk, len(chunkData), rowLen*colLen*4)
			}

			for rIdx := 0; rIdx < rowLen; rIdx++ {
				for cIdx := 0; cIdx < colLen; cIdx++ {
					off := (rIdx*colLen + cIdx) * 4
					u := uint32(chunkData[off]) |
						uint32(chunkData[off+1])<<8 |
						uint32(chunkData[off+2])<<16 |
						uint32(chunkData[off+3])<<24
					v := int32(u)

					globalRow := rowStart + rIdx
					globalCol := colStart + cIdx
					if globalCol == 0 {
						binCoords[globalRow][0] = v
					} else if globalCol == 1 {
						binCoords[globalRow][1] = v
					}
				}
			}
		}
	}

	// Build bins
	bins := make([]Bin, nBins)
	for i := 0; i < nBins; i++ {
		bins[i] = Bin{
			X:         binCoords[i][0],
			Y:         binCoords[i][1],
			CellCount: cellCounts[i],
		}
	}

	return bins, nil
}

// GetBinsForTile returns bins that fall within a specific tile.
func (r *Reader) GetBinsForTile(zoom, tileX, tileY int) ([]Bin, error) {
	allBins, err := r.GetBins(zoom, tileX, tileY)
	if err != nil {
		return nil, err
	}

	// For now, return all bins since we're not tiling yet
	// In a real implementation, we'd filter bins by tile coordinates
	return allBins, nil
}

// GetExpression returns expression values for a gene at a zoom level.
func (r *Reader) GetExpression(zoom int, geneIdx int, aggregation string) ([]float32, error) {
	zoomPath := filepath.Join(r.basePath, fmt.Sprintf("zoom_%d", zoom))

	var arrayName string
	switch aggregation {
	case "mean":
		arrayName = "expression_mean"
	case "max":
		arrayName = "expression_max"
	default:
		arrayName = "expression_mean"
	}

	// Load expression array
	exprPath := filepath.Join(zoomPath, arrayName)
	exprMeta, err := r.loadArrayMeta(exprPath)
	if err != nil {
		return nil, fmt.Errorf("failed to load expression metadata: %w", err)
	}

	// Expression is stored as [N_bins, N_genes]
	if len(exprMeta.Shape) != 2 {
		return nil, fmt.Errorf("unexpected expression shape: %v", exprMeta.Shape)
	}
	nBins := exprMeta.Shape[0]
	nGenes := exprMeta.Shape[1]
	if geneIdx < 0 || geneIdx >= nGenes {
		return nil, fmt.Errorf("gene index out of range: %d (n_genes=%d)", geneIdx, nGenes)
	}
	if len(exprMeta.ChunkGrid.Configuration.ChunkShape) != 2 {
		return nil, fmt.Errorf("unexpected expression chunk shape: %v", exprMeta.ChunkGrid.Configuration.ChunkShape)
	}

	rowChunk := exprMeta.ChunkGrid.Configuration.ChunkShape[0]
	colChunk := exprMeta.ChunkGrid.Configuration.ChunkShape[1]
	if rowChunk <= 0 || colChunk <= 0 {
		return nil, fmt.Errorf("invalid expression chunk shape: %v", exprMeta.ChunkGrid.Configuration.ChunkShape)
	}

	nRowChunks := ceilDiv(nBins, rowChunk)
	nColChunks := ceilDiv(nGenes, colChunk)
	geneColChunk := geneIdx / colChunk
	if geneColChunk >= nColChunks {
		return nil, fmt.Errorf("gene chunk out of range: %d (n_col_chunks=%d)", geneColChunk, nColChunks)
	}
	geneOffset := geneIdx % colChunk

	geneExpr := make([]float32, nBins)
	for rChunk := 0; rChunk < nRowChunks; rChunk++ {
		rowStart := rChunk * rowChunk
		rowLen := min(rowChunk, nBins-rowStart)
		colStart := geneColChunk * colChunk
		colLen := min(colChunk, nGenes-colStart)

		chunkData, err := r.readChunkAt(exprPath, exprMeta, []int{rChunk, geneColChunk})
		if err != nil {
			return nil, fmt.Errorf("failed to load expression chunk %d/%d: %w", rChunk, geneColChunk, err)
		}
		if len(chunkData) < rowLen*colLen*4 {
			return nil, fmt.Errorf("expression chunk %d/%d too short: got %d bytes, expected %d", rChunk, geneColChunk, len(chunkData), rowLen*colLen*4)
		}
		if geneOffset >= colLen {
			return nil, fmt.Errorf("gene offset out of chunk range: geneOffset=%d colLen=%d", geneOffset, colLen)
		}

		for i := 0; i < rowLen; i++ {
			off := (i*colLen + geneOffset) * 4
			bits := uint32(chunkData[off]) |
				uint32(chunkData[off+1])<<8 |
				uint32(chunkData[off+2])<<16 |
				uint32(chunkData[off+3])<<24
			geneExpr[rowStart+i] = float32frombits(bits)
		}
	}

	return geneExpr, nil
}

// GetCategoryDataForBins returns category indices for bins at a zoom level.
// Returns the dominant category index for each bin (argmax of category counts).
func (r *Reader) GetCategoryDataForBins(zoom int, categoryColumn string) ([]int, error) {
	zoomPath := filepath.Join(r.basePath, fmt.Sprintf("zoom_%d", zoom))
	catPath := filepath.Join(zoomPath, "category_counts", categoryColumn)

	// Check if category column exists
	if _, err := os.Stat(catPath); os.IsNotExist(err) {
		return nil, fmt.Errorf("category column not found: %s", categoryColumn)
	}

	// Get category info
	catInfo, ok := r.metadata.Categories[categoryColumn]
	if !ok {
		return nil, fmt.Errorf("category info not found in metadata: %s", categoryColumn)
	}
	nCategories := len(catInfo.Values)
	if nCategories == 0 {
		return nil, fmt.Errorf("category has no values: %s", categoryColumn)
	}

	// Read category counts array [N_bins, N_categories] as uint32 across chunks.
	catMeta, err := r.loadArrayMeta(catPath)
	if err != nil {
		return nil, fmt.Errorf("failed to load category metadata: %w", err)
	}
	if len(catMeta.Shape) != 2 || catMeta.Shape[1] != nCategories {
		return nil, fmt.Errorf("unexpected category_counts shape: %v (expected [N,%d])", catMeta.Shape, nCategories)
	}
	if len(catMeta.ChunkGrid.Configuration.ChunkShape) != 2 {
		return nil, fmt.Errorf("unexpected category_counts chunk shape: %v", catMeta.ChunkGrid.Configuration.ChunkShape)
	}

	nBins := catMeta.Shape[0]
	rowChunk := catMeta.ChunkGrid.Configuration.ChunkShape[0]
	colChunk := catMeta.ChunkGrid.Configuration.ChunkShape[1]
	if rowChunk <= 0 || colChunk <= 0 {
		return nil, fmt.Errorf("invalid category_counts chunk shape: %v", catMeta.ChunkGrid.Configuration.ChunkShape)
	}

	nRowChunks := ceilDiv(nBins, rowChunk)
	nColChunks := ceilDiv(nCategories, colChunk)

	dominantIdx := make([]int, nBins)
	maxCounts := make([]uint32, nBins)
	for rChunk := 0; rChunk < nRowChunks; rChunk++ {
		rowStart := rChunk * rowChunk
		rowLen := min(rowChunk, nBins-rowStart)

		for cChunk := 0; cChunk < nColChunks; cChunk++ {
			colStart := cChunk * colChunk
			colLen := min(colChunk, nCategories-colStart)

			chunkData, err := r.readChunkAt(catPath, catMeta, []int{rChunk, cChunk})
			if err != nil {
				return nil, fmt.Errorf("failed to load category_counts chunk %d/%d: %w", rChunk, cChunk, err)
			}
			if len(chunkData) < rowLen*colLen*4 {
				return nil, fmt.Errorf("category_counts chunk %d/%d too short: got %d bytes, expected %d", rChunk, cChunk, len(chunkData), rowLen*colLen*4)
			}

			for i := 0; i < rowLen; i++ {
				globalRow := rowStart + i
				maxCount := maxCounts[globalRow]
				maxIdx := dominantIdx[globalRow]

				for j := 0; j < colLen; j++ {
					off := (i*colLen + j) * 4
					count := uint32(chunkData[off]) |
						uint32(chunkData[off+1])<<8 |
						uint32(chunkData[off+2])<<16 |
						uint32(chunkData[off+3])<<24
					if count > maxCount {
						maxCount = count
						maxIdx = colStart + j
					}
				}

				dominantIdx[globalRow] = maxIdx
				maxCounts[globalRow] = maxCount
			}
		}
	}

	return dominantIdx, nil
}

// GetCategoryColors returns the color mapping for a category column.
// Returns a map of category value -> hex color string.
func (r *Reader) GetCategoryColors(categoryColumn string) (map[string]string, error) {
	catInfo, ok := r.metadata.Categories[categoryColumn]
	if !ok {
		return nil, fmt.Errorf("category not found: %s", categoryColumn)
	}

	// 20-color categorical palette (matches frontend)
	colors := []string{
		"#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
		"#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
		"#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
		"#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
	}

	result := make(map[string]string)
	for i, value := range catInfo.Values {
		result[value] = colors[i%len(colors)]
	}

	return result, nil
}

// GetAvailableCategories returns all available category column names.
func (r *Reader) GetAvailableCategories() []string {
	categories := make([]string, 0, len(r.metadata.Categories))
	for col := range r.metadata.Categories {
		categories = append(categories, col)
	}
	return categories
}

// BinExpressionInfo contains bin information with expression data.
type BinExpressionInfo struct {
	BinIndex   int               `json:"bin_index"`
	BinX       int32             `json:"bin_x"`
	BinY       int32             `json:"bin_y"`
	CellCount  uint32            `json:"cell_count"`
	Expression float32           `json:"expression"`
	Categories map[string]string `json:"categories,omitempty"`
}

// BinQueryResult contains paginated bin query results.
type BinQueryResult struct {
	Bins       []BinExpressionInfo `json:"bins"`
	TotalCount int                 `json:"total_count"`
	Offset     int                 `json:"offset"`
	Limit      int                 `json:"limit"`
	Gene       string              `json:"gene"`
	Threshold  float32             `json:"threshold"`
}

// GetBinsExpressingGene returns bins where gene expression is above threshold.
func (r *Reader) GetBinsExpressingGene(
	gene string,
	threshold float32,
	offset, limit int,
) (*BinQueryResult, error) {
	// Get gene index
	geneIdx, ok := r.metadata.GeneIndex[gene]
	if !ok {
		return nil, fmt.Errorf("gene not found: %s", gene)
	}

	// Use zoom level 0 for highest resolution
	zoom := 0

	// Load bins
	bins, err := r.GetBins(zoom, 0, 0)
	if err != nil {
		return nil, fmt.Errorf("failed to load bins: %w", err)
	}

	// Load expression
	expression, err := r.GetExpression(zoom, geneIdx, "mean")
	if err != nil {
		return nil, fmt.Errorf("failed to load expression: %w", err)
	}

	// Load category data for all categories
	categoryData := make(map[string][]int)
	for catName := range r.metadata.Categories {
		catIdx, err := r.GetCategoryDataForBins(zoom, catName)
		if err == nil {
			categoryData[catName] = catIdx
		}
	}

	// Filter bins by threshold
	var matchingBins []BinExpressionInfo
	for i, bin := range bins {
		if i >= len(expression) {
			break
		}

		exprValue := expression[i]
		if exprValue < threshold {
			continue
		}

		binInfo := BinExpressionInfo{
			BinIndex:   i,
			BinX:       bin.X,
			BinY:       bin.Y,
			CellCount:  bin.CellCount,
			Expression: exprValue,
			Categories: make(map[string]string),
		}

		// Add category info
		for catName, catIdx := range categoryData {
			if i < len(catIdx) {
				catInfo := r.metadata.Categories[catName]
				if catIdx[i] < len(catInfo.Values) {
					binInfo.Categories[catName] = catInfo.Values[catIdx[i]]
				}
			}
		}

		matchingBins = append(matchingBins, binInfo)
	}

	// Apply pagination
	totalCount := len(matchingBins)
	if offset >= totalCount {
		matchingBins = []BinExpressionInfo{}
	} else {
		end := offset + limit
		if end > totalCount {
			end = totalCount
		}
		matchingBins = matchingBins[offset:end]
	}

	return &BinQueryResult{
		Bins:       matchingBins,
		TotalCount: totalCount,
		Offset:     offset,
		Limit:      limit,
		Gene:       gene,
		Threshold:  threshold,
	}, nil
}

// GeneStats contains statistics for a gene.
type GeneStats struct {
	Gene           string  `json:"gene"`
	Index          int     `json:"index"`
	ExpressingBins int     `json:"expressing_bins"`
	TotalBins      int     `json:"total_bins"`
	TotalCells     int     `json:"total_cells"`
	MeanExpression float32 `json:"mean_expression"`
	MaxExpression  float32 `json:"max_expression"`
	P80Expression  float32 `json:"p80_expression"`
}

// GetGeneStats returns statistics for a gene at a specific zoom level.
func (r *Reader) GetGeneStats(gene string, zoom int) (*GeneStats, error) {
	geneIdx, ok := r.metadata.GeneIndex[gene]
	if !ok {
		return nil, fmt.Errorf("gene not found: %s", gene)
	}

	// Validate zoom level
	if zoom < 0 || zoom >= r.metadata.ZoomLevels {
		return nil, fmt.Errorf("invalid zoom level: %d (valid: 0-%d)", zoom, r.metadata.ZoomLevels-1)
	}

	// Load expression at specified zoom level
	expression, err := r.GetExpression(zoom, geneIdx, "mean")
	if err != nil {
		return nil, err
	}

	// Load bins for cell counts at specified zoom level
	bins, err := r.GetBins(zoom, 0, 0)
	if err != nil {
		return nil, err
	}

	// Calculate statistics
	var sum, maxExpr float32
	var expressing, totalCells int
	totalBins := len(expression)
	expressingValues := expression[:0]

	for i, v := range expression {
		if v > 0 {
			expressing++
			sum += v
			if v > maxExpr {
				maxExpr = v
			}
			expressingValues = append(expressingValues, v)
		}
		if i < len(bins) {
			totalCells += int(bins[i].CellCount)
		}
	}

	var mean float32
	if expressing > 0 {
		mean = sum / float32(expressing)
	}

	var p80 float32
	if len(expressingValues) > 0 {
		sort.Slice(expressingValues, func(i, j int) bool { return expressingValues[i] < expressingValues[j] })
		n := len(expressingValues)
		// idx = ceil(0.80*n) - 1, computed with integers.
		idx := (80*n+99)/100 - 1
		if idx < 0 {
			idx = 0
		} else if idx >= n {
			idx = n - 1
		}
		p80 = expressingValues[idx]
	}

	return &GeneStats{
		Gene:           gene,
		Index:          geneIdx,
		ExpressingBins: expressing,
		TotalBins:      totalBins,
		TotalCells:     totalCells,
		MeanExpression: mean,
		MaxExpression:  maxExpr,
		P80Expression:  p80,
	}, nil
}

// Close releases resources.
func (r *Reader) Close() {
	if r.decoder != nil {
		r.decoder.Close()
	}
}

func float32frombits(b uint32) float32 {
	return *(*float32)(unsafe.Pointer(&b))
}

func mathFloat32bits(f float32) uint32 {
	return *(*uint32)(unsafe.Pointer(&f))
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func ceilDiv(a, b int) int {
	if b <= 0 {
		return 0
	}
	return (a + b - 1) / b
}
