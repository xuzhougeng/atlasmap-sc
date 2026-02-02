//go:build !soma

package soma

import (
	"fmt"
	"os"
)

// Reader is a stub when built without "-tags soma".
type Reader struct {
	experimentURI string
}

// NewReader creates a SOMA reader (stub). It still resolves and validates the experiment path,
// so config issues can be caught early, but all read methods return ErrUnsupported.
func NewReader(somaPath string) (*Reader, error) {
	uri, err := ResolveExperimentURI(somaPath)
	if err != nil {
		return nil, err
	}
	if _, statErr := os.Stat(uri); statErr != nil {
		return nil, fmt.Errorf("soma experiment not found at %s: %w", uri, statErr)
	}
	return &Reader{experimentURI: uri}, nil
}

func (r *Reader) Supported() bool { return false }

func (r *Reader) ExperimentURI() string { return r.experimentURI }

func (r *Reader) GeneJoinID(gene string) (int64, error) {
	return 0, ErrUnsupported
}

// ExpressionByCellJoinID reads expression values for one gene at given cell joinids.
// Returns a sparse map: cell_joinid -> value (only non-zero entries).
func (r *Reader) ExpressionByCellJoinID(gene string, cellJoinIDs []int64) (geneJoinID int64, values map[int64]float32, err error) {
	return 0, nil, ErrUnsupported
}

// GeneRange returns the min and max soma_joinid of genes in the var DataFrame.
func (r *Reader) GeneRange() (minID, maxID int64, err error) {
	return 0, 0, ErrUnsupported
}

// AllGenes returns a map of gene_id -> soma_joinid for all genes.
func (r *Reader) AllGenes() (map[string]int64, error) {
	return nil, ErrUnsupported
}

// ScanXForCells streams through ms/RNA/X/data for the given cell joinids (all genes).
func (r *Reader) ScanXForCells(cellJoinIDs []int64, onRow func(cell, gene int64, val float32)) error {
	return ErrUnsupported
}

// ObsColumns returns the list of attribute names in the obs DataFrame.
func (r *Reader) ObsColumns() ([]string, error) {
	return nil, ErrUnsupported
}

// ObsColumnValues returns unique values for a string column in obs.
func (r *Reader) ObsColumnValues(column string) ([]string, error) {
	return nil, ErrUnsupported
}

// ObsGroupIndex returns a map of column value -> cell joinids for a string column.
func (r *Reader) ObsGroupIndex(column string) (map[string][]int64, error) {
	return nil, ErrUnsupported
}

// SpatialIndex holds pre-computed spatial index for fast bbox queries.
type SpatialIndex struct {
	RenderZoom      int
	BinsPerAxis     int
	BinSize         float64
	CoordinateRange float64
	NCells          int
	X               []float32
	Y               []float32
	Offsets         []int
	CellIDs         []int64
}

// CellWithCoords holds a cell's joinID and its coordinates.
type CellWithCoords struct {
	JoinID int64
	X      float32
	Y      float32
}

// LoadSpatialIndex is a stub that returns ErrUnsupported.
func (r *Reader) LoadSpatialIndex(renderZoom int, coordinateRange float64) (*SpatialIndex, error) {
	return nil, ErrUnsupported
}

// QueryCellsInBounds is a stub that returns nil.
func (idx *SpatialIndex) QueryCellsInBounds(minX, minY, maxX, maxY float64, limit int) []int64 {
	return nil
}

// QueryCellsInBoundsWithCoords is a stub that returns nil.
func (idx *SpatialIndex) QueryCellsInBoundsWithCoords(minX, minY, maxX, maxY float64, limit int) []CellWithCoords {
	return nil
}

// GetCellCoordinates is a stub that returns empty maps.
func (idx *SpatialIndex) GetCellCoordinates(cellJoinIDs []int64) (map[int64]float32, map[int64]float32) {
	return nil, nil
}

// ReadObsCategoryColumn is a stub that returns ErrUnsupported.
func (r *Reader) ReadObsCategoryColumn(column string, cellJoinIDs []int64) (map[int64]string, error) {
	return nil, ErrUnsupported
}

// ReadObsCategoryColumnCached is a stub that returns ErrUnsupported.
func (r *Reader) ReadObsCategoryColumnCached(column string, cellJoinIDs []int64) (map[int64]string, error) {
	return nil, ErrUnsupported
}
