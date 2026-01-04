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


