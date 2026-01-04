// Package soma provides minimal, read-only access to a TileDB-SOMA experiment using TileDB arrays.
//
// This is intentionally small: we only support what SOMA-Tiles needs today:
//   - map gene_id -> gene soma_joinid (from ms/RNA/var)
//   - read sparse X for (cells subset) x (one gene) (from ms/RNA/X/data)
package soma

import (
	"errors"
	"os"
	"path/filepath"
	"strings"
)

var (
	// ErrUnsupported indicates this binary was built without SOMA/TileDB support.
	ErrUnsupported = errors.New("soma support is not enabled in this build (build server with: go build -tags soma)")
)

// ResolveExperimentURI accepts either:
//   - /path/to/.../soma/experiment.soma
//   - /path/to/.../soma  (parent directory)
// and returns the experiment.soma path.
func ResolveExperimentURI(somaPath string) (string, error) {
	p := strings.TrimSpace(somaPath)
	if p == "" {
		return "", errors.New("empty soma_path")
	}
	p = os.ExpandEnv(p)
	p = filepath.Clean(p)

	// If user points directly to experiment.soma
	if strings.HasSuffix(p, ".soma") {
		return p, nil
	}
	// If user points to parent "soma/" dir
	return filepath.Join(p, "experiment.soma"), nil
}


