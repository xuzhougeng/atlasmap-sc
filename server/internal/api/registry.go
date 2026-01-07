package api

import (
	"os"

	"github.com/soma-tiles/server/internal/service"
)

// DatasetInfo contains information about a dataset for the API response.
type DatasetInfo struct {
	ID        string `json:"id"`
	Name      string `json:"name"`
	HasH5AD   bool   `json:"has_h5ad,omitempty"`
	HasBlastP bool   `json:"has_blastp,omitempty"`
}

// DatasetRegistry holds tile services for all configured datasets.
type DatasetRegistry struct {
	datasets       map[string]*datasetServices
	defaultDataset string
	datasetOrder   []string
	title          string
	h5adPaths      map[string]string
	blastpPaths    map[string]string
}

type datasetServices struct {
	defaultCoord string
	services     map[string]*service.TileService // coordKey -> service
}

// NewDatasetRegistry creates a new dataset registry.
func NewDatasetRegistry(defaultDataset string, order []string, title string) *DatasetRegistry {
	return &DatasetRegistry{
		datasets:       make(map[string]*datasetServices),
		defaultDataset: defaultDataset,
		datasetOrder:   order,
		title:          title,
		h5adPaths:      make(map[string]string),
		blastpPaths:    make(map[string]string),
	}
}

// Register adds a tile service for a dataset.
func (r *DatasetRegistry) Register(datasetID string, svc *service.TileService) {
	r.RegisterCoord(datasetID, "", svc)
	r.SetDefaultCoord(datasetID, "")
}

// RegisterCoord adds a tile service for a dataset + coordinate system.
func (r *DatasetRegistry) RegisterCoord(datasetID, coordKey string, svc *service.TileService) {
	ds := r.datasets[datasetID]
	if ds == nil {
		ds = &datasetServices{
			services: make(map[string]*service.TileService),
		}
		r.datasets[datasetID] = ds
	}
	ds.services[coordKey] = svc
	if ds.defaultCoord == "" {
		ds.defaultCoord = coordKey
	}
}

// SetDefaultCoord sets the default coordinate system for a dataset.
func (r *DatasetRegistry) SetDefaultCoord(datasetID, coordKey string) {
	ds := r.datasets[datasetID]
	if ds == nil {
		return
	}
	if _, ok := ds.services[coordKey]; ok {
		ds.defaultCoord = coordKey
	}
}

// SetH5ADPath sets the path to a dataset's downloadable H5AD file.
// An empty path disables H5AD download for the dataset.
func (r *DatasetRegistry) SetH5ADPath(datasetID string, path string) {
	if path == "" {
		delete(r.h5adPaths, datasetID)
		return
	}
	r.h5adPaths[datasetID] = path
}

// H5ADPath returns the configured H5AD file path for a dataset (may be empty).
func (r *DatasetRegistry) H5ADPath(datasetID string) string {
	return r.h5adPaths[datasetID]
}

// HasH5AD reports whether the dataset has a configured and existing H5AD file.
func (r *DatasetRegistry) HasH5AD(datasetID string) bool {
	path := r.h5adPaths[datasetID]
	if path == "" {
		return false
	}
	info, err := os.Stat(path)
	if err != nil {
		return false
	}
	return !info.IsDir()
}

// SetBlastPPath sets the BLASTP database path for a dataset.
// An empty path disables BLASTP for the dataset.
func (r *DatasetRegistry) SetBlastPPath(datasetID string, path string) {
	if path == "" {
		delete(r.blastpPaths, datasetID)
		return
	}
	r.blastpPaths[datasetID] = path
}

// BlastPPath returns the configured BLASTP database path for a dataset (may be empty).
func (r *DatasetRegistry) BlastPPath(datasetID string) string {
	return r.blastpPaths[datasetID]
}

// HasBlastP reports whether the dataset has a configured BLASTP database.
// Checks if at least one of the expected BLAST DB files exists (e.g., .phr, .pin, .psq).
func (r *DatasetRegistry) HasBlastP(datasetID string) bool {
	path := r.blastpPaths[datasetID]
	if path == "" {
		return false
	}
	// Check for common BLAST DB file extensions
	for _, ext := range []string{".phr", ".pin", ".psq", ".pal"} {
		if _, err := os.Stat(path + ext); err == nil {
			return true
		}
	}
	return false
}

// BlastPPaths returns a map of all configured BLASTP paths (datasetID -> path).
func (r *DatasetRegistry) BlastPPaths() map[string]string {
	result := make(map[string]string, len(r.blastpPaths))
	for k, v := range r.blastpPaths {
		result[k] = v
	}
	return result
}

// Get returns the tile service for a dataset, or nil if not found.
func (r *DatasetRegistry) Get(datasetID string) *service.TileService {
	ds := r.datasets[datasetID]
	if ds == nil {
		return nil
	}
	if svc, ok := ds.services[ds.defaultCoord]; ok {
		return svc
	}
	// Fallback: return any registered service.
	for _, svc := range ds.services {
		return svc
	}
	return nil
}

// GetCoord returns the tile service for a dataset + coordinate system.
// If coordKey is empty or not found, it falls back to the dataset's default coordinate.
func (r *DatasetRegistry) GetCoord(datasetID, coordKey string) *service.TileService {
	ds := r.datasets[datasetID]
	if ds == nil {
		return nil
	}
	if coordKey != "" {
		if svc, ok := ds.services[coordKey]; ok {
			return svc
		}
	}
	if svc, ok := ds.services[ds.defaultCoord]; ok {
		return svc
	}
	for _, svc := range ds.services {
		return svc
	}
	return nil
}

// Default returns the default dataset's tile service.
func (r *DatasetRegistry) Default() *service.TileService {
	return r.Get(r.defaultDataset)
}

// DefaultDatasetID returns the default dataset ID.
func (r *DatasetRegistry) DefaultDatasetID() string {
	return r.defaultDataset
}

// DatasetIDs returns all dataset IDs in config order.
func (r *DatasetRegistry) DatasetIDs() []string {
	return r.datasetOrder
}

// Title returns the configured site title.
func (r *DatasetRegistry) Title() string {
	if r.title != "" {
		return r.title
	}
	return "SOMA-Tiles"
}

// Datasets returns dataset info for all registered datasets.
func (r *DatasetRegistry) Datasets() []DatasetInfo {
	infos := make([]DatasetInfo, 0, len(r.datasetOrder))
	for _, id := range r.datasetOrder {
		// Use the config ID as the display name (user-defined in server.yaml)
		infos = append(infos, DatasetInfo{
			ID:        id,
			Name:      id,
			HasH5AD:   r.HasH5AD(id),
			HasBlastP: r.HasBlastP(id),
		})
	}
	return infos
}
