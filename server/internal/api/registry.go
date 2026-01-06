package api

import (
	"os"

	"github.com/soma-tiles/server/internal/service"
)

// DatasetInfo contains information about a dataset for the API response.
type DatasetInfo struct {
	ID      string `json:"id"`
	Name    string `json:"name"`
	HasH5AD bool   `json:"has_h5ad,omitempty"`
}

// DatasetRegistry holds tile services for all configured datasets.
type DatasetRegistry struct {
	services       map[string]*service.TileService
	defaultDataset string
	datasetOrder   []string
	title          string
	h5adPaths      map[string]string
}

// NewDatasetRegistry creates a new dataset registry.
func NewDatasetRegistry(defaultDataset string, order []string, title string) *DatasetRegistry {
	return &DatasetRegistry{
		services:       make(map[string]*service.TileService),
		defaultDataset: defaultDataset,
		datasetOrder:   order,
		title:          title,
		h5adPaths:      make(map[string]string),
	}
}

// Register adds a tile service for a dataset.
func (r *DatasetRegistry) Register(datasetID string, svc *service.TileService) {
	r.services[datasetID] = svc
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

// Get returns the tile service for a dataset, or nil if not found.
func (r *DatasetRegistry) Get(datasetID string) *service.TileService {
	return r.services[datasetID]
}

// Default returns the default dataset's tile service.
func (r *DatasetRegistry) Default() *service.TileService {
	return r.services[r.defaultDataset]
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
			ID:      id,
			Name:    id,
			HasH5AD: r.HasH5AD(id),
		})
	}
	return infos
}
