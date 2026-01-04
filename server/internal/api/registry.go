package api

import (
	"github.com/soma-tiles/server/internal/service"
)

// DatasetInfo contains information about a dataset for the API response.
type DatasetInfo struct {
	ID   string `json:"id"`
	Name string `json:"name"`
}

// DatasetRegistry holds tile services for all configured datasets.
type DatasetRegistry struct {
	services       map[string]*service.TileService
	defaultDataset string
	datasetOrder   []string
	title          string
}

// NewDatasetRegistry creates a new dataset registry.
func NewDatasetRegistry(defaultDataset string, order []string, title string) *DatasetRegistry {
	return &DatasetRegistry{
		services:       make(map[string]*service.TileService),
		defaultDataset: defaultDataset,
		datasetOrder:   order,
		title:          title,
	}
}

// Register adds a tile service for a dataset.
func (r *DatasetRegistry) Register(datasetID string, svc *service.TileService) {
	r.services[datasetID] = svc
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
			ID:   id,
			Name: id,
		})
	}
	return infos
}
