// Package config handles configuration loading for the SOMA-Tiles server.
package config

import (
	"fmt"
	"os"
	"sort"

	"gopkg.in/yaml.v3"
)

// Config represents the server configuration.
type Config struct {
	Server ServerConfig `yaml:"server"`
	Data   DataConfig   `yaml:"-"` // Custom unmarshal
	Cache  CacheConfig  `yaml:"cache"`
	Render RenderConfig `yaml:"render"`
}

// ServerConfig contains HTTP server settings.
type ServerConfig struct {
	Port        int      `yaml:"port"`
	CORSOrigins []string `yaml:"cors_origins"`
}

// DatasetConfig contains paths for a single dataset.
type DatasetConfig struct {
	ZarrPath string `yaml:"zarr_path"`
	SomaPath string `yaml:"soma_path"`
}

// DataConfig contains data source settings (supports legacy single or multi-dataset).
type DataConfig struct {
	Datasets       map[string]DatasetConfig
	DefaultDataset string
	DatasetOrder   []string // preserves YAML key order
}

// CacheConfig contains caching settings.
type CacheConfig struct {
	TileSizeMB     int `yaml:"tile_size_mb"`
	TileTTLMinutes int `yaml:"tile_ttl_minutes"`
}

// RenderConfig contains rendering settings.
type RenderConfig struct {
	TileSize        int    `yaml:"tile_size"`
	DefaultColormap string `yaml:"default_colormap"`
}

// rawConfig is used for initial YAML parsing.
type rawConfig struct {
	Server ServerConfig   `yaml:"server"`
	Data   yaml.Node      `yaml:"data"`
	Cache  CacheConfig    `yaml:"cache"`
	Render RenderConfig   `yaml:"render"`
}

// Load reads configuration from a YAML file.
func Load(path string) (*Config, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		// Return default config if file doesn't exist
		return DefaultConfig(), nil
	}

	var raw rawConfig
	if err := yaml.Unmarshal(data, &raw); err != nil {
		return nil, err
	}

	dataConfig, err := parseDataConfig(&raw.Data)
	if err != nil {
		return nil, fmt.Errorf("failed to parse data config: %w", err)
	}

	cfg := &Config{
		Server: raw.Server,
		Data:   *dataConfig,
		Cache:  raw.Cache,
		Render: raw.Render,
	}

	applyDefaults(cfg)
	return cfg, nil
}

// parseDataConfig handles both legacy (zarr_path/soma_path) and multi-dataset formats.
func parseDataConfig(node *yaml.Node) (*DataConfig, error) {
	if node == nil || node.Kind == 0 {
		// No data section; return default single dataset
		return &DataConfig{
			Datasets: map[string]DatasetConfig{
				"default": {
					ZarrPath: "./data/preprocessed/zarr/bins.zarr",
					SomaPath: "./data/soma",
				},
			},
			DefaultDataset: "default",
			DatasetOrder:   []string{"default"},
		}, nil
	}

	if node.Kind != yaml.MappingNode {
		return nil, fmt.Errorf("data section must be a mapping")
	}

	// Check if this is legacy format (has zarr_path or soma_path as direct keys)
	isLegacy := false
	for i := 0; i < len(node.Content); i += 2 {
		key := node.Content[i].Value
		if key == "zarr_path" || key == "soma_path" {
			isLegacy = true
			break
		}
	}

	if isLegacy {
		var ds DatasetConfig
		if err := node.Decode(&ds); err != nil {
			return nil, fmt.Errorf("failed to parse legacy data config: %w", err)
		}
		return &DataConfig{
			Datasets:       map[string]DatasetConfig{"default": ds},
			DefaultDataset: "default",
			DatasetOrder:   []string{"default"},
		}, nil
	}

	// Multi-dataset format: each key is a dataset ID
	datasets := make(map[string]DatasetConfig)
	var order []string

	for i := 0; i < len(node.Content); i += 2 {
		keyNode := node.Content[i]
		valNode := node.Content[i+1]

		datasetID := keyNode.Value
		if datasetID == "" {
			continue
		}

		var ds DatasetConfig
		if err := valNode.Decode(&ds); err != nil {
			return nil, fmt.Errorf("failed to parse dataset %q: %w", datasetID, err)
		}

		datasets[datasetID] = ds
		order = append(order, datasetID)
	}

	if len(datasets) == 0 {
		return nil, fmt.Errorf("no datasets defined in data section")
	}

	return &DataConfig{
		Datasets:       datasets,
		DefaultDataset: order[0], // First dataset in YAML order is default
		DatasetOrder:   order,
	}, nil
}

// DefaultConfig returns the default configuration.
func DefaultConfig() *Config {
	return &Config{
		Server: ServerConfig{
			Port:        8080,
			CORSOrigins: []string{"http://localhost:3000", "http://localhost:5173"},
		},
		Data: DataConfig{
			Datasets: map[string]DatasetConfig{
				"default": {
					ZarrPath: "./data/preprocessed/zarr/bins.zarr",
					SomaPath: "./data/soma",
				},
			},
			DefaultDataset: "default",
			DatasetOrder:   []string{"default"},
		},
		Cache: CacheConfig{
			TileSizeMB:     512,
			TileTTLMinutes: 10,
		},
		Render: RenderConfig{
			TileSize:        256,
			DefaultColormap: "viridis",
		},
	}
}

func applyDefaults(cfg *Config) {
	defaults := DefaultConfig()

	if cfg.Server.Port == 0 {
		cfg.Server.Port = defaults.Server.Port
	}
	if len(cfg.Server.CORSOrigins) == 0 {
		cfg.Server.CORSOrigins = defaults.Server.CORSOrigins
	}
	if cfg.Cache.TileSizeMB == 0 {
		cfg.Cache.TileSizeMB = defaults.Cache.TileSizeMB
	}
	if cfg.Cache.TileTTLMinutes == 0 {
		cfg.Cache.TileTTLMinutes = defaults.Cache.TileTTLMinutes
	}
	if cfg.Render.TileSize == 0 {
		cfg.Render.TileSize = defaults.Render.TileSize
	}
	if cfg.Render.DefaultColormap == "" {
		cfg.Render.DefaultColormap = defaults.Render.DefaultColormap
	}

	// Apply defaults to datasets with missing paths
	for id, ds := range cfg.Data.Datasets {
		if ds.ZarrPath == "" {
			ds.ZarrPath = defaults.Data.Datasets["default"].ZarrPath
			cfg.Data.Datasets[id] = ds
		}
	}
}

// DatasetIDs returns dataset IDs in config order.
func (c *DataConfig) DatasetIDs() []string {
	if len(c.DatasetOrder) > 0 {
		return c.DatasetOrder
	}
	// Fallback: sorted keys
	ids := make([]string, 0, len(c.Datasets))
	for id := range c.Datasets {
		ids = append(ids, id)
	}
	sort.Strings(ids)
	return ids
}
