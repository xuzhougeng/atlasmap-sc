// Package config handles configuration loading for the SOMA-Tiles server.
package config

import (
	"os"

	"gopkg.in/yaml.v3"
)

// Config represents the server configuration.
type Config struct {
	Server ServerConfig `yaml:"server"`
	Data   DataConfig   `yaml:"data"`
	Cache  CacheConfig  `yaml:"cache"`
	Render RenderConfig `yaml:"render"`
}

// ServerConfig contains HTTP server settings.
type ServerConfig struct {
	Port        int      `yaml:"port"`
	CORSOrigins []string `yaml:"cors_origins"`
}

// DataConfig contains data source settings.
type DataConfig struct {
	ZarrPath string `yaml:"zarr_path"`
	SomaPath string `yaml:"soma_path"`
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

// Load reads configuration from a YAML file.
func Load(path string) (*Config, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		// Return default config if file doesn't exist
		return DefaultConfig(), nil
	}

	var cfg Config
	if err := yaml.Unmarshal(data, &cfg); err != nil {
		return nil, err
	}

	// Apply defaults for missing values
	applyDefaults(&cfg)

	return &cfg, nil
}

// DefaultConfig returns the default configuration.
func DefaultConfig() *Config {
	return &Config{
		Server: ServerConfig{
			Port:        8080,
			CORSOrigins: []string{"http://localhost:3000", "http://localhost:5173"},
		},
		Data: DataConfig{
			ZarrPath: "./data/preprocessed/zarr/bins.zarr",
			SomaPath: "./data/soma",
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
	if cfg.Data.ZarrPath == "" {
		cfg.Data.ZarrPath = defaults.Data.ZarrPath
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
}
