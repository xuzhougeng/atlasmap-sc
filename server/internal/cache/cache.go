// Package cache provides caching for tiles and query results.
package cache

import (
	"context"
	"crypto/sha256"
	"encoding/hex"
	"fmt"
	"time"

	"github.com/allegro/bigcache/v3"
	lru "github.com/hashicorp/golang-lru/v2"
)

// Config contains cache configuration.
type Config struct {
	TileCacheSizeMB int
	TileTTL         time.Duration
	QueryCacheSize  int
}

// Manager manages tile and query caches.
type Manager struct {
	tileCache  *bigcache.BigCache
	queryCache *lru.Cache[string, []byte]
}

// NewManager creates a new cache manager.
func NewManager(cfg Config) (*Manager, error) {
	// Configure tile cache
	tileCacheConfig := bigcache.Config{
		Shards:             1024,
		LifeWindow:         cfg.TileTTL,
		CleanWindow:        cfg.TileTTL / 2,
		MaxEntriesInWindow: 100000,
		MaxEntrySize:       100 * 1024, // 100KB per tile
		HardMaxCacheSize:   cfg.TileCacheSizeMB,
		Verbose:            false,
	}

	tileCache, err := bigcache.New(context.Background(), tileCacheConfig)
	if err != nil {
		return nil, fmt.Errorf("failed to create tile cache: %w", err)
	}

	// Create query cache
	queryCache, err := lru.New[string, []byte](cfg.QueryCacheSize)
	if err != nil {
		return nil, fmt.Errorf("failed to create query cache: %w", err)
	}

	return &Manager{
		tileCache:  tileCache,
		queryCache: queryCache,
	}, nil
}

// GetTile retrieves a tile from cache.
func (m *Manager) GetTile(key string) ([]byte, bool) {
	data, err := m.tileCache.Get(key)
	if err != nil {
		return nil, false
	}
	return data, true
}

// SetTile stores a tile in cache.
func (m *Manager) SetTile(key string, data []byte) error {
	return m.tileCache.Set(key, data)
}

// GetQuery retrieves a query result from cache.
func (m *Manager) GetQuery(key string) ([]byte, bool) {
	return m.queryCache.Get(key)
}

// SetQuery stores a query result in cache.
func (m *Manager) SetQuery(key string, data []byte) {
	m.queryCache.Add(key, data)
}

// TileKey generates a cache key for a tile.
func TileKey(z, x, y int, filters map[string]interface{}) string {
	base := fmt.Sprintf("tile:%d/%d/%d", z, x, y)
	if len(filters) == 0 {
		return base
	}

	// Hash filters for cache key
	h := sha256.New()
	h.Write([]byte(base))
	for k, v := range filters {
		h.Write([]byte(fmt.Sprintf("%s=%v", k, v)))
	}
	return base + ":" + hex.EncodeToString(h.Sum(nil))[:16]
}

// ExpressionTileKey generates a cache key for an expression tile.
func ExpressionTileKey(z, x, y int, gene, colormap string) string {
	return fmt.Sprintf("expr:%d/%d/%d:%s:%s", z, x, y, gene, colormap)
}

// Stats returns cache statistics.
func (m *Manager) Stats() map[string]interface{} {
	return map[string]interface{}{
		"tile_cache_len":     m.tileCache.Len(),
		"tile_cache_cap":     m.tileCache.Capacity(),
		"query_cache_len":    m.queryCache.Len(),
	}
}

// Close closes the cache manager.
func (m *Manager) Close() error {
	return m.tileCache.Close()
}
