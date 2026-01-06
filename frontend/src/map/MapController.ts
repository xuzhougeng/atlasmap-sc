// Map Controller - Wraps Leaflet for UMAP visualization

import L from 'leaflet';
import { CategoryPostTileLayer } from './CategoryPostTileLayer';
import type { CategoryCentroidItem } from '../api/client';

export type ColorMode = 'expression' | 'category' | 'default';

export interface MapConfig {
    apiUrl: string;
    tileSize: number;
    // Leaflet zoom range (allows extra "magnification" zooms beyond native data)
    maxZoom: number;
    // Highest tile zoom the backend can serve; above this Leaflet scales tiles.
    maxNativeZoom: number;
    // Lowest Leaflet zoom level.
    minZoom?: number;
    // Initial Leaflet zoom level (defaults to fitBounds result).
    initialZoom?: number;
    bounds: {
        min_x: number;
        max_x: number;
        min_y: number;
        max_y: number;
    };
}

export interface TileRequest {
    url: string;
    init?: RequestInit;
}

export interface ExpressionColorRange {
    min: number;
    max: number;
}

export class MapController {
    private map: L.Map;
    private config: MapConfig;
    private expressionLayer: L.TileLayer | null = null;
    private categoryLayer: L.TileLayer | null = null;
    private categoryLabelLayer: L.LayerGroup;
    private bounds: L.LatLngBounds;
    private currentColorMode: ColorMode = 'default';
    private currentCategory: string | null = null;
    private currentExpressionGene: string | null = null;
    private currentExpressionColorScale: string | null = null;
    // null => auto range (use gene min/max on backend)
    private currentExpressionRange: ExpressionColorRange | null = null;
    // null => no filter (show all); [] => filter to none (show none)
    private selectedCategories: string[] | null = null;
    private categoryLabelItems: CategoryCentroidItem[] | null = null;
    private categoryLabelFilter: string[] | null = null;

    constructor(container: HTMLElement, config: MapConfig) {
        this.config = config;

        // Calculate bounds from config (CRS.Simple expects [lat, lng] = [y, x])
        this.bounds = L.latLngBounds(
            [config.bounds.min_y, config.bounds.min_x],
            [config.bounds.max_y, config.bounds.max_x]
        );

        // Initialize map with a CRS that keeps Y increasing upward (Cartesian/UMAP-like),
        // while shifting the origin so tile indices stay non-negative.
        const crs = L.Util.extend({}, L.CRS.Simple, {
            transformation: new L.Transformation(
                1,
                -config.bounds.min_x,
                -1,
                config.bounds.max_y
            ),
        }) as L.CRS;

        this.map = L.map(container, {
            crs,
            minZoom: config.minZoom ?? 0,
            maxZoom: config.maxZoom,
            zoomSnap: 1,
            zoomDelta: 1,
            attributionControl: false,
            zoomControl: true,
        });

        // Set initial view
        this.map.fitBounds(this.bounds);
        if (typeof config.initialZoom === 'number') {
            this.map.setZoom(config.initialZoom);
        }
        this.map.setMaxBounds(this.bounds.pad(0.1));

        // Note: No base tile layer is created here.
        // Tiles are loaded on-demand via setCategoryColumn or setExpressionGene.

        const labelPane = this.map.createPane('categoryLabels');
        labelPane.style.zIndex = '650';
        labelPane.style.pointerEvents = 'none';

        this.categoryLabelLayer = L.layerGroup().addTo(this.map);

        // Add zoom info display
        this.addZoomDisplay();

        // Handle map events
        this.setupEvents();
    }

    private addZoomDisplay(): void {
        const ZoomInfo = L.Control.extend({
            onAdd: (map: L.Map) => {
                const div = L.DomUtil.create('div', 'zoom-info');
                div.innerHTML = `Zoom: ${map.getZoom()}`;

                map.on('zoomend', () => {
                    div.innerHTML = `Zoom: ${map.getZoom()}`;
                });

                return div;
            },
        });

        new ZoomInfo({ position: 'bottomright' }).addTo(this.map);
    }

    private setupEvents(): void {
        // Click handler for cell/bin info
        this.map.on('click', (e: L.LeafletMouseEvent) => {
            const { lat, lng } = e.latlng;
            console.log(`Clicked at: (${lng.toFixed(2)}, ${lat.toFixed(2)})`);
            // Could fetch bin info here
        });

        // Zoom change handler
        this.map.on('zoomend', () => {
            console.log(`Zoom level: ${this.map.getZoom()}`);
            this.renderCategoryLabels();
        });
    }

    setCategoryLabelItems(items: CategoryCentroidItem[] | null): void {
        this.categoryLabelItems = items;
        this.renderCategoryLabels();
    }

    setCategoryLabelFilter(values: string[] | null): void {
        this.categoryLabelFilter = values;
        this.renderCategoryLabels();
    }

    private renderCategoryLabels(): void {
        this.categoryLabelLayer.clearLayers();

        if (!this.categoryLabelItems) return;
        if (this.categoryLabelFilter && this.categoryLabelFilter.length === 0) return;

        const filterSet =
            this.categoryLabelFilter === null
                ? null
                : new Set(this.categoryLabelFilter);

        for (const item of this.categoryLabelItems) {
            if (filterSet && !filterSet.has(item.value)) continue;
            if (item.x === null || item.y === null) continue;

            const text = item.value;
            const color = item.color || '#ffffff';

            const marker = L.marker([item.y, item.x], {
                interactive: false,
                keyboard: false,
                pane: 'categoryLabels',
                icon: L.divIcon({
                    className: 'category-label',
                    html: `<span class="category-label__text" style="color:${color}">${escapeHtml(text)}</span>`,
                }),
            });
            this.categoryLabelLayer.addLayer(marker);
        }
    }

    /**
     * Set expression gene for coloring
     */
    setExpressionGene(
        gene: string,
        colorScale: string = 'viridis',
        range?: ExpressionColorRange | null
    ): void {
        // Remove existing coloring layers
        this.clearExpression();
        this.clearCategoryLayer();

        if (typeof range !== 'undefined') {
            this.currentExpressionRange = range;
        }

        const query = new URLSearchParams();
        query.set('colormap', colorScale);
        if (this.currentExpressionRange !== null) {
            query.set('min', String(this.currentExpressionRange.min));
            query.set('max', String(this.currentExpressionRange.max));
        }

        // Add expression-colored tile layer
        this.expressionLayer = L.tileLayer(
            `${this.config.apiUrl}/tiles/{z}/{x}/{y}/expression/${gene}.png?${query.toString()}`,
            {
                tileSize: this.config.tileSize,
                noWrap: true,
                bounds: this.bounds,
                maxZoom: this.config.maxZoom,
                maxNativeZoom: this.config.maxNativeZoom,
                minZoom: this.config.minZoom ?? 0,
            }
        ).addTo(this.map);

        // Bring expression layer to front
        this.expressionLayer.bringToFront();
        this.currentColorMode = 'expression';
        this.currentExpressionGene = gene;
        this.currentExpressionColorScale = colorScale;
    }

    /**
     * Clear expression coloring
     */
    clearExpression(): void {
        if (this.expressionLayer) {
            this.map.removeLayer(this.expressionLayer);
            this.expressionLayer = null;
        }
        this.currentExpressionGene = null;
        this.currentExpressionColorScale = null;
    }

    /**
     * Set category column for coloring
     */
    setCategoryColumn(column: string, categoryFilter?: string[] | null): void {
        // Clear other coloring layers
        this.clearExpression();
        this.clearCategoryLayer();

        if (typeof categoryFilter !== 'undefined') {
            this.selectedCategories = categoryFilter;
        }
        const tileUrl = `${this.config.apiUrl}/tiles/{z}/{x}/{y}/category/${column}.png`;

        // Add category-colored tile layer
        const layerOptions: L.TileLayerOptions = {
            tileSize: this.config.tileSize,
            noWrap: true,
            bounds: this.bounds,
            maxZoom: this.config.maxZoom,
            maxNativeZoom: this.config.maxNativeZoom,
            minZoom: this.config.minZoom ?? 0,
        };

        this.categoryLayer =
            this.selectedCategories === null
                ? L.tileLayer(tileUrl, layerOptions)
                : new CategoryPostTileLayer(tileUrl, {
                      ...layerOptions,
                      categories: this.selectedCategories,
                  });
        this.categoryLayer.addTo(this.map);

        // Bring category layer to front
        this.categoryLayer.bringToFront();
        this.currentCategory = column;
        this.currentColorMode = 'category';
    }

    /**
     * Clear category coloring layer
     */
    clearCategoryLayer(): void {
        if (this.categoryLayer) {
            this.map.removeLayer(this.categoryLayer);
            this.categoryLayer = null;
        }
        this.currentCategory = null;
    }

    /**
     * Get current color mode
     */
    getColorMode(): ColorMode {
        return this.currentColorMode;
    }

    /**
     * Get current category column
     */
    getCurrentCategory(): string | null {
        return this.currentCategory;
    }

    /**
     * Get current expression gene (if in expression mode)
     */
    getCurrentExpressionGene(): string | null {
        return this.currentExpressionGene;
    }

    /**
     * Get current expression colormap (if in expression mode)
     */
    getCurrentExpressionColorScale(): string | null {
        return this.currentExpressionColorScale;
    }

    /**
     * Get configured tile size
     */
    getTileSize(): number {
        return this.config.tileSize;
    }

    /**
     * Get max native zoom supported by backend tiles
     */
    getMaxNativeZoom(): number {
        return this.config.maxNativeZoom;
    }

    /**
     * Get configured data bounds (numeric)
     */
    getDataBounds(): MapConfig['bounds'] {
        return this.config.bounds;
    }

    /**
     * Build a request for a concrete tile for the current display mode.
     *
     * Note: when category filtering is active, category tiles are fetched via POST and
     * the filter is sent in the request body (not in the URL).
     * Returns null when no tile layer is active (e.g., after resetView()).
     */
    getTileRequest(z: number, x: number, y: number): TileRequest | null {
        const apiUrl = this.config.apiUrl.replace(/\/$/, '');
        if (this.currentColorMode === 'expression' && this.currentExpressionGene) {
            const gene = encodeURIComponent(this.currentExpressionGene);
            const query = new URLSearchParams();
            query.set('colormap', this.currentExpressionColorScale ?? 'viridis');
            if (this.currentExpressionRange !== null) {
                query.set('min', String(this.currentExpressionRange.min));
                query.set('max', String(this.currentExpressionRange.max));
            }
            return {
                url: `${apiUrl}/tiles/${z}/${x}/${y}/expression/${gene}.png?${query.toString()}`,
            };
        }

        if (this.currentColorMode === 'category' && this.currentCategory) {
            const column = encodeURIComponent(this.currentCategory);
            const url = `${apiUrl}/tiles/${z}/${x}/${y}/category/${column}.png`;
            if (this.selectedCategories !== null) {
                return {
                    url,
                    init: {
                        method: 'POST',
                        body: JSON.stringify(this.selectedCategories),
                    },
                };
            }
            return { url };
        }

        return null;
    }

    /**
     * Build a concrete tile URL for the current display mode.
     * Returns null when no tile layer is active (e.g., after resetView()).
     */
    getTileUrl(z: number, x: number, y: number): string | null {
        return this.getTileRequest(z, x, y)?.url ?? null;
    }

    /**
     * Reset view to initial state
     */
    resetView(): void {
        this.clearExpression();
        this.clearCategoryLayer();
        this.map.fitBounds(this.bounds);
        if (typeof this.config.initialZoom === 'number') {
            this.map.setZoom(this.config.initialZoom);
        }
        this.currentColorMode = 'default';
    }

    /**
     * Pan to specific coordinates
     */
    panTo(x: number, y: number, zoom?: number): void {
        const targetZoom = zoom ?? this.map.getZoom();
        this.map.setView([y, x], targetZoom);
    }

    /**
     * Get current zoom level
     */
    getZoom(): number {
        return this.map.getZoom();
    }

    /**
     * Get current center
     */
    getCenter(): [number, number] {
        const center = this.map.getCenter();
        return [center.lat, center.lng];
    }

    /**
     * Refresh tiles (e.g., after filter change)
     */
    refreshTiles(): void {
        if (this.expressionLayer) {
            this.expressionLayer.redraw();
        }
        if (this.categoryLayer) {
            this.categoryLayer.redraw();
        }
    }

    /**
     * Update category filter and refresh tiles
     */
    updateCategoryFilter(categories: string[] | null): void {
        this.selectedCategories = categories;

        if (this.currentColorMode !== 'category' || !this.currentCategory) return;
        if (!this.categoryLayer) {
            this.setCategoryColumn(this.currentCategory);
            return;
        }

        const isPostLayer = this.categoryLayer instanceof CategoryPostTileLayer;
        const needsPostLayer = categories !== null;

        if (needsPostLayer) {
            if (isPostLayer) {
                (this.categoryLayer as CategoryPostTileLayer).setCategories(categories);
            } else {
                this.setCategoryColumn(this.currentCategory, categories);
            }
            return;
        }

        if (isPostLayer) {
            this.setCategoryColumn(this.currentCategory, null);
        } else {
            this.categoryLayer.redraw();
        }
    }

    /**
     * Get the underlying Leaflet map instance
     */
    getMap(): L.Map {
        return this.map;
    }
}

function escapeHtml(value: string): string {
    return value.replace(/[&<>"']/g, (ch) => {
        switch (ch) {
            case '&':
                return '&amp;';
            case '<':
                return '&lt;';
            case '>':
                return '&gt;';
            case '"':
                return '&quot;';
            case "'":
                return '&#39;';
            default:
                return ch;
        }
    });
}
