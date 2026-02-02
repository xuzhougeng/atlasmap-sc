// Map Controller - Wraps Leaflet for UMAP visualization

import L from 'leaflet';
import { CategoryPostTileLayer } from './CategoryPostTileLayer';
import type { CategoryCentroidItem } from '../api/client';

export type ColorMode = 'expression' | 'category' | 'default';

export interface MapConfig {
    apiUrl: string;
    coord?: string;
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
    private container: HTMLElement;
    private config: MapConfig;
    private expressionLayer: L.TileLayer | null = null;
    private categoryLayer: L.TileLayer | null = null;
    private categoryLabelLayer: L.LayerGroup;
    private categoryLabelPane: HTMLElement;
    private bounds: L.LatLngBounds;
    private maxBounds: L.LatLngBounds;
    private panLocked: boolean = true;
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
    private categoryLabelsVisible: boolean = true;
    private pointSize: number = 1.0;

    constructor(container: HTMLElement, config: MapConfig) {
        this.container = container;
        this.config = config;

        // Calculate bounds from config (CRS.Simple expects [lat, lng] = [y, x])
        this.bounds = L.latLngBounds(
            [config.bounds.min_y, config.bounds.min_x],
            [config.bounds.max_y, config.bounds.max_x]
        );
        // Use a slightly padded bounds to allow a small amount of slack when panning is locked.
        this.maxBounds = L.latLngBounds(this.bounds.getSouthWest(), this.bounds.getNorthEast()).pad(0.1);

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
        this.setPanLocked(true);

        // Note: No base tile layer is created here.
        // Tiles are loaded on-demand via setCategoryColumn or setExpressionGene.

        const labelPane = this.map.createPane('categoryLabels');
        labelPane.style.zIndex = '650';
        labelPane.style.pointerEvents = 'none';
        this.categoryLabelPane = labelPane;

        this.categoryLabelLayer = L.layerGroup().addTo(this.map);
        this.updateCategoryLabelStyle();

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
            this.updateCategoryLabelStyle();
            this.renderCategoryLabels();
        });
    }

    setCategoryLabelsVisible(visible: boolean): void {
        this.categoryLabelsVisible = visible;
        this.categoryLabelPane.style.display = visible ? '' : 'none';
        if (!visible) {
            this.categoryLabelLayer.clearLayers();
            return;
        }
        this.updateCategoryLabelStyle();
        this.renderCategoryLabels();
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

        if (!this.categoryLabelsVisible) return;
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

    private updateCategoryLabelStyle(): void {
        const zoom = this.map.getZoom();
        const baseSizePx = 12;
        const zoom0 = this.config.initialZoom ?? 0;
        const scale = Math.pow(2, (zoom - zoom0) * 0.15);
        const sizePx = Math.max(11, Math.min(26, baseSizePx * scale));
        this.container.style.setProperty('--category-label-font-size', `${sizePx.toFixed(2)}px`);
    }

    private formatPointSize(): string {
        return this.pointSize.toFixed(3);
    }

    setPointSize(pointSize: number): void {
        const clamped = Math.min(5.0, Math.max(0.1, pointSize));
        const quantized = Math.round(clamped * 1000) / 1000;
        if (quantized === this.pointSize) return;

        this.pointSize = quantized;

        if (this.expressionLayer && this.currentExpressionGene) {
            this.expressionLayer.setUrl(this.buildExpressionTileUrl(this.currentExpressionGene));
        }

        if (this.categoryLayer && this.currentCategory) {
            this.categoryLayer.setUrl(this.buildCategoryTileUrl(this.currentCategory));
        }
    }

    getPointSize(): number {
        return this.pointSize;
    }

    private buildExpressionTileUrl(gene: string): string {
        const query = new URLSearchParams();
        if (this.config.coord) query.set('coord', this.config.coord);
        query.set('colormap', this.currentExpressionColorScale ?? 'viridis');
        query.set('point_size', this.formatPointSize());
        if (this.currentExpressionRange !== null) {
            query.set('min', String(this.currentExpressionRange.min));
            query.set('max', String(this.currentExpressionRange.max));
        }
        const geneParam = encodeURIComponent(gene);
        return `${this.config.apiUrl}/tiles/{z}/{x}/{y}/expression/${geneParam}.png?${query.toString()}`;
    }

    private buildCategoryTileUrl(column: string): string {
        const query = new URLSearchParams();
        if (this.config.coord) query.set('coord', this.config.coord);
        query.set('point_size', this.formatPointSize());
        const columnParam = encodeURIComponent(column);
        return `${this.config.apiUrl}/tiles/{z}/{x}/{y}/category/${columnParam}.png?${query.toString()}`;
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
        this.currentExpressionGene = gene;
        this.currentExpressionColorScale = colorScale;

        // Add expression-colored tile layer
        // Keep tiles visible beyond maxNativeZoom, but avoid requesting higher-zoom tiles.
        // Leaflet will reuse tiles at maxNativeZoom and scale them.
        // (cell canvas layer will take over visually at higher zooms)
        this.expressionLayer = L.tileLayer(
            this.buildExpressionTileUrl(gene),
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
        const tileUrl = this.buildCategoryTileUrl(column);

        // Add category-colored tile layer
        // Keep tiles visible beyond maxNativeZoom, but avoid requesting higher-zoom tiles.
        // Leaflet will reuse tiles at maxNativeZoom and scale them.
        // (cell canvas layer will take over visually at higher zooms)
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
            if (this.config.coord) query.set('coord', this.config.coord);
            query.set('colormap', this.currentExpressionColorScale ?? 'viridis');
            query.set('point_size', this.formatPointSize());
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
            const query = new URLSearchParams();
            if (this.config.coord) query.set('coord', this.config.coord);
            query.set('point_size', this.formatPointSize());
            const url = `${apiUrl}/tiles/${z}/${x}/${y}/category/${column}.png?${query.toString()}`;
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

    /**
     * Whether panning is constrained to the dataset bounds.
     */
    isPanLocked(): boolean {
        return this.panLocked;
    }

    /**
     * Lock/unlock panning to the dataset bounds.
     * When unlocked, users can freely drag the map view beyond the data extent.
     */
    setPanLocked(locked: boolean): void {
        this.panLocked = locked;
        if (locked) {
            this.map.setMaxBounds(this.maxBounds);
        } else {
            // Calling without arguments clears the maxBounds constraint.
            this.map.setMaxBounds();
        }
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
