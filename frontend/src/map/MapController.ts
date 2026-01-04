// Map Controller - Wraps Leaflet for UMAP visualization

import L from 'leaflet';

export type ColorMode = 'expression' | 'category' | 'default';

export interface MapConfig {
    apiUrl: string;
    tileSize: number;
    maxZoom: number;
    bounds: {
        min_x: number;
        max_x: number;
        min_y: number;
        max_y: number;
    };
}

export class MapController {
    private map: L.Map;
    private config: MapConfig;
    private expressionLayer: L.TileLayer | null = null;
    private categoryLayer: L.TileLayer | null = null;
    private bounds: L.LatLngBounds;
    private currentColorMode: ColorMode = 'default';
    private currentCategory: string | null = null;

    constructor(container: HTMLElement, config: MapConfig) {
        this.config = config;

        // Calculate bounds from config (CRS.Simple expects [lat, lng] = [y, x])
        this.bounds = L.latLngBounds(
            [config.bounds.min_y, config.bounds.min_x],
            [config.bounds.max_y, config.bounds.max_x]
        );

        // Initialize map with a CRS that keeps Y increasing downward (image-like coordinates),
        // which also keeps tile y indices non-negative.
        const crs = L.Util.extend({}, L.CRS.Simple, {
            transformation: new L.Transformation(1, 0, 1, 0),
        }) as L.CRS;

        this.map = L.map(container, {
            crs,
            minZoom: 0,
            maxZoom: config.maxZoom,
            zoomSnap: 1,
            zoomDelta: 1,
            attributionControl: false,
            zoomControl: true,
        });

        // Set initial view
        this.map.fitBounds(this.bounds);
        this.map.setMaxBounds(this.bounds.pad(0.1));

        // Note: No base tile layer is created here.
        // Tiles are loaded on-demand via setCategoryColumn or setExpressionGene.

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
        });
    }

    /**
     * Set expression gene for coloring
     */
    setExpressionGene(gene: string, colorScale: string = 'viridis'): void {
        // Remove existing coloring layers
        this.clearExpression();
        this.clearCategoryLayer();

        // Add expression-colored tile layer
        this.expressionLayer = L.tileLayer(
            `${this.config.apiUrl}/tiles/{z}/{x}/{y}/expression/${gene}.png?colormap=${colorScale}`,
            {
                tileSize: this.config.tileSize,
                noWrap: true,
                bounds: this.bounds,
                maxZoom: this.config.maxZoom,
                minZoom: 0,
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
    }

    /**
     * Set category column for coloring
     */
    setCategoryColumn(column: string): void {
        // Clear other coloring layers
        this.clearExpression();
        this.clearCategoryLayer();

        // Add category-colored tile layer
        this.categoryLayer = L.tileLayer(
            `${this.config.apiUrl}/tiles/{z}/{x}/{y}/category/${column}.png`,
            {
                tileSize: this.config.tileSize,
                noWrap: true,
                bounds: this.bounds,
                maxZoom: this.config.maxZoom,
                minZoom: 0,
            }
        ).addTo(this.map);

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
     * Reset view to initial state
     */
    resetView(): void {
        this.clearExpression();
        this.clearCategoryLayer();
        this.map.fitBounds(this.bounds);
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
     * Get the underlying Leaflet map instance
     */
    getMap(): L.Map {
        return this.map;
    }
}
