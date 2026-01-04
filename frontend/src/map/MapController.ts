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
    private baseOverlay: L.ImageOverlay;
    private expressionOverlay: L.ImageOverlay | null = null;
    private categoryOverlay: L.ImageOverlay | null = null;
    private bounds: L.LatLngBounds;
    private currentColorMode: ColorMode = 'default';
    private currentCategory: string | null = null;
    private currentExpressionGene: string | null = null;
    private currentExpressionColormap: string = 'viridis';

    constructor(container: HTMLElement, config: MapConfig) {
        this.config = config;

        // Calculate bounds from config (CRS.Simple expects [lat, lng] = [y, x])
        this.bounds = L.latLngBounds(
            [config.bounds.min_y, config.bounds.min_x],
            [config.bounds.max_y, config.bounds.max_x]
        );

        // Initialize map with a CRS that keeps Y increasing downward (more natural for image-like coordinates)
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

        // Create base overlay image (server currently renders full-view tiles at x=0,y=0)
        const initialZoom = this.getTileZoom();
        this.baseOverlay = L.imageOverlay(
            this.getBaseUrl(initialZoom),
            this.bounds,
            { interactive: false, zIndex: 0 }
        ).addTo(this.map);

        // Add zoom info display
        this.addZoomDisplay();

        // Handle map events
        this.setupEvents();
    }

    private getTileZoom(): number {
        const zoom = this.map.getZoom();
        // Leaflet zoom is integer here (zoomSnap=1), but clamp anyway
        return Math.max(0, Math.min(this.config.maxZoom, Math.round(zoom)));
    }

    private getBaseUrl(zoom: number): string {
        return `${this.config.apiUrl}/tiles/${zoom}/0/0.png`;
    }

    private getCategoryUrl(zoom: number, column: string): string {
        return `${this.config.apiUrl}/tiles/${zoom}/0/0/category/${column}.png`;
    }

    private getExpressionUrl(zoom: number, gene: string, colormap: string): string {
        return `${this.config.apiUrl}/tiles/${zoom}/0/0/expression/${gene}.png?colormap=${colormap}`;
    }

    private updateOverlaysForZoom(): void {
        const zoom = this.getTileZoom();
        this.baseOverlay.setUrl(this.getBaseUrl(zoom));

        if (this.categoryOverlay && this.currentCategory) {
            this.categoryOverlay.setUrl(this.getCategoryUrl(zoom, this.currentCategory));
        }

        if (this.expressionOverlay && this.currentExpressionGene) {
            this.expressionOverlay.setUrl(
                this.getExpressionUrl(zoom, this.currentExpressionGene, this.currentExpressionColormap)
            );
        }
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
            this.updateOverlaysForZoom();
        });
    }

    /**
     * Set expression gene for coloring
     */
    setExpressionGene(gene: string, colorScale: string = 'viridis'): void {
        // Remove existing coloring layers
        this.clearExpression();
        this.clearCategoryLayer();

        // Add expression-colored overlay
        const zoom = this.getTileZoom();
        this.currentExpressionGene = gene;
        this.currentExpressionColormap = colorScale;
        this.expressionOverlay = L.imageOverlay(
            this.getExpressionUrl(zoom, gene, colorScale),
            this.bounds,
            { interactive: false, zIndex: 2 }
        ).addTo(this.map);

        // Bring expression overlay to front
        this.expressionOverlay.bringToFront();
        this.currentColorMode = 'expression';
    }

    /**
     * Clear expression coloring
     */
    clearExpression(): void {
        if (this.expressionOverlay) {
            this.map.removeLayer(this.expressionOverlay);
            this.expressionOverlay = null;
        }
        this.currentExpressionGene = null;
    }

    /**
     * Set category column for coloring
     */
    setCategoryColumn(column: string): void {
        // Clear other coloring layers
        this.clearExpression();
        this.clearCategoryLayer();

        // Add category-colored overlay
        const zoom = this.getTileZoom();
        this.currentCategory = column;
        this.categoryOverlay = L.imageOverlay(
            this.getCategoryUrl(zoom, column),
            this.bounds,
            { interactive: false, zIndex: 1 }
        ).addTo(this.map);

        // Bring category overlay to front
        this.categoryOverlay.bringToFront();
        this.currentCategory = column;
        this.currentColorMode = 'category';
    }

    /**
     * Clear category coloring layer
     */
    clearCategoryLayer(): void {
        if (this.categoryOverlay) {
            this.map.removeLayer(this.categoryOverlay);
            this.categoryOverlay = null;
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
        this.updateOverlaysForZoom();
    }

    /**
     * Get the underlying Leaflet map instance
     */
    getMap(): L.Map {
        return this.map;
    }
}
