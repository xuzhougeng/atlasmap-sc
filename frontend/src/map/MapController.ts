// Map Controller - Wraps Leaflet for UMAP visualization

import L from 'leaflet';

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
    private tileLayer: L.TileLayer;
    private expressionLayer: L.TileLayer | null = null;
    private bounds: L.LatLngBounds;

    constructor(container: HTMLElement, config: MapConfig) {
        this.config = config;

        // Calculate bounds from config
        const coordMax = config.bounds.max_x;
        this.bounds = L.latLngBounds(
            [0, 0],
            [coordMax, coordMax]
        );

        // Initialize map with CRS.Simple for UMAP coordinate space
        this.map = L.map(container, {
            crs: L.CRS.Simple,
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

        // Create base tile layer
        this.tileLayer = L.tileLayer(
            `${config.apiUrl}/tiles/{z}/{x}/{y}.png`,
            {
                tileSize: config.tileSize,
                noWrap: true,
                bounds: this.bounds,
                maxZoom: config.maxZoom,
                minZoom: 0,
            }
        ).addTo(this.map);

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
        // Remove existing expression layer
        this.clearExpression();

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
     * Reset view to initial state
     */
    resetView(): void {
        this.clearExpression();
        this.map.fitBounds(this.bounds);
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
        this.tileLayer.redraw();
        if (this.expressionLayer) {
            this.expressionLayer.redraw();
        }
    }

    /**
     * Get the underlying Leaflet map instance
     */
    getMap(): L.Map {
        return this.map;
    }
}
