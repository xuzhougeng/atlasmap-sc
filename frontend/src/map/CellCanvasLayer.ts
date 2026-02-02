// Cell Canvas Layer - Renders individual cells when zoomed beyond maxNativeZoom

import L from 'leaflet';
import type { ApiClient, CellInfo } from '../api/client';
import { COLORMAP_STOPS, type RgbTuple } from '../utils/colormap';

export type CellColorMode = 'expression' | 'category';

export interface CellCanvasLayerConfig {
    api: ApiClient;
    map: L.Map;
    maxNativeZoom: number;
    bounds: {
        min_x: number;
        max_x: number;
        min_y: number;
        max_y: number;
    };
}

export interface CellColorOptions {
    mode: CellColorMode;
    gene?: string | null;
    categoryColumn?: string | null;
    categoryFilter?: string[] | null;
    categoryColors?: Record<string, string>;
    expressionColormap?: string;
    expressionMin?: number;
    expressionMax?: number;
}

export class CellCanvasLayer {
    private config: CellCanvasLayerConfig;
    private canvas: HTMLCanvasElement;
    private ctx: CanvasRenderingContext2D;
    private cells: CellInfo[] = [];
    private colorOptions: CellColorOptions = { mode: 'category' };

    // Request management
    private abortController: AbortController | null = null;
    private requestId = 0;
    private lastRequestedBounds: { minX: number; minY: number; maxX: number; maxY: number } | null = null;
    private isActive = false;
    private hasFetchedOnce = false;
    private pointSize = 3;
    // Global enable flag: when false, the layer never activates regardless of zoom
    private enabled = true;

    // Debounce
    private moveEndTimer: number | null = null;
    private fadeTimer: number | null = null;
    private dpr = 1;

    constructor(config: CellCanvasLayerConfig) {
        this.config = config;

        // Create canvas element
        this.canvas = document.createElement('canvas');
        this.canvas.style.position = 'absolute';
        this.canvas.style.top = '0';
        this.canvas.style.left = '0';
        this.canvas.style.pointerEvents = 'none';
        this.canvas.style.zIndex = '400'; // Above tiles but below controls
        this.canvas.style.opacity = '0';
        this.canvas.style.transition = 'opacity 140ms ease';
        this.canvas.style.willChange = 'opacity';
        this.canvas.className = 'cell-canvas-layer';

        const ctx = this.canvas.getContext('2d');
        if (!ctx) {
            throw new Error('Failed to get 2d context');
        }
        this.ctx = ctx;

        // Add canvas to map pane
        const pane = config.map.getPane('overlayPane');
        if (pane) {
            pane.appendChild(this.canvas);
        }

        // Set up event listeners
        this.setupEvents();

        // Initial check
        this.checkZoomLevel();
    }

    private updateCanvasPosition(): void {
        const topLeft = this.config.map.containerPointToLayerPoint([0, 0]);
        L.DomUtil.setPosition(this.canvas, topLeft);
    }

    private setupEvents(): void {
        // Use moveend for both pan and zoom changes
        this.config.map.on('moveend', this.onMoveEnd, this);
        this.config.map.on('resize', this.onResize, this);

        // Initial resize
        this.onResize();
    }

    private onMoveEnd = (): void => {
        // Clear any pending timer
        if (this.moveEndTimer !== null) {
            clearTimeout(this.moveEndTimer);
        }

        // Debounce to avoid multiple rapid requests
        this.moveEndTimer = window.setTimeout(() => {
            this.moveEndTimer = null;
            this.updateCanvasPosition();
            this.checkZoomLevel();
        }, 100);
    };

    private onResize = (): void => {
        const size = this.config.map.getSize();
        // Render at device pixel ratio for sharp points.
        this.dpr = Math.max(1, Math.min(3, window.devicePixelRatio || 1));
        this.canvas.style.width = `${size.x}px`;
        this.canvas.style.height = `${size.y}px`;
        this.canvas.width = Math.round(size.x * this.dpr);
        this.canvas.height = Math.round(size.y * this.dpr);
        this.ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);

        // Update position
        this.updateCanvasPosition();

        // Redraw if active
        if (this.isActive) {
            this.draw();
        }
    };

    private checkZoomLevel(): void {
        const zoom = this.config.map.getZoom();
        // Only activate if enabled AND zoom exceeds maxNativeZoom
        const shouldBeActive = this.enabled && zoom > this.config.maxNativeZoom;

        if (shouldBeActive && !this.isActive) {
            this.activate();
        } else if (!shouldBeActive && this.isActive) {
            this.deactivate();
        } else if (shouldBeActive) {
            // Already active, check if we need to fetch new data
            this.fetchCellsIfNeeded();
        }
    }

    private activate(): void {
        this.isActive = true;
        this.updateCanvasPosition();
        this.canvas.style.display = 'block';
        // Keep hidden until first successful fetch+draw, so tiles remain visible (no blank flash).
        this.canvas.style.opacity = '0';
        this.hasFetchedOnce = false;
        this.fetchCellsIfNeeded();
    }

    private deactivate(): void {
        this.isActive = false;
        // Fade out, then hide.
        this.canvas.style.opacity = '0';
        this.cells = [];
        this.lastRequestedBounds = null;
        this.hasFetchedOnce = false;

        // Cancel any pending request
        if (this.abortController) {
            this.abortController.abort();
            this.abortController = null;
        }

        if (this.fadeTimer !== null) {
            clearTimeout(this.fadeTimer);
        }
        this.fadeTimer = window.setTimeout(() => {
            this.fadeTimer = null;
            if (this.isActive) return;
            this.canvas.style.display = 'none';
            // Clear canvas (in CSS pixel space; ctx is already scaled by dpr).
            const size = this.config.map.getSize();
            this.ctx.clearRect(0, 0, size.x, size.y);
        }, 150);
    }

    private getCurrentBounds(): { minX: number; minY: number; maxX: number; maxY: number } {
        const leafletBounds = this.config.map.getBounds();
        // In Leaflet CRS.Simple, lat = y, lng = x
        return {
            minX: leafletBounds.getWest(),
            minY: leafletBounds.getSouth(),
            maxX: leafletBounds.getEast(),
            maxY: leafletBounds.getNorth(),
        };
    }

    private padBounds(
        bounds: { minX: number; minY: number; maxX: number; maxY: number },
        padding = 0.2
    ): { minX: number; minY: number; maxX: number; maxY: number } {
        const width = bounds.maxX - bounds.minX;
        const height = bounds.maxY - bounds.minY;
        const padX = width * padding;
        const padY = height * padding;
        return {
            minX: bounds.minX - padX,
            minY: bounds.minY - padY,
            maxX: bounds.maxX + padX,
            maxY: bounds.maxY + padY,
        };
    }

    private boundsContainedIn(
        inner: { minX: number; minY: number; maxX: number; maxY: number },
        outer: { minX: number; minY: number; maxX: number; maxY: number }
    ): boolean {
        return (
            inner.minX >= outer.minX &&
            inner.minY >= outer.minY &&
            inner.maxX <= outer.maxX &&
            inner.maxY <= outer.maxY
        );
    }

    private async fetchCellsIfNeeded(): Promise<void> {
        if (!this.isActive) return;

        const currentBounds = this.getCurrentBounds();

        // Check if current bounds are within last requested (padded) bounds
        if (this.lastRequestedBounds && this.boundsContainedIn(currentBounds, this.lastRequestedBounds)) {
            // No need to fetch, just redraw
            this.draw();
            return;
        }

        // If the very first fetch is still in-flight (often triggers backend spatial index build),
        // don't spam new requests; aborting won't cancel the backend warm-up anyway.
        if (this.abortController && !this.hasFetchedOnce) {
            return;
        }

        // Cancel previous request
        if (this.abortController) {
            this.abortController.abort();
        }
        this.abortController = new AbortController();

        const paddedBounds = this.padBounds(currentBounds, 0.3);
        const thisRequestId = ++this.requestId;

        // Calculate limit based on zoom level:
        // At lower zooms (closer to maxNativeZoom), we need fewer points since they're denser.
        // At higher zooms, we can show more points since they're sparser on screen.
        const zoom = this.config.map.getZoom();
        const zoomDiff = Math.max(0, zoom - this.config.maxNativeZoom);
        // Base limit of 5000, scaling up with zoom difference (up to 20000)
        const limit = Math.min(20000, 5000 * Math.pow(2, Math.min(2, zoomDiff)));

        try {
            const result = await this.config.api.getSomaCellsInBounds(
                paddedBounds,
                {
                    gene: this.colorOptions.mode === 'expression' ? this.colorOptions.gene ?? undefined : undefined,
                    category: this.colorOptions.mode === 'category' ? this.colorOptions.categoryColumn ?? undefined : undefined,
                    categories: this.colorOptions.categoryFilter ?? undefined,
                    limit,
                    seed: 0, // Fixed seed for deterministic downsampling (stable renders)
                },
                this.abortController.signal
            );

            // Check if this is still the current request
            if (thisRequestId !== this.requestId) return;

            // Log downsampling info for debugging
            if (result.truncated) {
                console.log(
                    `[CellCanvas] Downsampling applied: returned ${result.cells.length} cells ` +
                    `(limit=${limit}, seed=0, zoom=${zoom.toFixed(1)}, zoomDiff=${zoomDiff.toFixed(1)})`
                );
            } else {
                console.log(
                    `[CellCanvas] Fetched ${result.cells.length} cells (no downsampling needed, limit=${limit})`
                );
            }

            this.cells = result.cells;
            this.lastRequestedBounds = paddedBounds;
            this.hasFetchedOnce = true;
            this.draw();
            this.abortController = null;
        } catch (err) {
            if (err instanceof Error && err.name === 'AbortError') {
                // Request was cancelled, ignore
                this.abortController = null;
                return;
            }
            console.error('Failed to fetch cells:', err);
            this.abortController = null;
        }
    }

    private draw(): void {
        if (!this.isActive) return;

        // Keep canvas aligned to the current map transform (pan/zoom).
        this.updateCanvasPosition();

        const { ctx, canvas, cells } = this;
        const size = this.config.map.getSize();
        ctx.clearRect(0, 0, size.x, size.y);
        // White background to match expected map look.
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, size.x, size.y);

        // If we haven't fetched yet, keep hidden so scaled tiles show through.
        if (!this.hasFetchedOnce) return;

        // Fade in once we have at least one successful fetch+draw.
        if (this.canvas.style.opacity !== '1') {
            // Let the browser apply display/position before starting transition.
            requestAnimationFrame(() => {
                if (!this.isActive) return;
                this.canvas.style.opacity = '1';
            });
        }

        if (cells.length === 0) return;

        const map = this.config.map;

        // Compute point size based on zoom
        const zoom = map.getZoom();
        const zoomDiff = zoom - this.config.maxNativeZoom;
        // Scale point size with zoom: larger at higher zooms
        const baseSize = this.pointSize;
        const scaledSize = baseSize * Math.pow(1.5, zoomDiff);

        // Get color function based on mode
        const getColor = this.getColorFunction();

        // Draw cells
        for (const cell of cells) {
            // Convert data coordinates to screen coordinates
            // In Leaflet CRS.Simple: lat = y, lng = x
            const point = map.latLngToContainerPoint([cell.y, cell.x]);

            // Skip if outside canvas
            if (point.x < -scaledSize || point.x > canvas.width + scaledSize ||
                point.y < -scaledSize || point.y > canvas.height + scaledSize) {
                continue;
            }

            const color = getColor(cell);
            ctx.fillStyle = color;
            // Small points are much faster as squares; larger ones can be circles.
            if (scaledSize <= 1.25) {
                const s = Math.max(1, Math.round(scaledSize * 2));
                ctx.fillRect(Math.round(point.x), Math.round(point.y), s, s);
            } else {
                ctx.beginPath();
                ctx.arc(point.x, point.y, scaledSize, 0, Math.PI * 2);
                ctx.fill();
            }
        }
    }

    private getColorFunction(): (cell: CellInfo) => string {
        if (this.colorOptions.mode === 'expression') {
            return this.getExpressionColorFunction();
        } else {
            return this.getCategoryColorFunction();
        }
    }

    private getExpressionColorFunction(): (cell: CellInfo) => string {
        const colormap = this.colorOptions.expressionColormap ?? 'viridis';
        const stops = COLORMAP_STOPS[colormap] ?? COLORMAP_STOPS['viridis'];
        const min = this.colorOptions.expressionMin ?? 0;
        const max = this.colorOptions.expressionMax ?? 1;
        const range = max - min || 1;

        return (cell: CellInfo) => {
            const expr = cell.expression ?? 0;
            const t = Math.max(0, Math.min(1, (expr - min) / range));
            return this.interpolateColormap(stops, t);
        };
    }

    private getCategoryColorFunction(): (cell: CellInfo) => string {
        const colors = this.colorOptions.categoryColors ?? {};
        const defaultColor = '#808080';

        return (cell: CellInfo) => {
            const cat = cell.category;
            if (cat && colors[cat]) {
                return colors[cat];
            }
            return defaultColor;
        };
    }

    private interpolateColormap(stops: readonly RgbTuple[], t: number): string {
        if (stops.length === 0) return '#808080';
        if (stops.length === 1) return this.rgbToHex(stops[0]);
        if (t <= 0) return this.rgbToHex(stops[0]);
        if (t >= 1) return this.rgbToHex(stops[stops.length - 1]);

        const scaledT = t * (stops.length - 1);
        const i = Math.floor(scaledT);
        const frac = scaledT - i;

        const c0 = stops[i];
        const c1 = stops[Math.min(i + 1, stops.length - 1)];

        const r = Math.round(c0[0] + (c1[0] - c0[0]) * frac);
        const g = Math.round(c0[1] + (c1[1] - c0[1]) * frac);
        const b = Math.round(c0[2] + (c1[2] - c0[2]) * frac);

        return this.rgbToHex([r, g, b]);
    }

    private rgbToHex(rgb: RgbTuple | readonly [number, number, number]): string {
        const toHex = (n: number) => Math.max(0, Math.min(255, Math.round(n))).toString(16).padStart(2, '0');
        return `#${toHex(rgb[0])}${toHex(rgb[1])}${toHex(rgb[2])}`;
    }

    /**
     * Update color options and refresh display.
     */
    setColorOptions(options: Partial<CellColorOptions>): void {
        const modeChanged = options.mode !== undefined && options.mode !== this.colorOptions.mode;
        const geneChanged = options.gene !== undefined && options.gene !== this.colorOptions.gene;
        const categoryChanged = options.categoryColumn !== undefined && options.categoryColumn !== this.colorOptions.categoryColumn;
        const filterChanged = options.categoryFilter !== undefined;

        this.colorOptions = { ...this.colorOptions, ...options };

        if (this.isActive) {
            // If mode, gene, or category changed, we need to refetch
            if (modeChanged || geneChanged || categoryChanged || filterChanged) {
                this.lastRequestedBounds = null; // Force refetch
                void this.fetchCellsIfNeeded();
            } else {
                // Just color options changed, redraw
                this.draw();
            }
        }
    }

    /**
     * Set point size and redraw.
     */
    setPointSize(size: number): void {
        this.pointSize = Math.max(1, Math.min(10, size));
        if (this.isActive) {
            this.draw();
        }
    }

    /**
     * Enable or disable the cell canvas layer globally.
     * When disabled, the layer will deactivate and never activate regardless of zoom level.
     * When enabled, the layer will activate/deactivate based on zoom level as normal.
     */
    setEnabled(enabled: boolean): void {
        if (this.enabled === enabled) return;

        this.enabled = enabled;

        if (!enabled && this.isActive) {
            // Disable: force deactivation
            this.deactivate();
        } else if (enabled) {
            // Enable: check if we should activate based on current zoom
            this.checkZoomLevel();
        }
    }

    /**
     * Check if the layer is enabled.
     */
    isEnabled(): boolean {
        return this.enabled;
    }

    /**
     * Update the maxNativeZoom threshold.
     * This determines when the cell canvas layer activates (zoom > maxNativeZoom).
     */
    setMaxNativeZoom(maxNativeZoom: number): void {
        if (this.config.maxNativeZoom === maxNativeZoom) return;

        this.config.maxNativeZoom = maxNativeZoom;

        // Recheck zoom level with the new threshold
        this.checkZoomLevel();
    }

    /**
     * Get the current maxNativeZoom threshold.
     */
    getMaxNativeZoom(): number {
        return this.config.maxNativeZoom;
    }

    /**
     * Force refresh (refetch and redraw).
     */
    refresh(): void {
        if (this.isActive) {
            this.lastRequestedBounds = null;
            void this.fetchCellsIfNeeded();
        }
    }

    /**
     * Check if the layer is currently active (showing cells).
     */
    isLayerActive(): boolean {
        return this.isActive;
    }

    /**
     * Clean up resources.
     */
    destroy(): void {
        this.config.map.off('moveend', this.onMoveEnd, this);
        this.config.map.off('resize', this.onResize, this);

        if (this.moveEndTimer !== null) {
            clearTimeout(this.moveEndTimer);
        }
        if (this.fadeTimer !== null) {
            clearTimeout(this.fadeTimer);
        }

        if (this.abortController) {
            this.abortController.abort();
        }

        if (this.canvas.parentNode) {
            this.canvas.parentNode.removeChild(this.canvas);
        }
    }
}
