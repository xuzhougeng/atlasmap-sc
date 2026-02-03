import L from 'leaflet';
import type { ApiClient } from '../api/client';
import type { MapController } from '../map/MapController';

export interface LassoSelectionPanelConfig {
    getCategoryColumn: () => string | null;
    getCategoryFilter: () => string[] | null;
    onEnabledChange?: (enabled: boolean) => void;
}

interface SelectionStats {
    categoryColumn: string | null;
    totalSelected: number;
    truncated: boolean;
    categoryCounts: Array<{ value: string; count: number }>;
    selectedJoinIDs: number[];
}

export class LassoSelectionPanel {
    private container: HTMLElement;
    private mapController: MapController;
    private api: ApiClient;
    private config: LassoSelectionPanelConfig;
    private enabled = false;
    private drawing = false;
    private points: L.LatLng[] = [];
    private drawingLine: L.Polyline | null = null;
    private selectionPolygon: L.Polygon | null = null;
    private selectionStats: SelectionStats | null = null;
    private status: 'idle' | 'loading' | 'error' = 'idle';
    private statusMessage: string | null = null;
    private abortController: AbortController | null = null;

    private readonly onMapMouseDown = (e: L.LeafletEvent): void => {
        if (!this.enabled) return;
        if (this.drawing) return;
        const ev = e as unknown as L.LeafletMouseEvent;
        if (!ev.latlng) return;
        if (typeof (ev.originalEvent as MouseEvent | undefined)?.button === 'number') {
            const btn = (ev.originalEvent as MouseEvent).button;
            if (btn !== 0) return;
        }
        this.startDrawing(ev.latlng);
    };

    private readonly onMapMouseMove = (e: L.LeafletEvent): void => {
        if (!this.enabled) return;
        if (!this.drawing) return;
        const ev = e as unknown as L.LeafletMouseEvent;
        if (!ev.latlng) return;
        this.addPoint(ev.latlng);
    };

    private readonly onMapMouseUp = (): void => {
        if (!this.enabled) return;
        if (!this.drawing) return;
        void this.finishDrawing();
    };

    private readonly onWindowPointerUp = (): void => {
        if (!this.enabled) return;
        if (!this.drawing) return;
        void this.finishDrawing();
    };

    private readonly onKeyDown = (e: KeyboardEvent): void => {
        if (!this.enabled) return;
        if (e.key === 'Escape') {
            e.preventDefault();
            this.setEnabled(false);
        }
    };

    constructor(
        container: HTMLElement,
        mapController: MapController,
        api: ApiClient,
        config: LassoSelectionPanelConfig
    ) {
        this.container = container;
        this.mapController = mapController;
        this.api = api;
        this.config = config;

        const map = this.mapController.getMap();
        map.on('mousedown', this.onMapMouseDown);
        map.on('mousemove', this.onMapMouseMove);
        map.on('mouseup', this.onMapMouseUp);

        // Touch support (Leaflet maps touch to similar events, but ensure we capture end)
        map.on('touchstart', this.onMapMouseDown);
        map.on('touchmove', this.onMapMouseMove);
        map.on('touchend', this.onMapMouseUp);

        this.render();
    }

    isEnabled(): boolean {
        return this.enabled;
    }

    setEnabled(enabled: boolean): void {
        if (this.enabled === enabled) return;

        this.enabled = enabled;
        if (enabled) {
            this.mapController.setInteractionsEnabled(false);
            window.addEventListener('keydown', this.onKeyDown);
        } else {
            this.cancelInFlight();
            this.cancelDrawing();
            this.mapController.setInteractionsEnabled(true);
            window.removeEventListener('keydown', this.onKeyDown);
        }

        this.config.onEnabledChange?.(this.enabled);
        this.render();
    }

    toggle(): void {
        this.setEnabled(!this.enabled);
    }

    clearSelection(): void {
        this.cancelInFlight();
        this.cancelDrawing();
        if (this.selectionPolygon) {
            this.mapController.getMap().removeLayer(this.selectionPolygon);
            this.selectionPolygon = null;
        }
        this.selectionStats = null;
        this.status = 'idle';
        this.statusMessage = null;
        this.render();
    }

    private startDrawing(start: L.LatLng): void {
        this.cancelInFlight();
        this.cancelDrawing();

        this.drawing = true;
        this.points = [start];

        window.addEventListener('mouseup', this.onWindowPointerUp);
        window.addEventListener('touchend', this.onWindowPointerUp);

        const map = this.mapController.getMap();
        this.drawingLine = L.polyline(this.points, {
            color: '#4f9cff',
            weight: 2,
            opacity: 0.9,
        }).addTo(map);

        this.status = 'idle';
        this.statusMessage = null;
        this.render();
    }

    private addPoint(point: L.LatLng): void {
        if (this.points.length > 2000) return;
        const map = this.mapController.getMap();
        const last = this.points[this.points.length - 1];
        const p1 = map.latLngToContainerPoint(last);
        const p2 = map.latLngToContainerPoint(point);
        if (p1.distanceTo(p2) < 2) return;

        this.points.push(point);
        this.drawingLine?.setLatLngs(this.points);
    }

    private async finishDrawing(): Promise<void> {
        this.drawing = false;
        window.removeEventListener('mouseup', this.onWindowPointerUp);
        window.removeEventListener('touchend', this.onWindowPointerUp);

        const map = this.mapController.getMap();
        if (this.drawingLine) {
            map.removeLayer(this.drawingLine);
            this.drawingLine = null;
        }

        if (this.points.length < 3) {
            this.points = [];
            this.statusMessage = 'Selection cancelled (too small).';
            this.render();
            return;
        }

        if (this.selectionPolygon) {
            map.removeLayer(this.selectionPolygon);
            this.selectionPolygon = null;
        }

        this.selectionPolygon = L.polygon(this.points, {
            color: '#4f9cff',
            weight: 2,
            opacity: 0.9,
            fillColor: '#4f9cff',
            fillOpacity: 0.12,
        }).addTo(map);

        this.status = 'loading';
        this.statusMessage = 'Querying cells…';
        this.render();

        const { bounds, polygon } = this.getPolygonBoundsAndPoints(this.points);
        const categoryColumn = this.config.getCategoryColumn();
        const categoryFilter = this.config.getCategoryFilter();

        try {
            this.cancelInFlight();
            this.abortController = new AbortController();
            const result = await this.api.getSomaCellsInBounds(
                bounds,
                {
                    category: categoryColumn ?? undefined,
                    categories: categoryFilter ?? null,
                    limit: 50000,
                    seed: 0,
                },
                this.abortController.signal
            );

            const selected = result.cells.filter((cell) =>
                LassoSelectionPanel.pointInPolygon(cell.x, cell.y, polygon)
            );

            const counts = new Map<string, number>();
            const joinIDs: number[] = [];
            for (const cell of selected) {
                joinIDs.push(cell.joinid);
                const value = (cell.category ?? '(unknown)').trim() || '(unknown)';
                counts.set(value, (counts.get(value) ?? 0) + 1);
            }

            const categoryCounts = Array.from(counts.entries())
                .map(([value, count]) => ({ value, count }))
                .sort((a, b) => b.count - a.count || a.value.localeCompare(b.value));

            this.selectionStats = {
                categoryColumn,
                totalSelected: selected.length,
                truncated: result.truncated,
                categoryCounts,
                selectedJoinIDs: joinIDs,
            };

            if (result.truncated) {
                this.statusMessage = 'Selection is based on a sampled subset (hit cell query limit).';
            } else {
                this.statusMessage = null;
            }
            this.status = 'idle';
            this.render();
        } catch (error) {
            if ((error as any)?.name === 'AbortError') return;
            this.status = 'error';
            this.statusMessage = error instanceof Error ? error.message : 'Failed to query selection.';
            this.render();
        } finally {
            this.abortController = null;
        }
    }

    private cancelInFlight(): void {
        if (this.abortController) {
            this.abortController.abort();
            this.abortController = null;
        }
    }

    private cancelDrawing(): void {
        this.drawing = false;
        window.removeEventListener('mouseup', this.onWindowPointerUp);
        window.removeEventListener('touchend', this.onWindowPointerUp);
        this.points = [];
        const map = this.mapController.getMap();
        if (this.drawingLine) {
            map.removeLayer(this.drawingLine);
            this.drawingLine = null;
        }
    }

    private render(): void {
        const modeText = this.enabled ? 'ON' : 'OFF';
        const selected = this.selectionStats;
        const categoryColumn = selected?.categoryColumn ?? this.config.getCategoryColumn();

        const cellsText = selected
            ? `${selected.totalSelected.toLocaleString()}${selected.truncated ? ' (sampled)' : ''}`
            : '-';

        const statusText =
            this.status === 'loading' ? 'Loading…' : this.status === 'error' ? 'Error' : '';

        const hint = this.enabled
            ? 'Drag on the map to draw a lasso. Release to finish. Press Esc to exit.'
            : 'Click the lasso tool to start selecting a region.';

        const listHtml = selected
            ? `
                <div class="selection-list" role="list">
                    ${selected.categoryCounts.length === 0
                        ? `<div class="selection-item"><span class="selection-item__label">(no category data)</span><span class="selection-item__count">-</span></div>`
                        : selected.categoryCounts
                              .map(
                                  (row) => `
                                    <div class="selection-item" role="listitem">
                                        <span class="selection-item__label" title="${escapeHtml(row.value)}">${escapeHtml(row.value)}</span>
                                        <span class="selection-item__count">${row.count.toLocaleString()}</span>
                                    </div>
                                `
                              )
                              .join('')}
                </div>
            `
            : '';

        this.container.innerHTML = `
            <div class="selection-panel">
                <div class="selection-hint">${escapeHtml(hint)}</div>
                <div class="selection-meta">
                    <div><strong>Mode</strong>: ${escapeHtml(modeText)} ${statusText ? ` · ${escapeHtml(statusText)}` : ''}</div>
                    <div><strong>Category</strong>: ${escapeHtml(categoryColumn ?? '-')}</div>
                    <div><strong>Cells</strong>: ${escapeHtml(cellsText)}</div>
                </div>
                <div class="selection-actions">
                    <button id="btn-selection-clear" class="btn-primary" ${selected ? '' : 'disabled'}>
                        Clear
                    </button>
                    <button
                        id="btn-selection-download"
                        class="btn-link"
                        type="button"
                        disabled
                        title="TODO: allow download of selected cells"
                    >
                        Download selected cells (TODO)
                    </button>
                </div>
                ${this.statusMessage ? `<div class="selection-hint">${escapeHtml(this.statusMessage)}</div>` : ''}
                ${listHtml}
            </div>
        `;

        const clearBtn = this.container.querySelector('#btn-selection-clear') as HTMLButtonElement | null;
        clearBtn?.addEventListener('click', () => this.clearSelection());
    }

    private getPolygonBoundsAndPoints(latlngs: L.LatLng[]): {
        bounds: { minX: number; minY: number; maxX: number; maxY: number };
        polygon: Array<{ x: number; y: number }>;
    } {
        let minX = Number.POSITIVE_INFINITY;
        let minY = Number.POSITIVE_INFINITY;
        let maxX = Number.NEGATIVE_INFINITY;
        let maxY = Number.NEGATIVE_INFINITY;

        const polygon: Array<{ x: number; y: number }> = [];
        for (const ll of latlngs) {
            const x = ll.lng;
            const y = ll.lat;
            polygon.push({ x, y });
            if (x < minX) minX = x;
            if (x > maxX) maxX = x;
            if (y < minY) minY = y;
            if (y > maxY) maxY = y;
        }

        return {
            bounds: { minX, minY, maxX, maxY },
            polygon,
        };
    }

    private static pointInPolygon(
        x: number,
        y: number,
        polygon: Array<{ x: number; y: number }>
    ): boolean {
        let inside = false;
        for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
            const xi = polygon[i].x;
            const yi = polygon[i].y;
            const xj = polygon[j].x;
            const yj = polygon[j].y;

            const intersects = (yi > y) !== (yj > y) && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi;
            if (intersects) inside = !inside;
        }
        return inside;
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
