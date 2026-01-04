import L from 'leaflet';

export interface CategoryPostTileLayerOptions extends L.TileLayerOptions {
    categories: string[];
}

/**
 * Leaflet tile layer that fetches category tiles via POST to avoid large query strings.
 *
 * The request body is a JSON array of category values, e.g. ["T","B"].
 */
export class CategoryPostTileLayer extends L.TileLayer {
    private categories: string[];
    private abortControllers = new WeakMap<HTMLImageElement, AbortController>();

    constructor(urlTemplate: string, options: CategoryPostTileLayerOptions) {
        const { categories, ...tileOptions } = options;
        super(urlTemplate, tileOptions);
        this.categories = categories;

        this.on('tileunload', (e: L.TileEvent) => {
            const tile = e.tile as HTMLImageElement | undefined;
            if (!tile) return;
            const controller = this.abortControllers.get(tile);
            controller?.abort();
            this.abortControllers.delete(tile);
        });
    }

    setCategories(categories: string[]): void {
        this.categories = categories;
        this.redraw();
    }

    createTile(coords: L.Coords, done: L.DoneCallback): HTMLElement {
        const img = L.DomUtil.create('img', 'leaflet-tile') as HTMLImageElement;
        img.alt = '';
        img.setAttribute('role', 'presentation');
        img.decoding = 'async';

        const tileSize = this.getTileSize();
        img.width = tileSize.x;
        img.height = tileSize.y;

        const url = this.getTileUrl(coords);
        const controller = new AbortController();
        this.abortControllers.set(img, controller);

        fetch(url, {
            method: 'POST',
            body: JSON.stringify(this.categories),
            signal: controller.signal,
        })
            .then((resp) => {
                if (!resp.ok) {
                    throw new Error(`Tile request failed: ${resp.status} ${resp.statusText}`);
                }
                return resp.blob();
            })
            .then((blob) => {
                const objectUrl = URL.createObjectURL(blob);
                img.onload = () => {
                    URL.revokeObjectURL(objectUrl);
                    done(undefined, img);
                };
                img.onerror = () => {
                    URL.revokeObjectURL(objectUrl);
                    done(new Error('Failed to decode tile image'), img);
                };
                img.src = objectUrl;
            })
            .catch((err: unknown) => {
                if ((err as { name?: string }).name === 'AbortError') return;
                done(err as Error, img);
            });

        return img;
    }
}
