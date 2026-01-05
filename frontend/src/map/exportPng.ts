import L from 'leaflet';
import { MapController } from './MapController';

export interface DownloadPngProgress {
    done: number;
    total: number;
}

export interface DownloadPngOptions {
    maxDim?: number;
    concurrency?: number;
    onProgress?: (progress: DownloadPngProgress) => void;
}

export interface DownloadPngResult {
    filename: string;
    zoom: number;
    width: number;
    height: number;
    scale: number;
}

type DecodedImage = ImageBitmap | HTMLImageElement;

async function decodePng(blob: Blob): Promise<DecodedImage> {
    if ('createImageBitmap' in window) {
        return await createImageBitmap(blob);
    }

    return await new Promise((resolve, reject) => {
        const url = URL.createObjectURL(blob);
        const img = new Image();
        img.onload = () => {
            URL.revokeObjectURL(url);
            resolve(img);
        };
        img.onerror = () => {
            URL.revokeObjectURL(url);
            reject(new Error('Failed to decode PNG'));
        };
        img.src = url;
    });
}

function closeDecodedImage(img: DecodedImage): void {
    if (typeof (img as ImageBitmap).close === 'function') {
        (img as ImageBitmap).close();
    }
}

function canvasToBlob(canvas: HTMLCanvasElement): Promise<Blob> {
    return new Promise((resolve, reject) => {
        canvas.toBlob(
            (blob) => {
                if (!blob) {
                    reject(new Error('Failed to create PNG'));
                    return;
                }
                resolve(blob);
            },
            'image/png',
            1.0
        );
    });
}

function sanitizeFilenamePart(value: string): string {
    return value.replace(/[^a-zA-Z0-9._-]+/g, '_').replace(/^_+|_+$/g, '');
}

function buildFilename(mapController: MapController, zoom: number): string {
    const mode = mapController.getColorMode();
    if (mode === 'expression') {
        const gene = mapController.getCurrentExpressionGene() ?? 'gene';
        return `expression_${sanitizeFilenamePart(gene)}_z${zoom}.png`;
    }
    if (mode === 'category') {
        const col = mapController.getCurrentCategory() ?? 'category';
        return `category_${sanitizeFilenamePart(col)}_z${zoom}.png`;
    }
    return `tiles_z${zoom}.png`;
}

function triggerDownload(blob: Blob, filename: string): void {
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    window.setTimeout(() => URL.revokeObjectURL(url), 1000);
}

async function runWithConcurrency<T>(
    items: readonly T[],
    concurrency: number,
    fn: (item: T, idx: number) => Promise<void>
): Promise<void> {
    const limit = Math.max(1, Math.floor(concurrency));
    let nextIndex = 0;

    const workers = Array.from({ length: Math.min(limit, items.length) }, async () => {
        while (true) {
            const idx = nextIndex;
            nextIndex += 1;
            if (idx >= items.length) return;
            await fn(items[idx], idx);
        }
    });

    await Promise.all(workers);
}

function nextFrame(): Promise<void> {
    return new Promise((resolve) => window.requestAnimationFrame(() => resolve()));
}

/**
 * Download a stitched PNG for the current layer at the current (native) zoom.
 * The stitched image covers the full dataset bounds.
 */
export async function downloadPngAtCurrentZoom(
    mapController: MapController,
    options: DownloadPngOptions = {}
): Promise<DownloadPngResult> {
    const map = mapController.getMap();

    const tileSize = mapController.getTileSize();
    const maxNativeZoom = mapController.getMaxNativeZoom();
    const requestedZoom = map.getZoom();
    const zoom = Math.min(requestedZoom, maxNativeZoom);

    const sampleReq = mapController.getTileRequest(zoom, 0, 0);
    if (!sampleReq) {
        throw new Error('No active tile layer to export');
    }

    const dataBounds = mapController.getDataBounds();
    const topLeftPx = map.project(
        L.latLng(dataBounds.max_y, dataBounds.min_x),
        zoom
    );
    const bottomRightPx = map.project(
        L.latLng(dataBounds.min_y, dataBounds.max_x),
        zoom
    );

    const fullWidth = bottomRightPx.x - topLeftPx.x;
    const fullHeight = bottomRightPx.y - topLeftPx.y;
    if (fullWidth <= 0 || fullHeight <= 0) {
        throw new Error('Invalid bounds for export');
    }

    const maxDim = options.maxDim ?? 8192;
    const scale =
        maxDim > 0 ? Math.min(1, maxDim / Math.max(fullWidth, fullHeight)) : 1;

    const outWidth = Math.max(1, Math.round(fullWidth * scale));
    const outHeight = Math.max(1, Math.round(fullHeight * scale));

    const canvas = document.createElement('canvas');
    canvas.width = outWidth;
    canvas.height = outHeight;

    const ctx = canvas.getContext('2d');
    if (!ctx) {
        throw new Error('Canvas 2D context not available');
    }
    ctx.imageSmoothingEnabled = scale !== 1;

    const x0 = Math.max(0, Math.floor(topLeftPx.x / tileSize));
    const y0 = Math.max(0, Math.floor(topLeftPx.y / tileSize));
    const x1 = Math.max(0, Math.floor((bottomRightPx.x - 1) / tileSize));
    const y1 = Math.max(0, Math.floor((bottomRightPx.y - 1) / tileSize));

    const tiles: Array<{ x: number; y: number }> = [];
    for (let y = y0; y <= y1; y += 1) {
        for (let x = x0; x <= x1; x += 1) {
            tiles.push({ x, y });
        }
    }

    const total = tiles.length;
    let done = 0;
    const concurrency = options.concurrency ?? 12;

    const drawSize = Math.max(1, Math.round(tileSize * scale));

    await runWithConcurrency(tiles, concurrency, async ({ x, y }) => {
        const req = mapController.getTileRequest(zoom, x, y);
        if (!req) return;

        const resp = await fetch(req.url, req.init);
        if (!resp.ok) {
            throw new Error(`Failed to fetch tile: ${resp.status} ${resp.statusText}`);
        }

        const blob = await resp.blob();
        const img = await decodePng(blob);
        try {
            const dx = Math.round((x * tileSize - topLeftPx.x) * scale);
            const dy = Math.round((y * tileSize - topLeftPx.y) * scale);
            ctx.drawImage(img, dx, dy, drawSize, drawSize);
        } finally {
            closeDecodedImage(img);
        }

        done += 1;
        options.onProgress?.({ done, total });
        if (done % 200 === 0) {
            await nextFrame();
        }
    });

    const png = await canvasToBlob(canvas);
    const filename = buildFilename(mapController, zoom);
    triggerDownload(png, filename);

    return {
        filename,
        zoom,
        width: outWidth,
        height: outHeight,
        scale,
    };
}
