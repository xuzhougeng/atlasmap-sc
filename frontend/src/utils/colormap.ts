export type RgbTuple = readonly [number, number, number];

export const COLORMAP_STOPS: Record<string, readonly RgbTuple[]> = {
    // Keep these in sync with `server/pkg/colormap/colormap.go`
    viridis: [
        [68, 1, 84],
        [72, 35, 116],
        [64, 67, 135],
        [52, 94, 141],
        [41, 120, 142],
        [32, 144, 140],
        [34, 167, 132],
        [68, 190, 112],
        [121, 209, 81],
        [189, 222, 38],
        [253, 231, 37],
    ],
    plasma: [
        [13, 8, 135],
        [75, 3, 161],
        [125, 3, 168],
        [168, 34, 150],
        [203, 70, 121],
        [229, 107, 93],
        [248, 148, 65],
        [253, 195, 40],
        [240, 249, 33],
    ],
    inferno: [
        [0, 0, 4],
        [40, 11, 84],
        [101, 21, 110],
        [159, 42, 99],
        [212, 72, 66],
        [245, 125, 21],
        [250, 193, 39],
        [252, 255, 164],
    ],
    magma: [
        [0, 0, 4],
        [28, 16, 68],
        [79, 18, 123],
        [129, 37, 129],
        [181, 54, 122],
        [229, 80, 100],
        [251, 135, 97],
        [254, 194, 135],
        [252, 253, 191],
    ],
    seurat: [
        [211, 211, 211],
        [255, 0, 0],
    ],
};

function toHex(n: number): string {
    const clamped = Math.max(0, Math.min(255, Math.round(n)));
    return clamped.toString(16).padStart(2, '0');
}

export function rgbToHex(rgb: RgbTuple): string {
    return `#${toHex(rgb[0])}${toHex(rgb[1])}${toHex(rgb[2])}`;
}

export function getColormapCssGradient(colormap: string): string | null {
    const stops = COLORMAP_STOPS[colormap];
    if (!stops || stops.length === 0) return null;
    if (stops.length === 1) return rgbToHex(stops[0]);

    const parts = stops.map((rgb, idx) => {
        const t = (idx / (stops.length - 1)) * 100;
        return `${rgbToHex(rgb)} ${t.toFixed(2)}%`;
    });
    return `linear-gradient(to top, ${parts.join(', ')})`;
}

