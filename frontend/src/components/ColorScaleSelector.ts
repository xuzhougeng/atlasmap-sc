// Color Scale Selector Component
// Dropdown for selecting expression colormap

export interface ColorScaleSelectorConfig {
    onScaleChange: (scale: string) => void;
}

export class ColorScaleSelector {
    private container: HTMLElement;
    private config: ColorScaleSelectorConfig;
    private availableScales: string[] = [];
    private selectedScale: string = 'viridis';

    constructor(container: HTMLElement, config: ColorScaleSelectorConfig) {
        this.container = container;
        this.config = config;
        this.render();
    }

    private render(): void {
        this.container.innerHTML = `
            <div class="colormap-selector">
                <label>Colormap:</label>
                <select id="colormap-select">
                    ${this.availableScales
                        .map(scale => {
                            const label = scale.charAt(0).toUpperCase() + scale.slice(1);
                            return `<option value="${scale}" ${scale === this.selectedScale ? 'selected' : ''}>${label}</option>`;
                        })
                        .join('')}
                </select>
            </div>
        `;

        this.setupEvents();
    }

    private setupEvents(): void {
        const select = this.container.querySelector('#colormap-select') as HTMLSelectElement | null;
        if (!select) return;

        select.addEventListener('change', () => {
            this.selectedScale = select.value;
            this.config.onScaleChange(this.selectedScale);
        });
    }

    setScales(scales: string[]): void {
        this.availableScales = scales;
        if (scales.length > 0 && !scales.includes(this.selectedScale)) {
            this.selectedScale = scales[0];
        }
        this.render();
    }

    getSelectedScale(): string {
        return this.selectedScale;
    }

    setSelectedScale(scale: string): void {
        if (this.availableScales.includes(scale)) {
            this.selectedScale = scale;
            const select = this.container.querySelector('#colormap-select') as HTMLSelectElement | null;
            if (select) {
                select.value = scale;
            }
        }
    }
}

