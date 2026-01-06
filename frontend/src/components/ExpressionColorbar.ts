import { getColormapCssGradient } from '../utils/colormap';

type ExpressionRangeMode = 'auto' | 'manual';

export class ExpressionColorbar {
    private container: HTMLElement;
    private root: HTMLElement;
    private titleEl: HTMLElement;
    private subtitleEl: HTMLElement;
    private minEl: HTMLElement;
    private maxEl: HTMLElement;
    private barEl: HTMLElement;

    constructor(container: HTMLElement) {
        this.container = container;
        this.container.classList.add('expression-colorbar-container');

        this.container.innerHTML = `
            <div class="expression-colorbar hidden">
                <div class="expression-colorbar-header">
                    <div class="expression-colorbar-title"></div>
                    <div class="expression-colorbar-subtitle"></div>
                </div>
                <div class="expression-colorbar-scale">
                    <div class="expression-colorbar-label expression-colorbar-max">-</div>
                    <div class="expression-colorbar-bar"></div>
                    <div class="expression-colorbar-label expression-colorbar-min">-</div>
                </div>
            </div>
        `;

        const root = this.container.querySelector('.expression-colorbar') as HTMLElement | null;
        const titleEl = this.container.querySelector('.expression-colorbar-title') as HTMLElement | null;
        const subtitleEl = this.container.querySelector('.expression-colorbar-subtitle') as HTMLElement | null;
        const minEl = this.container.querySelector('.expression-colorbar-min') as HTMLElement | null;
        const maxEl = this.container.querySelector('.expression-colorbar-max') as HTMLElement | null;
        const barEl = this.container.querySelector('.expression-colorbar-bar') as HTMLElement | null;
        if (!root || !titleEl || !subtitleEl || !minEl || !maxEl || !barEl) {
            throw new Error('ExpressionColorbar: missing required elements');
        }

        this.root = root;
        this.titleEl = titleEl;
        this.subtitleEl = subtitleEl;
        this.minEl = minEl;
        this.maxEl = maxEl;
        this.barEl = barEl;
    }

    show(): void {
        this.root.classList.remove('hidden');
    }

    hide(): void {
        this.root.classList.add('hidden');
    }

    setGene(gene: string): void {
        this.titleEl.textContent = gene;
        this.titleEl.title = gene;
    }

    setColormap(colormap: string): void {
        const gradient = getColormapCssGradient(colormap) ?? getColormapCssGradient('viridis');
        if (gradient) {
            this.barEl.style.background = gradient;
        }
    }

    setRange(min: number | null, max: number | null, mode: ExpressionRangeMode): void {
        const modeLabel = mode === 'auto' ? 'Auto (p80)' : 'Manual';
        this.subtitleEl.textContent = modeLabel;
        this.subtitleEl.title = modeLabel;

        this.minEl.textContent = min === null ? '…' : this.formatValue(min);
        this.maxEl.textContent = max === null ? '…' : this.formatValue(max);
    }

    private formatValue(value: number): string {
        if (!Number.isFinite(value)) return '-';
        if (value === 0) return '0';

        const abs = Math.abs(value);
        let formatted: string;
        if (abs >= 1000) {
            formatted = value.toFixed(0);
        } else if (abs >= 100) {
            formatted = value.toFixed(1);
        } else if (abs >= 10) {
            formatted = value.toFixed(2);
        } else {
            formatted = value.toFixed(3);
        }

        // Trim trailing zeros and dot.
        formatted = formatted.replace(/\.?0+$/, '');
        return formatted;
    }
}
